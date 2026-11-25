/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "onnx_model.hpp"

#include <onnx/shape_inference/implementation.h>

/// \cond INTERNAL
namespace casadi {

  // Read structural integer constants without passing their payload through DM.
  static bool integer_tensor(const onnx::TensorProto& tensor, std::vector<int64_t>& values) {
    if (tensor.data_type() != onnx::TensorProto::INT64 &&
        tensor.data_type() != onnx::TensorProto::INT32) return false;
    values.clear();
    if (tensor.has_raw_data()) {
      const auto& raw = tensor.raw_data();
      size_t width = tensor.data_type() == onnx::TensorProto::INT64 ? 8 : 4;
      casadi_assert(raw.size() % width == 0, "Invalid integer tensor raw_data size");
      for (size_t i = 0; i < raw.size(); i += width) {
        uint64_t bits = 0;
        for (size_t j = 0; j < width; ++j) {
          bits |= static_cast<uint64_t>(static_cast<unsigned char>(raw[i + j])) << (8 * j);
        }
        // ONNX raw_data is little-endian two's complement.
        if (width == 4) {
          values.push_back(static_cast<int64_t>(bits) - static_cast<int64_t>((bits >> 31) << 32));
        } else {
          values.push_back(bits >> 63 ? -1 - static_cast<int64_t>(~bits) :
            static_cast<int64_t>(bits));
        }
      }
    } else if (tensor.data_type() == onnx::TensorProto::INT64) {
      values.assign(tensor.int64_data().begin(), tensor.int64_data().end());
    } else {
      values.assign(tensor.int32_data().begin(), tensor.int32_data().end());
    }
    return true;
  }

  // Whether two stored shapes broadcast against each other (equal or 1 along each axis)
  static bool broadcastable(const MX& a, const MX& b) {
    return (a.size1() == b.size1() || a.size1() == 1 || b.size1() == 1)
        && (a.size2() == b.size2() || a.size2() == 1 || b.size2() == 1);
  }

  // Orient a vector operand v towards the other operand o of a binary op
  static void orient_vector(MX* v, const MX& o) {
    if (!v->is_vector()) return;
    // Two equal-length vectors in different orientations are the same rank-1 tensor;
    // settle on the column form, the convention for rank-1 initializers and inputs
    if (o.is_vector() && v->numel() == o.numel()) {
      if (v->size1() != o.size1()) *v = reshape(*v, v->numel(), 1);
      return;
    }
    if (broadcastable(*v, o)) return;
    if (v->numel() == o.size1()) *v = reshape(*v, o.size1(), 1);
    else if (v->numel() == o.size2()) *v = reshape(*v, 1, o.size2());
  }

  // ONNX (numpy-style) broadcasting of two operands in the 2-D stored representation.
  // A rank-1 tensor may be stored as a row or a column depending on where it came from
  // (Gather/Slice results are rows; initializers and graph inputs are columns), so a vector
  // that is not already broadcast-compatible is first oriented along the other operand's
  // matching dimension, e.g. a Gemm/Add bias [N] against a stored (N,M) result.
  static void broadcast_pair(MX* x, MX* y) {
    if (x->is_scalar() || y->is_scalar()) return;
    orient_vector(x, *y);
    orient_vector(y, *x);
    casadi_int nr = std::max(x->size1(), y->size1()), nc = std::max(x->size2(), y->size2());
    for (MX* v : {x, y}) {
      casadi_assert(broadcastable(*v, MX::zeros(nr, nc)), "ONNX import: cannot broadcast "
        + x->dim() + " with " + y->dim());
      if (v->size1() != nr || v->size2() != nc) *v = repmat(*v, nr / v->size1(), nc / v->size2());
    }
  }

  // Tensors of any rank (see Onnx::shapes_): an ONNX shape [d0, d1, ..., dk] is stored as a
  // (d1*...*dk) x d0 matrix, whose column-major order is the row-major ONNX element order.
  static std::pair<casadi_int, casadi_int> stored_size(const std::vector<casadi_int>& s) {
    if (s.empty()) return {1, 1};
    if (s.size() == 1) return {s[0], 1};
    casadi_int r = 1;
    for (size_t k = 1; k < s.size(); ++k) r *= s[k];
    return {r, s[0]};
  }

  static casadi_int shape_numel(const std::vector<casadi_int>& s) {
    casadi_int n = 1;
    for (casadi_int d : s) n *= d;
    return n;
  }

  static std::string shape_str(const std::vector<casadi_int>& s) {
    std::string r = "[";
    for (size_t k = 0; k < s.size(); ++k) r += (k ? ", " : "") + str(s[k]);
    return r + "]";
  }

  // Row-major strides of a shape
  static std::vector<casadi_int> shape_strides(const std::vector<casadi_int>& s) {
    std::vector<casadi_int> st(s.size(), 1);
    for (casadi_int k = static_cast<casadi_int>(s.size()) - 2; k >= 0; --k) {
      st[k] = st[k + 1] * s[k + 1];
    }
    return st;
  }

  // x in the stored form of shape s
  static MX as_shape(const MX& x, const std::vector<casadi_int>& s) {
    auto sz = stored_size(s);
    casadi_assert(x.numel() == sz.first * sz.second, "ONNX import: " + x.dim()
      + " value cannot have shape " + shape_str(s));
    return x.size1() == sz.first ? x : reshape(x, sz.first, sz.second);
  }

  // The elements of x at the given row-major positions, as a tensor of shape s
  static MX take_elements(const MX& x, const std::vector<casadi_int>& pos,
                          const std::vector<casadi_int>& s) {
    MX r;
    vec(densify(x)).get_nz(r, false, Matrix<casadi_int>(pos));
    return as_shape(r, s);
  }

  // Numpy broadcasting of two shapes
  static bool broadcast_shapes(const std::vector<casadi_int>& a, const std::vector<casadi_int>& b,
                               std::vector<casadi_int>& r) {
    r.resize(std::max(a.size(), b.size()));
    for (size_t k = 0; k < r.size(); ++k) {
      casadi_int da = k < r.size() - a.size() ? 1 : a[k - (r.size() - a.size())];
      casadi_int db = k < r.size() - b.size() ? 1 : b[k - (r.size() - b.size())];
      if (da != db && da != 1 && db != 1) return false;
      r[k] = da == 1 ? db : da;
    }
    return true;
  }

  static std::vector<casadi_int> broadcast_shapes(const std::vector<casadi_int>& a,
                                                  const std::vector<casadi_int>& b) {
    std::vector<casadi_int> r;
    casadi_assert(broadcast_shapes(a, b, r), "ONNX import: cannot broadcast " + shape_str(a)
      + " with " + shape_str(b));
    return r;
  }

  // Row-major positions of the elements selected along each axis
  static std::vector<casadi_int> selected_positions(
      const std::vector<std::vector<casadi_int>>& sel, const std::vector<casadi_int>& strides) {
    std::vector<casadi_int> pos = {0};
    for (size_t k = 0; k < sel.size(); ++k) {
      std::vector<casadi_int> next;
      for (casadi_int p : pos) {
        for (casadi_int i : sel[k]) next.push_back(p + i * strides[k]);
      }
      pos = next;
    }
    return pos;
  }

  // x of shape s broadcast to shape t
  static MX broadcast_to(const MX& x, const std::vector<casadi_int>& s,
                         const std::vector<casadi_int>& t) {
    if (s == t) return as_shape(x, t);
    std::vector<casadi_int> ts = shape_strides(t), ss = shape_strides(s);
    size_t off = t.size() - s.size();
    std::vector<casadi_int> pos(shape_numel(t));
    for (casadi_int j = 0; j < static_cast<casadi_int>(pos.size()); ++j) {
      casadi_int p = 0;
      for (size_t k = off; k < t.size(); ++k) {
        casadi_int i = (j / ts[k]) % t[k];
        if (s[k - off] != 1) p += i * ss[k - off];
      }
      pos[j] = p;
    }
    return take_elements(x, pos, t);
  }

  // Names a graph uses from outer scopes, including those its subgraphs use
  static void free_names(const onnx::GraphProto& g, std::set<std::string>& names) {
    std::set<std::string> defined;
    for (const auto& i : g.input()) defined.insert(i.name());
    for (const auto& t : g.initializer()) defined.insert(t.name());
    for (const auto& t : g.sparse_initializer()) defined.insert(t.values().name());
    for (const auto& n : g.node()) {
      for (const auto& in : n.input()) {
        if (!in.empty() && !defined.count(in)) names.insert(in);
      }
      for (const auto& a : n.attribute()) {
        std::set<std::string> inner;
        if (a.has_g()) free_names(a.g(), inner);
        for (const auto& sub : a.graphs()) free_names(sub, inner);
        for (const auto& x : inner) {
          if (!defined.count(x)) names.insert(x);
        }
      }
      for (const auto& o : n.output()) defined.insert(o);
    }
  }

  static bool is_activation(const std::string& op_type);

  // Elementwise ops: the result has the broadcast shape of the operands
  static bool is_elementwise(const std::string& op_type) {
    if (is_activation(op_type)) return true;
    const OpMapping* m = get_op_mapping_by_name(op_type);
    if (m) {
      return m->arity == 2 || (m->casadi_op != OP_TRANSPOSE && m->casadi_op != OP_NORM1
        && m->casadi_op != OP_NORM2 && m->casadi_op != OP_NORMF && m->casadi_op != OP_MMIN
        && m->casadi_op != OP_MMAX);
    }
    static const std::set<std::string> ops = {"Identity", "Cast", "Not", "And", "Or", "Less",
      "LessOrEqual", "Greater", "GreaterOrEqual", "Equal", "Min", "Max", "Mod", "Where"};
    return ops.count(op_type) > 0;
  }

  // ONNX initializers are pre-loaded constants
  void Onnx::process_graph_initializers(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      bool verbose) const {

    for (int i = 0; i < graph.initializer_size(); ++i) {
      const onnx::TensorProto& tensor = graph.initializer(i);
      std::string tensor_name = tensor.name();

      if (verbose) {
        uout() << "  Processing initializer: " << tensor_name << std::endl;
      }

      value_map[tensor_name] = MX(tensor_to_dm(tensor));
    }
  }

  // The constant for a graph input fixed with GraphBuilder::set. The value is in ONNX row-major
  // order, which is the column-major order of the transposed stored form sp. All values are
  // doubles here, so a Boolean or integer input must be given a value that its type can hold
  // exactly: the numeric backend would convert it, and the two paths would then disagree.
  static MX baked_input(const std::string& name, int32_t elem_type, const Sparsity& sp,
                        const std::vector<double>& value) {
    casadi_assert(static_cast<casadi_int>(value.size()) == sp.numel(),
                  "Baked value for '" + name + "' has " + str(value.size())
                  + " elements, expected " + str(sp.numel()));
    bool is_bool = elem_type == onnx::TensorProto::BOOL;
    bool is_unsigned = elem_type == onnx::TensorProto::UINT8
      || elem_type == onnx::TensorProto::UINT16 || elem_type == onnx::TensorProto::UINT32
      || elem_type == onnx::TensorProto::UINT64;
    bool is_integer = is_unsigned || elem_type == onnx::TensorProto::INT8
      || elem_type == onnx::TensorProto::INT16 || elem_type == onnx::TensorProto::INT32
      || elem_type == onnx::TensorProto::INT64;
    for (double v : value) {
      if (is_bool) {
        casadi_assert(v == 0 || v == 1, "Baked value for Boolean input '" + name
                      + "' must be 0 or 1, got " + str(v));
      } else if (is_integer) {
        casadi_assert(std::isfinite(v) && v == std::floor(v) && (!is_unsigned || v >= 0),
                      "Baked value for integer input '" + name + "' must be "
                      + (is_unsigned ? "a nonnegative" : "an") + " integer, got " + str(v));
      }
    }
    DM d = reshape(DM(value), sp.size1(), sp.size2());
    if (!sp.is_dense()) {
      // A sparse input pattern: entries outside it are structural zeros of the model
      casadi_assert(norm_inf(d - project(d, sp)).scalar() == 0,
                    "Baked value for '" + name + "' is nonzero outside the input's sparsity");
      d = project(d, sp);
    }
    return MX(d);
  }

  // Create MX symbols for graph inputs
  void Onnx::process_graph_inputs(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      std::vector<MX>& func_inputs,
      std::vector<std::string>& input_names,
      const BakedValues& baked,
      bool verbose) const {

    for (int i = 0; i < graph.input_size(); ++i) {
      const onnx::ValueInfoProto& input = graph.input(i);
      std::string input_name = input.name();

      // Skip if already in value_map (it's an initializer, not a variable)
      if (value_map.count(input_name)) {
        if (verbose) {
          uout() << "  Skipping input '" << input_name
                 << "' (it's an initializer)" << std::endl;
        }
        continue;
      }

      // Extract shape
      const onnx::TensorShapeProto& shape =
          input.type().tensor_type().shape();
      casadi_int rows = get_dimension(shape, 0);
      casadi_int cols = get_dimension(shape, 1);
      // Rank > 2: the trailing dimensions are flattened (see shapes_)
      for (int k = 2; k < shape.dim_size(); ++k) cols *= get_dimension(shape, k);

      if (verbose) {
        uout() << "  Creating input: " << input_name
               << " [" << rows << ", " << cols << "]" << std::endl;
      }

      // Sparsity overlay: if a sparse_initializer shares this input's name it carries a non-dense
      // input pattern -- make a SPARSE input symbol and feed it AS-IS, so CasADi propagates the
      // true sparsity through the body (no densification on entry; densify only where an ONNX op
      // requires
      // the dense layout, e.g. Reshape). Transpose-rep: a rank-2 ONNX input (r,c) is stored as
      // (c,r); a rank<2 input [n] is a column (n,1), matching tensor_to_dm for initializers so
      // that e.g. a 1-D Gemm bias broadcasts the same whether it is an initializer or an input.
      Sparsity ov = input_pattern(graph, input_name);
      Sparsity sp = !ov.is_null() ? ov
                  : shape.dim_size() < 2 ? Sparsity::dense(rows, cols)
                  : Sparsity::dense(cols, rows);

      // An input fixed with GraphBuilder::set is a constant, not a Function input
      auto bv = baked.find(input_name);
      if (bv != baked.end()) {
        value_map[input_name] = baked_input(input_name, input.type().tensor_type().elem_type(),
                                            sp, bv->second);
        if (verbose) uout() << "  Input '" << input_name << "' is baked" << std::endl;
        continue;
      }

      MX mx_input = MX::sym(input_name, sp);
      value_map[input_name] = mx_input;
      func_inputs.push_back(mx_input);
      input_names.push_back(input_name);
    }
  }

  // ONNX ranks the graph declares (inputs, intermediate value_info, outputs, initializers)
  static std::map<std::string, casadi_int> declared_ranks(const onnx::GraphProto& graph) {
    std::map<std::string, casadi_int> ranks;
    for (const auto* infos : {&graph.input(), &graph.value_info(), &graph.output()}) {
      for (const auto& info : *infos) {
        if (info.type().tensor_type().has_shape()) {
          ranks[info.name()] = info.type().tensor_type().shape().dim_size();
        }
      }
    }
    for (const auto& tensor : graph.initializer()) ranks[tensor.name()] = tensor.dims_size();
    return ranks;
  }

  void Onnx::declare_shapes(const onnx::GraphProto& graph) {
    for (const auto* infos : {&graph.input(), &graph.value_info(), &graph.output()}) {
      for (const auto& info : *infos) {
        if (!info.type().tensor_type().has_shape()) continue;
        const onnx::TensorShapeProto& shape = info.type().tensor_type().shape();
        std::vector<casadi_int> s;
        for (const auto& d : shape.dim()) {
          auto it = d.has_dim_param() ? dimension_overrides_.find(d.dim_param())
                                      : dimension_overrides_.end();
          s.push_back(d.has_dim_value() ? static_cast<casadi_int>(d.dim_value())
                      : it != dimension_overrides_.end() ? it->second : -1);
        }
        shapes_[info.name()] = s;
      }
    }
    for (const auto& tensor : graph.initializer()) {
      shapes_[tensor.name()] = std::vector<casadi_int>(tensor.dims().begin(), tensor.dims().end());
    }
  }

  std::vector<casadi_int> Onnx::onnx_shape(const std::string& name, const MX& v) const {
    // A recorded shape that does not fit the value is ignored: exporters declare wrong shapes
    // for some outputs, and subgraphs may reuse an outer name
    auto it = shapes_.find(name);
    if (it != shapes_.end()) {
      const std::vector<casadi_int>& s = it->second;
      if (s.size() > 2) {
        // The trailing dimensions of a declared shape may be unknown; d0 is the stored columns
        std::vector<casadi_int> r = s;
        r[0] = v.size2();
        casadi_int known = 1, unknown = -1, n_unknown = 0;
        for (size_t k = 1; k < r.size(); ++k) {
          if (r[k] >= 0) {
            known *= r[k];
          } else {
            unknown = static_cast<casadi_int>(k);
            n_unknown++;
          }
        }
        if (n_unknown == 1 && known > 0) r[unknown] = v.size1() / known;
        if (n_unknown <= 1 && shape_numel(r) == v.numel()
            && stored_size(r).first == v.size1()) return r;
      } else if (s.size() == 2) {
        return {v.size2(), v.size1()};
      } else if (s.size() == 1 && v.is_vector()) {
        return {v.numel()};
      } else if (s.empty() && v.is_scalar()) {
        return {};
      }
    }
    if (v.is_scalar()) return {};
    if (v.is_vector()) return {v.numel()};
    return {v.size2(), v.size1()};
  }

  // ONNX counts a negative "axis" attribute from the end, relative to the rank of input 0
  // (OneHot, relative to the output rank, aside). CasADi values are 2-D and do not keep the ONNX
  // rank, e.g. [1, n] and [n] are both stored as (n, 1), so resolve the axis against the rank
  // the graph declares. Returns node itself if there is nothing to resolve, else *copy.
  static const onnx::NodeProto& resolve_negative_axis(const onnx::NodeProto& node,
      const std::map<std::string, casadi_int>& ranks, onnx::NodeProto* copy) {
    if (node.op_type() == "OneHot" || node.input_size() == 0) return node;
    for (int a = 0; a < node.attribute_size(); ++a) {
      if (node.attribute(a).name() != "axis" || node.attribute(a).i() >= 0) continue;
      auto it = ranks.find(node.input(0));
      if (it == ranks.end()) return node;  // undeclared rank: leave it to the op
      *copy = node;
      copy->mutable_attribute(a)->set_i(node.attribute(a).i() + it->second);
      return *copy;
    }
    return node;
  }

  // Process all nodes in the graph
  void Onnx::process_graph_nodes(
      const onnx::GraphProto& graph,
      std::map<std::string, MX>& value_map,
      bool verbose) {

    IntegerConstants integer_constants;
    for (const auto& tensor : graph.initializer()) {
      std::vector<int64_t> values;
      if (integer_tensor(tensor, values)) integer_constants[tensor.name()] = std::move(values);
    }

    std::map<std::string, casadi_int> ranks = declared_ranks(graph);
    declare_shapes(graph);

    // Per-graph capture of kron operands, keyed by the exporter's node-name tag "kron<k>".
    std::map<std::string, MX> kron_operand_a_, kron_operand_b_;

    for (int i = 0; i < graph.node_size(); ++i) {
      onnx::NodeProto resolved;
      const onnx::NodeProto& node = resolve_negative_axis(graph.node(i), ranks, &resolved);
      std::string op_type = node.op_type();

      if (verbose) {
        uout() << "  Processing node " << i << ": " << op_type << std::endl;
      }

      // Only direct Constant nodes and initializers retain an integer payload.
      if (op_type == "Constant" && node.output_size() == 1) {
        for (const auto& attr : node.attribute()) {
          if (attr.name() == "value") {
            shapes_[node.output(0)] = std::vector<casadi_int>(attr.t().dims().begin(),
                                                              attr.t().dims().end());
          }
          std::vector<int64_t> values;
          if (attr.name() == "value_int") {
            values = {attr.i()};
          } else if (attr.name() == "value_ints") {
            values.assign(attr.ints().begin(), attr.ints().end());
          } else if (attr.name() != "value" || !integer_tensor(attr.t(), values)) {
            continue;
          }
          integer_constants[node.output(0)] = std::move(values);
        }
      }

      // ========== Special Handling: Kronecker product group ==========
      // The exporter tags its kron node group with node name "kron<k>" (see OP_KRON in
      // onnx_operations.cpp). The 4-D intermediate tensors cannot be MX (MX is 2-D), so we do NOT
      // process the group generically (the helper tensors kron<k>_A4/_B4/_P are never added to the
      // value_map). We capture the two logical operands at the leading Reshapes (whose single data
      // input is the imported, transposed-back A / B) and reconstruct kron(A_imp, B_imp) at the
      // final Reshape. CasADi kron propagates sparsity natively, so the pattern round-trips for
      // free
      // (kron commutes with transpose). This is handled BEFORE the generic input gathering because
      // the helper tensors are intentionally absent from value_map.
      if (node.name().rfind("kron", 0) == 0) {
        const std::string& nm = node.name();  // "kron<k>_<role>"
        std::size_t us = nm.rfind('_');
        std::string tag = nm.substr(0, us), role = nm.substr(us + 1);
        if (role == "A4") {
          kron_operand_a_[tag] = value_map.at(node.input(0));  // data input = imported A
          continue;
        }
        if (role == "B4") {
          kron_operand_b_[tag] = value_map.at(node.input(0));  // data input = imported B
          continue;
        }
        if (role == "P") continue;  // intermediate broadcast product; nothing to record
        if (role == "R") {
          // Final Reshape: emit kron of the captured operands and skip generic 4-D handling.
          casadi_assert(kron_operand_a_.count(tag) && kron_operand_b_.count(tag),
                        "ONNX import: incomplete kron node group '" + tag + "'");
          value_map[node.output(0)] = MX::kron(kron_operand_a_.at(tag), kron_operand_b_.at(tag));
          continue;
        }
      }

      // Gather input tensors
      std::vector<MX> node_inputs;
      for (int j = 0; j < node.input_size(); ++j) {
        std::string input_name = node.input(j);

        // Empty string means optional input not provided
        if (input_name.empty()) {
          node_inputs.push_back(MX());
          continue;
        }

        casadi_assert(value_map.count(input_name),
                      "Unknown input tensor '" + input_name +
                      "' required by node " + std::to_string(i) +
                      " (op_type: " + op_type + ")");
        node_inputs.push_back(value_map[input_name]);
      }

      // Compute output based on operation type
      MX output;
      std::vector<casadi_int> out_shape;  // its ONNX shape, if shaped
      bool shaped = false;

      // ========== Special Handling: Multi-output operations ==========
      if (op_type == "Split") {
        casadi_assert(node_inputs.size() >= 1, "Split requires 1 input");
        casadi_int axis = get_int_attribute(node, "axis", 0);
        casadi_assert(axis == 0 || axis == 1, "Split: only axis 0 and 1 supported");

        // Split sizes: 2nd input (opset 13+), else 'split' attribute (opset <13), else equal
        std::vector<casadi_int> split_sizes;
        if (node_inputs.size() >= 2) {
          DM split_dm = DM(node_inputs[1]);
          for (casadi_int k = 0; k < split_dm.nnz(); ++k) {
            split_sizes.push_back(static_cast<casadi_int>(static_cast<double>(split_dm.nz(k))));
          }
        } else {
          for (int a = 0; a < node.attribute_size(); ++a) {
            if (node.attribute(a).name() == "split") {
              for (int k = 0; k < node.attribute(a).ints_size(); ++k) {
                split_sizes.push_back(node.attribute(a).ints(k));
              }
              break;
            }
          }
        }
        if (split_sizes.empty()) {  // equal split across the outputs
          // Transpose-rep: ONNX axis 0 is CasADi columns, axis 1 is CasADi rows
          casadi_int total = (axis == 0) ? node_inputs[0].size2() : node_inputs[0].size1();
          split_sizes.assign(node.output_size(), total / node.output_size());
        }

        // Offsets [0, s0, s0+s1, ...]; ONNX axis 0 -> CasADi horzsplit, axis 1 -> vertsplit
        std::vector<casadi_int> offset = {0};
        for (casadi_int sz : split_sizes) offset.push_back(offset.back() + sz);
        std::vector<MX> outputs = (axis == 0) ? horzsplit(node_inputs[0], offset)
                                              : vertsplit(node_inputs[0], offset);

        for (casadi_int j = 0; j < outputs.size(); ++j) value_map[node.output(j)] = outputs[j];
        continue;  // Don't use standard output handling

      // ========== Scan: iterate a body over columns -> CasADi Map ==========
      } else if (op_type == "Scan") {
        const onnx::GraphProto* body = get_graph_attribute(node, "body");
        casadi_assert(body != nullptr, "Scan node requires a 'body' subgraph");
        casadi_int num_scan_inputs = get_int_attribute(node, "num_scan_inputs", node.input_size());
        casadi_int M = node.input_size() - num_scan_inputs;  // state variables (reduce_out count)

        if (M == 0) {
          // No state variables: the body is evaluated for each slice along ONNX axis 0 of the
          // scan inputs, a CasADi map. The outer-scope values the body uses (captures) become
          // extra inputs of the body Function that every evaluation gets unchanged (reduce_in);
          // constant ones are baked into it, so that they can still serve as indices.
          for (const char* a : {"scan_input_axes", "scan_input_directions", "scan_output_axes",
                                "scan_output_directions"}) {
            for (const auto& attr : node.attribute()) {
              if (attr.name() != a) continue;
              for (auto v : attr.ints()) {
                casadi_assert(v == 0, std::string("ONNX import: Scan with nonzero ") + a
                  + " is not supported");
              }
            }
          }
          onnx::GraphProto g = *body;
          // A scan input slice without a declared shape: the trailing dimensions of the input
          for (casadi_int j = 0; j < num_scan_inputs; ++j) {
            onnx::ValueInfoProto* vi = g.mutable_input(j);
            if (vi->type().tensor_type().has_shape()) continue;
            std::vector<casadi_int> sj = onnx_shape(node.input(M + j), node_inputs[M + j]);
            vi->mutable_type()->mutable_tensor_type()->set_elem_type(real_type());
            auto* shape = vi->mutable_type()->mutable_tensor_type()->mutable_shape();
            for (size_t k = 1; k < sj.size(); ++k) shape->add_dim()->set_dim_value(sj[k]);
          }
          std::set<std::string> captures;
          free_names(*body, captures);
          BakedValues baked;
          std::vector<MX> args;
          for (casadi_int j = 0; j < num_scan_inputs; ++j) args.push_back(node_inputs[M + j]);
          std::vector<bool> reduce_in(num_scan_inputs, false);
          for (const std::string& c : captures) {
            casadi_assert(value_map.count(c), "Unknown tensor '" + c + "' used by a Scan body");
            const MX& v = value_map.at(c);
            onnx::ValueInfoProto* vi = g.add_input();
            vi->set_name(c);
            vi->mutable_type()->mutable_tensor_type()->set_elem_type(real_type());
            auto* shape = vi->mutable_type()->mutable_tensor_type()->mutable_shape();
            for (casadi_int d : onnx_shape(c, v)) shape->add_dim()->set_dim_value(d);
            if (v.is_constant()) {
              DM dv = densify(static_cast<DM>(v));
              baked[c] = std::vector<double>(dv->begin(), dv->end());
            } else {
              args.push_back(v);
              reduce_in.push_back(true);
            }
          }
          Function base = function_from_graph(g, op_type + "_body", baked);
          // The slices of a scan input side by side: column-major, a stored [n, ...] tensor is
          // the slices' row-major elements one after the other (see shapes_)
          casadi_int n = -1;
          for (casadi_int j = 0; j < num_scan_inputs; ++j) {
            casadi_int r = base.size1_in(j), c = base.size2_in(j);
            casadi_int nj = r * c == 0 ? 0 : args[j].numel() / (r * c);
            casadi_assert(nj * r * c == args[j].numel() && (n < 0 || nj == n),
              "ONNX import: Scan input " + str(j) + " (" + args[j].dim()
              + ") does not match the body input (" + base.sparsity_in(j).dim() + ")");
            n = nj;
            args[j] = reshape(args[j], r, c * n);
          }
          bool any_reduce = std::find(reduce_in.begin(), reduce_in.end(), true) != reduce_in.end();
          std::vector<MX> outs = (any_reduce ? base.map(n, reduce_in) : base.map(n))(args);
          for (int j = 0; j < node.output_size(); ++j) {
            std::vector<casadi_int> sj = onnx_shape(body->output(j).name(),
                                                    MX::zeros(base.sparsity_out(j)));
            casadi_int r = base.size1_out(j) * base.size2_out(j);
            value_map[node.output(j)] = r == 1 ? outs[j] : reshape(outs[j], r, n);
            sj.insert(sj.begin(), n);
            shapes_[node.output(j)] = sj;
          }
          continue;
        }

        std::set<std::string> body_in_names;
        for (int b = 0; b < body->input_size(); ++b) body_in_names.insert(body->input(b).name());

        // Reduce-map: the body wraps a single base call; recover base + masks and rebuild
        // base.map(n, reduce_in, reduce_out).
        const onnx::NodeProto* call = nullptr;
        for (int nn = 0; nn < body->node_size(); ++nn)
          if (!body->node(nn).domain().empty()) { call = &body->node(nn); break; }
        casadi_assert(call, "reduce-map Scan body must contain a base function call");
        casadi_int nin = call->input_size(), nout = call->output_size();

        // reduce_in: base arg is a capture; reduce_out: base result feeds an accumulator Add
        std::vector<bool> reduce_in(nin), reduce_out(nout, false);
        for (casadi_int j = 0; j < nin; ++j) reduce_in[j] = !body_in_names.count(call->input(j));
        for (int nn = 0; nn < body->node_size(); ++nn) {
          if (body->node(nn).op_type() != "Add") continue;
          for (int ii = 0; ii < body->node(nn).input_size(); ++ii)
            for (casadi_int j = 0; j < nout; ++j)
              if (body->node(nn).input(ii) == call->output(j)) reduce_out[j] = true;
        }

        // Base input shapes: captures from the outer tensor, scanned ones from the body input
        std::vector<std::pair<casadi_int, casadi_int>> in_shapes(nin);
        for (casadi_int j = 0; j < nin; ++j) {
          if (reduce_in[j]) {
            MX cap = value_map.at(call->input(j));
            in_shapes[j] = {cap.size1(), cap.size2()};
          } else {
            for (int b = 0; b < body->input_size(); ++b)
              if (body->input(b).name() == call->input(j)) {
                const onnx::TensorShapeProto& sh = body->input(b).type().tensor_type().shape();
                in_shapes[j] = {get_dimension(sh, 0), get_dimension(sh, 1)};
                break;
              }
          }
        }

        const onnx::FunctionProto* fp = find_function(call->op_type(), call->domain());
        casadi_assert(fp, "reduce-map base function '" + call->op_type() + "' not found");
        Function base = function_from_function_proto(*fp, in_shapes, op_type + "_base");

        // base.map args (base order): captures broadcast, repeated come from scan node inputs
        std::vector<MX> args(nin);
        casadi_int n = 0, r = 0;
        for (casadi_int j = 0; j < nin; ++j) {
          if (reduce_in[j]) {
            args[j] = value_map.at(call->input(j));
          } else {
            // The slices side by side, as in the map above
            casadi_int rj = base.size1_in(j), cj = base.size2_in(j);
            if (n == 0 && rj * cj > 0) n = node_inputs[M + r].numel() / (rj * cj);
            args[j] = reshape(node_inputs[M + r], rj, node_inputs[M + r].numel() / rj);
            ++r;
          }
        }
        casadi_assert(n > 0, "reduce-map import: could not infer map size");
        std::vector<MX> outs = base.map(n, reduce_in, reduce_out)(args);

        // Node outputs are [state accumulators (reduce_out), scan outputs (repeated)], base order
        casadi_int si = 0, ci = 0;
        for (casadi_int j = 0; j < nout; ++j) {
          if (reduce_out[j]) value_map[node.output(si++)] = outs[j];
          else
            value_map[node.output(M + ci++)] = outs[j];
        }
        continue;

      // ========== If: two captured branches combined with if_else ==========
      } else if (op_type == "If") {
        const onnx::GraphProto* then_b = get_graph_attribute(node, "then_branch");
        const onnx::GraphProto* else_b = get_graph_attribute(node, "else_branch");
        casadi_assert(then_b && else_b, "If node requires then_branch and else_branch");
        MX cond = node_inputs[0];
        // A condition known at import (from initializers, Constant nodes or inputs fixed with
        // GraphBuilder::set) selects a branch: the other one is not built
        if (cond.is_constant()) {
          casadi_assert(cond.numel() == 1, "If: the condition must have one element, got "
                        + cond.dim());
          bool taken = densify(static_cast<DM>(cond)).scalar() != 0;
          const onnx::GraphProto& branch = taken ? *then_b : *else_b;
          std::vector<MX> r = eval_captured_subgraph(branch, value_map);
          for (int j = 0; j < node.output_size(); ++j) {
            value_map[node.output(j)] = r[j];
            copy_shape(branch.output(j).name(), node.output(j));
          }
          continue;
        }
        std::vector<MX> t = eval_captured_subgraph(*then_b, value_map);
        std::vector<MX> e = eval_captured_subgraph(*else_b, value_map);
        for (int j = 0; j < node.output_size(); ++j) {
          value_map[node.output(j)] = if_else(cond, t[j], e[j]);
          copy_shape(then_b->output(j).name(), node.output(j));
        }
        continue;

      // ========== Loop: a counted loop, unrolled ==========
      } else if (op_type == "Loop") {
        const onnx::GraphProto* body = get_graph_attribute(node, "body");
        casadi_assert(body != nullptr, "Loop node requires a 'body' subgraph");
        casadi_assert(node_inputs.size() >= 2 && node_inputs[0].is_constant()
          && !node_inputs[0].is_empty(), "ONNX import: only a Loop with a constant trip count "
          "is supported.");
        casadi_assert(node_inputs[1].is_empty() || (node_inputs[1].is_constant()
          && static_cast<DM>(node_inputs[1]).scalar() != 0),
          "ONNX import: only a Loop with a constant true initial condition is supported.");
        casadi_int trip_count = static_cast<casadi_int>(static_cast<DM>(node_inputs[0]).scalar());

        // Body: (iteration, condition, carried...) -> (condition, carried..., scan outputs...)
        std::vector<MX> carried(node_inputs.begin() + 2, node_inputs.end());
        casadi_int n_scan = node.output_size() - static_cast<casadi_int>(carried.size());
        casadi_assert(body->input_size() == 2 + static_cast<casadi_int>(carried.size())
          && n_scan >= 0 && body->output_size() == 1 + node.output_size(),
          "Loop body inputs and outputs do not match the node");
        std::vector<std::vector<MX>> scan(n_scan);
        std::vector<std::vector<casadi_int>> scan_shape(n_scan);  // per iteration
        for (casadi_int it = 0; it < trip_count; ++it) {
          std::map<std::string, MX> scope = value_map;
          scope[body->input(0).name()] = MX(static_cast<double>(it));
          scope[body->input(1).name()] = MX(1);
          for (casadi_int j = 0; j < static_cast<casadi_int>(carried.size()); ++j) {
            scope[body->input(2 + j).name()] = carried[j];
          }
          std::vector<MX> out = eval_captured_subgraph(*body, scope);
          casadi_assert(out[0].is_constant() && static_cast<DM>(out[0]).scalar() != 0,
            "ONNX import: a Loop that can exit early is not supported.");
          for (casadi_int j = 0; j < static_cast<casadi_int>(carried.size()); ++j) {
            carried[j] = out[1 + j];
          }
          for (casadi_int j = 0; j < n_scan; ++j) {
            const MX& e = out[1 + carried.size() + j];
            std::vector<casadi_int> s = onnx_shape(body->output(1 + carried.size() + j).name(), e);
            casadi_assert(it == 0 || s == scan_shape[j], "ONNX import: Loop scan output "
              + str(j) + " changes shape from " + shape_str(scan_shape[j]) + " to " + shape_str(s));
            scan_shape[j] = s;
            scan[j].push_back(vec(e));
          }
        }
        for (casadi_int j = 0; j < static_cast<casadi_int>(carried.size()); ++j) {
          value_map[node.output(j)] = carried[j];
          if (trip_count > 0) copy_shape(body->output(1 + j).name(), node.output(j));
        }
        // Scan outputs stack along a new ONNX axis 0, which is CasADi columns as in Concat: the
        // stacked [n, ...] tensor is the n columns of the elements' row-major vec (see shapes_)
        for (casadi_int j = 0; j < n_scan; ++j) {
          value_map[node.output(carried.size() + j)] = horzcat(scan[j]);
          if (trip_count > 0) {
            scan_shape[j].insert(scan_shape[j].begin(), trip_count);
            shapes_[node.output(carried.size() + j)] = scan_shape[j];
          }
        }
        continue;

      // ========== Standard Operations (delegated to helper) ==========
      } else {
        // A non-empty domain marks a function call.
        std::string node_domain = node.domain();
        if (!node_domain.empty()) {
          const onnx::FunctionProto* func_proto = find_function(op_type, node_domain);

          if (func_proto != nullptr) {
            if (verbose) {
              uout() << "    Function call to: " << node_domain << "." << op_type << std::endl;
            }

            // Inline the function body: only its nodes are needed (process_graph_nodes reads
            // graph.node() only); formal inputs/outputs are wired through func_value_map below.
            onnx::GraphProto func_graph;
            func_graph.set_name(func_proto->name());
            for (int n = 0; n < func_proto->node_size(); ++n) {
              *func_graph.add_node() = func_proto->node(n);
            }

            // Map the function's formal inputs to the call-site values. The body has a name
            // scope of its own, so shapes recorded in it are dropped afterwards.
            std::map<std::string, MX> func_value_map;
            std::map<std::string, std::vector<casadi_int>> outer_shapes = shapes_;
            for (size_t n = 0; n < node_inputs.size() && n < func_proto->input_size(); ++n) {
              func_value_map[func_proto->input(n)] = node_inputs[n];
              shapes_[func_proto->input(n)] = onnx_shape(node.input(n), node_inputs[n]);
            }

            process_graph_nodes(func_graph, func_value_map, false);

            // Store the function outputs in the outer value_map.
            std::vector<std::vector<casadi_int>> out_shapes;
            for (int n = 0; n < node.output_size() && n < func_proto->output_size(); ++n) {
              const std::string& o = func_proto->output(n);
              out_shapes.push_back(func_value_map.count(o) ? onnx_shape(o, func_value_map[o])
                                                           : std::vector<casadi_int>());
            }
            shapes_ = outer_shapes;
            for (int n = 0; n < node.output_size() && n < func_proto->output_size(); ++n) {
              std::string func_output_name = func_proto->output(n);
              casadi_assert(func_value_map.count(func_output_name),
                           "Function output '" + func_output_name + "' not found");
              value_map[node.output(n)] = func_value_map[func_output_name];
              shapes_[node.output(n)] = out_shapes[n];

              if (verbose) {
                uout() << "      -> " << node.output(n) << std::endl;
              }
            }

            // Skip the normal single-output handling
            continue;
          }
          // Fall through to error if function not found
        }

        // All other (single-output) operations; those that need ONNX shapes first
        if (process_shaped_operation(op_type, node, node_inputs, integer_constants,
                                     output, out_shape)) {
          shaped = true;
        } else {
          output = process_node_operation(op_type, node, node_inputs, integer_constants);
          // An elementwise op keeps the broadcast shape of its operands. If an operand has a
          // recorded shape, this overrides a declared one, which exporters get wrong at times.
          if (is_elementwise(op_type)) {
            bool recorded = false, ok = true;
            for (size_t j = 0; j < node_inputs.size() && ok; ++j) {
              if (node.input(j).empty() || (op_type == "Clip" && j > 0)) continue;
              recorded = recorded || shapes_.count(node.input(j));
              std::vector<casadi_int> sj = onnx_shape(node.input(j), node_inputs[j]);
              ok = j == 0 ? (out_shape = sj, true) : broadcast_shapes(out_shape, sj, out_shape);
            }
            auto sz = stored_size(out_shape);
            shaped = recorded && ok && sz.first * sz.second == output.numel()
              && (sz.first == output.size1() || out_shape.size() < 2);
          }
        }
      }

      // Store output (assume single output for now)
      casadi_assert(node.output_size() >= 1,
                    "Node must have at least one output");

      std::string output_name = node.output(0);
      value_map[output_name] = output;
      if (shaped) shapes_[output_name] = out_shape;

      if (verbose) {
        uout() << "    -> " << output_name << " " << output.dim() << std::endl;
      }
    }
  }

  // Collect graph outputs from value_map
  void Onnx::collect_graph_outputs(
      const onnx::GraphProto& graph,
      const std::map<std::string, MX>& value_map,
      std::vector<MX>& func_outputs,
      std::vector<std::string>& output_names,
      bool verbose) const {

    for (int i = 0; i < graph.output_size(); ++i) {
      const onnx::ValueInfoProto& output = graph.output(i);
      std::string output_name = output.name();

      casadi_assert(value_map.count(output_name),
                    "Unknown output tensor: " + output_name +
                    ". This usually means the ONNX graph contains unsupported operations.");

      // No output overlay: each op restores its own output sparsity from seeds during export, so
      // the value already carries the exact CasADi pattern (re-propagated natively from the seeds).
      func_outputs.push_back(value_map.at(output_name));
      output_names.push_back(output_name);

      if (verbose) {
        uout() << "  Graph output: " << output_name << std::endl;
      }
    }
  }

  const onnx::FunctionProto* Onnx::find_function(const std::string& name,
                                                 const std::string& domain) const {
    for (int f = 0; f < model_.functions_size(); ++f)
      if (model_.functions(f).name() == name && model_.functions(f).domain() == domain)
        return &model_.functions(f);
    return nullptr;
  }

  // Whether a node has a negative axis that resolve_negative_axis cannot resolve
  static bool has_unresolved_negative_axis(const onnx::GraphProto& graph) {
    std::map<std::string, casadi_int> ranks = declared_ranks(graph);
    for (const auto& node : graph.node()) {
      if (node.op_type() == "OneHot" || node.input_size() == 0) continue;
      if (get_int_attribute(node, "axis", 0) < 0 && !ranks.count(node.input(0))) return true;
    }
    return false;
  }

  Function Onnx::create(const std::string& name, const BakedValues& baked) {
    casadi_assert(has_model_, "No ONNX model loaded. Call load() first.");
    if (!has_unresolved_negative_axis(model_.graph())) {
      return function_from_graph(model_.graph(), name, baked);
    }
    // Exporters need not declare the shapes of intermediate tensors (PyTorch's TorchScript
    // exporter doesn't). Infer them on a copy, leaving model_ as loaded, so that negative axes
    // can be resolved. If inference fails, import as is: the affected op reports the axis.
    onnx::ModelProto inferred = model_;
    try {
      onnx::shape_inference::InferShapes(inferred, onnx::OpSchemaRegistry::Instance(),
        onnx::ShapeInferenceOptions(false, 0, true));
    } catch (const std::exception& e) {
      if (verbose_) uout() << "ONNX shape inference failed: " << e.what() << std::endl;
      return function_from_graph(model_.graph(), name, baked);
    }
    return function_from_graph(inferred.graph(), name, baked);
  }

  Function Onnx::function_from_graph(const onnx::GraphProto& graph,
                                              const std::string& name,
                                              const BakedValues& baked) {
    std::map<std::string, MX> value_map;
    std::vector<MX> inputs, outputs;
    std::vector<std::string> input_names, output_names;

    if (verbose_) {
      uout() << "Building CasADi Function '" << name << "' from graph '" << graph.name()
             << "' (" << graph.initializer_size() << " initializers, "
             << graph.input_size() << " inputs, " << graph.output_size() << " outputs, "
             << graph.node_size() << " nodes)" << std::endl;
    }

    process_graph_initializers(graph, value_map, verbose_);
    process_graph_inputs(graph, value_map, inputs, input_names, baked, verbose_);
    process_graph_nodes(graph, value_map, verbose_);
    collect_graph_outputs(graph, value_map, outputs, output_names, verbose_);

    return Function(name, inputs, outputs, input_names, output_names);
  }

  Function Onnx::function_from_function_proto(
      const onnx::FunctionProto& fp,
      const std::vector<std::pair<casadi_int, casadi_int>>& in_shapes,
      const std::string& name) {
    // A FunctionProto carries no shapes; wrap it as a graph with caller-supplied input shapes
    onnx::GraphProto g;
    g.set_name(fp.name());
    for (int j = 0; j < fp.input_size(); ++j) {
      onnx::ValueInfoProto* vi = g.add_input();
      vi->set_name(fp.input(j));
      onnx::TypeProto::Tensor* tt = vi->mutable_type()->mutable_tensor_type();
      tt->set_elem_type(real_type());
      tt->mutable_shape()->add_dim()->set_dim_value(in_shapes[j].first);
      tt->mutable_shape()->add_dim()->set_dim_value(in_shapes[j].second);
    }
    for (int k = 0; k < fp.node_size(); ++k) *g.add_node() = fp.node(k);
    for (int j = 0; j < fp.output_size(); ++j) g.add_output()->set_name(fp.output(j));
    return function_from_graph(g, name);
  }

  std::vector<MX> Onnx::eval_captured_subgraph(const onnx::GraphProto& graph,
                                                        std::map<std::string, MX> scope) {
    // No formal inputs: the branch captures outer tensors, already present in `scope`
    process_graph_initializers(graph, scope, verbose_);
    process_graph_nodes(graph, scope, verbose_);
    std::vector<MX> outputs;
    for (int i = 0; i < graph.output_size(); ++i)
      outputs.push_back(scope.at(graph.output(i).name()));
    return outputs;
  }


  // Convert a constant tensor input to a vector of integers
  // Whether a value is known at import: a constant, or an expression of constants that is not
  // folded, such as an index computed in an unrolled Loop
  static bool is_known(const MX& m) {
    return m.is_constant() || symvar(m).empty();
  }

  static std::vector<casadi_int> constant_ints(const MX& m) {
    casadi_assert(is_known(m), "Expected a constant integer tensor");
    DM dm = m.is_constant() ? static_cast<DM>(m) : evalf(m);
    std::vector<casadi_int> v;
    for (casadi_int k = 0; k < dm.numel(); ++k)
      v.push_back(static_cast<casadi_int>(dm(k).scalar()));
    return v;
  }

  // Numerically stable softplus log(1 + exp(x)). Above x = 40, log1p(exp(x)) == x in double;
  // the clamp keeps exp finite in the discarded branch, so reverse mode gets no 0*inf = NaN.
  static MX softplus(const MX& x) {
    return if_else(x > 40, x, log1p(exp(fmin(x, 50))));
  }

  // Equal to 1/(1 + exp(-x)), but with no overflow for large |x|
  static MX sigmoid(const MX& x) {
    return 0.5 * (1 + tanh(0.5 * x));
  }

  // alpha * (exp(x) - 1) for x <= 0, x otherwise. The clamp keeps expm1 finite for large x.
  static MX elu(const MX& x, double alpha) {
    return if_else(x > 0, x, alpha * expm1(fmin(x, 1)));
  }

  static bool is_activation(const std::string& op_type) {
    static const std::set<std::string> ops = {"Relu", "LeakyRelu", "PRelu", "Sigmoid",
      "HardSigmoid", "HardSwish", "Softplus", "Softsign", "Elu", "Selu", "Celu", "Mish", "Gelu",
      "Clip", "ThresholdedRelu", "Shrink", "Swish"};
    return ops.count(op_type) > 0;
  }

  // Element-wise activation functions, composed of CasADi primitives so that AD and code
  // generation work unchanged. Non-smooth ones get CasADi's usual subgradients at the kinks.
  // Attribute defaults are float literals: ONNX float attributes are float32.
  static MX activation(const std::string& op_type, const onnx::NodeProto& node,
                       const std::vector<MX>& node_inputs) {
    casadi_assert(node_inputs.size() >= 1, op_type + " requires 1 input");
    const MX& x = node_inputs[0];

    if (op_type == "Relu") {
      return fmax(x, 0);
    } else if (op_type == "LeakyRelu") {
      return if_else(x >= 0, x, get_float_attribute(node, "alpha", 0.01f) * x);
    } else if (op_type == "PRelu") {
      casadi_assert(node_inputs.size() >= 2, "PRelu requires 2 inputs");
      MX bx = x, slope = node_inputs[1];
      broadcast_pair(&bx, &slope);
      return if_else(bx >= 0, bx, slope * bx);
    } else if (op_type == "ThresholdedRelu") {
      return if_else(x > get_float_attribute(node, "alpha", 1.0), x, 0);
    } else if (op_type == "Shrink") {
      double lambd = get_float_attribute(node, "lambd", 0.5);
      double bias = get_float_attribute(node, "bias", 0.0);
      return if_else(x < -lambd, x + bias, if_else(x > lambd, x - bias, 0));
    } else if (op_type == "Sigmoid") {
      return sigmoid(x);
    } else if (op_type == "Swish") {
      return x * sigmoid(get_float_attribute(node, "alpha", 1.0) * x);
    } else if (op_type == "HardSigmoid") {
      double alpha = get_float_attribute(node, "alpha", 0.2f);
      double beta = get_float_attribute(node, "beta", 0.5f);
      return fmax(0, fmin(1, alpha * x + beta));
    } else if (op_type == "HardSwish") {
      return x * fmax(0, fmin(1, x / 6 + 0.5));
    } else if (op_type == "Softplus") {
      return softplus(x);
    } else if (op_type == "Softsign") {
      return x / (1 + fabs(x));
    } else if (op_type == "Elu") {
      return elu(x, get_float_attribute(node, "alpha", 1.0));
    } else if (op_type == "Selu") {
      double alpha = get_float_attribute(node, "alpha", 1.67326319217681884765625f);
      double gamma = get_float_attribute(node, "gamma", 1.05070102214813232421875f);
      return gamma * elu(x, alpha);
    } else if (op_type == "Celu") {
      // max(0, x) + min(0, alpha * (exp(x / alpha) - 1)) is, for either sign of alpha, x for
      // x > 0 and alpha * expm1(x / alpha) otherwise. Written so that it stays C1 at 0.
      double alpha = get_float_attribute(node, "alpha", 1.0);
      casadi_assert(alpha != 0, "Celu requires a nonzero alpha");
      return if_else(x > 0, x, alpha * expm1(fmin(x, std::fabs(alpha)) / alpha));
    } else if (op_type == "Mish") {
      return x * tanh(softplus(x));
    } else if (op_type == "Gelu") {
      std::string approximate = get_string_attribute(node, "approximate");
      if (approximate.empty() || approximate == "none") {
        return 0.5 * x * (1 + erf(x / std::sqrt(2.0)));
      }
      casadi_assert(approximate == "tanh",
                    "Gelu: unsupported approximate '" + approximate + "'");
      return 0.5 * x * (1 + tanh(std::sqrt(2 / pi) * (x + 0.044715 * pow(x, 3))));
    } else if (op_type == "Clip") {
      // Opset >= 11: optional inputs min and max (an empty name is absent).
      // Opset < 11: attributes min and max. An absent bound is not applied.
      MX output = x;
      for (size_t k = 1; k <= 2; ++k) {
        MX bound;
        if (node_inputs.size() > k) {
          bound = node_inputs[k];
        } else {
          double inf = std::numeric_limits<double>::infinity();
          double b = get_float_attribute(node, k == 1 ? "min" : "max", k == 1 ? -inf : inf);
          if (!std::isinf(b)) bound = b;
        }
        if (bound.is_empty()) continue;
        output = k == 1 ? fmax(output, bound) : fmin(output, bound);
      }
      return output;
    }
    casadi_error("Unsupported activation '" + op_type + "'");
  }

  bool Onnx::process_shaped_operation(const std::string& op_type, const onnx::NodeProto& node,
                                      const std::vector<MX>& node_inputs,
                                      const IntegerConstants& integer_constants,
                                      MX& output, std::vector<casadi_int>& out_shape) {
    // ONNX shapes of the inputs
    std::vector<std::vector<casadi_int>> s(node_inputs.size());
    bool nd = false;  // an input of rank > 2
    for (size_t j = 0; j < node_inputs.size(); ++j) {
      if (node.input(j).empty()) continue;
      s[j] = onnx_shape(node.input(j), node_inputs[j]);
      nd = nd || s[j].size() > 2;
    }
    auto ints = [&](size_t j) -> std::vector<casadi_int> {
      auto it = integer_constants.find(node.input(j));
      if (it == integer_constants.end()) return constant_ints(node_inputs[j]);
      return std::vector<casadi_int>(it->second.begin(), it->second.end());
    };
    // The axes of Squeeze/Unsqueeze: an input from opset 13, an attribute before
    auto axes_of = [&](std::vector<casadi_int>& axes) -> bool {
      if (node_inputs.size() > 1 && !node.input(1).empty()) {
        axes = ints(1);
        return true;
      }
      for (const auto& attr : node.attribute()) {
        if (attr.name() == "axes") {
          axes.assign(attr.ints().begin(), attr.ints().end());
          return true;
        }
      }
      return false;
    };
    // An axis counted from the end
    auto axis_in = [&](casadi_int a, size_t rank) -> casadi_int {
      casadi_int r = static_cast<casadi_int>(rank);
      casadi_assert(a >= -r && a < r, op_type + ": axis " + str(a) + " out of range for rank "
        + str(r));
      return a < 0 ? a + r : a;
    };

    if (op_type == "Unsqueeze") {
      // Insert dimensions of size 1; the elements stay in order
      std::vector<casadi_int> axes;
      casadi_assert(axes_of(axes), "Unsqueeze requires axes");
      size_t rank = s[0].size() + axes.size();
      std::vector<bool> inserted(rank, false);
      for (casadi_int a : axes) inserted[axis_in(a, rank)] = true;
      auto in = s[0].begin();
      for (size_t k = 0; k < rank; ++k) out_shape.push_back(inserted[k] ? 1 : *in++);
      output = as_shape(node_inputs[0], out_shape);
      return true;

    } else if (op_type == "Squeeze") {
      // Remove dimensions of size 1 (all of them, if no axes are given)
      std::vector<casadi_int> axes;
      std::vector<casadi_int> in = s[0];
      bool has_axes = axes_of(axes);
      // A vector of unknown rank may be a 2-D tensor with one dimension 1
      if (has_axes && in.size() == 1 && !shapes_.count(node.input(0))) {
        for (casadi_int a : axes) {
          if ((a != 0 && a != -1) || in[0] != 1) {
            in = {node_inputs[0].size2(), node_inputs[0].size1()};
          }
        }
      }
      std::vector<bool> removed(in.size(), !has_axes);
      for (casadi_int a : axes) {
        casadi_int k = axis_in(a, in.size());
        casadi_assert(in[k] == 1, "Squeeze: dimension " + str(k) + " of " + shape_str(in)
          + " is not 1");
        removed[k] = true;
      }
      for (size_t k = 0; k < in.size(); ++k) {
        if (!(removed[k] && in[k] == 1)) out_shape.push_back(in[k]);
      }
      output = as_shape(node_inputs[0], out_shape);
      return true;
    }

    if (op_type == "Reshape") {
      // The elements stay in order. A target with rank <= 2, from a value of rank <= 2, is
      // left to process_node_operation.
      casadi_assert(node_inputs.size() >= 2 && node_inputs[1].is_constant(),
                    "Reshape shape must be a constant");
      std::vector<casadi_int> t = ints(1);
      if (!nd && t.size() <= 2) return false;
      bool allowzero = get_int_attribute(node, "allowzero", 0) != 0;
      casadi_int known = 1, infer = -1;
      for (size_t k = 0; k < t.size(); ++k) {
        if (t[k] == 0 && !allowzero) {
          casadi_assert(k < s[0].size(), "Reshape: no dimension " + str(k) + " to copy");
          t[k] = s[0][k];
        }
        if (t[k] == -1) {
          infer = static_cast<casadi_int>(k);
        } else {
          known *= t[k];
        }
      }
      if (infer >= 0) t[infer] = known ? node_inputs[0].numel() / known : 0;
      out_shape = t;
      output = as_shape(densify(node_inputs[0]), t);
      return true;
    }

    if (op_type == "Slice" && s[0].size() != 2) {
      // A rank-1 or rank > 2 tensor: the elements in the ranges along the axes (rank 2 is
      // left to process_node_operation)
      const std::vector<casadi_int>& d = s[0];
      casadi_assert(!d.empty(), "Slice of a scalar");
      std::vector<casadi_int> starts = ints(1), ends = ints(2), axes, steps;
      if (node_inputs.size() > 3 && !node.input(3).empty()) {
        axes = ints(3);
      } else {
        for (size_t k = 0; k < starts.size(); ++k) axes.push_back(k);
      }
      if (node_inputs.size() > 4 && !node.input(4).empty()) {
        steps = ints(4);
      } else {
        steps.assign(starts.size(), 1);
      }
      casadi_assert(ends.size() == starts.size() && axes.size() == starts.size()
        && steps.size() == starts.size(), "Slice: starts, ends, axes and steps must have "
        "equal lengths");
      std::vector<std::vector<casadi_int>> sel(d.size());
      for (size_t k = 0; k < d.size(); ++k) {
        for (casadi_int i = 0; i < d[k]; ++i) sel[k].push_back(i);
      }
      for (size_t a = 0; a < axes.size(); ++a) {
        casadi_int k = axis_in(axes[a], d.size()), n = d[k], st = steps[a];
        casadi_assert(st != 0, "Slice: step must not be zero");
        // Negative positions count from the end; then clamp as ONNX does
        casadi_int b = starts[a] < 0 ? starts[a] + n : starts[a];
        casadi_int e = ends[a] < 0 ? ends[a] + n : ends[a];
        b = st > 0 ? std::min(std::max(b, casadi_int(0)), n)
                   : std::min(std::max(b, casadi_int(-1)), n - 1);
        e = st > 0 ? std::min(std::max(e, casadi_int(0)), n)
                   : std::min(std::max(e, casadi_int(-1)), n - 1);
        sel[k].clear();
        for (casadi_int i = b; st > 0 ? i < e : i > e; i += st) sel[k].push_back(i);
      }
      for (const auto& v : sel) out_shape.push_back(v.size());
      output = take_elements(node_inputs[0], selected_positions(sel, shape_strides(d)),
                             out_shape);
      // A rank-1 row stays a row
      if (d.size() == 1 && node_inputs[0].size1() == 1) output = reshape(output, 1, out_shape[0]);
      return true;
    }

    if (op_type == "Concat" && s[0].size() != 2) {
      // Rank 1 or > 2: the inputs' elements in order along the axis (rank 2 is left to
      // process_node_operation)
      casadi_int axis = get_int_attribute(node, "axis", 0);
      if (!nd && axis != 0 && axis != -1) return false;  // a vector of unknown rank
      axis = axis_in(axis, s[0].size());
      out_shape = s[0];
      out_shape[axis] = 0;
      std::vector<MX> flat;
      std::vector<casadi_int> offset = {0};
      bool rows = true;
      for (size_t j = 0; j < node_inputs.size(); ++j) {
        casadi_assert(s[j].size() == s[0].size(), "Concat: inputs of different ranks");
        out_shape[axis] += s[j][axis];
        flat.push_back(vec(densify(node_inputs[j])));
        offset.push_back(offset.back() + node_inputs[j].numel());
        rows = rows && node_inputs[j].size1() == 1;
      }
      casadi_int outer = 1, inner = 1;
      for (casadi_int k = 0; k < axis; ++k) outer *= out_shape[k];
      for (size_t k = axis + 1; k < out_shape.size(); ++k) inner *= out_shape[k];
      std::vector<casadi_int> pos;
      for (casadi_int o = 0; o < outer; ++o) {
        for (size_t j = 0; j < node_inputs.size(); ++j) {
          for (casadi_int t = 0; t < s[j][axis] * inner; ++t) {
            pos.push_back(offset[j] + o * s[j][axis] * inner + t);
          }
        }
      }
      output = take_elements(vertcat(flat), pos, out_shape);
      // Rank-1 rows stay a row
      if (out_shape.size() == 1 && rows) output = reshape(output, 1, out_shape[0]);
      return true;
    }

    if (op_type == "Gather") {
      // out[a.., i.., b..] = data[a.., indices[i..], b..] along axis
      const std::vector<casadi_int>& d = s[0];
      const std::vector<casadi_int>& is = s[1];
      if (d.size() < 3 && d.size() + is.size() < 4) return false;  // 2-D in and out
      casadi_int axis = axis_in(get_int_attribute(node, "axis", 0), d.size());
      casadi_assert(is_known(node_inputs[1]), "Gather indices must be constant");
      std::vector<casadi_int> idx = ints(1);
      casadi_int outer = 1, inner = 1;
      for (casadi_int k = 0; k < axis; ++k) outer *= d[k];
      for (size_t k = axis + 1; k < d.size(); ++k) inner *= d[k];
      out_shape.assign(d.begin(), d.begin() + axis);
      out_shape.insert(out_shape.end(), is.begin(), is.end());
      out_shape.insert(out_shape.end(), d.begin() + axis + 1, d.end());
      std::vector<casadi_int> pos;
      for (casadi_int o = 0; o < outer; ++o) {
        for (casadi_int t : idx) {
          if (t < 0) t += d[axis];
          casadi_assert(t >= 0 && t < d[axis], "Gather: index " + str(t) + " out of range");
          for (casadi_int i = 0; i < inner; ++i) pos.push_back((o * d[axis] + t) * inner + i);
        }
      }
      output = take_elements(node_inputs[0], pos, out_shape);
      return true;
    }

    if (op_type == "GatherND" || op_type == "ScatterND") {
      // Index tuples: the last axis of indices holds q indices into the leading axes of data,
      // each selecting a block of the trailing axes. The selected blocks of data, in row-major
      // positions, are the output of GatherND and the elements ScatterND overwrites.
      std::vector<casadi_int> d = s[0];
      const std::vector<casadi_int>& is = s[1];
      casadi_assert(get_int_attribute(node, "batch_dims", 0) == 0,
        "ONNX import: GatherND with batch_dims is not supported");
      casadi_assert(!is.empty() && is_known(node_inputs[1]),
        "ONNX import: " + op_type + " requires constant indices of rank >= 1");
      casadi_int q = is.back();
      // A vector of unknown rank may be a 2-D tensor with one dimension 1
      if (static_cast<casadi_int>(d.size()) < q && d.size() <= 1) {
        d = {node_inputs[0].size2(), node_inputs[0].size1()};
      }
      casadi_assert(q >= 1 && q <= static_cast<casadi_int>(d.size()), op_type + ": "
        + str(q) + " indices per tuple for data of shape " + shape_str(d));
      std::vector<casadi_int> idx = ints(1), st = shape_strides(d);
      casadi_int block = 1;
      for (size_t k = q; k < d.size(); ++k) block *= d[k];
      std::vector<casadi_int> pos;
      for (size_t t = 0; t + q <= idx.size(); t += q) {
        casadi_int p = 0;
        for (casadi_int k = 0; k < q; ++k) {
          casadi_int i = idx[t + k] < 0 ? idx[t + k] + d[k] : idx[t + k];
          casadi_assert(i >= 0 && i < d[k], op_type + ": index " + str(idx[t + k])
            + " out of range for dimension " + str(d[k]));
          p += i * st[k];
        }
        for (casadi_int i = 0; i < block; ++i) pos.push_back(p + i);
      }
      if (op_type == "GatherND") {
        out_shape.assign(is.begin(), is.end() - 1);
        out_shape.insert(out_shape.end(), d.begin() + q, d.end());
        output = take_elements(node_inputs[0], pos, out_shape);
        return true;
      }
      // ScatterND: data with the updates (of shape indices[:-1] + data[q:]) written at pos.
      // Element positions are nonzero positions of the dense stored form.
      casadi_assert(static_cast<casadi_int>(pos.size()) == node_inputs[2].numel(),
        "ScatterND: " + str(node_inputs[2].numel()) + " updates for " + str(pos.size())
        + " positions");
      MX data = densify(node_inputs[0]), upd = vec(densify(node_inputs[2]));
      if (get_string_attribute(node, "reduction") == "add") {
        // ADDNONZEROS: out = data with out.nz[pos] += updates, accumulating DUPLICATE pos. An
        // addnonzeros MX cannot be built directly; it is the reverse-mode adjoint of a
        // getnonzeros. Build g = w.nz[pos] on a dense numel-vector w, then
        // jtimes(g, w, updates, true) is exactly the scatter-add of updates into a zero frame
        // (duplicates summed); add it to data.
        MX w = MX::sym("addnz_w", data.numel(), 1), g;
        w.get_nz(g, false, Matrix<casadi_int>(pos));
        output = data + reshape(jtimes(g, w, upd, true), data.size1(), data.size2());
      } else {
        casadi_assert(get_string_attribute(node, "reduction").empty()
          || get_string_attribute(node, "reduction") == "none",
          "ScatterND: reduction '" + get_string_attribute(node, "reduction")
          + "' is not supported");
        output = data;
        output.set_nz(upd, false, Matrix<casadi_int>(pos));
      }
      out_shape = d;
      return true;
    }

    if (!nd) return false;

    if (is_elementwise(op_type)) {
      // Broadcast the tensor operands to a common shape, then the 2-D op applies as is
      std::vector<MX> args = node_inputs;
      for (size_t j = 0; j < s.size(); ++j) {
        if (node.input(j).empty() || (op_type == "Clip" && j > 0)) continue;
        out_shape = out_shape.empty() && j == 0 ? s[j] : broadcast_shapes(out_shape, s[j]);
      }
      for (size_t j = 0; j < s.size(); ++j) {
        if (node.input(j).empty() || (op_type == "Clip" && j > 0)) continue;
        args[j] = broadcast_to(node_inputs[j], s[j], out_shape);
      }
      output = process_node_operation(op_type, node, args, integer_constants);
      return true;

    } else if (op_type == "Transpose") {
      // Permute the axes (default: reverse them)
      const std::vector<casadi_int>& d = s[0];
      std::vector<casadi_int> perm;
      for (const auto& attr : node.attribute()) {
        if (attr.name() == "perm") perm.assign(attr.ints().begin(), attr.ints().end());
      }
      if (perm.empty()) {
        for (casadi_int k = static_cast<casadi_int>(d.size()) - 1; k >= 0; --k) perm.push_back(k);
      }
      casadi_assert(perm.size() == d.size(), "Transpose: perm does not match the rank");
      for (casadi_int p : perm) out_shape.push_back(d[p]);
      std::vector<casadi_int> ds = shape_strides(d), os = shape_strides(out_shape);
      std::vector<casadi_int> pos(shape_numel(d));
      for (casadi_int j = 0; j < static_cast<casadi_int>(pos.size()); ++j) {
        casadi_int p = 0;
        for (size_t k = 0; k < perm.size(); ++k) p += ((j / os[k]) % out_shape[k]) * ds[perm[k]];
        pos[j] = p;
      }
      output = take_elements(node_inputs[0], pos, out_shape);
      return true;
    }

    casadi_error("ONNX import: " + op_type + " is not supported for tensors of rank > 2");
    return false;
  }

  MX Onnx::process_node_operation(
      const std::string& op_type,
      const onnx::NodeProto& node,
      const std::vector<MX>& node_inputs,
      const IntegerConstants& integer_constants) {

    auto input_ints = [&](size_t i) -> std::vector<casadi_int> {
      auto it = integer_constants.find(node.input(i));
      if (it == integer_constants.end()) return constant_ints(node_inputs[i]);
      std::vector<casadi_int> values;
      for (int64_t v : it->second) {
        casadi_assert(v >= std::numeric_limits<casadi_int>::min() &&
          v <= std::numeric_limits<casadi_int>::max(), "Integer constant out of range: " + str(v));
        values.push_back(static_cast<casadi_int>(v));
      }
      return values;
    };

    MX output;

    // Simple ops via the centralized lookup table
    const OpMapping* mapping = get_op_mapping_by_name(op_type);
    if (mapping) {
      const MX& x = node_inputs[0];
      if (mapping->arity == 1) {
        casadi_assert(node_inputs.size() >= 1, op_type + " requires 1 input");
        switch (mapping->casadi_op) {
          case OP_SIN: return sin(x);
          case OP_COS: return cos(x);
          case OP_TAN: return tan(x);
          case OP_ASIN: return asin(x);
          case OP_ACOS: return acos(x);
          case OP_ATAN: return atan(x);
          case OP_SINH: return sinh(x);
          case OP_COSH: return cosh(x);
          case OP_TANH: return tanh(x);
          case OP_ASINH: return asinh(x);
          case OP_ACOSH: return acosh(x);
          case OP_ATANH: return atanh(x);
          case OP_EXP: return exp(x);
          case OP_LOG: return log(x);
          case OP_SQRT: return sqrt(x);
          case OP_NEG: return -x;
          case OP_FABS: return fabs(x);
          case OP_CEIL: return ceil(x);
          case OP_FLOOR: return floor(x);
          case OP_SIGN: return sign(x);
          case OP_ERF: return erf(x);
          case OP_INV: return 1.0 / x;
          case OP_TRANSPOSE: return x.T();
          case OP_NORM1: return norm_1(x);
          case OP_NORM2: return norm_2(x);
          case OP_NORMF: return norm_fro(x);
          case OP_MMIN: return mmin(x);
          case OP_MMAX: return mmax(x);
          default: break;
        }
      } else if (mapping->arity == 2) {
        casadi_assert(node_inputs.size() >= 2, op_type + " requires 2 inputs");
        MX bx = x, by = node_inputs[1];
        broadcast_pair(&bx, &by);
        switch (mapping->casadi_op) {
          case OP_ADD: return bx + by;
          case OP_SUB: return bx - by;
          case OP_MUL: return bx * by;
          case OP_DIV: return bx / by;
          case OP_POW: return pow(bx, by);
          default: break;
        }
      }
    }

    // Everything else needs attributes, multiple inputs, etc.

    if (op_type == "Less") {
      casadi_assert(node_inputs.size() >= 2, "Less requires 2 inputs");
      output = if_else(node_inputs[0] < node_inputs[1], MX(1.0), MX(0.0));

    } else if (op_type == "Equal") {
      casadi_assert(node_inputs.size() >= 2, "Equal requires 2 inputs");
      output = !ne(node_inputs[0], node_inputs[1]);  // CasADi has ne() but no eq()

    } else if (op_type == "LessOrEqual") {
      casadi_assert(node_inputs.size() >= 2, "LessOrEqual requires 2 inputs");
      output = if_else(node_inputs[0] <= node_inputs[1], MX(1.0), MX(0.0));

    } else if (op_type == "Greater") {
      casadi_assert(node_inputs.size() >= 2, "Greater requires 2 inputs");
      output = if_else(node_inputs[0] > node_inputs[1], MX(1.0), MX(0.0));

    } else if (op_type == "GreaterOrEqual") {
      casadi_assert(node_inputs.size() >= 2, "GreaterOrEqual requires 2 inputs");
      output = if_else(node_inputs[0] >= node_inputs[1], MX(1.0), MX(0.0));

    } else if (op_type == "Min") {
      casadi_assert(node_inputs.size() >= 2, "Min requires 2 inputs");
      output = fmin(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Max") {
      casadi_assert(node_inputs.size() >= 2, "Max requires 2 inputs");
      output = fmax(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Mod") {
      casadi_assert(node_inputs.size() >= 2, "Mod requires 2 inputs");
      output = fmod(node_inputs[0], node_inputs[1]);

    } else if (op_type == "ReduceSum") {
      // Sum all elements (ReduceMin/Max/L1/L2 are in the op_map table)
      casadi_assert(node_inputs.size() >= 1, "ReduceSum requires 1 input");
      output = sum1(sum2(node_inputs[0]));

    } else if (op_type == "Not") {
      casadi_assert(node_inputs.size() >= 1, "Not requires 1 input");
      output = logic_not(node_inputs[0]);

    } else if (op_type == "And") {
      casadi_assert(node_inputs.size() >= 2, "And requires 2 inputs");
      output = logic_and(node_inputs[0], node_inputs[1]);

    } else if (op_type == "Or") {
      casadi_assert(node_inputs.size() >= 2, "Or requires 2 inputs");
      output = logic_or(node_inputs[0], node_inputs[1]);
      // Or with a constant true operand is known (as And with a constant false one already is),
      // so that an If on it resolves to one branch
      for (const MX& a : {node_inputs[0], node_inputs[1]}) {
        if (!output.is_constant() && a.is_constant() && a.is_scalar(true)
            && densify(static_cast<DM>(a)).scalar() != 0) {
          output = MX::ones(output.size());
        }
      }

    } else if (op_type == "Where") {
      casadi_assert(node_inputs.size() >= 3, "Where requires 3 inputs");
      output = if_else(node_inputs[0], node_inputs[1], node_inputs[2]);

    } else if (op_type == "Identity") {
      casadi_assert(node_inputs.size() >= 1, "Identity requires 1 input");
      output = node_inputs[0];

    } else if (op_type == "Cast") {
      // Everything is double in CasADi (converted at the tensor_to_dm boundary), so Cast is an
      // identity during symbolic computation; the 'to' attribute is irrelevant here.
      casadi_assert(node_inputs.size() >= 1, "Cast requires 1 input");
      output = node_inputs[0];

    } else if (op_type == "MatMul") {
      // Transpose-rep: ONNX MatMul(P,Q)=P*Q is the transpose-rep of q*p, so swap operands.
      casadi_assert(node_inputs.size() >= 2, "MatMul requires 2 inputs");
      output = mtimes(node_inputs[1], node_inputs[0]);

    } else if (op_type == "Gemm") {
      // Transpose-rep: the stored result is (opQ(Q)^T)*(opP(P)^T); recover the CasADi value by
      // swapping operands and the transA/transB roles.
      casadi_assert(node_inputs.size() >= 2, "Gemm requires at least 2 inputs");
      MX A = get_int_attribute(node, "transB", 0) ? node_inputs[1].T() : node_inputs[1];
      MX B = get_int_attribute(node, "transA", 0) ? node_inputs[0].T() : node_inputs[0];
      output = get_float_attribute(node, "alpha", 1.0) * mtimes(A, B);
      if (node_inputs.size() >= 3) {
        MX C = node_inputs[2];  // unidirectionally broadcastable to the result, e.g. a bias [N]
        broadcast_pair(&output, &C);
        output = output + get_float_attribute(node, "beta", 1.0) * C;
      }

    } else if (op_type == "Sum") {
      // Variadic elementwise sum
      casadi_assert(node_inputs.size() >= 1, "Sum requires at least 1 input");
      output = node_inputs[0];
      for (casadi_int idx = 1; idx < node_inputs.size(); ++idx) output = output + node_inputs[idx];

    } else if (op_type == "Pad") {
      // Inverse of the block-diagonal Pad: drop the (dense) block into a zero frame at its offset.
      // Transpose-rep: pads are [axis0_begin, axis1_begin, axis0_end, axis1_end] on the stored
      // (C x R) tensor, i.e. [col_off, row_off, col_end, row_end].
      casadi_assert(node_inputs.size() >= 2, "Pad requires data and pads");
      MX block = node_inputs[0];
      std::vector<casadi_int> pads = constant_ints(node_inputs[1]);
      casadi_int col_off = pads[0], row_off = pads[1];
      casadi_int br = block.size1(), bc = block.size2();
      casadi_int C = col_off + bc + pads[2];
      // Pad with SPARSE zero blocks via horz/vertcat: the block's pattern is preserved and no
      // project/densify is introduced (a triangular block stays triangular), so the block matrix's
      // sparsity emerges natively from CasADi's own propagation.
      MX padded = block;
      if (col_off > 0 || pads[2] > 0) {
        padded = horzcat(MX(Sparsity(br, col_off)), padded, MX(Sparsity(br, pads[2])));
      }
      if (row_off > 0 || pads[3] > 0) {
        padded = vertcat(MX(Sparsity(row_off, C)), padded, MX(Sparsity(pads[3], C)));
      }
      output = padded;

    } else if (op_type == "Einsum") {
      // Reconstruct a CasADi einstein contraction (inverse of the export envelope). The
      // operands arrive as the 2-D reshaped tensors; their column-major vec is the original.
      casadi_assert(node_inputs.size() >= 2, "Einsum requires 2 inputs");
      std::string eq;
      for (char ch : get_string_attribute(node, "equation")) if (ch != ' ') eq += ch;
      size_t comma = eq.find(','), arrow = eq.find("->");
      casadi_assert(comma != std::string::npos && arrow != std::string::npos,
                    "ONNX import: only binary Einsum 'a,b->c' is supported");
      std::string sa = eq.substr(0, comma);
      std::string sb = eq.substr(comma + 1, arrow - comma - 1);
      std::string sc = eq.substr(arrow + 2);

      // Transpose-rep: the exporter reversed each subscript to label the axis-reversed ONNX
      // tensors. node_inputs here are the imported CasADi values (transposed back to LOGICAL axis
      // order), so reverse the subscripts back to natural order to pair them with the logical axes
      // (size1,size2).
      std::reverse(sa.begin(), sa.end());
      std::reverse(sb.begin(), sb.end());
      std::reverse(sc.begin(), sc.end());

      // Each letter's size, read from the operand shapes (axes follow the subscript order). A
      // single-index (vector) operand is stored as a rank-1 ONNX tensor that imports with an
      // ambiguous 1xN/Nx1 orientation, so take its length from numel() (orientation-independent);
      // a 2-index operand reads size1/size2 directly.
      std::map<char, casadi_int> lsize;
      for (int t = 0; t < 2; ++t) {
        const std::string& s = (t == 0) ? sa : sb;
        const MX& m = node_inputs[t];
        if (s.size() == 1) {
          lsize[s[0]] = m.numel();
        } else if (s.size() >= 2) {
          lsize[s[0]] = m.size1();
          lsize[s[1]] = m.size2();
        }
      }
      // Assign a label to each distinct letter
      std::map<char, casadi_int> lab;
      casadi_int next_label = -1;
      for (char ch : sa + sb + sc) if (!lab.count(ch)) lab[ch] = next_label--;

      // Labels in natural order; operands' column-major vec is the original einstein input
      std::vector<casadi_int> da, db, dc, La, Lb, Lc;
      for (char ch : sa) { da.push_back(lsize[ch]); La.push_back(lab[ch]); }
      for (char ch : sb) { db.push_back(lsize[ch]); Lb.push_back(lab[ch]); }
      for (char ch : sc) { dc.push_back(lsize[ch]); Lc.push_back(lab[ch]); }
      output = einstein(vec(node_inputs[0]), vec(node_inputs[1]), da, db, dc, La, Lb, Lc);

    } else if (op_type == "Det") {
      casadi_assert(node_inputs.size() >= 1, "Det requires 1 input");
      output = det(node_inputs[0]);

    } else if (op_type == "ReduceLogSumExp") {
      // log(sum(exp(x))) over all elements
      casadi_assert(node_inputs.size() >= 1, "ReduceLogSumExp requires 1 input");
      output = log(sum1(sum2(exp(node_inputs[0]))));

    } else if (op_type == "Constant") {
      // ONNX permits tensor, sparse tensor, and scalar/vector numeric attributes.
      const onnx::AttributeProto* value_attr = nullptr;
      for (int a = 0; a < node.attribute_size(); ++a) {
        const auto& attr = node.attribute(a);
        const std::string& an = attr.name();
        if (an == "value" || an == "sparse_value" || an == "value_int" ||
            an == "value_ints" || an == "value_float" || an == "value_floats") {
          casadi_assert(value_attr == nullptr, "Constant node has multiple value attributes");
          value_attr = &attr;
        }
      }
      casadi_assert(value_attr != nullptr,
        "Constant node requires a numeric value or sparse_value attribute");
      const auto& attr = *value_attr;
      if (attr.name() == "sparse_value") {
        output = MX(sparse_tensor_to_dm(attr.sparse_tensor()));
      } else if (attr.name() == "value") {
        output = MX(tensor_to_dm(attr.t()));
      } else if (attr.name() == "value_int") {
        output = MX(static_cast<double>(attr.i()));
      } else if (attr.name() == "value_float") {
        output = MX(static_cast<double>(attr.f()));
      } else {
        std::vector<double> values;
        if (attr.name() == "value_ints") {
          for (auto v : attr.ints()) values.push_back(static_cast<double>(v));
        } else {
          for (auto v : attr.floats()) values.push_back(static_cast<double>(v));
        }
        output = MX(DM(values));
      }

    // Complex tensor operations (Transpose is handled by the op_map table)
    } else if (op_type == "Reshape") {
      casadi_assert(node_inputs.size() >= 2,
                    "Reshape operation requires 2 inputs (data and shape)");
      // Second input is the target shape - should be a constant
      casadi_assert(node_inputs[1].is_constant(),
                    "Reshape shape must be a constant");
      std::vector<casadi_int> shape = input_ints(1);
      // ONNX Reshape is PSEUDO-DENSE (it operates on the full numel), so densify the operand here:
      // this is the ONE place the imported graph must drop sparsity, because CasADi's own reshape
      // PRESERVES it (a sparse operand would reshape to a flat with only nnz entries, shifting
      // every
      // downstream Gather/Scatter index). True sparsity propagates everywhere else.
      // Transpose-rep: the ONNX target (s0,s1) is the reverse of the CasADi target. A target of
      // rank > 2 is handled by process_shaped_operation.
      casadi_int s0 = shape.empty() ? 1 : shape[0];
      casadi_int s1 = shape.size() > 1 ? shape[1] : 1;
      output = reshape(densify(node_inputs[0]), s1, s0);

    } else if (op_type == "Concat") {
      // Transpose-rep: ONNX axis 0 is CasADi columns (horzcat), axis 1 is CasADi rows (vertcat)
      casadi_int axis = get_int_attribute(node, "axis", 0);
      if (axis == 0) {
        output = horzcat(node_inputs);
      } else if (axis == 1) {
        output = vertcat(node_inputs);
      } else {
        casadi_error("Concat with axis=" + std::to_string(axis) +
                     " not supported. Only axis=0 (vertcat) and axis=1 (horzcat) are supported.");
      }

    } else if (op_type == "Slice") {
      // Inputs: data, starts, ends, [axes], [steps].
      casadi_assert(node_inputs.size() >= 3,
                    "Slice requires at least 3 inputs (data, starts, ends)");

      MX data = node_inputs[0];
      casadi_assert(node_inputs[1].is_constant() && node_inputs[2].is_constant(),
                    "Slice starts and ends must be constants");

      DM starts_dm = static_cast<DM>(node_inputs[1]);
      DM ends_dm = static_cast<DM>(node_inputs[2]);

      // Axes (default [0, 1, ...]) and steps (default all 1s)
      std::vector<casadi_int> axes, steps;
      if (node_inputs.size() >= 4 && !node_inputs[3].is_empty()) {
        axes = input_ints(3);
      } else {
        for (casadi_int k = 0; k < starts_dm.numel(); ++k) axes.push_back(k);
      }
      if (node_inputs.size() >= 5 && !node_inputs[4].is_empty()) {
        steps = input_ints(4);
      } else {
        steps.assign(starts_dm.numel(), 1);
      }

      casadi_assert(axes.size() <= 2, "Slice: only up to 2D slicing supported");
      casadi_assert(starts_dm.numel() == axes.size() && ends_dm.numel() == axes.size() &&
        steps.size() == axes.size(), "Slice: starts, ends, axes and steps must have equal lengths");

      // Transpose-rep: ONNX axis 0 is CasADi columns, axis 1 is CasADi rows.
      Slice row_slice, col_slice;
      casadi_int nrow = data.size1(), ncol = data.size2();
      for (casadi_int a = 0; a < static_cast<casadi_int>(axes.size()); ++a) {
        casadi_int axis = axes[a] < 0 ? axes[a] + 2 : axes[a];
        casadi_assert(axis == 0 || axis == 1, "Slice: axis out of range");
        casadi_int n = axis == 0 ? ncol : nrow;
        casadi_int sp = steps[a];
        casadi_assert(sp != 0, "Slice: step must not be zero");
        sp = std::max(-std::max(n, casadi_int(1)), std::min(sp, std::max(n, casadi_int(1))));
        auto bound = [&](size_t input, const DM& dm) -> casadi_int {
          auto it = integer_constants.find(node.input(input));
          if (it != integer_constants.end()) {
            int64_t v = it->second.at(a);
            if (v < 0) v += n;
            return static_cast<casadi_int>(std::max<int64_t>(sp > 0 ? 0 : -1,
              std::min<int64_t>(v, sp > 0 ? n : n - 1)));
          }
          // Preserve the existing path for computed constants; clip before narrowing.
          double v = dm(a).scalar();
          if (v < 0) v += n;
          return static_cast<casadi_int>(std::max(sp > 0 ? 0. : -1.,
            std::min(v, static_cast<double>(sp > 0 ? n : n - 1))));
        };
        casadi_int st = bound(1, starts_dm);
        casadi_int en = bound(2, ends_dm);
        // CasADi uses an unbounded stop to represent -1 for a backwards slice.
        Slice indices = (sp > 0 ? st >= en : st <= en) ? Slice(0, 0) :
          Slice(st, en == -1 ? std::numeric_limits<casadi_int>::max() : en, sp);
        if (axis == 0) {
          col_slice = indices;
        } else {
          row_slice = indices;
        }
      }
      output = data(row_slice, col_slice);

    } else if (op_type == "Gather") {
      // Extract element(s) at the given indices along an axis.
      casadi_assert(node_inputs.size() >= 2, "Gather requires data and indices");

      // Gather indexes the pseudo-dense column-major layout -> densify a sparse operand.
      MX data = densify(node_inputs[0]);
      MX indices_mx = node_inputs[1];
      casadi_int axis = get_int_attribute(node, "axis", 0);

      casadi_assert(is_known(indices_mx), "Gather indices must be constant");
      std::vector<casadi_int> indices = input_ints(1);

      // Keep even scalar indices in an integer vector: Matrix's scalar constructor takes double.
      if (data.size2() == 1 || axis == 1) {
        output = data(indices, Slice());
      } else if (axis == 0) {
        output = data(Slice(), indices);
      } else {
        casadi_error("Gather: only axis 0 and 1 supported for 2D tensors");
      }

    } else if (op_type == "ScatterElements") {
      // data with data.nz[indices] {=,+=} updates -- nz indexing is pseudo-dense, so densify the
      // target. Constant indices -> (set/add)nonzeros; runtime indices -> the _PARAM variants.
      casadi_assert(node_inputs.size() >= 3, "ScatterElements requires data, indices, updates");
      MX se_data = densify(node_inputs[0]);
      MX se_idx = node_inputs[1];
      MX se_upd = vec(node_inputs[2]);
      if (get_string_attribute(node, "reduction") == "add") {
        // {add,addparam}nonzeros: out = data with data.nz[idx] += updates, accumulating DUPLICATE
        // idx. No direct MX builder; it is the reverse-mode adjoint of a (parametric)
        // getnonzeros --
        // jtimes(w.nz[idx], w, updates, true) is the scatter-add (dups summed) into a zero frame.
        MX w = MX::sym("scel_w", se_data.numel(), 1), g;
        if (se_idx.is_constant()) {
          w.get_nz(g, false, Matrix<casadi_int>(constant_ints(se_idx)));
        } else {
          w.get_nz(g, false, vec(se_idx));   // runtime indices -> getnonzeros_param
        }
        output = se_data + reshape(jtimes(g, w, se_upd, true), se_data.size1(), se_data.size2());
      } else {
        output = se_data;
        if (se_idx.is_constant()) {
          output.set_nz(se_upd, false, Matrix<casadi_int>(constant_ints(se_idx)));
        } else {
          output.set_nz(se_upd, false, vec(se_idx));   // runtime indices -> setnonzeros_param
        }
      }

    } else if (op_type == "Tile") {
      // Repeat the tensor along each dimension by repeats = [rows, cols].
      casadi_assert(node_inputs.size() >= 2, "Tile requires data and repeats inputs");
      MX data = node_inputs[0];
      MX repeats_mx = node_inputs[1];
      casadi_assert(repeats_mx.is_constant(), "Tile repeats must be constant");
      DM repeats_dm = static_cast<DM>(repeats_mx);

      casadi_int rows_repeat = 1, cols_repeat = 1;
      if (repeats_dm.numel() >= 1) rows_repeat = static_cast<casadi_int>(repeats_dm(0).scalar());
      if (repeats_dm.numel() >= 2) cols_repeat = static_cast<casadi_int>(repeats_dm(1).scalar());

      // Transpose-rep: ONNX repeats [rows,cols] map to CasADi repmat(data, cols, rows)
      output = repmat(data, cols_repeat, rows_repeat);

    } else if (op_type == "GatherElements") {
      // out[..,j,..] = data[.., idx[..,j,..], ..] along `axis`. The export of OP_GETNONZEROS_PARAM
      // emits this on a 1 x N flat row with a 1 x M index row (-> imported as N x 1 / M x 1 cols),
      // which is exactly the parametric getnonzeros y = data.nz[idx]. Reconstruct via get_nz so the
      // sparse/runtime composition (loc-mapping gather then data gather) rebuilds the right MX.
      casadi_assert(node_inputs.size() >= 2, "GatherElements requires data and indices");
      MX data = densify(node_inputs[0]);  // nz indexing is pseudo-dense -> densify the operand
      MX indices_mx = node_inputs[1];
      // Orient indices as a column vector of nz-indices (a 1 x M ONNX row imports as M x 1).
      MX idx_vec = vec(indices_mx);
      MX flat = vec(data);  // column-major flat of data: flat.nz[k] = data's k-th nonzero
      if (indices_mx.is_constant()) {
        // Constant indices -> the plain constant getnonzeros (OP_GETNONZEROS) path.
        flat.get_nz(output, false, Matrix<casadi_int>(constant_ints(indices_mx)));
      } else {
        // Runtime indices -> parametric getnonzeros (OP_GETNONZEROS_PARAM).
        flat.get_nz(output, false, idx_vec);
      }

    } else if (is_activation(op_type)) {
      output = activation(op_type, node, node_inputs);

    } else {
      casadi_error("Unsupported operation '" + op_type + "'");
    }

    return output;
  }

} // namespace casadi
/// \endcond
