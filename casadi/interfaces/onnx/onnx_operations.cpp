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

/// \cond INTERNAL
namespace casadi {

  // ========== Operation Mapping ==========
  // Single source of truth for simple CasADi <-> ONNX ops; shared by export and import.
  // (OpMapping is declared in onnx_model.hpp.)

  static const OpMapping op_map[] = {
    // Unary operations
    {OP_SIN, "Sin", 1},
    {OP_COS, "Cos", 1},
    {OP_TAN, "Tan", 1},
    {OP_ASIN, "Asin", 1},
    {OP_ACOS, "Acos", 1},
    {OP_ATAN, "Atan", 1},
    {OP_SINH, "Sinh", 1},
    {OP_COSH, "Cosh", 1},
    {OP_TANH, "Tanh", 1},
    {OP_ASINH, "Asinh", 1},
    {OP_ACOSH, "Acosh", 1},
    {OP_ATANH, "Atanh", 1},
    {OP_EXP, "Exp", 1},
    {OP_LOG, "Log", 1},
    {OP_SQRT, "Sqrt", 1},
    {OP_NEG, "Neg", 1},
    {OP_FABS, "Abs", 1},
    {OP_CEIL, "Ceil", 1},
    {OP_FLOOR, "Floor", 1},
    {OP_SIGN, "Sign", 1},
    {OP_ERF, "Erf", 1},
    {OP_INV, "Reciprocal", 1},
    {OP_TRANSPOSE, "Transpose", 1},
    {OP_NORM1, "ReduceL1", 1},
    {OP_NORMF, "ReduceL2", 1},  // NORMF first: import uses norm_fro (works for vectors and matrices)
    {OP_NORM2, "ReduceL2", 1},
    {OP_MMIN, "ReduceMin", 1},
    {OP_MMAX, "ReduceMax", 1},
    // Binary operations
    {OP_ADD, "Add", 2},
    {OP_SUB, "Sub", 2},
    {OP_MUL, "Mul", 2},
    {OP_DIV, "Div", 2},
    {OP_POW, "Pow", 2},
    {OP_CONSTPOW, "Pow", 2},
  };

  // Lookup by CasADi opcode (for export)
  const OpMapping* get_op_mapping(casadi_int op) {
    for (const auto& entry : op_map) {
      if (entry.casadi_op == op) return &entry;
    }
    return nullptr;
  }

  // Lookup by ONNX name (for import)
  const OpMapping* get_op_mapping_by_name(const std::string& onnx_name) {
    for (const auto& entry : op_map) {
      if (entry.onnx_name == onnx_name) return &entry;
    }
    return nullptr;
  }

  // ========== Export Helpers ==========

  // Create a Constant node with a tensor value. Setting the attribute type explicitly is
  // required for ONNX Runtime compatibility.
  onnx::TensorProto* create_constant_tensor(
      AddNodeFn add_node,
      const std::string& output_name,
      onnx::TensorProto::DataType data_type) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Constant");
    node->add_output(output_name);
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name("value");
    attr->set_type(onnx::AttributeProto::TENSOR);  // Required for ONNX Runtime
    onnx::TensorProto* tensor = attr->mutable_t();
    tensor->set_data_type(data_type);
    return tensor;
  }

  // Add a Constant node holding an INT64 tensor (default shape: 1-D of data.size())
  void add_int_constant(AddNodeFn add_node, const std::string& name,
                        const std::vector<casadi_int>& data,
                        std::vector<casadi_int> dims = {}) {
    onnx::TensorProto* t = create_constant_tensor(add_node, name, onnx::TensorProto::INT64);
    if (dims.empty()) dims = {static_cast<casadi_int>(data.size())};
    for (casadi_int d : dims) t->add_dims(d);
    for (casadi_int v : data) t->add_int64_data(v);
  }

  // Add a Constant node holding a real tensor in the configured type (empty dims => scalar)
  void Onnx::add_real_constant(AddNodeFn add_node, const std::string& name,
                                         const std::vector<double>& data,
                                         const std::vector<casadi_int>& dims) {
    onnx::TensorProto* t = create_constant_tensor(add_node, name, real_type());
    for (casadi_int d : dims) t->add_dims(d);
    if (real_type() == onnx::TensorProto::FLOAT) {
      for (double v : data) t->add_float_data(static_cast<float>(v));
    } else {
      for (double v : data) t->add_double_data(v);
    }
  }

  // Add a Constant node holding a SPARSE real tensor (ONNX sparse_value attribute, COO format).
  void Onnx::add_sparse_constant(AddNodeFn add_node, const std::string& name, const DM& dm) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Constant");
    node->add_output(name);
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name("sparse_value");
    attr->set_type(onnx::AttributeProto::SPARSE_TENSOR);
    fill_sparse_tensor(attr->mutable_sparse_tensor(), name, dm);
  }

  void Onnx::fill_sparse_tensor(onnx::SparseTensorProto* st, const std::string& name,
                                const DM& dm) const {
    // Transpose-rep: store V^T (shape c x r). V's NATIVE column-major nonzero order is exactly
    // V^T's ascending row-major order, so emit V's (col,row) coordinates directly -- no sort,
    // no transpose.
    st->add_dims(dm.size2());
    st->add_dims(dm.size1());
    casadi_int nnz = dm.nnz();
    std::vector<casadi_int> row = dm.sparsity().get_row();  // V's rows (column-major order)
    std::vector<casadi_int> col = dm.sparsity().get_col();  // V's cols (column-major order)
    const std::vector<double>& vals = dm.nonzeros();
    // values [NNZ]
    onnx::TensorProto* vt = st->mutable_values();
    vt->set_name(name);  // a SparseTensorProto's name is its values' name
    vt->set_data_type(real_type());
    vt->add_dims(nnz);
    if (real_type() == onnx::TensorProto::FLOAT) {
      for (double v : vals) vt->add_float_data(static_cast<float>(v));
    } else {
      for (double v : vals) vt->add_double_data(v);
    }
    // indices [NNZ, 2] = V^T coordinates (col, row), 0-based, ascending row-major of V^T
    onnx::TensorProto* it = st->mutable_indices();
    it->set_data_type(onnx::TensorProto::INT64);
    it->add_dims(nnz);
    it->add_dims(2);
    for (casadi_int k = 0; k < nnz; ++k) {
      it->add_int64_data(col[k]);  // row in V^T (= column of V)
      it->add_int64_data(row[k]);  // column in V^T (= row of V)
    }
  }

  // Add an integer attribute (e.g. axis) to a node
  void add_int_attribute(onnx::NodeProto* node, const std::string& name, casadi_int value) {
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name(name);
    attr->set_type(onnx::AttributeProto::INT);
    attr->set_i(value);
  }

  // Add an integer-list attribute (e.g. scan_input_axes) to a node
  void add_ints_attribute(onnx::NodeProto* node, const std::string& name,
                          const std::vector<casadi_int>& values) {
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name(name);
    attr->set_type(onnx::AttributeProto::INTS);
    for (casadi_int v : values) attr->add_ints(v);
  }

  // GraphProto convenience wrapper for add_int_constant
  void add_int_constant(onnx::GraphProto* graph, const std::string& name,
                        const std::vector<casadi_int>& data, std::vector<casadi_int> dims) {
    add_int_constant([graph]() { return graph->add_node(); }, name, data, dims);
  }

  // Gather(data, indices) along an axis -> output
  onnx::NodeProto* create_gather_node(AddNodeFn add_node, const std::string& data,
                                      const std::string& indices, const std::string& output,
                                      casadi_int axis = 0) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Gather");
    node->add_input(data);
    node->add_input(indices);
    node->add_output(output);
    add_int_attribute(node, "axis", axis);
    return node;
  }

  // Slice(data) on the given axes from starts to ends (step 1), emitting the index constants
  onnx::NodeProto* create_slice_node(AddNodeFn add_node, const std::string& uniq,
      const std::string& data, const std::vector<casadi_int>& starts,
      const std::vector<casadi_int>& ends, const std::vector<casadi_int>& axes,
      const std::vector<casadi_int>& steps, const std::string& output) {
    add_int_constant(add_node, uniq + "_starts", starts);
    add_int_constant(add_node, uniq + "_ends", ends);
    add_int_constant(add_node, uniq + "_axes", axes);
    if (!steps.empty()) add_int_constant(add_node, uniq + "_steps", steps);  // before the Slice node
    onnx::NodeProto* node = add_node();
    node->set_op_type("Slice");
    node->add_input(data);
    node->add_input(uniq + "_starts");
    node->add_input(uniq + "_ends");
    node->add_input(uniq + "_axes");
    if (!steps.empty()) node->add_input(uniq + "_steps");  // omit -> ONNX defaults steps to 1
    node->add_output(output);
    return node;
  }

  // Reshape `data` (which, under the transpose-representation invariant, holds the transpose of a
  // CasADi value) to a CasADi target shape `dims`. Because every ONNX tensor stores its CasADi
  // value's column-major bytes with the shape declared reversed, a CasADi column-major reshape to
  // `dims` is just a plain row-major Reshape to the REVERSED dims -- no Transpose nodes needed.
  void emit_colmajor_reshape(AddNodeFn add_node, const std::string& data,
                             const std::vector<casadi_int>& dims, const std::string& output,
                             const std::string& uniq) {
    std::vector<casadi_int> rev(dims.rbegin(), dims.rend());
    add_int_constant(add_node, uniq + "_s", rev);
    create_unary_node(add_node, "Reshape", data, output)->add_input(uniq + "_s");
  }

  // Callback-based implementations (main logic)
  onnx::NodeProto* create_binary_node(
      AddNodeFn add_node,
      const std::string& op_type,
      const std::string& input1,
      const std::string& input2,
      const std::string& output) {
    onnx::NodeProto* node = add_node();
    node->set_op_type(op_type);
    node->add_input(input1);
    node->add_input(input2);
    node->add_output(output);
    return node;
  }

  onnx::NodeProto* create_unary_node(
      AddNodeFn add_node,
      const std::string& op_type,
      const std::string& input,
      const std::string& output) {
    onnx::NodeProto* node = add_node();
    node->set_op_type(op_type);
    node->add_input(input);
    node->add_output(output);
    return node;
  }

  // Cast(input) -> output of the given ONNX data type
  onnx::NodeProto* create_cast_node(
      AddNodeFn add_node,
      const std::string& input,
      const std::string& output,
      onnx::TensorProto::DataType to_type) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Cast");
    node->add_input(input);
    node->add_output(output);
    onnx::AttributeProto* attr = node->add_attribute();
    attr->set_name("to");
    attr->set_type(onnx::AttributeProto::INT);
    attr->set_i(to_type);
    return node;
  }

  // Gemm: the CALLER intends CasADi-level output = (A or A')*(B or B') [+ C]. Under the
  // transpose-rep invariant the stored output must be out^T = (A*B+C)^T = B^T*A^T + C^T, so we
  // emit Gemm(stored_B, stored_A [, stored_C]) with the transA/transB flags swapped too.
  // C empty -> 2-input form (alpha=beta=1).
  onnx::NodeProto* create_gemm_node(
      AddNodeFn add_node,
      const std::string& A, const std::string& B, const std::string& C,
      const std::string& output, bool transA = false, bool transB = false) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Gemm");
    node->add_input(B);
    node->add_input(A);
    if (!C.empty()) node->add_input(C);
    node->add_output(output);
    if (transB) add_int_attribute(node, "transA", 1);
    if (transA) add_int_attribute(node, "transB", 1);
    return node;
  }

  // Where(condition, if_true, if_false) -> output
  onnx::NodeProto* create_where_node(
      AddNodeFn add_node,
      const std::string& cond,
      const std::string& if_true,
      const std::string& if_false,
      const std::string& output) {
    onnx::NodeProto* node = add_node();
    node->set_op_type("Where");
    node->add_input(cond);
    node->add_input(if_true);
    node->add_input(if_false);
    node->add_output(output);
    return node;
  }

  // GraphProto convenience wrappers
  onnx::NodeProto* create_binary_node(
      onnx::GraphProto* graph,
      const std::string& op_type,
      const std::string& input1,
      const std::string& input2,
      const std::string& output) {
    return create_binary_node([graph]() { return graph->add_node(); },
                               op_type, input1, input2, output);
  }

  onnx::NodeProto* create_unary_node(
      onnx::GraphProto* graph,
      const std::string& op_type,
      const std::string& input,
      const std::string& output) {
    return create_unary_node([graph]() { return graph->add_node(); },
                              op_type, input, output);
  }

  // Emit the general "rearrange nonzeros between two sparsities" envelope (the poorest-form
  // primitive shared by OP_GETNONZEROS's fallback and OP_PROJECT/densify). `idx[k]` = input
  // nonzero index feeding output nonzero k, or -1 = fill-with-zero.
  void Onnx::emit_nonzero_remap(AddNodeFn add_node, const std::string& data,
                                const Sparsity& sp_in, const Sparsity& out_sp,
                                std::vector<casadi_int> idx, const std::string& uniq,
                                const std::string& node_output) {
    casadi_int numel_in = sp_in.size1() * sp_in.size2();
    std::vector<casadi_int> loc = sp_in.find();
    bool has_fill = false;
    for (casadi_int& v : idx) {
      if (v < 0) { v = numel_in; has_fill = true; }  // -> the appended zero (below)
      else v = loc[v];
    }

    // Flatten the input to a 1 x numel ROW vector (skip when it is already a row)
    std::string flat = data;
    if (sp_in.size2() != 1) {
      flat = "gnz_flat_" + uniq;
      emit_colmajor_reshape(add_node, data, {numel_in, 1}, flat, "gnz_flat_" + uniq);
    }
    // Zero-fill (-1) entries gather position `numel_in`: append a single 0 to the flat row.
    if (has_fill) {
      std::string zc = "gnz_fz_" + uniq, fl2 = "gnz_flz_" + uniq;
      add_real_constant(add_node, zc, {0.0}, {1, 1});
      onnx::NodeProto* cc = add_node();
      cc->set_op_type("Concat");
      cc->add_input(flat); cc->add_input(zc); cc->add_output(fl2);
      add_int_attribute(cc, "axis", 1);
      flat = fl2;
    }

    // Gather the nnz_out values along axis 1 -> a (1, nnz_out) row.
    casadi_int numel_out = out_sp.size1() * out_sp.size2();
    bool dense_out = (static_cast<casadi_int>(idx.size()) == numel_out);
    // A dense column-vector output is already in stored form -> write straight to the output.
    std::string gathered = (dense_out && out_sp.size2() == 1) ? node_output : ("gnz_g_" + uniq);
    add_int_constant(add_node, "indices_" + uniq, idx);
    create_gather_node(add_node, flat, "indices_" + uniq, gathered, 1);

    if (dense_out) {
      if (out_sp.size2() != 1) {
        emit_colmajor_reshape(add_node, gathered, {out_sp.size1(), out_sp.size2()}, node_output,
                              "gnz_out_" + uniq);
      }
    } else {
      // SPARSE output: scatter the gathered values into a dense (1, numel_out) zero frame at the
      // output's dense column-major positions, then reshape to a DENSE (out_sp shape) value with
      // exact numerics (zero outside out_sp). The pattern is then restored to out_sp natively via
      // the seed recipe -- import re-propagates out_sp for free (no overlay).
      std::vector<casadi_int> ofind = out_sp.find();
      std::string zname = "gnz_zeros_" + uniq, fname = "gnz_ofind_" + uniq, scat = "gnz_scat_" + uniq;
      add_real_constant(add_node, zname, std::vector<double>(numel_out, 0.0), {1, numel_out});
      add_int_constant(add_node, fname, ofind, {1, static_cast<casadi_int>(ofind.size())});
      onnx::NodeProto* sc = add_node();
      sc->set_op_type("ScatterElements");
      sc->add_input(zname); sc->add_input(fname); sc->add_input(gathered);
      sc->add_output(scat);
      add_int_attribute(sc, "axis", 1);
      std::string dframe = "gnz_dense_" + uniq;
      emit_colmajor_reshape(add_node, scat, {out_sp.size1(), out_sp.size2()}, dframe,
                            "gnz_out_" + uniq);
      emit_sparsity_restore(add_node, dframe, Sparsity::dense(out_sp.size1(),
                            out_sp.size2()), out_sp, "gnzr_" + uniq, node_output);
    }
  }

  // Plant a sparsity seed so a dense ONNX value imports to exactly target_sp (see header).
  std::string Onnx::emit_sparsity_restore(AddNodeFn add_node, const std::string& value,
                                          const Sparsity& value_sp, const Sparsity& target_sp,
                                          const std::string& uniq,
                                          const std::string& final_output) {
    // Already imports to the right pattern -> nothing to do (0-node passthrough).
    if (value_sp == target_sp) return value;

    if (target_sp.is_dense()) {
      // value is non-dense here (a dense value_sp of this shape would equal target_sp above).
      // Add a dense real zero of target_sp's shape: dense + dense -> dense.
      std::string zc = "sr_dz_" + uniq;
      add_real_constant(add_node, zc,
                        std::vector<double>(target_sp.size1() * target_sp.size2(), 0.0),
                        {target_sp.size2(), target_sp.size1()});
      create_binary_node(add_node, "Add", value, zc, final_output);
      return final_output;
    }

    // Non-dense target. Classify the value pattern against the target:
    //  - value_sp subset of target_sp: the value only DROPS entries vs target; a single
    //    Add(value, zeros_target) suffices, since pattern(value) U target = target.
    //  - value_sp superset of target_sp (e.g. dense -> sparse narrowing project): the value only
    //    carries EXTRA entries; a single Mul(value, ones_target) suffices, since
    //    pattern(value) ^ target = target.
    //  - neither: the value both adds and drops entries, so mask with ones_target (intersection)
    //    then union with zeros_target -- exactly project(value, target).
    bool subset = value_sp.is_subset(target_sp);
    bool superset = target_sp.is_subset(value_sp);
    if (subset) {
      std::string zeros = "sr_zeros_" + uniq;
      add_sparse_constant(add_node, zeros, DM(target_sp, 0.0));
      create_binary_node(add_node, "Add", value, zeros, final_output);
      return final_output;
    }
    if (superset) {
      std::string ones = "sr_ones_" + uniq;
      add_sparse_constant(add_node, ones, DM(target_sp, 1.0));
      create_binary_node(add_node, "Mul", value, ones, final_output);
      return final_output;
    }
    // General: Mul(value, ones_target) then Add(_, zeros_target).
    std::string ones = "sr_ones_" + uniq;
    add_sparse_constant(add_node, ones, DM(target_sp, 1.0));
    std::string mul = "sr_mul_" + uniq;
    create_binary_node(add_node, "Mul", value, ones, mul);
    std::string zeros = "sr_zeros_" + uniq;
    add_sparse_constant(add_node, zeros, DM(target_sp, 0.0));
    create_binary_node(add_node, "Add", mul, zeros, final_output);
    return final_output;
  }

  // Assemble a block matrix pseudo-dense: Pad each DENSE block to (R,C) at its offset, then Sum.
  void Onnx::emit_blockdiag(AddNodeFn add_node, const std::vector<std::string>& names,
                            const std::vector<casadi_int>& row_off,
                            const std::vector<casadi_int>& col_off,
                            const std::vector<casadi_int>& br, const std::vector<casadi_int>& bc,
                            casadi_int R, casadi_int C, const std::string& output,
                            const std::string& uniq) {
    std::vector<std::string> padded;
    for (casadi_int b = 0; b < static_cast<casadi_int>(names.size()); ++b) {
      // Stored block is (bc x br); stored output is (C x R). ONNX Pad's `pads` is
      // [axis0_begin, axis1_begin, axis0_end, axis1_end] = [col_off, row_off, end_col, end_row].
      std::string pn = uniq + "_pp" + std::to_string(b);
      add_int_constant(add_node, pn, {col_off[b], row_off[b],
                                      C - col_off[b] - bc[b], R - row_off[b] - br[b]});
      std::string pd = uniq + "_pad" + std::to_string(b);
      onnx::NodeProto* pad = add_node();
      pad->set_op_type("Pad");
      pad->add_input(names[b]);
      pad->add_input(pn);  // mode defaults to "constant", value to 0
      pad->add_output(pd);
      padded.push_back(pd);
    }
    if (padded.size() == 1) {
      create_unary_node(add_node, "Identity", padded[0], output);
    } else {
      onnx::NodeProto* sum = add_node();
      sum->set_op_type("Sum");
      for (const auto& p : padded) sum->add_input(p);
      sum->add_output(output);
    }
  }

  // Callback-based process_operation (main implementation)
  bool Onnx::process_operation(
      AddNodeFn add_node,
      const Function& f,
      casadi_int op,
      casadi_int k,
      const std::vector<casadi_int>& i_vec,
      const std::vector<casadi_int>& o_vec,
      std::map<casadi_int, std::string>& work_to_onnx,
      const std::string& node_output) {

    onnx::NodeProto* node = nullptr;

    // Simple ops via the centralized lookup table
    const OpMapping* mapping = get_op_mapping(op);
    if (mapping) {
      if (mapping->arity == 1 && i_vec.size() >= 1 && o_vec.size() == 1) {
        create_unary_node(add_node, mapping->onnx_name, work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      } else if (mapping->arity == 2 && i_vec.size() >= 2 && o_vec.size() == 1) {
        create_binary_node(add_node, mapping->onnx_name, work_to_onnx[i_vec[0]],
                           work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }
    }

    // Everything else needs attributes, multiple nodes, or shape logic.
    switch (op) {
      case OP_INPUT: {
        std::string input_name = onnx_input_name(f, i_vec[0]);

        // A sub-slice of a larger input has offset/numel narrower than the whole input.
        MX mx_input = f.instruction_MX(k);
        Dict info = mx_input.info();
        casadi_int offset = info["offset"];
        casadi_int input_numel = f.numel_in(i_vec[0]);
        casadi_int output_numel = mx_input.numel();

        if (output_numel < input_numel) {
          // Extract output_numel consecutive column-major elements starting at `offset`.
          // Transpose-rep: flatten the stored input to a 1 x numel row, Gather along axis 1.
          std::string uniq = std::to_string(k);
          std::string flat = input_name;
          if (f.size2_in(i_vec[0]) != 1) {  // not already stored as a row
            flat = "in_flat_" + uniq;
            emit_colmajor_reshape(add_node, input_name, {input_numel, 1}, flat, "in_flat_" + uniq);
          }
          std::vector<casadi_int> idx;
          for (casadi_int t = 0; t < output_numel; ++t) idx.push_back(offset + t);
          std::string index_name = "input_idx_" + uniq;
          add_int_constant(add_node, index_name, idx);
          bool out_is_col = (mx_input.size2() == 1);
          std::string gathered = out_is_col ? node_output : ("in_g_" + uniq);
          create_gather_node(add_node, flat, index_name, gathered, 1);
          if (!out_is_col) {
            emit_colmajor_reshape(add_node, gathered, {mx_input.size1(), mx_input.size2()},
                                  node_output, "in_out_" + uniq);
          }
        } else {
          // Full input access (output size matches input size): no rename node -- downstream nodes
          // reference the graph input tensor directly. If this input is returned DIRECTLY as an
          // output, the OP_OUTPUT path emits the single Identity needed to give it the output name
          // (a graph input cannot be renamed in place).
          work_to_onnx[o_vec[0]] = input_name;
          return true;
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_OUTPUT: {
        // Connect the work value to the output name; o_vec[0] is the output slot, not a work index
        create_unary_node(add_node, "Identity", work_to_onnx[i_vec[0]], onnx_output_name(f, o_vec[0]));
        return true;
      }

      case OP_CONST: {
        DM dm_const = static_cast<DM>(f.instruction_MX(k));
        if (dm_const.is_dense()) {
          // Transpose-rep: the column-major bytes go in as-is, shape declared REVERSED (c,r).
          std::vector<double> data(dm_const->begin(), dm_const->end());
          add_real_constant(add_node, node_output, data, {dm_const.size2(), dm_const.size1()});
        } else {
          // Store the sparsity pattern faithfully as an ONNX sparse constant (COO)
          add_sparse_constant(add_node, node_output, dm_const);
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_MTIMES: {
        // Fused multiply-add input[1]*input[2] + input[0] -> a single Gemm
        create_gemm_node(add_node, work_to_onnx[i_vec[1]], work_to_onnx[i_vec[2]],
                         work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_SQ: {
        // x^2 = Mul(x, x)
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[0]],
                          work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_TWICE: {
        // 2*x
        std::string const_name = "const_2_" + std::to_string(k);
        add_real_constant(add_node, const_name, {2.0});
        create_binary_node(add_node, "Mul", const_name, work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_DETERMINANT:
        create_unary_node(add_node, "Det", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_LOGSUMEXP:
        create_unary_node(add_node, "ReduceLogSumExp", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_LOG1P: {
        // log(1 + x)
        std::string one = "one_" + std::to_string(k), s = "log1p_" + std::to_string(k);
        add_real_constant(add_node, one, {1.0});
        create_binary_node(add_node, "Add", one, work_to_onnx[i_vec[0]], s);
        create_unary_node(add_node, "Log", s, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_EXPM1: {
        // exp(x) - 1
        std::string e = "exp_" + std::to_string(k), one = "one_" + std::to_string(k);
        create_unary_node(add_node, "Exp", work_to_onnx[i_vec[0]], e);
        add_real_constant(add_node, one, {1.0});
        create_binary_node(add_node, "Sub", e, one, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HYPOT: {
        // sqrt(x^2 + y^2)
        std::string x2 = "hx_" + std::to_string(k), y2 = "hy_" + std::to_string(k);
        std::string s = "hs_" + std::to_string(k);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[0]], x2);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[1]], work_to_onnx[i_vec[1]], y2);
        create_binary_node(add_node, "Add", x2, y2, s);
        create_unary_node(add_node, "Sqrt", s, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Pure pass-throughs (drop CasADi-internal markers)
      case OP_ASSIGN:
      case OP_LIFT:
      case OP_ASSERTION:   // attachAssert: output = dep(0); the condition (dep(1)) is a side guard
      case OP_MONITOR:     // monitor/printme-style debug passthrough: output = dep(0)
        create_unary_node(add_node, "Identity", work_to_onnx[i_vec[0]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_SOLVE: {
        // ONNX has no native linear solver. tr selects A'*x = b and is reported for diagnostics.
        MX mx_solve = f.instruction_MX(k);
        Dict info = mx_solve.info();
        bool tr = info["tr"];
        casadi_error("ONNX export: OP_SOLVE (linear solver) is not supported. "
                    "ONNX does not provide native linear algebra solvers. "
                    "Consider using explicit matrix operations or iterative methods. "
                    "Transpose flag was: " + std::string(tr ? "true" : "false"));
        return false;
      }

      case OP_NORMINF: {
        // Infinity norm: max(abs(x))
        std::string abs_result = "abs_" + std::to_string(k);
        create_unary_node(add_node, "Abs", work_to_onnx[i_vec[0]], abs_result);
        create_unary_node(add_node, "ReduceMax", abs_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Comparisons output a bool tensor in ONNX; cast to double for the double-typed graph
      case OP_LT: {
        std::string b = "cmp_" + std::to_string(k);
        create_binary_node(add_node, "Less", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], b);
        create_cast_node(add_node, b, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_LE: {
        std::string b = "cmp_" + std::to_string(k);
        create_binary_node(add_node, "LessOrEqual", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], b);
        create_cast_node(add_node, b, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_EQ: {
        std::string b = "cmp_" + std::to_string(k);
        create_binary_node(add_node, "Equal", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], b);
        create_cast_node(add_node, b, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Logic ops need bool operands: cast 0/1 doubles to bool, apply, cast back
      case OP_NOT: {
        std::string b = "bin_" + std::to_string(k);
        create_cast_node(add_node, work_to_onnx[i_vec[0]], b, onnx::TensorProto::BOOL);
        std::string r = "bres_" + std::to_string(k);
        create_unary_node(add_node, "Not", b, r);
        create_cast_node(add_node, r, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_AND:
      case OP_OR: {
        std::string b0 = "bin0_" + std::to_string(k), b1 = "bin1_" + std::to_string(k);
        create_cast_node(add_node, work_to_onnx[i_vec[0]], b0, onnx::TensorProto::BOOL);
        create_cast_node(add_node, work_to_onnx[i_vec[1]], b1, onnx::TensorProto::BOOL);
        std::string r = "bres_" + std::to_string(k);
        create_binary_node(add_node, op == OP_AND ? "And" : "Or", b0, b1, r);
        create_cast_node(add_node, r, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_FMIN:
        create_binary_node(add_node, "Min", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_FMOD: {
        node = add_node();
        node->set_op_type("Mod");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(work_to_onnx[i_vec[1]]);
        node->add_output(node_output);
        add_int_attribute(node, "fmod", 1);  // fmod=1 required for floating point types
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_COPYSIGN: {
        // copysign(x, y) = sign(y) * abs(x)
        std::string sign_result = "sign_" + std::to_string(k);
        std::string abs_result = "abs_" + std::to_string(k);
        create_unary_node(add_node, "Sign", work_to_onnx[i_vec[1]], sign_result);
        create_unary_node(add_node, "Abs", work_to_onnx[i_vec[0]], abs_result);
        create_binary_node(add_node, "Mul", sign_result, abs_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_FMAX:
        create_binary_node(add_node, "Max", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;

      case OP_NE: {
        // Not equal: Not(Equal(...)), then cast the bool result to double
        std::string equal_result = "equal_" + std::to_string(k);
        create_binary_node(add_node, "Equal", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], equal_result);
        std::string r = "bres_" + std::to_string(k);
        create_unary_node(add_node, "Not", equal_result, r);
        create_cast_node(add_node, r, node_output, real_type());
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_IF_ELSE_ZERO: {
        // if (condition != 0) value else 0  =>  Where(bool(condition), value, 0)
        std::string cond = "cond_" + std::to_string(k);
        create_cast_node(add_node, work_to_onnx[i_vec[0]], cond, onnx::TensorProto::BOOL);
        std::string zero_name = "const_0_" + std::to_string(k);
        add_real_constant(add_node, zero_name, {0.0});
        create_where_node(add_node, cond, work_to_onnx[i_vec[1]], zero_name, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_DOT: {
        // dot(a, b) = ReduceSum(Mul(a, b))
        std::string mul_result = "mul_" + std::to_string(k);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[0]], work_to_onnx[i_vec[1]], mul_result);
        create_unary_node(add_node, "ReduceSum", mul_result, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_BILIN: {
        // Bilinear form input[1]' * input[0] * input[2]. Gemm(transA) avoids a separate
        // Transpose (which ORT would fuse into a float64-less FusedMatMul).
        std::string xa = "bilin_" + std::to_string(k);
        create_gemm_node(add_node, work_to_onnx[i_vec[1]], work_to_onnx[i_vec[0]], "", xa, true);
        // Transpose-rep: stored out = y_stored * xa_stored, so emit MatMul(y, xa)
        create_binary_node(add_node, "MatMul", work_to_onnx[i_vec[2]], xa, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_RANK1: {
        // Rank-1 update input[0] + input[1]*input[2]*input[3]'. Gemm(transB) forms x*y'
        // without a separate Transpose; input[1] is a scalar alpha.
        std::string xyt = "rank1_" + std::to_string(k);
        create_gemm_node(add_node, work_to_onnx[i_vec[2]], work_to_onnx[i_vec[3]], "", xyt, false, true);
        std::string scaled = "rank1s_" + std::to_string(k);
        create_binary_node(add_node, "Mul", work_to_onnx[i_vec[1]], xyt, scaled);
        create_binary_node(add_node, "Add", work_to_onnx[i_vec[0]], scaled, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_EINSTEIN: {
        // Tensor contraction C_c += A_a * B_b. deps [C, A, B], operands column-vec flattened.
        // Reshape each operand (column-major) to its logical dims and contract with Einsum
        // using the labels in natural order; reshape the result back to a column vector.
        MX mx_e = f.instruction_MX(k);
        Dict info = mx_e.info();
        std::vector<casadi_int> la = info["a"], lb = info["b"], lc = info["c"];
        std::vector<casadi_int> da = info["dim_a"], db = info["dim_b"], dc = info["dim_c"];
        casadi_assert(da.size() <= 2 && db.size() <= 2 && dc.size() <= 2,
          "ONNX export: einstein with >2-index operands is not supported.");

        // Assign a letter to each distinct label (natural order)
        std::map<casadi_int, char> letter;
        char nxt = 'a';
        for (const std::vector<casadi_int>& labs : {la, lb, lc}) {
          for (casadi_int L : labs) if (!letter.count(L)) letter[L] = nxt++;
        }
        std::string sa, sb, sc;
        for (casadi_int L : la) sa += letter[L];
        for (casadi_int L : lb) sb += letter[L];
        for (casadi_int L : lc) sc += letter[L];
        // Transpose-rep: emit_colmajor_reshape stores each operand with axes REVERSED (logical d0xd1
        // -> ONNX [d1,d0]). Reverse each subscript so the equation labels the stored axes correctly;
        // the reversed-axis output then flattens column-major back into C's vector.
        std::reverse(sa.begin(), sa.end());
        std::reverse(sb.begin(), sb.end());
        std::reverse(sc.begin(), sc.end());

        std::string a_re = "ein_a_" + std::to_string(k), b_re = "ein_b_" + std::to_string(k);
        emit_colmajor_reshape(add_node, work_to_onnx[i_vec[1]], da, a_re, a_re);
        emit_colmajor_reshape(add_node, work_to_onnx[i_vec[2]], db, b_re, b_re);

        std::string ein_out = "ein_o_" + std::to_string(k);
        node = add_node();
        node->set_op_type("Einsum");
        node->add_input(a_re);
        node->add_input(b_re);
        node->add_output(ein_out);
        onnx::AttributeProto* eq = node->add_attribute();
        eq->set_name("equation");
        eq->set_type(onnx::AttributeProto::STRING);
        eq->set_s(sa + "," + sb + "->" + sc);

        // Flatten the contraction result back to a column vector and add the C accumulator
        casadi_int prod_c = 1;
        for (casadi_int d : dc) prod_c *= d;
        std::string ein_vec = "ein_v_" + std::to_string(k);
        emit_colmajor_reshape(add_node, ein_out, {prod_c, 1}, ein_vec, ein_vec);
        create_binary_node(add_node, "Add", work_to_onnx[i_vec[0]], ein_vec, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_KRON: {
        // Kronecker product kron(A, B), A is ra x ca, B is rb x cb, result (ra*rb) x (ca*cb).
        // kron commutes with transpose: kron(A,B)^T == kron(A^T,B^T). The stored tensors are the
        // transposes (At=[ca,ra], Bt=[cb,rb]) and the stored output is the result's transpose, so in
        // stored space kron maps to kron of the stored operands with NO extra transposes. Emit
        //   A4 = Reshape(At, [ca,1,ra,1]);  B4 = Reshape(Bt, [1,cb,1,rb])
        //   P  = Mul(A4, B4)  -> broadcast [ca,cb,ra,rb]
        //   R  = Reshape(P, [ca*cb, ra*rb])  (= stored result)
        // These are plain row-major Reshapes (insert/collapse singleton axes), so emit the shape
        // constants directly (no axis reversal -> do NOT use emit_colmajor_reshape). The node group
        // is tagged with name "kron<k>" so the importer reconstructs kron(A_imp,B_imp) directly
        // (the 4-D intermediates cannot be represented as 2-D MX).
        MX mx_kron = f.instruction_MX(k);
        casadi_int ra = mx_kron.dep(0).size1(), ca = mx_kron.dep(0).size2();
        casadi_int rb = mx_kron.dep(1).size1(), cb = mx_kron.dep(1).size2();
        std::string uniq = std::to_string(k);
        std::string tag = "kron" + uniq;

        // ONNX node names must be unique; give each kron node a distinct name sharing the "kron<k>_"
        // prefix and a role suffix the importer dispatches on.
        auto emit_reshape = [&](const std::string& data, const std::vector<casadi_int>& shape,
                                const std::string& out, const std::string& role) {
          add_int_constant(add_node, out + "_s", shape);
          onnx::NodeProto* rn = create_unary_node(add_node, "Reshape", data, out);
          rn->add_input(out + "_s");
          rn->set_name(tag + "_" + role);
        };

        std::string a4 = tag + "_A4", b4 = tag + "_B4", p = tag + "_P";
        emit_reshape(work_to_onnx[i_vec[0]], {ca, 1, ra, 1}, a4, "A4");
        emit_reshape(work_to_onnx[i_vec[1]], {1, cb, 1, rb}, b4, "B4");
        onnx::NodeProto* mul = create_binary_node(add_node, "Mul", a4, b4, p);
        mul->set_name(tag + "_P");
        emit_reshape(p, {ca * cb, ra * rb}, node_output, "R");

        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Tensor operations
      case OP_RESHAPE: {
        // ONNX Reshape is pseudo-dense; the importer densifies the operand (it must, to keep flat
        // indices aligned), so a non-dense reshape result is restored to its CasADi pattern natively
        // from the seed (no overlay).
        auto sp = f.instruction_MX(k).sparsity();
        if (sp.is_dense()) {
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[0]], {sp.size1(), sp.size2()},
                                node_output, "rs_" + std::to_string(k));
        } else {
          std::string rs = "rs_d_" + std::to_string(k);
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[0]], {sp.size1(), sp.size2()},
                                rs, "rs_" + std::to_string(k));
          std::string r = emit_sparsity_restore(add_node, rs,
                                                Sparsity::dense(sp.size1(), sp.size2()), sp,
                                                "rsr_" + std::to_string(k), node_output);
          work_to_onnx[o_vec[0]] = r;
          return true;
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZCAT:
      case OP_VERTCAT: {
        // Transpose-rep: CasADi columns/rows are ONNX rows/columns, so axes swap:
        // horzcat (CasADi axis 1) -> ONNX axis 0, vertcat (CasADi axis 0) -> ONNX axis 1.
        node = add_node();
        node->set_op_type("Concat");
        for (casadi_int idx : i_vec) node->add_input(work_to_onnx[idx]);
        node->add_output(node_output);
        add_int_attribute(node, "axis", op == OP_HORZCAT ? 0 : 1);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_DIAGCAT: {
        // Block-diagonal, pseudo-dense: Pad each DENSE block to the output shape at its diagonal
        // offset, then Sum -- no flatten/gather (a triangular block is just a dense block with a
        // structural zero; the output pattern is restored from the seed).
        MX mx_dc = f.instruction_MX(k);
        std::string uniq = std::to_string(k);
        std::vector<std::string> names;
        std::vector<casadi_int> row_off, col_off, brs, bcs;
        casadi_int ro = 0, co = 0;
        for (casadi_int j = 0; j < static_cast<casadi_int>(i_vec.size()); ++j) {
          casadi_int br = mx_dc.dep(j).size1(), bc = mx_dc.dep(j).size2();
          names.push_back(work_to_onnx[i_vec[j]]);
          row_off.push_back(ro); col_off.push_back(co); brs.push_back(br); bcs.push_back(bc);
          ro += br; co += bc;  // diagonal placement
        }
        Sparsity dc_sp = mx_dc.sparsity();
        if (dc_sp.is_dense()) {
          emit_blockdiag(add_node, names, row_off, col_off, brs, bcs,
                         mx_dc.size1(), mx_dc.size2(), node_output, "dc" + uniq);
        } else {
          // The Pad+Sum assembly imports dense; restore the block-diagonal pattern from the seed.
          std::string tmp = "dc_pre_" + uniq;
          emit_blockdiag(add_node, names, row_off, col_off, brs, bcs,
                         mx_dc.size1(), mx_dc.size2(), tmp, "dc" + uniq);
          std::string r = emit_sparsity_restore(add_node, tmp,
                              Sparsity::dense(mx_dc.size1(), mx_dc.size2()), dc_sp, "dcr_" + uniq,
                              node_output);
          work_to_onnx[o_vec[0]] = r;
          return true;
        }
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZSPLIT:
      case OP_VERTSPLIT: {
        // Split the pseudo-dense tensor by COLUMN (horzsplit) / ROW (vertsplit) counts. The MX
        // "offset" is a NONZERO offset, so dividing by a stride only yields column/row counts for a
        // DENSE input; for a sparse input it is wrong (nonzeros are not uniform per column). The
        // segment dimensions are the authoritative sizes -- read them straight off the split
        // primitive's output Function. The importer splits the sparse seed natively, so each
        // segment recovers its own pattern (no overlay).
        MX mx_split = f.instruction_MX(k);
        bool horz = (op == OP_HORZSPLIT);
        Function split_out = mx_split.info()["output"];

        std::vector<casadi_int> split_sizes;
        for (casadi_int j = 0; j < split_out.n_out(); ++j) {
          split_sizes.push_back(horz ? split_out.size2_out(j) : split_out.size1_out(j));
        }

        std::string split_sizes_name = "split_sizes_" + std::to_string(k);
        add_int_constant(add_node, split_sizes_name, split_sizes);

        node = add_node();
        node->set_op_type("Split");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(split_sizes_name);
        // Transpose-rep: horzsplit (CasADi axis 1) -> ONNX axis 0, vertsplit -> ONNX axis 1
        add_int_attribute(node, "axis", horz ? 0 : 1);

        for (casadi_int j = 0; j < o_vec.size(); ++j) {
          std::string output_name = "n" + std::to_string(k) + "_out" + std::to_string(j);
          node->add_output(output_name);
          work_to_onnx[o_vec[j]] = output_name;
        }
        return true;
      }

      case OP_HORZREPMAT: {
        // repmat(x, 1, n) -> ONNX Tile with repeats [1, n]
        MX mx_repmat = f.instruction_MX(k);
        casadi_int n = mx_repmat.size2() / mx_repmat.dep(0).size2();

        std::string repeats_name = "repeats_" + std::to_string(k);
        add_int_constant(add_node, repeats_name, {n, 1});  // transpose-rep: tile over rows

        node = add_node();
        node->set_op_type("Tile");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(repeats_name);
        node->add_output(node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_HORZREPSUM: {
        // repsum(x, 1, n): inverse of repmat - split [rows, cols*n] into n parts and sum
        MX mx_repsum = f.instruction_MX(k);
        casadi_int output_cols = mx_repsum.size2();
        casadi_int n = mx_repsum.dep(0).size2() / output_cols;

        std::string split_sizes_name = "split_sizes_" + std::to_string(k);
        add_int_constant(add_node, split_sizes_name, std::vector<casadi_int>(n, output_cols));

        node = add_node();
        node->set_op_type("Split");
        node->add_input(work_to_onnx[i_vec[0]]);
        node->add_input(split_sizes_name);
        add_int_attribute(node, "axis", 0);  // transpose-rep: repmat tiled over rows

        // Sum all n parts in a single variadic Sum node
        onnx::NodeProto* sum_node = add_node();
        sum_node->set_op_type("Sum");
        for (casadi_int i = 0; i < n; ++i) {
          std::string out_name = "split_" + std::to_string(k) + "_" + std::to_string(i);
          node->add_output(out_name);
          sum_node->add_input(out_name);
        }
        sum_node->add_output(node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_GETNONZEROS: {
        // x.nz[indices], CasADi column-major. Transpose-rep: row-major-flatten the stored input
        // (its row-major buffer IS x's column-major order) to a 1 x numel ROW vector, Gather along
        // axis 1, reshape back. The nz indices address x's NONZERO array; map them to dense
        // column-major flat positions with find() (identity when x is dense).
        MX mx_getnonzeros = f.instruction_MX(k);
        Dict info = mx_getnonzeros.info();
        std::string data = work_to_onnx[i_vec[0]];
        std::string uniq = std::to_string(k);

        // Collect the explicit column-major nonzero indices
        std::vector<casadi_int> idx;
        if (info.count("nz")) {
          idx = info["nz"];
        } else if (info.count("slice")) {
          Dict s = info["slice"];
          casadi_int start = s["start"], stop = s["stop"], step = s["step"];
          for (casadi_int i = start; i < stop; i += step) idx.push_back(i);
        } else if (info.count("inner")) {
          Dict inner_info = info["inner"], outer_info = info["outer"];
          casadi_int is = inner_info["start"], ip = inner_info["stop"], ist = inner_info["step"];
          casadi_int os = outer_info["start"], op2 = outer_info["stop"], ost = outer_info["step"];
          for (casadi_int o = os; o < op2; o += ost) {
            for (casadi_int i = is; i < ip; i += ist) idx.push_back(o + i);
          }
        } else {
          return false;  // unknown variant
        }

        // Recognize a structured slice of the (pretend-dense) input and emit a clean ONNX Slice.
        //  - DENSE input: ANY regular grid A[r0:r1:rs, c0:c1:cs] (single row/col, range, block,
        //    strided) decoded geometrically from the column-major positions -> Slice on both axes.
        //  - SPARSE input: a single full row/col matched against the ACTUAL nonzeros, which
        //    round-trips WITH sparsity (the importer rebuilds A[r,:]/A[:,c] on the sparse value).
        // Transpose-rep: CasADi rows -> ONNX axis 1, CasADi cols -> ONNX axis 0.
        {
          const Sparsity& sp_in = mx_getnonzeros.dep(0).sparsity();
          casadi_int m_in = sp_in.size1(), n_in = sp_in.size2();
          if (sp_in.is_dense() && !idx.empty()) {
            // `idx` are dense column-major positions. Column-major layout means the first run of
            // equal-column entries is the row set; its repeat stride is the column stride.
            casadi_int N = static_cast<casadi_int>(idx.size());
            casadi_int c0 = idx[0] / m_in, r0 = idx[0] % m_in;
            casadi_int k = 1;
            while (k < N && idx[k] / m_in == c0) ++k;        // rows in the first column
            bool grid = (N % k == 0);
            casadi_int nc = grid ? N / k : 0;
            casadi_int rstep = (k > 1) ? (idx[1] % m_in) - r0 : 1;
            casadi_int cstep = (grid && nc > 1) ? idx[k] / m_in - c0 : 1;
            if (grid && rstep > 0 && cstep > 0) {
              for (casadi_int cj = 0; cj < nc && grid; ++cj)
                for (casadi_int ri = 0; ri < k && grid; ++ri)
                  if (idx[cj * k + ri] != (c0 + cj * cstep) * m_in + (r0 + ri * rstep)) grid = false;
              if (grid) {
                create_slice_node(add_node, "gnzs_" + uniq, data,
                    {c0, r0}, {c0 + (nc - 1) * cstep + 1, r0 + (k - 1) * rstep + 1},
                    {0, 1}, {cstep, rstep}, node_output);
                work_to_onnx[o_vec[0]] = node_output;
                return true;
              }
            }
          } else if (!idx.empty()) {
            std::vector<casadi_int> irow = sp_in.get_row(), icol = sp_in.get_col();
            casadi_int nnz_in = static_cast<casadi_int>(irow.size());
            if (mx_getnonzeros.size1() == 1 && mx_getnonzeros.size2() == n_in) {
              casadi_int r = irow[idx[0]];
              std::vector<casadi_int> exp;
              for (casadi_int p = 0; p < nnz_in; ++p) if (irow[p] == r) exp.push_back(p);
              if (exp == idx) {
                create_slice_node(add_node, "gnzs_" + uniq, data, {r}, {r + 1}, {1}, {}, node_output);
                work_to_onnx[o_vec[0]] = node_output;
                return true;
              }
            }
            if (mx_getnonzeros.size1() == m_in && mx_getnonzeros.size2() == 1) {
              casadi_int c = icol[idx[0]];
              std::vector<casadi_int> exp;
              for (casadi_int p = 0; p < nnz_in; ++p) if (icol[p] == c) exp.push_back(p);
              if (exp == idx) {
                create_slice_node(add_node, "gnzs_" + uniq, data, {c}, {c + 1}, {0}, {}, node_output);
                work_to_onnx[o_vec[0]] = node_output;
                return true;
              }
            }
          }
        }

        // General fallback: route the explicit nonzero map `idx` (output-nz -> input-nz, -1=fill)
        // through the shared remap helper -- the poorest-form primitive shared with OP_PROJECT.
        emit_nonzero_remap(add_node, data, mx_getnonzeros.dep(0).sparsity(),
                           mx_getnonzeros.sparsity(), idx, uniq, node_output);
        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_GETNONZEROS_PARAM: {
        // y = x.nz[p] with p a RUNTIME index vector (dep(1)), p[k] addressing x's column-major
        // NONZERO array. Mirrors the constant OP_GETNONZEROS but the indices flow as a value.
        // Transpose-rep: flatten the stored x to a 1 x numel ROW (its row-major buffer IS x's
        // column-major order); p (shape np x 1) is stored as a 1 x np ROW. Map nz->dense
        // col-major positions with find() (a compile-time GatherElements; identity if x dense),
        // then GatherElements the flat row by those positions -> a 1 x np ROW == stored output.
        // GatherElements(data[1,N], idx[1,M], axis=1) returns shape == idx [1,M] (no rank blowup).
        MX mx_gnzp = f.instruction_MX(k);
        const Sparsity& sp_in = mx_gnzp.dep(0).sparsity();
        casadi_int numel_in = sp_in.size1() * sp_in.size2();
        std::string uniq = std::to_string(k);

        // x's column-major flat as a 1 x numel ROW.
        std::string flat = work_to_onnx[i_vec[0]];
        if (sp_in.size2() != 1) {
          flat = "gnzp_flat_" + uniq;
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[0]], {numel_in, 1}, flat,
                                "gnzp_flat_" + uniq);
        }

        // Runtime indices p (stored as a 1 x np row of real values) -> INT64 nz-indices.
        std::string nzidx = "gnzp_pi_" + uniq;
        create_cast_node(add_node, work_to_onnx[i_vec[1]], nzidx, onnx::TensorProto::INT64);

        // Map nz-indices -> dense column-major positions. Identity if x dense; otherwise gather the
        // constant find() table by the runtime nz-indices (GatherElements out shape == idx [1,np]).
        std::string dpos = nzidx;
        if (!sp_in.is_dense()) {
          std::vector<casadi_int> loc = sp_in.find();  // nz-index -> dense col-major position
          std::string loc_c = "gnzp_loc_" + uniq;
          add_int_constant(add_node, loc_c, loc, {1, static_cast<casadi_int>(loc.size())});
          dpos = "gnzp_dpos_" + uniq;
          onnx::NodeProto* ge = add_node();
          ge->set_op_type("GatherElements");
          ge->add_input(loc_c); ge->add_input(nzidx); ge->add_output(dpos);
          add_int_attribute(ge, "axis", 1);
        }

        // Gather x's flat row by the dense positions: GatherElements(flat[1,numel], dpos[1,np]) ->
        // [1,np] == the stored output (a column vector shaped like p stored transposed).
        onnx::NodeProto* ge = add_node();
        ge->set_op_type("GatherElements");
        ge->add_input(flat); ge->add_input(dpos); ge->add_output(node_output);
        add_int_attribute(ge, "axis", 1);

        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      case OP_SETNONZEROS_PARAM:
      case OP_ADDNONZEROS_PARAM: {
        // y = x with y.nz[p] {=,+=} v, p a RUNTIME index vector (dep(2)) into x's column-major
        // NONZERO array; v = dep(1) the values, x = dep(0) the base. Pseudo-dense, on the flattened
        // x: ScatterElements(flat[1,numel], dpos[1,np], updates[1,np], axis=1, reduction). ADD uses
        // reduction="add" (adds onto x AND accumulates DUPLICATE p, as the getnonzeros-adjoint
        // requires); SET uses "none" (overwrite, last wins). Mirrors OP_GETNONZEROS_PARAM's
        // flatten + INT64 cast + find() nz->dense map. ADDNONZEROS_PARAM has no direct MX builder --
        // it is the reverse-mode adjoint of a parametric getnonzeros; the importer rebuilds it so.
        bool add = (op == OP_ADDNONZEROS_PARAM);
        MX mx_snzp = f.instruction_MX(k);
        const Sparsity& sp_in = mx_snzp.dep(0).sparsity();
        casadi_int numel_in = sp_in.size1() * sp_in.size2();
        casadi_int np_idx = mx_snzp.dep(2).nnz();
        std::string uniq = std::to_string(k);

        // x's column-major flat as a 1 x numel ROW.
        std::string flat = work_to_onnx[i_vec[0]];
        if (sp_in.size2() != 1) {
          flat = "snzp_flat_" + uniq;
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[0]], {numel_in, 1}, flat,
                                "snzp_flat_" + uniq);
        }

        // Runtime indices p -> INT64 nz-indices, then map nz -> dense col-major positions (identity
        // if x dense; else gather the constant find() table -- GatherElements out shape == idx).
        std::string nzidx = "snzp_pi_" + uniq;
        create_cast_node(add_node, work_to_onnx[i_vec[2]], nzidx, onnx::TensorProto::INT64);
        std::string dpos = nzidx;
        if (!sp_in.is_dense()) {
          std::vector<casadi_int> loc = sp_in.find();
          std::string loc_c = "snzp_loc_" + uniq;
          add_int_constant(add_node, loc_c, loc, {1, static_cast<casadi_int>(loc.size())});
          dpos = "snzp_dpos_" + uniq;
          onnx::NodeProto* ge = add_node();
          ge->set_op_type("GatherElements");
          ge->add_input(loc_c); ge->add_input(nzidx); ge->add_output(dpos);
          add_int_attribute(ge, "axis", 1);
        }

        // Values v as a 1 x np ROW (its column-major flat; ORT materializes any structural zeros).
        std::string updates = work_to_onnx[i_vec[1]];
        if (mx_snzp.dep(1).size2() != 1) {
          updates = "snzp_v_" + uniq;
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[1]], {np_idx, 1}, updates, "snzp_v_" + uniq);
        }

        // ScatterElements onto the flat row.
        std::string scat = "snzp_scat_" + uniq;
        onnx::NodeProto* se = add_node();
        se->set_op_type("ScatterElements");
        se->add_input(flat); se->add_input(dpos); se->add_input(updates); se->add_output(scat);
        add_int_attribute(se, "axis", 1);
        if (add) {
          onnx::AttributeProto* red = se->add_attribute();
          red->set_name("reduction");
          red->set_type(onnx::AttributeProto::STRING);
          red->set_s("add");
        }

        // Reshape the modified flat back to x's stored shape, then restore the output pattern.
        std::string back = "snzp_back_" + uniq;
        emit_colmajor_reshape(add_node, scat, {mx_snzp.size1(), mx_snzp.size2()}, back, "snzp_back_" + uniq);
        std::string r = emit_sparsity_restore(add_node, back,
                            Sparsity::dense(mx_snzp.size1(), mx_snzp.size2()),
                            mx_snzp.sparsity(), "snzpr_" + uniq, node_output);
        work_to_onnx[o_vec[0]] = r;
        return true;
      }

      case OP_PROJECT: {
        // project(A, sp_out) / densify keeps the matrix shape and only ZEROS the dense cells not in
        // sp_out (it never permutes data). Numerically that is project(A, sp_out); pseudo-dense, the
        // emit_sparsity_restore seed reproduces both the numerics and the exact import pattern
        // sp_out natively (no overlay): densify -> dense Add-zero; widening -> sparse Add-zeros;
        // narrowing -> sparse Mul-ones then Add-zeros.
        MX mx_proj = f.instruction_MX(k);
        Sparsity sp_out = mx_proj.sparsity();
        std::string r = emit_sparsity_restore(add_node, work_to_onnx[i_vec[0]],
                                              mx_proj.dep(0).sparsity(), sp_out,
                                              std::to_string(k), node_output);
        work_to_onnx[o_vec[0]] = r;
        return true;
      }

      case OP_SETNONZEROS:
      case OP_ADDNONZEROS: {
        // y = x with y.nz[idx] = z (SETNONZEROS) or y.nz[idx] += z (ADDNONZEROS). Emit ONNX
        // ScatterND: scatter z's values into the 2-D pseudo-dense x at the (row,col) COORDINATES of
        // the target cells -- the canonical subtensor write, no flatten/reshape. idx address x's
        // NONZERO array -> map to dense column-major positions with find() (identity if dense).
        // ADDNONZEROS uses reduction="add" (adds onto x AND accumulates DUPLICATE idx, as required by
        // the adjoint of a getnonzeros with repeated indices); SETNONZEROS overwrites (reduction
        // "none"). ADDNONZEROS cannot be built as an MX directly -- it arises from reverse-mode AD
        // (jtimes) of a getnonzeros; the importer reconstructs it that way.
        MX mx_set = f.instruction_MX(k);
        Dict info = mx_set.info();
        std::string uniq = std::to_string(k);
        bool add = info["add"];

        std::vector<casadi_int> idx;
        if (info.count("nz")) {
          idx = info["nz"];
        } else if (info.count("slice")) {
          Dict s = info["slice"];
          casadi_int start = s["start"], stop = s["stop"], step = s["step"];
          for (casadi_int i = start; i < stop; i += step) idx.push_back(i);
        } else if (info.count("inner")) {
          Dict ii = info["inner"], oo = info["outer"];
          casadi_int is = ii["start"], ip = ii["stop"], ist = ii["step"];
          casadi_int os = oo["start"], op2 = oo["stop"], ost = oo["step"];
          for (casadi_int o = os; o < op2; o += ost) {
            for (casadi_int i = is; i < ip; i += ist) idx.push_back(o + i);
          }
        } else {
          casadi_error("ONNX export: unsupported setnonzeros pattern.");
        }
        std::vector<casadi_int> loc = mx_set.dep(0).sparsity().find();
        for (casadi_int& v : idx) v = loc[v];  // dense column-major positions in x

        // ScatterND coordinates [nnz, 2]: x is CasADi (nrow,ncol) stored as x^T (ncol,nrow), so a
        // dense position p -> CasADi cell (p%nrow, p/nrow) -> ONNX coord (p/nrow, p%nrow).
        casadi_int nrow = mx_set.dep(0).size1();
        std::vector<casadi_int> coords;
        coords.reserve(2 * idx.size());
        for (casadi_int p : idx) { coords.push_back(p / nrow); coords.push_back(p % nrow); }
        std::string ind = "snz_idx_" + uniq;
        add_int_constant(add_node, ind, coords, {static_cast<casadi_int>(idx.size()), 2});

        // updates for ScatterND must be a 1-D [nnz] tensor. Build z's nonzeros while still 2-D (CasADi
        // has no true 1-D, so a 2-D Gather round-trips cleanly), then reshape to 1-D last.
        casadi_int numel_z = mx_set.dep(1).numel();
        std::string vals = work_to_onnx[i_vec[1]];
        if (mx_set.dep(1).size2() != 1) {  // flatten to a (1,numel_z) row (a column on import)
          vals = "snz_zr_" + uniq;
          emit_colmajor_reshape(add_node, work_to_onnx[i_vec[1]], {numel_z, 1}, vals, "snz_zr_" + uniq);
        }
        // The write positions are z's NONZEROS, so a SPARSE z contributes only its nonzeros.
        if (mx_set.dep(1).nnz() != numel_z) {
          std::vector<casadi_int> zfind = mx_set.dep(1).sparsity().find();
          std::string zi = "snz_zi_" + uniq, vg = "snz_zg_" + uniq;
          add_int_constant(add_node, zi, zfind);
          create_gather_node(add_node, vals, zi, vg, 1);  // axis 1 (column-flat on import)
          vals = vg;
        }
        std::string upd = "snz_upd_" + uniq, ush = "snz_ush_" + uniq;
        add_int_constant(add_node, ush, {static_cast<casadi_int>(idx.size())});  // 1-D shape [nnz]
        node = add_node();
        node->set_op_type("Reshape");
        node->add_input(vals); node->add_input(ush); node->add_output(upd);

        std::string scat = "snz_out_" + uniq;
        node = add_node();
        node->set_op_type("ScatterND");
        node->add_input(work_to_onnx[i_vec[0]]);  // 2-D pseudo-dense x (no flatten)
        node->add_input(ind);
        node->add_input(upd);
        node->add_output(scat);
        if (add) {  // accumulate onto x and over duplicate indices (opset 16 reduction)
          onnx::AttributeProto* red = node->add_attribute();
          red->set_name("reduction");
          red->set_type(onnx::AttributeProto::STRING);
          red->set_s("add");
        }
        // Restore the instruction's output pattern natively from the seed (no overlay). ScatterND
        // imports as a DENSE value (the importer densifies the target), so seed mx_set.sparsity().
        std::string r = emit_sparsity_restore(add_node, scat,
                                              Sparsity::dense(mx_set.size1(), mx_set.size2()),
                                              mx_set.sparsity(), uniq, node_output);
        work_to_onnx[o_vec[0]] = r;
        return true;
      }

      case OP_ATAN2: {
        // atan2(y, x) = atan(y/x) + correction for x<0, with origin→0
        std::string y = work_to_onnx[i_vec[0]], x = work_to_onnx[i_vec[1]];
        std::string p = "atan2_" + std::to_string(k) + "_";

        add_real_constant(add_node, p+"c0", {0.0});
        add_real_constant(add_node, p+"pi", {M_PI});
        create_binary_node(add_node, "Div", y, x, p+"r");
        create_unary_node(add_node, "Atan", p+"r", p+"b");
        create_binary_node(add_node, "Less", x, p+"c0", p+"xn");
        create_binary_node(add_node, "Less", y, p+"c0", p+"yn");
        create_unary_node(add_node, "Neg", p+"pi", p+"npi");
        create_where_node(add_node, p+"yn", p+"npi", p+"pi", p+"corr"); // x<0: -π if y<0 else +π
        create_where_node(add_node, p+"xn", p+"corr", p+"c0", p+"adj"); // 0 when x>=0
        create_binary_node(add_node, "Add", p+"b", p+"adj", p+"res");
        create_binary_node(add_node, "Equal", x, p+"c0", p+"xz");
        create_binary_node(add_node, "Equal", y, p+"c0", p+"yz");
        create_binary_node(add_node, "And", p+"xz", p+"yz", p+"oz");
        create_where_node(add_node, p+"oz", p+"c0", p+"res", node_output); // origin → 0

        work_to_onnx[o_vec[0]] = node_output;
        return true;
      }

      // Handled by the caller (OP_CALL: control flow / function call) or unsupported.
      case OP_CALL:
      case OP_SUBREF:
      default:
        return false;
    }
  }

  // GraphProto convenience wrapper for process_operation
  bool Onnx::process_operation(
      onnx::GraphProto* graph,
      const Function& f,
      casadi_int op,
      casadi_int k,
      const std::vector<casadi_int>& i_vec,
      const std::vector<casadi_int>& o_vec,
      std::map<casadi_int, std::string>& work_to_onnx,
      const std::string& node_output) {
    return process_operation([graph]() { return graph->add_node(); },
                              f, op, k, i_vec, o_vec, work_to_onnx, node_output);
  }

} // namespace casadi
/// \endcond
