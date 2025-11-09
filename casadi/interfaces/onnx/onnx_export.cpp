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
#include <casadi/core/casadi_meta.hpp>
#include <memory>
#include <set>

/// \cond INTERNAL
namespace casadi {

  std::string onnx_input_name(const Function& f, casadi_int i) {
    std::string n = f.name_in(i);
    return n.empty() ? "input_" + std::to_string(i) : n;
  }

  std::string onnx_output_name(const Function& f, casadi_int i) {
    std::string n = f.name_out(i);
    return n.empty() ? "output_" + std::to_string(i) : n;
  }

  void Onnx::set_real_tensor_type(onnx::ValueInfoProto* value, const Sparsity& sp) {
    // Transpose-representation invariant: an ONNX tensor stores its CasADi value's column-major
    // bytes with the shape declared REVERSED, so a CasADi (r,c) value is an ONNX (c,r) tensor.
    onnx::TypeProto::Tensor* tensor_type = value->mutable_type()->mutable_tensor_type();
    tensor_type->set_elem_type(real_type());
    onnx::TensorShapeProto* shape = tensor_type->mutable_shape();
    shape->add_dim()->set_dim_value(sp.size2());
    shape->add_dim()->set_dim_value(sp.size1());
  }

  // Add graph inputs. Subgraphs use a prefixed name for uniqueness; the main graph (empty
  // prefix) uses the function's own input names.
  void Onnx::add_graph_inputs(onnx::GraphProto* graph, const Function& f,
                                        const std::string& name_prefix) {
    for (casadi_int i = 0; i < f.n_in(); ++i) {
      onnx::ValueInfoProto* input = graph->add_input();
      input->set_name(name_prefix.empty() ? onnx_input_name(f, i)
                                           : name_prefix + "_i" + std::to_string(i));
      set_real_tensor_type(input, f.sparsity_in(i));
    }
  }

  // Add graph outputs (see add_graph_inputs for the naming convention)
  void Onnx::add_graph_outputs(onnx::GraphProto* graph, const Function& f,
                                         const std::string& name_prefix) {
    for (casadi_int i = 0; i < f.n_out(); ++i) {
      onnx::ValueInfoProto* output = graph->add_output();
      output->set_name(name_prefix.empty() ? onnx_output_name(f, i)
                                            : name_prefix + "_o" + std::to_string(i));
      set_real_tensor_type(output, f.sparsity_out(i));
    }
  }

  void Onnx::load(const Function& f) {
    model_.Clear();
    exported_functions_.clear();
    model_.set_ir_version(8);
    model_.set_producer_name("CasADi");
    model_.set_producer_version(CasadiMeta::version());

    onnx::OperatorSetIdProto* opset = model_.add_opset_import();
    opset->set_domain("");   // default ONNX domain
    opset->set_version(16);  // opset 16: ScatterND reduction

    onnx::GraphProto* graph = model_.mutable_graph();
    graph->set_name(f.name());
    add_graph_inputs(graph, f);

    // work vector index -> the ONNX tensor name that produced it
    std::map<casadi_int, std::string> work_to_onnx;

    casadi_int n_instr = f.n_instructions();

    // Pre-scan: count segments per output (multi-segment = horz/vert/diagcat output)
    std::map<casadi_int, casadi_int> output_segment_count;
    for (casadi_int k = 0; k < n_instr; ++k) {
      if (f.instruction_id(k) == OP_OUTPUT) {
        MX mx = f.instruction_MX(k);
        Dict info = mx.info();
        casadi_int output_idx = info["ind"];
        output_segment_count[output_idx]++;
      }
    }

    // Track segment values for multi-segment outputs: output_idx -> (offset -> onnx_name)
    std::map<casadi_int, std::map<casadi_int, std::string>> output_segment_values;
    // and each segment's source sparsity, to extract its nonzero array on reassembly
    std::map<casadi_int, std::map<casadi_int, Sparsity>> output_segment_sparsity;

    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i = f.instruction_input(k);

      // Unique result name keyed on the instruction index
      std::string node_output = "n" + std::to_string(k);

      // OP_OUTPUT is special: it supports multi-segment outputs (horz/vert/diagcat)
      if (op == OP_OUTPUT) {
        MX mx = f.instruction_MX(k);
        Dict info = mx.info();
        casadi_int output_idx = info["ind"];
        casadi_int offset = info["offset"];

        std::string output_name = onnx_output_name(f, output_idx);
        std::string input_onnx_name = work_to_onnx[i[0]];

        if (output_segment_count[output_idx] > 1) {
          // Multi-segment output (horzcat-like): each segment provides a CONTIGUOUS run of the
          // output's column-major nonzero array (offset = nnz offset). Record the source value and
          // its sparsity; reassemble by concatenating the segments' nonzero arrays below.
          output_segment_values[output_idx][offset] = input_onnx_name;
          output_segment_sparsity[output_idx][offset] = mx.dep(0).sparsity();
        } else {
          // Single-segment output. A reshape can fold into the output (no OP_RESHAPE
          // instruction): the source shape then differs from the declared output shape. A folded
          // reshape densifies on import, so when the output pattern is non-dense, route the assembly
          // through a temp and restore the exact pattern natively from the seed (no overlay).
          MX dep = mx.dep(0);
          Sparsity out_sp = f.sparsity_out(output_idx);
          auto add_node = [graph]() { return graph->add_node(); };
          bool folded = (dep.size1() != out_sp.size1() || dep.size2() != out_sp.size2());
          if (folded && !out_sp.is_dense()) {
            std::string tmp = "out_pre_" + std::to_string(k);
            emit_output_node(graph, input_onnx_name, dep.size1(), dep.size2(), out_sp,
                             tmp, "out_rs_" + std::to_string(k));
            emit_sparsity_restore(add_node, tmp,
                                Sparsity::dense(out_sp.size1(), out_sp.size2()), out_sp,
                                "outr_" + std::to_string(k), output_name);
          } else {
            emit_output_node(graph, input_onnx_name, dep.size1(), dep.size2(), out_sp,
                             output_name, "out_rs_" + std::to_string(k));
          }
        }
        continue;
      }

      if (process_operation(graph, f, op, k, i, o, work_to_onnx, node_output)) continue;

      // Operations not handled by process_operation
      if (op == OP_CALL) {
        Function called_func = f.instruction_MX(k).which_function();
        if (is_map_function(called_func)) {
          export_map(graph, called_func, i, o, work_to_onnx, node_output + "_out");
        } else if (is_reduce_map_function(called_func)) {
          export_reduce_map(graph, called_func, i, o, work_to_onnx, node_output + "_out");
        } else if (is_if_else_function(called_func)) {
          export_if(graph, called_func, i, o, work_to_onnx, node_output + "_out");
        } else {
          assert_not_control_flow(called_func);
          export_call(graph, called_func, i, o, work_to_onnx, node_output + "_out");
        }
        continue;
      }

      // Unknown/unsupported operation
      casadi_error("ONNX export: unsupported operation code " +
                  std::to_string(op) + " at instruction " + std::to_string(k) +
                  ". The CasADi Function contains operations that cannot be exported to ONNX.");
    }

    // Reassemble multi-segment outputs pseudo-dense. The block structure is one of:
    //  - horzcat (all blocks full-height, size1==R): columns concatenated -> a single ONNX Concat.
    //    Transpose-rep: CasADi horzcat (axis 1) is ONNX axis 0 (matches OP_HORZCAT dispatch).
    //  - vertcat (all blocks full-width, size2==C): rows concatenated -> a single ONNX Concat on
    //    the other axis. Transpose-rep: CasADi vertcat (axis 0) is ONNX axis 1 (matches OP_VERTCAT).
    //  - diagcat (neither, genuine block-diagonal): Pad each DENSE block to (R,C) at its offset, Sum.
    // The assembly imports as a DENSE block; when the output pattern is non-dense, restore it
    // natively from the seed (no overlay) by routing the assembly through a temp.
    auto seg_add_node = [graph]() { return graph->add_node(); };
    for (const auto& kv : output_segment_values) {
      casadi_int output_idx = kv.first;
      const auto& segments = kv.second;  // offset -> onnx_name, sorted by offset
      const auto& seg_sp = output_segment_sparsity[output_idx];
      casadi_int R = f.size1_out(output_idx), C = f.size2_out(output_idx);

      std::vector<std::string> names;
      std::vector<casadi_int> brs, bcs;
      bool all_full_height = true, all_full_width = true;
      for (const auto& seg : segments) {
        const Sparsity& sp = seg_sp.at(seg.first);
        names.push_back(seg.second);
        brs.push_back(sp.size1()); bcs.push_back(sp.size2());
        if (sp.size1() != R) all_full_height = false;  // not horzcat
        if (sp.size2() != C) all_full_width = false;   // not vertcat
      }
      Sparsity out_sp = f.sparsity_out(output_idx);
      std::string oname = onnx_output_name(f, output_idx);
      std::string uniq = "oseg" + std::to_string(output_idx);

      // Helper: write the assembly directly to oname when dense, else to a temp + restore.
      auto finish = [&](const std::function<void(const std::string&)>& emit) {
        if (out_sp.is_dense()) {
          emit(oname);
        } else {
          std::string tmp = oname + "_pre";
          emit(tmp);
          emit_sparsity_restore(seg_add_node, tmp, Sparsity::dense(R, C), out_sp, uniq + "r", oname);
        }
      };

      if (all_full_height) {
        // horzcat -> ONNX Concat axis 0
        finish([&](const std::string& dst) {
          onnx::NodeProto* cc = seg_add_node();
          cc->set_op_type("Concat");
          for (const auto& n : names) cc->add_input(n);
          cc->add_output(dst);
          add_int_attribute(cc, "axis", 0);
        });
      } else if (all_full_width) {
        // vertcat -> ONNX Concat axis 1
        finish([&](const std::string& dst) {
          onnx::NodeProto* cc = seg_add_node();
          cc->set_op_type("Concat");
          for (const auto& n : names) cc->add_input(n);
          cc->add_output(dst);
          add_int_attribute(cc, "axis", 1);
        });
      } else {
        // genuine block-diagonal (diagcat) -> Pad each block to (R,C) at its offset, then Sum
        std::vector<casadi_int> row_off, col_off;
        casadi_int ro = 0, co = 0;
        for (casadi_int b = 0; b < static_cast<casadi_int>(names.size()); ++b) {
          row_off.push_back(ro); col_off.push_back(co);
          ro += brs[b]; co += bcs[b];
        }
        finish([&](const std::string& dst) {
          emit_blockdiag(seg_add_node, names, row_off, col_off, brs, bcs, R, C, dst, uniq);
        });
      }
    }

    // Fuse trailing rename Identities: a node of the form Identity(src) -> oname, where oname is a
    // declared graph output and src is an internal node output used nowhere else, is a pure rename.
    // Make the producing node write straight to oname and drop the Identity. Conservative: only when
    // src is produced by exactly one node and consumed by exactly this Identity (so we never collapse
    // an output that aliases an input, or a value feeding two outputs / another consumer).
    {
      std::set<std::string> output_names;
      for (casadi_int i = 0; i < f.n_out(); ++i) output_names.insert(onnx_output_name(f, i));
      // How many places consume each tensor name (as a node input).
      std::map<std::string, int> consumers;
      for (const auto& nd : graph->node())
        for (const auto& in : nd.input()) consumers[in]++;
      // Which node produces each tensor name, and is the name produced more than once.
      std::map<std::string, int> producers;
      for (const auto& nd : graph->node())
        for (const auto& o : nd.output()) producers[o]++;

      auto* nodes = graph->mutable_node();
      // Plan rewrites: src -> oname for fusable trailing Identities.
      std::map<std::string, std::string> rename;     // src -> oname
      std::set<int> drop;                             // node indices (Identity) to remove
      for (int n = 0; n < nodes->size(); ++n) {
        const onnx::NodeProto& nd = nodes->Get(n);
        if (nd.op_type() != "Identity" || nd.input_size() != 1 || nd.output_size() != 1) continue;
        const std::string& src = nd.input(0);
        const std::string& dst = nd.output(0);
        if (!output_names.count(dst)) continue;        // only fuse into graph outputs
        if (output_names.count(src)) continue;         // src is itself an output (e.g. aliased)
        if (producers[src] != 1) continue;             // src must be produced by a single node
        if (consumers[src] != 1) continue;             // and consumed only by this Identity
        if (rename.count(src) || drop.count(n)) continue;  // src already claimed by another output
        rename[src] = dst;
        drop.insert(n);
      }
      if (!drop.empty()) {
        // Apply: rewrite producer outputs, then rebuild node list without the dropped Identities.
        for (int n = 0; n < nodes->size(); ++n) {
          if (drop.count(n)) continue;
          onnx::NodeProto* nd = nodes->Mutable(n);
          for (int o = 0; o < nd->output_size(); ++o) {
            auto it = rename.find(nd->output(o));
            if (it != rename.end()) nd->set_output(o, it->second);
          }
        }
        google::protobuf::RepeatedPtrField<onnx::NodeProto> kept;
        for (int n = 0; n < nodes->size(); ++n)
          if (!drop.count(n)) kept.Add()->CopyFrom(nodes->Get(n));
        nodes->Swap(&kept);
      }
    }

    // Add graph outputs (for main graph, use empty prefix to use function's output names)
    add_graph_outputs(graph, f, "");

    // An output not produced by any node (e.g. an all-structural-zero output, whose MXFunction has
    // no instructions) needs an explicit zero Constant so both ORT and the importer have a tensor.
    // A non-dense pattern is planted directly as a sparse_value (all-zero) Constant -- the seed that
    // imports to exactly that pattern; a dense output gets a plain dense zero Constant.
    {
      std::set<std::string> produced;
      for (const auto& nd : graph->node())
        for (const auto& o : nd.output()) produced.insert(o);
      auto add_node = [graph]() { return graph->add_node(); };
      for (casadi_int i = 0; i < f.n_out(); ++i) {
        std::string oname = onnx_output_name(f, i);
        if (!produced.count(oname)) {
          const Sparsity& sp = f.sparsity_out(i);
          if (sp.is_dense()) {
            add_real_constant(add_node, oname,
                              std::vector<double>(sp.size1() * sp.size2(), 0.0),
                              {sp.size2(), sp.size1()});
          } else {
            add_sparse_constant(add_node, oname, DM(sp, 0.0));
          }
        }
      }
    }

    // Input sparsity seed: the ONNX graph value-flow is dense, so a non-dense CasADi input pattern
    // is planted as a standard sparse_initializer (all-zero values) sharing the input name -- an
    // optional input with a sparse zero default, readable by any ONNX tool. Import picks it up as a
    // sparse MX and re-propagates sparsity natively. OUTPUT patterns need no overlay: each op
    // restores its own output sparsity from seeds (emit_sparsity_restore), so they recover for free.
    for (casadi_int i = 0; i < f.n_in(); ++i) {
      if (!f.sparsity_in(i).is_dense()) {
        fill_sparse_tensor(graph->add_sparse_initializer(), onnx_input_name(f, i),
                           DM(f.sparsity_in(i), 0.0));
      }
    }

    // If we exported any functions, add the casadi domain opset_import
    if (!exported_functions_.empty()) {
      onnx::OperatorSetIdProto* casadi_opset = model_.add_opset_import();
      casadi_opset->set_domain("casadi");
      casadi_opset->set_version(1);

      if (verbose_) {
        uout() << "  Exported " << exported_functions_.size()
               << " function(s) to casadi domain" << std::endl;
      }
    }

    has_model_ = true;

    if (verbose_) {
      uout() << "Converted CasADi Function to ONNX model: " << f.name() << std::endl;
      uout() << "  Instructions processed: " << n_instr << std::endl;
      uout() << "  ONNX nodes created: " << graph->node_size() << std::endl;
    }
  }

  // ========== Control Flow Support ==========

  bool Onnx::is_if_else_function(const Function& f) const {
    return f.class_name() == "Switch";
  }

  // map and mapaccum both have class_name "Map"/"OmpMap"; the name distinguishes them.
  bool Onnx::is_mapaccum_function(const Function& f) const {
    std::string class_name = f.class_name();
    if (class_name != "Map" && class_name != "OmpMap") return false;
    std::string fname = f.name();
    return fname.find("mapaccum") != std::string::npos ||
           fname.find("accum") != std::string::npos;
  }

  bool Onnx::is_map_function(const Function& f) const {
    std::string class_name = f.class_name();
    if (class_name != "Map" && class_name != "OmpMap") return false;
    std::string fname = f.name();
    return fname.find("mapaccum") == std::string::npos &&
           fname.find("accum") == std::string::npos;
  }

  bool Onnx::is_reduce_map_function(const Function& f) const {
    // A reduce-map is an MXFunction wrapping a single MapSum; spot the MapSum sub-function
    for (const std::string& nm : f.get_function()) {
      if (f.get_function(nm).class_name() == "MapSum") return true;
    }
    return false;
  }

  void Onnx::assert_not_control_flow(const Function& called_func) const {
    casadi_assert(!is_mapaccum_function(called_func), "ONNX export: mapaccum (Loop) is not supported.");
  }

  // Prefix the names a subgraph *defines* (graph inputs, initializers, node outputs) so they
  // can't collide with the enclosing graph. References to names defined elsewhere are left
  // untouched, which preserves an If branch's implicit captures of outer-scope tensors.
  static void prefix_graph_names(onnx::GraphProto* g, const std::string& prefix) {
    std::set<std::string> defined;
    for (int i = 0; i < g->input_size(); ++i) defined.insert(g->input(i).name());
    for (int i = 0; i < g->initializer_size(); ++i) defined.insert(g->initializer(i).name());
    for (int n = 0; n < g->node_size(); ++n) {
      for (int j = 0; j < g->node(n).output_size(); ++j) defined.insert(g->node(n).output(j));
    }
    for (int n = 0; n < g->node_size(); ++n) {
      onnx::NodeProto* nd = g->mutable_node(n);
      for (int i = 0; i < nd->input_size(); ++i) {
        if (defined.count(nd->input(i))) nd->set_input(i, prefix + nd->input(i));
      }
      for (int i = 0; i < nd->output_size(); ++i) nd->set_output(i, prefix + nd->output(i));
    }
    for (int i = 0; i < g->input_size(); ++i) {
      g->mutable_input(i)->set_name(prefix + g->input(i).name());
    }
    for (int i = 0; i < g->output_size(); ++i) {
      g->mutable_output(i)->set_name(prefix + g->output(i).name());
    }
  }

  template<typename Container>
  void Onnx::emit_reshape(Container* container, const std::string& data,
                                    const std::vector<casadi_int>& shape,
                                    const std::string& output, const std::string& shape_name) {
    onnx::NodeProto* sc = container->add_node();
    sc->set_op_type("Constant");
    sc->add_output(shape_name);
    onnx::AttributeProto* attr = sc->add_attribute();
    attr->set_name("value");
    attr->set_type(onnx::AttributeProto::TENSOR);
    onnx::TensorProto* t = attr->mutable_t();
    t->set_data_type(onnx::TensorProto::INT64);
    t->add_dims(static_cast<casadi_int>(shape.size()));
    for (casadi_int v : shape) t->add_int64_data(v);
    onnx::NodeProto* r = container->add_node();
    r->set_op_type("Reshape");
    r->add_input(data);
    r->add_input(shape_name);
    r->add_output(output);
  }

  template<typename Container>
  void Onnx::colmajor_reshape_into(Container* container, const std::string& data,
                                             const std::vector<casadi_int>& dims,
                                             const std::string& output, const std::string& uniq) {
    // Transpose-rep: a CasADi column-major reshape to `dims` is a plain row-major Reshape of the
    // stored (already transposed) tensor to the REVERSED dims -- no Transpose nodes.
    std::vector<casadi_int> rev(dims.rbegin(), dims.rend());
    emit_reshape(container, data, rev, output, uniq + "_s");
  }

  template<typename Container>
  void Onnx::emit_output_node(Container* container, const std::string& data,
                              casadi_int src_rows, casadi_int src_cols, const Sparsity& out_sp,
                              const std::string& output, const std::string& uniq) {
    // A reshape folded into the output shows up as a source/output shape mismatch
    if (src_rows != out_sp.size1() || src_cols != out_sp.size2()) {
      colmajor_reshape_into(container, data, {out_sp.size1(), out_sp.size2()}, output, uniq);
    } else {
      onnx::NodeProto* id = container->add_node();
      id->set_op_type("Identity");
      id->add_input(data);
      id->add_output(output);
    }
  }

  onnx::GraphProto Onnx::build_scan_body(const Function& base) {
    // Scan slices the iteration axis off the 3-D inputs, so the body runs the base on a
    // clean 2-D (rows, c) block per iteration -- no rank fixing needed.
    onnx::GraphProto body;
    body.set_name(base.name() + "_scan_body");
    std::map<casadi_int, std::string> work_to_onnx;

    for (casadi_int j = 0; j < base.n_in(); ++j) {
      onnx::ValueInfoProto* vi = body.add_input();
      vi->set_name("body_in_" + std::to_string(j));
      set_real_tensor_type(vi, base.sparsity_in(j));
    }

    casadi_int n_instr = base.n_instructions();
    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = base.instruction_id(k);
      std::vector<casadi_int> o = base.instruction_output(k);
      std::vector<casadi_int> i_vec = base.instruction_input(k);
      std::string node_output = "n" + std::to_string(k);

      if (op == OP_INPUT) {
        work_to_onnx[o[0]] = "body_in_" + std::to_string(i_vec[0]);
        continue;
      }

      if (op == OP_OUTPUT) {
        std::string out_name = "body_out_" + std::to_string(o[0]);
        MX dep = base.instruction_MX(k).dep(0);
        Sparsity out_sp = base.sparsity_out(o[0]);
        emit_output_node(&body, work_to_onnx[i_vec[0]], dep.size1(), dep.size2(), out_sp,
                         out_name, "body_out_rs_" + std::to_string(k));
        onnx::ValueInfoProto* vi = body.add_output();
        vi->set_name(out_name);
        set_real_tensor_type(vi, out_sp);
        continue;
      }

      if (op == OP_CALL) {
        Function called = base.instruction_MX(k).which_function();
        if (is_map_function(called)) {
          export_map(&body, called, i_vec, o, work_to_onnx, node_output + "_out");
        } else if (is_reduce_map_function(called)) {
          export_reduce_map(&body, called, i_vec, o, work_to_onnx, node_output + "_out");
        } else if (is_if_else_function(called)) {
          export_if(&body, called, i_vec, o, work_to_onnx, node_output + "_out");
        } else {
          assert_not_control_flow(called);
          export_call(&body, called, i_vec, o, work_to_onnx, node_output + "_out");
        }
        continue;
      }

      if (process_operation(&body, base, op, k, i_vec, o, work_to_onnx, node_output)) continue;

      casadi_error("ONNX export: unsupported operation code " + std::to_string(op) +
                  " in Map body of '" + base.name() + "'");
    }
    return body;
  }

  template<typename Container>
  void Onnx::export_map(Container* container, const Function& map_fn,
                                  const std::vector<casadi_int>& i_vec,
                                  const std::vector<casadi_int>& o,
                                  std::map<casadi_int, std::string>& work_to_onnx,
                                  const std::string& out_prefix) {
    Function base = map_fn.get_function(map_fn.get_function().at(0));
    casadi_int n = map_fn.size2_out(0) / base.size2_out(0);

    // Only the plain Map is handled: every input repeated n times, every output concatenated.
    // reduce_in (broadcast inputs) / reduce_out (summed outputs) are not yet supported.
    for (casadi_int j = 0; j < base.n_in(); ++j) {
      casadi_assert(map_fn.size2_in(j) == n * base.size2_in(j),
        "ONNX export: Map with non-repeated (reduce_in) inputs is not yet supported.");
    }
    for (casadi_int j = 0; j < base.n_out(); ++j) {
      casadi_assert(map_fn.size2_out(j) == n * base.size2_out(j),
        "ONNX export: Map with reduced (reduce_out) outputs is not yet supported.");
    }

    // Transpose-rep: the stored input is (n*c, rows); lift to 3-D (n, c, rows) so Scan slices the
    // leading n-axis into 2-D (c, rows) blocks (each the stored transpose of a (rows, c) block).
    std::vector<std::string> scan_inputs;
    for (casadi_int j = 0; j < base.n_in(); ++j) {
      std::string lifted = out_prefix + "_in" + std::to_string(j);
      emit_reshape(container, work_to_onnx[i_vec[j]],
                   {n, base.size2_in(j), base.size1_in(j)}, lifted,
                   out_prefix + "_insh" + std::to_string(j));
      scan_inputs.push_back(lifted);
    }

    onnx::NodeProto* scan = container->add_node();
    scan->set_op_type("Scan");
    for (const std::string& s : scan_inputs) scan->add_input(s);
    std::vector<std::string> scan_outputs;
    for (casadi_int j = 0; j < base.n_out(); ++j) {
      std::string so = out_prefix + "_s" + std::to_string(j);
      scan->add_output(so);
      scan_outputs.push_back(so);
    }
    add_int_attribute(scan, "num_scan_inputs", base.n_in());
    add_ints_attribute(scan, "scan_input_axes", std::vector<casadi_int>(base.n_in(), 0));
    add_ints_attribute(scan, "scan_output_axes", std::vector<casadi_int>(base.n_out(), 0));

    onnx::GraphProto body = build_scan_body(base);
    prefix_graph_names(&body, out_prefix + "_b_");  // keep body names out of the outer scope

    onnx::AttributeProto* body_attr = scan->add_attribute();
    body_attr->set_name("body");
    body_attr->set_type(onnx::AttributeProto::GRAPH);
    *body_attr->mutable_g() = body;

    // Transpose-rep: stacked scan output is (n, c, rows); flatten to (n*c, rows) = stored output
    for (casadi_int j = 0; j < base.n_out(); ++j) {
      std::string out = out_prefix + std::to_string(j);
      emit_reshape(container, scan_outputs[j],
                   {n * base.size2_out(j), base.size1_out(j)}, out,
                   out_prefix + "_outsh" + std::to_string(j));
      work_to_onnx[o[j]] = out;
    }
  }

  onnx::GraphProto Onnx::build_reduce_scan_body(const Function& base,
      const std::vector<bool>& reduce_in, const std::vector<bool>& reduce_out,
      const std::vector<std::string>& capture_names) {
    // The body calls the base once; reduce_in args are captured from the outer scope,
    // reduce_out results are accumulated into state variables (ONNX needs states declared first).
    onnx::GraphProto body;
    body.set_name(base.name() + "_redscan_body");
    std::map<casadi_int, std::string> w;

    for (casadi_int j = 0; j < base.n_out(); ++j) {
      if (!reduce_out[j]) continue;
      onnx::ValueInfoProto* vi = body.add_input();
      vi->set_name("acc_in_" + std::to_string(j));
      set_real_tensor_type(vi, base.sparsity_out(j));
    }

    std::vector<std::string> args(base.n_in());
    for (casadi_int j = 0; j < base.n_in(); ++j) {
      if (reduce_in[j]) {
        args[j] = capture_names[j];  // implicit outer-scope capture
      } else {
        args[j] = "scan_in_" + std::to_string(j);
        onnx::ValueInfoProto* vi = body.add_input();
        vi->set_name(args[j]);
        set_real_tensor_type(vi, base.sparsity_in(j));
      }
    }

    std::vector<casadi_int> ii(base.n_in()), oo(base.n_out());
    for (casadi_int j = 0; j < base.n_in(); ++j) { ii[j] = j; w[j] = args[j]; }
    for (casadi_int j = 0; j < base.n_out(); ++j) oo[j] = base.n_in() + j;
    export_call(&body, base, ii, oo, w, "bcall_");  // base outputs land in w[oo[j]]

    for (casadi_int j = 0; j < base.n_out(); ++j) {
      if (!reduce_out[j]) continue;
      std::string out = "acc_out_" + std::to_string(j);
      create_binary_node(&body, "Add", "acc_in_" + std::to_string(j), w[oo[j]], out);
      onnx::ValueInfoProto* vi = body.add_output();
      vi->set_name(out);
      set_real_tensor_type(vi, base.sparsity_out(j));
    }
    for (casadi_int j = 0; j < base.n_out(); ++j) {
      if (reduce_out[j]) continue;
      std::string out = "yscan_" + std::to_string(j);
      create_unary_node(&body, "Identity", w[oo[j]], out);
      onnx::ValueInfoProto* vi = body.add_output();
      vi->set_name(out);
      set_real_tensor_type(vi, base.sparsity_out(j));
    }
    return body;
  }

  template<typename Container>
  void Onnx::export_reduce_map(Container* container, const Function& wrapper,
                                         const std::vector<casadi_int>& i_vec,
                                         const std::vector<casadi_int>& o,
                                         std::map<casadi_int, std::string>& work_to_onnx,
                                         const std::string& out_prefix) {
    // Unwrap the MXFunction(wrapper) -> MapSum -> base function
    Function mapsum;
    for (const std::string& nm : wrapper.get_function()) {
      if (wrapper.get_function(nm).class_name() == "MapSum") { mapsum = wrapper.get_function(nm); break; }
    }
    Function base = mapsum.get_function(mapsum.get_function().at(0));

    // Repeated I/O is n*base wide; reduce_in/reduce_out keep base width
    std::vector<bool> reduce_in(base.n_in()), reduce_out(base.n_out());
    for (casadi_int j = 0; j < base.n_in(); ++j)
      reduce_in[j] = wrapper.size2_in(j) == base.size2_in(j);
    for (casadi_int j = 0; j < base.n_out(); ++j)
      reduce_out[j] = wrapper.size2_out(j) == base.size2_out(j);
    casadi_int n = 0;
    for (casadi_int j = 0; j < base.n_in(); ++j)
      if (!reduce_in[j]) { n = wrapper.size2_in(j) / base.size2_in(j); break; }
    if (n == 0) for (casadi_int j = 0; j < base.n_out(); ++j)
      if (!reduce_out[j]) { n = wrapper.size2_out(j) / base.size2_out(j); break; }
    casadi_assert(n > 0, "ONNX export: reduce-map with no repeated inputs or outputs.");

    AddNodeFn add_node = [container]() { return container->add_node(); };

    // Alias each reduce_in input so the body captures a uniquely-named outer tensor
    std::vector<std::string> capture_names(base.n_in());
    for (casadi_int j = 0; j < base.n_in(); ++j) {
      if (!reduce_in[j]) continue;
      capture_names[j] = out_prefix + "_cap" + std::to_string(j);
      create_unary_node(add_node, "Identity", work_to_onnx[i_vec[j]], capture_names[j]);
    }

    // Lift each repeated input (rows, n*c) -> 3-D (rows, n, c) so Scan slices the n-axis
    std::vector<std::string> scan_inputs;
    for (casadi_int j = 0; j < base.n_in(); ++j) {
      if (reduce_in[j]) continue;
      std::string lifted = out_prefix + "_in" + std::to_string(j);
      emit_reshape(container, work_to_onnx[i_vec[j]],
                   {n, base.size2_in(j), base.size1_in(j)}, lifted,
                   out_prefix + "_insh" + std::to_string(j));
      scan_inputs.push_back(lifted);
    }

    // Zero initial accumulators for each reduce_out output (transpose-rep shape (c,r))
    std::vector<std::string> acc_inits;
    for (casadi_int j = 0; j < base.n_out(); ++j) {
      if (!reduce_out[j]) continue;
      std::string nm = out_prefix + "_acc0_" + std::to_string(j);
      std::vector<double> zeros(base.size1_out(j) * base.size2_out(j), 0.0);
      add_real_constant(add_node, nm, zeros, {base.size2_out(j), base.size1_out(j)});
      acc_inits.push_back(nm);
    }

    onnx::NodeProto* scan = add_node();
    scan->set_op_type("Scan");
    for (const std::string& s : acc_inits) scan->add_input(s);
    for (const std::string& s : scan_inputs) scan->add_input(s);

    // Scan outputs: state accumulators (reduce_out) first, then concatenated scan outputs
    std::vector<std::string> state_out, scan_out;
    for (casadi_int j = 0; j < base.n_out(); ++j)
      if (reduce_out[j]) { state_out.push_back(out_prefix + "_acc" + std::to_string(j)); scan->add_output(state_out.back()); }
    for (casadi_int j = 0; j < base.n_out(); ++j)
      if (!reduce_out[j]) { scan_out.push_back(out_prefix + "_s" + std::to_string(j)); scan->add_output(scan_out.back()); }

    add_int_attribute(scan, "num_scan_inputs", static_cast<casadi_int>(scan_inputs.size()));
    add_ints_attribute(scan, "scan_input_axes", std::vector<casadi_int>(scan_inputs.size(), 0));
    add_ints_attribute(scan, "scan_output_axes", std::vector<casadi_int>(scan_out.size(), 0));

    onnx::GraphProto body = build_reduce_scan_body(base, reduce_in, reduce_out, capture_names);
    prefix_graph_names(&body, out_prefix + "_b_");  // captures stay unprefixed (defined outside)
    onnx::AttributeProto* battr = scan->add_attribute();
    battr->set_name("body");
    battr->set_type(onnx::AttributeProto::GRAPH);
    *battr->mutable_g() = body;

    // reduce_out -> accumulator output directly; repeated -> flatten (rows, n, c) -> (rows, n*c)
    casadi_int si = 0, ci = 0;
    for (casadi_int j = 0; j < base.n_out(); ++j) {
      if (reduce_out[j]) {
        work_to_onnx[o[j]] = state_out[si++];
      } else {
        std::string out = out_prefix + std::to_string(j);
        emit_reshape(container, scan_out[ci++], {n * base.size2_out(j), base.size1_out(j)}, out,
                     out_prefix + "_outsh" + std::to_string(j));
        work_to_onnx[o[j]] = out;
      }
    }
  }

  onnx::GraphProto Onnx::build_if_branch(const Function& f,
      const std::vector<std::string>& arg_names, const std::string& prefix) {
    // An ONNX If branch has no formal inputs: it captures the outer arg tensors by name.
    onnx::GraphProto g;
    g.set_name(f.name() + "_branch");
    std::map<casadi_int, std::string> work_to_onnx;

    casadi_int n_instr = f.n_instructions();
    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i_vec = f.instruction_input(k);
      std::string node_output = "n" + std::to_string(k);

      if (op == OP_INPUT) {
        work_to_onnx[o[0]] = arg_names.at(i_vec[0]);  // capture the outer tensor
        continue;
      }
      if (op == OP_OUTPUT) {
        std::string out_name = "branch_out_" + std::to_string(o[0]);
        MX dep = f.instruction_MX(k).dep(0);
        Sparsity out_sp = f.sparsity_out(o[0]);
        emit_output_node(&g, work_to_onnx[i_vec[0]], dep.size1(), dep.size2(), out_sp,
                         out_name, "branch_out_rs_" + std::to_string(k));
        onnx::ValueInfoProto* vi = g.add_output();
        vi->set_name(out_name);
        set_real_tensor_type(vi, out_sp);
        continue;
      }
      if (op == OP_CALL) {
        Function called = f.instruction_MX(k).which_function();
        if (is_map_function(called)) {
          export_map(&g, called, i_vec, o, work_to_onnx, node_output + "_out");
        } else if (is_reduce_map_function(called)) {
          export_reduce_map(&g, called, i_vec, o, work_to_onnx, node_output + "_out");
        } else if (is_if_else_function(called)) {
          export_if(&g, called, i_vec, o, work_to_onnx, node_output + "_out");
        } else {
          assert_not_control_flow(called);
          export_call(&g, called, i_vec, o, work_to_onnx, node_output + "_out");
        }
        continue;
      }
      if (process_operation(&g, f, op, k, i_vec, o, work_to_onnx, node_output)) continue;

      casadi_error("ONNX export: unsupported operation code " + std::to_string(op) +
                  " in if_else branch of '" + f.name() + "'");
    }
    prefix_graph_names(&g, prefix);  // unique internal names; captures of arg_names are kept
    return g;
  }

  template<typename Container>
  void Onnx::export_if(Container* container, const Function& switch_fn,
                                 const std::vector<casadi_int>& i_vec,
                                 const std::vector<casadi_int>& o,
                                 std::map<casadi_int, std::string>& work_to_onnx,
                                 const std::string& out_prefix) {
    // CasADi Switch: index 0 selects f[0], otherwise the default. if_else has a single case,
    // so f[0] is the else-branch and f_def the then-branch. ONNX If is strictly 2-way.
    Dict info = switch_fn.info();
    std::vector<Function> cases = info.at("f");
    Function f_then = info.at("f_def");
    casadi_assert(cases.size() == 1,
      "ONNX export: only 2-way if_else is supported (got a multi-case Switch).");
    Function f_else = cases[0];

    // i_vec[0] is the condition; i_vec[1..] are the data arguments captured by both branches.
    // Alias each through an Identity with a unique name so a branch's internal "n<k>" names
    // can never coincide with a captured outer tensor name.
    std::vector<std::string> arg_names;
    for (size_t j = 1; j < i_vec.size(); ++j) {
      std::string alias = out_prefix + "_arg" + std::to_string(j - 1);
      onnx::NodeProto* id = container->add_node();
      id->set_op_type("Identity");
      id->add_input(work_to_onnx[i_vec[j]]);
      id->add_output(alias);
      arg_names.push_back(alias);
    }

    std::string cond = out_prefix + "_cond";
    onnx::NodeProto* cast = container->add_node();
    cast->set_op_type("Cast");
    cast->add_input(work_to_onnx[i_vec[0]]);
    cast->add_output(cond);
    add_int_attribute(cast, "to", onnx::TensorProto::BOOL);

    onnx::NodeProto* if_node = container->add_node();
    if_node->set_op_type("If");
    if_node->add_input(cond);
    for (size_t j = 0; j < o.size(); ++j) {
      std::string out = out_prefix + std::to_string(j);
      if_node->add_output(out);
      work_to_onnx[o[j]] = out;
    }
    onnx::AttributeProto* then_attr = if_node->add_attribute();
    then_attr->set_name("then_branch");
    then_attr->set_type(onnx::AttributeProto::GRAPH);
    *then_attr->mutable_g() = build_if_branch(f_then, arg_names, out_prefix + "_t_");

    onnx::AttributeProto* else_attr = if_node->add_attribute();
    else_attr->set_name("else_branch");
    else_attr->set_type(onnx::AttributeProto::GRAPH);
    *else_attr->mutable_g() = build_if_branch(f_else, arg_names, out_prefix + "_e_");
  }

  template<typename Container>
  void Onnx::export_call(Container* container, const Function& called_func,
                                   const std::vector<casadi_int>& i_vec,
                                   const std::vector<casadi_int>& o,
                                   std::map<casadi_int, std::string>& work_to_onnx,
                                   const std::string& out_prefix) {
    std::string func_name = called_func.name();
    std::string domain = "casadi";

    // Export the called function as a local FunctionProto the first time we see it
    if (!exported_functions_.count(func_name)) {
      onnx::FunctionProto* func_proto = function_to_function_proto(called_func, domain);
      *model_.add_functions() = *func_proto;
      delete func_proto;
      exported_functions_.insert(func_name);
    }

    onnx::NodeProto* call_node = container->add_node();
    call_node->set_op_type(func_name);
    call_node->set_domain(domain);
    for (casadi_int idx : i_vec) call_node->add_input(work_to_onnx[idx]);
    for (size_t j = 0; j < o.size(); ++j) {
      std::string out = out_prefix + std::to_string(j);
      call_node->add_output(out);
      work_to_onnx[o[j]] = out;
    }
  }

  onnx::FunctionProto* Onnx::function_to_function_proto(
      const Function& f,
      const std::string& domain) {

    // Create a new FunctionProto (guarded so a casadi_assert/error mid-build frees it)
    std::unique_ptr<onnx::FunctionProto> func_guard(new onnx::FunctionProto());
    onnx::FunctionProto* func = func_guard.get();
    func->set_name(f.name());
    func->set_domain(domain);

    // Add opset imports for the function (required by ONNX Runtime)
    // Default ONNX opset
    onnx::OperatorSetIdProto* opset = func->add_opset_import();
    opset->set_domain("");  // Empty string = default ONNX domain
    opset->set_version(16);
    // CasADi domain for nested function calls
    onnx::OperatorSetIdProto* casadi_opset = func->add_opset_import();
    casadi_opset->set_domain("casadi");
    casadi_opset->set_version(1);

    if (verbose_) {
      uout() << "  Converting Function '" << f.name() << "' to ONNX FunctionProto"
             << " (domain: " << domain << ")" << std::endl;
      uout() << "    Inputs: " << f.n_in() << ", Outputs: " << f.n_out()
             << ", Instructions: " << f.n_instructions() << std::endl;
    }

    // Add function inputs (parameter names)
    for (casadi_int i = 0; i < f.n_in(); ++i) {
      func->add_input(onnx_input_name(f, i));
    }

    // Add function outputs (parameter names)
    for (casadi_int i = 0; i < f.n_out(); ++i) {
      func->add_output(onnx_output_name(f, i));
    }

    // Map from work vector index to ONNX node name
    std::map<casadi_int, std::string> work_to_onnx;

    // Multi-segment I/O (a slot read/written in several pieces, as mapaccum produces) needs
    // the Gather/Concat machinery the top-level graph has; reject it here instead of emitting
    // a wrong FunctionProto that reuses an output name across segments.
    std::set<casadi_int> written_outputs;

    // Process instructions and add nodes to the function
    casadi_int n_instr = f.n_instructions();
    for (casadi_int k = 0; k < n_instr; ++k) {
      casadi_int op = f.instruction_id(k);
      std::vector<casadi_int> o = f.instruction_output(k);
      std::vector<casadi_int> i_vec = f.instruction_input(k);

      std::string node_output = "n" + std::to_string(k);

      // OP_INPUT maps a work slot directly to the function parameter name (no node)
      if (op == OP_INPUT) {
        casadi_assert(f.instruction_MX(k).numel() == f.numel_in(i_vec[0]),
          "ONNX export: nested functions with multi-segment inputs (e.g. mapaccum) "
          "are not supported.");
        work_to_onnx[o[0]] = onnx_input_name(f, i_vec[0]);
        continue;
      }

      // OP_OUTPUT connects the internal result to the declared output name (a folded reshape
      // shows up as a source/output shape mismatch)
      if (op == OP_OUTPUT) {
        casadi_assert(written_outputs.insert(o[0]).second,
          "ONNX export: nested functions with multi-segment outputs (e.g. mapaccum) "
          "are not supported.");
        MX dep = f.instruction_MX(k).dep(0);
        Sparsity out_sp = f.sparsity_out(o[0]);
        emit_output_node(func, work_to_onnx[i_vec[0]], dep.size1(), dep.size2(), out_sp,
                         onnx_output_name(f, o[0]), "out_rs_" + std::to_string(k));
        continue;
      }

      // OP_CALL - recursive function call
      if (op == OP_CALL) {
        Function called_func = f.instruction_MX(k).which_function();
        std::string out_prefix = "call_" + called_func.name() + "_" + std::to_string(k) + "_out";
        if (is_map_function(called_func)) {
          export_map(func, called_func, i_vec, o, work_to_onnx, out_prefix);
        } else if (is_reduce_map_function(called_func)) {
          export_reduce_map(func, called_func, i_vec, o, work_to_onnx, out_prefix);
        } else if (is_if_else_function(called_func)) {
          export_if(func, called_func, i_vec, o, work_to_onnx, out_prefix);
        } else {
          assert_not_control_flow(called_func);
          export_call(func, called_func, i_vec, o, work_to_onnx, out_prefix);
        }
        continue;
      }

      // All other operations use shared implementation via callback
      if (process_operation([&]() { return func->add_node(); },
                            f, op, k, i_vec, o, work_to_onnx, node_output)) {
        continue;
      }

      // Unsupported operation
      casadi_error("ONNX export: Unsupported operation code " + std::to_string(op) +
                  " in function '" + f.name() + "' at instruction " + std::to_string(k));
    }

    if (verbose_) {
      uout() << "    Created " << func->node_size() << " nodes in FunctionProto" << std::endl;
    }

    return func_guard.release();
  }

} // namespace casadi
/// \endcond
