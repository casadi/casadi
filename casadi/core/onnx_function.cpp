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

#include "onnx_function_impl.hpp"
#include "graph_builder_internal.hpp"
#include "casadi_misc.hpp"

namespace casadi {

  bool has_onnxbackend(const std::string& solver) {
    return OnnxFunction::has_plugin(solver);
  }

  void load_onnxbackend(const std::string& solver) {
    OnnxFunction::load_plugin(solver);
  }

  std::vector<std::string> onnxbackend_solvers() {
    std::vector<std::string> ret;
    for (auto&& s : OnnxFunction::solvers_) ret.push_back(s.first);
    return ret;
  }

  std::string onnxbackend_doc(const std::string& solver) {
    return OnnxFunction::getPlugin(solver).doc;
  }

  std::map<std::string, OnnxFunction::Plugin> OnnxFunction::solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::mutex OnnxFunction::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

  const std::string OnnxFunction::infix_ = "onnx";

  std::string OnnxFunction::meta_doc = "";

  const Options OnnxFunction::options_
  = {{&FunctionInternal::options_},
     {{"provider",
       {OT_STRING, "Execution provider for the ONNX runtime backend"}},
      {"dim_bindings",
       {OT_DICT, "Sizes for symbolic/dynamic tensor dimensions (name -> size)"}},
      {"input_shapes",
       {OT_DICT, "Explicit shapes for inputs (name -> shape), overriding the model's"}},
      {"input_values",
       {OT_DICT, "Baked-in input values (name -> value); these inputs are not exposed"}},
      {"fwd_dim",
       {OT_STRING, "Symbolic dimension naming the forward seed count [nfwd]"}},
      {"adj_dim",
       {OT_STRING, "Symbolic dimension naming the adjoint seed count [nadj]"}}
     }
    };

  std::string onnx_dtype_name(casadi_int t) {
    switch (t) {
      case 1: return "FLOAT";       case 2: return "UINT8";    case 3: return "INT8";
      case 4: return "UINT16";      case 5: return "INT16";    case 6: return "INT32";
      case 7: return "INT64";       case 8: return "STRING";   case 9: return "BOOL";
      case 10: return "FLOAT16";    case 11: return "DOUBLE";  case 12: return "UINT32";
      case 13: return "UINT64";     case 16: return "BFLOAT16";
      default: return "TYPE" + str(t);
    }
  }

  casadi_int onnx_dtype_enum(const std::string& name) {
    static const std::map<std::string, casadi_int> m = {
      {"FLOAT", 1},  {"UINT8", 2},  {"INT8", 3},   {"UINT16", 4},   {"INT16", 5},
      {"INT32", 6},  {"INT64", 7},  {"STRING", 8}, {"BOOL", 9},     {"FLOAT16", 10},
      {"DOUBLE", 11}, {"UINT32", 12}, {"UINT64", 13}, {"BFLOAT16", 16}};
    auto it = m.find(name);
    return it != m.end() ? it->second : 0;
  }

  // Select a subset of tensors by name, preserving the requested order (empty request => all)
  static std::vector<OnnxTensorInfo> select_tensors(const std::vector<OnnxTensorInfo>& all,
                                                    const std::vector<std::string>& req) {
    if (req.empty()) return all;
    std::vector<OnnxTensorInfo> sel;
    for (const std::string& name : req) {
      bool found = false;
      for (const OnnxTensorInfo& t : all)
        if (t.name == name) { sel.push_back(t); found = true; break; }
      casadi_assert(found, "ONNX tensor '" + name + "' not found in model");
    }
    return sel;
  }

  OnnxFunction::OnnxFunction(const std::string& name, const GraphBuilderInternal* gb,
                             const std::vector<std::string>& inputs,
                             const std::vector<std::string>& outputs)
    : FunctionInternal(name) {
    // Freeze a snapshot: the builder may mutate afterwards without affecting this function.
    // Shapes are resolved by the builder (dim bindings + input-shape overrides applied).
    model_data_ = gb->model_data_;
    input_values_ = gb->input_values_;
    std::vector<OnnxTensorInfo> all_out;
    for (const Node& n : gb->node_list()) {
      OnnxTensorInfo t;
      t.name = n.name;
      t.elem_type = onnx_dtype_enum(n.dtype);
      t.shape = gb->resolved_shape(n);
      t.numel = 1;
      for (casadi_int d : t.shape) t.numel *= d;
      if (n.io == "input") {
        all_in_.push_back(t);
        model_inputs_.insert(t.name);
      } else {
        all_out.push_back(t);
        model_outputs_.insert(t.name);
      }
    }
    in_ = select_tensors(all_in_, inputs);    // exposed CasADi inputs (selection)
    out_ = select_tensors(all_out, outputs);  // exposed CasADi outputs (selection)
  }

  OnnxFunction::~OnnxFunction() {
    clear_mem();
  }

  // Pack/unpack a vector<OnnxTensorInfo> as parallel arrays
  static void pack_tensors(SerializingStream& s, const std::string& d,
                           const std::vector<OnnxTensorInfo>& v) {
    std::vector<std::string> names;
    std::vector<std::vector<casadi_int>> shapes;
    std::vector<casadi_int> elem_types, numels;
    for (const OnnxTensorInfo& t : v) {
      names.push_back(t.name);
      shapes.push_back(t.shape);
      elem_types.push_back(t.elem_type);
      numels.push_back(t.numel);
    }
    s.pack(d + "::names", names);
    s.pack(d + "::shapes", shapes);
    s.pack(d + "::elem_types", elem_types);
    s.pack(d + "::numels", numels);
  }
  static void unpack_tensors(DeserializingStream& s, const std::string& d,
                             std::vector<OnnxTensorInfo>& v) {
    std::vector<std::string> names;
    std::vector<std::vector<casadi_int>> shapes;
    std::vector<casadi_int> elem_types, numels;
    s.unpack(d + "::names", names);
    s.unpack(d + "::shapes", shapes);
    s.unpack(d + "::elem_types", elem_types);
    s.unpack(d + "::numels", numels);
    v.clear();
    for (size_t k = 0; k < names.size(); ++k)
      v.push_back(OnnxTensorInfo{names[k], shapes[k], elem_types[k], numels[k]});
  }

  void OnnxFunction::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    PluginInterface<OnnxFunction>::serialize_type(s);
  }

  void OnnxFunction::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);
    s.version("OnnxFunction", 1);
    s.pack("OnnxFunction::model_data", std::string(model_data_.begin(), model_data_.end()));
    pack_tensors(s, "OnnxFunction::in", in_);
    pack_tensors(s, "OnnxFunction::out", out_);
    pack_tensors(s, "OnnxFunction::all_in", all_in_);
    s.pack("OnnxFunction::in_src", in_src_);
    s.pack("OnnxFunction::in_val", in_val_);
    s.pack("OnnxFunction::model_inputs",
           std::vector<std::string>(model_inputs_.begin(), model_inputs_.end()));
    s.pack("OnnxFunction::model_outputs",
           std::vector<std::string>(model_outputs_.begin(), model_outputs_.end()));
    s.pack("OnnxFunction::fwd_dim", fwd_dim_);
    s.pack("OnnxFunction::adj_dim", adj_dim_);
    s.pack("OnnxFunction::input_values", input_values_);
  }

  OnnxFunction::OnnxFunction(DeserializingStream& s) : FunctionInternal(s) {
    s.version("OnnxFunction", 1);
    std::string bytes;
    s.unpack("OnnxFunction::model_data", bytes);
    model_data_.assign(bytes.begin(), bytes.end());
    unpack_tensors(s, "OnnxFunction::in", in_);
    unpack_tensors(s, "OnnxFunction::out", out_);
    unpack_tensors(s, "OnnxFunction::all_in", all_in_);
    s.unpack("OnnxFunction::in_src", in_src_);
    s.unpack("OnnxFunction::in_val", in_val_);
    std::vector<std::string> mi, mo;
    s.unpack("OnnxFunction::model_inputs", mi);
    model_inputs_ = std::set<std::string>(mi.begin(), mi.end());
    s.unpack("OnnxFunction::model_outputs", mo);
    model_outputs_ = std::set<std::string>(mo.begin(), mo.end());
    s.unpack("OnnxFunction::fwd_dim", fwd_dim_);
    s.unpack("OnnxFunction::adj_dim", adj_dim_);
    s.unpack("OnnxFunction::input_values", input_values_);
  }

  ProtoFunction* OnnxFunction::deserialize(DeserializingStream& s) {
    return PluginInterface<OnnxFunction>::deserialize(s);
  }

  void OnnxFunction::init(const Dict& opts) {
    for (auto&& op : opts) {
      if (op.first == "fwd_dim") fwd_dim_ = op.second.to_string();
      else if (op.first == "adj_dim") adj_dim_ = op.second.to_string();
    }
    // Baked inputs are fed a fixed value, not exposed as Function inputs
    std::vector<OnnxTensorInfo> exposed;
    for (const OnnxTensorInfo& t : in_) if (!input_values_.count(t.name)) exposed.push_back(t);
    in_ = exposed;
    build_io_map();  // baked-feed map over all model inputs
    FunctionInternal::init(opts);
  }

  void OnnxFunction::build_io_map() {
    // in_src: exposed-arg index (>=0), -2 baked value, or -1 unwired (fed a default later)
    in_src_.clear();
    in_val_.clear();
    for (const OnnxTensorInfo& t : all_in_) {
      casadi_int src = -1;
      for (casadi_int j = 0; j < static_cast<casadi_int>(in_.size()); ++j)
        if (in_[j].name == t.name) { src = j; break; }
      auto bv = input_values_.find(t.name);
      if (src < 0 && bv != input_values_.end()) {
        casadi_assert(static_cast<casadi_int>(bv->second.size()) == t.numel,
                      "Baked value for '" + t.name + "' has " + str(bv->second.size())
                      + " elements, expected " + str(t.numel));
        src = -2;
        for (double v : bv->second) in_val_.push_back(v);
      } else {
        in_val_.insert(in_val_.end(), t.numel, 0.0);  // placeholder block keeps voff aligned
      }
      in_src_.push_back(src);
    }
  }

  Sparsity OnnxFunction::tensor_sparsity(const std::vector<casadi_int>& shape) {
    // CasADi is 2-D: rank-0/1/2 map directly, higher ranks flatten to a column vector
    if (shape.empty()) return Sparsity::dense(1, 1);
    if (shape.size() == 1) return Sparsity::dense(shape[0], 1);
    if (shape.size() == 2) return Sparsity::dense(shape[0], shape[1]);
    casadi_int numel = 1;
    for (casadi_int d : shape) numel *= d;
    return Sparsity::dense(numel, 1);
  }

  Function OnnxFunction::wrap_derivative(const std::string& name,
      const std::vector<std::string>& inames, const std::vector<std::string>& onames,
      const std::vector<Sparsity>& in_sp, const std::vector<Sparsity>& out_sp,
      const Dict& dim_bind, const Dict& opts) const {
    // Re-create from the same model using only the derivative tensors it actually has
    std::vector<std::string> cin, con;
    for (const std::string& nm : inames) if (model_inputs_.count(nm)) cin.push_back(nm);
    for (const std::string& nm : onames) if (model_outputs_.count(nm)) con.push_back(nm);
    Dict o;
    if (!dim_bind.empty()) o["dim_bindings"] = dim_bind;
    Function g = from_model_data(plugin_name(), name + "_core", model_data_, cin, con, o);
    // Present CasADi's full derivative signature: feed present inputs, zero-fill absent outputs
    std::map<std::string, MX> m;
    std::vector<MX> args(inames.size());
    for (size_t i = 0; i < inames.size(); ++i) {
      args[i] = MX::sym(inames[i], in_sp[i]);
      m[inames[i]] = args[i];
    }
    std::vector<MX> gin;
    for (const std::string& nm : cin) gin.push_back(m[nm]);
    std::vector<MX> gout = g(gin);
    std::map<std::string, MX> mo;
    for (size_t k = 0; k < con.size(); ++k) mo[con[k]] = gout[k];
    std::vector<MX> outs(onames.size());
    for (size_t j = 0; j < onames.size(); ++j) {
      auto it = mo.find(onames[j]);
      outs[j] = it != mo.end() ? it->second : MX::zeros(out_sp[j]);
    }
    Dict wopts;
    auto it = opts.find("derivative_of");
    if (it != opts.end()) wopts["derivative_of"] = it->second;
    return Function(name, args, outs, inames, onames, wopts);
  }

  bool OnnxFunction::has_forward(casadi_int nfwd) const {
    // Need a fwd_<tensor> for every DIFFERENTIABLE input (seed) and output (sensitivity)
    if (in_.empty() || out_.empty()) return false;
    std::string pref = diff_prefix("fwd");
    bool any_in = false, any_out = false;
    for (size_t i = 0; i < in_.size(); ++i) {
      if (!diff_in(i)) continue;
      any_in = true;
      if (!model_inputs_.count(pref + in_[i].name)) return false;
    }
    for (size_t j = 0; j < out_.size(); ++j) {
      if (!diff_out(j)) continue;
      any_out = true;
      if (!model_outputs_.count(pref + out_[j].name)) return false;
    }
    return any_in && any_out;
  }

  Function OnnxFunction::get_forward(casadi_int nfwd, const std::string& name,
                                     const std::vector<std::string>& inames,
                                     const std::vector<std::string>& onames,
                                     const Dict& opts) const {
    std::vector<Sparsity> isp, osp;
    for (const OnnxTensorInfo& x : in_) isp.push_back(tensor_sparsity(x.shape));
    for (const OnnxTensorInfo& y : out_) isp.push_back(tensor_sparsity(y.shape));
    for (const OnnxTensorInfo& x : in_) {
      Sparsity s = tensor_sparsity(x.shape);
      isp.push_back(Sparsity::dense(s.size1(), nfwd * s.size2()));
    }
    for (const OnnxTensorInfo& y : out_) {
      Sparsity s = tensor_sparsity(y.shape);
      osp.push_back(Sparsity::dense(s.size1(), nfwd * s.size2()));
    }
    Dict db; db[fwd_dim_] = nfwd;
    return wrap_derivative(name, inames, onames, isp, osp, db, opts);
  }

  bool OnnxFunction::has_reverse(casadi_int nadj) const {
    // Need an adj_<tensor> output for every differentiable input and seed input per output
    if (in_.empty() || out_.empty()) return false;
    std::string pref = diff_prefix("adj");
    bool any_in = false, any_out = false;
    for (size_t i = 0; i < in_.size(); ++i) {
      if (!diff_in(i)) continue;
      any_in = true;
      if (!model_outputs_.count(pref + in_[i].name)) return false;
    }
    for (size_t j = 0; j < out_.size(); ++j) {
      if (!diff_out(j)) continue;
      any_out = true;
      if (!model_inputs_.count(pref + out_[j].name)) return false;
    }
    return any_in && any_out;
  }

  Function OnnxFunction::get_reverse(casadi_int nadj, const std::string& name,
                                     const std::vector<std::string>& inames,
                                     const std::vector<std::string>& onames,
                                     const Dict& opts) const {
    std::vector<Sparsity> isp, osp;
    for (const OnnxTensorInfo& x : in_) isp.push_back(tensor_sparsity(x.shape));
    for (const OnnxTensorInfo& y : out_) isp.push_back(tensor_sparsity(y.shape));
    for (const OnnxTensorInfo& y : out_) {
      Sparsity s = tensor_sparsity(y.shape);
      isp.push_back(Sparsity::dense(s.size1(), nadj * s.size2()));
    }
    for (const OnnxTensorInfo& x : in_) {
      Sparsity s = tensor_sparsity(x.shape);
      osp.push_back(Sparsity::dense(s.size1(), nadj * s.size2()));
    }
    Dict db; db[adj_dim_] = nadj;
    return wrap_derivative(name, inames, onames, isp, osp, db, opts);
  }

  bool OnnxFunction::has_jacobian() const {
    // Need a jac_<out>_<in> output for every differentiable output/input pair
    if (in_.empty() || out_.empty()) return false;
    bool any = false;
    for (size_t j = 0; j < out_.size(); ++j) {
      if (!diff_out(j)) continue;
      for (size_t i = 0; i < in_.size(); ++i) {
        if (!diff_in(i)) continue;
        any = true;
        if (!model_outputs_.count("jac_" + out_[j].name + "_" + in_[i].name)) return false;
      }
    }
    return any;
  }

  Function OnnxFunction::get_jacobian(const std::string& name,
                                      const std::vector<std::string>& inames,
                                      const std::vector<std::string>& onames,
                                      const Dict& opts) const {
    std::vector<Sparsity> isp, osp;
    for (const OnnxTensorInfo& x : in_) isp.push_back(tensor_sparsity(x.shape));
    for (const OnnxTensorInfo& y : out_) isp.push_back(tensor_sparsity(y.shape));
    for (const OnnxTensorInfo& y : out_) {
      Sparsity so = tensor_sparsity(y.shape);
      for (const OnnxTensorInfo& x : in_) {
        Sparsity si = tensor_sparsity(x.shape);
        osp.push_back(Sparsity::dense(so.size1() * so.size2(), si.size1() * si.size2()));
      }
    }
    return wrap_derivative(name, inames, onames, isp, osp, Dict(), opts);
  }

  Function OnnxFunction::create(const std::string& solver, const std::string& name,
                                const GraphBuilderInternal* gb,
                                const std::vector<std::string>& inputs,
                                const std::vector<std::string>& outputs,
                                const Dict& opts) {
    return Function::create(getPlugin(solver).creator(name, gb, inputs, outputs, opts), opts);
  }

  Function OnnxFunction::from_model_data(const std::string& solver, const std::string& name,
                                         const std::vector<uint8_t>& model_data,
                                         const std::vector<std::string>& inputs,
                                         const std::vector<std::string>& outputs,
                                         const Dict& opts) {
    // Stage the model + configuration in a GraphBuilder, then freeze an OnnxFunction snapshot
    GraphBuilder b(name, model_data, "onnx");
    GraphBuilderInternal* bi = b.get();
    Dict fopts;
    for (auto&& op : opts) {
      if (op.first == "dim_bindings") {
        for (auto&& d : static_cast<Dict>(op.second)) bi->dim_bindings_[d.first] = d.second;
      } else if (op.first == "input_shapes") {
        for (auto&& d : static_cast<Dict>(op.second))
          bi->input_shapes_[d.first] = d.second.to_int_vector();
      } else if (op.first == "input_values") {
        for (auto&& d : static_cast<Dict>(op.second))
          bi->input_values_[d.first] = d.second.to_double_vector();
      } else {
        fopts[op.first] = op.second;
      }
    }
    return OnnxFunction::create(solver, name, bi, inputs, outputs, fopts);
  }

} // namespace casadi
