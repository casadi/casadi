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

#include "ort_interface.hpp"

namespace casadi {

  extern "C"
  int CASADI_ONNX_ORT_EXPORT
  casadi_register_onnx_ort(OnnxFunction::Plugin* plugin) {
    plugin->creator = OnnxRuntimeInterface::creator;
    plugin->name = "ort";
    plugin->doc = OnnxRuntimeInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &OnnxRuntimeInterface::options_;
    plugin->deserialize = &OnnxRuntimeInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_ONNX_ORT_EXPORT casadi_load_onnx_ort() {
    OnnxFunction::registerPlugin(casadi_register_onnx_ort);
  }

  const std::string OnnxRuntimeInterface::meta_doc =
    "Black-box ONNX model evaluation through Microsoft's ONNX Runtime.\n";

  const Options OnnxRuntimeInterface::options_
  = {{&OnnxFunction::options_},
     {{"provider",
       {OT_STRING, "Execution provider ('CPU', 'CUDA') [CPU]"}}
     }
    };

  OnnxRuntimeInterface::OnnxRuntimeInterface(const std::string& name,
                                             const GraphBuilderInternal* gb,
                                             const std::vector<std::string>& inputs,
                                             const std::vector<std::string>& outputs)
    : OnnxFunction(name, gb, inputs, outputs),
      ort_api_(OrtGetApiBase()->GetApi(ORT_API_VERSION)), provider_("CPU") {
    casadi_assert(ort_api_ != nullptr, "Failed to obtain the ONNX Runtime API");
  }

  void OnnxRuntimeInterface::ort_check(OrtStatus* status, const std::string& what) const {
    if (!status) return;
    std::string msg = ort_api_->GetErrorMessage(status);
    ort_api_->ReleaseStatus(status);
    casadi_error("ONNX Runtime error (" + what + "): " + msg);
  }

  void OnnxRuntimeInterface::build_prob() {
    // Metadata (all_in_, in_, out_, in_src_, in_val_, resolved shapes) is frozen by the base
    in_names_c_.clear(); out_names_c_.clear();
    in_elem_.clear(); out_elem_.clear(); in_ndim_.clear(); out_ndim_.clear();
    in_dims_.clear(); out_dims_.clear(); in_numel_.clear(); out_numel_.clear();
    // Inputs: ALL model inputs (the base in_src_/in_val_ feed map covers them)
    for (const OnnxTensorInfo& t : all_in_) {
      in_names_c_.push_back(t.name.c_str());
      in_elem_.push_back(t.elem_type);
      in_ndim_.push_back(static_cast<casadi_int>(t.shape.size()));
      for (casadi_int d : t.shape) in_dims_.push_back(d);
      in_numel_.push_back(t.numel);
    }
    // Outputs: only the exposed selection (ORT runs just the needed subgraph)
    for (const OnnxTensorInfo& t : out_) {
      out_names_c_.push_back(t.name.c_str());
      out_elem_.push_back(t.elem_type);
      out_ndim_.push_back(static_cast<casadi_int>(t.shape.size()));
      for (casadi_int d : t.shape) out_dims_.push_back(d);
      out_numel_.push_back(t.numel);
    }
    prob_.n_in = all_in_.size();
    prob_.n_out = out_.size();
    prob_.input_names = in_names_c_.data();
    prob_.output_names = out_names_c_.data();
    prob_.in_src = in_src_.data();
    prob_.in_val = in_val_.data();
    prob_.in_elem_type = in_elem_.data();
    prob_.out_elem_type = out_elem_.data();
    prob_.in_ndim = in_ndim_.data();
    prob_.out_ndim = out_ndim_.data();
    prob_.in_dims = in_dims_.data();
    prob_.out_dims = out_dims_.data();
    prob_.in_numel = in_numel_.data();
    prob_.out_numel = out_numel_.data();
    prob_.model_data = model_data_.data();
    prob_.model_size = static_cast<casadi_int>(model_data_.size());
  }

  void OnnxRuntimeInterface::init(const Dict& opts) {
    OnnxFunction::init(opts);

    for (auto&& op : opts) {
      if (op.first == "provider") provider_ = op.second.to_string();
    }

    build_prob();
  }

  int OnnxRuntimeInterface::init_mem(void* mem) const {
    if (OnnxFunction::init_mem(mem)) return 1;
    auto* m = static_cast<OnnxRuntimeMemory*>(mem);
    casadi_assert(casadi_onnxruntime_init(&m->d, &prob_) == 0,
                  "Failed to create ONNX Runtime session for '" + name_ + "'");
    return 0;
  }

  void OnnxRuntimeInterface::free_mem(void* mem) const {
    auto* m = static_cast<OnnxRuntimeMemory*>(mem);
    casadi_onnxruntime_data& d = m->d;
    // Release the reusable scaffolding (built lazily by the runtime's prepare step) ...
    if (d.inv) {
      for (casadi_int i = 0; i < prob_.n_in; ++i) if (d.inv[i]) ort_api_->ReleaseValue(d.inv[i]);
      free(d.inv);
    }
    if (d.buf) {
      for (casadi_int i = 0; i < prob_.n_in; ++i) if (d.buf[i]) free(d.buf[i]);
      free(d.buf);
    }
    free(d.outv);
    free(d.row);
    // ... then the session handles created by init_mem
    if (d.mem) ort_api_->ReleaseMemoryInfo(d.mem);
    if (d.session) ort_api_->ReleaseSession(d.session);
    if (d.env) ort_api_->ReleaseEnv(d.env);
    delete m;
  }

  OnnxRuntimeInterface::~OnnxRuntimeInterface() {
    // Free the memory pool while this is still the dynamic type, so the virtual
    // free_mem() above runs (the base ~OnnxFunction::clear_mem() is too late).
    clear_mem();
  }

  void OnnxRuntimeInterface::serialize_body(SerializingStream &s) const {
    OnnxFunction::serialize_body(s);
    s.version("OnnxRuntimeInterface", 1);
    s.pack("OnnxRuntimeInterface::provider", provider_);
  }

  OnnxRuntimeInterface::OnnxRuntimeInterface(DeserializingStream& s)
    : OnnxFunction(s),
      ort_api_(OrtGetApiBase()->GetApi(ORT_API_VERSION)), provider_("CPU") {
    casadi_assert(ort_api_ != nullptr, "Failed to obtain the ONNX Runtime API");
    s.version("OnnxRuntimeInterface", 1);
    s.unpack("OnnxRuntimeInterface::provider", provider_);
    build_prob();
  }

  int OnnxRuntimeInterface::eval(const double** arg, double** res,
                                 casadi_int* iw, double* w, void* mem) const {
    auto* m = static_cast<OnnxRuntimeMemory*>(mem);
    return casadi_onnxruntime_solve(&m->d, &prob_, arg, res);
  }

  void OnnxRuntimeInterface::codegen_declarations(CodeGenerator& g) const {
    g.add_include("ort_runtime.h", false);
  }

  void OnnxRuntimeInterface::codegen_body(CodeGenerator& g) const {
    // Embed the model and its metadata (all inputs + the in_src map), then call the runtime
    std::vector<std::string> in_nm, out_nm;
    for (const OnnxTensorInfo& t : all_in_) in_nm.push_back(t.name);
    for (const OnnxTensorInfo& t : out_) out_nm.push_back(t.name);
    std::string model = g.constant(std::vector<char>(model_data_.begin(), model_data_.end()));
    std::string in_names = g.constant(in_nm), out_names = g.constant(out_nm);
    std::string in_src = g.constant(in_src_);
    std::string in_val = g.constant(in_val_.empty() ? std::vector<double>{0.0} : in_val_);
    std::string in_elem = g.constant(in_elem_), out_elem = g.constant(out_elem_);
    std::string in_ndim = g.constant(in_ndim_), out_ndim = g.constant(out_ndim_);
    std::string in_dims = g.constant(in_dims_.empty() ? std::vector<casadi_int>{0} : in_dims_);
    std::string out_dims = g.constant(out_dims_.empty() ? std::vector<casadi_int>{0} : out_dims_);
    std::string in_numel = g.constant(in_numel_), out_numel = g.constant(out_numel_);

    g << "static struct casadi_onnxruntime_prob prob = {"
      << all_in_.size() << ", " << out_.size() << ", "
      << in_names << ", " << out_names << ", "
      << in_src << ", " << in_val << ", "
      << in_elem << ", " << out_elem << ", "
      << in_ndim << ", " << out_ndim << ", "
      << in_dims << ", " << out_dims << ", "
      << in_numel << ", " << out_numel << ", "
      << "(const unsigned char*)" << model << ", " << model_data_.size() << "};\n";
    g << "static struct casadi_onnxruntime_data data = {0};\n";
    g << "if (casadi_onnxruntime_init(&data, &prob)) return 1;\n";
    g << "return casadi_onnxruntime_solve(&data, &prob, arg, res);\n";
  }

} // namespace casadi
