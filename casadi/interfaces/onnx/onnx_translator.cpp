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


#include "onnx_translator.hpp"
#include <casadi/core/casadi_misc.hpp>
#include <casadi/core/calculus.hpp>
#include <casadi/core/mx.hpp>
#include <casadi/config.h>
#include <fstream>

namespace casadi {
  extern "C"
  int CASADI_TRANSLATOR_ONNX_EXPORT
  casadi_register_translator_onnx(TranslatorInternal::Plugin* plugin) {
    plugin->creator = OnnxTranslator::creator;
    plugin->name = "onnx";
    plugin->doc = OnnxTranslator::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &OnnxTranslator::options_;
    return 0;
  }

  extern "C"
  void CASADI_TRANSLATOR_ONNX_EXPORT casadi_load_translator_onnx() {
    TranslatorInternal::registerPlugin(casadi_register_translator_onnx);
  }

  OnnxTranslator::OnnxTranslator()
    : TranslatorInternal("onnx"), has_model_(false) {
  }

  OnnxTranslator::~OnnxTranslator() {
  }

  const Options OnnxTranslator::options_
  = {{&TranslatorInternal::options_},
     {{"opset_version",
       {OT_INT,
        "ONNX opset version to use when exporting (default: 13)"}}}
  };

  void OnnxTranslator::init(const Dict& opts) {
    // Call base class
    TranslatorInternal::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="opset_version") {
        // Store for later use
      }
    }
  }

  void OnnxTranslator::load(const std::string& filename) {
    // Read ONNX file
    std::ifstream input(filename, std::ios::binary);
    casadi_assert(input.good(), "Failed to open ONNX file: " + filename);

    // Parse protobuf
    casadi_assert(model_.ParseFromIstream(&input),
                  "Failed to parse ONNX model from file: " + filename);

    input.close();
    has_model_ = true;

    if (verbose_) {
      uout() << "Loaded ONNX model: " << model_.graph().name() << std::endl;
      uout() << "  Inputs: " << model_.graph().input_size() << std::endl;
      uout() << "  Outputs: " << model_.graph().output_size() << std::endl;
      uout() << "  Nodes: " << model_.graph().node_size() << std::endl;
    }
  }

  void OnnxTranslator::save(const std::string& filename) {
    casadi_assert(has_model_, "No ONNX model loaded. Call load() first.");

    // Write ONNX file
    std::ofstream output(filename, std::ios::binary);
    casadi_assert(output.good(), "Failed to open output file: " + filename);

    casadi_assert(model_.SerializeToOstream(&output),
                  "Failed to serialize ONNX model to file: " + filename);

    output.close();

    if (verbose_) {
      uout() << "Saved ONNX model to: " << filename << std::endl;
    }
  }

} // namespace casadi
