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


#ifndef CASADI_ONNX_TRANSLATOR_HPP
#define CASADI_ONNX_TRANSLATOR_HPP

#include <casadi/core/translator_internal.hpp>
#include <casadi/interfaces/onnx/casadi_translator_onnx_export.h>

#define ONNX_ML 1
#define ONNX_NAMESPACE onnx
#include <onnx/onnx_pb.h>

#include <map>
#include <string>

/// \cond INTERNAL
namespace casadi {

  /** \brief ONNX Translator

      Imports and exports ONNX computational graphs

      \author Joris Gillis
      \date 2025

      \identifier{onnx_translator} */
  class CASADI_EXPORT OnnxTranslator : public TranslatorInternal {
  public:
    /// Constructor
    OnnxTranslator();

    /// Destructor
    ~OnnxTranslator() override;

    /** \brief Get type name

        \identifier{onnx_translator_type} */
    std::string class_name() const override { return "OnnxTranslator";}

    /// Query plugin name
    const char* plugin_name() const override { return "onnx";}

    ///@{
    /** \brief Options

        \identifier{onnx_translator_options} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /** \brief Initialize

        \identifier{onnx_translator_init} */
    void init(const Dict& opts) override;

    /** \brief Load a graph from ONNX file

        \identifier{onnx_translator_load_file} */
    void load(const std::string& filename) override;

    /** \brief Load a CasADi Function and convert to ONNX representation

        \identifier{onnx_translator_load_function} */
    void load(const Function& f) override;

    /** \brief Set dimension for a symbolic variable

        \identifier{onnx_translator_set_dimension} */
    void set_dimension(const std::string& name, casadi_int dim) override;

    /** \brief Create a CasADi Function from the loaded ONNX graph

        \identifier{onnx_translator_create} */
    Function create(const std::string& name) override;

    /** \brief Save the loaded graph to ONNX file

        \identifier{onnx_translator_save} */
    void save(const std::string& filename) override;

    /// Creator function for plugin
    static TranslatorInternal* creator() {
      return new OnnxTranslator();
    }

    /// A documentation string
    static const std::string meta_doc;

  protected:
    /// ONNX model protocol buffer
    onnx::ModelProto model_;

    /// Dimension overrides for symbolic dimensions
    std::map<std::string, casadi_int> dimension_overrides_;

    /// Whether a model has been loaded
    bool has_model_;

  private:
    /** \brief Get dimension value from ONNX shape

        Handles both concrete dimensions and symbolic dimensions
        using dimension_overrides_ map.

        \identifier{onnx_translator_get_dimension} */
    casadi_int get_dimension(const onnx::TensorShapeProto& shape, int idx) const;

    /** \brief Convert ONNX TensorProto to CasADi DM

        Extracts shape and data from ONNX tensor.
        Currently only supports DOUBLE data type.

        \identifier{onnx_translator_tensor_to_dm} */
    DM tensor_to_dm(const onnx::TensorProto& tensor) const;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_ONNX_TRANSLATOR_HPP
