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


#ifndef CASADI_GRAPH_BUILDER_HPP
#define CASADI_GRAPH_BUILDER_HPP

#include "function.hpp"

namespace casadi {

  class GraphBuilderInternal;

  /** \brief Metadata for one graph tensor (graph input or output)

      The format-neutral GraphBuilder analogue of DaeBuilder's Variable.

      */
  struct CASADI_EXPORT Node
    : public SWIG_IF_ELSE(PrintableCommon, Printable<Node>) {
    std::string name;                     ///< Tensor name
    std::string io;                       ///< "input" or "output"
    std::string dtype;                    ///< Element type name (FLOAT, INT64, ...)
    std::vector<casadi_int> shape;        ///< Declared shape, -1 for dynamic dimensions
    std::vector<std::string> dim_params;  ///< Symbolic name per axis ("" if static)
    std::vector<double> value;            ///< Baked value (optional)
    bool baked = false;                   ///< True if a fixed value was set() for this input

    std::string type_name() const { return "Node"; }
    void disp(std::ostream& stream, bool more=false) const;
    std::string get_str(bool more=false) const;
  };

  /** \brief A mutable, format-neutral interface to a computational-graph model

      Two-stage workflow: explore and configure a model (introspection, dynamic-dimension
      binding) with GraphBuilder, then freeze it into an immutable, evaluable Function.

      \verbatim
      GraphBuilder b("model.onnx");
      b.bind_dim("batch", 4);
      Function f = b.create("net");
      \endverbatim

      \date 2026

      */
  class CASADI_EXPORT GraphBuilder
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<GraphBuilder>) {
   public:
    /// Readable name of the class
    std::string type_name() const { return "GraphBuilder"; }

    /// Default constructor
    GraphBuilder();

    /// Construct from a model file (import lifecycle; format from the file suffix)
    explicit GraphBuilder(const std::string& model_path, const Dict& opts = Dict());

    /// Construct from a CasADi Function (export lifecycle)
    explicit GraphBuilder(const Function& f, const Dict& opts = Dict());

#ifndef SWIG
    /// Construct from model data in memory
    GraphBuilder(const std::string& name, const std::vector<uint8_t>& model_data,
                 const std::string& format, const Dict& opts = Dict());
#endif // SWIG

    /// Name of the model
    const std::string& name() const;

    ///@{
    /** @name Model introspection */
    casadi_int n_in() const;
    casadi_int n_out() const;
    std::vector<std::string> name_in() const;
    std::vector<std::string> name_out() const;
    std::vector<casadi_int> input_shape(const std::string& name) const;
    std::vector<casadi_int> output_shape(const std::string& name) const;
    /// Names of the symbolic/dynamic dimensions in the model
    std::vector<std::string> dynamic_params() const;
    /// Metadata for a single tensor (input or output) by name
    Node node(const std::string& name) const;
#ifndef SWIG
    /// All tensors (inputs followed by outputs)
    std::vector<Node> nodes() const;
#endif // SWIG
    ///@}

    ///@{
    /** @name Configuration */
    /// Bind a symbolic/dynamic dimension to a concrete size
    void bind_dim(const std::string& param, casadi_int value);
    /// Pin the full shape of an input
    void bind_shape(const std::string& input_name, const std::vector<casadi_int>& shape);
    /// Bake a fixed value into an input; it is fed at create() and not exposed as a Function input
    void set(const std::string& input_name, const std::vector<double>& value);
    /// Bake a scalar value into an input
    void set(const std::string& input_name, double value);
    ///@}

    /** \brief Freeze into an evaluable Function

        \param name      Name assigned to the resulting Function
        \param name_in   Names of the inputs to expose (empty = all model inputs)
        \param name_out  Names of the outputs to expose (empty = all model outputs)
        \param opts      "symbolic" (bool, default false) and "backend" (numeric backend,
                         default "ort"); any remaining options pass through to the backend.

        */
    Function create(const std::string& name,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out,
      const Dict& opts = Dict()) const;

    /** \brief Freeze into an evaluable Function, exposing all model inputs and outputs

        \param name  Name assigned to the resulting Function
        \param opts  See the full create() overload

        */
    Function create(const std::string& name, const Dict& opts = Dict()) const;

    /** \brief Freeze into an evaluable Function, default naming */
    Function create() const { return create(name() + "_graph"); }

    /** \brief Export to an ONNX model file

        */
    void export_onnx(const std::string& filename, const Dict& opts = Dict());

#ifndef SWIG
    ///@{
    /// Access functions of the node
    GraphBuilderInternal* operator->();
    const GraphBuilderInternal* operator->() const;
    GraphBuilderInternal* get() const;
    ///@}
#endif // SWIG
  };

} // namespace casadi

#endif // CASADI_GRAPH_BUILDER_HPP
