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

/// \cond INTERNAL
namespace casadi {

  void OnnxTranslator::set_dimension(const std::string& name, casadi_int dim) {
    dimension_overrides_[name] = dim;

    if (verbose_) {
      uout() << "Set dimension override: " << name << " = " << dim << std::endl;
    }
  }

  casadi_int OnnxTranslator::get_dimension(
      const onnx::TensorShapeProto& shape, int idx) const {
    // Out of range or no dimensions = scalar (1)
    if (idx >= shape.dim_size()) return 1;

    const auto& dim = shape.dim(idx);

    // Concrete dimension
    if (dim.has_dim_value()) {
      return static_cast<casadi_int>(dim.dim_value());
    }

    // Symbolic dimension - look up in overrides
    if (dim.has_dim_param()) {
      std::string param_name = dim.dim_param();

      auto it = dimension_overrides_.find(param_name);
      casadi_assert(it != dimension_overrides_.end(),
                    "Symbolic dimension '" + param_name + "' not specified. " +
                    "Call set_dimension(\"" + param_name + "\", value) before create().");

      return it->second;
    }

    // No dimension info - default to 1
    return 1;
  }

  DM OnnxTranslator::tensor_to_dm(const onnx::TensorProto& tensor) const {
    // Extract shape
    std::vector<casadi_int> dims;
    for (int i = 0; i < tensor.dims_size(); ++i) {
      dims.push_back(static_cast<casadi_int>(tensor.dims(i)));
    }

    if (verbose_) {
      uout() << "      tensor_to_dm: dims_size=" << tensor.dims_size()
             << ", double_data_size=" << tensor.double_data_size()
             << ", raw_data_size=" << tensor.raw_data().size()
             << ", data_type=" << tensor.data_type();
      if (dims.size() > 0) {
        uout() << ", dims[0]=" << dims[0];
      }
      uout() << std::endl;
    }

    // Determine matrix dimensions (assume 2D or less)
    casadi_int rows = dims.size() > 0 ? dims[0] : 1;
    casadi_int cols = dims.size() > 1 ? dims[1] : 1;

    if (verbose_) {
      uout() << "      rows=" << rows << ", cols=" << cols << std::endl;
    }

    // For higher dimensional tensors, flatten to 2D for now
    if (dims.size() > 2) {
      casadi_warning("ONNX tensor has " + std::to_string(dims.size()) +
                     " dimensions, flattening to 2D");
      cols = 1;
      for (size_t i = 1; i < dims.size(); ++i) {
        cols *= dims[i];
      }
    }

    // Create DM
    DM dm(rows, cols);

    // Fill with data (only handle DOUBLE for now)
    casadi_assert(tensor.data_type() == onnx::TensorProto::DOUBLE,
                  "Only DOUBLE tensors supported in import, got type " +
                  std::to_string(tensor.data_type()));

    // ONNX can store data in either double_data field or raw_data field
    // Check raw_data first since that's what our export uses
    if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
      // Data in raw bytes
      const std::string& raw = tensor.raw_data();
      casadi_assert(raw.size() == rows * cols * sizeof(double),
                    "Raw data size mismatch: expected " +
                    std::to_string(rows * cols * sizeof(double)) + ", got " +
                    std::to_string(raw.size()));

      const double* data_ptr = reinterpret_cast<const double*>(raw.data());
      for (casadi_int i = 0; i < rows * cols; ++i) {
        dm(i) = data_ptr[i];
      }
    } else if (tensor.double_data_size() > 0) {
      // Data in typed field
      casadi_assert(tensor.double_data_size() == rows * cols,
                    "Tensor data size mismatch: expected " +
                    std::to_string(rows * cols) + ", got " +
                    std::to_string(tensor.double_data_size()));

      for (casadi_int i = 0; i < tensor.double_data_size(); ++i) {
        dm(i) = tensor.double_data(i);
      }
    } else {
      casadi_error("Tensor has neither double_data nor raw_data");
    }

    return dm;
  }

} // namespace casadi
/// \endcond
