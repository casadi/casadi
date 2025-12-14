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

    // Fill with data based on type
    // CasADi DM is always double, so we convert from various ONNX types
    onnx::TensorProto::DataType dtype =
        static_cast<onnx::TensorProto::DataType>(tensor.data_type());

    switch (dtype) {
      case onnx::TensorProto::DOUBLE: {
        // ONNX can store data in either double_data field or raw_data field
        if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
          const std::string& raw = tensor.raw_data();
          casadi_assert(raw.size() == rows * cols * sizeof(double),
                        "Raw data size mismatch for DOUBLE");
          const double* data_ptr = reinterpret_cast<const double*>(raw.data());
          for (casadi_int i = 0; i < rows * cols; ++i) {
            dm(i) = data_ptr[i];
          }
        } else if (tensor.double_data_size() > 0) {
          casadi_assert(tensor.double_data_size() == rows * cols,
                        "Tensor data size mismatch for DOUBLE");
          for (casadi_int i = 0; i < tensor.double_data_size(); ++i) {
            dm(i) = tensor.double_data(i);
          }
        } else {
          casadi_error("DOUBLE tensor has neither double_data nor raw_data");
        }
        break;
      }

      case onnx::TensorProto::FLOAT: {
        // 32-bit float -> promoted to 64-bit double
        if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
          const std::string& raw = tensor.raw_data();
          casadi_assert(raw.size() == rows * cols * sizeof(float),
                        "Raw data size mismatch for FLOAT");
          const float* data_ptr = reinterpret_cast<const float*>(raw.data());
          for (casadi_int i = 0; i < rows * cols; ++i) {
            dm(i) = static_cast<double>(data_ptr[i]);
          }
        } else if (tensor.float_data_size() > 0) {
          casadi_assert(tensor.float_data_size() == rows * cols,
                        "Tensor data size mismatch for FLOAT");
          for (casadi_int i = 0; i < tensor.float_data_size(); ++i) {
            dm(i) = static_cast<double>(tensor.float_data(i));
          }
        } else {
          casadi_error("FLOAT tensor has neither float_data nor raw_data");
        }
        break;
      }

      case onnx::TensorProto::INT32: {
        // 32-bit integer -> converted to double
        if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
          const std::string& raw = tensor.raw_data();
          casadi_assert(raw.size() == rows * cols * sizeof(int32_t),
                        "Raw data size mismatch for INT32");
          const int32_t* data_ptr = reinterpret_cast<const int32_t*>(raw.data());
          for (casadi_int i = 0; i < rows * cols; ++i) {
            dm(i) = static_cast<double>(data_ptr[i]);
          }
        } else if (tensor.int32_data_size() > 0) {
          casadi_assert(tensor.int32_data_size() == rows * cols,
                        "Tensor data size mismatch for INT32");
          for (casadi_int i = 0; i < tensor.int32_data_size(); ++i) {
            dm(i) = static_cast<double>(tensor.int32_data(i));
          }
        } else {
          casadi_error("INT32 tensor has neither int32_data nor raw_data");
        }
        break;
      }

      case onnx::TensorProto::INT64: {
        // 64-bit integer -> converted to double (may lose precision for large values)
        bool has_precision_risk = false;
        if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
          const std::string& raw = tensor.raw_data();
          casadi_assert(raw.size() == rows * cols * sizeof(int64_t),
                        "Raw data size mismatch for INT64");
          const int64_t* data_ptr = reinterpret_cast<const int64_t*>(raw.data());
          for (casadi_int i = 0; i < rows * cols; ++i) {
            int64_t val = data_ptr[i];
            if (std::abs(val) > (1LL << 53)) has_precision_risk = true;
            dm(i) = static_cast<double>(val);
          }
        } else if (tensor.int64_data_size() > 0) {
          casadi_assert(tensor.int64_data_size() == rows * cols,
                        "Tensor data size mismatch for INT64");
          for (casadi_int i = 0; i < tensor.int64_data_size(); ++i) {
            int64_t val = tensor.int64_data(i);
            if (std::abs(val) > (1LL << 53)) has_precision_risk = true;
            dm(i) = static_cast<double>(val);
          }
        } else {
          casadi_error("INT64 tensor has neither int64_data nor raw_data");
        }
        if (has_precision_risk && verbose_) {
          uout() << "WARNING: INT64 tensor contains values exceeding double precision "
                 << "(2^53), precision loss may occur" << std::endl;
        }
        break;
      }

      case onnx::TensorProto::BOOL: {
        // Boolean -> converted to 0.0/1.0
        // ONNX stores BOOL in int32_data field (0 or 1)
        if (tensor.int32_data_size() > 0) {
          casadi_assert(tensor.int32_data_size() == rows * cols,
                        "Tensor data size mismatch for BOOL");
          for (casadi_int i = 0; i < tensor.int32_data_size(); ++i) {
            dm(i) = tensor.int32_data(i) ? 1.0 : 0.0;
          }
        } else if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
          // Some ONNX implementations may use raw_data for bool (1 byte per bool)
          const std::string& raw = tensor.raw_data();
          casadi_assert(raw.size() == rows * cols,
                        "Raw data size mismatch for BOOL");
          for (casadi_int i = 0; i < rows * cols; ++i) {
            dm(i) = raw[i] ? 1.0 : 0.0;
          }
        } else {
          casadi_error("BOOL tensor has neither int32_data nor raw_data");
        }
        break;
      }

      default:
        casadi_error("Unsupported ONNX tensor data type: " + std::to_string(dtype) +
                     ". Supported types: DOUBLE(11), FLOAT(1), INT32(6), INT64(7), BOOL(9)");
    }

    return dm;
  }

} // namespace casadi
/// \endcond
