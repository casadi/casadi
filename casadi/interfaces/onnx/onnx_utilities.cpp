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

  Sparsity Onnx::input_pattern(const onnx::GraphProto& graph, const std::string& name) const {
    for (int i = 0; i < graph.sparse_initializer_size(); ++i) {
      if (graph.sparse_initializer(i).values().name() == name) {
        return sparse_tensor_to_dm(graph.sparse_initializer(i)).sparsity();
      }
    }
    return Sparsity();  // null -> no overlay
  }

  casadi_int get_int_attribute(const onnx::NodeProto& node, const std::string& name,
                               casadi_int default_value) {
    for (int a = 0; a < node.attribute_size(); ++a) {
      if (node.attribute(a).name() == name) return node.attribute(a).i();
    }
    return default_value;
  }

  double get_float_attribute(const onnx::NodeProto& node, const std::string& name,
                             double default_value) {
    for (int a = 0; a < node.attribute_size(); ++a) {
      if (node.attribute(a).name() == name) return node.attribute(a).f();
    }
    return default_value;
  }

  std::string get_string_attribute(const onnx::NodeProto& node, const std::string& name) {
    for (int a = 0; a < node.attribute_size(); ++a) {
      if (node.attribute(a).name() == name) return node.attribute(a).s();
    }
    return "";
  }

  const onnx::GraphProto* get_graph_attribute(const onnx::NodeProto& node,
                                              const std::string& name) {
    for (int a = 0; a < node.attribute_size(); ++a) {
      if (node.attribute(a).name() == name) return &node.attribute(a).g();
    }
    return nullptr;
  }

  void Onnx::set_dimension(const std::string& name, casadi_int dim) {
    dimension_overrides_[name] = dim;

    if (verbose_) {
      uout() << "Set dimension override: " << name << " = " << dim << std::endl;
    }
  }

  casadi_int Onnx::get_dimension(
      const onnx::TensorShapeProto& shape, int idx) const {
    if (idx >= shape.dim_size()) return 1;  // out of range / no dims -> scalar
    const auto& dim = shape.dim(idx);

    if (dim.has_dim_value()) return static_cast<casadi_int>(dim.dim_value());

    if (dim.has_dim_param()) {  // symbolic dimension -> resolve via overrides
      std::string param_name = dim.dim_param();
      auto it = dimension_overrides_.find(param_name);
      casadi_assert(it != dimension_overrides_.end(),
                    "Symbolic dimension '" + param_name + "' not specified. " +
                    "Call set_dimension(\"" + param_name + "\", value) before create().");
      return it->second;
    }

    return 1;  // no dimension info
  }

  // Decode n raw_data values of type T into dm; return false if there's no raw_data
  template<typename T>
  static bool load_raw(DM& dm, casadi_int n, const onnx::TensorProto& tensor,
                       const char* type_name) {
    if (!(tensor.has_raw_data() && tensor.raw_data().size() > 0)) return false;
    const std::string& raw = tensor.raw_data();
    casadi_assert(raw.size() == static_cast<size_t>(n) * sizeof(T),
                  std::string("Raw data size mismatch for ") + type_name);
    const T* p = reinterpret_cast<const T*>(raw.data());
    for (casadi_int i = 0; i < n; ++i) dm(i) = static_cast<double>(p[i]);
    return true;
  }

  DM Onnx::tensor_to_dm(const onnx::TensorProto& tensor) const {
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

    // 2-D or less; higher-rank tensors are flattened to 2-D (columns absorb the extra axes).
    casadi_int rows = dims.size() > 0 ? dims[0] : 1;
    casadi_int cols = dims.size() > 1 ? dims[1] : 1;
    if (dims.size() > 2) {
      casadi_warning("ONNX tensor has " + std::to_string(dims.size()) +
                     " dimensions, flattening to 2D");
      cols = 1;
      for (size_t i = 1; i < dims.size(); ++i) {
        cols *= dims[i];
      }
    }

    // Transpose-rep: a 2-D ONNX tensor (c,r) holds a CasADi (r,c) value's column-major bytes, so
    // reading it back is just swapping the declared dims (the byte/fill order already matches).
    if (dims.size() == 2) std::swap(rows, cols);

    DM dm(rows, cols);  // DM is always double; ONNX integer/bool/float types are converted below

    onnx::TensorProto::DataType dtype =
        static_cast<onnx::TensorProto::DataType>(tensor.data_type());

    casadi_int n = rows * cols;
    switch (dtype) {
      case onnx::TensorProto::DOUBLE: {
        if (!load_raw<double>(dm, n, tensor, "DOUBLE")) {
          casadi_assert(tensor.double_data_size() == n, "Tensor data size mismatch for DOUBLE");
          for (casadi_int i = 0; i < n; ++i) dm(i) = tensor.double_data(i);
        }
        break;
      }

      case onnx::TensorProto::FLOAT: {  // 32-bit float -> promoted to double
        if (!load_raw<float>(dm, n, tensor, "FLOAT")) {
          casadi_assert(tensor.float_data_size() == n, "Tensor data size mismatch for FLOAT");
          for (casadi_int i = 0; i < n; ++i) dm(i) = static_cast<double>(tensor.float_data(i));
        }
        break;
      }

      case onnx::TensorProto::INT32: {  // 32-bit integer -> double
        if (!load_raw<int32_t>(dm, n, tensor, "INT32")) {
          casadi_assert(tensor.int32_data_size() == n, "Tensor data size mismatch for INT32");
          for (casadi_int i = 0; i < n; ++i) dm(i) = static_cast<double>(tensor.int32_data(i));
        }
        break;
      }

      case onnx::TensorProto::INT64: {  // 64-bit integer -> double
        if (!load_raw<int64_t>(dm, n, tensor, "INT64")) {
          casadi_assert(tensor.int64_data_size() == n, "Tensor data size mismatch for INT64");
          for (casadi_int i = 0; i < n; ++i) dm(i) = static_cast<double>(tensor.int64_data(i));
        }
        break;
      }

      case onnx::TensorProto::BOOL: {
        // Boolean -> 0.0/1.0. ONNX stores BOOL in int32_data, or 1 byte per bool in raw_data.
        if (tensor.int32_data_size() > 0) {
          casadi_assert(tensor.int32_data_size() == n, "Tensor data size mismatch for BOOL");
          for (casadi_int i = 0; i < n; ++i) dm(i) = tensor.int32_data(i) ? 1.0 : 0.0;
        } else if (tensor.has_raw_data() && tensor.raw_data().size() > 0) {
          const std::string& raw = tensor.raw_data();
          casadi_assert(raw.size() == static_cast<size_t>(n), "Raw data size mismatch for BOOL");
          for (casadi_int i = 0; i < n; ++i) dm(i) = raw[i] ? 1.0 : 0.0;
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

  DM Onnx::sparse_tensor_to_dm(const onnx::SparseTensorProto& st) const {
    // Dense shape (2-D; flatten any higher dims into the columns, like tensor_to_dm)
    casadi_int nrow = st.dims_size() > 0 ? static_cast<casadi_int>(st.dims(0)) : 1;
    casadi_int ncol = st.dims_size() > 1 ? static_cast<casadi_int>(st.dims(1)) : 1;
    for (int i = 2; i < st.dims_size(); ++i) ncol *= static_cast<casadi_int>(st.dims(i));

    // Values: a 1-D [NNZ] tensor -> reuse the dense reader, then take its entries
    std::vector<double> values = tensor_to_dm(st.values()).nonzeros();
    casadi_int nnz = static_cast<casadi_int>(values.size());

    // Indices: INT64, either [NNZ,2] explicit (row,col) or [NNZ] row-major linearized
    const onnx::TensorProto& it = st.indices();
    std::vector<int64_t> flat;
    if (it.int64_data_size() > 0) {
      flat.assign(it.int64_data().begin(), it.int64_data().end());
    } else if (it.has_raw_data() && !it.raw_data().empty()) {
      const std::string& raw = it.raw_data();
      const int64_t* p = reinterpret_cast<const int64_t*>(raw.data());
      flat.assign(p, p + raw.size() / sizeof(int64_t));
    }

    const bool coord = it.dims_size() == 2 && it.dims(1) == 2;  // [NNZ,2] coordinate form
    casadi_assert(static_cast<casadi_int>(flat.size()) == (coord ? 2 * nnz : nnz),
                  "Sparse tensor indices/values size mismatch");

    std::vector<casadi_int> rows(nnz), cols(nnz);
    for (casadi_int k = 0; k < nnz; ++k) {
      if (coord) {
        rows[k] = static_cast<casadi_int>(flat[2 * k]);
        cols[k] = static_cast<casadi_int>(flat[2 * k + 1]);
      } else {
        casadi_int lin = static_cast<casadi_int>(flat[k]);  // row-major linearized
        rows[k] = lin / ncol;
        cols[k] = lin % ncol;
      }
    }
    // Transpose-rep: the stored COO is V^T (dims (c,r), coords (col,row)); build it then
    // transpose to recover the CasADi value V. DM::triplet is order-agnostic.
    return DM::triplet(rows, cols, values, nrow, ncol).T();
  }

} // namespace casadi
/// \endcond
