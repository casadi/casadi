/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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


#ifndef CASADI_MATRIX_IMPL_HPP
#define CASADI_MATRIX_IMPL_HPP

// The declaration of the class is in a separate file
#include "matrix.hpp"
#include "matrix_tools.hpp"
#include "sparsity_tools.hpp"
#include "../std_vector_tools.hpp"

/// \cond INTERNAL

namespace casadi {
  // Implementations

  template<typename DataType>
  const DataType& Matrix<DataType>::elem(int rr, int cc) const {
    int ind = sparsity().getNZ(rr, cc);
    if (ind==-1)
      return casadi_limits<DataType>::zero;
    else
      return at(ind);
  }

  template<typename DataType>
  int Matrix<DataType>::stream_precision_ = 6;
  template<typename DataType>
  int Matrix<DataType>::stream_width_ = 0;
  template<typename DataType>
  bool Matrix<DataType>::stream_scientific_ = false;

  template<typename DataType>
  DataType& Matrix<DataType>::elem(int rr, int cc) {
    int oldsize = sparsity().size();
    int ind = sparsityRef().getNZ(rr, cc);
    if (oldsize != sparsity().size())
      data().insert(begin()+ind, DataType(0));
    return at(ind);
  }

  template<typename DataType>
  bool Matrix<DataType>::__nonzero__() const {
    if (numel()!=1) {casadi_error("Only scalar Matrix could have a truth value, but you "
                                  "provided a shape" << dimString());}
    return at(0)!=0;
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::sub(int rr, int cc) const {
    return elem(rr, cc);
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::sub(const std::vector<int>& jj,
                                               const std::vector<int>& ii) const {
    // Nonzero mapping from submatrix to full
    std::vector<int> mapping;

    // Get the sparsity pattern - does bounds checking
    Sparsity sp = sparsity().sub(jj, ii, mapping);

    // Create return object
    Matrix<DataType> ret(sp);

    // Copy nonzeros
    for (int k=0; k<mapping.size(); ++k)
      ret.data()[k] = data()[mapping[k]];

    // Return (RVO)
    return ret;
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::sub(const Matrix<int>& k,
                                               const std::vector<int>& ii) const {
    std::vector< int > rows = range(size1());
    std::vector< Matrix<DataType> > temp;

    if (!inBounds(ii, size2())) {
      casadi_error("Slicing [ii, k] out of bounds. Your ii contains "
                   << *std::min_element(ii.begin(), ii.end()) << " up to "
                   << *std::max_element(ii.begin(), ii.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    for (int i=0;i<ii.size();++i) {
      Matrix<DataType> m = k;
      for (int j=0;j<m.size();++j) {
        m.data()[j] = elem(k.at(j), ii.at(i));
      }
      temp.push_back(m);
    }

    return horzcat(temp);
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::sub(const std::vector<int>& jj,
                                               const Matrix<int>& k) const {
    std::vector< int > cols = range(size2());
    std::vector< Matrix<DataType> > temp;

    if (!inBounds(jj, size1())) {
      casadi_error("Slicing [ii, k] out of bounds. Your jj contains "
                   << *std::min_element(jj.begin(), jj.end()) << " up to "
                   << *std::max_element(jj.begin(), jj.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    for (int j=0;j<jj.size();++j) {
      Matrix<DataType> m = k;
      for (int i=0;i<m.size();++i) {
        m.data()[i] = elem(jj.at(j), k.at(i));
      }
      temp.push_back(m);
    }

    return vertcat(temp);
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::sub(const Matrix<int>& j, const Matrix<int>& i) const {
    casadi_assert_message(i.sparsity()==j.sparsity(),
                          "sub(Imatrix i, Imatrix j): sparsities must match. Got "
                          << i.dimString() << " and " << j.dimString() << ".");

    Matrix<DataType> ret(i.sparsity());
    for (int k=0;k<i.size();++k) {
      ret.data()[k] = elem(j.at(k), i.at(k));
    }

    return ret;
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::sub(const Sparsity& sp, int dummy) const {
    casadi_assert_message(
      size1()==sp.size1() && size2()==sp.size2(),
      "sub(Sparsity sp): shape mismatch. This matrix has shape "
      << size1() << " x " << size2()
      << ", but supplied sparsity index has shape "
      << sp.size1() << " x " << sp.size2() << ".");
    Matrix<DataType> ret(sp);

    std::vector<unsigned char> mapping; // Mapping that will be filled by patternunion
    sparsity().patternCombine(sp, false, true, mapping);

    int k = 0;     // Flat index into non-zeros of this matrix
    int j = 0;     // Flat index into non-zeros of the resultant matrix;
    for (int i=0;i<mapping.size();++i) {
      if (mapping[i] & 1) { // If the original matrix has a non-zero entry in the union
        if (!(mapping[i] & 4)) ret.at(j) = at(k); // If this non-zero entry appears in the
                                                  // intersection, add it to the mapping
        k++;                 // Increment the original matrix' non-zero index counter
      }
      if (mapping[i] & 2) j++;
    }

    return ret;
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, int j, int i) {
    if (m.isDense()) {
      elem(j, i) = m.toScalar();
    } else {
      setSub(m, std::vector<int>(1, j), std::vector<int>(1, i));
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, const std::vector<int>& rr,
                                const std::vector<int>& cc) {
    casadi_assert_message(m.numel()==1 || (cc.size() == m.size2() && rr.size() == m.size1()),
                          "Dimension mismatch." << std::endl << "lhs is " << rr.size() << " x "
                          << cc.size() << ", while rhs is " << m.dimString());

    if (!inBounds(rr, size1())) {
      casadi_error("setSub[., rr, cc] out of bounds. Your rr contains "
                   << *std::min_element(rr.begin(), rr.end()) << " up to "
                   << *std::max_element(rr.begin(), rr.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(cc, size2())) {
      casadi_error("setSub [., rr, cc] out of bounds. Your cc contains "
                   << *std::min_element(cc.begin(), cc.end()) << " up to "
                   << *std::max_element(cc.begin(), cc.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    // If m is scalar
    if (m.numel() != cc.size() * rr.size()) {
      setSub(Matrix<DataType>::repmat(m.toScalar(), rr.size(), cc.size()), rr, cc);
      return;
    }

    if (isDense() && m.isDense()) {
      // Dense mode
      for (int i=0; i<cc.size(); ++i) {
        for (int j=0; j<rr.size(); ++j) {
          data()[cc[i]*size1() + rr[j]]=m.data()[i*m.size1()+j];
        }
      }
    } else {
      // Sparse mode

      // Remove submatrix to be replaced
      erase(rr, cc);

      // Extend el to the same dimension as this
      Matrix<DataType> el_ext = m;
      el_ext.enlarge(size1(), size2(), rr, cc);

      // Unite the sparsity patterns
      *this = unite(*this, el_ext);
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, const std::vector<int>& jj,
                                const Matrix<int>& i) {
    // If el is scalar
    if (m.isScalar() && (jj.size() > 1 || i.size() > 1)) {
      setSub(repmat(Matrix<DataType>(i.sparsity(), m.toScalar()), jj.size(), 1), jj, i);
      return;
    }

    if (!inBounds(jj, size1())) {
      casadi_error(
        "setSub[., i, jj] out of bounds. Your jj contains "
        << *std::min_element(jj.begin(), jj.end()) << " up to "
        << *std::max_element(jj.begin(), jj.end())
        << ", which is outside of the matrix shape " << dimString() << ".");
    }

    Sparsity result_sparsity = repmat(i, jj.size(), 1).sparsity();

    casadi_assert_message(
      result_sparsity == m.sparsity(),
      "setSub(Imatrix" << i.dimString() << ", Ivector(length=" << jj.size()
      << "), Matrix<DataType>)::Dimension mismatch. The sparsity of repmat(IMatrix, "
      << jj.size() << ",1) = " << result_sparsity.dimString()
      << " must match the sparsity of Matrix<DataType> = "  << m.dimString()
      << ".");

    std::vector<int> slice_i = range(i.size2());

    for (int k=0; k<jj.size(); ++k) {
      Matrix<DataType> el_k = m(range(k*i.size1(), (k+1)*i.size1()), slice_i);
      for (int j=0;j<i.size();++j) {
        elem(jj[k], i.at(j))=el_k.at(j);
      }
    }

  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, const Matrix<int>& j,
                                const std::vector<int>& ii) {

    // If el is scalar
    if (m.isScalar() && (ii.size() > 1 || j.size() > 1)) {
      setSub(repmat(Matrix<DataType>(j.sparsity(), m.toScalar()), 1, ii.size()), j, ii);
      return;
    }

    if (!inBounds(ii, size2())) {
      casadi_error("setSub[., ii, j] out of bounds. Your ii contains "
                   << *std::min_element(ii.begin(), ii.end()) << " up to "
                   << *std::max_element(ii.begin(), ii.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    Sparsity result_sparsity = repmat(j, 1, ii.size()).sparsity();


    casadi_assert_message(
        result_sparsity == m.sparsity(),
        "setSub(Ivector(length=" << ii.size() << "), Imatrix" << j.dimString()
        << ", Matrix<DataType>)::Dimension mismatch. The sparsity of repmat(Imatrix,1, "
        << ii.size() << ") = " << result_sparsity.dimString()
        << " must match the sparsity of Matrix<DataType> = " << m.dimString() << ".");

    std::vector<int> slice_j = range(j.size1());

    for (int k=0; k<ii.size(); ++k) {
      Matrix<DataType> el_k = m(slice_j, range(k*j.size2(), (k+1)*j.size2()));
      for (int i=0;i<j.size();++i) {
        elem(j.at(i), ii[k])=el_k.at(i);
      }
    }

  }


  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, const Matrix<int>& j,
                                const Matrix<int>& i) {
    casadi_assert_message(i.sparsity()==j.sparsity(),
                          "setSub(., Imatrix i, Imatrix j): sparsities must match. Got "
                          << i.dimString() << " for i and " << j.dimString() << " for j.");

    // If m is scalar
    if (m.isScalar() && i.numel() > 1) {
      setSub(Matrix<DataType>(i.sparsity(), m.toScalar()), j, i);
      return;
    }

    casadi_assert_message(
      m.sparsity()==i.sparsity(),
      "setSub(Matrix m, Imatrix i, Imatrix j): sparsities must match. Got "
      << m.dimString() << " for m and " << j.dimString() << " for i and j.");

    for (int k=0; k<i.size(); ++k) {
      elem(j.at(k), i.at(k)) = m.at(k);
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setSub(const Matrix<DataType>& m, const Sparsity& sp, int dummy) {
    casadi_assert_message(
      size2()==sp.size2() && size1()==sp.size1(),
      "sub(Sparsity sp): shape mismatch. This matrix has shape "
      << size2() << " x " << size1()
      << ", but supplied sparsity index has shape "
      << sp.size2() << " x " << sp.size1() << ".");
    // TODO(Joel): optimize this for speed
    Matrix<DataType> elm;
    if (m.isScalar()) {
      elm = Matrix<DataType>(sp, m.at(0));
    } else {
      elm = m.sub(sp);
    }

    for (int i=0; i<sp.colind().size()-1; ++i) {
      for (int k=sp.colind()[i]; k<sp.colind()[i+1]; ++k) {
        int j=sp.row()[k];
        elem(j, i)=elm.data()[k];
      }
    }
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getNZ(const std::vector<int>& k) const {
    try {
      Matrix<DataType> ret = zeros(k.size());
      for (int el=0; el<k.size(); ++el)
        ret.data()[el] = data().at(k[el]);

      return ret;
    } catch(std::out_of_range& ex) {
      std::stringstream ss;
      ss << "Out of range error in Matrix<>::getNZ: " << k
         << " not all in range [0, " << size() << ")";
      throw CasadiException(ss.str());
    }
  }

  template<typename DataType>
  const Matrix<DataType> Matrix<DataType>::getNZ(const Matrix<int>& k) const {
    try {
      Matrix<DataType> ret = zeros(k.sparsity());
      for (int el=0; el<k.size(); ++el)
        ret.data()[el] = data().at(k.at(el));

      return ret;
    } catch(std::out_of_range& ex) {
      std::stringstream ss;
      ss << "Out of range error in Matrix<>::getNZ: " << k
         << " not all in range [0, " << size() << ")";
      throw CasadiException(ss.str());
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setNZ(int k, const Matrix<DataType>& m) {
    if (k<0) k+=size();
    at(k) = m.toScalar();
  }

  template<typename DataType>
  void Matrix<DataType>::setNZ(const std::vector<int>& kk, const Matrix<DataType>& m) {
    if (m.isScalar()) {
      // Assign all elements with the same scalar
      for (int k=0; k<kk.size(); ++k) {
        setNZ(kk[k], m);
      }
    } else {
      // Assignment elementwise
      casadi_assert_message(kk.size()==m.size(),
                            "Matrix<DataType>::setNZ: length of non-zero indices ("
                            << kk.size() << ") " << std::endl << "must match size of rhs ("
                            << m.size() << ").");
      for (int k=0; k<kk.size(); ++k) {
        setNZ(kk[k], m[k]);
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setNZ(const Matrix<int>& kk, const Matrix<DataType>& m) {
    if (m.isScalar()) {
      // Assign all elements with the same scalar
      for (int k=0; k<kk.size(); ++k) {
        setNZ(kk.at(k), m);
      }
    } else if (kk.isDense() && !m.isDense() && kk.size2()==m.size2() && kk.size1()==m.size1()) {
      const std::vector<int>& row = m.sparsity().row();
      const std::vector<int>& colind = m.sparsity().colind();
      for (int i=0; i<colind.size()-1; ++i) {
        for (int k=colind[i]; k<colind[i+1]; ++k) {
          int j=row[k];
          setNZ(kk.elem(j, i), m[k]);
        }
      }
    } else {
      casadi_assert_message(kk.sparsity()==m.sparsity(),
                            "Matrix<DataType>::setNZ: sparsity of IMatrix index "
                            << kk.dimString() << " " << std::endl
                            << "must match sparsity of rhs " << m.dimString() << ".");
      for (int k=0; k<kk.size(); ++k) {
        setNZ(kk.at(k), m[k]);
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::densify(const DataType& val) {
    // Quick return if possible
    if (isDense()) return;

    // Get sparsity pattern
    int nrow = size1();
    int ncol = size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();

    // Resize data and copy
    data_.resize(nrow*ncol, val);

    // Loop over the columns in reverse order
    for (int cc=ncol-1; cc>=0; --cc) {
      // Loop over nonzero elements of the column in reverse order
      for (int el=colind[cc+1]-1; el>=colind[cc]; --el) {
        int rr = row[el];
        int new_el = cc*nrow + rr;
        if (el==new_el) break; // Already done, the rest of the elements must be in the same place
        std::swap(data_[new_el], data_[el]);
      }
    }

    // Update the sparsity pattern
    sparsity_ = Sparsity::dense(shape());
  }

  template<typename DataType>
  void Matrix<DataType>::sparsify(double tol) {
    // Quick return if there are no entries to be removed
    bool remove_nothing = true;
    for (typename std::vector<DataType>::iterator it=begin(); it!=end() && remove_nothing; ++it) {
      remove_nothing = !casadi_limits<DataType>::isAlmostZero(*it, tol);
    }
    if (remove_nothing) return;

    // Get the current sparsity pattern
    int size1 = this->size1();
    int size2 = this->size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();

    // Construct the new sparsity pattern
    std::vector<int> new_colind(1, 0), new_row;

    // Loop over the columns
    for (int cc=0; cc<size2; ++cc) {
      // Loop over existing nonzeros
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // If it is not known to be a zero
        if (!casadi_limits<DataType>::isAlmostZero(data_[el], tol)) {
          // Save the nonzero in its new location
          data_[new_row.size()] = data_[el];

          // Add to pattern
          new_row.push_back(row[el]);
        }
      }
      // Save the new column offset
      new_colind.push_back(new_row.size());
    }

    // Trim the data vector
    data_.resize(new_row.size());

    // Update the sparsity pattern
    sparsity_ = Sparsity(size1, size2, new_colind, new_row);
  }

  template<typename DataType>
  Matrix<DataType>::Matrix() : sparsity_(Sparsity::sparse(0, 0)) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Matrix<DataType>& m) : sparsity_(m.sparsity_), data_(m.data_) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const std::vector<DataType>& x) :
      sparsity_(Sparsity::dense(x.size(), 1)), data_(x) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const std::vector<DataType>& x, int nrow, int ncol) :
      sparsity_(Sparsity::dense(nrow, ncol)), data_(x) {
    casadi_assert_message(x.size() == nrow*ncol, "Dimension mismatch." << std::endl
                          << "You supplied a vector of length " << x.size() << ", but " << nrow
                          << " x " << ncol << " = " << nrow*ncol);
  }

  template<typename DataType>
  Matrix<DataType>& Matrix<DataType>::operator=(const Matrix<DataType>& m) {
    sparsity_ = m.sparsity_;
    data_ = m.data_;
    return *this;
  }

  template<typename DataType>
  std::string Matrix<DataType>::className() { return matrixName<DataType>(); }

  template<typename DataType>
  void Matrix<DataType>::printScalar(std::ostream &stream, bool trailing_newline) const {
    casadi_assert_message(numel()==1, "Not a scalar");

    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();

    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }

    if (size()==0) {
      stream << "00";
    } else {
      stream << toScalar();
    }

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename DataType>
  void Matrix<DataType>::printVector(std::ostream &stream, bool trailing_newline) const {
    casadi_assert_message(isVector(), "Not a vector");

    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();

    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }

    // Access data structures
    const std::vector<int>& r = row();

    // Nonzero
    int el=0;

    // Loop over rows
    stream << "[";
    for (int rr=0; rr<size1(); ++rr) {
      // Add delimiter
      if (rr!=0) stream << ", ";

      // Check if nonzero
      if (el<r.size() && rr==r[el]) {
        stream << at(el++);
      } else {
        stream << "00";
      }
    }
    stream << "]";

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename DataType>
  void Matrix<DataType>::printDense(std::ostream &stream, bool trailing_newline) const {
    // Print as a single line
    bool oneliner=this->size1()<=1;

    std::streamsize precision = stream.precision();
    std::streamsize width = stream.width();
    std::ios_base::fmtflags flags = stream.flags();

    stream.precision(stream_precision_);
    stream.width(stream_width_);
    if (stream_scientific_) {
      stream.setf(std::ios::scientific);
    } else {
      stream.unsetf(std::ios::scientific);
    }

    // Index counter for each column
    std::vector<int> cind = colind();

    // Loop over rows
    for (int rr=0; rr<size1(); ++rr) {
      // Beginning of row
      if (rr==0) {
        if (!oneliner) stream << std::endl;
        stream << "[[";
      } else {
        stream << " [";
      }

      // Loop over columns
      for (int cc=0; cc<size2(); ++cc) {
        // Separating comma
        if (cc>0) stream << ", ";

        // Check if nonzero
        if (cind[cc]<colind(cc+1) && row(cind[cc])==rr) {
          stream << data().at(cind[cc]++);
        } else {
          stream << "00";
        }
      }

      // End of row
      if (rr<size1()-1) {
        stream << "], ";
        if (!oneliner) stream << std::endl;
      } else {
        stream << "]]";
      }
    }

    if (trailing_newline) stream << std::endl;
    stream << std::flush;
    stream.precision(precision);
    stream.width(width);
    stream.flags(flags);
  }

  template<typename DataType>
  void Matrix<DataType>::printSparse(std::ostream &stream, bool trailing_newline) const {
    if (size()==0) {
      stream << "all zero sparse: " << size1() << "-by-" << size2();
    } else {
      stream << "sparse: " << size1() << "-by-" << size2() << ", " << size() << " nnz";
      for (int cc=0; cc<size2(); ++cc) {
        for (int el=colind(cc); el<colind(cc+1); ++el) {
          int rr=row(el);
          stream << std::endl << " (" << rr << ", " << cc << ") -> " << at(el);
        }
      }
    }
    if (trailing_newline) stream << std::endl;
    stream << std::flush;
  }

  template<typename DataType>
  void Matrix<DataType>::print(std::ostream &stream, bool trailing_newline) const {
    if (isEmpty()) {
      stream << "[]";
    } else if (numel()==1) {
      printScalar(stream, false);
    } else if (isVector()) {
      printVector(stream, false);
    } else if (std::max(size1(), size2())<=10 || static_cast<double>(size())/numel()>=0.5) {
      // if "small" or "dense"
      printDense(stream, false);
    } else {
      printSparse(stream, false);
    }
    if (trailing_newline) stream << std::endl;
  }

  template<typename DataType>
  void Matrix<DataType>::repr(std::ostream &stream, bool trailing_newline) const {
    stream << className() << "(";
    print(stream, false);
    stream << ")";
    if (trailing_newline) stream << std::endl;
    stream << std::flush;
  }

  template<typename DataType>
  const std::vector<int>& Matrix<DataType>::row() const {
    return sparsity().row();
  }

  template<typename DataType>
  const std::vector<int>& Matrix<DataType>::colind() const {
    return sparsity_.colind();
  }

  template<typename DataType>
  int Matrix<DataType>::row(int el) const {
    return sparsity_.row(el);
  }

  template<typename DataType>
  int Matrix<DataType>::colind(int col) const {
    return sparsity_.colind(col);
  }

  template<typename DataType>
  void Matrix<DataType>::reserve(int nnz) {
    reserve(nnz, size2());
  }

  template<typename DataType>
  void Matrix<DataType>::reserve(int nnz, int ncol) {
    data().reserve(nnz);
    sparsity_.reserve(nnz, ncol);
  }

  template<typename DataType>
  void Matrix<DataType>::resize(int nrow, int ncol) {
    sparsity_.resize(nrow, ncol);
  }

  template<typename DataType>
  void Matrix<DataType>::clear() {
    sparsity_ = Sparsity::sparse(0, 0);
    data().clear();
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(double val) :
      sparsity_(Sparsity::dense(1, 1)), data_(std::vector<DataType>(1, val)) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const std::vector< std::vector<DataType> >& d) {
    // Get dimensions
    int nrow=d.size();
    int ncol=d.empty() ? 1 : d.front().size();

    // Assert consistency
    for (int rr=0; rr<nrow; ++rr) {
      casadi_assert_message(ncol==d[rr].size(),
        "Matrix<DataType>::Matrix(const std::vector< std::vector<DataType> >& d): "
        "shape mismatch" << std::endl
        << "Attempting to construct a matrix from a nested list." << std::endl
        << "I got convinced that the desired size is ("<< nrow << " x " << ncol
        << " ), but now I encounter a vector of size ("
        << d[rr].size() <<  " )" << std::endl);
    }

    // Form matrix
    sparsity_ = Sparsity::dense(nrow, ncol);
    data().resize(nrow*ncol);
    typename std::vector<DataType>::iterator it=begin();
    for (int cc=0; cc<ncol; ++cc) {
      for (int rr=0; rr<nrow; ++rr) {
        *it++ = d[rr][cc];
      }
    }
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Sparsity& sparsity, const DataType& val) :
      sparsity_(sparsity), data_(sparsity.size(), val) {
  }

  template<typename DataType>
  Matrix<DataType>::Matrix(const Sparsity& sparsity, const std::vector<DataType>& d) :
      sparsity_(sparsity), data_(d) {
    casadi_assert_message(sparsity.size()==d.size(), "Size mismatch." << std::endl
                          << "You supplied a sparsity of " << sparsity.dimString()
                          << ", but the supplied vector is of length " << d.size());
  }

  template<typename DataType>
  void Matrix<DataType>::setZero() {
    setAll(0);
  }

  template<typename DataType>
  void Matrix<DataType>::setAll(const DataType& val) {
    std::fill(begin(), end(), val);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::unary(int op, const Matrix<DataType> &x) {
    // Return value
    Matrix<DataType> ret(x.sparsity());

    // Nonzeros
    std::vector<DataType>& ret_data = ret.data();
    const std::vector<DataType>& x_data = x.data();

    // Do the operation on all non-zero elements
    for (int el=0; el<x.size(); ++el) {
      casadi_math<DataType>::fun(op, x_data[el], x_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.isDense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, 0, 0, fcn_0);
      if (!casadi_limits<DataType>::isZero(fcn_0)) { // Remove this if?
        ret.densify(fcn_0);
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::operator-() const {
    return unary(OP_NEG, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::operator+() const {
    return *this;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__add__(const Matrix<DataType> &y) const {
    return binary(OP_ADD, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__sub__(const Matrix<DataType> &y) const {
    return binary(OP_SUB, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__mul__(const Matrix<DataType> &y) const {
    return binary(OP_MUL, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__div__(const Matrix<DataType> &y) const {
    return binary(OP_DIV, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__lt__(const Matrix<DataType> &y) const {
    return binary(OP_LT, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__le__(const Matrix<DataType> &y) const {
    return binary(OP_LE, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__eq__(const Matrix<DataType> &y) const {
    return binary(OP_EQ, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__ne__(const Matrix<DataType> &y) const {
    return binary(OP_NE, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__mrdivide__(const Matrix<DataType>& b) const {
    if (b.numel()==1) return *this/b;
    throw CasadiException("mrdivide: Not implemented");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__mpower__(const Matrix<DataType>& b) const {
    if (b.numel()==1) return (*this).__pow__(b);
    throw CasadiException("mpower: Not implemented");
  }

  template<typename DataType>
  Sparsity& Matrix<DataType>::sparsityRef() {
    sparsity_.makeUnique(); // NOTE: Remove?
    return sparsity_;
  }

  template<typename DataType>
  void Matrix<DataType>::getBand(int kl, int ku, int ldres, DataType *res) const {
    // delete the content of the matrix
    for (int j=0; j<size1(); ++j) // loop over rows
      for (int s=0; s<kl+ku+1; ++s) // loop over the subdiagonals
        res[s + ldres*j] = 0;

    // loop over cols
    for (int i=0; i<size2(); ++i) {

      // loop over the non-zero elements
      for (int el=colind(i); el<colind(i+1); ++el) {
        int j=row(el);  // row

        // Check if we have not yet inside the band
        if (j<i-kl) continue;

        // Check if we are already outside the band
        if (j>i+ku) break;

        // Get the subdiagonal
        int s = i - j + ku;

        // Store the element
        res[s + ldres*j] = data()[el];
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::set(DataType val, SparsityType sp) {
    std::fill(data().begin(), data().end(), val);
  }

  template<typename DataType>
  void Matrix<DataType>::get(DataType& val, SparsityType sp) const {
    getArray(&val, 1, DENSE);
  }

  template<typename DataType>
  void Matrix<DataType>::set(const std::vector<DataType>& val, SparsityType sp) {
    setArray(val.empty() ? 0 : &val.front(), val.size(), sp);
  }

  template<typename DataType>
  void Matrix<DataType>::get(std::vector<DataType>& val, SparsityType sp) const {
    getArray(val.empty() ? 0 : &val.front(), val.size(), sp);
  }

  template<typename DataType>
  void Matrix<DataType>::set(const Matrix<DataType>& val, SparsityType sp) {
    sparsity().set(getPtr(data()), getPtr(val.data()), val.sparsity());
  }

  template<typename DataType>
  void Matrix<DataType>::setBV(const Matrix<DataType>& val) {
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    const bvec_t* bw_val = reinterpret_cast<const bvec_t*>(getPtr(val.data()));
    sparsity().set(bw_this, bw_val, val.sparsity());
  }

  template<typename DataType>
  void Matrix<DataType>::setZeroBV() {
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    std::fill(bw_this, bw_this+size(), bvec_t(0));
  }

  template<typename DataType>
  void Matrix<DataType>::borBV(const Matrix<DataType>& val) {
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    const bvec_t* bw_val = reinterpret_cast<const bvec_t*>(getPtr(val.data()));
    sparsity().bor(bw_this, bw_val, val.sparsity());
  }

  template<typename DataType>
  void Matrix<DataType>::getArrayBV(bvec_t* val, int len) const {
    casadi_assert(len==size());
    const bvec_t* bw_this = reinterpret_cast<const bvec_t*>(getPtr(data()));
    std::copy(bw_this, bw_this+len, val);
  }

  template<typename DataType>
  void Matrix<DataType>::setArrayBV(const bvec_t* val, int len) {
    casadi_assert(len==size());
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    std::copy(val, val+len, bw_this);
  }

  template<typename DataType>
  void Matrix<DataType>::borArrayBV(const bvec_t* val, int len) {
    casadi_assert(len==size());
    bvec_t* bw_this = reinterpret_cast<bvec_t*>(getPtr(data()));
    for (int i=0; i<len; ++i) *bw_this++ |= *val++;
  }

  template<typename DataType>
  void Matrix<DataType>::get(Matrix<DataType>& val, SparsityType sp) const {
    val.set(*this, sp);
  }

  template<typename DataType>
  void Matrix<DataType>::set(const DataType* val, SparsityType sp) {
    int len = sp==SPARSE ? size() :
        sp==DENSE || sp==DENSETRANS ? numel() :
        sp==SPARSESYM ? sizeU() : -1;
    setArray(val, len, sp);
  }

  template<typename DataType>
  void Matrix<DataType>::get(DataType* val, SparsityType sp) const {
    int len = sp==SPARSE ? size() :
        sp==DENSE || sp==DENSETRANS ? numel() :
        sp==SPARSESYM ? sizeU() : -1;
    getArray(val, len, sp);
  }

  template<typename DataType>
  void Matrix<DataType>::getArray(DataType* val, int len, SparsityType sp) const {
    // Get references to data for quick access
    const std::vector<DataType> &data = this->data();
    const int size1 = this->size1();
    const int size2 = this->size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();

    if (sp==SPARSE || (sp==DENSE && isDense())) {
      casadi_assert_message(len==size(),
                            "Matrix<DataType>::getArray: Dimension mismatch." << std::endl <<
                            "Trying to fetch " << len << " elements from a " << dimString()
                            << " matrix with " << size() << " non-zeros.");
      copy(data.begin(), data.end(), val);
    } else if (sp==DENSE) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::getArray: Dimension mismatch." << std::endl <<
                            "Trying to fetch " << len << " elements from a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Begin with all zeros
      std::fill(val, val+len, 0);

      // Add the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[rr+cc*size1] = data[el];
        }
      }
    } else if (sp==DENSETRANS) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::getArray: Dimension mismatch." << std::endl <<
                            "Trying to fetch " << len << " elements from a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Begin with all zeros
      std::fill(val, val+len, 0);

      // Add the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[cc+rr*size2] = data[el];
        }
      }
    } else if (sp==SPARSESYM) {
      // copy to the result vector
      int nz = 0;
      for (int cc=0; cc<size2; ++cc) {
        // Loop over the elements in the col
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          if (row[el] > cc) break; // break inner loop (only upper triangular part is used)
          val[nz++] = data[el];
        }
      }
    } else {
      casadi_error("Matrix<DataType>::getArray: not SPARSE, SPARSESYM, DENSE or DENSETRANS");
    }
  }

  /**
     Set stride to zero for unstrided acces
  */
  template<typename DataType>
  void Matrix<DataType>::getStridedArray(DataType* val, int len, int stride1, int stride2,
                                         SparsityType sp) const {
    if (stride1==0 || stride2==0 || (sp==DENSE && stride2==1 && stride1==size1()))
        return getArray(val, len, sp);

    // Get references to data for quick access
    const std::vector<DataType> &data = this->data();
    //const int size1 = this->size1();
    const int size2 = this->size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();

    if (sp==SPARSE) {
      throw CasadiException("Matrix<DataType>::getArray: strided SPARSE not implemented");
    } else if (sp==DENSE && isDense()) {
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[rr*stride2 + cc*stride1] = data[el];
        }
      }
    } else if (sp==DENSETRANS && isDense()) {
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          val[cc*stride2 + rr*stride1] = data[el];
        }
      }
    } else if (sp==DENSE) {
      throw CasadiException("Matrix<DataType>::getStridedArray: "
                            "strided sparse DMatrix->dense not implemented");
    } else if (sp==SPARSESYM) {
      throw CasadiException("Matrix<DataType>::getStridedArray: strided SPARSESYM not implemented");
    } else {
      throw CasadiException("Matrix<DataType>::getStridedArray: not SPARSE or DENSE");
    }

  }

  template<typename DataType>
  void Matrix<DataType>::setArray(const DataType* val, int len, SparsityType sp) {
    // Get references to data for quick access
    std::vector<DataType> &data = this->data();
    const int size1 = this->size1();
    const int size2 = this->size2();
    const std::vector<int>& colind = this->colind();
    const std::vector<int>& row = this->row();

    if (sp==SPARSE || (sp==DENSE && numel()==size())) {
      casadi_assert_message(len==size(),
                            "Matrix<DataType>::setArray: Dimension mismatch." << std::endl <<
                            "Trying to pass " << len << " elements to a " << dimString()
                            << " matrix with " << size() << " non-zeros.");
      copy(val, val+len, data.begin());
    } else if (sp==DENSE) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::setArray: Dimension mismatch." << std::endl <<
                            "Trying to pass " << len << " elements to a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Get the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          data[el] = val[rr+cc*size1];
        }
      }
    } else if (sp==DENSETRANS) {
      casadi_assert_message(len==numel(),
                            "Matrix<DataType>::setArray: Dimension mismatch." << std::endl <<
                            "Trying to pass " << len << " elements to a " << dimString()
                            << " matrix with " << numel() << " entries.");
      // Get the nonzeros
      for (int cc=0; cc<size2; ++cc) { // loop over columns
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          int rr=row[el];
          data[el] = val[cc+rr*size2];
        }
      }
    } else if (sp==SPARSESYM) {
      // NOTE: Has to be rewritten! sparsity().transpose(...) involves memory allocation
      // and is not threadsafe!!!
      // These routines is supposed to be used inside threadsafe code.
      std::vector<int> mapping;
      sparsity().transpose(mapping, false);
      // copy to the result vector
      int nz = 0;
      for (int cc=0; cc<size2; ++cc) {
        // Loop over the elements in the col
        for (int el=colind[cc]; el<colind[cc+1]; ++el) { // loop over the non-zero elements
          if (row[el] > cc) break; // break inner loop (only lower triangular part is used)
          data[mapping[el]] = data[el] = val[nz++];
        }
      }
    } else {
      throw CasadiException("Matrix<DataType>::setArray: "
                            "not SPARSE, SPARSESYM, DENSE or DENSETRANS");
    }
  }

  template<typename DataType>
  void Matrix<DataType>::getArray(DataType* val) const {
    getArray(val, size(), SPARSE);
  }

  template<typename DataType>
  void Matrix<DataType>::setArray(const DataType* val) {
    setArray(val, size(), SPARSE);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__pow__(const Matrix<DataType>& y) const {
    return binary(OP_POW, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__constpow__(const Matrix<DataType>& y) const {
    return binary(OP_CONSTPOW, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::sin() const {
    return unary(OP_SIN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::cos() const {
    return unary(OP_COS, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::tan() const {
    return unary(OP_TAN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::erf() const {
    return unary(OP_ERF, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::arcsin() const {
    return unary(OP_ASIN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::arccos() const {
    return unary(OP_ACOS, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::arctan() const {
    return unary(OP_ATAN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::sinh() const {
    return unary(OP_SINH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::cosh() const {
    return unary(OP_COSH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::tanh() const {
    return unary(OP_TANH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::arcsinh() const {
    return unary(OP_ASINH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::arccosh() const {
    return unary(OP_ACOSH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::arctanh() const {
    return unary(OP_ATANH, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::exp() const {
    return unary(OP_EXP, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::log() const {
    return unary(OP_LOG, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::log10() const {
    return log()*(1/std::log(10.));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::sqrt() const {
    return unary(OP_SQRT, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::floor() const {
    return unary(OP_FLOOR, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::ceil() const {
    return unary(OP_CEIL, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::fmod(const Matrix<DataType>& y) const {
    return binary(OP_FMOD, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::fabs() const {
    return unary(OP_FABS, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::sign() const {
    return unary(OP_SIGN, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::__copysign__(const Matrix<DataType>& y) const {
    return binary(OP_COPYSIGN, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::erfinv() const {
    return unary(OP_ERFINV, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::fmin(const Matrix<DataType>& y) const {
    return binary(OP_FMIN, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::arctan2(const Matrix<DataType>& y) const {
    return binary(OP_ATAN2, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::fmax(const Matrix<DataType>& y) const {
    return binary(OP_FMAX, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::printme(const Matrix<DataType>& y) const {
    return binary(OP_PRINTME, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::logic_not() const {
    return unary(OP_NOT, *this);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::logic_and(const Matrix<DataType>& y) const {
    return binary(OP_AND, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::logic_or(const Matrix<DataType>& y) const {
    return binary(OP_OR, *this, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::if_else_zero(const Matrix<DataType>& y) const {
    return binary(OP_IF_ELSE_ZERO, *this, y);
  }

  template<typename DataType>
  std::vector<DataType>& Matrix<DataType>::data() {
    return data_;
  }

  template<typename DataType>
  const std::vector<DataType>& Matrix<DataType>::data() const {
    return data_;
  }

  template<typename DataType>
  void Matrix<DataType>::erase(const std::vector<int>& rr, const std::vector<int>& cc) {
    // Erase from sparsity pattern
    std::vector<int> mapping = sparsityRef().erase(rr, cc);

    // Update non-zero entries
    for (int k=0; k<mapping.size(); ++k)
      data()[k] = data()[mapping[k]];

    // Truncate nonzero vector
    data().resize(mapping.size());
  }


  template<typename DataType>
  void Matrix<DataType>::remove(const std::vector<int>& rr, const std::vector<int>& cc) {
    if (!inBounds(rr, size1())) {
      casadi_error("Remove(rr, cc) out of bounds. Your rr contains "
                   << *std::min_element(rr.begin(), rr.end()) << " up to "
                   << *std::max_element(rr.begin(), rr.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(cc, size2())) {
      casadi_error("Remove(rr, cc) out of bounds. Your cc contains "
                   << *std::min_element(cc.begin(), cc.end()) << " up to "
                   << *std::max_element(cc.begin(), cc.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    // Remove by performing a complementary slice
    std::vector<int> rrc = complement(rr, size1());
    std::vector<int> ccc = complement(cc, size2());

    Matrix<DataType> ret = (*this)(rrc, ccc);

    operator=(ret);

  }

  template<typename DataType>
  void Matrix<DataType>::enlarge(int nrow, int ncol, const std::vector<int>& rr,
                                 const std::vector<int>& cc) {
    sparsityRef().enlarge(nrow, ncol, rr, cc);
  }

  template<typename DataType>
  void Matrix<DataType>::sanityCheck(bool complete) const {
    sparsity_.sanityCheck(complete);

    if (data_.size()!=sparsity_.row().size()) {
      std::stringstream s;
      s << "Matrix:Compressed Col Storage is not sane. The following must hold:" << std::endl;
      s << "  data.size() = nrow, but got   row.size()  = " << data_.size()
        << "   and   nrow = "  << sparsity_.row().size() << std::endl;
      s << "  Note that the signature is as follows: DMatrix (ncol, nrow, row, colind, data)."
        << std::endl;
      casadi_error(s.str());
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::mul(const Matrix<DataType> &y, const Sparsity& sp_z) const {
    return this->mul_smart(y, sp_z);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::mul_full(const Matrix<DataType> &y,
                                              const Sparsity& sp_z) const {
    // First factor
    const Matrix<DataType>& x = *this;

    // Return object (assure RVO)
    Matrix<DataType> ret;
    if (sp_z.isNull()) {
      // Create the sparsity pattern for the matrix-matrix product
      Sparsity spres = x.sparsity().patternProductNew(y.sparsity());

      // Create the return object
      ret = Matrix<DataType>::zeros(spres);
    } else {
      ret = Matrix<DataType>::zeros(sp_z);
    }

    // Carry out the matrix product
    mul_no_alloc_nn(x, y, ret);

    return ret;
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_nn(const Matrix<DataType> &x, const Matrix<DataType> &y,
                                         Matrix<DataType>& z, std::vector<DataType>& work) {
    // Dimensions of the result
    int d1 = x.size1();
    int d2 = y.size2();

    // Assert dimensions
    casadi_assert_message(d1==z.size1(), "Dimension error. Got x=" << x.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(d2==z.size2(), "Dimension error. Got y=" << y.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==x.size2(), "Dimension error. Got y=" << y.dimString()
                          << " and x=" << x.dimString() << ".");

    // Assert work vector large enough
    casadi_assert_message(work.size()>=d1, "Work vector too small. Got length "
                          << work.size() << " < " << x.size1());

    // Direct access to the arrays
    const std::vector<int> &y_colind = y.colind();
    const std::vector<int> &y_row = y.row();
    const std::vector<DataType> &y_data = y.data();
    const std::vector<int> &x_colind = x.colind();
    const std::vector<int> &x_row = x.row();
    const std::vector<DataType> &x_data = x.data();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &z_row = z.row();
    std::vector<DataType> &z_data = z.data();

    // Loop over the columns of y and z
    for (int cc=0; cc<d2; ++cc) {
      // Get the dense column of z
      for (int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        work[z_row[kk]] = z_data[kk];
      }

      // Loop over the nonzeros of y
      for (int kk=y_colind[cc]; kk<y_colind[cc+1]; ++kk) {
        int rr = y_row[kk];
        DataType yy = y_data[kk];

        // Loop over corresponding columns of x
        for (int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
          work[x_row[kk1]] += x_data[kk1]*yy;
        }
      }

      // Get the sparse column of z
      for (int kk=z_colind[cc]; kk<z_colind[cc+1]; ++kk) {
        z_data[kk] = work[z_row[kk]];
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_nn(const Matrix<DataType> &x, const Matrix<DataType> &y,
                                         Matrix<DataType>& z) {
    // Assert dimensions
    casadi_assert_message(x.size1()==z.size1(), "Dimension error. Got x=" << x.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size2()==z.size2(), "Dimension error. Got y=" << y.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==x.size2(), "Dimension error. Got y=" << y.dimString()
                          << " and x=" << x.dimString() << ".");

    // Direct access to the arrays
    const std::vector<int> &y_colind = y.colind();
    const std::vector<int> &y_row = y.row();
    const std::vector<DataType> &y_data = y.data();
    const std::vector<int> &x_colind = x.colind();
    const std::vector<int> &x_row = x.row();
    const std::vector<DataType> &x_data = x.data();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &z_row = z.row();
    std::vector<DataType> &z_data = z.data();

    // loop over the cols of the first argument
    for (int i=0; i<y_colind.size()-1; ++i) {
      // loop over the non-zeros of the first argument
      for (int el=y_colind[i]; el<y_colind[i+1]; ++el) {
        int j = y_row[el];
        int el1 = z_colind[i];
        int el2 = x_colind[j];
        while (el1 < z_colind[i+1] && el2 < x_colind[j+1]) { // loop over matching non-zero elements
          int j1 = z_row[el1];
          int i2 = x_row[el2];
          if (j1==i2) {
            z_data[el1++] += y_data[el]*x_data[el2++];
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_tn(const Matrix<DataType> &x_trans,
                                         const std::vector<DataType> &y,
                                         std::vector<DataType>& z) {
    // Assert dimensions
    casadi_assert_message(x_trans.size2()==z.size(), "Dimension error. Got x_trans="
                          << x_trans.dimString() << " and z=" << z.size() << ".");
    casadi_assert_message(x_trans.size1()==y.size(), "Dimension error. Got x_trans="
                          << x_trans.dimString() << " and y=" << y.size() << ".");

    // Direct access to the arrays
    const std::vector<int> &x_rowind = x_trans.colind();
    const std::vector<int> &x_col = x_trans.row();
    const std::vector<DataType> &x_trans_data = x_trans.data();

    // loop over the columns of the matrix
    for (int i=0; i<x_rowind.size()-1; ++i) {
      for (int el=x_rowind[i]; el<x_rowind[i+1]; ++el) { // loop over the non-zeros of the matrix
        int j = x_col[el];

        // Perform operation
        z[i] += x_trans_data[el] * y[j];
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_nn(const Matrix<DataType>& x,
                                         const std::vector<DataType> &y,
                                         std::vector<DataType> &z) {
    // Assert dimensions
    casadi_assert_message(x.size1()==z.size(), "Dimension error. Got x=" << x.dimString()
                          << " and z=" << z.size() << ".");
    casadi_assert_message(x.size2()==y.size(), "Dimension error. Got x=" << x.dimString()
                          << " and y=" << y.size() << ".");

    // Direct access to the arrays
    const std::vector<int> &x_colind = x.colind();
    const std::vector<int> &x_row = x.row();
    const std::vector<DataType> &x_data = x.data();

    // loop over the rows of the matrix
    for (int i=0; i<x_colind.size()-1; ++i) {
      for (int el=x_colind[i]; el<x_colind[i+1]; ++el) { // loop over the non-zeros of the matrix
        int j = x_row[el];
        z[j] += x_data[el] * y[i];
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_nt(const Matrix<DataType> &x,
                                         const Matrix<DataType>& y_trans,
                                         Matrix<DataType> &z) {
    // Assert dimensions
    casadi_assert_message(y_trans.size1()==z.size2(), "Dimension error. Got y_trans="
                          << y_trans.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(x.size1()==z.size1(), "Dimension error. Got x=" << x.dimString()
                          << " and z=" << z.dimString() << ".");
    casadi_assert_message(y_trans.size2()==x.size2(), "Dimension error. Got y_trans="
                          << y_trans.dimString() << " and x=" << x.dimString() << ".");

    // Direct access to the arrays
    const std::vector<int> &y_rowind = y_trans.colind();
    const std::vector<int> &y_col = y_trans.row();
    const std::vector<DataType> &y_trans_data = y_trans.data();
    const std::vector<int> &x_colind = x.colind();
    const std::vector<int> &x_row = x.row();
    const std::vector<DataType> &x_data = x.data();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &z_row = z.row();
    std::vector<DataType> &z_data = z.data();

    // loop over the rows of the first argument
    for (int i=0; i<y_rowind.size()-1; ++i) {
      // loop over the non-zeros of the first argument
      for (int el=y_rowind[i]; el<y_rowind[i+1]; ++el) {
        int j = y_col[el];
        int el1 = x_colind[i];
        int el2 = z_colind[j];
        while (el1 < x_colind[i+1] && el2 < z_colind[j+1]) { // loop over matching non-zero elements
          int j1 = x_row[el1];
          int i2 = z_row[el2];
          if (j1==i2) {
            z_data[el2++] += y_trans_data[el] * x_data[el1++];
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<typename DataType>
  void Matrix<DataType>::mul_no_alloc_tn(const Matrix<DataType> &x_trans,
                                         const Matrix<DataType> &y,
                                         Matrix<DataType>& z) {
    // Assert dimensions
    casadi_assert_message(y.size2()==z.size2(), "Dimension error. Got y="
                          << y.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(x_trans.size2()==z.size1(), "Dimension error. Got x_trans="
                          << x_trans.dimString() << " and z=" << z.dimString() << ".");
    casadi_assert_message(y.size1()==x_trans.size1(), "Dimension error. Got y="
                          << y.dimString() << " and x_trans=" << x_trans.dimString() << ".");

    // Direct access to the arrays
    const std::vector<int> &y_colind = y.colind();
    const std::vector<int> &y_row = y.row();
    const std::vector<DataType> &y_data = y.data();
    const std::vector<int> &x_rowind = x_trans.colind();
    const std::vector<int> &x_col = x_trans.row();
    const std::vector<DataType> &x_trans_data = x_trans.data();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &z_row = z.row();
    std::vector<DataType> &z_data = z.data();

    // loop over the cols of the resulting matrix
    for (int i=0; i<z_colind.size()-1; ++i) {
      // loop over the non-zeros of the resulting matrix
      for (int el=z_colind[i]; el<z_colind[i+1]; ++el) {
        int j = z_row[el];
        int el1 = y_colind[i];
        int el2 = x_rowind[j];
        while (el1 < y_colind[i+1] && el2 < x_rowind[j+1]) { // loop over non-zero elements
          int j1 = y_row[el1];
          int i2 = x_col[el2];
          if (j1==i2) {
            z_data[el] += y_data[el1++] * x_trans_data[el2++];
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<typename DataType>
  template<bool Fwd>
  void Matrix<DataType>::mul_sparsity(Matrix<DataType> &x_trans,
                                      Matrix<DataType> &y,
                                      Matrix<DataType>& z) {
    // Direct access to the arrays
    const std::vector<int> &z_row = z.row();
    const std::vector<int> &z_colind = z.colind();
    const std::vector<int> &y_row = y.row();
    const std::vector<int> &x_col = x_trans.row();
    const std::vector<int> &y_colind = y.colind();
    const std::vector<int> &x_rowind = x_trans.colind();

    // Convert data array to arrays of integers
    bvec_t *y_data = get_bvec_t(y.data());
    bvec_t *x_trans_data = get_bvec_t(x_trans.data());
    bvec_t *z_data = get_bvec_t(z.data());

    // loop over the cols of the resulting matrix)
    for (int i=0; i<z_colind.size()-1; ++i) {
      // loop over the non-zeros of the resulting matrix
      for (int el=z_colind[i]; el<z_colind[i+1]; ++el) {
        int j = z_row[el];
        int el1 = y_colind[i];
        int el2 = x_rowind[j];
        while (el1 < y_colind[i+1] && el2 < x_rowind[j+1]) { // loop over non-zero elements
          int j1 = y_row[el1];
          int i2 = x_col[el2];
          if (j1==i2) {
            // | and not & since we are propagating dependencies
            if (Fwd) {
              z_data[el] |= y_data[el1] | x_trans_data[el2];
            } else {
              y_data[el1] |= z_data[el];
              x_trans_data[el2] |= z_data[el];
            }
            el1++;
            el2++;
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
      }
    }
  }

  template<typename DataType>
  DataType Matrix<DataType>::quad_form(const std::vector<DataType>& x, const Matrix<DataType>& A) {
    // Assert dimensions
    casadi_assert_message(x.size()==A.size2() && x.size()==A.size1(),
                          "Dimension mismatch. Got x=" << x.size() << " and A=" << A.dimString());

    // Access the internal data of A
    const std::vector<int> &A_colind = A.colind();
    const std::vector<int> &A_row = A.row();
    const std::vector<DataType> &A_data = A.data();

    // Return value
    DataType ret=0;

    // Loop over the cols of A
    for (int i=0; i<x.size(); ++i) {
      // Loop over the nonzeros of A
      for (int el=A_colind[i]; el<A_colind[i+1]; ++el) {
        // Get row
        int j = A_row[el];

        // Add contribution
        ret += x[i]*A_data[el]*x[j];
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::trans() const {
    // quick return if empty or scalar
    if ((size1()==0 && size2()==0) || isScalar()) return *this;

    // Create the new sparsity pattern and the mapping
    std::vector<int> mapping;
    Sparsity s = sparsity().transpose(mapping);

    // create the return matrix
    Matrix<DataType> ret(s);

    // Copy the content
    for (int i=0; i<mapping.size(); ++i)
      ret.at(i) = at(mapping[i]);

    return ret;
  }

  template<typename DataType>
  const DataType Matrix<DataType>::toScalar() const {
    // Make sure that the matrix is 1-by-1
    casadi_assert_message(isScalar(), "Can only convert 1-by-1 matrices to scalars");

    // return zero or the nonzero element
    if (size()==1)
      return data()[0];
    else
      return casadi_limits<DataType>::zero;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::binary(int op,
                                            const Matrix<DataType> &x,
                                            const Matrix<DataType> &y) {
    if (x.numel()==1)
      return scalar_matrix(op, x, y);
    else if (y.numel()==1)
      return matrix_scalar(op, x, y);
    else
      return matrix_matrix(op, x, y);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::scalar_matrix(int op,
                                                   const Matrix<DataType> &x,
                                                   const Matrix<DataType> &y) {
    // Return value
    Matrix<DataType> ret(y.sparsity());

    // Nonzeros
    std::vector<DataType>& ret_data = ret.data();
    const std::vector<DataType>& x_data = x.data();
    const DataType& x_val = x_data.empty() ? casadi_limits<DataType>::zero : x.front();
    const std::vector<DataType>& y_data = y.data();

    // Do the operation on all non-zero elements
    for (int el=0; el<y.size(); ++el) {
      casadi_math<DataType>::fun(op, x_val, y_data[el], ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!y.isDense() && !operation_checker<Function0Checker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, x_val, casadi_limits<DataType>::zero, fcn_0);
      if (!casadi_limits<DataType>::isZero(fcn_0)) { // Remove this if?
        ret.densify(fcn_0);
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::matrix_scalar(int op,
                                                   const Matrix<DataType> &x,
                                                   const Matrix<DataType> &y) {
    // Return value
    Matrix<DataType> ret(x.sparsity());

    // Nonzeros
    std::vector<DataType>& ret_data = ret.data();
    const std::vector<DataType>& x_data = x.data();
    const std::vector<DataType>& y_data = y.data();
    const DataType& y_val = y_data.empty() ? casadi_limits<DataType>::zero : y.front();

    // Do the operation on all non-zero elements
    for (int el=0; el<x.size(); ++el) {
      casadi_math<DataType>::fun(op, x_data[el], y_val, ret_data[el]);
    }

    // Check the value of the structural zero-entries, if there are any
    if (!x.isDense() && !operation_checker<F0XChecker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, casadi_limits<DataType>::zero, y_val, fcn_0);
      if (!casadi_limits<DataType>::isZero(fcn_0)) { // Remove this if?
        ret.densify(fcn_0);
      }
    }

    return ret;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::matrix_matrix(int op,
                                                   const Matrix<DataType> &x,
                                                   const Matrix<DataType> &y) {

    if (!(x.size2() == y.size2() && x.size1() == y.size1())) {
      std::stringstream ss;
      casadi_math<DataType>::print(op, ss, "lhs", "rhs");
      casadi_error("matrix_matrix: dimension mismatch in element-wise matrix operation "
                   << ss.str() <<"." << std::endl << "Left argument has shape " << x.dimString()
                   << ", right has shape " << y.dimString() << ". They should be equal.");
    }

    // Get the sparsity pattern of the result
    // (ignoring structural zeros giving rise to nonzero result)
    const Sparsity& x_sp = x.sparsity();
    const Sparsity& y_sp = y.sparsity();
    Sparsity r_sp = x_sp.patternCombine(y_sp,
                                        operation_checker<F0XChecker>(op),
                                        operation_checker<Function0Checker>(op));

    // Return value
    Matrix<DataType> r(r_sp);

    // Perform the operations elementwise
    if (x_sp==y_sp) {
      // Matching sparsities
      casadi_math<DataType>::fun(op, getPtr(x.data()), getPtr(y.data()),
                                 getPtr(r.data()), r_sp.size());
    } else if (y_sp==r_sp) {
      // Project first argument
      Matrix<DataType> x_mod = x(r_sp);
      casadi_math<DataType>::fun(op, getPtr(x_mod.data()), getPtr(y.data()),
                                 getPtr(r.data()), r_sp.size());
    } else if (x_sp==r_sp) {
      // Project second argument
      Matrix<DataType> y_mod = y(r_sp);
      casadi_math<DataType>::fun(op, getPtr(x.data()),
                                 getPtr(y_mod.data()), getPtr(r.data()), r_sp.size());
    } else {
      // Project both arguments
      Matrix<DataType> x_mod = x(r_sp);
      Matrix<DataType> y_mod = y(r_sp);
      casadi_math<DataType>::fun(op, getPtr(x_mod.data()), getPtr(y_mod.data()),
                                 getPtr(r.data()), r_sp.size());
    }

    // Handle structural zeros giving rise to nonzero result, e.g. cos(0) == 1
    if (!r.isDense() && !operation_checker<F00Checker>(op)) {
      // Get the value for the structural zeros
      DataType fcn_0;
      casadi_math<DataType>::fun(op, casadi_limits<DataType>::zero,
                                 casadi_limits<DataType>::zero, fcn_0);
      r.densify(fcn_0);
    }

    return r;
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const std::vector<DataType>& d) {
    return triplet(row, col, d, *std::max_element(row.begin(), row.end()),
                   *std::max_element(col.begin(), col.end()));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const std::vector<DataType>& d,
                                             const std::pair<int, int>& rc) {
    return triplet(row, col, d, rc.first, rc.second);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::triplet(const std::vector<int>& row,
                                             const std::vector<int>& col,
                                             const std::vector<DataType>& d,
                                             int nrow, int ncol) {
    casadi_assert_message(col.size()==row.size() && col.size()==d.size(),
                          "Argument error in Matrix<DataType>::sparse(row, col, d): "
                          "supplied lists must all be of equal length, but got: "
                          << row.size() << ", " << col.size()  << " and " << d.size());
    std::vector<int> mapping;
    Sparsity sp = Sparsity::triplet(nrow, ncol, row, col, mapping);
    std::vector<DataType> v(mapping.size());
    for (int k=0; k<v.size(); ++k) v[k] = d[mapping[k]];
    return Matrix<DataType>(sp, v);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::repmat(const DataType& x, const Sparsity& sp) {
    return Matrix<DataType>(sp, x);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::repmat(const Matrix<DataType>& x, const Sparsity& sp) {
    casadi_assert_message(x.isScalar(),
                          "repmat(Matrix<DataType> x, Sparsity sp) only defined for scalar x");
    return Matrix<DataType>(sp, x.toScalar());
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::repmat(const Matrix<DataType>& x,
                                            const std::pair<int, int>& rc) {
    return repmat(x, rc.first, rc.second);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::repmat(const Matrix<DataType>& x, int nrow, int ncol) {
    if (x.isScalar()) {
      if (x.isDense()) {
        return Matrix<DataType>(Sparsity::dense(nrow, ncol), x.toScalar());
      } else {
        return sparse(nrow, ncol);
      }
    } else {
      return vertcat(
        std::vector< Matrix<DataType> >(nrow, horzcat(std::vector< Matrix<DataType> >(ncol, x))));
    }
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::eye(int n) {
    return Matrix<DataType>(Sparsity::diag(n), 1);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::inf(const Sparsity& sp) {
    casadi_assert_message(std::numeric_limits<DataType>::has_infinity,
                          "Datatype cannot represent infinity");
    return Matrix<DataType>(sp, std::numeric_limits<DataType>::infinity());
  }


  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::inf(const std::pair<int, int>& rc) {
    return inf(rc.first, rc.second);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::inf(int nrow, int ncol) {
    return inf(Sparsity::dense(nrow, ncol));
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::nan(const Sparsity& sp) {
    casadi_assert_message(std::numeric_limits<DataType>::has_quiet_NaN,
                          "Datatype cannot represent not-a-number");
    return Matrix<DataType>(sp, std::numeric_limits<DataType>::quiet_NaN());
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::nan(const std::pair<int, int>& rc) {
    return nan(rc.first, rc.second);
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::nan(int nrow, int ncol) {
    return nan(Sparsity::dense(nrow, ncol));
  }

  template<typename DataType>
  void Matrix<DataType>::append(const Matrix<DataType>& y) {
    // Quick return if expr is empty
    if (size2()==0 && size1()==0) {
      *this=y;
      return;
    }

    // Quick return if empty
    if (y.size2()==0 && y.size1()==0) return;

    // Appending can be done efficiently if vectors
    if (isVector()) {
      // Append the sparsity pattern vertically
      sparsityRef().append(y.sparsity());

      // Add the non-zeros at the end
      data().insert(end(), y.begin(), y.end());
    } else {
      // Fall back on vertical concatenation
      *this = vertcat(*this, y);
    }
  }

  template<typename DataType>
  void Matrix<DataType>::appendColumns(const Matrix<DataType>& y) {

    // Quick return if expr is empty
    if (size2()==0 && size1()==0) {
      *this=y;
      return;
    }

    // Quick return if empty
    if (y.size2()==0 && y.size1()==0) return;

    // Append the sparsity pattern
    sparsityRef().appendColumns(y.sparsity());

    // Add the non-zeros at the end
    data().insert(end(), y.begin(), y.end());
  }

  template<typename DataType>
  NonZeroIterator<DataType>::NonZeroIterator(const Matrix<DataType> & m)
    : m_(m) {
    nz.i  = 0;
    nz.j  = 0;
    nz.k  = 0;
  }

  template<typename DataType>
  bool NonZeroIterator<DataType>::operator==(const NonZeroIterator<DataType>& rhs)
  { return (m_ == rhs.m_) && (nz.k==rhs.nz.k); }

  template<typename DataType>
  NonZero<DataType>& NonZeroIterator<DataType>::operator*() {
    return nz;
  }

  template<typename DataType>
  NonZeroIterator<DataType>& NonZeroIterator<DataType>::operator++() {
    nz.k ++;

    if (nz.k < m_.size()) {
      nz.j = m_.row()[nz.k];
      nz.el = m_.data()[nz.k];
      while (nz.k>=m_.colind(nz.i)) {nz.i++; }
    }
    return *this;
  }

  template<typename DataType>
  NonZeroIterator<DataType> NonZeroIterator<DataType>::begin() {
    NonZeroIterator<DataType> it = NonZeroIterator<DataType>(m_);
    return it;

  }
  template<typename DataType>
  NonZeroIterator<DataType> NonZeroIterator<DataType>::end() {
    NonZeroIterator<DataType> it = NonZeroIterator<DataType>(m_);
    it.nz.k = m_.size()-1;
    return it;
  }

  template<typename DataType>
  bool Matrix<DataType>::isRegular() const {
    return casadi::isRegular(data_);
  }

  template<typename DataType>
  bool Matrix<DataType>::isSmooth() const {
    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isSymbolic() const {
    return false;
  }

  template<typename DataType>
  bool Matrix<DataType>::isSymbolicSparse() const {
    return false;
  }

  template<typename DataType>
  bool Matrix<DataType>::isInteger() const {
    // loop over non-zero elements
    for (int k=0; k<size(); ++k)
      if (!casadi_limits<DataType>::isInteger(at(k))) // if an element is not integer
        return false;

    // Integer if reached this point
    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isConstant() const {
    // loop over non-zero elements
    for (int k=0; k<size(); ++k)
      if (!casadi_limits<DataType>::isConstant(at(k))) // if an element is not constant
        return false;

    // Constant if we reach this point
    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isOne() const {
    if (!isDense()) {
      return false;
    }

    // loop over non-zero elements
    for (int el=0; el<size(); ++el)
      if (!casadi_limits<DataType>::isOne(at(el)))
        return false;

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isMinusOne() const {
    if (!isDense()) {
      return false;
    }

    // loop over non-zero elements
    for (int el=0; el<size(); ++el)
      if (!casadi_limits<DataType>::isMinusOne(at(el)))
        return false;

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isZero() const {

    // loop over (potentially) non-zero elements
    for (int el=0; el<size(); ++el)
      if (!casadi_limits<DataType>::isZero(at(el)))
        return false;

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isIdentity() const {

    // Make sure that the matrix is diagonal
    if (!sparsity().isDiagonal())
      return false;

    // Make sure that all entries are one
    for (typename Matrix<DataType>::const_iterator it=begin(); it!=end(); ++it) {
      if (!casadi_limits<DataType>::isOne(*it))
        return false;
    }

    return true;
  }

  template<typename DataType>
  bool Matrix<DataType>::isEqual(const Matrix<DataType> &ex2) const {
    // TODO(Joel): Very inefficient, refactor
    if ((size()!=0 || ex2.size()!=0) && shape()!=ex2.shape()) return false;
    Matrix<DataType> difference = *this - ex2;
    return difference.isZero();
  }

  template<typename DataType>
  bool Matrix<DataType>::hasNonStructuralZeros() const {
    // Check if the structural nonzero is known to be zero
    for (int el=0; el<size(); ++el) {
      if (casadi_limits<DataType>::isZero(at(el)))
        return true;
    }

    // No known zeros amongst the structurally nonzero entries
    return false;
  }

  template<typename DataType>
  double Matrix<DataType>::getValue() const {
    return static_cast<double>(toScalar());
  }

  template<typename DataType>
  std::string Matrix<DataType>::getName() const {
    throw CasadiException("\"getName\" not defined for instantiation");
  }

  template<typename DataType>
  Matrix<DataType> Matrix<DataType>::setSparse(const Sparsity& sp, bool intersect) const {
    if (intersect) {
      return setSparse(sp.patternIntersection(sparsity()), false);
    } else {
      Matrix<DataType> ret(sp);
      ret.set(*this);
      return ret;
    }
  }

  template<typename DataType>
  void Matrix<DataType>::setMaxNumCallsInPrint(long num) {
    throw CasadiException("\"setMaxNumCallsInPrint\" not defined for instantiation");
  }

  template<typename DataType>
  long Matrix<DataType>::getMaxNumCallsInPrint() {
    throw CasadiException("\"getMaxNumCallsInPrint\" not defined for instantiation");
  }

  template<typename DataType>
  void Matrix<DataType>::setEqualityCheckingDepth(int eq_depth) {
    throw CasadiException("\"setEqualityCheckingDepth\" not defined for instantiation");
  }

  template<typename DataType>
  int Matrix<DataType>::getEqualityCheckingDepth() {
    throw CasadiException("\"getEqualityCheckingDepth\" not defined for instantiation");
  }


} // namespace casadi

/// \endcond

#endif // CASADI_MATRIX_IMPL_HPP

