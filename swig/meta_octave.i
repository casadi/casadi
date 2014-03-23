/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

%inline %{

/// std::vector<double>
template<> char meta< std::vector<double> >::expected_message[] = "Expecting (1xn) array(number)";

template <>
int meta< std::vector<double> >::as(const octave_value& p, std::vector<double> &m) {
  NATIVERETURN(std::vector<double>, m);
  if(p.is_real_matrix() && p.is_numeric_type()){
    const Matrix &mat = p.matrix_value();
    if (mat.cols()==1) {
      m.resize(mat.rows());
      for(int i=0; i<mat.cols(); ++i) m[i] = mat(i,0);
    } else if (mat.rows()==1) {
      m.resize(mat.cols());
      for(int j=0; j<mat.cols(); ++j) m[j] = mat(0,j);
    } else {
      return false;
    }
    return true;
  }
}

template <> bool meta< std::vector<double> >::couldbe(const octave_value& p) { 
  if (meta< std::vector<double> >::isa(p)) return true; 
  if(p.is_real_matrix() && p.is_numeric_type()){
    const Matrix &mat = p.matrix_value();
    return (mat.rows()==1 || mat.cols()==1);
  } else {
    return false;
  }
}

/// std::vector<int>
template<> char meta< std::vector<int> >::expected_message[] = "Expecting (1xn) array(number)";

template <>
int meta< std::vector<int> >::as(const octave_value& p, std::vector<int> &m) {
  NATIVERETURN(std::vector<int>, m);
  if(p.is_real_matrix()  && p.is_numeric_type()){
    const Matrix &mat = p.matrix_value();
    if (mat.cols()==1) {
      m.resize(mat.rows());
      for(int i=0; i<mat.cols(); ++i) m[i] = mat(i,0);
    } else if (mat.rows()==1) {
      m.resize(mat.cols());
      for(int j=0; j<mat.cols(); ++j) m[j] = mat(0,j);
    } else {
      return false;
    }
    return true;
  }
}

template <> bool meta< std::vector<int> >::couldbe(const octave_value& p) { 
  if (meta< std::vector<int> >::isa(p)) return true; 
  if(p.is_real_matrix() && p.is_numeric_type()) {
    const Matrix &mat = p.matrix_value();
    return (mat.rows()==1 || mat.cols()==1);
  } else {
    return false;
  }
}


/// CasADi::GenericType
template<> char meta< CasADi::GenericType >::expected_message[] = "Expecting number, string, vector(number)";

template <>
int meta< CasADi::GenericType >::as(const octave_value& p,CasADi::GenericType &s) {
  NATIVERETURN(CasADi::GenericType, s)
  if (p.is_real_scalar()) {
    s=CasADi::GenericType(p.double_value());
  } else if (meta< std::vector<double> >::couldbe(p)) {
    std::vector<double> temp;
    int ret = meta< std::vector<double> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (meta< std::vector<int> >::couldbe(p)) {
    std::vector<int> temp;
    int ret = meta< std::vector<int> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (p.is_string()) {
    s = CasADi::GenericType(p.string_value());
  } else {
    return false;
  }
  return true;
}

template <>
bool meta< CasADi::GenericType >::couldbe(const octave_value& p) {
  return p.is_real_scalar() || meta< std::vector<int> >::couldbe(p) || meta< std::vector<double> >::couldbe(p) || p.is_string() ;
}

/// CasADi::Matrix<double>
template<> char meta< CasADi::Matrix<double> >::expected_message[] = "Expecting numpy.array2D, numpy.matrix, csc_matrix, DMatrix";

template <>
int meta< CasADi::Matrix<double> >::as(const octave_value& p,CasADi::Matrix<double> &m) {
  NATIVERETURN(CasADi::Matrix<double>,m)
  if((p.is_real_matrix() && p.is_numeric_type() && p.is_sparse_type())){
    SparseMatrix mat = p.sparse_matrix_value();
    
    int size = mat.nnz();
    
    std::vector<double> data(size);
    for (int k=0;k<data.size();k++) data[k]=mat.data(k);

    std::vector<int> cidx(mat.cols()+1);
    std::vector<int> ridx(size);
    for (int k=0;k<cidx.size();k++) cidx[k]=mat.cidx(k);
    for (int k=0;k<ridx.size();k++) ridx[k]=mat.ridx(k);
    
    CasADi::Sparsity A = CasADi::Sparsity(mat.rows(),mat.cols(),cidx,ridx);
    m = CasADi::Matrix<double>(A,data);
    
    return true;
  }
  if((p.is_real_matrix() && p.is_numeric_type())){
    Matrix mat = p.matrix_value();
    m = CasADi::DMatrix::zeros(mat.rows(),mat.cols());
    for(int i=0; i<mat.rows(); ++i){
      for(int j=0; j<mat.cols(); ++j){
        m(i,j) = mat(i,j);
      }
    }
    return true;
  }
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m = CasADi::DMatrix(p.double_value());
    return true;
  } 
  return false;
}

// Disallow 1D numpy arrays. Allowing them may introduce conflicts with other typemaps or overloaded methods
template <>
bool meta< CasADi::Matrix<double> >::couldbe(const octave_value& p) { return meta< CasADi::Matrix<double> >::isa(p) || (p.is_real_matrix() && p.is_numeric_type()) || (p.is_real_scalar() && p.is_numeric_type());}


// Explicit intialization of these two member functions, so we can use them in meta< CasADi::SXElement >
template<> int meta< CasADi::Matrix<CasADi::SXElement> >::as(GUESTOBJECT,CasADi::Matrix<CasADi::SXElement> &);
template<> bool meta< CasADi::Matrix<CasADi::SXElement> >::couldbe(GUESTOBJECT);

/// CasADi::SX
template<> char meta< CasADi::SXElement >::expected_message[] = "Expecting SXElement or number";

template <>
int meta< CasADi::SXElement >::as(const octave_value& p,CasADi::SXElement &s) {
  NATIVERETURN(CasADi::SXElement, s)
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    s=CasADi::SXElement(p.double_value());
    return true;
  } else if (meta< CasADi::Matrix< CasADi::SXElement > >::isa(p)) {
    CasADi::Matrix< CasADi::SXElement > m;
    meta< CasADi::Matrix< CasADi::SXElement > >::as(p,m);
    if (m.numel()==1 && m.size()==1) {
      s = m.at(0);
      return true;
    }
  }
  return false;
}

template <>
bool meta< CasADi::SXElement >::couldbe(const octave_value& p) {
  if (meta< CasADi::Matrix< CasADi::SXElement > >::isa(p)) {
    CasADi::Matrix< CasADi::SXElement > m;
    meta<CasADi::Matrix< CasADi::SXElement > >::as(p,m);
    if (m.numel()==1 && m.size()==1)
      return true;
  }
  return (meta< CasADi::SXElement >::isa(p) || (p.is_real_scalar() && p.is_numeric_type()) );
}


/// CasADi::Matrix<CasADi::SXElement>
template<> char meta< CasADi::Matrix<CasADi::SXElement> >::expected_message[] = "Expecting one of: numpy.ndarray(SX/number) , SX, SX, number, sequence(SX/number)";

template <>
int meta< CasADi::Matrix<CasADi::SXElement> >::as(const octave_value& p,CasADi::Matrix<CasADi::SXElement> &m) {
  NATIVERETURN(CasADi::Matrix<CasADi::SXElement>, m)
  NATIVERETURN(CasADi::Matrix<double>, m)
  NATIVERETURN(CasADi::SXElement, m)
  if((p.is_real_matrix() && p.is_numeric_type())){
    Matrix mat = p.matrix_value();
    m = CasADi::SX::zeros(mat.rows(),mat.cols());
    for(int i=0; i<mat.rows(); ++i){
      for(int j=0; j<mat.cols(); ++j){
        m(i,j) = mat(i,j);
      }
    }
    return true;
  } 
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m = CasADi::SXElement(p.double_value());
    return true;
  }
  return false;
}

template <> bool meta< CasADi::Matrix<CasADi::SXElement> >::couldbe(const octave_value& p) { return meta< CasADi::Matrix<CasADi::SXElement> >::isa(p) || meta< CasADi::SXElement >::isa(p) || meta< CasADi::Matrix<double> >::isa(p)  || (p.is_real_matrix() && p.is_numeric_type()) || (p.is_real_scalar() && p.is_numeric_type());}


meta_vector(CasADi::Matrix<CasADi::SXElement>);
meta_vector(std::vector< CasADi::Matrix<CasADi::SXElement> >);
meta_vector(std::vector<CasADi::SXElement>);
meta_vector(CasADi::SXElement);


/// CasADi::MX
template<> char meta< CasADi::MX >::expected_message[] = "Expecting (MX, numberarray)";

template <>
bool meta< CasADi::MX >::couldbe(const octave_value& p) {
  return (meta< CasADi::MX >::isa(p) || meta< CasADi::Matrix<double> >::couldbe(p) );
}

template <>
int meta< CasADi::MX >::as(const octave_value& p,CasADi::MX &m) {
  NATIVERETURN(CasADi::MX,m)
  NATIVERETURN(CasADi::Matrix<double>,m)
  if(meta< CasADi::Matrix<double> >::couldbe(p)) {
    CasADi::DMatrix mt;
    bool result=meta< CasADi::Matrix<double> >::as(p,mt);
    if (!result)
      return false;
    m = CasADi::MX(mt);
    return true;
  }
  return false;
}

meta_vector(CasADi::MX);
%}

