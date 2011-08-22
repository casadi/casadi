/// std::vector<double>
%inline %{
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
  NATIVECOULDBE(std::vector<double>)
  if(p.is_real_matrix() && p.is_numeric_type()){
    const Matrix &mat = p.matrix_value();
    return (mat.rows()==1 || mat.cols()==1);
  } else {
    return false;
  }
}

%}

/// std::vector<int>
%inline %{
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
  NATIVECOULDBE(std::vector<int>)
  if(p.is_real_matrix() && p.is_numeric_type()) {
    const Matrix &mat = p.matrix_value();
    return (mat.rows()==1 || mat.cols()==1);
  } else {
    return false;
  }
}

%}


/// CasADi::GenericType
%inline %{

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

%}

/// CasADi::Matrix<double>
%inline %{
template<> char meta< CasADi::Matrix<double> >::expected_message[] = "Expecting numpy.array2D, numpy.matrix, csr_matrix, DMatrix";

template <>
int meta< CasADi::Matrix<double> >::as(const octave_value& p,CasADi::Matrix<double> &m) {
  NATIVERETURN(CasADi::Matrix<double>,m)
  if((p.is_real_matrix() && p.is_numeric_type() && p.is_sparse_type())){
    // Note: octave uses column-major storage
    SparseMatrix mat = p.sparse_matrix_value();
    
    int size = mat.nnz();
    
    std::vector<double> data(size);
    for (int k=0;k<data.size();k++) data[k]=mat.data(k);

    std::vector<int> cidx(mat.cols()+1);
    std::vector<int> ridx(size);
    for (int k=0;k<cidx.size();k++) cidx[k]=mat.cidx(k);
    for (int k=0;k<ridx.size();k++) ridx[k]=mat.ridx(k);
    
    CasADi::CRSSparsity A = CasADi::CRSSparsity(mat.cols(),mat.rows(),ridx,cidx);
    CasADi::Matrix<double> ret = CasADi::Matrix<double>(A,data);
    
    m = ret.trans();
    
    return true;
  }
  if((p.is_real_matrix() && p.is_numeric_type())){
    Matrix mat = p.matrix_value();
    m = CasADi::DMatrix(mat.rows(),mat.cols(),0);
    for(int i=0; i<mat.rows(); ++i){
      for(int j=0; j<mat.cols(); ++j){
        m(i,j) = mat(i,j);
      }
    }
    return true;
  }
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m = CasADi::DMatrix(1,1,p.double_value());
    return true;
  } 
  return false;
}

// Disallow 1D numpy arrays. Allowing them may introduce conflicts with other typemaps or overloaded methods
template <>
bool meta< CasADi::Matrix<double> >::couldbe(const octave_value& p) { return meta< CasADi::Matrix<double> >::isa(p) || (p.is_real_matrix() && p.is_numeric_type()) || (p.is_real_scalar() && p.is_numeric_type());}

%}

%inline %{
// Explicit intialization of these two member functions, so we can use them in meta< CasADi::SX >
template<> int meta< CasADi::Matrix<CasADi::SX> >::as(GUESTOBJECT,CasADi::Matrix<CasADi::SX> &);
template<> bool meta< CasADi::Matrix<CasADi::SX> >::couldbe(GUESTOBJECT);
%}

/// CasADi::SX
%inline %{
template<> char meta< CasADi::SX >::expected_message[] = "Expecting SX or number";

template <>
int meta< CasADi::SX >::as(const octave_value& p,CasADi::SX &s) {
  NATIVERETURN(CasADi::SX, s)
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    s=CasADi::SX(p.double_value());
    return true;
  } else if (meta< CasADi::Matrix< CasADi::SX > >::isa(p)) {
    CasADi::Matrix< CasADi::SX > m;
    meta< CasADi::Matrix< CasADi::SX > >::as(p,m);
    if (m.numel()==1 && m.size()==1) {
      s = m.at(0);
      return true;
    }
  }
  return false;
}

template <>
bool meta< CasADi::SX >::couldbe(const octave_value& p) {
  if (meta< CasADi::Matrix< CasADi::SX > >::isa(p)) {
    CasADi::Matrix< CasADi::SX > m;
    meta<CasADi::Matrix< CasADi::SX > >::as(p,m);
    if (m.numel()==1 && m.size()==1)
      return true;
  }
  return (meta< CasADi::SX >::isa(p) || (p.is_real_scalar() && p.is_numeric_type()) );
}

%}


/// CasADi::Matrix<CasADi::SX>
%inline %{
template<> char meta< CasADi::Matrix<CasADi::SX> >::expected_message[] = "Expecting one of: numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number)";

template <>
int meta< CasADi::Matrix<CasADi::SX> >::as(const octave_value& p,CasADi::Matrix<CasADi::SX> &m) {
  NATIVERETURN(CasADi::Matrix<CasADi::SX>, m)
  NATIVERETURN(CasADi::Matrix<double>, m)
  NATIVERETURN(CasADi::SX, m)
  if((p.is_real_matrix() && p.is_numeric_type())){
    Matrix mat = p.matrix_value();
    m = CasADi::SXMatrix(mat.rows(),mat.cols(),0);
    for(int i=0; i<mat.rows(); ++i){
      for(int j=0; j<mat.cols(); ++j){
        m(i,j) = mat(i,j);
      }
    }
    return true;
  } 
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m = CasADi::SX(p.double_value());
    return true;
  }
  return false;
}

template <> bool meta< CasADi::Matrix<CasADi::SX> >::couldbe(const octave_value& p) { return meta< CasADi::Matrix<CasADi::SX> >::isa(p) || meta< CasADi::SX >::isa(p) || meta< CasADi::Matrix<double> >::isa(p)  || (p.is_real_matrix() && p.is_numeric_type()) || (p.is_real_scalar() && p.is_numeric_type());}
%}


/// CasADi::Slice
%inline %{
template<> char meta< CasADi::Slice >::expected_message[] = "Expecting Slice or number";

template <>
int meta< CasADi::Slice >::as(const octave_value& p,CasADi::Slice &m) {
  if (p.is_range()) {
    Range r = p.range_value();
    m.start_ = r.base()-1;
    m.stop_ = r.limit();
    m.step_ = r.inc();
  } else if (p.is_magic_colon()) {
    m.start_ = 0;
    m.stop_ = std::numeric_limits<int>::max();
  } else if (p.is_numeric_type()) {
    m.start_ = p.int_value()-1;
    m.stop_ = m.start_+1;
  } else {
    return false;
  }
  return true;
}

template <>
bool meta<  CasADi::Slice >::couldbe(const octave_value& p) {
  return p.is_range() || p.is_magic_colon()|| (p.is_real_scalar() && p.is_numeric_type());
}

%}


/// CasADi::IndexList
%inline %{
template<> char meta< CasADi::IndexList >::expected_message[] = "Expecting Slice or number or list of ints";

template <>
int meta< CasADi::IndexList >::as(const octave_value& p,CasADi::IndexList &m) {
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m.type = CasADi::IndexList::INT;
    m.i = p.int_value()-1;
  } else if (meta< std::vector<int> >::couldbe(p)) {
    m.type = CasADi::IndexList::IVECTOR;
    bool result = meta< std::vector<int> >::as(p,m.iv);
    if (!result) return false;
    for (int k=0; k < m.iv.size();k++) m.iv[k]--;
  } else if (meta< CasADi::Slice>::couldbe(p)) {
    m.type = CasADi::IndexList::SLICE;
    return meta< CasADi::Slice >::as(p,m.slice);
  } else {
    return false;
  }
  return true;
}


template <>
bool meta<  CasADi::IndexList >::couldbe(const octave_value& p) {
  return meta< CasADi::Slice >::couldbe(p) || meta< std::vector<int> >::couldbe(p) || (p.is_real_scalar() && p.is_numeric_type());
}
%}


/// CasADi::MX
%inline %{
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
%}
