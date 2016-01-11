%exception  casadi::Function::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::derivative(const DMVector &arg, DMVector &output_res, const DMVectorVector &fseed, DMVectorVector &output_fsens, const DMVectorVector &aseed, DMVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::derivative(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::derivative(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_arg() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_iw() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_res() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_w() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< SX >::sym(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::finalize_metadata(Index n, const StringMetaDataMapType &var_string_md, const IntegerMetaDataMapType &var_integer_md, const NumericMetaDataMapType &var_numeric_md, Index m, const StringMetaDataMapType &con_string_md, const IntegerMetaDataMapType &con_integer_md, const NumericMetaDataMapType &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_number_of_nonlinear_variables() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_var_con_metadata(Index n, StringMetaDataMapType &var_string_md, IntegerMetaDataMapType &var_integer_md, NumericMetaDataMapType &var_numeric_md, Index m, StringMetaDataMapType &con_string_md, IntegerMetaDataMapType &con_integer_md, NumericMetaDataMapType &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LibInfo< Compiler >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LibInfo< Compiler >::get(FcnPtr &fcnPtr, const std::string &sym) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LibInfo< std::string >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LibInfo< std::string >::get(FcnPtr &fcnPtr, const std::string &sym) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getTemp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::has_duplicates() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::resetInput() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::adj(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::all(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::any(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::binary(int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::chol(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::cofactor(const Matrix< Scalar > &x, int i, int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::eig_symbolic(const Matrix< Scalar > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::expand(const Matrix< Scalar > &ex, Matrix< Scalar > &weights, Matrix< Scalar > &terms) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, int order, const Matrix< Scalar > &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, int order=5) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::getMinor(const Matrix< Scalar > &x, int i, int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::has_duplicates() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::heaviside(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::matrix_matrix(int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::matrix_scalar(int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mtaylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, int order, const std::vector< int > &order_contributions) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mtaylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, int order=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::norm_inf_mul(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::poly_coeff(const Matrix< Scalar > &f, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::poly_roots(const Matrix< Scalar > &p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::pw_const(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::pw_lin(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::qr(const Matrix< Scalar > &A, Matrix< Scalar > &Q, Matrix< Scalar > &R) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::ramp(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::rectangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::resetInput() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::scalar_matrix(int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::sparsify(const Matrix< Scalar > &A, double tol=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::taylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, int order=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::unary(int op, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< double >::pinv(const Matrix< double > &A, const std::string &lsolver, const Dict &dict) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< double >::solve(const Matrix< double > &a, const Matrix< double > &b, const std::string &lsolver, const Dict &dict) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::optionAllowedIndex(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::optionEnumValue(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::setOptionByAllowedIndex(const std::string &name, int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::setOptionByEnumValue(const std::string &name, int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::printPtr(std::ostream &stream=casadi::userOut()) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::data() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::data() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::hasNZ(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::resize(int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::scalar() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::alive() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::shared() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::check_exposed(T t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::diffTimers(const Timer t1, const Timer t0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::diffToDict(const DiffTime &diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getTimerTime(void) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_regular(N_Vector v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::iszero(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName< SXElem >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(const std::vector< std::vector< T > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ptrVec(std::vector< std::vector< T > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::replaceMat(const M &arg, const Sparsity &inp, bool hcat=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::timerPlusEq(DiffTime &t, const DiffTime diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::zip(const std::vector< std::string > &id, const std::vector< T > &mat) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IpoptUserClass::IpoptUserClass(const IpoptInterface &ipoptInterface, IpoptMemory &mem) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< Compiler >::LibInfo() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< Compiler >::LibInfo(const Compiler &compiler) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< std::string >::LibInfo() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< std::string >::LibInfo(const std::string &bin_name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const Scalar &val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const std::vector< Scalar > &d, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::vector< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const DataType &val=DataType(0)) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(SharedObject shared) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(int dummy=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}