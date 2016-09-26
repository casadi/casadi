%exception  casadi::BonminUserClass::branchingInfo() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BonminUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, TNLP::IndexStyleEnum &index_style) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BonminUserClass::get_number_of_nonlinear_variables() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BonminUserClass::sosConstraints() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::compile(const std::string &compiler="gcc -fPIC -O2") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::tic() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::toc() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_input(const std::string &s, const MatType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_output(const std::string &s, const MatType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::get_input(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::get_output(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::has_in(const std::string &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::has_out(const std::string &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::name_in() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::name_out() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::request_input(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::request_output(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
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
%exception  casadi::Matrix< Scalar >::T() const  {
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
%exception  casadi::Matrix< Scalar >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::cofactor(const Matrix< Scalar > &x, int i, int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::eig_symbolic(const Matrix< Scalar > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::enlarge(int nrow, int ncol, const std::vector< int > &rr, const std::vector< int > &cc, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::erase(const std::vector< int > &rr, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::erase(const std::vector< int > &rr, const std::vector< int > &cc, bool ind1=false) {
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
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &output_m, bool ind1, const Matrix< int > &rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &output_m, bool ind1, const Matrix< int > &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &output_m, bool ind1, const Matrix< int > &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &output_m, bool ind1, const Slice &rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &output_m, bool ind1, const Slice &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &output_m, bool ind1, const Slice &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &output_m, bool ind1, const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::getMinor(const Matrix< Scalar > &x, int i, int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_nonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_nz(Matrix< Scalar > &output_m, bool ind1, const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_nz(Matrix< Scalar > &output_m, bool ind1, const Slice &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::has_zeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::heaviside(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inf(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inf(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inf(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_constant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_identity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_integer() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_minus_one() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_one() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_zero() const  {
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
%exception  casadi::Matrix< Scalar >::nan(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nan(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nan(int nrow=1, int ncol=1) {
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
%exception  casadi::Matrix< Scalar >::print(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_dense(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_scalar(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_sparse(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_vector(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::printme(const Matrix< Scalar > &y) const  {
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
%exception  casadi::Matrix< Scalar >::remove(const std::vector< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::repr(std::ostream &stream=casadi::userOut(), bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reserve(int nnz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reserve(int nnz, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::resize(int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::sanity_check(bool complete=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::scalar_matrix(int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< int > &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< int > &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< int > &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Slice &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Slice &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Slice &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set_nz(const Matrix< Scalar > &m, bool ind1, const Matrix< int > &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set_nz(const Matrix< Scalar > &m, bool ind1, const Slice &k) {
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
%exception  casadi::Matrix< Scalar >::triplet(const std::vector< int > &row, const std::vector< int > &col, const Matrix< Scalar > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triplet(const std::vector< int > &row, const std::vector< int > &col, const Matrix< Scalar > &d, const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triplet(const std::vector< int > &row, const std::vector< int > &col, const Matrix< Scalar > &d, int nrow, int ncol) {
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
%exception  casadi::SX::dep(int ch=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::element_hash() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::grad(const Function &f, const std::string &iname, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::grad(const Function &f, const std::string &iname, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::grad(const Function &f, int iind, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::grad(const Function &f, int iind=0, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::has_duplicates() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::hess(const Function &f, const std::string &iname, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::hess(const Function &f, const std::string &iname, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::hess(const Function &f, int iind, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::hess(const Function &f, int iind=0, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_commutative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_leaf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_regular() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_smooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_symbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_valid_input() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::jac(const Function &f, const std::string &iname, const std::string &oname, bool compact=false, bool symmetric=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::jac(const Function &f, const std::string &iname, int oind=0, bool compact=false, bool symmetric=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::jac(const Function &f, int iind, const std::string &oname, bool compact=false, bool symmetric=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::jac(const Function &f, int iind=0, int oind=0, bool compact=false, bool symmetric=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::n_dep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::print_split(std::vector< std::string > &output_nz, std::vector< std::string > &output_inter) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::resetInput() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::tang(const Function &f, const std::string &iname, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::tang(const Function &f, const std::string &iname, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::tang(const Function &f, int iind, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::tang(const Function &f, int iind=0, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::printPtr(std::ostream &stream=casadi::userOut()) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::has_nz(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::nonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::nonzeros() {
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
%exception  casadi::SparseStorage< DataType >::scalar() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sparsityRef() {
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
%exception  casadi::combine(const Dict &first, const Dict &second) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_regular(N_Vector v) {
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
%exception  casadi::zip(const std::vector< std::string > &id, const std::vector< T > &mat) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CodeGenerator::CodeGenerator(const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FStats::FStats() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Factory< MatType >::Factory(const Function::AuxOut &aux) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Matrix< Scalar > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const Matrix< Scalar > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::vector< std::vector< double > > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< T >::Matrix(const Sparsity &sp, const Scalar &val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< T >::Matrix(const Sparsity &sp, const std::vector< Scalar > &d, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< T >::Matrix(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< T >::Matrix(const std::vector< Scalar > &x) {
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