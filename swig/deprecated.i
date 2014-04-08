%exception  casadi::Function::call(const MX &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::call(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::call(const std::vector< std::vector< MX > > &arg, const Dictionary &paropt=Dictionary()) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::eval(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::eval(const SX &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::eval(const SXVector &arg, std::vector< SX > &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::eval(const std::vector< MX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::eval(const std::vector< SX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::evalMX(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::evalMX(const std::vector< MX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::evalSX(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::evalSX(const std::vector< SX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::indexed_one_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Function::indexed_zero_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::SXFunction::indexed_one_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::SXFunction::indexed_zero_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Sparsity::createDiagonal(int nrow) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Sparsity::createDiagonal(int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Sparsity::getSparsity(std::vector< int > &output_row, std::vector< int > &output_col) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Sparsity::getSparsityCCS(std::vector< int > &output_colind, std::vector< int > &output_row) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Sparsity::getSparsityCRS(std::vector< int > &output_rowind, std::vector< int > &output_col) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Sparsity::lowerNZ() const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::Sparsity::upperNZ() const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::densify(const MX &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::densify(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::getIntValue(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::getNZDense(const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::getName(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::getValue(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::hasNonStructuralZeros(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isConstant(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isDense(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isDense(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isEmpty(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isIdentity(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isIdentity(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isInteger(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isMinusOne(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isMinusOne(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isOne(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isOne(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isRegular(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isRegular(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isRegular(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isRegular(const SXElement &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isScalar(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isSingular(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isSmooth(const SX &ex) {
 START DEPRECATED_MSG("use ex.isSmooth()") $action STOP { $action } 
}
%exception  casadi::isSymbolic(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isSymbolic(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isSymbolicSparse(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isSymbolicSparse(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isTranspose(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isTril(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isTriu(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isVector(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isZero(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::isZero(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::lowerNZ(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::lowerSparsity(const Sparsity &a, bool includeDiagonal=true) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::makeDense(MX &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::makeDense(Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::makeSparse(Matrix< DataType > &A, double tol=0) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const Matrix< double > &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const std::string &name, const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const std::string &name, const Sparsity &sp, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const std::string &name, const Sparsity &sp, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const std::string &name, const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const std::string &name, int nrow, int ncol, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const std::string &name, int nrow, int ncol, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::msym(const std::string &name, int nrow=1, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::nnz(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_band(int n, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_banded(int n, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_compress(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_compress(const int *v) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_compress(const std::vector< int > &v) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_dense(const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_dense(int nrow, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_diag(int n) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_rowcol(const std::vector< int > &row, const std::vector< int > &col, int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_sparse(const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_sparse(int nrow, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_tril(int n) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_triplet(int nrow, int ncol, const std::vector< int > &row, const std::vector< int > &col) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_triplet(int nrow, int ncol, const std::vector< int > &row, const std::vector< int > &col, std::vector< int > &mapping, bool invert_mapping=false) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_triu(int n) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::sp_unit(int n, int el) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const Matrix< double > &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const std::string &name, const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const std::string &name, const Sparsity &sp, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const std::string &name, const Sparsity &sp, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const std::string &name, const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const std::string &name, int nrow, int ncol, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const std::string &name, int nrow, int ncol, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::ssym(const std::string &name, int nrow=1, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::trans(const MX &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::trans(const Matrix< DataType > &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::trans(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::upperNZ(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  casadi::upperSparsity(const Sparsity &a, bool includeDiagonal=true) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::IpoptSolver::IpoptSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::KnitroSolver::KnitroSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::MX::MX(const std::string &name, const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::MX::MX(const std::string &name, const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::MX::MX(const std::string &name, int nrow=1, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::MX::MX(int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::MX::MX(int nrow, int ncol, const MX &val) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(int nrow, int ncol, const DataType &val) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(int nrow, int ncol, const std::vector< int > &colind, const std::vector< int > &row, const std::vector< DataType > &d=std::vector< DataType >()) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::SCPgen::SCPgen(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::SQPMethod::SQPMethod(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::SXElement::SXElement(const std::string &name) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::SnoptSolver::SnoptSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::Sparsity::Sparsity(int nrow, int ncol, bool dense=false) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::StabilizedSQPMethod::StabilizedSQPMethod(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception casadi::WorhpSolver::WorhpSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}