%exception  CasADi::Function::call(const MX &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::call(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::call(const std::vector< std::vector< MX > > &arg, const Dictionary &paropt=Dictionary()) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::eval(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::eval(const SX &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::eval(const SXVector &arg, std::vector< SX > &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::eval(const std::vector< MX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::eval(const std::vector< SX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::evalMX(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::evalMX(const std::vector< MX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::evalSX(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::evalSX(const std::vector< SX > &arg) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::indexed_one_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Function::indexed_zero_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::SXFunction::indexed_one_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::SXFunction::indexed_zero_based(int k) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Sparsity::createDiagonal(int nrow) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Sparsity::createDiagonal(int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Sparsity::getSparsity(std::vector< int > &output_row, std::vector< int > &output_col) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Sparsity::getSparsityCCS(std::vector< int > &output_colind, std::vector< int > &output_row) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Sparsity::getSparsityCRS(std::vector< int > &output_rowind, std::vector< int > &output_col) const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Sparsity::lowerNZ() const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::Sparsity::upperNZ() const  {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::densify(const MX &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::densify(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::getIntValue(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::getNZDense(const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::getName(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::getValue(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::hasNonStructuralZeros(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isConstant(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isDense(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isDense(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isEmpty(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isIdentity(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isIdentity(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isInteger(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isMinusOne(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isMinusOne(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isOne(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isOne(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isRegular(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isRegular(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isRegular(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isRegular(const SXElement &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isScalar(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isSingular(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isSmooth(const SX &ex) {
 START DEPRECATED_MSG("use ex.isSmooth()") $action STOP { $action } 
}
%exception  CasADi::isSymbolic(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isSymbolic(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isSymbolicSparse(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isSymbolicSparse(const SX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isTranspose(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isTril(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isTriu(const Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isVector(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isZero(const MX &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::isZero(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::lowerNZ(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::lowerSparsity(const Sparsity &a, bool includeDiagonal=true) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::makeDense(MX &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::makeDense(Matrix< DataType > &A) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::makeSparse(Matrix< DataType > &A, double tol=0) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const Matrix< double > &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const std::string &name, const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const std::string &name, const Sparsity &sp, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const std::string &name, const Sparsity &sp, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const std::string &name, const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const std::string &name, int nrow, int ncol, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const std::string &name, int nrow, int ncol, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::msym(const std::string &name, int nrow=1, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::nnz(const Matrix< DataType > &ex) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_band(int n, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_banded(int n, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_compress(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_compress(const int *v) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_compress(const std::vector< int > &v) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_dense(const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_dense(int nrow, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_diag(int n) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_rowcol(const std::vector< int > &row, const std::vector< int > &col, int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_sparse(const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_sparse(int nrow, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_tril(int n) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_triplet(int nrow, int ncol, const std::vector< int > &row, const std::vector< int > &col) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_triplet(int nrow, int ncol, const std::vector< int > &row, const std::vector< int > &col, std::vector< int > &mapping, bool invert_mapping=false) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_triu(int n) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::sp_unit(int n, int el) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const Matrix< double > &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const std::string &name, const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const std::string &name, const Sparsity &sp, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const std::string &name, const Sparsity &sp, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const std::string &name, const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const std::string &name, int nrow, int ncol, int p) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const std::string &name, int nrow, int ncol, int p, int r) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::ssym(const std::string &name, int nrow=1, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::trans(const MX &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::trans(const Matrix< DataType > &x) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::trans(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::upperNZ(const Sparsity &a) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception  CasADi::upperSparsity(const Sparsity &a, bool includeDiagonal=true) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::IpoptSolver::IpoptSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::KnitroSolver::KnitroSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::MX::MX(const std::string &name, const Sparsity &sp) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::MX::MX(const std::string &name, const std::pair< int, int > &rc) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::MX::MX(const std::string &name, int nrow=1, int ncol=1) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::MX::MX(int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::MX::MX(int nrow, int ncol, const MX &val) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(int nrow, int ncol) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(int nrow, int ncol, const DataType &val) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(int nrow, int ncol, const std::vector< int > &colind, const std::vector< int > &row, const std::vector< DataType > &d=std::vector< DataType >()) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::SCPgen::SCPgen(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::SQPMethod::SQPMethod(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::SXElement::SXElement(const std::string &name) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::SnoptSolver::SnoptSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::Sparsity::Sparsity(int nrow, int ncol, bool dense=false) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::StabilizedSQPMethod::StabilizedSQPMethod(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}
%exception CasADi::WorhpSolver::WorhpSolver(const Function &F, const Function &G) {
 START DEPRECATED_MSG("") $action STOP { $action } 
}