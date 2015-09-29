%exception  casadi::Adaptor< Derived, Solver >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::inputD(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::inputD(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::outputD(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::outputD(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::readInputs() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseIO< Derived >::writeOutputs() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::callDerivative(const DMatrixVector &arg, DMatrixVector &output_res, const DMatrixVectorVector &fseed, DMatrixVectorVector &output_fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::callDerivative(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::callDerivative(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spInit(bool fwd) {
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
%exception  casadi::IOInterface< Function  >::getInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::getInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, int oind=0) {
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
%exception  casadi::LibInfo< Compiler >::name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LibInfo< std::string >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LibInfo< std::string >::get(FcnPtr &fcnPtr, const std::string &sym) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LibInfo< std::string >::name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getTemp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::hasDuplicates() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::resetInput() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXFunction::generateLiftingFunctions(MXFunction &output_vdef_fcn, MXFunction &output_vinit_fcn) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::binary(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::hasDuplicates() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_scalar(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::resetInput() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::scalar_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::unary(int op, const Matrix< DataType > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::getOptionAllowedIndex(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::getOptionEnumValue(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::setOptionByAllowedIndex(const std::string &name, int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionality::setOptionByEnumValue(const std::string &name, int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_ENTRY >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_EXIT >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_IO >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_NAME >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_SOURCE >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProfilingType< ProfilingData_TIMELINE >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::assertInit() const  {
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
%exception  casadi::SparseStorage< DataType >::toScalar() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::reCache() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_tril(bool includeDiagonal=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_triu(bool includeDiagonal=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::alive() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::shared() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< Derived >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< Derived >::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::check_exposed(T t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::diffTimers(const timer t1, const timer t0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::diffToDict(const diffTime &diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getRealTime() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getTimerTime(void) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_combine(std::size_t &seed, T v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_combine(std::size_t &seed, const std::vector< int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_sparsity(int nrow, int ncol, const std::vector< int > &colind, const std::vector< int > &row) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::hash_value(T v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::iszero(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::operation_checker(unsigned int op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::profileWrite(std::ofstream &f, const T &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::profileWriteBare(std::ofstream &f, const T &s) {
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
%exception  casadi::slicot_periodic_schur(const std::vector< Matrix< double > > &a, std::vector< Matrix< double > > &t, std::vector< Matrix< double > > &z, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &dwork, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::timerPlusEq(diffTime &t, const diffTime diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblem::setLog(isnLog snLog, isnLog2 snLog2, isqLog sqLog) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblem::setSTOP(isnSTOP snSTOP) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblem::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::computeJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setNeA(int neA) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setNeG(int neG) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setObjective(int ObjRow, double ObjAdd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setProblemSize(int n, int neF) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setUserFun(snFunA usrfun) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::setWorkspace() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemA::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemB::setFuncon(snConB funcon) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemB::setFunobj(snObjB funobj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemB::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setObjective(int iObj, double ObjAdd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setProblemSize(int m, int n, int nnCon, int nnJac, int nnObj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setUserFun(snFunC usrfun) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::setWorkspace() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  snoptProblemC::solve(int starttype) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< Compiler >::LibInfo() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< Compiler >::LibInfo(const Compiler &compiler, const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< std::string >::LibInfo() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LibInfo< std::string >::LibInfo(const std::string &bin_name, const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinearSolver::LinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sp, const DataType &val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sp, const std::vector< DataType > &d, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x) {
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
%exception snoptProblemA::snoptProblemA() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception snoptProblemB::snoptProblemB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception snoptProblemC::snoptProblemC() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}