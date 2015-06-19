%exception  casadi::Adaptor< Derived, Solver >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addAuxiliary(Auxiliary f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addDependency(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addInclude(const std::string &new_include, bool relative_path=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::addSparsity(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::casadi_dot(int n, const std::string &x, int inc_x, const std::string &y, int inc_y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::copyVector(std::ostream &s, const std::string &arg, std::size_t n, const std::string &res, const std::string &it="i", bool only_if_exists=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::flush(std::ostream &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< double > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< int > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getDependency(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::getSparsity(const Sparsity &sp) const  {
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
%exception  casadi::Function::checkAdjSeed(const std::vector< std::vector< M > > &aseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkArg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkFwdSeed(const std::vector< std::vector< M > > &fseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkRes(const std::vector< M > &res) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::inputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::inputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::input_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::input_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingAdjSeed(const std::vector< std::vector< M > > &aseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingArg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingFwdSeed(const std::vector< std::vector< M > > &fseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::matchingRes(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::outputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::outputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::output_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::output_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceAdjSeed(const std::vector< std::vector< M > > &aseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceArg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceFwdSeed(const std::vector< std::vector< M > > &fseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::replaceRes(const std::vector< M > &res) const  {
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
%exception  casadi::GenericExpression< SXElement  >::zz_ge(const SXElement &y) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElement  >::zz_gt(const SXElement &y) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_a() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDictionary() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDoubleVector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVectorVector() {
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
%exception  casadi::IOInterface< Function  >::inputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::inputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::inputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::outputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Function  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOScheme::print(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOScheme::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::vector() const {
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
%exception  casadi::LinearSolver::spSolve(DMatrix &X, const DMatrix &B, bool transpose=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::at(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::at(int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getTemp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXFunction::algorithm() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXFunction::generateLiftingFunctions(MXFunction &vdef_fcn, MXFunction &vinit_fcn) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::at(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::at(int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::back() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::back() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::begin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::begin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::binary(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::borBV(const Matrix< DataType > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::data() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::data() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::end() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::end() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::front() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::front() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(double &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getBV(Matrix< DataType > &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getElement(int rr, int cc=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNZ(double &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::makeDense(const DataType &val=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_scalar(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::ptr() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::ptr() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::rbegin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::rbegin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::rend() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::rend() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::scalar_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setAll(const DataType &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setBV(const Matrix< DataType > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setZeroBV() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::toScalar() const  {
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
%exception  casadi::PrintableObject< SXElement  >::getDescription() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PrintableObject< SXElement  >::getRepresentation() const {
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
%exception  casadi::SXElement::assignIfDuplicate(const SXElement &scalar, int depth=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::assignNoDelete(const SXElement &scalar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::constpow(const Matrix< SXElement > &n) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::constpow(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::get() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getDep(int ch=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getNdeps() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getTemp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::inv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isCommutative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isDoubled() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isLeaf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isMinusInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isNan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isNonNegative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isNull() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isOp(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isRegular() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::mark() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::marked() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::printme(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::sq() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_abs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_acos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_acosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_and(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_asin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_asinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atan2(const Matrix< SXElement > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atan2(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_atanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_ceil() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_cos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_cosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_eq(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_erf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_erfinv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_exp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_floor() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_if_else_zero(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_isEqual(const SXElement &scalar, int depth=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_le(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_log() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_log10() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_lt(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_max(const Matrix< SXElement > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_max(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_min(const Matrix< SXElement > &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_min(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_minus(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_mod(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_mpower(const SXElement &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_mul(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_ne(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_not() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_or(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_plus(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_power(const SXElement &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_rdivide(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sign() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_simplify() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_sqrt() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_tan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_tanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::zz_times(const SXElement &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXFunction::algorithm() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::assertInit() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::get() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::get() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::getCount() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::printPtr(std::ostream &stream=CASADI_COUT) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::swap(SharedObject &other) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::weak() {
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
%exception  casadi::Sparsity::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::find(std::vector< int > &loc, bool ind1=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::reCache() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_blockcat(const std::vector< std::vector< Sparsity > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_diagcat(const std::vector< Sparsity > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_diagsplit(const std::vector< int > &offset1, const std::vector< int > &offset2) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_horzcat(const std::vector< Sparsity > &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_horzsplit(const std::vector< int > &output_offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_kron(const Sparsity &b) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_mtimes(const Sparsity &Y, const Sparsity &Z) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_mtimes(const Sparsity &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_norm_0_mul(const Sparsity &B) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_reshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_reshape(int nrow, int ncol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_sprank() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_tril(bool includeDiagonal=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_triu(bool includeDiagonal=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_vecNZ() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_vertcat(const std::vector< Sparsity > &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::zz_vertsplit(const std::vector< int > &output_offset) const  {
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
%exception  casadi::XmlFile::parse(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::acosh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::asinh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atanh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::check_exposed(T t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocationInterpolators(const std::vector< double > &tau_root, std::vector< std::vector< double > > &C, std::vector< double > &D) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocationPointsL(int order, const std::string &scheme="radau") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::constpow(const T &x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::copysign(const T &x, const T &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::copysign(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum0(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::deepcopy(const A &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erf(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erfinv(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erfinv(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fabs(int x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(int x, int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(int x, int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getDescription(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getPtr(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getPtr(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getRealTime() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getRepresentation(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(const std::vector< double > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(std::vector< double > &v) {
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
%exception  casadi::if_else_zero(const T &x, const T &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::if_else_zero(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::inner_prod(const std::vector< T > &a, const std::vector< T > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isEqual(double x, double y, int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_a(const SharedObject &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isinf(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isnan(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::iszero(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::linspace(std::vector< T > &v, const F &first, const L &last) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::makeVector(int size, int ind0=-1, const T &val0=T(), int ind1=-1, const T &val1=T(), int ind2=-1, const T &val2=T(), int ind3=-1, const T &val3=T(), int ind4=-1, const T &val4=T(), int ind5=-1, const T &val5=T(), int ind6=-1, const T &val6=T(), int ind7=-1, const T &val7=T(), int ind8=-1, const T &val8=T(), int ind9=-1, const T &val9=T(), int ind10=-1, const T &val10=T(), int ind11=-1, const T &val11=T(), int ind12=-1, const T &val12=T(), int ind13=-1, const T &val13=T(), int ind14=-1, const T &val14=T(), int ind15=-1, const T &val15=T(), int ind16=-1, const T &val16=T(), int ind17=-1, const T &val17=T(), int ind18=-1, const T &val18=T(), int ind19=-1, const T &val19=T()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName< SXElement >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::norm_1(const std::vector< T > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::norm_2(const std::vector< T > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::norm_inf(const std::vector< T > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::operation_checker(unsigned int op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::print(const std::vector< T > &v, std::ostream &stream=CASADI_COUT) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::printme(const T &x, const T &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::printme(double x, double y) {
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
%exception  casadi::range(int start, int stop, int step=1, int len=std::numeric_limits< int >::max()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::range(int stop) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::repr(const std::vector< T > &v, std::ostream &stream=CASADI_COUT) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::reverse(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::shared_cast(SharedObject &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::shared_cast(const SharedObject &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sign(double x) {
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
%exception  casadi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< int > &indices, bool invert_indices=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sq(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::toVector(const T &v0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::toVector(const T &v0, const T &v1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::toVector(const T &v0, const T &v1, const T &v2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::twice(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_slice(const std::vector< T > &v, const std::vector< int > &i) {
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
%exception casadi::CLEInputIOSchemeVector< M >::CLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CLEOutputIOSchemeVector< M >::CLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CleStructIOSchemeVector< M >::CleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ControlSimulatorInputIOSchemeVector< M >::ControlSimulatorInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ControlledDAEInputIOSchemeVector< M >::ControlledDAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DAEInputIOSchemeVector< M >::DAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DAEOutputIOSchemeVector< M >::DAEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DLEInputIOSchemeVector< M >::DLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DLEOutputIOSchemeVector< M >::DLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DPLEInputIOSchemeVector< M >::DPLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DPLEOutputIOSchemeVector< M >::DPLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DleStructIOSchemeVector< M >::DleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DpleVecStructIOSchemeVector< M >::DpleVecStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GradFInputIOSchemeVector< M >::GradFInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GradFOutputIOSchemeVector< M >::GradFOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::HNLPInputIOSchemeVector< M >::HNLPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::HessLagInputIOSchemeVector< M >::HessLagInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::HessLagOutputIOSchemeVector< M >::HessLagOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IntegratorInputIOSchemeVector< M >::IntegratorInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IntegratorOutputIOSchemeVector< M >::IntegratorOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::JacGInputIOSchemeVector< M >::JacGInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::JacGOutputIOSchemeVector< M >::JacGOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LPStructIOSchemeVector< M >::LPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DLEInputIOSchemeVector< M >::LR_DLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DLEOutputIOSchemeVector< M >::LR_DLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DPLEInputIOSchemeVector< M >::LR_DPLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LR_DPLEOutputIOSchemeVector< M >::LR_DPLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinearSolver::LinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolInputIOSchemeVector< M >::LinsolInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolOutputIOSchemeVector< M >::LinsolOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LpSolverInputIOSchemeVector< M >::LpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LpSolverOutputIOSchemeVector< M >::LpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDleStructIOSchemeVector< M >::LrDleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDpleVecStructIOSchemeVector< M >::LrDpleVecStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp, double val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp, int val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const MX &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const std::vector< MX > &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MXFunction::MXFunction(const std::vector< MX > &input, const MX &output) {
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
%exception casadi::NLPInputIOSchemeVector< M >::NLPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NLPOutputIOSchemeVector< M >::NLPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlpSolverInputIOSchemeVector< M >::NlpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlpSolverOutputIOSchemeVector< M >::NlpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QCQPStructIOSchemeVector< M >::QCQPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QPStructIOSchemeVector< M >::QPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpSolverInputIOSchemeVector< M >::QcqpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpSolverOutputIOSchemeVector< M >::QcqpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpSolverInputIOSchemeVector< M >::QpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpSolverOutputIOSchemeVector< M >::QpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RDAEInputIOSchemeVector< M >::RDAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RDAEOutputIOSchemeVector< M >::RDAEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPInputIOSchemeVector< M >::SDPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPOutputIOSchemeVector< M >::SDPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPStructIOSchemeVector< M >::SDPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPInputIOSchemeVector< M >::SDQPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPOutputIOSchemeVector< M >::SDQPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPStructIOSchemeVector< M >::SDQPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPInputIOSchemeVector< M >::SOCPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPOutputIOSchemeVector< M >::SOCPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPStructIOSchemeVector< M >::SOCPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElement::SXElement() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElement::SXElement(const SXElement &scalar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElement::SXElement(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const SX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< SX > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< std::vector< SXElement > > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const std::vector< SX > &arg, const SX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const SX &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const std::vector< std::vector< SXElement > > &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject(const SharedObject &ref) {
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
%exception casadi::Sparsity::Sparsity(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedQpSolverInputIOSchemeVector< M >::StabilizedQpSolverInputIOSchemeVector(const std::vector< M > &t) {
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