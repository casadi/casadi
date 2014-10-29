%exception  casadi::Adaptor< Derived, Solver >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< DleToLrDle , DleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< DpleToDle , DpleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< DpleToLrDple , DpleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< LpToQp , QpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< LrDleToDle , LrDleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< LrDpleToDple , LrDpleInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QcqpToSocp , SocpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QpToImplicit , NlpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QpToNlp , NlpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< QpToQcqp , QcqpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< SdqpToSdp , SdpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Adaptor< SocpToSdp , SdpSolverInternal  >::addOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Assertion::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::isBinaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinaryMX< ScX, ScY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::dep(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::dep(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BinarySX::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getFunctionInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getFunctionOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getNumOutputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CallFunction::sparsity(int oind) const  {
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
%exception  casadi::CollocationIntegrator::calculateInitialConditions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::calculateInitialConditionsB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CollocationIntegrator::setupFG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Concat::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getSetSparse(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isValue(double val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Constant< Value >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantDMatrix::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getInnerProd(const MX &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantSX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantSX::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::freeCplex() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CplexInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CsparseInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CsparseInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CsparseInterface::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::freeCVodes() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJacB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJacGen() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::getJacGenB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::initAdj() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesInterface::setStopTime(double tf) {
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
%exception  casadi::DenseMultiplication< TrX, TrY >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseMultiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DenseTranspose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Determinant::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagcat::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::getDiagcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Diagsplit::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::create(const LrDleStructure &st, const std::vector< int > &Hs) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::getDerivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DleToLrDle::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToDle::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToDle::create(const DleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToDle::create(const DleStructure &st, const std::vector< int > &Hs) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToDle::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToDle::getDerivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToDle::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToDle::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::create(const LrDpleStructure &st, const std::vector< std::vector< int > > &Hs) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::getDerivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DpleToLrDple::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DsdpInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DsdpInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DsdpInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditionsB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::getExplicit() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::getExplicitB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::setupFG() {
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
%exception  casadi::Function::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::spInit(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call(const DMatrixVector &arg, DMatrixVector &res, const DMatrixVectorVector &fseed, DMatrixVectorVector &fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &asens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call(const MXVector &arg, MXVector &res, const MXVectorVector &fseed, MXVectorVector &fsens, const MXVectorVector &aseed, MXVectorVector &asens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call(const SXVector &arg, SXVector &res, const SXVectorVector &fseed, SXVectorVector &fsens, const SXVectorVector &aseed, SXVectorVector &asens, bool always_inline, bool never_inline) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::callSelf(const std::vector< MX > &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::createCall(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::createCallDerivative(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::derivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr, std::string compiler) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::evalMX(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::evalSX(const std::vector< SX > &arg, std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::evalSXsparse(const std::vector< SX > &arg, std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::fullJacobian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateCode(std::ostream &cfile, bool generate_main) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateFunction(std::ostream &stream, const std::string &fname, const std::string &input_type, const std::string &output_type, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generateIO(CodeGenerator &gen) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getDerivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getDerivativeViaJac(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getFullJacobian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getGradient(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getHessian(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsity(int iind, int oind, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchical(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityPlain(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumInputElements() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumInputNonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumOutputElements() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumOutputNonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getNumericJacobian(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getPartition(int iind, int oind, Sparsity &D1, Sparsity &D2, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getSanitizedName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getStat(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getStats() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getTangent(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::gradient(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::hessian(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::inputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::input_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::input_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacSparsity(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacobian(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::log(const std::string &fcn, const std::string &msg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::log(const std::string &msg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::monitored(const std::string &mod) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputScheme() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::outputScheme() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::output_struct() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::output_struct() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::setDerivative(const Function &fcn, int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::setJacSparsity(const Sparsity &sp, int iind, int oind, bool compact) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::setJacobian(const Function &jac, int iind, int oind, bool compact) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spEvaluateViaJacSparsity(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::spInit(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicInputSX() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicOutput(const std::vector< MX > &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::tangent(int iind, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::verbose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::wrapMXFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::shape() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::shape() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::shape() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_a() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDictionary() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDictionary() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDoubleVector() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toDoubleVector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVector() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::toIntVector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzeros::mapping() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice2::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosSlice::simplifyMe(MX &ex) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GetNonzerosVector::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzcat::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Horzsplit::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::inputSchemeEntry(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::outputSchemeEntry(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< Derived >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const  {
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
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(int iind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInputScheme() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumInputs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumOutputs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(int oind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutputScheme() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputSchemeEntry(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, const std::string &iname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, int iind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInputScheme(const casadi::IOScheme &scheme) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumInputs(int num_in) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumOutputs(int num_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, const std::string &oname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, int oind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutputScheme(const casadi::IOScheme &scheme) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOScheme::print(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOScheme::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::print(std::ostream &stream=std::cout, bool trailing_newline=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::repr(std::ostream &stream=std::cout, bool trailing_newline=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IOSchemeVector< M  >::vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::correctInitialConditions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::freeIDAS() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJacB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJacGen() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::getJacGenB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initAdj() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initBandedLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initBandedLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initDenseLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initDenseLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initIterativeLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initIterativeLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initTaping() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolver() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolverB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IdasInterface::setStopTime(double tf) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicit() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicitB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFunctionInternal::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IndexList::getAll(int len) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InfSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InfSX::isInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::InnerProd::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegerSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegerSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegerSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getAugOffset(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getAugmented(int nfwd, int nadj, AugOffset &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getDerivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::p() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::qf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rqf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rx0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rxf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rz0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::rzf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::setDerivativeOptions(Integrator &integrator, const AugOffset &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::setStopTime(double tf) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spCanEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spEvaluate(bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spJacF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::spJacG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::x0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::xf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::z0() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IntegratorInternal::zf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Inverse::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::finalize_metadata(int n, const std::map< std::string, std::vector< std::string > > &var_string_md, const std::map< std::string, std::vector< int > > &var_integer_md, const std::map< std::string, std::vector< double > > &var_numeric_md, int m, const std::map< std::string, std::vector< std::string > > &con_string_md, const std::map< std::string, std::vector< int > > &con_integer_md, const std::map< std::string, std::vector< double > > &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::freeIpopt() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::getReducedHessian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::get_nlp_info(int &n, int &m, int &nnz_jac_g, int &nnz_h_lag) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::get_number_of_nonlinear_variables() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::get_var_con_metadata(int n, std::map< std::string, std::vector< std::string > > &var_string_md, std::map< std::string, std::vector< int > > &var_integer_md, std::map< std::string, std::vector< double > > &var_numeric_md, int m, std::map< std::string, std::vector< std::string > > &con_string_md, std::map< std::string, std::vector< int > > &con_integer_md, std::map< std::string, std::vector< double > > &con_numeric_md) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptInterface::setQPOptions() {
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
%exception  casadi::KinsolInterface::bjac(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::create(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::djac(long N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::func(N_Vector u, N_Vector fval) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::kinsol_error(const std::string &module, int flag, bool fatal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::lsetup(KINMem kin_mem) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::psetup(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector tmp1, N_Vector tmp2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::psolve(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v, N_Vector tmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KinsolInterface::solveNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KnitroInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KnitroInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::KnitroInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLuDense::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLuDense::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLuDense::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQrDense::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQrDense::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQrDense::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolver::spSolve(DMatrix &X, const DMatrix &B, bool transpose=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::evaluateDGen(const DMatrixPtrV &input, DMatrixPtrV &output, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::evaluateMXGen(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::evaluateSXGen(const SXPtrV &input, SXPtrV &output, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::getFactorization(bool transpose) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::getFactorizationSparsity(bool transpose) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::ncol() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::nnz() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::nrow() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::propagateSparsityGen(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd, bool transpose) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::solve(bool transpose) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::solve(const MX &A, const MX &B, bool transpose) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearSolverInternal::spSolve(DMatrix &X, const DMatrix &B, bool transpose) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpToQp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpToQp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LpToQp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::create(const DleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::getDerivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDleToDle::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::create(const DpleStructure &st) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::getDerivative(int nfwd, int nadj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LrDpleToDple::printStats(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::T() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::at(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::at(int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getNZ(const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getNZ(const Slice &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getNZ(const std::vector< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getNZ(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::getTemp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setNZ(const Matrix< int > &k, const MX &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setNZ(const Slice &k, const MX &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setNZ(const std::vector< int > &k, const MX &el) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setNZ(int k, const MX &el) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const Slice &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const Sparsity &sp, int dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, Slice cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, int rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, int rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setSub(const MX &m, int rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Slice &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Slice &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Slice &rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const Sparsity &sp, int dummy=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(int rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(int rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(int rr, const std::vector< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sub(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXFunction::algorithm() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXFunction::generateLiftingFunctions(MXFunction &vdef_fcn, MXFunction &vinit_fcn) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::addDependency(const MX &dep) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, bool add=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, const std::vector< int > &onz, bool add=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::dep(int ind=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::dep(int ind=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getAddNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getAssertion(const MX &y, const std::string &fail_message="") const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getBinarySwitch(int op, const MX &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getDeterminant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getDiagcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getDiagsplit(const std::vector< int > &offset1, const std::vector< int > &offset2) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunction() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunctionInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getFunctionOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getHorzsplit(const std::vector< int > &output_offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getInnerProd(const MX &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getInverse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getMultiplication(const MX &y, const Sparsity &sp_z=Sparsity()) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNorm1() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNorm2() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNormF() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNormInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getNumOutputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getOutput(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSetSparse(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSolve(const MX &r, bool tr, const LinearSolver &linear_solver) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSubAssign(const MX &y, const Slice &i, const Slice &j) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getSubRef(const Slice &i, const Slice &j) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::getVertsplit(const std::vector< int > &output_offset) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isBinaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isMultipleOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isOutputNode() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isUnaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isValue(double val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::mapping() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::numel() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const MX &dep) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2, const MX &dep3) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setDependencies(const std::vector< MX > &dep) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::setSparsity(const Sparsity &sparsity) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::shape() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::simplifyMe(MX &ex) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::size() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::size1() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::size2() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MXNode::sparsity(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::append(const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::appendColumns(const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::arccos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::arccosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::arcsin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::arcsinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::arctan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::arctan2(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::arctanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::binary(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::borBV(const Matrix< DataType > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::ceil() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::colind(int col) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::cos() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::cosh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::data() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::data() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::densify(const DataType &val=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::enlarge(int nrow, int ncol, const std::vector< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::erase(const std::vector< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::erf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::erfinv() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::exp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::fabs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::floor() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::fmax(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::fmin(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::fmod(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(DataType &val, SparsityType sp=SPARSE) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &val, SparsityType sp=SPARSE) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::get(std::vector< DataType > &val, SparsityType sp=SPARSE) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNZ(const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getNZ(const std::vector< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::hasNonStructuralZeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::if_else_zero(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::inf(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::inf(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::inf(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isEqual(const Matrix< DataType > &ex2) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isIdentity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isRegular() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isSymbolicSparse() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::log() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::log10() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::logic_and(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::logic_not() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::logic_or(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::matrix_scalar(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::mul(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::mul_full(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nan(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nan(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::nan(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::print(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printDense(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printScalar(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printSparse(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printVector(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::printme(const Matrix< DataType > &y) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::remove(const std::vector< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::repmat(const DataType &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, int nrow, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::repr(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::resize(int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::row(int el) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sanityCheck(bool complete=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::scalar_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(DataType val, SparsityType sp=SPARSE) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &val, SparsityType sp=SPARSE) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::set(const std::vector< DataType > &val, SparsityType sp=SPARSE) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setAll(const DataType &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setBV(const Matrix< DataType > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setNZ(const Matrix< int > &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setNZ(const std::vector< int > &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setNZ(int k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSparse(const Sparsity &sp, bool intersect=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Sparsity &sp, int dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, int rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setZero() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::setZeroBV() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sign() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sinh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sparsify(double tol=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sparsityRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sqrt() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sub(const Sparsity &sp, int dummy=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::sub(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::tan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::tanh() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::toScalar() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::trans() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, const std::pair< int, int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< DataType >::unary(int op, const Matrix< DataType > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::T() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::addSub(const Matrix< DataType > &m, RR rr, CC cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::at(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::at(int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::back() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::back() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::begin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::begin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::end() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::end() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::front() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::front() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::getBV(Matrix< DataType > &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::getSub(Matrix< DataType > &m, RR rr, CC cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr, const IndexList &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const IndexList &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const IndexList &cc, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const IndexList &cc, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Slice &cc, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Slice &cc, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Sparsity &sp, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_one_based(const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_one_based(int rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_one_based(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(int rr, const DataType &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(int rr, int cc, const DataType &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_zero_based(const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_zero_based(int rr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_zero_based(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(int rr, const DataType &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(int rr, int cc, const DataType &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed(const IndexList &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed(const Slice &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_assignment(const IndexList &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_assignment(const Slice &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based(const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based_assignment(int k, const DataType &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based(const Matrix< int > &k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based_assignment(int k, const DataType &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::ptr() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::ptr() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rbegin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rbegin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rend() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::rend() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const int rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, int rr, const Matrix< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, int rr, const std::vector< int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const Matrix< int > &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const Matrix< int > &rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const std::vector< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const std::vector< int > &rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(const std::vector< int > &rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(int rr, const Matrix< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(int rr, const Slice &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::sub(int rr, const std::vector< int > &cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusInfSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusInfSX::isMinusInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MinusOneSX::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::getNumOutputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::getOutput(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::isMultipleOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MultipleOutput::sparsity(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Multiplication< TrX, TrY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NanSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NanSX::isNan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::create(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::solveNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::checkInitialBounds() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getGradF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getGradLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getHessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getJacF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getJacG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getReducedHessian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::getSpHessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::gradF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::gradLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::hessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::jacF() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::jacG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::reportConstraints(std::ostream &stream=std::cout) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpSolverInternal::spHessLag() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NonZeroIterator< DataType >::begin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NonZeroIterator< DataType >::end() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm1::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm1::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm1::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm2::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm2::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Norm2::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormF::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormInf::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormInf::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NormInf::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::integrate(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::integrateB(double t_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OldCollocationIntegrator::resetB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OneSX::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OoqpInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OoqpInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OoqpInterface::init() {
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
%exception  casadi::OptionsFunctionalityNode::addOption(const std::string &str, const opt_type &type, const GenericType &def_val, const std::string &desc, const std::string &allowed_vals, bool inherit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::addOption(const std::string &str, const opt_type &type, const GenericType &def_val=GenericType(), const std::string &desc="n/a", const std::vector< GenericType > &allowed_vals=std::vector< GenericType >(), bool inherit=false, std::vector< int > enum_values=std::vector< int >(), std::vector< std::string > enum_descr=std::vector< std::string >()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::copyOptions(const OptionsFunctionality &obj, bool skipUnknown=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::dictionary() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getBestMatches(const std::string &name, std::vector< std::string > &suggestions, int amount=5) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOption(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowed(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowedIndex(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDefault(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDescription(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionEnumValue(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionNames() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionType(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::getOptionTypeName(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::hasOption(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::hasSetOption(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::printOption(const std::string &name, std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::printOptions(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const Dictionary &dict, bool skipUnknown=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const std::string &str, const GenericType &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByAllowedIndex(const std::string &name, int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByEnumValue(const std::string &name, int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getFunctionInput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getFunctionOutput() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getHorzcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::isNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::isOutputNode() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OutputNode::printPart(std::ostream &stream, int part) const  {
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
%exception  casadi::QcqpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpSolverInternal::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpToSocp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpToSocp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QcqpToSocp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::setLPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::create(const Function &f) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToImplicit::solveNonLinear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToQcqp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToQcqp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToQcqp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpoasesInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpoasesInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpoasesInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RealtypeSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RealtypeSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RealtypeSX::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Reshape::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::create(const Function &f, const Function &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RkIntegrator::setupFG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::assignIfDuplicate(const SXElement &scalar, int depth=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::assignNoDelete(const SXElement &scalar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::get() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::getTemp() const  {
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
%exception  casadi::SXElement::print(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::repr(std::ostream &stream=std::cout, bool trailing_newline=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElement::setTemp(int t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXFunction::algorithm() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::dep(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::dep(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isConstant() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isMinusInf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isMinusOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isNan() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isOne() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::mark() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::marked() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXNode::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::dualInfeasibility() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_exp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_mat() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_res() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::eval_vec() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::getQpSolver() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::line_search(int &ls_iter, bool &ls_success) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::primalInfeasibility() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double reg, int ls_trials, bool ls_success) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::regularize() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::solve_qp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::setSOCPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpSolverInternal::setSOCQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpToSdp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpToSdp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SdqpToSdp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::mapping() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzeros< Add >::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice2< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::isAssignment() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosSlice< Add >::simplifyMe(MX &ex) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::getAll() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetNonzerosVector< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SetSparse::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
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
%exception  casadi::SharedObject::printPtr(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::swap(SharedObject &other) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::weak() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::assertInit() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::getCount() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::isInit() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::repr(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectNode::weak() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SimpleHomotopyNlp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SimpleHomotopyNlp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SimpleHomotopyNlp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::formatStatus(int status) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SnoptInterface::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpToSdp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpToSdp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SocpToSdp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::getFunction() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Solve< Tr >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::at(int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::at(int k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::back() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::back() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::begin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::begin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::colind(int col) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::data() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::data() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::end() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::end() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::front() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::front() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::getElement(int rr, int cc=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::hasNZ(int rr, int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::rbegin() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::rbegin() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::rend() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::rend() {
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
%exception  casadi::SparseStorage< DataType >::row() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::row(int el) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sanityCheck(bool complete=false) const  {
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
%exception  casadi::Sparsity::T() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::colindRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::reCache() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::rowRef() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::shape() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::getNumOutputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Split::sparsity(int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SqicInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_f(const std::vector< double > &x, double &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_g(const std::vector< double > &x, std::vector< double > &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::getQpSolver() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::getRegularization(const Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, int ls_trials, bool ls_success) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::regularize(Matrix< double > &H, double reg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::reset_h() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpSolverInternal::checkInputs() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpSolverInternal::setLPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpSolverInternal::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedQpToQp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::generateNativeCode(std::ostream &file) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqicInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_f(const std::vector< double > &x, double &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_g(const std::vector< double > &x, std::vector< double > &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::getRegularization(const Matrix< double > &H) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::getStabilizedQpSolver() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::mat_vec(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::mat_vectran(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::meritfg() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::norm1matrix(const DMatrix &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, double TRdelta, int ls_trials, bool ls_success, char info) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::regularize(Matrix< double > &H, double reg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::reset_h() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StabilizedSqp::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt, double muR, const std::vector< double > &mu, const std::vector< double > &muE) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubAssign::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SubRef::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getBandwidth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getBandwidthB() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getJac() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::getJacB() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsInterface::setStopTime(double tf) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::evaluateSXGen(const SXPtrV &input, SXPtrV &output, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicSX::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicSX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicSX::isSymbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::TinyXmlInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::TinyXmlInterface::parse(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::nTmp(size_t &ni, size_t &nr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Transpose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::isUnaryOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::numInplace() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnaryMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::dep(int i) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::dep(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::hasDep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::isSmooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::ndep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::UnarySX::print(std::ostream &stream, long &remaining_calls) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertcat::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::getOp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::getVertcat(const std::vector< MX > &x) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Vertsplit::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::alive() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakRef::shared() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::formatStatus(int status) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::passOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WorhpInterface::setQPOptions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< Derived >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< Derived >::evaluate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< DleToLrDle  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< DpleToLrDple  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< LrDleToDle  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Wrapper< LrDpleToDple  >::checkDimensions() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlFile::parse(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlFileInternal::print(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::checkName(const std::string &str) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::dump(std::ostream &stream, int indent=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getAttribute(const std::string &attribute_name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getText() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::getText(T &val) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::hasAttribute(const std::string &attribute_name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::hasChild(const std::string &childname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::readAttribute(const std::string &attribute_name, T &val, bool assert_existance=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::setAttribute(const std::string &attribute_name, const std::string &attribute) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::setName(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlNode::size() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::clone() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getMatrixValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getReshape(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getSetSparse(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getTranspose() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getUnary(int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroByZero::printPart(std::ostream &stream, int part) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::getIntValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::getValue() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::isAlmostZero(double tol) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::isInteger() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ZeroSX::isZero() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::abs(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::acos(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::acosh(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::acosh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::asin(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::asinh(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::asinh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atan(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atan2(const T &x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atan2(const T &x, double n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atan2(double x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atanh(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::atanh(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::blkdiag(const MX &A, const MX &B) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::blkdiag(const Sparsity &a, const Sparsity &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::blockcat(const MX &A, const MX &B, const MX &C, const MX &D) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::blockcat(const Matrix< DataType > &A, const Matrix< DataType > &B, const Matrix< DataType > &C, const Matrix< DataType > &D) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::blockmatrix(SX array[n]) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::blockmatrix(SX array[n][m]) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::casadi_load_linearsolver_csparsecholesky() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ceil(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocationInterpolators(const std::vector< double > &tau_root, std::vector< std::vector< double > > &C, std::vector< double > &D) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocationPointsL(int order, const std::string &scheme) {
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
%exception  casadi::cos(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cosh(const T &x) {
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
%exception  casadi::diagsplitNative(const MX &x, const std::vector< int > &output_offset1, const std::vector< int > &output_offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::diagsplitNative(const Matrix< DataType > &x, const std::vector< int > &output_offset1, const std::vector< int > &output_offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erf(const T &x) {
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
%exception  casadi::exp(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fabs(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::floor(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(const T &x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(const T &x, double n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(double x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmax(int x, int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(const T &x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(const T &x, double n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(double x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmin(int x, int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmod(const T &x, const T &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getDescription(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getPtr(Matrix< DataType > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::getPtr(const Matrix< DataType > &v) {
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
%exception  casadi::horzcat(const MX &a, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::horzcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::horzcat(const Sparsity &a, const Sparsity &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::if_else(const SXElement &cond, const T &if_true, const T &if_false) {
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
%exception  casadi::is_a(const SharedObject &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isinf(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isnan(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::linspace(std::vector< T > &v, const F &first, const L &last) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::log(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::log10(const T &x) {
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
%exception  casadi::pow(const T &x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::pow(const T &x, double n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::pow(double x, const T &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::print(const std::vector< T > &v, std::ostream &stream=std::cout) {
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
%exception  casadi::qr(const Matrix< DataType > &A, Matrix< DataType > &Q, Matrix< DataType > &R) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::range(int start, int stop, int step=1, int len=std::numeric_limits< int >::max()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::range(int stop) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::repr(const std::vector< T > &v, std::ostream &stream=std::cout) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::reshape(const MX &x, int nrow, int ncol) {
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
%exception  casadi::sign(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sign(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::simplify(SXElement &ex) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sin(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sinh(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &dwork, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< int > &indices, bool invert_indices=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sq(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sqrt(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, bool reverse=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, std::vector< MX > &ex, bool reverse=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::tan(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::tanh(const T &x) {
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
%exception  casadi::vertcat(const MX &a, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vertcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vertcat(const Sparsity &a, const Sparsity &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  sqicDestroy() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Assertion::Assertion(const MX &x, const MX &y, const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::BinaryMX< ScX, ScY >::BinaryMX(Operation op, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CLEInputIOSchemeVector< M >::CLEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CLEOutputIOSchemeVector< M >::CLEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CallFunction::CallFunction(const Function &fcn, std::vector< MX > arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CleStructIOSchemeVector< T >::CleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CollocationIntegrator::CollocationIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Concat::Concat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Constant< Value >::Constant(const Sparsity &sp, Value v=Value()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ConstantDMatrix::ConstantDMatrix(const Matrix< double > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ConstantMX::ConstantMX(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ControlSimulatorInputIOSchemeVector< M >::ControlSimulatorInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ControlledDAEInputIOSchemeVector< M >::ControlledDAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CplexInterface::CplexInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CplexInterface::CplexInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CsparseInterface::CsparseInterface(const CsparseInterface &linsol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CsparseInterface::CsparseInterface(const Sparsity &sp, int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CvodesInterface::CvodesInterface(const Function &f, const Function &g) {
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
%exception casadi::DenseMultiplication< TrX, TrY >::DenseMultiplication(const MX &z, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DenseTranspose::DenseTranspose(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Determinant::Determinant(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Diagcat::Diagcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Diagsplit::Diagsplit(const MX &x, const std::vector< int > &offset1, const std::vector< int > &offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DleStructIOSchemeVector< T >::DleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DleToLrDle::DleToLrDle(const LrDleStructure &st, const std::vector< int > &Hs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DpleToDle::DpleToDle(const DleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DpleToDle::DpleToDle(const DleStructure &st, const std::vector< int > &Hs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DpleToLrDple::DpleToLrDple(const LrDpleStructure &st, const std::vector< std::vector< int > > &Hs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DpleVecStructIOSchemeVector< T >::DpleVecStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DsdpInterface::DsdpInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::EmptySparsity::EmptySparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FixedStepIntegrator::FixedStepIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzeros::GetNonzeros(const Sparsity &sp, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzerosSlice2::GetNonzerosSlice2(const Sparsity &sp, const MX &x, const Slice &inner, const Slice &outer) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzerosSlice::GetNonzerosSlice(const Sparsity &sp, const MX &x, const Slice &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GetNonzerosVector::GetNonzerosVector(const Sparsity &sp, const MX &x, const std::vector< int > &nz) {
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
%exception casadi::Horzcat::Horzcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Horzsplit::Horzsplit(const MX &x, const std::vector< int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IdasInterface::IdasInterface(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IndexList::IndexList() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IndexList::IndexList(const Slice &i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IndexList::IndexList(const std::vector< int > &i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IndexList::IndexList(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::InfSX::InfSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::InnerProd::InnerProd(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IntegratorInputIOSchemeVector< M >::IntegratorInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IntegratorOutputIOSchemeVector< M >::IntegratorOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Inverse::Inverse(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::IpoptInterface::IpoptInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::JacGInputIOSchemeVector< M >::JacGInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::JacGOutputIOSchemeVector< M >::JacGOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::KinsolInterface::KinsolInterface(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::KnitroInterface::KnitroInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LPStructIOSchemeVector< T >::LPStructIOSchemeVector(const std::vector< M > &t) {
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
%exception casadi::LapackLuDense::LapackLuDense(const Sparsity &sparsity, int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LapackQrDense::LapackQrDense(const Sparsity &sparsity, int nrhs) {
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
%exception casadi::LpToQp::LpToQp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDleStructIOSchemeVector< T >::LrDleStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDleToDle::LrDleToDle(const DleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDpleToDple::LrDpleToDple(const DpleStructure &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LrDpleVecStructIOSchemeVector< T >::LrDpleVecStructIOSchemeVector(const std::vector< M > &t) {
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
%exception casadi::MXNode::MXNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const DataType &val=DataType(0)) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const std::vector< DataType > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x, int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< std::vector< DataType > > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< DataType >::Matrix(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MinusInfSX::MinusInfSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MinusOneSX::MinusOneSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MultipleOutput::MultipleOutput() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Multiplication< TrX, TrY >::Multiplication(const MX &z, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NLPInputIOSchemeVector< M >::NLPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NLPOutputIOSchemeVector< M >::NLPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NanSX::NanSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Newton::Newton(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlpSolverInputIOSchemeVector< M >::NlpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlpSolverOutputIOSchemeVector< M >::NlpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NonZeroIterator< DataType >::NonZeroIterator(const Matrix< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Norm1::Norm1(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Norm2::Norm2(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Norm::Norm(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NormF::NormF(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NormInf::NormInf(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OldCollocationIntegrator::OldCollocationIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OneSX::OneSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OoqpInterface::OoqpInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OoqpInterface::OoqpInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OptionsFunctionalityNode::OptionsFunctionalityNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OutputNode::OutputNode(const MX &parent, int oind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QCQPStructIOSchemeVector< T >::QCQPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QPStructIOSchemeVector< T >::QPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpSolverInputIOSchemeVector< M >::QcqpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpSolverOutputIOSchemeVector< M >::QcqpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QcqpToSocp::QcqpToSocp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpSolverInputIOSchemeVector< M >::QpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpSolverOutputIOSchemeVector< M >::QpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToImplicit::QpToImplicit(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToNlp::QpToNlp() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToNlp::QpToNlp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToQcqp::QpToQcqp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpoasesInterface::QpoasesInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpoasesInterface::QpoasesInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RDAEInputIOSchemeVector< M >::RDAEInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RDAEOutputIOSchemeVector< M >::RDAEOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Reshape::Reshape(const MX &x, Sparsity sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RkIntegrator::RkIntegrator(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RuntimeConst< T >::RuntimeConst() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RuntimeConst< T >::RuntimeConst(T v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPInputIOSchemeVector< M >::SDPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPOutputIOSchemeVector< M >::SDPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDPStructIOSchemeVector< T >::SDPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPInputIOSchemeVector< M >::SDQPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPOutputIOSchemeVector< M >::SDQPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SDQPStructIOSchemeVector< T >::SDQPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPInputIOSchemeVector< M >::SOCPInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPOutputIOSchemeVector< M >::SOCPOutputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SOCPStructIOSchemeVector< T >::SOCPStructIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElement::SXElement(const SXElement &scalar) {
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
%exception casadi::SXNode::SXNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ScalarSparseSparsity::ScalarSparseSparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ScalarSparsity::ScalarSparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Scpgen::Scpgen(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SdqpToSdp::SdqpToSdp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzeros< Add >::SetNonzeros(const MX &y, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzerosSlice2< Add >::SetNonzerosSlice2(const MX &y, const MX &x, const Slice &inner, const Slice &outer) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzerosSlice< Add >::SetNonzerosSlice(const MX &y, const MX &x, const Slice &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetNonzerosVector< Add >::SetNonzerosVector(const MX &y, const MX &x, const std::vector< int > &nz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SetSparse::SetSparse(const MX &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject(const SharedObject &ref) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObjectNode::SharedObjectNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObjectNode::SharedObjectNode(const SharedObjectNode &node) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SimpleHomotopyNlp::SimpleHomotopyNlp(const Function &hnlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SnoptInterface::SnoptInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SocpToSdp::SocpToSdp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Solve< Tr >::Solve(const MX &r, const MX &A, const LinearSolver &linear_solver) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< A > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const DataType &val=DataType(0)) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const std::vector< DataType > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< A > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< A > &x, int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< DataType > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< DataType > &x, int nrow, int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< std::vector< DataType > > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Split::Split(const MX &x, const std::vector< int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SqicInterface::SqicInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SqicInterface::SqicInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Sqpmethod::Sqpmethod(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedQpSolverInputIOSchemeVector< M >::StabilizedQpSolverInputIOSchemeVector(const std::vector< M > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedQpToQp::StabilizedQpToQp(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface(const std::vector< Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StabilizedSqp::StabilizedSqp(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubAssign::SubAssign(const MX &x, const MX &y, const Slice &i, const Slice &j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubRef::SubRef(const MX &x, const Slice &i, const Slice &j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SundialsInterface::SundialsInterface(const Function &f, const Function &g) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, int nrow=1, int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicQr::SymbolicQr(const Sparsity &sparsity, int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicSX::SymbolicSX(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::TinyXmlInterface::TinyXmlInterface() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Transpose::Transpose(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::UnaryMX::UnaryMX(Operation op, MX x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Vertcat::Vertcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Vertsplit::Vertsplit(const MX &x, const std::vector< int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(SharedObject shared) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(int dummy=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WorhpInterface::WorhpInterface(const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::XmlNode::XmlNode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ZeroSX::ZeroSX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}