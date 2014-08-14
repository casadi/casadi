%exception  casadi::Assertion::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Assertion::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Assertion::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Assertion::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Assertion::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Assertion::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Assertion::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::getBinary(int op, const MX &y, bool scX, bool scY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::isBinaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinaryMX< ScX, ScY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinarySX::dep(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinarySX::dep(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinarySX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinarySX::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinarySX::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinarySX::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::BinarySX::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::getFunction() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::getFunctionInput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::getFunctionOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CallFunction::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::addAuxiliary(Auxiliary f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::addDependency(const Function &f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::addInclude(const std::string &new_include, bool relative_path=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::addSparsity(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::casadi_dot(int n, const std::string &x, int inc_x, const std::string &y, int inc_y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::copyVector(std::ostream &s, const std::string &arg, std::size_t n, const std::string &res, const std::string &it="i", bool only_if_exists=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::flush(std::ostream &s) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< double > &v, bool allow_adding=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< int > &v, bool allow_adding=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::getDependency(const Function &f) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CodeGenerator::getSparsity(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CollocationIntegrator::calculateInitialConditions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CollocationIntegrator::calculateInitialConditionsB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CollocationIntegrator::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CollocationIntegrator::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CollocationIntegrator::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CollocationIntegrator::setupFG() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Concat::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Concat::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Concat::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Concat::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Concat::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Concat::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getSetSparse(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::isValue(double val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Constant< Value >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantDMatrix::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::getInnerProd(const MX &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantSX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ConstantSX::isConstant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CplexInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CplexInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CplexInterface::freeCplex() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CplexInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CsparseInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CsparseInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CsparseInterface::prepare() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::freeCVodes() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::getJac() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::getJacB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::getJacGen() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::getJacGenB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::initAdj() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::integrate(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::integrateB(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::printStats(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::reset() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::resetB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::CvodesInterface::setStopTime(double tf) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseMultiplication< TrX, TrY >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseMultiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseTranspose::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseTranspose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseTranspose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseTranspose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseTranspose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseTranspose::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DenseTranspose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Determinant::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Determinant::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Determinant::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Determinant::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DsdpInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DsdpInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::DsdpInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditionsB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::getExplicit() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::getExplicitB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::integrate(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::integrateB(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::reset() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::resetB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FixedStepIntegrator::setupFG() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::callDerivative(const DMatrixVector &arg, DMatrixVector &output_res, const DMatrixVectorVector &fseed, DMatrixVectorVector &output_fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::callDerivative(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::callDerivative(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::checkNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::inputScheme() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::inputScheme() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::input_struct() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::input_struct() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::outputScheme() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::outputScheme() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::output_struct() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::output_struct() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::spCanEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::spEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Function::spInit(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::call(const DMatrixVector &arg, DMatrixVector &res, const DMatrixVectorVector &fseed, DMatrixVectorVector &fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &asens, bool always_inline, bool never_inline) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::call(const MXVector &arg, MXVector &res, const MXVectorVector &fseed, MXVectorVector &fsens, const MXVectorVector &aseed, MXVectorVector &asens, bool always_inline, bool never_inline) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::call(const SXVector &arg, SXVector &res, const SXVectorVector &fseed, SXVectorVector &fsens, const SXVectorVector &aseed, SXVectorVector &asens, bool always_inline, bool never_inline) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::callSelf(const std::vector< MX > &arg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::createCall(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::createCallDerivative(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::derivative(int nfwd, int nadj) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr, std::string compiler) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::evalMX(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::evalSX(const std::vector< SX > &arg, std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::evalSXsparse(const std::vector< SX > &arg, std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::fullJacobian() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::generateCode(std::ostream &cfile, bool generate_main) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::generateFunction(std::ostream &stream, const std::string &fname, const std::string &input_type, const std::string &output_type, const std::string &type, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::generateIO(CodeGenerator &gen) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getDerivative(int nfwd, int nadj) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getDerivativeViaJac(int nfwd, int nadj) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getFullJacobian() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getGradient(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getHessian(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getJacSparsity(int iind, int oind, bool symmetric) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchical(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getJacSparsityPlain(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getNumInputElements() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getNumInputNonzeros() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getNumOutputElements() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getNumOutputNonzeros() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getNumericJacobian(int iind, int oind, bool compact, bool symmetric) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getPartition(int iind, int oind, Sparsity &D1, Sparsity &D2, bool compact, bool symmetric) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getSanitizedName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getStat(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getStats() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::getTangent(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::gradient(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::hessian(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::inputScheme() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::inputScheme() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::input_struct() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::input_struct() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::jacSparsity(int iind, int oind, bool compact, bool symmetric) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::jacobian(int iind, int oind, bool compact, bool symmetric) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::log(const std::string &fcn, const std::string &msg) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::log(const std::string &msg) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::monitored(const std::string &mod) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::outputScheme() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::outputScheme() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::output_struct() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::output_struct() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::setDerivative(const Function &fcn, int nfwd, int nadj) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::setJacSparsity(const Sparsity &sp, int iind, int oind, bool compact) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::setJacobian(const Function &jac, int iind, int oind, bool compact) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::spCanEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::spEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::spEvaluateViaJacSparsity(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::spInit(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::symbolicInput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::symbolicInputSX() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::symbolicOutput(const std::vector< MX > &arg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::tangent(int iind, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::verbose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::FunctionInternal::wrapMXFunction() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericMatrix< MX  >::shape() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericMatrix< MatType >::shape() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::shape() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericType::is_a() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericType::toDictionary() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericType::toDictionary() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericType::toDoubleVector() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericType::toDoubleVector() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericType::toIntVector() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GenericType::toIntVector() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzeros::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzeros::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzeros::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzeros::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzeros::mapping() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice2::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosSlice::simplifyMe(MX &ex) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::GetNonzerosVector::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzcat::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzcat::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzcat::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzsplit::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzsplit::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzsplit::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Horzsplit::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, const std::string &iname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, const std::string &oname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::inputSchemeEntry(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::outputSchemeEntry(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Derived >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::getInput(T val, const std::string &iname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::getInput(T val, int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, const std::string &oname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::inputS(int i) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::inputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::inputSchemeEntry(const std::string &name) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::outputSchemeEntry(const std::string &name) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< Function  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, const std::string &iname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(const std::string &iname) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(int iind=0) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInputScheme() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumInputs() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumOutputs() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, const std::string &oname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(const std::string &oname) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(int oind=0) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutputScheme() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputSchemeEntry(const std::string &name) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputSchemeEntry(const std::string &name) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, const std::string &iname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInputScheme(const casadi::IOScheme &scheme) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumInputs(int num_in) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumOutputs(int num_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, const std::string &oname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutputScheme(const casadi::IOScheme &scheme) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOScheme::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOScheme::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOSchemeVector< M  >::print(std::ostream &stream=std::cout) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOSchemeVector< M  >::repr(std::ostream &stream=std::cout) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOSchemeVector< M  >::vector() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOSchemeVector< T >::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IOSchemeVector< T >::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::correctInitialConditions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::freeIDAS() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::getJac() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::getJacB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::getJacGen() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::getJacGenB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initAdj() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initBandedLinearSolver() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initBandedLinearSolverB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initDenseLinearSolver() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initDenseLinearSolverB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initIterativeLinearSolver() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initIterativeLinearSolverB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initTaping() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolver() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolverB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::integrate(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::integrateB(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::printStats(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::reset() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::resetB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IdasInterface::setStopTime(double tf) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFixedStepIntegrator::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFixedStepIntegrator::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicit() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicitB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFixedStepIntegrator::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFunctionInternal::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFunctionInternal::spCanEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ImplicitFunctionInternal::spEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IndexList::getAll(int len) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InfSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InfSX::isInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::InnerProd::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegerSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegerSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegerSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::getAugOffset(int nfwd, int nadj) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::getAugmented(int nfwd, int nadj, AugOffset &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::getDerivative(int nfwd, int nadj) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::integrate(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::integrateB(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::p() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::printStats(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::qf() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::resetB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::rp() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::rqf() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::rx0() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::rxf() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::rz0() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::rzf() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::setDerivativeOptions(Integrator &integrator, const AugOffset &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::setStopTime(double tf) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::spCanEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::spEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::spJacF() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::spJacG() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::x0() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::xf() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::z0() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IntegratorInternal::zf() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Inverse::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Inverse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Inverse::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Inverse::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::finalize_metadata(int n, const std::map< std::string, std::vector< std::string > > &var_string_md, const std::map< std::string, std::vector< int > > &var_integer_md, const std::map< std::string, std::vector< double > > &var_numeric_md, int m, const std::map< std::string, std::vector< std::string > > &con_string_md, const std::map< std::string, std::vector< int > > &con_integer_md, const std::map< std::string, std::vector< double > > &con_numeric_md) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::freeIpopt() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::getReducedHessian() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::get_nlp_info(int &n, int &m, int &nnz_jac_g, int &nnz_h_lag) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::get_number_of_nonlinear_variables() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::get_var_con_metadata(int n, std::map< std::string, std::vector< std::string > > &var_string_md, std::map< std::string, std::vector< int > > &var_integer_md, std::map< std::string, std::vector< double > > &var_numeric_md, int m, std::map< std::string, std::vector< std::string > > &con_string_md, std::map< std::string, std::vector< int > > &con_integer_md, std::map< std::string, std::vector< double > > &con_numeric_md) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptInterface::setQPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptUserClass::finalize_metadata(Index n, const StringMetaDataMapType &var_string_md, const IntegerMetaDataMapType &var_integer_md, const NumericMetaDataMapType &var_numeric_md, Index m, const StringMetaDataMapType &con_string_md, const IntegerMetaDataMapType &con_integer_md, const NumericMetaDataMapType &con_numeric_md) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptUserClass::get_number_of_nonlinear_variables() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::IpoptUserClass::get_var_con_metadata(Index n, StringMetaDataMapType &var_string_md, IntegerMetaDataMapType &var_integer_md, NumericMetaDataMapType &var_numeric_md, Index m, StringMetaDataMapType &con_string_md, IntegerMetaDataMapType &con_integer_md, NumericMetaDataMapType &con_numeric_md) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::bjac(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::create(const Function &f, const Function &jac, const LinearSolver &linsol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::djac(long N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::func(N_Vector u, N_Vector fval) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::kinsol_error(const std::string &module, int flag, bool fatal=true) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::lsetup(KINMem kin_mem) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::psetup(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector tmp1, N_Vector tmp2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::psolve(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v, N_Vector tmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KinsolInterface::solveNonLinear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KnitroInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KnitroInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::KnitroInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LapackLuDense::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LapackLuDense::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LapackLuDense::prepare() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LapackQrDense::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LapackQrDense::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LapackQrDense::prepare() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolver::spSolve(DMatrix &X, const DMatrix &B, bool transpose=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::colind() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::evaluateDGen(const DMatrixPtrV &input, DMatrixPtrV &output, bool tr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::evaluateMXGen(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given, bool tr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::evaluateSXGen(const SXPtrV &input, SXPtrV &output, bool tr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::getFactorization(bool transpose) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::getFactorizationSparsity(bool transpose) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::ncol() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::nnz() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::nrow() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::propagateSparsityGen(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd, bool transpose) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::row() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::solve(bool transpose) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::solve(const MX &A, const MX &B, bool transpose) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LinearSolverInternal::spSolve(DMatrix &X, const DMatrix &B, bool transpose) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LpSolverInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LpSolverInternal::solve() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LpToQp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LpToQp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::LpToQp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::T() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::at(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::at(int k) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::getNZ(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::getNZ(const Slice &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::getNZ(const std::vector< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::getNZ(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::getTemp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setNZ(const Matrix< int > &k, const MX &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setNZ(const Slice &k, const MX &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setNZ(const std::vector< int > &k, const MX &el) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setNZ(int k, const MX &el) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const Slice &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const Sparsity &sp, int dummy) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, Slice cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, int rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, int rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setSub(const MX &m, int rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::setTemp(int t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Slice &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Slice &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Slice &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const Sparsity &sp, int dummy=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(int rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(int rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(int rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MX::sub(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXFunction::algorithm() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXFunction::generateLiftingFunctions(MXFunction &vdef_fcn, MXFunction &vinit_fcn) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::addDependency(const MX &dep) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, bool add=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, const std::vector< int > &onz, bool add=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::dep(int ind=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::dep(int ind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getAddNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getAssertion(const MX &y, const std::string &fail_message="") const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getBinary(int op, const MX &y, bool scX, bool scY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getBinarySwitch(int op, const MX &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getDeterminant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getFunction() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getFunction() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getFunctionInput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getFunctionOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getHorzsplit(const std::vector< int > &output_offset) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getInnerProd(const MX &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getInverse() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getMultiplication(const MX &y, const Sparsity &sp_z=Sparsity()) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getNorm1() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getNorm2() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getNormF() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getNormInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getOutput(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getSetSparse(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getSolve(const MX &r, bool tr, const LinearSolver &linear_solver) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getSubAssign(const MX &y, const Slice &i, const Slice &j) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getSubRef(const Slice &i, const Slice &j) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::getVertsplit(const std::vector< int > &output_offset) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isBinaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isMultipleOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isNonLinear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isOutputNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isUnaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isValue(double val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::mapping() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::numel() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::setDependencies(const MX &dep) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2, const MX &dep3) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::setDependencies(const std::vector< MX > &dep) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::setSparsity(const Sparsity &sparsity) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::shape() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::simplifyMe(MX &ex) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::size() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::size1() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::size2() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::sparsity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MXNode::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::append(const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::appendColumns(const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::arccos() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::arccosh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::arcsin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::arcsinh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::arctan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::arctan2(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::arctanh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::binary(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::borBV(const Matrix< DataType > &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::ceil() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::className() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::clear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::colind() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::colind(int col) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::cos() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::cosh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::data() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::data() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::densify(const DataType &val=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::enlarge(int nrow, int ncol, const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::erase(const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::erf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::erfinv() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::exp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::fabs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::floor() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::fmax(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::fmin(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::fmod(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::get(DataType &val, SparsityType sp=SPARSE) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &val, SparsityType sp=SPARSE) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::get(std::vector< DataType > &val, SparsityType sp=SPARSE) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::getNZ(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::getNZ(const std::vector< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::hasNonStructuralZeros() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::if_else_zero(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::inf(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::inf(const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::inf(int nrow=1, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isConstant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isEqual(const Matrix< DataType > &ex2) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isRegular() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isSymbolic() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isSymbolicSparse() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::log() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::log10() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::logic_and(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::logic_not() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::logic_or(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::matrix_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::matrix_scalar(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::mul(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::mul_full(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::nan(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::nan(const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::nan(int nrow=1, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::printDense(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::printScalar(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::printSparse(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::printVector(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::printme(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::remove(const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::repmat(const DataType &x, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, int nrow, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::resize(int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::row() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::row(int el) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sanityCheck(bool complete=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::scalar_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::set(DataType val, SparsityType sp=SPARSE) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &val, SparsityType sp=SPARSE) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::set(const std::vector< DataType > &val, SparsityType sp=SPARSE) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setAll(const DataType &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setBV(const Matrix< DataType > &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setNZ(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setNZ(const std::vector< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setNZ(int k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setSparse(const Sparsity &sp, bool intersect=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Sparsity &sp, int dummy) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, int rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setZero() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::setZeroBV() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sign() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sinh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sparsify(double tol=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sparsityRef() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sqrt() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sub(const Sparsity &sp, int dummy=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::sub(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::tan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::tanh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::toScalar() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::trans() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< DataType >::unary(int op, const Matrix< DataType > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::T() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::addSub(const Matrix< DataType > &m, RR rr, CC cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::at(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::at(int k) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::back() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::back() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::begin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::begin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::end() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::end() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::front() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::front() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::getBV(Matrix< DataType > &val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::getSub(Matrix< DataType > &m, RR rr, CC cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr, const IndexList &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const IndexList &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const IndexList &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const IndexList &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Slice &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Slice &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Sparsity &sp, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_one_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_one_based(int rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_one_based(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(int rr, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(int rr, int cc, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_zero_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_zero_based(int rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_zero_based(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(int rr, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(int rr, int cc, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed(const IndexList &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed(const Slice &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_assignment(const IndexList &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_assignment(const Slice &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based_assignment(int k, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based_assignment(int k, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::ptr() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::ptr() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::rbegin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::rbegin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::rend() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::rend() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const int rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, int rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, int rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const Matrix< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const Matrix< int > &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const std::vector< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(const std::vector< int > &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(int rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(int rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Matrix< T >::sub(int rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MinusInfSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MinusInfSX::isMinusInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MinusOneSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MinusOneSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MinusOneSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MinusOneSX::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MultipleOutput::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MultipleOutput::getOutput(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MultipleOutput::isMultipleOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::MultipleOutput::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Multiplication< TrX, TrY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NanSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NanSX::isNan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Newton::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Newton::create(const Function &f, const Function &jac, const LinearSolver &linsol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Newton::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Newton::solveNonLinear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::checkInitialBounds() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::getGradF() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::getGradLag() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::getHessLag() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::getJacF() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::getJacG() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::getReducedHessian() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::getSpHessLag() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::gradF() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::gradLag() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::hessLag() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::jacF() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::jacG() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::reportConstraints(std::ostream &stream=std::cout) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::setOptionsFromFile(const std::string &file) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::setQPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NlpSolverInternal::spHessLag() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NonZeroIterator< DataType >::begin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NonZeroIterator< DataType >::end() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Norm1::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Norm1::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Norm1::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Norm2::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Norm2::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Norm2::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormF::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormInf::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormInf::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::NormInf::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OldCollocationIntegrator::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OldCollocationIntegrator::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OldCollocationIntegrator::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OldCollocationIntegrator::integrate(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OldCollocationIntegrator::integrateB(double t_out) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OldCollocationIntegrator::reset() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OldCollocationIntegrator::resetB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OneSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OneSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OneSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OneSX::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OoqpInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OoqpInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OoqpInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionality::checkNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionality::getOptionAllowedIndex(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionality::getOptionEnumValue(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionality::setOptionByAllowedIndex(const std::string &name, int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionality::setOptionByEnumValue(const std::string &name, int v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::copyOptions(const OptionsFunctionality &obj, bool skipUnknown=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::dictionary() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getBestMatches(const std::string &name, std::vector< std::string > &suggestions, int amount=5) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOption(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowed(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowedIndex(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDefault(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDescription(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionEnumValue(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionNames() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionType(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::getOptionTypeName(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::hasOption(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::hasSetOption(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::printOption(const std::string &name, std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::printOptions(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const Dictionary &dict, bool skipUnknown=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const std::string &str, const GenericType &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByAllowedIndex(const std::string &name, int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByEnumValue(const std::string &name, int v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::getFunctionInput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::getFunctionOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::isNonLinear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::isOutputNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::OutputNode::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Polynomial::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Polynomial::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::PrintableObject::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::PrintableObject::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ProfilingType() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ProfilingType< ProfilingData_ENTRY >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ProfilingType< ProfilingData_EXIT >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ProfilingType< ProfilingData_IO >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ProfilingType< ProfilingData_NAME >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ProfilingType< ProfilingData_SOURCE >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ProfilingType< ProfilingData_TIMELINE >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QcqpSolverInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QcqpSolverInternal::setQPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QcqpSolverInternal::solve() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QcqpToSocp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QcqpToSocp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QcqpToSocp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpSolverInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpSolverInternal::generateNativeCode(std::ostream &file) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpSolverInternal::setLPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpSolverInternal::solve() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToImplicit::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToImplicit::create(const Function &f, const Function &jac, const LinearSolver &linsol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToImplicit::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToImplicit::solveNonLinear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToNlp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToNlp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToNlp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToQcqp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToQcqp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpToQcqp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpoasesInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpoasesInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::QpoasesInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::RealtypeSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::RealtypeSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::RealtypeSX::isAlmostZero(double tol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Reshape::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::RkIntegrator::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::RkIntegrator::create(const Function &f, const Function &g) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::RkIntegrator::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::RkIntegrator::setupFG() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::assignIfDuplicate(const SXElement &scalar, int depth=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::assignNoDelete(const SXElement &scalar) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::get() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::getTemp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::mark() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::marked() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::setTemp(int t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXElement::toString() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXFunction::algorithm() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::dep(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::dep(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isAlmostZero(double tol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isConstant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isMinusInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isNan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isSymbolic() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::mark() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::marked() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SXNode::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::dualInfeasibility() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::eval_exp() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::eval_mat() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::eval_res() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::eval_vec() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::getQpSolver() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::line_search(int &ls_iter, bool &ls_success) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::primalInfeasibility() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double reg, int ls_trials, bool ls_success) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::regularize() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Scpgen::solve_qp() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdpSolverInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdpSolverInternal::setSOCPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdpSolverInternal::solve() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdqpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdqpSolverInternal::setSOCQPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdqpSolverInternal::solve() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdqpToSdp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdqpToSdp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SdqpToSdp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzeros< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzeros< Add >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzeros< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzeros< Add >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzeros< Add >::mapping() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzeros< Add >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice2< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::isAssignment() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosSlice< Add >::simplifyMe(MX &ex) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetNonzerosVector< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SetSparse::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::assertInit() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::checkNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::get() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::get() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::getCount() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::printPtr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::swap(SharedObject &other) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObject::weak() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::assertInit() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::getCount() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::isInit() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SharedObjectNode::weak() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SimpleHomotopyNlp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SimpleHomotopyNlp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SimpleHomotopyNlp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Slice::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::formatStatus(int status) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::passOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::reset() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::setOptionsFromFile(const std::string &file) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SnoptInterface::setQPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SocpSolverInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SocpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SocpSolverInternal::solve() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SocpToSdp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SocpToSdp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SocpToSdp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::getFunction() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Solve< Tr >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::at(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::at(int k) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::back() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::back() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::begin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::begin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::clear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::colind() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::colind(int col) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::data() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::data() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::end() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::end() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::front() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::front() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::getElement(int rr, int cc=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::hasNZ(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::rbegin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::rbegin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::rend() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::rend() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::resize(int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::row() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::row(int el) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::sanityCheck(bool complete=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::sparsity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::sparsityRef() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SparseStorage< DataType >::toScalar() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sparsity::T() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sparsity::colindRef() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sparsity::reCache() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sparsity::rowRef() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sparsity::shape() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Split::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Split::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Split::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Split::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Split::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Split::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Split::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SqicInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SqicInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SqicInterface::generateNativeCode(std::ostream &file) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SqicInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::eval_f(const std::vector< double > &x, double &f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::eval_g(const std::vector< double > &x, std::vector< double > &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::getQpSolver() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::getRegularization(const Matrix< double > &H) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, int ls_trials, bool ls_success) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::regularize(Matrix< double > &H, double reg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::reset_h() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Sqpmethod::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedQpSolverInternal::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedQpSolverInternal::setLPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedQpSolverInternal::solve() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedQpToQp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedQpToQp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedQpToQp::generateNativeCode(std::ostream &file) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedQpToQp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqicInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqicInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqicInterface::generateNativeCode(std::ostream &file) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqicInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::eval_f(const std::vector< double > &x, double &f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::eval_g(const std::vector< double > &x, std::vector< double > &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::getRegularization(const Matrix< double > &H) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::getStabilizedQpSolver() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::mat_vec(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::mat_vectran(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::meritfg() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::norm1matrix(const DMatrix &A) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, double TRdelta, int ls_trials, bool ls_success, char info) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::regularize(Matrix< double > &H, double reg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::reset_h() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::StabilizedSqp::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt, double muR, const std::vector< double > &mu, const std::vector< double > &muE) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubAssign::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SubRef::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SundialsInterface::getJac() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SundialsInterface::getJacB() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SundialsInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SundialsInterface::reset() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SundialsInterface::setStopTime(double tf) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicNLP::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicNLP::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicOCP::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicOCP::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicQr::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicQr::evaluateSXGen(const SXPtrV &input, SXPtrV &output, bool tr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicQr::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicQr::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicQr::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicQr::prepare() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicSX::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicSX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicSX::isSymbolic() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::TinyXmlInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::TinyXmlInterface::parse(const std::string &filename) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Transpose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::getBinary(int op, const MX &y, bool scX, bool scY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::isUnaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnaryMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnarySX::dep(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnarySX::dep(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnarySX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnarySX::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnarySX::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnarySX::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::UnarySX::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertcat::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertcat::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertcat::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertsplit::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertsplit::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertsplit::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::Vertsplit::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WeakRef::alive() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WeakRef::shared() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::evaluate() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::formatStatus(int status) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::passOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::reset() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::setOptionsFromFile(const std::string &file) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::WorhpInterface::setQPOptions() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlFile::parse(const std::string &filename) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlFileInternal::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::checkName(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::dump(std::ostream &stream, int indent=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::getAttribute(const std::string &attribute_name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::getText() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::getText(T &val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::hasAttribute(const std::string &attribute_name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::hasChild(const std::string &childname) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::readAttribute(const std::string &attribute_name, T &val, bool assert_existance=true) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::setAttribute(const std::string &attribute_name, const std::string &attribute) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::setName(const std::string &name) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XmlNode::size() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getSetSparse(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroByZero::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroSX::isAlmostZero(double tol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ZeroSX::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::abs(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::acos(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::acosh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::acosh(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::asin(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::asinh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::asinh(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::atan(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::atan2(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::atan2(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::atan2(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::atanh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::atanh(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::blkdiag(const MX &A, const MX &B) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::blkdiag(const Sparsity &a, const Sparsity &b) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::blockcat(const MX &A, const MX &B, const MX &C, const MX &D) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::blockcat(const Matrix< DataType > &A, const Matrix< DataType > &B, const Matrix< DataType > &C, const Matrix< DataType > &D) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::blockmatrix(SX array[n]) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::blockmatrix(SX array[n][m]) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::casadi_load_linearsolver_csparsecholesky() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ceil(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::collocationInterpolators(const std::vector< double > &tau_root, std::vector< std::vector< double > > &C, std::vector< double > &D) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::collocationPointsL(int order, const std::string &scheme) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::constpow(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::copysign(const T &x, const T &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::copysign(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::cos(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::cosh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::deepcopy(const A &a) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::erf(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::erf(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::erfinv(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::erfinv(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::exp(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fabs(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::floor(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmax(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmax(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmax(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmax(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmax(int x, int y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmin(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmin(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmin(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmin(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmin(int x, int y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::fmod(const T &x, const T &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::getDescription(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::getPtr(Matrix< DataType > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::getPtr(const Matrix< DataType > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::getPtr(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::getPtr(std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::getRealTime() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::getRepresentation(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::get_bvec_t(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::get_bvec_t(const std::vector< double > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::get_bvec_t(std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::get_bvec_t(std::vector< double > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::hash_combine(std::size_t &seed, T v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::hash_combine(std::size_t &seed, const std::vector< int > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::hash_sparsity(int nrow, int ncol, const std::vector< int > &colind, const std::vector< int > &row) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::hash_value(T v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::horzcat(const MX &a, const MX &b) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::horzcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::horzcat(const Sparsity &a, const Sparsity &b) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::if_else(const SXElement &cond, const T &if_true, const T &if_false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::if_else_zero(const T &x, const T &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::if_else_zero(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::inner_prod(const std::vector< T > &a, const std::vector< T > &b) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::is_a(const SharedObject &A) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::isinf(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::isnan(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::linspace(std::vector< T > &v, const F &first, const L &last) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::log(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::log10(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::makeVector(int size, int ind0=-1, const T &val0=T(), int ind1=-1, const T &val1=T(), int ind2=-1, const T &val2=T(), int ind3=-1, const T &val3=T(), int ind4=-1, const T &val4=T(), int ind5=-1, const T &val5=T(), int ind6=-1, const T &val6=T(), int ind7=-1, const T &val7=T(), int ind8=-1, const T &val8=T(), int ind9=-1, const T &val9=T(), int ind10=-1, const T &val10=T(), int ind11=-1, const T &val11=T(), int ind12=-1, const T &val12=T(), int ind13=-1, const T &val13=T(), int ind14=-1, const T &val14=T(), int ind15=-1, const T &val15=T(), int ind16=-1, const T &val16=T(), int ind17=-1, const T &val17=T(), int ind18=-1, const T &val18=T(), int ind19=-1, const T &val19=T()) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::matrixName< SXElement >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::norm_1(const std::vector< T > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::norm_2(const std::vector< T > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::norm_inf(const std::vector< T > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::operation_checker(unsigned int op) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::pow(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::pow(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::pow(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::print(const std::vector< T > &v, std::ostream &stream=std::cout) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::printme(const T &x, const T &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::printme(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::profileWrite(std::ofstream &f, const T &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::profileWriteBare(std::ofstream &f, const T &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ptrVec(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ptrVec(const std::vector< std::vector< T > > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ptrVec(std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::ptrVec(std::vector< std::vector< T > > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::qr(const Matrix< DataType > &A, Matrix< DataType > &Q, Matrix< DataType > &R) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::range(int start, int stop, int step, int len) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::range(int stop) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::repr(const std::vector< T > &v, std::ostream &stream=std::cout) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::reshape(const MX &x, int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::shared_cast(SharedObject &A) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::shared_cast(const SharedObject &A) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::sign(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::sign(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::simplify(SXElement &ex) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::sin(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::sinh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &dwork, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< int > &indices, bool invert_indices=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::sq(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::sqrt(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, bool reverse=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, std::vector< MX > &ex, bool reverse=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::tan(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::tanh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::toVector(const T &v0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::toVector(const T &v0, const T &v1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::toVector(const T &v0, const T &v1, const T &v2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::twice(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::vertcat(const MX &a, const MX &b) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::vertcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::vertcat(const Sparsity &a, const Sparsity &b) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  snlog2_() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  snlog_() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  sqicDestroy() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  sqlog_() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Assertion::Assertion(const MX &x, const MX &y, const std::string &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::BinaryMX< ScX, ScY >::BinaryMX(Operation op, const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::CallFunction::CallFunction(const Function &fcn, std::vector< MX > arg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::CollocationIntegrator::CollocationIntegrator(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Concat::Concat(const std::vector< MX > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Constant< Value >::Constant(const Sparsity &sp, Value v=Value()) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ConstantDMatrix::ConstantDMatrix(const Matrix< double > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ConstantMX::ConstantMX(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ControlSimulatorInputIOSchemeVector< M >::ControlSimulatorInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ControlledDAEInputIOSchemeVector< M >::ControlledDAEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::CplexInterface::CplexInterface() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::CplexInterface::CplexInterface(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::CsparseInterface::CsparseInterface(const CsparseInterface &linsol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::CsparseInterface::CsparseInterface(const Sparsity &sp, int nrhs) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::CvodesInterface::CvodesInterface(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::DAEInputIOSchemeVector< M >::DAEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::DAEOutputIOSchemeVector< M >::DAEOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::DPLEInputIOSchemeVector< M >::DPLEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::DPLEOutputIOSchemeVector< M >::DPLEOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::DenseMultiplication< TrX, TrY >::DenseMultiplication(const MX &z, const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::DenseTranspose::DenseTranspose(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Determinant::Determinant(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::DsdpInterface::DsdpInterface(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::EmptySparsity::EmptySparsity() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::FixedStepIntegrator::FixedStepIntegrator(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::GetNonzeros::GetNonzeros(const Sparsity &sp, const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::GetNonzerosSlice2::GetNonzerosSlice2(const Sparsity &sp, const MX &x, const Slice &inner, const Slice &outer) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::GetNonzerosSlice::GetNonzerosSlice(const Sparsity &sp, const MX &x, const Slice &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::GetNonzerosVector::GetNonzerosVector(const Sparsity &sp, const MX &x, const std::vector< int > &nz) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::GradFInputIOSchemeVector< M >::GradFInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::GradFOutputIOSchemeVector< M >::GradFOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::HNLPInputIOSchemeVector< M >::HNLPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::HessLagInputIOSchemeVector< M >::HessLagInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::HessLagOutputIOSchemeVector< M >::HessLagOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Horzcat::Horzcat(const std::vector< MX > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Horzsplit::Horzsplit(const MX &x, const std::vector< int > &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IdasInterface::IdasInterface(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IndexList::IndexList() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IndexList::IndexList(const Slice &i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IndexList::IndexList(const std::vector< int > &i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IndexList::IndexList(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::InfSX::InfSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::InnerProd::InnerProd(const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IntegratorInputIOSchemeVector< M >::IntegratorInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IntegratorOutputIOSchemeVector< M >::IntegratorOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Inverse::Inverse(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::IpoptInterface::IpoptInterface(const Function &nlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::JacGInputIOSchemeVector< M >::JacGInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::JacGOutputIOSchemeVector< M >::JacGOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::KinsolInterface::KinsolInterface(const Function &f, const Function &jac, const LinearSolver &linsol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::KnitroInterface::KnitroInterface(const Function &nlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LPStructIOSchemeVector< T >::LPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LapackLuDense::LapackLuDense(const Sparsity &sparsity, int nrhs) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LapackQrDense::LapackQrDense(const Sparsity &sparsity, int nrhs) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LinearSolver::LinearSolver() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LinsolInputIOSchemeVector< M >::LinsolInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LinsolOutputIOSchemeVector< M >::LinsolOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LpSolverInputIOSchemeVector< M >::LpSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LpSolverOutputIOSchemeVector< M >::LpSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LpToQp::LpToQp(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const MX &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const std::vector< MX > &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::MXFunction::MXFunction(const std::vector< MX > &input, const MX &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::MXNode::MXNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const DataType &val=DataType(0)) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const std::vector< DataType > &d) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x, int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< std::vector< DataType > > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Matrix< DataType >::Matrix(double val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::MinusInfSX::MinusInfSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::MinusOneSX::MinusOneSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::MultipleOutput::MultipleOutput() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Multiplication< TrX, TrY >::Multiplication(const MX &z, const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NLPInputIOSchemeVector< M >::NLPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NLPOutputIOSchemeVector< M >::NLPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NanSX::NanSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Newton::Newton(const Function &f, const Function &jac, const LinearSolver &linsol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NlpSolverInputIOSchemeVector< M >::NlpSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NlpSolverOutputIOSchemeVector< M >::NlpSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NonZeroIterator< DataType >::NonZeroIterator(const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Norm1::Norm1(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Norm2::Norm2(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Norm::Norm(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NormF::NormF(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NormInf::NormInf(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OldCollocationIntegrator::OldCollocationIntegrator(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OneSX::OneSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OoqpInterface::OoqpInterface() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OoqpInterface::OoqpInterface(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OptionsFunctionalityNode::OptionsFunctionalityNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OutputNode::OutputNode(const MX &parent, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QCQPStructIOSchemeVector< T >::QCQPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QPStructIOSchemeVector< T >::QPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QcqpSolverInputIOSchemeVector< M >::QcqpSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QcqpSolverOutputIOSchemeVector< M >::QcqpSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QcqpToSocp::QcqpToSocp(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpSolverInputIOSchemeVector< M >::QpSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpSolverOutputIOSchemeVector< M >::QpSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpToImplicit::QpToImplicit(const Function &f, const Function &jac, const LinearSolver &linsol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpToNlp::QpToNlp() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpToNlp::QpToNlp(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpToQcqp::QpToQcqp(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpoasesInterface::QpoasesInterface() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QpoasesInterface::QpoasesInterface(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::RDAEInputIOSchemeVector< M >::RDAEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::RDAEOutputIOSchemeVector< M >::RDAEOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Reshape::Reshape(const MX &x, Sparsity sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::RkIntegrator::RkIntegrator(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::RuntimeConst< T >::RuntimeConst() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::RuntimeConst< T >::RuntimeConst(T v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SDPInputIOSchemeVector< M >::SDPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SDPOutputIOSchemeVector< M >::SDPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SDPStructIOSchemeVector< T >::SDPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SDQPInputIOSchemeVector< M >::SDQPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SDQPOutputIOSchemeVector< M >::SDQPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SDQPStructIOSchemeVector< T >::SDQPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SOCPInputIOSchemeVector< M >::SOCPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SOCPOutputIOSchemeVector< M >::SOCPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SOCPStructIOSchemeVector< T >::SOCPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXElement::SXElement(const SXElement &scalar) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const SX &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< SX > &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< std::vector< SXElement > > &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXFunction::SXFunction(const std::vector< SX > &arg, const SX &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const SX &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const std::vector< std::vector< SXElement > > &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SXNode::SXNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ScalarSparseSparsity::ScalarSparseSparsity() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ScalarSparsity::ScalarSparsity() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Scpgen::Scpgen(const Function &nlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SdqpToSdp::SdqpToSdp(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SetNonzeros< Add >::SetNonzeros(const MX &y, const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SetNonzerosSlice2< Add >::SetNonzerosSlice2(const MX &y, const MX &x, const Slice &inner, const Slice &outer) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SetNonzerosSlice< Add >::SetNonzerosSlice(const MX &y, const MX &x, const Slice &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SetNonzerosVector< Add >::SetNonzerosVector(const MX &y, const MX &x, const std::vector< int > &nz) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SetSparse::SetSparse(const MX &x, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SharedObject::SharedObject() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SharedObject::SharedObject(const SharedObject &ref) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SharedObjectNode::SharedObjectNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SharedObjectNode::SharedObjectNode(const SharedObjectNode &node) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SimpleHomotopyNlp::SimpleHomotopyNlp(const Function &hnlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SnoptInterface::SnoptInterface(const Function &nlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SocpToSdp::SocpToSdp(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Solve< Tr >::Solve(const MX &r, const MX &A, const LinearSolver &linear_solver) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< A > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const DataType &val=DataType(0)) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const std::vector< DataType > &d) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< A > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< A > &x, int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< DataType > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< DataType > &x, int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< std::vector< DataType > > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Split::Split(const MX &x, const std::vector< int > &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SqicInterface::SqicInterface() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SqicInterface::SqicInterface(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Sqpmethod::Sqpmethod(const Function &nlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::StabilizedQpSolverInputIOSchemeVector< M >::StabilizedQpSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::StabilizedQpToQp::StabilizedQpToQp(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface(const std::vector< Sparsity > &st) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::StabilizedSqp::StabilizedSqp(const Function &nlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SubAssign::SubAssign(const MX &x, const MX &y, const Slice &i, const Slice &j) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SubRef::SubRef(const MX &x, const Slice &i, const Slice &j) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SundialsInterface::SundialsInterface(const Function &f, const Function &g) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, int nrow=1, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SymbolicQr::SymbolicQr(const Sparsity &sparsity, int nrhs) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SymbolicSX::SymbolicSX(const std::string &name) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::TinyXmlInterface::TinyXmlInterface() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Transpose::Transpose(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::UnaryMX::UnaryMX(Operation op, MX x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Vertcat::Vertcat(const std::vector< MX > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Vertsplit::Vertsplit(const MX &x, const std::vector< int > &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::WeakRef::WeakRef(SharedObject shared) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::WeakRef::WeakRef(int dummy=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::WorhpInterface::WorhpInterface(const Function &nlp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::XmlNode::XmlNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ZeroSX::ZeroSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}