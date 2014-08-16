%exception  casadi::Assertion::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Assertion::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Assertion::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Assertion::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Assertion::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Assertion::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Assertion::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::getUnary(int op) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::isBinaryOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::numInplace() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinaryMX< ScX, ScY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinarySX::dep(int i) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinarySX::dep(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinarySX::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinarySX::hasDep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinarySX::isSmooth() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinarySX::ndep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::BinarySX::print(std::ostream &stream, long &remaining_calls) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::getFunction() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::getFunctionInput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::getFunctionOutput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::getNumOutputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::nTmp(size_t &ni, size_t &nr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CallFunction::sparsity(int oind) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::addAuxiliary(Auxiliary f) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::addDependency(const Function &f) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::addInclude(const std::string &new_include, bool relative_path=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::addSparsity(const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::casadi_dot(int n, const std::string &x, int inc_x, const std::string &y, int inc_y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::copyVector(std::ostream &s, const std::string &arg, std::size_t n, const std::string &res, const std::string &it="i", bool only_if_exists=false) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::flush(std::ostream &s) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< double > &v, bool allow_adding=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::getConstant(const std::vector< int > &v, bool allow_adding=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::getDependency(const Function &f) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CodeGenerator::getSparsity(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CollocationIntegrator::calculateInitialConditions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CollocationIntegrator::calculateInitialConditionsB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CollocationIntegrator::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CollocationIntegrator::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CollocationIntegrator::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CollocationIntegrator::setupFG() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Concat::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Concat::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Concat::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Concat::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Concat::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Concat::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getHorzcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getMatrixValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getReshape(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getSetSparse(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getTranspose() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getUnary(int op) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::getVertcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::isIdentity() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::isOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::isValue(double val) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::isZero() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Constant< Value >::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::getMatrixValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::isIdentity() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::isMinusOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::isOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::isZero() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantDMatrix::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::getInnerProd(const MX &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::getMatrixValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantSX::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ConstantSX::isConstant() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CplexInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CplexInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CplexInterface::freeCplex() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CplexInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CsparseInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CsparseInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CsparseInterface::prepare() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::freeCVodes() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::getJac() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::getJacB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::getJacGen() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::getJacGenB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::initAdj() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::integrate(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::integrateB(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::printStats(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::reset() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::resetB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::CvodesInterface::setStopTime(double tf) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseMultiplication< TrX, TrY >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseMultiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseTranspose::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseTranspose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseTranspose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseTranspose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseTranspose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseTranspose::nTmp(size_t &ni, size_t &nr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DenseTranspose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Determinant::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Determinant::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Determinant::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Determinant::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DsdpInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DsdpInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::DsdpInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::calculateInitialConditionsB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::getExplicit() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::getExplicitB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::integrate(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::integrateB(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::reset() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::resetB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FixedStepIntegrator::setupFG() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::callDerivative(const DMatrixVector &arg, DMatrixVector &output_res, const DMatrixVectorVector &fseed, DMatrixVectorVector &output_fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::callDerivative(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::callDerivative(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::checkNode() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::inputScheme() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::inputScheme() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::input_struct() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::input_struct() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::outputScheme() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::outputScheme() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::output_struct() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::output_struct() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::spCanEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::spEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Function::spInit(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::call(const DMatrixVector &arg, DMatrixVector &res, const DMatrixVectorVector &fseed, DMatrixVectorVector &fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &asens, bool always_inline, bool never_inline) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::call(const MXVector &arg, MXVector &res, const MXVectorVector &fseed, MXVectorVector &fsens, const MXVectorVector &aseed, MXVectorVector &asens, bool always_inline, bool never_inline) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::call(const SXVector &arg, SXVector &res, const SXVectorVector &fseed, SXVectorVector &fsens, const SXVectorVector &aseed, SXVectorVector &asens, bool always_inline, bool never_inline) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::callSelf(const std::vector< MX > &arg) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::createCall(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::createCallDerivative(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::derivative(int nfwd, int nadj) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::dynamicCompilation(Function f, std::string fname, std::string fdescr, std::string compiler) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::evalMX(const std::vector< MX > &arg, std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::evalSX(const std::vector< SX > &arg, std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::evalSXsparse(const std::vector< SX > &arg, std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::fullJacobian() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::generateCode(std::ostream &cfile, bool generate_main) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::generateFunction(std::ostream &stream, const std::string &fname, const std::string &input_type, const std::string &output_type, const std::string &type, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::generateIO(CodeGenerator &gen) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getDerivative(int nfwd, int nadj) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getDerivativeViaJac(int nfwd, int nadj) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getFullJacobian() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getGradient(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getHessian(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getJacSparsity(int iind, int oind, bool symmetric) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchical(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchicalSymm(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getJacSparsityPlain(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getNumInputElements() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getNumInputNonzeros() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getNumOutputElements() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getNumOutputNonzeros() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getNumericJacobian(int iind, int oind, bool compact, bool symmetric) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getPartition(int iind, int oind, Sparsity &D1, Sparsity &D2, bool compact, bool symmetric) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getSanitizedName() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getStat(const std::string &name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getStats() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::getTangent(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::gradient(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::hessian(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::inputNoCheck(int iind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::inputScheme() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::inputScheme() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::input_struct() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::input_struct() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::jacSparsity(int iind, int oind, bool compact, bool symmetric) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::jacobian(int iind, int oind, bool compact, bool symmetric) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::log(const std::string &fcn, const std::string &msg) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::log(const std::string &msg) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::monitored(const std::string &mod) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::outputNoCheck(int oind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::outputScheme() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::outputScheme() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::output_struct() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::output_struct() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::print(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::repr(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::setDerivative(const Function &fcn, int nfwd, int nadj) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::setJacSparsity(const Sparsity &sp, int iind, int oind, bool compact) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::setJacobian(const Function &jac, int iind, int oind, bool compact) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::spCanEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::spEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::spEvaluateViaJacSparsity(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::spInit(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::symbolicInput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::symbolicInputSX() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::symbolicOutput(const std::vector< MX > &arg) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::tangent(int iind, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::verbose() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::FunctionInternal::wrapMXFunction() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericMatrix< MX  >::shape() const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericMatrix< MatType >::shape() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::shape() const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericType::is_a() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericType::toDictionary() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericType::toDictionary() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericType::toDoubleVector() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericType::toDoubleVector() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericType::toIntVector() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GenericType::toIntVector() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzeros::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzeros::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzeros::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzeros::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzeros::mapping() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice2::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::isIdentity() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosSlice::simplifyMe(MX &ex) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::GetNonzerosVector::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzcat::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzcat::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzcat::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzsplit::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzsplit::getHorzcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzsplit::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Horzsplit::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, const std::string &iname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::getInput(T val, int iind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, const std::string &oname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::getOutput(T val, int oind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::inputS(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::inputSchemeEntry(const std::string &name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::outputS(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::outputSchemeEntry(const std::string &name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Derived >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::getInput(T val, const std::string &iname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::getInput(T val, int iind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, const std::string &oname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::getOutput(T val, int oind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::inputS(int i) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::inputS(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::inputSchemeEntry(const std::string &name) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::outputS(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::outputSchemeEntry(const std::string &name) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< Function  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, const std::string &iname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(T val, int iind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(const std::string &iname) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInput(int iind=0) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getInputScheme() const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumInputs() const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getNumOutputs() const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, const std::string &oname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(T val, int oind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(const std::string &oname) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutput(int oind=0) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::getOutputScheme() const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(const std::string &iname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::input(int iind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputS(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::inputSchemeEntry(const std::string &name) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(const std::string &oname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::output(int oind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputS(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::outputSchemeEntry(const std::string &name) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, const std::string &iname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInput(T val, int iind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setInputScheme(const casadi::IOScheme &scheme) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumInputs(int num_in) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setNumOutputs(int num_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, const std::string &oname) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutput(T val, int oind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOInterface< FunctionInternal  >::setOutputScheme(const casadi::IOScheme &scheme) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOScheme::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOScheme::repr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOSchemeVector< M  >::print(std::ostream &stream=std::cout) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOSchemeVector< M  >::repr(std::ostream &stream=std::cout) const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOSchemeVector< M  >::vector() const {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOSchemeVector< T >::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IOSchemeVector< T >::repr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::correctInitialConditions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::freeIDAS() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::getJac() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::getJacB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::getJacGen() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::getJacGenB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initAdj() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initBandedLinearSolver() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initBandedLinearSolverB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initDenseLinearSolver() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initDenseLinearSolverB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initIterativeLinearSolver() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initIterativeLinearSolverB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initTaping() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolver() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::initUserDefinedLinearSolverB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::integrate(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::integrateB(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::printStats(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::reset() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::resetB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IdasInterface::setStopTime(double tf) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFixedStepIntegrator::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFixedStepIntegrator::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicit() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicitB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFixedStepIntegrator::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFunctionInternal::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFunctionInternal::spCanEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ImplicitFunctionInternal::spEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IndexList::getAll(int len) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InfSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InfSX::isInf() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::InnerProd::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegerSX::getIntValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegerSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegerSX::isInteger() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::getAugOffset(int nfwd, int nadj) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::getAugmented(int nfwd, int nadj, AugOffset &offset) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::getDerivative(int nfwd, int nadj) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::integrate(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::integrateB(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::p() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::printStats(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::qf() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::resetB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::rp() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::rqf() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::rx0() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::rxf() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::rz0() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::rzf() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::setDerivativeOptions(Integrator &integrator, const AugOffset &offset) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::setStopTime(double tf) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::spCanEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::spEvaluate(bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::spJacF() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::spJacG() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::x0() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::xf() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::z0() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IntegratorInternal::zf() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Inverse::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Inverse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Inverse::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Inverse::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::finalize_metadata(int n, const std::map< std::string, std::vector< std::string > > &var_string_md, const std::map< std::string, std::vector< int > > &var_integer_md, const std::map< std::string, std::vector< double > > &var_numeric_md, int m, const std::map< std::string, std::vector< std::string > > &con_string_md, const std::map< std::string, std::vector< int > > &con_integer_md, const std::map< std::string, std::vector< double > > &con_numeric_md) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::freeIpopt() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::getReducedHessian() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::get_nlp_info(int &n, int &m, int &nnz_jac_g, int &nnz_h_lag) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::get_number_of_nonlinear_variables() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::get_var_con_metadata(int n, std::map< std::string, std::vector< std::string > > &var_string_md, std::map< std::string, std::vector< int > > &var_integer_md, std::map< std::string, std::vector< double > > &var_numeric_md, int m, std::map< std::string, std::vector< std::string > > &con_string_md, std::map< std::string, std::vector< int > > &con_integer_md, std::map< std::string, std::vector< double > > &con_numeric_md) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptInterface::setQPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptUserClass::finalize_metadata(Index n, const StringMetaDataMapType &var_string_md, const IntegerMetaDataMapType &var_integer_md, const NumericMetaDataMapType &var_numeric_md, Index m, const StringMetaDataMapType &con_string_md, const IntegerMetaDataMapType &con_integer_md, const NumericMetaDataMapType &con_numeric_md) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptUserClass::get_number_of_nonlinear_variables() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::IpoptUserClass::get_var_con_metadata(Index n, StringMetaDataMapType &var_string_md, IntegerMetaDataMapType &var_integer_md, NumericMetaDataMapType &var_numeric_md, Index m, StringMetaDataMapType &con_string_md, IntegerMetaDataMapType &con_integer_md, NumericMetaDataMapType &con_numeric_md) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::bjac(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::create(const Function &f, const Function &jac, const LinearSolver &linsol) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::djac(long N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::func(N_Vector u, N_Vector fval) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::kinsol_error(const std::string &module, int flag, bool fatal=true) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::lsetup(KINMem kin_mem) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::psetup(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector tmp1, N_Vector tmp2) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::psolve(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v, N_Vector tmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KinsolInterface::solveNonLinear() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KnitroInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KnitroInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::KnitroInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LapackLuDense::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LapackLuDense::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LapackLuDense::prepare() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LapackQrDense::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LapackQrDense::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LapackQrDense::prepare() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolver::spSolve(DMatrix &X, const DMatrix &B, bool transpose=false) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::colind() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::evaluateDGen(const DMatrixPtrV &input, DMatrixPtrV &output, bool tr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::evaluateMXGen(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given, bool tr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::evaluateSXGen(const SXPtrV &input, SXPtrV &output, bool tr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::getFactorization(bool transpose) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::getFactorizationSparsity(bool transpose) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::ncol() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::nnz() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::nrow() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::propagateSparsityGen(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd, bool transpose) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::row() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::solve(bool transpose) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::solve(const MX &A, const MX &B, bool transpose) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LinearSolverInternal::spSolve(DMatrix &X, const DMatrix &B, bool transpose) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LpSolverInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LpSolverInternal::solve() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LpToQp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LpToQp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::LpToQp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::T() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::at(int k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::at(int k) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::getNZ(const Matrix< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::getNZ(const Slice &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::getNZ(const std::vector< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::getNZ(int k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::getTemp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setNZ(const Matrix< int > &k, const MX &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setNZ(const Slice &k, const MX &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setNZ(const std::vector< int > &k, const MX &el) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setNZ(int k, const MX &el) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, const Matrix< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const Matrix< int > &rr, int cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const Slice &rr, const Slice &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const Sparsity &sp, int dummy) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, Slice cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, const Matrix< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, const std::vector< int > &rr, int cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, int rr, const Matrix< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, int rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setSub(const MX &m, int rr, int cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::setTemp(int t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Matrix< int > &rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Slice &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Slice &rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Slice &rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const Sparsity &sp, int dummy=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(const std::vector< int > &rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(int rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(int rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(int rr, const std::vector< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MX::sub(int rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXFunction::algorithm() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXFunction::generateLiftingFunctions(MXFunction &vdef_fcn, MXFunction &vinit_fcn) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::addDependency(const MX &dep) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, bool add=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::assign(const MX &d, const std::vector< int > &inz, const std::vector< int > &onz, bool add=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::dep(int ind=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::dep(int ind=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getAddNonzeros(const MX &y, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getAssertion(const MX &y, const std::string &fail_message="") const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getBinarySwitch(int op, const MX &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getDeterminant() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getFunction() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getFunction() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getFunctionInput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getFunctionOutput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getHorzcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getHorzsplit(const std::vector< int > &output_offset) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getInnerProd(const MX &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getInverse() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getMatrixValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getMultiplication(const MX &y, const Sparsity &sp_z=Sparsity()) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getName() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getNorm1() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getNorm2() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getNormF() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getNormInf() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getNumOutputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getOutput(int oind) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getReshape(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getSetSparse(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getSolve(const MX &r, bool tr, const LinearSolver &linear_solver) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getSubAssign(const MX &y, const Slice &i, const Slice &j) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getSubRef(const Slice &i, const Slice &j) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getTranspose() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getUnary(int op) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getVertcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::getVertsplit(const std::vector< int > &output_offset) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::hasDep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isBinaryOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isIdentity() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isMultipleOutput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isNonLinear() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isOutputNode() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isUnaryOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isValue(double val) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::isZero() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::mapping() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::nTmp(size_t &ni, size_t &nr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::ndep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::numInplace() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::numel() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::print(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::print(std::ostream &stream, long &remaining_calls) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::repr(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::setDependencies(const MX &dep) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::setDependencies(const MX &dep1, const MX &dep2, const MX &dep3) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::setDependencies(const std::vector< MX > &dep) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::setSparsity(const Sparsity &sparsity) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::shape() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::simplifyMe(MX &ex) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::size() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::size1() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::size2() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::sparsity() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MXNode::sparsity(int oind) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::append(const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::appendColumns(const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::arccos() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::arccosh() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::arcsin() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::arcsinh() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::arctan() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::arctan2(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::arctanh() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::binary(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::borBV(const Matrix< DataType > &val) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::ceil() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::className() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::clear() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::colind() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::colind(int col) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::cos() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::cosh() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::data() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::data() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::densify(const DataType &val=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::elem(int rr, int cc=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::enlarge(int nrow, int ncol, const std::vector< int > &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::erase(const std::vector< int > &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::erf() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::erfinv() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::exp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::fabs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::floor() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::fmax(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::fmin(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::fmod(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::get(DataType &val, SparsityType sp=SPARSE) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::get(Matrix< DataType > &val, SparsityType sp=SPARSE) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::get(std::vector< DataType > &val, SparsityType sp=SPARSE) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::getNZ(const Matrix< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::getNZ(const std::vector< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::getName() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::hasNonStructuralZeros() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::if_else_zero(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::inf(const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::inf(const std::pair< int, int > &rc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::inf(int nrow=1, int ncol=1) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isConstant() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isEqual(const Matrix< DataType > &ex2) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isIdentity() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isInteger() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isMinusOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isRegular() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isSmooth() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isSymbolic() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isSymbolicSparse() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::isZero() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::log() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::log10() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::logic_and(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::logic_not() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::logic_or(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::matrix_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::matrix_scalar(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::mul(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::mul_full(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::nan(const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::nan(const std::pair< int, int > &rc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::nan(int nrow=1, int ncol=1) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::printDense(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::printScalar(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::printSparse(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::printVector(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::printme(const Matrix< DataType > &y) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::remove(const std::vector< int > &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::repmat(const DataType &x, const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const std::pair< int, int > &rc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::repmat(const Matrix< DataType > &x, int nrow, int ncol=1) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::repr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::reserve(int nnz, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::resize(int nrow, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::row() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::row(int el) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sanityCheck(bool complete=false) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::scalar_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::set(DataType val, SparsityType sp=SPARSE) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::set(const Matrix< DataType > &val, SparsityType sp=SPARSE) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::set(const std::vector< DataType > &val, SparsityType sp=SPARSE) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setAll(const DataType &val) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setBV(const Matrix< DataType > &val) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setNZ(const Matrix< int > &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setNZ(const std::vector< int > &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setNZ(int k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setSparse(const Sparsity &sp, bool intersect=false) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Matrix< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Sparsity &sp, int dummy) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Matrix< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setSub(const Matrix< DataType > &m, int rr, int cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setZero() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::setZeroBV() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sign() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sin() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sinh() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sparsify(double tol=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sparsityRef() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sqrt() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sub(const Sparsity &sp, int dummy=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::sub(int rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::tan() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::tanh() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::toScalar() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::trans() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, const std::pair< int, int > &rc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, int nrow, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< DataType >::unary(int op, const Matrix< DataType > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::T() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::addSub(const Matrix< DataType > &m, RR rr, CC cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::at(int k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::at(int k) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::back() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::back() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::begin() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::begin() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::end() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::end() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::front() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::front() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::getBV(Matrix< DataType > &val) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::getSub(Matrix< DataType > &m, RR rr, CC cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr, const IndexList &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const IndexList &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const IndexList &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const Matrix< int > &rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const Slice &rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const IndexList &cc, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const IndexList &cc, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Slice &cc, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Slice &rr, const Slice &cc, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_assignment(const Sparsity &sp, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_one_based(const Matrix< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_one_based(int rr) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_one_based(int rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(int rr, const DataType &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_one_based_assignment(int rr, int cc, const DataType &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_zero_based(const Matrix< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_zero_based(int rr) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_zero_based(int rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(int rr, const DataType &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::indexed_zero_based_assignment(int rr, int cc, const DataType &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed(const IndexList &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed(const Slice &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_assignment(const IndexList &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_assignment(const Slice &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based(const Matrix< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based(int k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_one_based_assignment(int k, const DataType &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based(const Matrix< int > &k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based(int k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::nz_indexed_zero_based_assignment(int k, const DataType &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::ptr() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::ptr() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::rbegin() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::rbegin() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::rend() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::rend() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Slice &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, int cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Matrix< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Slice &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, int cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const int rr, const Slice &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Slice &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, int cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, int rr, const Matrix< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::setSub(const Matrix< DataType > &m, int rr, const std::vector< int > &cc) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const Matrix< int > &rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const Matrix< int > &rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, const std::vector< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const Slice &rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const std::vector< int > &rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(const std::vector< int > &rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(int rr, const Matrix< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(int rr, const Slice &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Matrix< T >::sub(int rr, const std::vector< int > &cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MinusInfSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MinusInfSX::isMinusInf() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MinusOneSX::getIntValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MinusOneSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MinusOneSX::isInteger() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MinusOneSX::isMinusOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MultipleOutput::getNumOutputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MultipleOutput::getOutput(int oind) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MultipleOutput::isMultipleOutput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::MultipleOutput::sparsity(int oind) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::numInplace() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Multiplication< TrX, TrY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NanSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NanSX::isNan() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Newton::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Newton::create(const Function &f, const Function &jac, const LinearSolver &linsol) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Newton::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Newton::solveNonLinear() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::checkInitialBounds() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::getGradF() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::getGradLag() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::getHessLag() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::getJacF() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::getJacG() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::getReducedHessian() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::getSpHessLag() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::gradF() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::gradLag() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::hessLag() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::jacF() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::jacG() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::reportConstraints(std::ostream &stream=std::cout) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::setOptionsFromFile(const std::string &file) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::setQPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NlpSolverInternal::spHessLag() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NonZeroIterator< DataType >::begin() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NonZeroIterator< DataType >::end() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Norm1::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Norm1::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Norm1::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Norm2::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Norm2::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Norm2::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormF::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormInf::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormInf::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::NormInf::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OldCollocationIntegrator::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OldCollocationIntegrator::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OldCollocationIntegrator::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OldCollocationIntegrator::integrate(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OldCollocationIntegrator::integrateB(double t_out) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OldCollocationIntegrator::reset() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OldCollocationIntegrator::resetB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OneSX::getIntValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OneSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OneSX::isInteger() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OneSX::isOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OoqpInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OoqpInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OoqpInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionality::checkNode() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionality::getOptionAllowedIndex(const std::string &name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionality::getOptionEnumValue(const std::string &name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionality::setOptionByAllowedIndex(const std::string &name, int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionality::setOptionByEnumValue(const std::string &name, int v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::copyOptions(const OptionsFunctionality &obj, bool skipUnknown=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::dictionary() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getBestMatches(const std::string &name, std::vector< std::string > &suggestions, int amount=5) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOption(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowed(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionAllowedIndex(const std::string &name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDefault(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionDescription(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionEnumValue(const std::string &name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionNames() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionType(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::getOptionTypeName(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::hasOption(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::hasSetOption(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::print(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::printOption(const std::string &name, std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::printOptions(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::repr(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const Dictionary &dict, bool skipUnknown=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::setOption(const std::string &str, const GenericType &val) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByAllowedIndex(const std::string &name, int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OptionsFunctionalityNode::setOptionByEnumValue(const std::string &name, int v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::getFunctionInput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::getFunctionOutput() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::getHorzcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::getVertcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::isNonLinear() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::isOutputNode() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::OutputNode::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Polynomial::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Polynomial::repr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::PrintableObject::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::PrintableObject::repr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ProfilingType() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ProfilingType< ProfilingData_ENTRY >() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ProfilingType< ProfilingData_EXIT >() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ProfilingType< ProfilingData_IO >() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ProfilingType< ProfilingData_NAME >() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ProfilingType< ProfilingData_SOURCE >() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ProfilingType< ProfilingData_TIMELINE >() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QcqpSolverInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QcqpSolverInternal::setQPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QcqpSolverInternal::solve() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QcqpToSocp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QcqpToSocp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QcqpToSocp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpSolverInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpSolverInternal::generateNativeCode(std::ostream &file) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpSolverInternal::setLPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpSolverInternal::solve() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToImplicit::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToImplicit::create(const Function &f, const Function &jac, const LinearSolver &linsol) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToImplicit::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToImplicit::solveNonLinear() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToNlp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToNlp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToNlp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToQcqp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToQcqp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpToQcqp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpoasesInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpoasesInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::QpoasesInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::RealtypeSX::getIntValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::RealtypeSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::RealtypeSX::isAlmostZero(double tol) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::getReshape(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::numInplace() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Reshape::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::RkIntegrator::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::RkIntegrator::create(const Function &f, const Function &g) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::RkIntegrator::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::RkIntegrator::setupFG() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::assignIfDuplicate(const SXElement &scalar, int depth=1) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::assignNoDelete(const SXElement &scalar) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::get() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::getTemp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::mark() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::marked() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::print(std::ostream &stream, long &remaining_calls) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::setTemp(int t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXElement::toString() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXFunction::algorithm() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::dep(int i) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::dep(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::getIntValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::getName() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::hasDep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isAlmostZero(double tol) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isConstant() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isInf() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isInteger() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isMinusInf() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isMinusOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isNan() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isOne() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isSmooth() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isSymbolic() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::isZero() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::mark() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::marked() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::ndep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::print(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SXNode::print(std::ostream &stream, long &remaining_calls) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::dualInfeasibility() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::eval_exp() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::eval_mat() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::eval_res() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::eval_vec() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::getQpSolver() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::line_search(int &ls_iter, bool &ls_success) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::primalInfeasibility() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double reg, int ls_trials, bool ls_success) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::regularize() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Scpgen::solve_qp() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdpSolverInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdpSolverInternal::setSOCPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdpSolverInternal::solve() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdqpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdqpSolverInternal::setSOCQPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdqpSolverInternal::solve() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdqpToSdp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdqpToSdp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SdqpToSdp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzeros< Add >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzeros< Add >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzeros< Add >::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzeros< Add >::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzeros< Add >::mapping() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzeros< Add >::numInplace() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice2< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::isAssignment() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosSlice< Add >::simplifyMe(MX &ex) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::getAll() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetNonzerosVector< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SetSparse::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::assertInit() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::checkNode() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::get() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::get() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::getCount() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::printPtr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::repr(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::swap(SharedObject &other) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObject::weak() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::assertInit() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::getCount() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::isInit() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::print(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::repr(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SharedObjectNode::weak() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SimpleHomotopyNlp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SimpleHomotopyNlp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SimpleHomotopyNlp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Slice::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::formatStatus(int status) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::passOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::reset() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::setOptionsFromFile(const std::string &file) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SnoptInterface::setQPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SocpSolverInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SocpSolverInternal::printProblem(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SocpSolverInternal::solve() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SocpToSdp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SocpToSdp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SocpToSdp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::getFunction() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::nTmp(size_t &ni, size_t &nr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::numInplace() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::print(std::ostream &stream, long &remaining_calls) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Solve< Tr >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::at(int k) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::at(int k) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::back() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::back() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::begin() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::begin() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::clear() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::colind() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::colind(int col) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::data() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::data() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::elem(int rr, int cc=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::end() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::end() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::front() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::front() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::getElement(int rr, int cc=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::hasNZ(int rr, int cc) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::rbegin() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::rbegin() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::rend() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::rend() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::reserve(int nnz, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::resize(int nrow, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::row() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::row(int el) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::sanityCheck(bool complete=false) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::sparsity() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::sparsityRef() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SparseStorage< DataType >::toScalar() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sparsity::T() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sparsity::colindRef() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sparsity::reCache() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sparsity::rowRef() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sparsity::shape() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Split::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Split::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Split::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Split::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Split::getNumOutputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Split::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Split::sparsity(int oind) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SqicInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SqicInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SqicInterface::generateNativeCode(std::ostream &file) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SqicInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::eval_f(const std::vector< double > &x, double &f) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::eval_g(const std::vector< double > &x, std::vector< double > &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::getQpSolver() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::getRegularization(const Matrix< double > &H) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, int ls_trials, bool ls_success) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::regularize(Matrix< double > &H, double reg) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::reset_h() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Sqpmethod::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedQpSolverInternal::checkInputs() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedQpSolverInternal::setLPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedQpSolverInternal::solve() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedQpToQp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedQpToQp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedQpToQp::generateNativeCode(std::ostream &file) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedQpToQp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqicInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqicInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqicInterface::generateNativeCode(std::ostream &file) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqicInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::eval_f(const std::vector< double > &x, double &f) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::eval_g(const std::vector< double > &x, std::vector< double > &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::eval_grad_f(const std::vector< double > &x, double &f, std::vector< double > &grad_f) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::eval_h(const std::vector< double > &x, const std::vector< double > &lambda, double sigma, Matrix< double > &H) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::eval_jac_g(const std::vector< double > &x, std::vector< double > &g, Matrix< double > &J) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::getRegularization(const Matrix< double > &H) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::getStabilizedQpSolver() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::mat_vec(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::mat_vectran(const std::vector< double > &x, const DMatrix &A, std::vector< double > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::meritfg() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::norm1matrix(const DMatrix &A) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::primalInfeasibility(const std::vector< double > &x, const std::vector< double > &lbx, const std::vector< double > &ubx, const std::vector< double > &g, const std::vector< double > &lbg, const std::vector< double > &ubg) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, double TRdelta, int ls_trials, bool ls_success, char info) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::regularize(Matrix< double > &H, double reg) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::reset_h() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::StabilizedSqp::solve_QP(const Matrix< double > &H, const std::vector< double > &g, const std::vector< double > &lbx, const std::vector< double > &ubx, const Matrix< double > &A, const std::vector< double > &lbA, const std::vector< double > &ubA, std::vector< double > &x_opt, std::vector< double > &lambda_x_opt, std::vector< double > &lambda_A_opt, double muR, const std::vector< double > &mu, const std::vector< double > &muE) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubAssign::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SubRef::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SundialsInterface::getJac() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SundialsInterface::getJacB() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SundialsInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SundialsInterface::reset() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SundialsInterface::setStopTime(double tf) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::getName() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicNLP::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicNLP::repr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicOCP::print(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicOCP::repr(std::ostream &stream=std::cout) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicQr::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicQr::evaluateSXGen(const SXPtrV &input, SXPtrV &output, bool tr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicQr::generateBody(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicQr::generateDeclarations(std::ostream &stream, const std::string &type, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicQr::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicQr::prepare() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicSX::getName() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicSX::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::SymbolicSX::isSymbolic() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::TinyXmlInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::TinyXmlInterface::parse(const std::string &filename) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::getTranspose() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::nTmp(size_t &ni, size_t &nr) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Transpose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::getBinary(int op, const MX &y, bool scX, bool scY) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::getUnary(int op) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::isUnaryOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::numInplace() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnaryMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnarySX::dep(int i) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnarySX::dep(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnarySX::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnarySX::hasDep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnarySX::isSmooth() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnarySX::ndep() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::UnarySX::print(std::ostream &stream, long &remaining_calls) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertcat::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertcat::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertcat::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertsplit::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertsplit::getOp() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertsplit::getVertcat(const std::vector< MX > &x) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::Vertsplit::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WeakRef::alive() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WeakRef::shared() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::evaluate() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::formatStatus(int status) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::init() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::passOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::reset() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::setOptionsFromFile(const std::string &file) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::WorhpInterface::setQPOptions() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlFile::parse(const std::string &filename) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlFileInternal::print(std::ostream &stream) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::checkName(const std::string &str) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::dump(std::ostream &stream, int indent=0) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::getAttribute(const std::string &attribute_name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::getName() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::getText() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::getText(T &val) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::hasAttribute(const std::string &attribute_name) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::hasChild(const std::string &childname) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::readAttribute(const std::string &attribute_name, T &val, bool assert_existance=true) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::setAttribute(const std::string &attribute_name, const std::string &attribute) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::setName(const std::string &name) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::XmlNode::size() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::clone() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getMatrixValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getReshape(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getSetSparse(const Sparsity &sp) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getTranspose() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getUnary(int op) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroByZero::printPart(std::ostream &stream, int part) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroSX::getIntValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroSX::getValue() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroSX::isAlmostZero(double tol) const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroSX::isInteger() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ZeroSX::isZero() const  {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::abs(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::acos(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::acosh(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::acosh(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::asin(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::asinh(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::asinh(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::atan(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::atan2(const T &x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::atan2(const T &x, double n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::atan2(double x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::atanh(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::atanh(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::blkdiag(const MX &A, const MX &B) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::blkdiag(const Sparsity &a, const Sparsity &b) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::blockcat(const MX &A, const MX &B, const MX &C, const MX &D) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::blockcat(const Matrix< DataType > &A, const Matrix< DataType > &B, const Matrix< DataType > &C, const Matrix< DataType > &D) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::blockmatrix(SX array[n]) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::blockmatrix(SX array[n][m]) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::casadi_load_linearsolver_csparsecholesky() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ceil(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::collocationInterpolators(const std::vector< double > &tau_root, std::vector< std::vector< double > > &C, std::vector< double > &D) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::collocationPointsL(int order, const std::string &scheme) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::constpow(const T &x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::copysign(const T &x, const T &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::copysign(double x, double y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::cos(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::cosh(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::deepcopy(const A &a) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::erf(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::erf(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::erfinv(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::erfinv(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::exp(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fabs(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::floor(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmax(const T &x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmax(const T &x, double n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmax(double x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmax(double x, double y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmax(int x, int y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmin(const T &x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmin(const T &x, double n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmin(double x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmin(double x, double y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmin(int x, int y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::fmod(const T &x, const T &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::getDescription(const std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::getPtr(Matrix< DataType > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::getPtr(const Matrix< DataType > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::getPtr(const std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::getPtr(std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::getRealTime() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::getRepresentation(const std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::get_bvec_t(const std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::get_bvec_t(const std::vector< double > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::get_bvec_t(std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::get_bvec_t(std::vector< double > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::hash_combine(std::size_t &seed, T v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::hash_combine(std::size_t &seed, const std::vector< int > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::hash_sparsity(int nrow, int ncol, const std::vector< int > &colind, const std::vector< int > &row) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::hash_value(T v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::horzcat(const MX &a, const MX &b) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::horzcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::horzcat(const Sparsity &a, const Sparsity &b) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::if_else(const SXElement &cond, const T &if_true, const T &if_false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::if_else_zero(const T &x, const T &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::if_else_zero(double x, double y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::inner_prod(const std::vector< T > &a, const std::vector< T > &b) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::is_a(const SharedObject &A) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::isinf(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::isnan(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::linspace(std::vector< T > &v, const F &first, const L &last) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::log(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::log10(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::makeVector(int size, int ind0=-1, const T &val0=T(), int ind1=-1, const T &val1=T(), int ind2=-1, const T &val2=T(), int ind3=-1, const T &val3=T(), int ind4=-1, const T &val4=T(), int ind5=-1, const T &val5=T(), int ind6=-1, const T &val6=T(), int ind7=-1, const T &val7=T(), int ind8=-1, const T &val8=T(), int ind9=-1, const T &val9=T(), int ind10=-1, const T &val10=T(), int ind11=-1, const T &val11=T(), int ind12=-1, const T &val12=T(), int ind13=-1, const T &val13=T(), int ind14=-1, const T &val14=T(), int ind15=-1, const T &val15=T(), int ind16=-1, const T &val16=T(), int ind17=-1, const T &val17=T(), int ind18=-1, const T &val18=T(), int ind19=-1, const T &val19=T()) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::matrixName< SXElement >() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::norm_1(const std::vector< T > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::norm_2(const std::vector< T > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::norm_inf(const std::vector< T > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::operation_checker(unsigned int op) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::pow(const T &x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::pow(const T &x, double n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::pow(double x, const T &n) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::print(const std::vector< T > &v, std::ostream &stream=std::cout) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::printme(const T &x, const T &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::printme(double x, double y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::profileWrite(std::ofstream &f, const T &s) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::profileWriteBare(std::ofstream &f, const T &s) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ptrVec(const std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ptrVec(const std::vector< std::vector< T > > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ptrVec(std::vector< T > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::ptrVec(std::vector< std::vector< T > > &v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::qr(const Matrix< DataType > &A, Matrix< DataType > &Q, Matrix< DataType > &R) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::range(int start, int stop, int step, int len) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::range(int stop) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::repr(const std::vector< T > &v, std::ostream &stream=std::cout) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::reshape(const MX &x, int nrow, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::shared_cast(SharedObject &A) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::shared_cast(const SharedObject &A) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::sign(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::sign(double x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::simplify(SXElement &ex) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::sin(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::sinh(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &dwork, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &eig_real, std::vector< double > &eig_imag, double num_zero) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< int > &indices, bool invert_indices=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::sq(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::sqrt(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, bool reverse=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, std::vector< MX > &ex, bool reverse=false) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::tan(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::tanh(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::toVector(const T &v0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::toVector(const T &v0, const T &v1) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::toVector(const T &v0, const T &v1, const T &v2) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::twice(const T &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::vertcat(const MX &a, const MX &b) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::vertcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  casadi::vertcat(const Sparsity &a, const Sparsity &b) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  snlog2_() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  snlog_() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  sqicDestroy() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception  sqlog_() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Assertion::Assertion(const MX &x, const MX &y, const std::string &s) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::BinaryMX< ScX, ScY >::BinaryMX(Operation op, const MX &x, const MX &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::CallFunction::CallFunction(const Function &fcn, std::vector< MX > arg) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::CollocationIntegrator::CollocationIntegrator(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Concat::Concat(const std::vector< MX > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Constant< Value >::Constant(const Sparsity &sp, Value v=Value()) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ConstantDMatrix::ConstantDMatrix(const Matrix< double > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ConstantMX::ConstantMX(const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ControlSimulatorInputIOSchemeVector< M >::ControlSimulatorInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ControlledDAEInputIOSchemeVector< M >::ControlledDAEInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::CplexInterface::CplexInterface() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::CplexInterface::CplexInterface(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::CsparseInterface::CsparseInterface(const CsparseInterface &linsol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::CsparseInterface::CsparseInterface(const Sparsity &sp, int nrhs) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::CvodesInterface::CvodesInterface(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::DAEInputIOSchemeVector< M >::DAEInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::DAEOutputIOSchemeVector< M >::DAEOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::DPLEInputIOSchemeVector< M >::DPLEInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::DPLEOutputIOSchemeVector< M >::DPLEOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::DenseMultiplication< TrX, TrY >::DenseMultiplication(const MX &z, const MX &x, const MX &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::DenseTranspose::DenseTranspose(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Determinant::Determinant(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::DsdpInterface::DsdpInterface(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::EmptySparsity::EmptySparsity() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::FixedStepIntegrator::FixedStepIntegrator(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::GetNonzeros::GetNonzeros(const Sparsity &sp, const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::GetNonzerosSlice2::GetNonzerosSlice2(const Sparsity &sp, const MX &x, const Slice &inner, const Slice &outer) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::GetNonzerosSlice::GetNonzerosSlice(const Sparsity &sp, const MX &x, const Slice &s) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::GetNonzerosVector::GetNonzerosVector(const Sparsity &sp, const MX &x, const std::vector< int > &nz) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::GradFInputIOSchemeVector< M >::GradFInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::GradFOutputIOSchemeVector< M >::GradFOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::HNLPInputIOSchemeVector< M >::HNLPInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::HessLagInputIOSchemeVector< M >::HessLagInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::HessLagOutputIOSchemeVector< M >::HessLagOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Horzcat::Horzcat(const std::vector< MX > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Horzsplit::Horzsplit(const MX &x, const std::vector< int > &offset) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IdasInterface::IdasInterface(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IndexList::IndexList() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IndexList::IndexList(const Slice &i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IndexList::IndexList(const std::vector< int > &i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IndexList::IndexList(int i) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::InfSX::InfSX() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::InnerProd::InnerProd(const MX &x, const MX &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IntegratorInputIOSchemeVector< M >::IntegratorInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IntegratorOutputIOSchemeVector< M >::IntegratorOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Inverse::Inverse(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::IpoptInterface::IpoptInterface(const Function &nlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::JacGInputIOSchemeVector< M >::JacGInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::JacGOutputIOSchemeVector< M >::JacGOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::KinsolInterface::KinsolInterface(const Function &f, const Function &jac, const LinearSolver &linsol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::KnitroInterface::KnitroInterface(const Function &nlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LPStructIOSchemeVector< T >::LPStructIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LapackLuDense::LapackLuDense(const Sparsity &sparsity, int nrhs) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LapackQrDense::LapackQrDense(const Sparsity &sparsity, int nrhs) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LinearSolver::LinearSolver() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LinsolInputIOSchemeVector< M >::LinsolInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LinsolOutputIOSchemeVector< M >::LinsolOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LpSolverInputIOSchemeVector< M >::LpSolverInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LpSolverOutputIOSchemeVector< M >::LpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::LpToQp::LpToQp(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const MX &output) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::MXFunction::MXFunction(const MX &input, const std::vector< MX > &output) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::MXFunction::MXFunction(const std::vector< MX > &input, const MX &output) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::MXNode::MXNode() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix(const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const DataType &val=DataType(0)) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const std::vector< DataType > &d) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< DataType > &x, int nrow, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix(const std::vector< std::vector< DataType > > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Matrix< DataType >::Matrix(double val) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::MinusInfSX::MinusInfSX() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::MinusOneSX::MinusOneSX() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::MultipleOutput::MultipleOutput() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Multiplication< TrX, TrY >::Multiplication(const MX &z, const MX &x, const MX &y) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NLPInputIOSchemeVector< M >::NLPInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NLPOutputIOSchemeVector< M >::NLPOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NanSX::NanSX() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Newton::Newton(const Function &f, const Function &jac, const LinearSolver &linsol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NlpSolverInputIOSchemeVector< M >::NlpSolverInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NlpSolverOutputIOSchemeVector< M >::NlpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NonZeroIterator< DataType >::NonZeroIterator(const Matrix< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Norm1::Norm1(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Norm2::Norm2(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Norm::Norm(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NormF::NormF(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::NormInf::NormInf(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::OldCollocationIntegrator::OldCollocationIntegrator(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::OneSX::OneSX() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::OoqpInterface::OoqpInterface() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::OoqpInterface::OoqpInterface(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::OptionsFunctionalityNode::OptionsFunctionalityNode() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::OutputNode::OutputNode(const MX &parent, int oind) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QCQPStructIOSchemeVector< T >::QCQPStructIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QPStructIOSchemeVector< T >::QPStructIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QcqpSolverInputIOSchemeVector< M >::QcqpSolverInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QcqpSolverOutputIOSchemeVector< M >::QcqpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QcqpToSocp::QcqpToSocp(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpSolverInputIOSchemeVector< M >::QpSolverInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpSolverOutputIOSchemeVector< M >::QpSolverOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpToImplicit::QpToImplicit(const Function &f, const Function &jac, const LinearSolver &linsol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpToNlp::QpToNlp() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpToNlp::QpToNlp(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpToQcqp::QpToQcqp(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpoasesInterface::QpoasesInterface() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::QpoasesInterface::QpoasesInterface(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::RDAEInputIOSchemeVector< M >::RDAEInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::RDAEOutputIOSchemeVector< M >::RDAEOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Reshape::Reshape(const MX &x, Sparsity sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::RkIntegrator::RkIntegrator(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::RuntimeConst< T >::RuntimeConst() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::RuntimeConst< T >::RuntimeConst(T v) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SDPInputIOSchemeVector< M >::SDPInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SDPOutputIOSchemeVector< M >::SDPOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SDPStructIOSchemeVector< T >::SDPStructIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SDQPInputIOSchemeVector< M >::SDQPInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SDQPOutputIOSchemeVector< M >::SDQPOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SDQPStructIOSchemeVector< T >::SDQPStructIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SOCPInputIOSchemeVector< M >::SOCPInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SOCPOutputIOSchemeVector< M >::SOCPOutputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SOCPStructIOSchemeVector< T >::SOCPStructIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXElement::SXElement(const SXElement &scalar) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const SX &res) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< SX > &res) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXFunction::SXFunction(const SX &arg, const std::vector< std::vector< SXElement > > &res) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXFunction::SXFunction(const std::vector< SX > &arg, const SX &res) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const SX &res) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const std::vector< std::vector< SXElement > > &res) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SXNode::SXNode() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ScalarSparseSparsity::ScalarSparseSparsity() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ScalarSparsity::ScalarSparsity() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Scpgen::Scpgen(const Function &nlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SdqpToSdp::SdqpToSdp(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SetNonzeros< Add >::SetNonzeros(const MX &y, const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SetNonzerosSlice2< Add >::SetNonzerosSlice2(const MX &y, const MX &x, const Slice &inner, const Slice &outer) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SetNonzerosSlice< Add >::SetNonzerosSlice(const MX &y, const MX &x, const Slice &s) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SetNonzerosVector< Add >::SetNonzerosVector(const MX &y, const MX &x, const std::vector< int > &nz) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SetSparse::SetSparse(const MX &x, const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SharedObject::SharedObject() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SharedObject::SharedObject(const SharedObject &ref) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SharedObjectNode::SharedObjectNode() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SharedObjectNode::SharedObjectNode(const SharedObjectNode &node) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SimpleHomotopyNlp::SimpleHomotopyNlp(const Function &hnlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SnoptInterface::SnoptInterface(const Function &nlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SocpToSdp::SocpToSdp(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Solve< Tr >::Solve(const MX &r, const MX &A, const LinearSolver &linear_solver) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< A > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< DataType > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const DataType &val=DataType(0)) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const std::vector< DataType > &d) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< A > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< A > &x, int nrow, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< DataType > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< DataType > &x, int nrow, int ncol) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const std::vector< std::vector< DataType > > &m) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Split::Split(const MX &x, const std::vector< int > &offset) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SqicInterface::SqicInterface() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SqicInterface::SqicInterface(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Sqpmethod::Sqpmethod(const Function &nlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::StabilizedQpSolverInputIOSchemeVector< M >::StabilizedQpSolverInputIOSchemeVector(const std::vector< M > &t) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::StabilizedQpToQp::StabilizedQpToQp(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::StabilizedSqicInterface::StabilizedSqicInterface(const std::vector< Sparsity > &st) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::StabilizedSqp::StabilizedSqp(const Function &nlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SubAssign::SubAssign(const MX &x, const MX &y, const Slice &i, const Slice &j) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SubRef::SubRef(const MX &x, const Slice &i, const Slice &j) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SundialsInterface::SundialsInterface(const Function &f, const Function &g) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, const Sparsity &sp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, int nrow=1, int ncol=1) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SymbolicQr::SymbolicQr(const Sparsity &sparsity, int nrhs) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::SymbolicSX::SymbolicSX(const std::string &name) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::TinyXmlInterface::TinyXmlInterface() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Transpose::Transpose(const MX &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::UnaryMX::UnaryMX(Operation op, MX x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Vertcat::Vertcat(const std::vector< MX > &x) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::Vertsplit::Vertsplit(const MX &x, const std::vector< int > &offset) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::WeakRef::WeakRef(SharedObject shared) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::WeakRef::WeakRef(int dummy=0) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::WorhpInterface::WorhpInterface(const Function &nlp) {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::XmlNode::XmlNode() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}
%exception casadi::ZeroSX::ZeroSX() {
 try { INTERNAL_MSG() $action } CATCH_OR_RETHROW 
}