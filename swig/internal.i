%exception  CasADi::Assertion::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Assertion::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Assertion::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Assertion::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Assertion::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Assertion::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Assertion::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::getBinary(int op, const MX &y, bool scX, bool scY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::isBinaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinaryMX< ScX, ScY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinarySX::dep(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinarySX::dep(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinarySX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinarySX::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinarySX::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinarySX::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::BinarySX::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::getFunction() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::getFunctionInput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::getFunctionOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CallFunction::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::addAuxiliary(Auxiliary f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::addDependency(const Function &f) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::addInclude(const std::string &new_include, bool relative_path=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::addSparsity(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::casadi_dot(int n, const std::string &x, int inc_x, const std::string &y, int inc_y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::copyVector(std::ostream &s, const std::string &arg, std::size_t n, const std::string &res, const std::string &it="i", bool only_if_exists=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::flush(std::ostream &s) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::getConstant(const std::vector< double > &v, bool allow_adding=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::getConstant(const std::vector< int > &v, bool allow_adding=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::getDependency(const Function &f) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::CodeGenerator::getSparsity(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Concat::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Concat::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Concat::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Concat::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Concat::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Concat::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getSetSparse(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::isValue(double val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Constant< Value >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantDMatrix::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::getInnerProd(const MX &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantSX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ConstantSX::isConstant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseMultiplication< TrX, TrY >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseMultiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseTranspose::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseTranspose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseTranspose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseTranspose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseTranspose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseTranspose::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::DenseTranspose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Determinant::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Determinant::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Determinant::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Determinant::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::callDerivative(const DMatrixVector &arg, DMatrixVector &output_res, const DMatrixVectorVector &fseed, DMatrixVectorVector &output_fsens, const DMatrixVectorVector &aseed, DMatrixVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::callDerivative(const MXVector &arg, MXVector &output_res, const MXVectorVector &fseed, MXVectorVector &output_fsens, const MXVectorVector &aseed, MXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::callDerivative(const SXVector &arg, SXVector &output_res, const SXVectorVector &fseed, SXVectorVector &output_fsens, const SXVectorVector &aseed, SXVectorVector &output_asens, bool always_inline=false, bool never_inline=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::checkInputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::checkNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::inputScheme() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::inputScheme() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::input_struct() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::input_struct() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::outputScheme() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::outputScheme() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::output_struct() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::output_struct() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::spCanEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::spEvaluate(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Function::spInit(bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericMatrix< MX  >::shape() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericMatrix< MatType >::shape() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericMatrix< Matrix< DataType >  >::shape() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericType::is_a() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericType::toDictionary() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericType::toDictionary() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericType::toDoubleVector() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericType::toDoubleVector() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericType::toIntVector() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GenericType::toIntVector() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzeros::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzeros::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzeros::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzeros::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzeros::mapping() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice2::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosSlice::simplifyMe(MX &ex) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::GetNonzerosVector::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzcat::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzcat::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzcat::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzsplit::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzsplit::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzsplit::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Horzsplit::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::getInput(T val, const std::string &iname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::getInput(T val, int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::getOutput(T val, const std::string &oname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::getOutput(T val, int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::inputS(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::inputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::inputSchemeEntry(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::outputS(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::outputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::outputSchemeEntry(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Derived >::schemeEntry(const CasADi::IOScheme &scheme, const std::string &name, bool input) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::getInput(T val, const std::string &iname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::getInput(T val, int iind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::getOutput(T val, const std::string &oname) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::getOutput(T val, int oind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::inputS(int i) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::inputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::inputSchemeEntry(const std::string &name) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::outputS(int i) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::outputS(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::outputSchemeEntry(const std::string &name) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOInterface< Function  >::schemeEntry(const CasADi::IOScheme &scheme, const std::string &name, bool input) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOScheme::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOScheme::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOSchemeVector< M  >::print(std::ostream &stream=std::cout) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOSchemeVector< M  >::repr(std::ostream &stream=std::cout) const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOSchemeVector< M  >::vector() const {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOSchemeVector< T >::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IOSchemeVector< T >::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IndexList::getAll(int len) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InfSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InfSX::isInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::InnerProd::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IntegerSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IntegerSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IntegerSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Inverse::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Inverse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Inverse::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Inverse::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IpoptUserClass::finalize_metadata(Index n, const StringMetaDataMapType &var_string_md, const IntegerMetaDataMapType &var_integer_md, const NumericMetaDataMapType &var_numeric_md, Index m, const StringMetaDataMapType &con_string_md, const IntegerMetaDataMapType &con_integer_md, const NumericMetaDataMapType &con_numeric_md) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IpoptUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IpoptUserClass::get_number_of_nonlinear_variables() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::IpoptUserClass::get_var_con_metadata(Index n, StringMetaDataMapType &var_string_md, IntegerMetaDataMapType &var_integer_md, NumericMetaDataMapType &var_numeric_md, Index m, StringMetaDataMapType &con_string_md, IntegerMetaDataMapType &con_integer_md, NumericMetaDataMapType &con_numeric_md) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::LinearSolver::spSolve(DMatrix &X, const DMatrix &B, bool transpose=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::T() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::at(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::at(int k) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::getNZ(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::getNZ(const Slice &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::getNZ(const std::vector< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::getNZ(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::getTemp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setNZ(const Matrix< int > &k, const MX &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setNZ(const Slice &k, const MX &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setNZ(const std::vector< int > &k, const MX &el) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setNZ(int k, const MX &el) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const Matrix< int > &k) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const Matrix< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const Matrix< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const Slice &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const Sparsity &sp, int dummy) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const std::vector< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, const std::vector< int > &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, int rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setSub(const MX &m, int rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::setTemp(int t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Matrix< int > &k, int dummy=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Matrix< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Slice &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Slice &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Slice &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const Sparsity &sp, int dummy=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(const std::vector< int > &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(int rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(int rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MX::sub(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXFunction::algorithm() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXFunction::generateLiftingFunctions(MXFunction &vdef_fcn, MXFunction &vinit_fcn) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::addDependency(const MX &dep) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::assign(const MX &d, const std::vector< int > &inz, bool add=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::assign(const MX &d, const std::vector< int > &inz, const std::vector< int > &onz, bool add=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::dep(int ind=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::dep(int ind=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getAddNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getAssertion(const MX &y, const std::string &fail_message="") const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getBinary(int op, const MX &y, bool scX, bool scY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getBinarySwitch(int op, const MX &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getDeterminant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getFunction() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getFunction() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getFunctionInput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getFunctionOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getHorzsplit(const std::vector< int > &output_offset) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getInnerProd(const MX &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getInverse() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getMultiplication(const MX &y, const Sparsity &sp_z=Sparsity()) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getNorm1() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getNorm2() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getNormF() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getNormInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getOutput(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getSetSparse(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getSolve(const MX &r, bool tr, const LinearSolver &linear_solver) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getSubAssign(const MX &y, const Slice &i, const Slice &j) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getSubRef(const Slice &i, const Slice &j) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::getVertsplit(const std::vector< int > &output_offset) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isBinaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isMultipleOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isNonLinear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isOutputNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isUnaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isValue(double val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::mapping() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::numel() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::setDependencies(const MX &dep) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::setDependencies(const MX &dep1, const MX &dep2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::setDependencies(const MX &dep1, const MX &dep2, const MX &dep3) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::setDependencies(const std::vector< MX > &dep) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::setSparsity(const Sparsity &sparsity) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::shape() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::simplifyMe(MX &ex) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::size() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::size1() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::size2() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::sparsity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MXNode::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::append(const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::appendColumns(const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::arccos() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::arccosh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::arcsin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::arcsinh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::arctan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::arctan2(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::arctanh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::binary(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::borBV(const Matrix< DataType > &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::ceil() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::className() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::clear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::colind() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::colind(int col) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::cos() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::cosh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::data() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::data() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::densify(const DataType &val=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::elem(int rr, int cc=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::elem(int rr, int cc=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::enlarge(int nrow, int ncol, const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::erase(const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::erf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::erfinv() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::exp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::fabs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::floor() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::fmax(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::fmin(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::get(DataType &val, SparsityType sp=SPARSE) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::get(Matrix< DataType > &val, SparsityType sp=SPARSE) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::get(std::vector< DataType > &val, SparsityType sp=SPARSE) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::getNZ(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::getNZ(const std::vector< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::hasNonStructuralZeros() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::if_else_zero(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::inf(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::inf(const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::inf(int nrow=1, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isConstant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isEqual(const Matrix< DataType > &ex2) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isIdentity() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isRegular() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isSymbolic() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isSymbolicSparse() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::log() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::log10() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::logic_and(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::logic_not() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::logic_or(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::matrix_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::matrix_scalar(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::mul(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::mul_full(const Matrix< DataType > &y, const Sparsity &sp_z=Sparsity()) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::nan(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::nan(const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::nan(int nrow=1, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::printDense(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::printScalar(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::printSparse(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::printVector(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::printme(const Matrix< DataType > &y) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::remove(const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::repmat(const DataType &x, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::repmat(const Matrix< DataType > &x, const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::repmat(const Matrix< DataType > &x, int nrow, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::reserve(int nnz) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::reserve(int nnz, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::resize(int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::row() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::row(int el) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sanityCheck(bool complete=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::scalar_matrix(int op, const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::set(DataType val, SparsityType sp=SPARSE) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::set(const Matrix< DataType > &val, SparsityType sp=SPARSE) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::set(const std::vector< DataType > &val, SparsityType sp=SPARSE) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setAll(const DataType &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setBV(const Matrix< DataType > &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setNZ(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setNZ(const std::vector< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setNZ(int k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setSparse(const Sparsity &sp, bool intersect=false) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const Sparsity &sp, int dummy) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setSub(const Matrix< DataType > &m, int rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setZero() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::setZeroBV() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sign() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sinh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sparsify(double tol=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sparsityRef() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sqrt() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sub(const Matrix< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sub(const Matrix< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sub(const Sparsity &sp, int dummy=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sub(const std::vector< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sub(const std::vector< int > &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::sub(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::tan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::tanh() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::toScalar() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::trans() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, const std::pair< int, int > &rc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::triplet(const std::vector< int > &row, const std::vector< int > &col, const std::vector< DataType > &d, int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< DataType >::unary(int op, const Matrix< DataType > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::T() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::addSub(const Matrix< DataType > &m, RR rr, CC cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::at(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::at(int k) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::back() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::back() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::begin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::begin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::end() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::end() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::front() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::front() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::getBV(Matrix< DataType > &val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::getSub(Matrix< DataType > &m, RR rr, CC cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const IndexList &rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const IndexList &rr, const IndexList &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const IndexList &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const Matrix< int > &rr, const IndexList &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const Matrix< int > &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const Matrix< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const Slice &rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const Slice &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const Slice &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const IndexList &rr, const IndexList &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const IndexList &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const IndexList &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const Matrix< int > &rr, const Slice &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const Slice &rr, const Matrix< int > &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const Slice &rr, const Slice &cc, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_assignment(const Sparsity &sp, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_one_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_one_based(int rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_one_based(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_one_based_assignment(int rr, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_one_based_assignment(int rr, int cc, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_zero_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_zero_based(int rr) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_zero_based(int rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_zero_based_assignment(int rr, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::indexed_zero_based_assignment(int rr, int cc, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed(const IndexList &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed(const Slice &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_assignment(const IndexList &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_assignment(const Slice &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_one_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_one_based(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_one_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_one_based_assignment(int k, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_zero_based(const Matrix< int > &k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_zero_based(int k) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_zero_based_assignment(const Matrix< int > &k, const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::nz_indexed_zero_based_assignment(int k, const DataType &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::ptr() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::ptr() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::rbegin() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::rbegin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::rend() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::rend() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const Matrix< int > &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Matrix< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const Slice &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const int rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, const Slice &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, const std::vector< int > &rr, int cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::setSub(const Matrix< DataType > &m, int rr, const std::vector< int > &cc) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(const Matrix< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(const Slice &rr, const Matrix< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(const Slice &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(const Slice &rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(const Slice &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(const std::vector< int > &rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(const std::vector< int > &rr, int cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(int rr, const Slice &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Matrix< T >::sub(int rr, const std::vector< int > &cc) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MinusInfSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MinusInfSX::isMinusInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MinusOneSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MinusOneSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MinusOneSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MinusOneSX::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MultipleOutput::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MultipleOutput::getOutput(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MultipleOutput::isMultipleOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::MultipleOutput::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Multiplication< TrX, TrY >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NanSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NanSX::isNan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NonZeroIterator< DataType >::begin() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NonZeroIterator< DataType >::end() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Norm1::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Norm1::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Norm1::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Norm2::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Norm2::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Norm2::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormF::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormInf::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormInf::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::NormInf::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OneSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OneSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OneSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OneSX::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionality::checkNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionality::getOptionAllowedIndex(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionality::getOptionEnumValue(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionality::setOptionByAllowedIndex(const std::string &name, int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionality::setOptionByEnumValue(const std::string &name, int v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::copyOptions(const OptionsFunctionality &obj, bool skipUnknown=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::dictionary() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getBestMatches(const std::string &name, std::vector< std::string > &suggestions, int amount=5) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOption(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionAllowed(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionAllowedIndex(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionDefault(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionDescription(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionEnumValue(const std::string &name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionNames() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionType(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::getOptionTypeName(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::hasOption(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::hasSetOption(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::printOption(const std::string &name, std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::printOptions(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::setOption(const Dictionary &dict, bool skipUnknown=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::setOption(const std::string &str, const GenericType &val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::setOptionByAllowedIndex(const std::string &name, int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OptionsFunctionalityNode::setOptionByEnumValue(const std::string &name, int v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::getFunctionInput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::getFunctionOutput() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::getHorzcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::isNonLinear() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::isOutputNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::OutputNode::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Polynomial::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Polynomial::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::PrintableObject::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::PrintableObject::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::RealtypeSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::RealtypeSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::RealtypeSX::isAlmostZero(double tol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Reshape::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::assignIfDuplicate(const SXElement &scalar, int depth=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::assignNoDelete(const SXElement &scalar) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::get() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::getTemp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::mark() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::marked() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::setTemp(int t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXElement::toString() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXFunction::algorithm() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::dep(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::dep(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isAlmostZero(double tol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isConstant() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isMinusInf() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isMinusOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isNan() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isOne() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isSymbolic() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::mark() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::marked() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SXNode::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzeros< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzeros< Add >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzeros< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzeros< Add >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzeros< Add >::mapping() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzeros< Add >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice2< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::isAssignment() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosSlice< Add >::simplifyMe(MX &ex) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::getAll() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetNonzerosVector< Add >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SetSparse::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::assertInit() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::checkNode() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::get() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::get() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::getCount() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::printPtr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::swap(SharedObject &other) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObject::weak() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::assertInit() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::getCount() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::init() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::isInit() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::print(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::repr(std::ostream &stream) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SharedObjectNode::weak() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Slice::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::getFunction() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Solve< Tr >::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Sparsity::T() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Sparsity::colindRef() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Sparsity::reCache() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Sparsity::rowRef() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Sparsity::shape() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Split::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Split::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Split::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Split::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Split::getNumOutputs() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Split::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Split::sparsity(int oind) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubAssign::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SubRef::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicNLP::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicNLP::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicOCP::print(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicOCP::repr(std::ostream &stream=std::cout) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicSX::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicSX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::SymbolicSX::isSymbolic() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::evaluateGen(const MatV &input, MatV &output, std::vector< int > &itmp, std::vector< T > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::nTmp(size_t &ni, size_t &nr) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Transpose::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::getBinary(int op, const MX &y, bool scX, bool scY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::isUnaryOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::numInplace() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnaryMX::propagateSparsity(DMatrixPtrV &input, DMatrixPtrV &output, bool fwd) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnarySX::dep(int i) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnarySX::dep(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnarySX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnarySX::hasDep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnarySX::isSmooth() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnarySX::ndep() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::UnarySX::print(std::ostream &stream, long &remaining_calls) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertcat::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertcat::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertcat::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertcat::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertsplit::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertsplit::evaluateMX(const MXPtrV &input, MXPtrV &output, const MXPtrVV &fwdSeed, MXPtrVV &fwdSens, const MXPtrVV &adjSeed, MXPtrVV &adjSens, bool output_given) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertsplit::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertsplit::getVertcat(const std::vector< MX > &x) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::Vertsplit::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::WeakRef::alive() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::WeakRef::shared() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::checkName(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::dump(std::ostream &stream, int indent=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::getAttribute(const std::string &attribute_name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::getText() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::getText(T &val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::hasAttribute(const std::string &attribute_name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::hasChild(const std::string &childname) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::readAttribute(const std::string &attribute_name, T &val, bool assert_existance=true) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::setAttribute(const std::string &attribute_name, const std::string &attribute) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::setName(const std::string &name) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::XMLNode::size() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::clone() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::evaluateD(const DMatrixPtrV &input, DMatrixPtrV &output, std::vector< int > &itmp, std::vector< double > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::evaluateSX(const SXPtrV &input, SXPtrV &output, std::vector< int > &itmp, std::vector< SXElement > &rtmp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::generateOperation(std::ostream &stream, const std::vector< std::string > &arg, const std::vector< std::string > &res, CodeGenerator &gen) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getBinary(int op, const MX &y, bool ScX, bool ScY) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getGetNonzeros(const Sparsity &sp, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getMatrixValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getReshape(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getSetNonzeros(const MX &y, const std::vector< int > &nz) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getSetSparse(const Sparsity &sp) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getTranspose() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getUnary(int op) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroByZero::printPart(std::ostream &stream, int part) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroSX::getIntValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroSX::getValue() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroSX::isAlmostZero(double tol) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroSX::isInteger() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ZeroSX::isZero() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::abs(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::acos(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::acosh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::acosh(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::asin(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::asinh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::asinh(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::atan(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::atan2(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::atan2(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::atan2(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::atanh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::atanh(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::blockcat(const Matrix< DataType > &A, const Matrix< DataType > &B, const Matrix< DataType > &C, const Matrix< DataType > &D) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::blockmatrix(SX array[n]) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::blockmatrix(SX array[n][m]) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ceil(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::collocationInterpolators(const std::vector< double > &tau_root, std::vector< std::vector< double > > &C, std::vector< double > &D) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::constpow(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::copysign(const T &x, const T &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::copysign(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::cos(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::cosh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::deepcopy(const A &a) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::erf(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::erf(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::erfinv(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::erfinv(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::exp(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fabs(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::floor(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmax(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmax(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmax(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmax(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmax(int x, int y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmin(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmin(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmin(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmin(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::fmin(int x, int y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::getDescription(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::getPtr(Matrix< DataType > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::getPtr(const Matrix< DataType > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::getPtr(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::getPtr(std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::getRepresentation(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::get_bvec_t(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::get_bvec_t(std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::hash_combine(std::size_t &seed, T v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::hash_combine(std::size_t &seed, const std::vector< int > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::hash_value(T v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::horzcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::if_else(const SXElement &cond, const T &if_true, const T &if_false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::if_else_zero(const T &x, const T &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::if_else_zero(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::inner_prod(const std::vector< T > &a, const std::vector< T > &b) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::isEqual(const MX &ex1, const MX &ex2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::is_a(const SharedObject &A) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::isinf(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::isnan(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::linspace(std::vector< T > &v, const F &first, const L &last) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::log(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::log10(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::makeVector(int size, int ind0=-1, const T &val0=T(), int ind1=-1, const T &val1=T(), int ind2=-1, const T &val2=T(), int ind3=-1, const T &val3=T(), int ind4=-1, const T &val4=T(), int ind5=-1, const T &val5=T(), int ind6=-1, const T &val6=T(), int ind7=-1, const T &val7=T(), int ind8=-1, const T &val8=T(), int ind9=-1, const T &val9=T(), int ind10=-1, const T &val10=T(), int ind11=-1, const T &val11=T(), int ind12=-1, const T &val12=T(), int ind13=-1, const T &val13=T(), int ind14=-1, const T &val14=T(), int ind15=-1, const T &val15=T(), int ind16=-1, const T &val16=T(), int ind17=-1, const T &val17=T(), int ind18=-1, const T &val18=T(), int ind19=-1, const T &val19=T()) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::matrixName< SXElement >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::norm_1(const std::vector< T > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::norm_2(const std::vector< T > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::norm_inf(const std::vector< T > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::operation_checker(unsigned int op) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::pow(const T &x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::pow(const T &x, double n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::pow(double x, const T &n) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::print(const std::vector< T > &v, std::ostream &stream=std::cout) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::printme(const T &x, const T &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::printme(double x, double y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ptrVec(const std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ptrVec(const std::vector< std::vector< T > > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ptrVec(std::vector< T > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::ptrVec(std::vector< std::vector< T > > &v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::qr(const Matrix< DataType > &A, Matrix< DataType > &Q, Matrix< DataType > &R) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::repr(const std::vector< T > &v, std::ostream &stream=std::cout) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::shared_cast(SharedObject &A) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::shared_cast(const SharedObject &A) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::sign(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::sign(double x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::sin(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::sinh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &dwork, std::vector< double > &eig_real, std::vector< double > &eig_imag) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &eig_real, std::vector< double > &eig_imag) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< int > &indices, bool invert_indices=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::sq(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::sqrt(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, bool reverse=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::substituteInPlace(const std::vector< MX > &v, std::vector< MX > &vdef, std::vector< MX > &ex, bool reverse=false) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::tan(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::tanh(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::toVector(const T &v0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::toVector(const T &v0, const T &v1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::toVector(const T &v0, const T &v1, const T &v2) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::twice(const T &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  CasADi::vertcat(const Matrix< DataType > &x, const Matrix< DataType > &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  ProfilingType() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  ProfilingType< ProfilingData_ENTRY >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  ProfilingType< ProfilingData_EXIT >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  ProfilingType< ProfilingData_IO >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  ProfilingType< ProfilingData_NAME >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  ProfilingType< ProfilingData_SOURCE >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  ProfilingType< ProfilingData_TIMELINE >() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  profileWrite(std::ofstream &f, const T &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  profileWriteBare(std::ofstream &f, const T &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  sqicDestroy() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Assertion::Assertion(const MX &x, const MX &y, const std::string &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::BinaryMX< ScX, ScY >::BinaryMX(Operation op, const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::CallFunction::CallFunction(const Function &fcn, std::vector< MX > arg) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Concat::Concat(const std::vector< MX > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Constant< Value >::Constant(const Sparsity &sp, Value v=Value()) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::ConstantDMatrix::ConstantDMatrix(const Matrix< double > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::ConstantMX::ConstantMX(const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::ControlSimulatorInputIOSchemeVector< M >::ControlSimulatorInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::ControlledDAEInputIOSchemeVector< M >::ControlledDAEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::DAEInputIOSchemeVector< M >::DAEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::DAEOutputIOSchemeVector< M >::DAEOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::DPLEInputIOSchemeVector< M >::DPLEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::DPLEOutputIOSchemeVector< M >::DPLEOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::DenseMultiplication< TrX, TrY >::DenseMultiplication(const MX &z, const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::DenseTranspose::DenseTranspose(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Determinant::Determinant(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::EmptySparsity::EmptySparsity() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::GetNonzeros::GetNonzeros(const Sparsity &sp, const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::GetNonzerosSlice2::GetNonzerosSlice2(const Sparsity &sp, const MX &x, const Slice &inner, const Slice &outer) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::GetNonzerosSlice::GetNonzerosSlice(const Sparsity &sp, const MX &x, const Slice &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::GetNonzerosVector::GetNonzerosVector(const Sparsity &sp, const MX &x, const std::vector< int > &nz) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::GradFInputIOSchemeVector< M >::GradFInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::GradFOutputIOSchemeVector< M >::GradFOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::HNLPInputIOSchemeVector< M >::HNLPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::HessLagInputIOSchemeVector< M >::HessLagInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::HessLagOutputIOSchemeVector< M >::HessLagOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Horzcat::Horzcat(const std::vector< MX > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Horzsplit::Horzsplit(const MX &x, const std::vector< int > &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::IndexList::IndexList() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::IndexList::IndexList(const Slice &i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::IndexList::IndexList(const std::vector< int > &i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::IndexList::IndexList(int i) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::InfSX::InfSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::InnerProd::InnerProd(const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::IntegratorInputIOSchemeVector< M >::IntegratorInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::IntegratorOutputIOSchemeVector< M >::IntegratorOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Inverse::Inverse(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::JacGInputIOSchemeVector< M >::JacGInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::JacGOutputIOSchemeVector< M >::JacGOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::LPSolverInputIOSchemeVector< M >::LPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::LPSolverOutputIOSchemeVector< M >::LPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::LPStructIOSchemeVector< T >::LPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::LinearSolver::LinearSolver() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::LinsolInputIOSchemeVector< M >::LinsolInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::LinsolOutputIOSchemeVector< M >::LinsolOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MXFunction::MXFunction(const MX &input, const MX &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MXFunction::MXFunction(const MX &input, const std::vector< MX > &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MXFunction::MXFunction(const std::vector< MX > &input, const MX &output) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MXNode::MXNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const DataType &val=DataType(0)) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(const Sparsity &sparsity, const std::vector< DataType > &d) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(const std::vector< DataType > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(const std::vector< DataType > &x, int nrow, int ncol) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(const std::vector< std::vector< DataType > > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Matrix< DataType >::Matrix(double val) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MayerInputIOSchemeVector< M >::MayerInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MinusInfSX::MinusInfSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MinusOneSX::MinusOneSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::MultipleOutput::MultipleOutput() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Multiplication< TrX, TrY >::Multiplication(const MX &z, const MX &x, const MX &y) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NLPInputIOSchemeVector< M >::NLPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NLPOutputIOSchemeVector< M >::NLPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NLPSolverInputIOSchemeVector< M >::NLPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NLPSolverOutputIOSchemeVector< M >::NLPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NanSX::NanSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NonZeroIterator< DataType >::NonZeroIterator(const Matrix< DataType > &m) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Norm1::Norm1(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Norm2::Norm2(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Norm::Norm(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NormF::NormF(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::NormInf::NormInf(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::OCPInputIOSchemeVector< M >::OCPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::OCPOutputIOSchemeVector< M >::OCPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::OneSX::OneSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::OptionsFunctionalityNode::OptionsFunctionalityNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::OutputNode::OutputNode(const MX &parent, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::QCQPSolverInputIOSchemeVector< M >::QCQPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::QCQPSolverOutputIOSchemeVector< M >::QCQPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::QCQPStructIOSchemeVector< T >::QCQPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::QPSolverInputIOSchemeVector< M >::QPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::QPSolverOutputIOSchemeVector< M >::QPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::QPStructIOSchemeVector< T >::QPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::RDAEInputIOSchemeVector< M >::RDAEInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::RDAEOutputIOSchemeVector< M >::RDAEOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Reshape::Reshape(const MX &x, Sparsity sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::RuntimeConst< T >::RuntimeConst() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::RuntimeConst< T >::RuntimeConst(T v) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SDPInputIOSchemeVector< M >::SDPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SDPOutputIOSchemeVector< M >::SDPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SDPStructIOSchemeVector< T >::SDPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SDQPInputIOSchemeVector< M >::SDQPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SDQPOutputIOSchemeVector< M >::SDQPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SDQPStructIOSchemeVector< T >::SDQPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SOCPInputIOSchemeVector< M >::SOCPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SOCPOutputIOSchemeVector< M >::SOCPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SOCPStructIOSchemeVector< T >::SOCPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXElement::SXElement(const SXElement &scalar) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXFunction::SXFunction(const SX &arg, const SX &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXFunction::SXFunction(const SX &arg, const std::vector< SX > &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXFunction::SXFunction(const SX &arg, const std::vector< std::vector< SXElement > > &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXFunction::SXFunction(const std::vector< SX > &arg, const SX &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const SX &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXFunction::SXFunction(const std::vector< std::vector< SXElement > > &arg, const std::vector< std::vector< SXElement > > &res) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SXNode::SXNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::ScalarSparseSparsity::ScalarSparseSparsity() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::ScalarSparsity::ScalarSparsity() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SetNonzeros< Add >::SetNonzeros(const MX &y, const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SetNonzerosSlice2< Add >::SetNonzerosSlice2(const MX &y, const MX &x, const Slice &inner, const Slice &outer) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SetNonzerosSlice< Add >::SetNonzerosSlice(const MX &y, const MX &x, const Slice &s) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SetNonzerosVector< Add >::SetNonzerosVector(const MX &y, const MX &x, const std::vector< int > &nz) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SetSparse::SetSparse(const MX &x, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SharedObject::SharedObject() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SharedObject::SharedObject(const SharedObject &ref) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SharedObjectNode::SharedObjectNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SharedObjectNode::SharedObjectNode(const SharedObjectNode &node) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Solve< Tr >::Solve(const MX &r, const MX &A, const LinearSolver &linear_solver) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Split::Split(const MX &x, const std::vector< int > &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::StabilizedQPSolverInputIOSchemeVector< M >::StabilizedQPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SubAssign::SubAssign(const MX &x, const MX &y, const Slice &i, const Slice &j) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SubRef::SubRef(const MX &x, const Slice &i, const Slice &j) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SymbolicMX::SymbolicMX(const std::string &name, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SymbolicMX::SymbolicMX(const std::string &name, int nrow=1, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::SymbolicSX::SymbolicSX(const std::string &name) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Transpose::Transpose(const MX &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::UnaryMX::UnaryMX(Operation op, MX x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Vertcat::Vertcat(const std::vector< MX > &x) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::Vertsplit::Vertsplit(const MX &x, const std::vector< int > &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::WeakRef::WeakRef(SharedObject shared) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::WeakRef::WeakRef(int dummy=0) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::XMLNode::XMLNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception CasADi::ZeroSX::ZeroSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}