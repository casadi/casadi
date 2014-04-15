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
%exception  casadi::LinearSolver::spSolve(DMatrix &X, const DMatrix &B, bool transpose=false) const  {
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
%exception  casadi::Slice::print(std::ostream &stream=std::cout) const  {
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
%exception  casadi::SymbolicSX::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicSX::getOp() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::SymbolicSX::isSymbolic() const  {
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
%exception  casadi::XMLNode::checkName(const std::string &str) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::dump(std::ostream &stream, int indent=0) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::getAttribute(const std::string &attribute_name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::getName() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::getText() const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::getText(T &val) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::hasAttribute(const std::string &attribute_name) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::hasChild(const std::string &childname) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::readAttribute(const std::string &attribute_name, T &val, bool assert_existance=true) const  {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::setAttribute(const std::string &attribute_name, const std::string &attribute) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::setName(const std::string &name) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::XMLNode::size() const  {
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
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &dwork, std::vector< double > &eig_real, std::vector< double > &eig_imag) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception  casadi::slicot_periodic_schur(int n, int K, const std::vector< double > &a, std::vector< double > &t, std::vector< double > &z, std::vector< double > &eig_real, std::vector< double > &eig_imag) {
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
%exception  sqicDestroy() {
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
%exception casadi::EmptySparsity::EmptySparsity() {
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
%exception casadi::JacGInputIOSchemeVector< M >::JacGInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::JacGOutputIOSchemeVector< M >::JacGOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LPSolverInputIOSchemeVector< M >::LPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LPSolverOutputIOSchemeVector< M >::LPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::LPStructIOSchemeVector< T >::LPStructIOSchemeVector(const std::vector< M > &t) {
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
%exception casadi::MayerInputIOSchemeVector< M >::MayerInputIOSchemeVector(const std::vector< M > &t) {
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
%exception casadi::NLPSolverInputIOSchemeVector< M >::NLPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NLPSolverOutputIOSchemeVector< M >::NLPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::NanSX::NanSX() {
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
%exception casadi::OCPInputIOSchemeVector< M >::OCPInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OCPOutputIOSchemeVector< M >::OCPOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OneSX::OneSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OptionsFunctionalityNode::OptionsFunctionalityNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::OutputNode::OutputNode(const MX &parent, int oind) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QCQPSolverInputIOSchemeVector< M >::QCQPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QCQPSolverOutputIOSchemeVector< M >::QCQPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QCQPStructIOSchemeVector< T >::QCQPStructIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QPSolverInputIOSchemeVector< M >::QPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QPSolverOutputIOSchemeVector< M >::QPSolverOutputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::QPStructIOSchemeVector< T >::QPStructIOSchemeVector(const std::vector< M > &t) {
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
%exception casadi::Solve< Tr >::Solve(const MX &r, const MX &A, const LinearSolver &linear_solver) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::Split::Split(const MX &x, const std::vector< int > &offset) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::StabilizedQPSolverInputIOSchemeVector< M >::StabilizedQPSolverInputIOSchemeVector(const std::vector< M > &t) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SubAssign::SubAssign(const MX &x, const MX &y, const Slice &i, const Slice &j) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SubRef::SubRef(const MX &x, const Slice &i, const Slice &j) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, const Sparsity &sp) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SymbolicMX::SymbolicMX(const std::string &name, int nrow=1, int ncol=1) {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::SymbolicSX::SymbolicSX(const std::string &name) {
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
%exception casadi::XMLNode::XMLNode() {
  START INTERNAL_MSG() $action STOP { $action } 
}
%exception casadi::ZeroSX::ZeroSX() {
  START INTERNAL_MSG() $action STOP { $action } 
}