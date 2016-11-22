%exception  casadi::Function::gradient(const std::string &iind, const std::string &oind) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::gradient(const std::string &iind, int oind=0) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::gradient(int iind, const std::string &oind) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::gradient(int iind=0, int oind=0) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::setJacobian(const Function &jac, int iind=0, int oind=0, bool compact=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::set_jac_sparsity(const Sparsity &sp, const std::string &iind, const std::string &oind, bool compact=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::set_jac_sparsity(const Sparsity &sp, const std::string &iind, int oind, bool compact=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::set_jac_sparsity(const Sparsity &sp, int iind, const std::string &oind, bool compact=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::set_jac_sparsity(const Sparsity &sp, int iind, int oind, bool compact=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::tangent(const std::string &iind, const std::string &oind) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::tangent(const std::string &iind, int oind=0) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::tangent(int iind, const std::string &oind) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::tangent(int iind=0, int oind=0) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}