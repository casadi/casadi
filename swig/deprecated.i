%exception  casadi::DaeBuilder::alg() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::derivatives() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::ode() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::outputs() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::quad() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::set_y(const std::vector< std::string > &name) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::ydef() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::zero() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::hessian_old(casadi_int iind, casadi_int oind) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::jacobian_old(casadi_int iind, casadi_int oind) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::sparsity_jac(casadi_int iind, casadi_int oind, bool compact=false, bool symmetric=false) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::sparsity_jac(casadi_int iind, const std::string &oind, bool compact=false, bool symmetric=false) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::sparsity_jac(const std::string &iind, casadi_int oind=0, bool compact=false, bool symmetric=false) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::sparsity_jac(const std::string &iind, const std::string &oind, bool compact=false, bool symmetric=false) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< MatType >::linearize(const MatType &f, const MatType &x, const MatType &x0, const Dict &opts=Dict()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const SXDict &dae, const Dict &opts=Dict()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}