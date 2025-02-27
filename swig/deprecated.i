%exception  casadi::DaeBuilder::add_c(const std::string &name, const MX &new_cdef) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_d(const std::string &name, const MX &new_ddef) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_init(const MX &lhs, const MX &rhs) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_p(const std::string &name=std::string()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_q(const std::string &name=std::string()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_t(const std::string &name="t") {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_u(const std::string &name=std::string()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_variable(const std::string &name, casadi_int n=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_variable(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_variable_new(const MX &new_v) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_variable_new(const std::string &name, casadi_int n=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_variable_new(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_w(const std::string &name, const MX &new_wdef) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_x(const std::string &name=std::string()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_y(const std::string &name, const MX &new_ydef) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::add_z(const std::string &name=std::string()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::clear_all(const std::string &v) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::e() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::eliminate_d() {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::eliminate_quad() {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::eliminate_w() {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::find(const std::string &name) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::find(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::has_variable(const std::string &name) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::name(const std::vector< size_t > &ind) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::name(size_t ind) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::ne() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::new_variable(const std::string &name, casadi_int numel=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::set_alg(const std::string &name, const MX &alg_rhs) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::set_all(const std::string &v, const std::vector< std::string > &name) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::set_beq(const std::string &name, const MX &val) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::set_ode(const std::string &name, const MX &ode_rhs) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::sort_d() {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::sort_w() {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::sort_z(const std::vector< std::string > &z_order) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::t() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::var(const std::vector< size_t > &ind) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::var(size_t ind) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::variable(const std::string &name) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::variable(const std::string &name) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::variable(size_t ind) const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::DaeBuilder::variable(size_t ind) {
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
%exception  casadi::integrator(const std::string &name, const std::string &solver, const SXDict &dae, const Dict &opts=Dict()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}