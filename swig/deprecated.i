%exception  casadi::Callback::get_forward(const std::string &name, int nfwd, Dict &opts) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Callback::get_reverse(const std::string &name, int nadj, Dict &opts) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::derivative(int nfwd, int nadj) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::forward(const std::vector< DM > &arg, const std::vector< DM > &res, const std::vector< std::vector< DM > > &fseed, std::vector< std::vector< DM > > &output_fsens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::forward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &output_fsens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::forward(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &output_fsens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::forward(int nfwd) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::integrator_dae() {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::kernel_sum(const std::string &name, const std::pair< int, int > &size, double r, int n, const Dict &opts=Dict()) const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::map(const std::map< std::string, MX > &arg, const std::string &parallelization="serial") {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::map(const std::string &name, const std::string &parallelization, int n, const Dict &opts=Dict()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::map(const std::vector< MX > &arg, const std::string &parallelization="serial") {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::nl_var(const std::string &s_in, const std::vector< std::string > &s_out) const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::printDimensions(std::ostream &stream=casadi::userOut()) const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::printOption(const std::string &name, std::ostream &stream=casadi::userOut()) const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::printOptions(std::ostream &stream=casadi::userOut()) const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::reverse(const std::vector< DM > &arg, const std::vector< DM > &res, const std::vector< std::vector< DM > > &aseed, std::vector< std::vector< DM > > &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::reverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::reverse(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &output_asens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::reverse(int nfwd) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::collocationInterpolators(const std::vector< double > &tau_root, std::vector< std::vector< double > > &output_C, std::vector< double > &output_D) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::doc_linsol(const std::string &name) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::has_linsol(const std::string &name) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::linsol_new(const std::string &name, const std::string &solver, const Sparsity &sp, int nrhs, const Dict &opts=Dict()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::load_linsol(const std::string &name) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  nl_var(const MatType &expr, const MatType &var) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  norm_F(const MatType &x) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}