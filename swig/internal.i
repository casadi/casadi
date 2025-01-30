%exception  casadi::AlpaqaProblem::eval_f(crvec x) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_f_grad_f(crvec x, rvec grad_fx) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_g(crvec x, rvec g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_grad_f(crvec x, rvec grad_fx) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_grad_gi(crvec x, index_t i, rvec grad_i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_grad_ps(crvec x, crvec y, crvec S, rvec grad_ps, rvec work_n, rvec work_m) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_hess_L(crvec x, crvec y, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_hess_ps(crvec x, crvec y, crvec S, real_t scale, rindexvec inner_idx, rindexvec outer_ptr, rvec H_values) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_hess_ps_prod(crvec x, crvec y, crvec S, real_t scale, crvec v, rvec Hv) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_jac_g(crvec x, rindexvec inner_idx, rindexvec outer_ptr, rvec J_values) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_ps(crvec x, crvec y, crvec S, rvec y) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::eval_ps_grad_ps(crvec x, crvec y, crvec S, rvec grad_ps, rvec work_n, rvec work_m) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::get_hess_L_num_nonzeros() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::get_hess_ps_num_nonzeros() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::get_jac_g_num_nonzeros() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::provides_eval_grad_gi() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::provides_eval_hess_L() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::provides_eval_hess_L_prod() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::provides_eval_hess_ps() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::provides_eval_hess_ps_prod() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::AlpaqaProblem::provides_eval_jac_g() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::construct(const std::string &name, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::eval(const std::vector< DM > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::finalize() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_name_in(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_name_out(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_sparsity_in(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::get_sparsity_out(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::has_eval_buffer() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::has_forward(casadi_int nfwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::has_jac_sparsity(casadi_int oind, casadi_int iind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::has_jacobian() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::has_reverse(casadi_int nadj) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Callback::uses_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CasadiException::what() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add(const Function &f, bool with_jac_sparsity=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add_auxiliary(Auxiliary f, const std::vector< std::string > &inst={"casadi_real"}) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add_dependency(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add_external(const std::string &new_external) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add_include(const std::string &new_include, bool relative_path=false, const std::string &use_ifdef=std::string()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add_io_sparsities(const std::string &name, const std::vector< Sparsity > &sp_in, const std::vector< Sparsity > &sp_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add_sparsity(const Sparsity &sp, bool canonical=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::arg(casadi_int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::avoid_stack() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::axpy(casadi_int n, const std::string &a, const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::bilin(const std::string &A, const Sparsity &sp_A, const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::bound_consistency(casadi_int n, const std::string &x, const std::string &lam, const std::string &lbx, const std::string &ubx) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::cache_check(const std::string &key, const std::string &cache, const std::string &loc, casadi_int stride, casadi_int sz, casadi_int key_sz, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::clear(const std::string &res, std::size_t n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::clip_max(const std::string &x, casadi_int n, const std::string &min, const std::string &mask) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::clip_min(const std::string &x, casadi_int n, const std::string &min, const std::string &mask) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::comment(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(casadi_int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(char v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::string &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::vector< casadi_int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::vector< char > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::vector< double > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::vector< int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::vector< std::string > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(double v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant_copy(const std::string &var_name, const std::vector< casadi_int > &v, const std::string &type="casadi_int") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::convexify_eval(const ConvexifyData &d, const std::string &Hin, const std::string &Hout, const std::string &iw, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::copy(const std::string &arg, std::size_t n, const std::string &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::copy_check(const std::string &arg, std::size_t n, const std::string &res, bool check_lhs=true, bool check_rhs=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::copy_default(const std::string &arg, std::size_t n, const std::string &res, const std::string &def, bool check_rhs=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::declare(std::string s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::define_pool_double(const std::string &name, const std::vector< double > &def) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::densify(const std::string &arg, const Sparsity &sp_arg, const std::string &res, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::dot(casadi_int n, const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::dump() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::dump(std::ostream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::elide_copy(casadi_int sz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::file_slurp(const std::string &fname, casadi_int n, const std::string &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::fill(const std::string &res, std::size_t n, const std::string &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::flush(std::ostream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::fmax(const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::fmin(const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::format_padded(casadi_int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::from_mex(std::string &arg, const std::string &res, std::size_t res_off, const Sparsity &sp_res, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::generate(const std::string &prefix="") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::get_constant(const std::vector< casadi_int > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::get_constant(const std::vector< char > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::get_constant(const std::vector< double > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::get_constant(const std::vector< std::string > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::get_sparsity(const Sparsity &sp) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::indent() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::init_local(const std::string &name, const std::string &def) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::initializer(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::interpn(const std::string &res, casadi_int ndim, const std::string &grid, const std::string &offset, const std::string &values, const std::string &x, const std::string &lookup_mode, casadi_int m, const std::string &iw, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::interpn_grad(const std::string &grad, casadi_int ndim, const std::string &grid, const std::string &offset, const std::string &values, const std::string &x, const std::string &lookup_mode, casadi_int m, const std::string &iw, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::lb_eig(const Sparsity &sp_h, const std::string &h) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::ldl(const std::string &sp_a, const std::string &a, const std::string &sp_lt, const std::string &lt, const std::string &d, const std::string &p, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::ldl_solve(const std::string &x, casadi_int nrhs, const std::string &sp_lt, const std::string &lt, const std::string &d, const std::string &p, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::local(const std::string &name, const std::string &type, const std::string &ref="") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::logsumexp(const std::string &A, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::low(const std::string &x, const std::string &grid, casadi_int ng, casadi_int lookup_mode) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::lsqr_solve(const std::string &A, const std::string &x, casadi_int nrhs, bool tr, const std::string &sp, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::masked_norm_inf(casadi_int n, const std::string &x, const std::string &mask) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::max(const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::max_viol(casadi_int n, const std::string &x, const std::string &lb, const std::string &ub) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::mem(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::min(const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::mmax(const std::string &x, casadi_int n, bool is_dense) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::mmin(const std::string &x, casadi_int n, bool is_dense) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::mtimes(const std::string &x, const Sparsity &sp_x, const std::string &y, const Sparsity &sp_y, const std::string &z, const Sparsity &sp_z, const std::string &w, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::mv(const std::string &x, casadi_int nrow_x, casadi_int ncol_x, const std::string &y, const std::string &z, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::mv(const std::string &x, const Sparsity &sp_x, const std::string &y, const std::string &z, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::norm_1(casadi_int n, const std::string &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::norm_2(casadi_int n, const std::string &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::norm_inf(casadi_int n, const std::string &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::ones(casadi_int sz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::pool_double(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::print_formatted(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::print_op(casadi_int op, const std::string &a0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::print_op(casadi_int op, const std::string &a0, const std::string &a1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::print_vector(std::ostream &s, const std::string &name, const std::vector< casadi_int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::print_vector(std::ostream &s, const std::string &name, const std::vector< char > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::print_vector(std::ostream &s, const std::string &name, const std::vector< double > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::print_vector(std::ostream &s, const std::string &name, const std::vector< std::string > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::printf(const std::string &str, const std::string &arg1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::printf(const std::string &str, const std::string &arg1, const std::string &arg2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::printf(const std::string &str, const std::string &arg1, const std::string &arg2, const std::string &arg3) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::printf(const std::string &str, const std::vector< std::string > &arg=std::vector< std::string >()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::project(const std::string &arg, const Sparsity &sp_arg, const std::string &res, const Sparsity &sp_res, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::qr(const std::string &sp, const std::string &A, const std::string &w, const std::string &sp_v, const std::string &v, const std::string &sp_r, const std::string &r, const std::string &beta, const std::string &prinv, const std::string &pc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::qr_solve(const std::string &x, casadi_int nrhs, bool tr, const std::string &sp_v, const std::string &v, const std::string &sp_r, const std::string &r, const std::string &beta, const std::string &prinv, const std::string &pc, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::rank1(const std::string &A, const Sparsity &sp_A, const std::string &alpha, const std::string &x, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::regularize(const Sparsity &sp_h, const std::string &h, const std::string &reg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::res(casadi_int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::reserve_work(casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sanitize_source(const std::string &src, const std::vector< std::string > &inst, bool add_shorthand=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::scal(casadi_int n, const std::string &alpha, const std::string &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::scope_enter() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::scope_exit() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::setup_callback(const std::string &s, const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::shorthand(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::shorthand(const std::string &name, bool allow_adding=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sparsify(const std::string &arg, const std::string &res, const Sparsity &sp_res, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sparsity(const Sparsity &sp, bool canonical=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sum_viol(casadi_int n, const std::string &x, const std::string &lb, const std::string &ub) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sx_work(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sz_work(size_t &sz_arg, size_t &sz_res, size_t &sz_iw, size_t &sz_w) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::to_mex(const Sparsity &sp, const std::string &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::trans(const std::string &x, const Sparsity &sp_x, const std::string &y, const Sparsity &sp_y, const std::string &iw) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::tri_project(const std::string &arg, const Sparsity &sp_arg, const std::string &res, bool lower) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::trilsolve(const Sparsity &sp_x, const std::string &x, const std::string &y, bool tr, bool unity, casadi_int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::triusolve(const Sparsity &sp_x, const std::string &x, const std::string &y, bool tr, bool unity, casadi_int nrhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::unindent() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::vector_fmax(casadi_int n, const std::string &x, const std::string &y, const std::string &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::vector_fmin(casadi_int n, const std::string &x, const std::string &y, const std::string &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::vfmax(const std::string &x, casadi_int n, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::vfmax(const std::string &x, const std::string &n, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::vfmin(const std::string &x, casadi_int n, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::vfmin(const std::string &x, const std::string &n, const std::string &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::work(casadi_int n, casadi_int sz, bool is_ref) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::workel(casadi_int n) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::zeros(casadi_int sz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Contraction(const T &a, const T &b, T &r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Contraction(const bvec_t &a, const bvec_t &b, bvec_t &r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_c(const std::string &name, const MX &new_cdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_d(const std::string &name, const MX &new_ddef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_e(const std::string &name, const MX &new_edef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_fun(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_fun(const std::string &name, const Importer &compiler, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_fun(const std::string &name, const std::vector< std::string > &arg, const std::vector< std::string > &res, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_lc(const std::string &name, const std::vector< std::string > &f_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_p(const std::string &name=std::string()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_q(const std::string &name=std::string()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_t(const std::string &name="t") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_u(const std::string &name=std::string()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_variable(const MX &new_v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_variable(const std::string &name, casadi_int n=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_variable(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_variable_new(const MX &new_v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_variable_new(const std::string &name, casadi_int n=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_variable_new(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_w(const std::string &name, const MX &new_wdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_when(const MX &cond, const MX &lhs, const MX &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_x(const std::string &name=std::string()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_y(const std::string &name, const MX &new_ydef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_z(const std::string &name=std::string()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::alg() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::all_variables() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::attribute(const std::string &a, const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::attribute(const std::string &a, const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::beq(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::c() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::causality(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::cdef() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::clear_all(const std::string &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::clear_when() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::create(const std::string &fname, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, bool sx, bool lifted_calls=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::create(const std::string &name, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::create(const std::string &name, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::d() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::ddef() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::dependent_fun(const std::string &fname, const std::vector< std::string > &s_in, const std::vector< std::string > &s_out) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::der(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::der(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::derivatives() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::description(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::dimension(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::display_unit(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::e() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::eliminate_d() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::eliminate_quad() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::eliminate_w() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::event_transition(const std::string &fname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::event_transition(const std::string &fname, casadi_int index) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::export_fmu(const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::find(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::find(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::fun() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::fun(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::gather_fun(casadi_int max_depth=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::get(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::get(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::has_fun(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::has_t() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::has_variable(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::init_lhs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::init_rhs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::initial(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::initial_unknowns() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::jac_sparsity(const std::vector< std::string > &onames, const std::vector< std::string > &inames) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::lift(bool lift_shared=true, bool lift_calls=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::load_fmi_description(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::max(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::max(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::min(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::min(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::name(const std::vector< size_t > &ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::name(size_t ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nc() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nd() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::ne() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::new_variable(const std::string &name, casadi_int numel=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nominal(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nominal(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::np() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nq() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nu() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::numel(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nw() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::ny() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::nz() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::ode() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::oracle(bool sx=false, bool elim_w=false, bool lifted_calls=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::outputs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::p() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::parse_fmi(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::provides_directional_derivative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::provides_directional_derivatives() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::prune(bool prune_p=true, bool prune_u=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::q() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::quad() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_c(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_d(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_e(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_p(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_q(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_t(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_u(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_w(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_x(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_y(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::register_z(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::sanity_check() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set(const std::string &name, double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set(const std::vector< std::string > &name, const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set(const std::vector< std::string > &name, const std::vector< std::string > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_alg(const std::string &name, const MX &alg_rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_all(const std::string &v, const std::vector< std::string > &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_attribute(const std::string &a, const std::string &name, double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_attribute(const std::string &a, const std::vector< std::string > &name, const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_beq(const std::string &name, const MX &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_causality(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_description(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_display_unit(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_init(const std::string &name, const MX &init_rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_initial(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_max(const std::string &name, double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_max(const std::vector< std::string > &name, const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_min(const std::string &name, double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_min(const std::vector< std::string > &name, const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_nominal(const std::string &name, double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_nominal(const std::vector< std::string > &name, const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_ode(const std::string &name, const MX &ode_rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_start(const std::string &name, double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_start(const std::vector< std::string > &name, const std::vector< double > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_start_time(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_step_size(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_stop_time(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_tolerance(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_type(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_unit(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_value_reference(const std::string &name, casadi_int val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::set_variability(const std::string &name, const std::string &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::sort_d() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::sort_w() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::sort_z(const std::vector< std::string > &z_order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::start(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::start(const std::vector< std::string > &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::start_time() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::step_size() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::stop_time() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::t() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::tear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::tolerance() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::type(const std::string &name, casadi_int fmi_version=3) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::type_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::u() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::unit(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::value_reference(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::var(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::var(const std::vector< size_t > &ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::var(size_t ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::variability(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::variable(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::variable(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::variable(size_t ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::variable(size_t ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::w() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::wdef() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::when_cond() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::when_lhs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::when_rhs() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::x() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::y() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::ydef() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::z() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::zero() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_dm() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_dm_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_double() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_double_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_function() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_function_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_generictype() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_generictype_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_int() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_int_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_linsol() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_linsol_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_mx() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_mx_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_sparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_sparsity_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_string() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_string_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_sx() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::blind_unpack_sx_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::connect(SerializerBase &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::pop_type() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_dm() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_dm_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_double() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_double_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_function() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_function_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_generictype() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_generictype_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_int() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_int_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_linsol() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_linsol_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_mx() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_mx_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_sparsity() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_sparsity_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_string() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_string_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_sx() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializerBase::unpack_sx_vector() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::connect(SerializingStream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(Fmu &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(Function &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(GenericType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(Importer &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(Linsol &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(MX &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(Matrix< T > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(SXElem &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(Slice &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(Sparsity &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(bool &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(casadi_int &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(char &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(const std::string &descr, T &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(double &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(int &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(size_t &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(std::map< K, V > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(std::ostream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(std::pair< A, B > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(std::string &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::unpack(std::vector< T > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::version(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::version(const std::string &name, int min, int max) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DeserializingStream::version(const std::string &name, int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::assert_size_in(casadi_int i, casadi_int nrow, casadi_int ncol) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::assert_size_out(casadi_int i, casadi_int nrow, casadi_int ncol) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::assert_sparsity_out(casadi_int i, const Sparsity &sp, casadi_int n=1, bool allow_all_zero_sparse=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::buf_in(MapArg arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::buf_in(VecArg arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::buf_out(MPrRes res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::buf_out(MapRes res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::buf_out(VPrRes res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::buf_out(VecRes res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::cache() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::call(const DMDict &arg, DMDict &res, bool always_inline=false, bool never_inline=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::call(const MXDict &arg, MXDict &res, bool always_inline=false, bool never_inline=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::call(const SXDict &arg, SXDict &res, bool always_inline=false, bool never_inline=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::call(const std::vector< DM > &arg, std::vector< DM > &res, bool always_inline=false, bool never_inline=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::call(const std::vector< MX > &arg, std::vector< MX > &res, bool always_inline=false, bool never_inline=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::call(const std::vector< SX > &arg, std::vector< SX > &res, bool always_inline=false, bool never_inline=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::call_gen(const std::map< std::string, M > &arg, std::map< std::string, M > &res, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::change_option(const std::string &option_name, const GenericType &option_value) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::checkout() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_in(const DMDict &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_in(const MXDict &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_in(const SXDict &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_in(const std::vector< DM > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_in(const std::vector< MX > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_in(const std::vector< SX > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_out(const DMDict &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_out(const MXDict &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_out(const SXDict &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_out(const std::vector< DM > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_out(const std::vector< MX > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::convert_out(const std::vector< SX > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::default_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::expand() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::expand(const std::string &name, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::export_code(const std::string &lang, const Dict &options=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::export_code(const std::string &lang, const std::string &fname, const Dict &options=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::export_code(const std::string &lang, std::ostream &stream, const Dict &options=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::factory(const std::string &name, const std::vector< std::string > &s_in, const std::vector< std::string > &s_out, const AuxOut &aux=AuxOut(), const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::find_function(const std::string &name, casadi_int max_depth=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::find_functions(casadi_int max_depth=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::fold(casadi_int N, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::forward(casadi_int nfwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::free_mx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::free_sx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate(const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate(const std::string &fname, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate_dependencies(const std::string &fname, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate_in(const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate_in(const std::string &fname, const std::vector< DM > &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate_lifted(Function &vdef_fcn, Function &vinit_fcn) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate_out(const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::generate_out(const std::string &fname, const std::vector< DM > &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::get_free() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::get_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::get_function(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::has_free() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::has_function(const std::string &fname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::has_option(const std::string &option_name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::has_spfwd() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::has_sprev() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::index_in(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::index_out(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::info() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::instruction_MX(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::instruction_constant(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::instruction_id(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::instruction_input(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::instruction_output(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::instructions_sx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::is_a(const std::string &type, bool recursive=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::is_diff_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::is_diff_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::is_diff_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::is_diff_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::jac_sparsity(bool compact=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::jac_sparsity(casadi_int oind, casadi_int iind, bool compact=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::jacobian() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::jit(const std::string &name, const std::string &body, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::jit(const std::string &name, const std::string &body, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const std::vector< Sparsity > &sparsity_in, const std::vector< Sparsity > &sparsity_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::map(casadi_int n, const std::string &parallelization, casadi_int max_num_threads) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::map(casadi_int n, const std::string &parallelization="serial") const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::map(casadi_int n, const std::vector< bool > &reduce_in, const std::vector< bool > &reduce_out=std::vector< bool >(), const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::map(const std::string &name, const std::string &parallelization, casadi_int n, const std::vector< casadi_int > &reduce_in, const std::vector< casadi_int > &reduce_out, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::map(const std::string &name, const std::string &parallelization, casadi_int n, const std::vector< std::string > &reduce_in, const std::vector< std::string > &reduce_out, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mapaccum(casadi_int N, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mapaccum(const std::string &name, casadi_int N, casadi_int n_accum, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mapaccum(const std::string &name, casadi_int N, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mapaccum(const std::string &name, casadi_int n, const std::vector< casadi_int > &accum_in, const std::vector< casadi_int > &accum_out, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mapaccum(const std::string &name, casadi_int n, const std::vector< std::string > &accum_in, const std::vector< std::string > &accum_out, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mapsum(const std::vector< MX > &x, const std::string &parallelization="serial") const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::max_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::memory(int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::merge(const std::vector< MX > &arg, std::vector< MX > &subs_from, std::vector< MX > &subs_to) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::min_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mx_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mx_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mx_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mx_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mx_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::mx_out(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::n_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::n_instructions() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::n_nodes() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::n_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::name_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::name_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::name_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::name_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nnz_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nnz_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nnz_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nnz_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nnz_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nnz_out(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nominal_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nominal_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::numel_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::numel_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::numel_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::numel_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::numel_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::numel_out(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nz_from_in(const std::vector< DM > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nz_from_out(const std::vector< DM > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nz_to_in(const std::vector< double > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::nz_to_out(const std::vector< double > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::oracle() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::print_dimensions(std::ostream &stream=casadi::uout()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::print_option(const std::string &name, std::ostream &stream=casadi::uout()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::print_options(std::ostream &stream=casadi::uout()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::release(int mem) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::reverse(casadi_int nadj) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::save(const std::string &fname, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::serialize(const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::serialize(std::ostream &stream, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size1_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size1_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size1_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size1_out(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size2_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size2_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size2_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size2_out(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::size_out(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::slice(const std::string &name, const std::vector< casadi_int > &order_in, const std::vector< casadi_int > &order_out, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sparsity_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sparsity_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sparsity_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sparsity_out(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::stats(int mem=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sx_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sx_in(casadi_int iind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sx_in(const std::string &iname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sx_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sx_out(casadi_int oind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sx_out(const std::string &oname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_arg() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_iw() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_res() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_w() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::sz_work(size_t &sz_arg, size_t &sz_res, size_t &sz_iw, size_t &sz_w) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::uses_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::which_depends(const std::string &s_in, const std::vector< std::string > &s_out, casadi_int order=1, bool tr=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::wrap() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::wrap_as_needed(const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionBuffer::ret() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionBuffer::stats() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::abs(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::acos(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::acosh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::asin(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::asinh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::atan(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::atan2(const ExType &y, const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::atanh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::ceil(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::constpow(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::copysign(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::cos(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::cosh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::eq(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::erf(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::erfinv(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::exp(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::expm1(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::fabs(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::floor(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::fmax(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::fmin(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::fmod(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::ge(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::gt(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::hypot(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::if_else_zero(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::is_equal(const ExType &x, const ExType &y, casadi_int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::le(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::log(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::log10(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::log1p(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::logic_and(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::logic_not(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::logic_or(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::lt(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::minus(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::mod(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::ne(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::plus(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::pow(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::printme(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::rdivide(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::remainder(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::sign(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::sin(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::sinh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::sq(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::sqrt(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::tan(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::tanh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression::times(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::abs(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::acos(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::acosh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::asin(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::asinh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::atan(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::atan2(const ExType &y, const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::atanh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::ceil(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::constpow(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::copysign(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::cos(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::cosh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::eq(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::erf(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::erfinv(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::exp(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::expm1(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::floor(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::fmax(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::fmin(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::ge(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::gt(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::hypot(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::if_else_zero(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::le(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::log(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::log10(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::log1p(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::logic_and(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::logic_not(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::logic_or(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::lt(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::minus(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::mod(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::ne(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::pow(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::printme(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::rdivide(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::remainder(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::sign(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::sin(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::sinh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::sq(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::sqrt(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::tan(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::tanh(const ExType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< ExType >::times(const ExType &x, const ExType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< MX  >::printme(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< Matrix< Scalar >  >::printme(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::abs(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::acos(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::acosh(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::asin(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::asinh(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::atan(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::atan2(const SXElem &y, const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::atanh(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::ceil(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::constpow(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::copysign(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::cos(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::cosh(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::eq(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::erf(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::erfinv(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::exp(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::expm1(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::floor(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::fmax(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::fmin(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::ge(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::gt(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::hypot(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::if_else_zero(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::le(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::log(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::log10(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::log1p(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::logic_and(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::logic_not(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::logic_or(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::lt(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::minus(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::mod(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::ne(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::pow(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::printme(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::rdivide(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::remainder(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::sign(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::sin(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::sinh(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::sq(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::sqrt(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::tan(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::tanh(const SXElem &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExpression< SXElem  >::times(const SXElem &x, const SXElem &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::bilin(const MatType &A, const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::conditional(const MatType &ind, const std::vector< MatType > &x, const MatType &x_default, bool short_circuit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::cross(const MatType &a, const MatType &b, casadi_int dim=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::cse(const MatType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::cse(const std::vector< MatType > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::cumsum(const MatType &x, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::densify(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::densify(const MatType &x, const MatType &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::depends_on(const MatType &f, const MatType &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::det(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::diag(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::diff(const MatType &x, casadi_int n=1, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::dot(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::einstein(const MatType &A, const MatType &B, const MatType &C, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::einstein(const MatType &A, const MatType &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::expm(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::expm_const(const MatType &A, const MatType &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::extract(std::vector< MatType > &ex, std::vector< MatType > &v, std::vector< MatType > &vdef, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::extract_parametric(const MatType &expr, const MatType &par, MatType &expr_ret, std::vector< MatType > &symbols, std::vector< MatType > &parametric, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::extract_parametric(const MatType &expr, const std::vector< MatType > &par, MatType &expr_ret, std::vector< MatType > &symbols, std::vector< MatType > &parametric, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::extract_parametric(const std::vector< MatType > &expr, const MatType &par, std::vector< MatType > &expr_ret, std::vector< MatType > &symbols, std::vector< MatType > &parametric, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::extract_parametric(const std::vector< MatType > &expr, const std::vector< MatType > &par, std::vector< MatType > &expr_ret, std::vector< MatType > &symbols, std::vector< MatType > &parametric, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::forward(const std::vector< MatType > &ex, const std::vector< MatType > &arg, const std::vector< std::vector< MatType > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::hessian(const MatType &ex, const MatType &arg, MatType &output_g, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::hessian(const MatType &ex, const MatType &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::if_else(const MatType &cond, const MatType &if_true, const MatType &if_false, bool short_circuit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::interp1d(const std::vector< double > &x, const MatType &v, const std::vector< double > &xq, const std::string &mode, bool equidistant=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::inv(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::inv(const MatType &A, const std::string &lsolver, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::inv_minor(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::inv_skew(const MatType &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::is_linear(const MatType &expr, const MatType &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::is_quadratic(const MatType &expr, const MatType &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::jacobian(const MatType &ex, const MatType &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::jacobian_sparsity(const MatType &f, const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::linear_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, bool check=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::linspace(const MatType &a, const MatType &b, casadi_int nsteps) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::logsumexp(const MatType &x, const MatType &margin) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::mldivide(const MatType &x, const MatType &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::mmax(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::mmin(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::mrdivide(const MatType &x, const MatType &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::n_nodes(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::norm_1(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::norm_2(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::norm_fro(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::norm_inf(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::nullspace(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::pinv(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::pinv(const MatType &A, const std::string &lsolver, const Dict &dict=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::polyval(const MatType &p, const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::print_operator(const MatType &xb, const std::vector< std::string > &args) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::project(const MatType &A, const Sparsity &sp, bool intersect=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::quadratic_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, MatType &c, bool check=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::repsum(const MatType &A, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::reverse(const std::vector< MatType > &ex, const std::vector< MatType > &arg, const std::vector< std::vector< MatType > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::separate_linear(const MatType &expr, const MatType &sym_lin, const MatType &sym_const, MatType &expr_const, MatType &expr_lin, MatType &expr_nonlin) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::shared(std::vector< MatType > &ex, std::vector< MatType > &v, std::vector< MatType > &vdef, const std::string &v_prefix="v_", const std::string &v_suffix="") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::simplify(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::skew(const MatType &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::solve(const MatType &A, const MatType &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::solve(const MatType &A, const MatType &b, const std::string &lsolver, const Dict &dict=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::substitute(const MatType &ex, const MatType &v, const MatType &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::substitute(const std::vector< MatType > &ex, const std::vector< MatType > &v, const std::vector< MatType > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::substitute_inplace(const std::vector< MatType > &v, std::vector< MatType > &inout_vdef, std::vector< MatType > &inout_ex, bool reverse=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::sumsqr(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::symvar(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::trace(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::tril2symm(const MatType &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::triu2symm(const MatType &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::unite(const MatType &A, const MatType &B) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::which_depends(const MatType &expr, const MatType &var, casadi_int order, bool tr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::bilin(const MX &A, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::colind(casadi_int col) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::columns() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::cross(const MX &a, const MX &b, casadi_int dim=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::diff(const MX &x, casadi_int n=1, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::dim(bool with_nz=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::get_colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::get_row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::gradient(const MX &ex, const MX &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::interp1d(const std::vector< double > &x, const MX &v, const std::vector< double > &xq, const std::string &mode, bool equidistant) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::inv_skew(const MX &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_column() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_dense() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_empty(bool both=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_linear(const MX &expr, const MX &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_quadratic(const MX &expr, const MX &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_scalar(bool scalar_and_dense=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_square() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_tril() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_triu() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::jtimes(const MX &ex, const MX &arg, const MX &v, bool tr=false, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::linear_coeff(const MX &expr, const MX &var, MX &A, MX &b, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::linearize(const MX &f, const MX &x, const MX &x0, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::linspace(const MX &a, const MX &b, casadi_int nsteps) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::mpower(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nnz() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nnz_diag() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nnz_lower() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nnz_upper() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::norm_0_mul(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::numel() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nz(const K &k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nz(const K &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::ones(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::ones(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::ones(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::quadratic_coeff(const MX &expr, const MX &var, MX &A, MX &b, MX &c, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::rank1(const MX &A, const MX &alpha, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::repsum(const MX &x, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::row(casadi_int el) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::rows() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::size() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::size(casadi_int axis) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::size1() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::size2() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::skew(const MX &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::soc(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sprank(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sumsqr(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sym(const std::string &name, casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sym(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sym(const std::string &name, const Sparsity &sp, casadi_int p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sym(const std::string &name, const Sparsity &sp, casadi_int p, casadi_int r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sym(const std::string &name, const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::tangent(const MX &ex, const MX &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::tril(const MX &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::tril2symm(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::triu(const MX &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::triu2symm(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::zeros(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::zeros(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::zeros(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::bilin(const MatType &A, const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::colind(casadi_int col) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::columns() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::cross(const MatType &a, const MatType &b, casadi_int dim=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::diff(const MatType &x, casadi_int n=1, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::dim(bool with_nz=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::get_colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::get_row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::gradient(const MatType &ex, const MatType &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::interp1d(const std::vector< double > &x, const MatType &v, const std::vector< double > &xq, const std::string &mode, bool equidistant) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::inv_skew(const MatType &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_column() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_dense() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_empty(bool both=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_linear(const MatType &expr, const MatType &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_quadratic(const MatType &expr, const MatType &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_scalar(bool scalar_and_dense=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_square() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_tril() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_triu() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::is_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::jtimes(const MatType &ex, const MatType &arg, const MatType &v, bool tr=false, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::linear_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::linearize(const MatType &f, const MatType &x, const MatType &x0, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::linspace(const MatType &a, const MatType &b, casadi_int nsteps) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::logsumexp(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::mpower(const MatType &x, const MatType &n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::mpower(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nnz() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nnz_diag() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nnz_lower() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nnz_upper() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::norm_0_mul(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::numel() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nz(const K &k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nz(const K &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::ones(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::ones(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::ones(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::quadratic_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, MatType &c, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::rank1(const MatType &A, const MatType &alpha, const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::repsum(const MatType &x, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::row(casadi_int el) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::rows() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::size() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::size(casadi_int axis) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::size1() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::size2() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::skew(const MatType &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::soc(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sparsity() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sprank(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sumsqr(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sym(const std::string &name, casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sym(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sym(const std::string &name, const Sparsity &sp, casadi_int p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sym(const std::string &name, const Sparsity &sp, casadi_int p, casadi_int r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::sym(const std::string &name, const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::tangent(const MatType &ex, const MatType &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::tril(const MatType &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::tril2symm(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::triu(const MatType &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::triu2symm(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::zeros(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::zeros(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::zeros(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::bilin(const Matrix< Scalar > &A, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::colind(casadi_int col) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::columns() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::cross(const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int dim=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::diff(const Matrix< Scalar > &x, casadi_int n=1, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::dim(bool with_nz=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::get_colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::get_row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::gradient(const Matrix< Scalar > &ex, const Matrix< Scalar > &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::interp1d(const std::vector< double > &x, const Matrix< Scalar > &v, const std::vector< double > &xq, const std::string &mode, bool equidistant) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::inv_skew(const Matrix< Scalar > &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_column() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_dense() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_empty(bool both=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_linear(const Matrix< Scalar > &expr, const Matrix< Scalar > &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_quadratic(const Matrix< Scalar > &expr, const Matrix< Scalar > &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_scalar(bool scalar_and_dense=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_square() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_tril() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_triu() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::is_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::jtimes(const Matrix< Scalar > &ex, const Matrix< Scalar > &arg, const Matrix< Scalar > &v, bool tr=false, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::linear_coeff(const Matrix< Scalar > &expr, const Matrix< Scalar > &var, Matrix< Scalar > &A, Matrix< Scalar > &b, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::linearize(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &x0, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::linspace(const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int nsteps) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::mpower(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::nnz() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::nnz_diag() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::nnz_lower() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::nnz_upper() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::norm_0_mul(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::numel() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::nz(const K &k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::nz(const K &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::ones(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::ones(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::ones(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::quadratic_coeff(const Matrix< Scalar > &expr, const Matrix< Scalar > &var, Matrix< Scalar > &A, Matrix< Scalar > &b, Matrix< Scalar > &c, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::rank1(const Matrix< Scalar > &A, const Matrix< Scalar > &alpha, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::repsum(const Matrix< Scalar > &x, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::row(casadi_int el) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::rows() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::size() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::size(casadi_int axis) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::size1() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::size2() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::skew(const Matrix< Scalar > &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::soc(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sprank(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sumsqr(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, casadi_int nrow, casadi_int ncol, casadi_int p, casadi_int r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, const Sparsity &sp, casadi_int p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, const Sparsity &sp, casadi_int p, casadi_int r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::sym(const std::string &name, const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::tangent(const Matrix< Scalar > &ex, const Matrix< Scalar > &arg, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::tril(const Matrix< Scalar > &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::tril2symm(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::triu(const Matrix< Scalar > &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::triu2symm(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::zeros(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::zeros(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< Matrix< Scalar >  >::zeros(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericShared< Shared, Internal >::debug_repr() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericShared< Shared, Internal >::is_null() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericShared< SharedObject , SharedObjectInternal  >::debug_repr() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericShared< SharedObject , SharedObjectInternal  >::is_null() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_bool() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_bool_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_dict() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_dict_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_double() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_double_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_double_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_function_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_int() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_int_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_int_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_string() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_string_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_string_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_void_pointer() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::can_cast_to(TypeID other) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::can_cast_to(const GenericType &other) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::getType() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::get_description() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_bool() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_bool_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_dict() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_dict_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_double() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_double_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_double_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_empty_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_function_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_int() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_int_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_int_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_string() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_string_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_string_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::is_void_pointer() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_bool() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_bool_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_dict() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_dict_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_double() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_double_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_double_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_function_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_int() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_int_type_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_int_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_int_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_string() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_string_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_string_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_vector_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::to_void_pointer() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericWeakRef< Shared, Internal >::alive() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericWeakRef< Shared, Internal >::shared() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericWeakRef< Shared, Internal >::shared_if_alive(Shared &shared) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericWeakRef< SharedObject , SharedObjectInternal  >::alive() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericWeakRef< SharedObject , SharedObjectInternal  >::shared() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericWeakRef< SharedObject , SharedObjectInternal  >::shared_if_alive(SharedObject &shared) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::body(const std::string &symname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::get_function(const std::string &symname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::get_meta(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::has_function(const std::string &symname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::has_meta(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::inlined(const std::string &symname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::library() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::meta_int(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::meta_set(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::meta_string(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::meta_vector(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::to(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::checkout() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::neig(const DM &A) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::nfact(const DM &A) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::rank(const DM &A) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::release(int mem) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::sfact(const DM &A) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::solve(const DM &A, const DM &B, bool tr=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::solve(const MX &A, const MX &B, bool tr=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::sparsity() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::stats(int mem=1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::T() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::ad_forward(const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::ad_reverse(const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::attachAssert(const MX &y, const std::string &fail_message="") const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::binary(casadi_int op, const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::blockcat(const std::vector< std::vector< MX > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::bspline(const MX &x, const DM &coeffs, const std::vector< std::vector< double > > &knots, const std::vector< casadi_int > &degree, casadi_int m, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::bspline(const MX &x, const MX &coeffs, const std::vector< std::vector< double > > &knots, const std::vector< casadi_int > &degree, casadi_int m, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::bspline_dual(const std::vector< double > &x, const std::vector< std::vector< double > > &knots, const std::vector< casadi_int > &degree, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::conditional(const MX &ind, const std::vector< MX > &x, const MX &x_default, bool short_circuit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::convexify(const MX &H, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::cse(const std::vector< MX > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::cumsum(const MX &x, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::densify(const MX &x, const MX &val=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::dep(casadi_int ch=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::depends_on(const MX &x, const MX &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::det(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::diag(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::diagcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::diagsplit(const MX &x, const std::vector< casadi_int > &offset1, const std::vector< casadi_int > &offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::difference(const std::vector< MX > &a, const std::vector< MX > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::dot(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::einstein(const MX &A, const MX &B, const MX &C, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::einstein(const MX &A, const MX &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::enlarge(casadi_int nrow, casadi_int ncol, const std::vector< casadi_int > &rr, const std::vector< casadi_int > &cc, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::erase(const std::vector< casadi_int > &rr, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::erase(const std::vector< casadi_int > &rr, const std::vector< casadi_int > &cc, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::eval_mx(const std::vector< MX > &arg, std::vector< MX > &res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::evalf(const MX &expr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::evalf(const MX &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::expm(const MX &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::expm_const(const MX &A, const MX &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::extract(std::vector< MX > &ex, std::vector< MX > &v, std::vector< MX > &vdef, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::extract_parametric(const MX &expr, const MX &par, MX &expr_ret, std::vector< MX > &symbols, std::vector< MX > &parametric, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::find(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::forward(const std::vector< MX > &ex, const std::vector< MX > &arg, const std::vector< std::vector< MX > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, casadi_int rr, casadi_int cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, casadi_int rr, const Slice &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const MX &rr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const MX &rr, const MX &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const MX &rr, const Slice &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Matrix< casadi_int > &rr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Matrix< casadi_int > &rr, const Slice &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Slice &rr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Slice &rr, casadi_int cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Slice &rr, const MX &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Slice &rr, const Matrix< casadi_int > &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Slice &rr, const Slice &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const Sparsity &sp) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get(MX &m, bool ind1, const casadi_int rr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nonzeros() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nz(MX &m, bool ind1, casadi_int kk) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nz(MX &m, bool ind1, const MX &inner, const MX &outer) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nz(MX &m, bool ind1, const MX &inner, const Slice &outer) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nz(MX &m, bool ind1, const MX &kk) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nz(MX &m, bool ind1, const Matrix< casadi_int > &kk) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nz(MX &m, bool ind1, const Slice &inner, const MX &outer) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_nz(MX &m, bool ind1, const Slice &kk) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_output(casadi_int oind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_sparsity() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const MX &ex, const std::vector< MX > &v, const std::vector< MX > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const MX &ex, const std::vector< MX > &v, const std::vector< MX > &vdef, bool &updated) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const MX &x, const std::vector< MX > &v, const std::vector< MX > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const MX &x, const std::vector< MX > &v, const std::vector< MX > &vdef, bool &updated) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef, bool &updated) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::has_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::hessian(const MX &f, const MX &x, MX &g, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::hessian(const MX &f, const MX &x, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::horzcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::horzsplit(const MX &x, const std::vector< casadi_int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::if_else(const MX &cond, const MX &if_true, const MX &if_false, bool short_circuit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inf(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inf(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inf(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::info() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inv(const MX &A, const std::string &lsolver="qr", const Dict &dict=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inv_minor(const MX &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inv_node(const MX &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inv_node(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_binary() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_call() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_commutative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_constant() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_equal(const MX &x, const MX &y, casadi_int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_eye() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_minus_one() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_multiplication() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_norm() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_one() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_op(casadi_int op) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_regular() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_symbolic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_transpose() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_unary() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_valid_input() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_zero() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::jacobian(const MX &f, const MX &x, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::jacobian_sparsity(const MX &f, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::join_primitives(const std::vector< DM > &v) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::join_primitives(const std::vector< MX > &v) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::join_primitives(const std::vector< SX > &v) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::kron(const MX &x, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::lift(const MX &x, const MX &x_guess) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::low(const MX &v, const MX &p, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mac(const MX &x, const MX &y, const MX &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mapping() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const MX &e, const std::vector< MX > &boundary, const Dict &options) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const MX &e, const std::vector< MX > &boundary=std::vector< MX >(), const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const std::vector< MX > &e, const std::vector< MX > &boundary, const Dict &options) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const std::vector< MX > &e, const std::vector< MX > &boundary=std::vector< MX >(), const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mldivide(const MX &a, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mmax(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mmin(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::monitor(const std::string &comment) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mrdivide(const MX &a, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mtimes(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::n_dep() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::n_nodes(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::n_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::n_primitives() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::nan(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::nan(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::nan(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::no_grad(const MX &expr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::no_hess(const MX &expr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::norm_1(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::norm_2(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::norm_fro(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::norm_inf(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::nullspace(const MX &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::op() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::pinv(const MX &A, const std::string &lsolver="qr", const Dict &dict=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::polyval(const MX &p, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::primitives() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::print_operator(const MX &x, const std::vector< std::string > &args) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::printme(const MX &b) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::project(const MX &x, const Sparsity &sp, bool intersect=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::repmat(const MX &x, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::repsum(const MX &x, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::reshape(const MX &x, casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::reshape(const MX &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::reverse(const std::vector< MX > &ex, const std::vector< MX > &arg, const std::vector< std::vector< MX > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::separate_linear(const MX &expr, const MX &sym_lin, const MX &sym_const, MX &expr_const, MX &expr_lin, MX &expr_nonlin) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set(const MX &m, bool ind1, const Matrix< casadi_int > &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set(const MX &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set(const MX &m, bool ind1, const Matrix< casadi_int > &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set(const MX &m, bool ind1, const Slice &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set(const MX &m, bool ind1, const Slice &rr, const Matrix< casadi_int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set(const MX &m, bool ind1, const Slice &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set(const MX &m, bool ind1, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set_nz(const MX &m, bool ind1, casadi_int kk) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set_nz(const MX &m, bool ind1, const MX &kk) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set_nz(const MX &m, bool ind1, const Matrix< casadi_int > &kk) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set_nz(const MX &m, bool ind1, const Slice &kk) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::shared(std::vector< MX > &ex, std::vector< MX > &v, std::vector< MX > &vdef, const std::string &v_prefix, const std::string &v_suffix) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::simplify(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::solve(const MX &a, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::solve(const MX &a, const MX &b, const std::string &lsolver, const Dict &dict=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sparsity() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sparsity_cast(const MX &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::split_primitives(const DM &x) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::split_primitives(const MX &x) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::split_primitives(const SX &x) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::stop_diff(const MX &expr, casadi_int order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::stop_diff(const MX &expr, const MX &var, casadi_int order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::substitute(const MX &ex, const MX &v, const MX &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::substitute_inplace(const std::vector< MX > &v, std::vector< MX > &vdef, std::vector< MX > &ex, bool reverse) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sum1(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::sum2(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::symvar(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::trace(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::unary(casadi_int op, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::unite(const MX &A, const MX &B) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::vertcat(const std::vector< MX > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::vertsplit(const MX &x, const std::vector< casadi_int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::which_depends(const MX &expr, const MX &var, casadi_int order=1, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::which_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::which_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::adj(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::all(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::any(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::chol(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::cofactor(const Matrix< Scalar > &x, casadi_int i, casadi_int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::eig_symbolic(const Matrix< Scalar > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::evalf(const Matrix< Scalar > &expr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::expand(const Matrix< Scalar > &ex, Matrix< Scalar > &weights, Matrix< Scalar > &terms) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order, const Matrix< Scalar > &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order=5) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::get_ptr(Matrix< Scalar > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::get_ptr(const Matrix< Scalar > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::heaviside(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::ldl(const Matrix< Scalar > &A, Matrix< Scalar > &D, Matrix< Scalar > &LT, std::vector< casadi_int > &p, bool amd=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::ldl_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &D, const Matrix< Scalar > &LT, const std::vector< casadi_int > &p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::minor(const Matrix< Scalar > &x, casadi_int i, casadi_int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::mtaylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, casadi_int order, const std::vector< casadi_int > &order_contributions) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::mtaylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, casadi_int order=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::norm_inf_mul(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::poly_coeff(const Matrix< Scalar > &f, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::poly_roots(const Matrix< Scalar > &p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::pw_const(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::pw_lin(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::qr(const Matrix< Scalar > &A, Matrix< Scalar > &Q, Matrix< Scalar > &R) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::qr_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &v, const Matrix< Scalar > &r, const Matrix< Scalar > &beta, const std::vector< casadi_int > &prinv, const std::vector< casadi_int > &pc, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::qr_sparse(const Matrix< Scalar > &A, Matrix< Scalar > &V, Matrix< Scalar > &R, Matrix< Scalar > &beta, std::vector< casadi_int > &prinv, std::vector< casadi_int > &pc, bool amd=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::ramp(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::rectangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::sparsify(const Matrix< Scalar > &A, double tol=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::taylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::taylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, casadi_int order=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix::triangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::T() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::adj(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::all(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::any(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::blockcat(const std::vector< std::vector< Matrix< Scalar > > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::chol(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::cofactor(const Matrix< Scalar > &A, casadi_int i, casadi_int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::conditional(const Matrix< Scalar > &ind, const std::vector< Matrix< Scalar > > &x, const Matrix< Scalar > &x_default, bool short_circuit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::cse(const std::vector< Matrix< Scalar > > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::cumsum(const Matrix< Scalar > &x, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::densify(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::densify(const Matrix< Scalar > &x, const Matrix< Scalar > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::depends_on(const Matrix< Scalar > &x, const Matrix< Scalar > &arg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::det(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::diag(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::diagcat(const std::vector< Matrix< Scalar > > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::diagsplit(const Matrix< Scalar > &x, const std::vector< casadi_int > &offset1, const std::vector< casadi_int > &offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::dot(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::eig_symbolic(const Matrix< Scalar > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::einstein(const Matrix< Scalar > &A, const Matrix< Scalar > &B, const Matrix< Scalar > &C, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::einstein(const Matrix< Scalar > &A, const Matrix< Scalar > &B, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::enlarge(casadi_int nrow, casadi_int ncol, const std::vector< casadi_int > &rr, const std::vector< casadi_int > &cc, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::erase(const std::vector< casadi_int > &rr, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::erase(const std::vector< casadi_int > &rr, const std::vector< casadi_int > &cc, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::evalf(const Matrix< Scalar > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::expand(const Matrix< Scalar > &x, Matrix< Scalar > &weights, Matrix< Scalar > &terms) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::expm(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::expm_const(const Matrix< Scalar > &A, const Matrix< Scalar > &t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::export_code(const std::string &lang, std::ostream &stream=casadi::uout(), const Dict &options=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::extract(std::vector< Matrix< Scalar >> &ex, std::vector< Matrix< Scalar >> &v, std::vector< Matrix< Scalar >> &vdef, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::extract_parametric(const Matrix< Scalar > &expr, const Matrix< Scalar > &par, Matrix< Scalar > &expr_ret, std::vector< Matrix< Scalar > > &symbols, std::vector< Matrix< Scalar >> &parametric, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::forward(const std::vector< Matrix< Scalar > > &ex, const std::vector< Matrix< Scalar > > &arg, const std::vector< std::vector< Matrix< Scalar > > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order, const Matrix< Scalar > &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order=5) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Slice &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Slice &rr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Slice &rr, const Matrix< casadi_int > &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Slice &rr, const Slice &cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get(Matrix< Scalar > &m, bool ind1, const Sparsity &sp) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_elements() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_nonzeros() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_nz(Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_nz(Matrix< Scalar > &m, bool ind1, const Slice &k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_sparsity() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::get_str(bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::has_nz(casadi_int rr, casadi_int cc) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::has_zeros() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::heaviside(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::hessian(const Matrix< Scalar > &f, const Matrix< Scalar > &x, Matrix< Scalar > &g, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::hessian(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::horzcat(const std::vector< Matrix< Scalar > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::horzsplit(const Matrix< Scalar > &x, const std::vector< casadi_int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::if_else(const Matrix< Scalar > &x, const Matrix< Scalar > &if_true, const Matrix< Scalar > &if_false, bool short_circuit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inf(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inf(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inf(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inv(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inv(const Matrix< Scalar > &A, const std::string &lsolver, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inv_minor(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_constant() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_equal(const Matrix< Scalar > &x, const Matrix< Scalar > &y, casadi_int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_eye() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_integer() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_minus_one() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_one() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_zero() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::jacobian(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::jacobian_sparsity(const Matrix< Scalar > &f, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::kron(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::ldl(const Matrix< Scalar > &A, Matrix< Scalar > &D, Matrix< Scalar > &LT, std::vector< casadi_int > &p, bool amd=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::ldl_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &D, const Matrix< Scalar > &LT, const std::vector< casadi_int > &p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mac(const Matrix< Scalar > &x, const Matrix< Scalar > &y, const Matrix< Scalar > &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::minor(const Matrix< Scalar > &x, casadi_int i, casadi_int j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mldivide(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mmax(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mmin(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mrdivide(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mtaylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, casadi_int order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mtaylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, casadi_int order, const std::vector< casadi_int > &order_contributions) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::mtimes(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::n_nodes(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nan(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nan(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nan(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nonzeros() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nonzeros() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::norm_1(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::norm_2(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::norm_fro(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::norm_inf(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::norm_inf_mul(const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::nullspace(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::pinv(const Matrix< Scalar > &A, const std::string &lsolver, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::pinv(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::poly_coeff(const Matrix< Scalar > &ex, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::poly_roots(const Matrix< Scalar > &p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::polyval(const Matrix< Scalar > &p, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_dense(std::ostream &stream, bool truncate=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_operator(const Matrix< Scalar > &x, const std::vector< std::string > &args) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_scalar(std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_sparse(std::ostream &stream, bool truncate=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_split(std::vector< std::string > &nz, std::vector< std::string > &inter) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::print_vector(std::ostream &stream, bool truncate=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::printme(const Matrix< Scalar > &y) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::project(const Matrix< Scalar > &x, const Sparsity &sp, bool intersect=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::ptr() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::ptr() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::pw_const(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::pw_lin(const Matrix< Scalar > &t, const Matrix< Scalar > &tval, const Matrix< Scalar > &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::qr(const Matrix< Scalar > &A, Matrix< Scalar > &Q, Matrix< Scalar > &R) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::qr_solve(const Matrix< Scalar > &b, const Matrix< Scalar > &v, const Matrix< Scalar > &r, const Matrix< Scalar > &beta, const std::vector< casadi_int > &prinv, const std::vector< casadi_int > &pc, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::qr_sparse(const Matrix< Scalar > &A, Matrix< Scalar > &V, Matrix< Scalar > &R, Matrix< Scalar > &beta, std::vector< casadi_int > &prinv, std::vector< casadi_int > &pc, bool amd=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::ramp(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::rand(casadi_int nrow=1, casadi_int ncol=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::rand(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::rand(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::rectangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::remove(const std::vector< casadi_int > &rr, const std::vector< casadi_int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reserve(casadi_int nnz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reserve(casadi_int nnz, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reshape(const Matrix< Scalar > &x, casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reshape(const Matrix< Scalar > &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::resize(casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reverse(const std::vector< Matrix< Scalar > > &ex, const std::vector< Matrix< Scalar > > &arg, const std::vector< std::vector< Matrix< Scalar > > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::scalar() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::separate_linear(const Matrix< Scalar > &expr, const Matrix< Scalar > &sym_lin, const Matrix< Scalar > &sym_const, Matrix< Scalar > &expr_const, Matrix< Scalar > &expr_lin, Matrix< Scalar > &expr_nonlin) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::serialize() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::serialize(std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Matrix< casadi_int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Slice &rr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Slice &rr, const Matrix< casadi_int > &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Slice &rr, const Slice &cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set(const Matrix< Scalar > &m, bool ind1, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set_nz(const Matrix< Scalar > &m, bool ind1, const Matrix< casadi_int > &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::set_nz(const Matrix< Scalar > &m, bool ind1, const Slice &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::shared(std::vector< Matrix< Scalar > > &ex, std::vector< Matrix< Scalar > > &v, std::vector< Matrix< Scalar > > &vdef, const std::string &v_prefix, const std::string &v_suffix) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::simplify(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::solve(const Matrix< Scalar > &A, const Matrix< Scalar > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::solve(const Matrix< Scalar > &A, const Matrix< Scalar > &b, const std::string &lsolver, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::sparsify(const Matrix< Scalar > &x, double tol=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::sparsity() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::sparsity_cast(const Matrix< Scalar > &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::substitute(const Matrix< Scalar > &ex, const Matrix< Scalar > &v, const Matrix< Scalar > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::substitute(const std::vector< Matrix< Scalar > > &ex, const std::vector< Matrix< Scalar > > &v, const std::vector< Matrix< Scalar > > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::substitute_inplace(const std::vector< Matrix< Scalar > > &v, std::vector< Matrix< Scalar > > &vdef, std::vector< Matrix< Scalar > > &ex, bool revers) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::sum1(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::sum2(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::symvar(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::taylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x, const Matrix< Scalar > &a, casadi_int order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::to_file(const std::string &filename, const std::string &format="") const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::trace(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triplet(const std::vector< casadi_int > &row, const std::vector< casadi_int > &col, const Matrix< Scalar > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triplet(const std::vector< casadi_int > &row, const std::vector< casadi_int > &col, const Matrix< Scalar > &d, casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triplet(const std::vector< casadi_int > &row, const std::vector< casadi_int > &col, const Matrix< Scalar > &d, const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::unite(const Matrix< Scalar > &A, const Matrix< Scalar > &B) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::vertcat(const std::vector< Matrix< Scalar > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::vertsplit(const Matrix< Scalar > &x, const std::vector< casadi_int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::which_depends(const Matrix< Scalar > &expr, const Matrix< Scalar > &var, casadi_int order=1, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpBuilder::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpBuilder::get_str(bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpBuilder::import_nl(const std::string &filename, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::NlpBuilder::type_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::advanced() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::callback_class() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::copy() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::debug() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::dual(const MX &m) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::f() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::g() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::get_str(bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::initial() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::lam_g() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::lbg() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::minimize(const MX &f, double linear_scale=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::ng() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::np() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::nx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::p() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::parameter(casadi_int n=1, casadi_int m=1, const std::string &attribute="full") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::parameter(const MX &symbol, const std::string &attribute="full") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::parameter(const Sparsity &sp, const std::string &attribute="full") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::return_status() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::set_domain(const MX &x, const std::string &domain) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::set_initial(const MX &x, const DM &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::set_initial(const std::vector< MX > &assignments) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::set_linear_scale(const MX &x, const DM &scale, const DM &offset=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::set_value(const MX &x, const DM &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::set_value(const std::vector< MX > &assignments) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::solve() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::solve_limited() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::solver(const std::string &solver, const Dict &plugin_options=Dict(), const Dict &solver_options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::stats() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::subject_to() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::subject_to(const MX &g, const DM &linear_scale, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::subject_to(const MX &g, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::subject_to(const std::vector< MX > &g, const DM &linear_scale, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::subject_to(const std::vector< MX > &g, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::to_function(const std::string &name, const std::map< std::string, MX > &dict, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::to_function(const std::string &name, const std::vector< MX > &args, const std::vector< MX > &res, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::to_function(const std::string &name, const std::vector< MX > &args, const std::vector< MX > &res, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::type_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::ubg() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::update_user_dict(const MX &m, const Dict &meta) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::update_user_dict(const std::vector< MX > &m, const Dict &meta) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::user_dict(const MX &m) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::value(const DM &x, const std::vector< MX > &values=std::vector< MX >()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::value(const MX &x, const std::vector< MX > &values=std::vector< MX >()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::value_parameters() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::value_variables() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::variable(casadi_int n=1, casadi_int m=1, const std::string &attribute="full") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::variable(const MX &symbol, const std::string &attribute="full") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::variable(const Sparsity &sp, const std::string &attribute="full") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Opti::x() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::active_symvar(VariableType type) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::active_values(VariableType type) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::arg() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::assert_active_symbol(const MX &m) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::assert_baked() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::assert_empty() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::assert_solved() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::bake() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::baked_copy() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::canon_expr(const MX &expr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::casadi_solver() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::constraints() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::describe(const MX &x, casadi_index indent=0, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::g_describe(casadi_index i, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::g_index_reduce_g(casadi_index i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::g_index_reduce_x(casadi_index i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::g_index_unreduce_g(casadi_index i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::g_lookup(casadi_index i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::get_meta(const MX &m) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::get_meta_con(const MX &m) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::instance_number() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::is_parametric(const MX &expr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::mark_problem_dirty(bool flag=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::mark_solved(bool flag=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::mark_solver_dirty(bool flag=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::objective() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::problem_dirty() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::res() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::res(const DMDict &res) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::set_meta(const MX &m, const MetaVar &meta) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::set_meta_con(const MX &m, const MetaCon &meta) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::show_infeasibilities(double tol=0, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::solve_actual(const DMDict &args) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::solve_prepare() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::solved() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::solver_dirty() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::symvar() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::symvar(const MX &expr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::symvar(const MX &expr, VariableType type) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::x_describe(casadi_index i, const Dict &opts=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiAdvanced::x_lookup(casadi_index i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiCallback::call(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::get_str(bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::opti() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::stats() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::type_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::value(const DM &x, const std::vector< MX > &values=std::vector< MX >()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::value(const MX &x, const std::vector< MX > &values=std::vector< MX >()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::value(const SX &x, const std::vector< MX > &values=std::vector< MX >()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::value_parameters() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OptiSol::value_variables() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::anti_derivative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::coeff() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::degree() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::derivative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::scalar() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::trim() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Polynomial::type_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Printable::repr(const Derived &obj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::dep(casadi_int ch=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::element_hash() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::get_output(casadi_int oind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::has_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::info() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_call() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_commutative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_leaf() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_op(casadi_int op) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_regular() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_smooth() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_symbolic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_valid_input() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::n_dep() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::op() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::which_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::which_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::dep(casadi_int ch=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::get_output(casadi_int oind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::has_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::if_else(const SXElem &x, const SXElem &y, const SXElem &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::inv() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_almost_zero(double tol) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_call() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_commutative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_constant() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_doubled() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_inf() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_integer() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_leaf() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_minus_inf() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_minus_one() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_nan() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_nonnegative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_null() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_one() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_op(casadi_int op) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_regular() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_symbolic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_zero() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::n_dep() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::op() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::which_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::which_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::connect(DeserializerBase &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const Function &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const GenericType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const Linsol &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const MX &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const Matrix< SXElem > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const Matrix< double > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const Sparsity &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const casadi_int &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const double &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::string &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< Function > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< GenericType > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< Linsol > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< MX > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< Matrix< SXElem > > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< Matrix< double > > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< Sparsity > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< casadi_int > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< double > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::pack(const std::vector< std::string > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializerBase::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::connect(DeserializingStream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(bool e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(casadi_int e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(char e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const Fmu &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const Function &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const GenericType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const Importer &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const Linsol &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const MX &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const Matrix< T > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const SXElem &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const Slice &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const Sparsity &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const std::map< K, V > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const std::pair< A, B > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const std::string &descr, T &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const std::string &descr, const T &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const std::string &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(const std::vector< T > &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(double e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(int e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(size_t e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::pack(std::istream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SerializingStream::version(const std::string &name, int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::class_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::get_str(bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::all() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::all(casadi_int len, bool ind1=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::all(const Slice &outer, casadi_int len) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::apply(casadi_int len, bool ind1=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::get_str(bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::info() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::is_empty() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::is_scalar(casadi_int len) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::scalar(casadi_int len) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::size() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Slice::type_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::blockcat(const MatType &A, const MatType &B, const MatType &C, const MatType &D) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::blockcat(const std::vector< std::vector< MatType > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::blocksplit(const MatType &x, casadi_int vert_incr=1, casadi_int horz_incr=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::blocksplit(const MatType &x, const std::vector< casadi_int > &vert_offset, const std::vector< casadi_int > &horz_offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagcat(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagcat(const MatType &x, const MatType &y, const MatType &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagcat(const std::vector< MatType > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagsplit(const MatType &x, casadi_int incr1, casadi_int incr2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagsplit(const MatType &x, casadi_int incr=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagsplit(const MatType &x, const std::vector< casadi_int > &output_offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::diagsplit(const MatType &x, const std::vector< casadi_int > &output_offset1, const std::vector< casadi_int > &output_offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzcat(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzcat(const MatType &x, const MatType &y, const MatType &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzcat(const std::vector< MatType > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzsplit(const MatType &x, casadi_int incr=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzsplit(const MatType &x, const std::vector< casadi_int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::horzsplit_n(const MatType &x, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::kron(const MatType &a, const MatType &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::mac(const MatType &x, const MatType &y, const MatType &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::mtimes(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::mtimes(const std::vector< MatType > &args) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::norm_0_mul(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::offset(const std::vector< MatType > &v, bool vert=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::repmat(const MatType &A, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::repmat(const MatType &A, const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::reshape(const MatType &x, casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::reshape(const MatType &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::reshape(const MatType &x, std::pair< casadi_int, casadi_int > rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::sparsity_cast(const MatType &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::sprank(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::sum1(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::sum2(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::transpose(const MatType &X) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::tril(const MatType &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::triu(const MatType &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vec(const MatType &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::veccat(const std::vector< MatType > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertcat(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertcat(const MatType &x, const MatType &y, const MatType &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertcat(const MatType &x, const MatType &y, const MatType &z, const MatType &w, const MatType &v, const MatType &u) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertcat(const std::vector< MatType > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertsplit(const MatType &x, casadi_int incr=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertsplit(const MatType &x, const std::vector< casadi_int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparsityInterface::vertsplit_n(const MatType &x, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StringDeserializer::decode(const std::string &string) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::StringSerializer::encode() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakCache< K, T >::cache(std::vector< K > &keys, std::vector< T > &entries) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakCache< K, T >::incache(const K &key, T &f, bool needs_lock=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakCache< K, T >::tocache(const K &key, const T &f, bool needs_lock=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::WeakCache< K, T >::tocache_if_missing(const K &key, T &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlFile::dump(const std::string &filename, const XmlNode &node) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlFile::parse(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::all(const std::vector< bool > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::any(const std::vector< bool > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::assign_vector(const std::vector< S > &s, std::vector< D > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::blazing_spline(const std::string &name, const std::vector< std::vector< double > > &knots, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::boolvec_and(const std::vector< bool > &lhs, const std::vector< bool > &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::boolvec_not(const std::vector< bool > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::boolvec_or(const std::vector< bool > &lhs, const std::vector< bool > &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::boolvec_to_index(const std::vector< bool > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::check_sos(casadi_int nx, const std::vector< std::vector< T > > &groups, std::vector< std::vector< double > > &weights, std::vector< casadi_int > &types) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocation_coeff(const std::vector< double > &tau, DM &C, DM &D, DM &B) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocation_interpolators(const std::vector< double > &tau, std::vector< std::vector< double > > &C, std::vector< double > &D) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocation_points(casadi_int order, const std::string &scheme="radau") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::collocation_pointsL(casadi_int order, const std::string &scheme) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::combine(const Dict &first, const Dict &second, bool recurse=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::complement(const std::vector< casadi_int > &v, casadi_int size) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic(const std::string &name, const std::string &solver, const SpDict &qp, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_debug(const Function &f, const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_debug(const Function &f, std::ostream &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_in(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_option_info(const std::string &name, const std::string &op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_option_type(const std::string &name, const std::string &op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_options(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::conic_out(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::copy_vector(const std::vector< S > &s, std::vector< D > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum0(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dae_init_gen(const MXDict &dae, const MXDict &dae_red, const std::string &init_solver, const DMDict &init_strength=DMDict(), const Dict &init_solver_options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dae_init_gen(const SXDict &dae, const SXDict &dae_red, const std::string &init_solver, const DMDict &init_strength, const Dict &init_solver_options) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dae_map_semi_expl(const MXDict &dae, const MXDict &dae_red, Function &state_to_orig, Function &phi) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dae_map_semi_expl(const SXDict &dae, const SXDict &dae_red, Function &state_to_orig, Function &phi) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dae_reduce_index(const MXDict &dae, Dict &stats, const Dict &opts={}) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dae_reduce_index(const SXDict &dae, Dict &stats, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::detect_simple_bounds(const MX &x, const MX &p, const MX &g, const MX &lbg, const MX &ubg, std::vector< casadi_int > &gi, MX &lbx, MX &ubx, Function &lam_forward, Function &lam_backward) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::detect_simple_bounds(const SX &xX, const SX &p, const SX &g, const SX &lbg, const SX &ubg, std::vector< casadi_int > &gi, SX &lbx, SX &ubx, Function &lam_forward, Function &lam_backward) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::diff(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_conic(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_dple(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_expm(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_integrator(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_interpolant(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_linsol(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_nlpsol(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::doc_rootfinder(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dot(const std::vector< T > &a, const std::vector< T > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dple_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dple_in(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dple_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dple_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dple_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dple_out(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dplesol(const DMVector &A, const DMVector &V, const std::string &solver, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dplesol(const MX &A, const MX &V, const std::string &solver, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dplesol(const MXVector &A, const MXVector &V, const std::string &solver, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dplesol(const std::string &name, const std::string &solver, const SpDict &st, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dyn_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dyn_in(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dyn_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dyn_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dyn_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dyn_out(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::einstein_process(const T &A, const T &B, const T &C, const std::vector< casadi_int > &dim_a, const std::vector< casadi_int > &dim_b, const std::vector< casadi_int > &dim_c, const std::vector< casadi_int > &a, const std::vector< casadi_int > &b, const std::vector< casadi_int > &c, std::vector< casadi_int > &iter_dims, std::vector< casadi_int > &strides_a, std::vector< casadi_int > &strides_b, std::vector< casadi_int > &strides_c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::event_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::event_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::expm_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::expm_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::expmsol(const std::string &name, const std::string &solver, const Sparsity &A, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::external(const std::string &name, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::external(const std::string &name, const Importer &li, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::external(const std::string &name, const std::string &bin_name, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::external_transform(const std::string &name, const std::string &op, const Function &f, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::extract_from_dict(const Dict &d, const std::string &key, T &value) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::extract_from_dict_inplace(Dict &d, const std::string &key, T &value) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::find(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::flatten_nested_vector(const std::vector< std::vector< T > > &nested, std::vector< S > &flat) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::flatten_nested_vector(const std::vector< std::vector< T > > &nested, std::vector< S > &flat, std::vector< I > &indices) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::fmtstr(const std::string &fmt, const std::vector< std::string > &args) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_bvec_t(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_from_dict(const Dict &d, const std::string &key, const T &default_value) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_from_dict(const std::map< std::string, T > &d, const std::string &key, const T &default_value) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_ptr(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::get_ptr(std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_conic(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_dple(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_expm(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_integrator(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_interpolant(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_linsol(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_negative(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_nlpsol(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::has_rootfinder(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::in_range(const std::vector< T > &v, casadi_int lower, casadi_int upper) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::in_range(const std::vector< T > &v, casadi_int upper) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::index_interp1d(const std::vector< double > &x, double xq, bool equidistant) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::init_vector(std::vector< S > &d, const std::vector< D > &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, double t0, const std::vector< double > &tout, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, double t0, double tf, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const MXDict &dae, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const MXDict &dae, double t0, const std::vector< double > &tout, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const MXDict &dae, double t0, double tf, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const SXDict &dae, double t0, const std::vector< double > &tout, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const SXDict &dae, double t0, double tf, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator_in(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator_out(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::interpolant(const std::string &name, const std::string &solver, const std::vector< casadi_int > &grid_dims, casadi_int m=1, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::interpolant(const std::string &name, const std::string &solver, const std::vector< casadi_int > &grid_dims, const std::vector< double > &values, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::interpolant(const std::string &name, const std::string &solver, const std::vector< std::vector< double > > &grid, casadi_int m=1, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::interpolant(const std::string &name, const std::string &solver, const std::vector< std::vector< double > > &grid, const std::vector< double > &values, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::invert_permutation(const std::vector< casadi_int > &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isUnique(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_decreasing(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_equally_spaced(const std::vector< double > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_increasing(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_monotone(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_nondecreasing(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_nonincreasing(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_permutation(const std::vector< casadi_int > &order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_range(const std::vector< casadi_int > &v, casadi_int start, casadi_int stop, casadi_int step=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_regular(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_slice(const IM &x, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_slice(const std::vector< casadi_int > &v, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_slice2(const std::vector< casadi_int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_strictly_monotone(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_zero(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::join(const std::vector< T > &a, const std::vector< T > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::join(const std::vector< T > &a, const std::vector< T > &b, const std::vector< T > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::join(const std::vector< std::string > &l, const std::string &delim) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::linspace(std::vector< T > &v, const F &first, const L &last) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_conic(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_dple(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_expm(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_integrator(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_interpolant(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_linsol(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_nlpsol(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::load_rootfinder(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::lookupvector(const std::vector< casadi_int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::lookupvector(const std::vector< casadi_int > &v, casadi_int size) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName< SXElem >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName< casadi_int >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName< double >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::message_prefix(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol(const std::string &name, const std::string &solver, const Function &nlp, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol(const std::string &name, const std::string &solver, const Importer &compiler, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol(const std::string &name, const std::string &solver, const MXDict &nlp, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol(const std::string &name, const std::string &solver, const NlpBuilder &nl, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol(const std::string &name, const std::string &solver, const SXDict &nlp, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol(const std::string &name, const std::string &solver, const std::string &fname, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_default_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_default_in(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_in(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_option_info(const std::string &name, const std::string &op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_option_type(const std::string &name, const std::string &op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_options(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol_out(casadi_int ind) {
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
%exception  casadi::normalized_in(std::istream &stream, double &ret) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::normalized_out(std::ostream &stream, double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::normalized_setup(std::istream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::normalized_setup(std::ostream &stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::permute(const std::vector< T > &a, const std::vector< casadi_int > &order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::product(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::qpsol(const std::string &name, const std::string &solver, const MXDict &qp, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::qpsol(const std::string &name, const std::string &solver, const SXDict &qp, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::range(casadi_int start, casadi_int stop, casadi_int step=1, casadi_int len=std::numeric_limits< casadi_int >::max()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::range(casadi_int stop) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::read_matlab(std::ifstream &file, std::vector< std::vector< T > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::read_matlab(std::istream &stream, std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::replace(const std::string &s, const std::string &p, const std::string &r) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::reverse(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder(const std::string &name, const std::string &solver, const Function &f, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder(const std::string &name, const std::string &solver, const MXDict &rfp, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder(const std::string &name, const std::string &solver, const SXDict &rfp, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_in(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_option_info(const std::string &name, const std::string &op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_option_type(const std::string &name, const std::string &op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_options(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::rootfinder_out(casadi_int ind) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::shared_cast(S &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::shared_cast(const S &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::simpleIRK(Function f, casadi_int N=10, casadi_int order=4, const std::string &scheme="radau", const std::string &solver="newton", const Dict &solver_options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::simpleIntegrator(Function f, const std::string &integrator="cvodes", const Dict &integrator_options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::simpleRK(Function f, casadi_int N=10, casadi_int order=4) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< casadi_int > &indices, bool invert_indices=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::startswith(const std::string &s, const std::string &p) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const T &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const T &v, bool more) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const std::array< T, N > &p, bool more=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const std::map< T1, T2 > &p, bool more=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const std::map< std::string, T2 > &p, bool more=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const std::pair< T1, T2 > &p, bool more=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const std::set< T > &v, bool more=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const std::vector< T > &v, bool more=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str_bvec(bvec_t v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::strvec() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::strvec(const T1 &t1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::strvec(const T1 &t1, const T2 &t2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::strvec(const T1 &t1, const T2 &t2, const T3 &t3) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::strvec(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::strvec(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::strvec(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sum(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::temporary_file(const std::string &prefix, const std::string &suffix) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::tensor_permute_mapping(const std::vector< casadi_int > &dims, const std::vector< casadi_int > &order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::text2set(const std::string &text) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::text2type(const std::string &text) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::text2vector(const std::string &text) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_int(casadi_int rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_int(const std::vector< casadi_int > &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_int(const std::vector< std::vector< casadi_int > > &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_slice(const IM &x, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_slice(const std::vector< casadi_int > &v, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_slice2(const std::vector< casadi_int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_string(DynIn v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_string(DynOut v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_string(EventIn v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_string(EventOut v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::trim_path(const std::string &full_path) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::uerr() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::uout() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::update_dict(Dict &target, const Dict &source, bool recurse=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_init(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_select(const std::vector< T > &v, const std::vector< bool > &s, bool invert=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_slice(const std::vector< T > &v, const std::vector< casadi_int > &i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_static_cast(const std::vector< S > &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_tail(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::write_matlab(std::ostream &stream, const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::write_matlab(std::ostream &stream, const std::vector< std::vector< T > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  std::mul_overflows(const T &a, const T &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Callback::Callback() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Callback::Callback(const Callback &obj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CasadiException::CasadiException() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CasadiException::CasadiException(const std::string &msg) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CodeGenerator::CodeGenerator(const std::string &name, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DaeBuilder::DaeBuilder() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DaeBuilder::DaeBuilder(const std::string &name, const std::string &path="", const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DeserializerBase::DeserializerBase(std::unique_ptr< std::istream > stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DeserializingStream::DeserializingStream(const DeserializingStream &)=delete {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DeserializingStream::DeserializingStream(std::istream &in_s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FileDeserializer::FileDeserializer(const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FileSerializer::FileSerializer(const std::string &fname, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, MXIList ex_in, MXIList ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, MXIList ex_in, MXIList ex_out, const StringVector &name_in, const StringVector &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, MXIList ex_in, const MXVector &ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, MXIList ex_in, const MXVector &ex_out, const StringVector &name_in, const StringVector &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, SXIList ex_in, SXIList ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, SXIList ex_in, SXIList ex_out, const StringVector &name_in, const StringVector &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, SXIList ex_in, const SXVector &ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, SXIList ex_in, const SXVector &ex_out, const StringVector &name_in, const StringVector &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const MXVector &ex_in, MXIList ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const MXVector &ex_in, MXIList ex_out, const StringVector &name_in, const StringVector &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const SXVector &ex_in, SXIList ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const SXVector &ex_in, SXIList ex_out, const StringVector &name_in, const StringVector &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const std::map< std::string, MX > &dict, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const std::map< std::string, SX > &dict, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const std::vector< MX > &ex_in, const std::vector< MX > &ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const std::vector< MX > &ex_in, const std::vector< MX > &ex_out, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const std::vector< SX > &ex_in, const std::vector< SX > &ex_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Function::Function(const std::string &name, const std::vector< SX > &ex_in, const std::vector< SX > &ex_out, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FunctionBuffer::FunctionBuffer(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FunctionBuffer::FunctionBuffer(const FunctionBuffer &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericShared< Shared, Internal >::GenericShared() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericShared< Shared, Internal >::GenericShared(const GenericShared &ref) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(bool b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const Dict &dict) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const char s[]) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< Dict > &dictv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< Function > &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< GenericType > &gv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< bool > &iv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< casadi_int > &iv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< double > &dv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< int > &iv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< std::string > &sv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< std::vector< GenericType > > &gvv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< std::vector< casadi_int > > &ivv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< std::vector< double > > &dv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< std::vector< std::string > > &sv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(double d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericWeakRef< Shared, Internal >::GenericWeakRef(Shared shared) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericWeakRef< Shared, Internal >::GenericWeakRef(int dummy=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Importer::Importer() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Importer::Importer(const std::string &name, const std::string &compiler, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::KeyboardInterruptException::KeyboardInterruptException() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Linsol::Linsol() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Linsol::Linsol(const std::string &name, const std::string &solver, const Sparsity &sp, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Logger::Stream< Err >::Stream() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Logger::Streambuf< Err >::Streambuf() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Matrix< double > &val, const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Matrix< double > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp, const MX &val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp, const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const Sparsity &sp, double val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(const std::vector< double > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::MX::MX(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Matrix< A > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Matrix< Scalar > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const Matrix< Scalar > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const Scalar &val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const std::vector< Scalar > &d, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::vector< A > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::vector< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::vector< std::vector< double > > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(std::initializer_list< Scalar > x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlImporter::NlImporter(NlpBuilder &nlp, const std::string &filename, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NonZeros< M, K >::NonZeros(M &mat, const K &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NonZeros< M, K >::NonZeros(const NonZeros< M, K > &y)=default {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Opti::Opti(const Opti &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Opti::Opti(const std::string &problem_type="nlp") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OptiAdvanced::OptiAdvanced(const Opti &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OptiCallback::OptiCallback() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OptiCallback::OptiCallback(const OptiCallback &obj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Polynomial::Polynomial(const std::vector< T > &coeff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Polynomial::Polynomial(double p0, double p1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Polynomial::Polynomial(double p0, double p1, double p2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Polynomial::Polynomial(double p0, double p1, double p2, double p3) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Polynomial::Polynomial(double scalar=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElem::SXElem() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElem::SXElem(const SXElem &scalar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SXElem::SXElem(double val) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SerializerBase::SerializerBase(std::unique_ptr< std::ostream > stream, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SerializingStream::SerializingStream(std::ostream &out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SerializingStream::SerializingStream(std::ostream &out, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Slice::Slice() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Slice::Slice(casadi_int i, bool ind1=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Slice::Slice(casadi_int start, casadi_int stop, casadi_int step=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Slice::Slice(casadi_int start, int stop, int step=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Slice::Slice(int start, casadi_int stop, int step=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Slice::Slice(int start, int stop, int step=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StringDeserializer::StringDeserializer(const std::string &string) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::StringSerializer::StringSerializer(const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubIndex< M, I >::SubIndex(M &mat, const I &i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubIndex< M, I >::SubIndex(const SubIndex< M, I > &y)=default {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubMatrix< M, I, J >::SubMatrix(M &mat, const I &i, const J &j) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SubMatrix< M, I, J >::SubMatrix(const SubMatrix< M, I, J > &y)=default {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(SharedObject shared) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::WeakRef::WeakRef(int dummy=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::XmlFile::XmlFile() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::XmlFile::XmlFile(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::conditional_lock_guard< _Mutex >::conditional_lock_guard(const conditional_lock_guard &)=delete {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::conditional_lock_guard< _Mutex >::conditional_lock_guard(mutex_type &m, bool condition) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::null_ptr_on_copy< T >::null_ptr_on_copy() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::null_ptr_on_copy< T >::null_ptr_on_copy(const null_ptr_on_copy &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::scoped_checkout< T >::scoped_checkout(const T &proto) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::scoped_checkout< T >::scoped_checkout(const scoped_checkout &that)=delete {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::scoped_checkout< T >::scoped_checkout(scoped_checkout &&that) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}