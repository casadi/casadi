%exception  casadi::BSplineInterpolant::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::get_diff_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::has_jacobian() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BackwardDiff::calc_stepsize(double abstol) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BackwardDiff::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Blocksqp::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Blocksqp::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Blocksqp::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Blocksqp::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Blocksqp::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Blocksqp::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BonminUserClass::branchingInfo() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BonminUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, TNLP::IndexStyleEnum &index_style) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BonminUserClass::get_number_of_nonlinear_variables() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BonminUserClass::sosConstraints() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::calc_fd() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::calc_stepsize(double abstol) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::get_abstol() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::has_err() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::n_pert() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::pert(casadi_int k, double h) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CentralDiff::pert(const std::string &k) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangCompiler::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangCompiler::get_function(const std::string &symname) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangCompiler::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangCompiler::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangCompiler::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangXCompiler::clone() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangXCompiler::getFunction(const std::string &symname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangXCompiler::init() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ClangXCompiler::plugin_name() const override {
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
%exception  casadi::CodeGenerator::add_io_sparsities(const std::string &name, const std::vector< Sparsity > &sp_in, const std::vector< Sparsity > &sp_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::add_sparsity(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::arg(casadi_int i) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::avoid_stack() {
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
%exception  casadi::CodeGenerator::comment(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(casadi_int v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::vector< casadi_int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::constant(const std::vector< double > &v) {
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
%exception  casadi::CodeGenerator::from_mex(std::string &arg, const std::string &res, std::size_t res_off, const Sparsity &sp_res, const std::string &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::get_constant(const std::vector< casadi_int > &v, bool allow_adding=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::get_constant(const std::vector< double > &v, bool allow_adding=false) {
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
%exception  casadi::CodeGenerator::initializer(const std::vector< casadi_int > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::initializer(const std::vector< double > &v) {
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
%exception  casadi::CodeGenerator::norm_inf(casadi_int n, const std::string &x) {
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
%exception  casadi::CodeGenerator::print_vector(std::ostream &s, const std::string &name, const std::vector< double > &v) {
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
%exception  casadi::CodeGenerator::shorthand(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::shorthand(const std::string &name, bool allow_adding=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sparsify(const std::string &arg, const std::string &res, const Sparsity &sp_res, bool tr=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sparsity(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sum_viol(casadi_int n, const std::string &x, const std::string &lb, const std::string &ub) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::sx_work(casadi_int i) {
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
%exception  casadi::CodeGenerator::work(casadi_int n, casadi_int sz) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CodeGenerator::workel(casadi_int n) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::algebraic_state_init(const MX &x0, const MX &z0) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::algebraic_state_output(const MX &Z) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Collocation::setupFG() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::generateNativeCode(std::ostream &file) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_default_in(casadi_int ind) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::integer_support() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::is_a(const std::string &type, bool recursive) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::psd_support() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConstantSX_deserialize(DeserializingStream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::codegen_qp_solve(CodeGenerator &cg, const std::string &H, const std::string &g, const std::string &lbdz, const std::string &ubdz, const std::string &A, const std::string &x_opt, const std::string &dlam) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::getConic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::print_iteration() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf, double dx_norm, double rg, casadi_int ls_trials, bool ls_success) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CustomNlpsol::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesSimulator::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesSimulator::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesSimulator::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesSimulator::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::CvodesSimulator::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::add_variable(const std::string &name, const Variable &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::find(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::find(const std::vector< std::string > &name) const {
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
%exception  casadi::DaeBuilder::var(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::var(const std::vector< size_t > &ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DaeBuilder::var(size_t ind) const {
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
%exception  casadi::DllLibrary::can_have_meta() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::finalize() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::get_function(const std::string &symname) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::init_handle() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::library() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Dple::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_diff_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::has_jac_sparsity(casadi_int oind, casadi_int iind) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::any_symbol_found() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_alloc_mem(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_checkout(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_decref(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_free_mem(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_incref(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_init_mem(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_mem_type() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_release(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::factory(const std::string &name, const std::vector< std::string > &s_in, const std::vector< std::string > &s_out, const Function::AuxOut &aux, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_default_in(casadi_int i) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::has_jacobian() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::init_external() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::join(FStats &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::tic() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::toc() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_dual(const Function::AuxOut &aux) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_input(const std::string &s, const MatType &e, bool is_diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_output(const std::string &s, const MatType &e, bool is_diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::block(const std::string &s1, const std::string &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate(const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate_adj(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate_fwd(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate_grad(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate_hess(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate_hess(const Dict &opts, size_t f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate_jac(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::find_hess(size_t f, size_t x1, size_t x2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::find_jac(size_t f, size_t x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::get_input(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::get_output(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::has_in(const std::string &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::has_out(const std::string &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::hblock(const std::string &s1, const std::string &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::imap(const std::string &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::iname() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::iname(const std::vector< size_t > &ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::omap(const std::string &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::oname() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::oname(const std::vector< size_t > &ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::request_input(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::request_output(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FastNewton::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_default_in(casadi_int ind) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FiniteDiff::uses_output() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::create_advanced(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::getExplicit() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::getExplicitB() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::setupFG() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::calc_fd() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::calc_stepsize(double abstol) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::get_abstol() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::has_err() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::n_pert() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::pert(casadi_int k, double h) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ForwardDiff::pert(const std::string &k) const override {
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
%exception  casadi::Function::call_gen(const std::map< std::string, M > &arg, std::map< std::string, M > &res, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::export_code(const std::string &lang, std::ostream &stream, const Dict &options=Dict()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::get() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Function::memory(int ind) const {
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
%exception  casadi::Function::sz_work(size_t &sz_arg, size_t &sz_res, size_t &sz_iw, size_t &sz_w) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::ad_weight() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::adjViaJac(casadi_int nadj) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::all_scalar() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::alloc(const Function &f, bool persistent=false, int num_threads=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::alloc_arg(size_t sz_arg, bool persistent=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::alloc_iw(size_t sz_iw, bool persistent=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::alloc_res(size_t sz_res, bool persistent=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::alloc_w(size_t sz_w, bool persistent=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call(const std::vector< M > &arg, std::vector< M > &res, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_forward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_forward(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_gen(const MXVector &arg, MXVector &res, casadi_int npar, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_gen(const std::vector< Matrix< D > > &arg, std::vector< Matrix< D > > &res, casadi_int npar, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_reverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_reverse(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::change_option(const std::string &option_name, const GenericType &option_value) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_arg(const std::vector< M > &arg, casadi_int &npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_mat(const Sparsity &arg, const Sparsity &inp, casadi_int &npar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_res(const std::vector< M > &res, casadi_int &npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen(CodeGenerator &g, const std::string &fname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_alloc_mem(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_body(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_checkout(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_declarations(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_decref(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_free_mem(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_incref(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_init_mem(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_mem(CodeGenerator &g, const std::string &index="mem") const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_mem_type() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_meta(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_name(const CodeGenerator &g, bool ns=true) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_release(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_sparsities(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_arg(const std::map< std::string, M > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_arg(const std::vector< M > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_res(const std::map< std::string, M > &res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_res(const std::vector< M > &res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::definition() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::disp(std::ostream &stream, bool more) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::disp_more(std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dm_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dm_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dm_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dm_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::eval_dm(const std::vector< DM > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::eval_mx(const MXVector &arg, MXVector &res, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::export_code(const std::string &lang, std::ostream &stream, const Dict &options) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::factory(const std::string &name, const std::vector< std::string > &s_in, const std::vector< std::string > &s_out, const Function::AuxOut &aux, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::finalize() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::forward(casadi_int nfwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::free_mx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::free_sx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::from_compact(casadi_int oind, casadi_int iind, const Sparsity &sp) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::fwdViaJac(casadi_int nfwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::fwd_seed(casadi_int nfwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generate_dependencies(const std::string &fname, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generate_lifted(Function &vdef_fcn, Function &vinit_fcn) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generate_options(bool is_temp=false) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getAdaptorSolverName() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_abstol() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_default_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_diff_in(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_diff_out(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_free() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_function(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jac_sparsity_gen(casadi_int oind, casadi_int iind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jac_sparsity_hierarchical(casadi_int oind, casadi_int iind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jac_sparsity_hierarchical_symm(casadi_int oind, casadi_int iind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_max_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_min_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_n_in() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_n_out() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_name_in(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_name_out(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_nominal_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_nominal_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_partition(casadi_int iind, casadi_int oind, Sparsity &D1, Sparsity &D2, bool compact, bool symmetric, bool allow_forward, bool allow_reverse) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_reltol() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_sparsity_in(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_sparsity_out(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_codegen() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_derivative() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_eval_dm() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_forward(casadi_int nfwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_free() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_function(const std::string &fname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_jac_sparsity(casadi_int oind, casadi_int iind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_jacobian() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_reverse(casadi_int nadj) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_spfwd() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_sprev() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::incache(const std::string &fname, Function &f, const std::string &suffix="") const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::index_in(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::index_out(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::info() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_MX(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_constant(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_id(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_input(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_output(casadi_int k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instructions_sx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::is_a(const std::string &type, bool recursive) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jac_is_symm(casadi_int oind, casadi_int iind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jac_sparsity(casadi_int oind, casadi_int iind, bool compact, bool symmetric) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacobian() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jit_dependencies(const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::map(casadi_int n, const std::string &parallelization) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mapsum_mx(const std::vector< MX > &arg, const std::string &parallelization) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::matching_arg(const std::vector< M > &arg, casadi_int &npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::matching_res(const std::vector< M > &arg, casadi_int &npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::n_instructions() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::n_nodes() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nz_in(const std::vector< DM > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nz_in(const std::vector< double > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nz_out(const std::vector< DM > &res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nz_out(const std::vector< double > &res) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::oracle() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::print_dimensions(std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::project_arg(const std::vector< M > &arg, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::project_res(const std::vector< M > &arg, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_arg(const std::vector< M > &arg, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_aseed(const std::vector< std::vector< M > > &aseed, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_aseed(const std::vector< std::vector< M >> &aseed, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_fseed(const std::vector< std::vector< M > > &fseed, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_fseed(const std::vector< std::vector< M >> &fseed, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_res(const std::vector< M > &res, casadi_int npar) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::reverse(casadi_int nadj) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::self() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::signature(const std::string &fname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size1_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size1_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size2_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size2_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::slice(const std::string &name, const std::vector< casadi_int > &order_in, const std::vector< casadi_int > &order_out, const Dict &opts) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sp_weight() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sparsity_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sparsity_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_in() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_in(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_out() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_out(casadi_int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicAdjSeed(casadi_int nadj, const std::vector< MatType > &v) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolic_output(const std::vector< MX > &arg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_arg() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_iw() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_res() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_w() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_work(size_t &sz_arg, size_t &sz_res, size_t &sz_iw, size_t &sz_w) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::to_compact(casadi_int oind, casadi_int iind, const Sparsity &sp) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::tocache(const Function &f, const std::string &suffix="") const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::uses_output() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::which_depends(const std::string &s_in, const std::vector< std::string > &s_out, casadi_int order, bool tr=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::wrap() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::wrap_as_needed(const Dict &opts) const {
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
%exception  casadi::GenericExpression::atan2(const ExType &x, const ExType &y) {
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
%exception  casadi::GenericExpression< ExType >::atan2(const ExType &x, const ExType &y) {
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
%exception  casadi::GenericExpression< SXElem  >::atan2(const SXElem &x, const SXElem &y) {
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
%exception  casadi::GenericExternal::any_symbol_found() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::get_diff_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::get_diff_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::get_jac_sparsity(casadi_int oind, casadi_int iind, bool symmetric) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::has_jac_sparsity(casadi_int oind, casadi_int iind) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::init_external() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::serialize_type(SerializingStream &s) const override {
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
%exception  casadi::GenericMatrix< MX  >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::cross(const MX &a, const MX &b, casadi_int dim=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::diff(const MX &x, casadi_int n=1, casadi_int axis=-1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::interp1d(const std::vector< double > &x, const MX &v, const std::vector< double > &xq, const std::string &mode, bool equidistant) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::inv_skew(const MX &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_linear(const MX &expr, const MX &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::is_quadratic(const MX &expr, const MX &var) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::linear_coeff(const MX &expr, const MX &var, MX &A, MX &b, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::linspace(const MX &a, const MX &b, casadi_int nsteps) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::norm_0_mul(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nz(const K &k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::nz(const K &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::quadratic_coeff(const MX &expr, const MX &var, MX &A, MX &b, MX &c, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::skew(const MX &a) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sparsity() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sprank(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MX  >::sumsqr(const MX &x) {
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
%exception  casadi::GenericMatrix< MatType >::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::interp1d(const std::vector< double > &x, const MatType &v, const std::vector< double > &xq, const std::string &mode, bool equidistant) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::linear_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::norm_0_mul(const MatType &x, const MatType &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nz(const K &k) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::nz(const K &k) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::quadratic_coeff(const MatType &expr, const MatType &var, MatType &A, MatType &b, MatType &c, bool check) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::repsum(const MatType &x, casadi_int n, casadi_int m=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix< MatType >::row() const {
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
%exception  casadi::GenericType::as_bool() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_bool_vector() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericType::as_dict() const {
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
%exception  casadi::GenericType::is_void_pointer() const {
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
%exception  casadi::GenericType::to_void_pointer() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericTypeBase::getType() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericTypeBase::serialize(SerializingStream &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicit() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::getExplicitB() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitFixedStepIntegrator::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitToNlp::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitToNlp::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitToNlp::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitToNlp::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImplicitToNlp::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Importer::get_function(const std::string &symname) {
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
%exception  casadi::Importer::to(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::body(const std::string &symname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::can_have_meta() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::construct(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::disp(std::ostream &stream, bool more) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::finalize() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::get_meta(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::get_options() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::has_function(const std::string &symname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::has_meta(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::init(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::inlined(const std::string &symname) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::library() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::read_external(const std::string &sym, bool inlined, std::istream &file, casadi_int &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::read_meta(std::istream &file, casadi_int &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::serialize_body(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::serialize_type(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::to_text(const std::string &cmd, casadi_int ind=-1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::algebraic_state_init(const MX &x0, const MX &z0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::algebraic_state_output(const MX &Z) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::aug_adj(casadi_int nadj) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::aug_fwd(casadi_int nfwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::create_advanced(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::getDerivativeOptions(bool fwd) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::has_spfwd() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::has_sprev() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::p() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::q() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rp() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rq() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rx() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rz() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::sp_jac_dae() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::sp_jac_rdae() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::t() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::x() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::z() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::arg_grid() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::arg_values() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::coeff_size() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_diff_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::has_parametric_grid() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::has_parametric_values() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::finalize_metadata(Index n, const StringMetaDataMapType &var_string_md, const IntegerMetaDataMapType &var_integer_md, const NumericMetaDataMapType &var_numeric_md, Index m, const StringMetaDataMapType &con_string_md, const IntegerMetaDataMapType &con_integer_md, const NumericMetaDataMapType &con_numeric_md) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_number_of_nonlinear_variables() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::IpoptUserClass::get_var_con_metadata(Index n, StringMetaDataMapType &var_string_md, IntegerMetaDataMapType &var_integer_md, NumericMetaDataMapType &var_numeric_md, Index m, StringMetaDataMapType &con_string_md, IntegerMetaDataMapType &con_integer_md, NumericMetaDataMapType &con_numeric_md) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Ipqp::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Ipqp::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Ipqp::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Ipqp::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Ipqp::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Ipqp::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::has_jacobian() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::JitFunction::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLu::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLu::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLu::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackLu::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQr::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQr::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQr::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQr::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQr::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LapackQr::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::has_jacobian() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolant::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::has_jacobian() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::has_parametric_grid() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::has_parametric_values() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::checkout() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Linsol::nfact(const DM &A) const {
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
%exception  casadi::LinsolInternal::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::disp(std::ostream &stream, bool more) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::disp_more(std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::generate(CodeGenerator &g, const std::string &A, const std::string &x, casadi_int nrhs, bool tr) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::ncol() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::nnz() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::nrow() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::generate(CodeGenerator &g, const std::string &A, const std::string &x, casadi_int nrhs, bool tr) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::finalize() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::generate(CodeGenerator &g, const std::string &A, const std::string &x, casadi_int nrhs, bool tr) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolTridiag::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolTridiag::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolTridiag::generate(CodeGenerator &g, const std::string &A, const std::string &x, casadi_int nrhs, bool tr) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolTridiag::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolTridiag::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Lsqr::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Lsqr::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Lsqr::generate(CodeGenerator &g, const std::string &A, const std::string &x, casadi_int nrhs, bool tr) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Lsqr::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::ad_forward(const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::ad_reverse(const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::blockcat(const std::vector< std::vector< MX > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::conditional(const MX &ind, const std::vector< MX > &x, const MX &x_default, bool short_circuit=false) {
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
%exception  casadi::MX::dot(const MX &x, const MX &y) {
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
%exception  casadi::MX::forward(const std::vector< MX > &ex, const std::vector< MX > &arg, const std::vector< std::vector< MX > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_temp() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const MX &x, const std::vector< MX > &v, const std::vector< MX > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const std::vector< MX > &ex, const std::vector< MX > &expr, const std::vector< MX > &exprs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::has_duplicates() const {
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
%exception  casadi::MX::inv(const MX &A, const std::string &lsolver="qr", const Dict &dict=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inv_minor(const MX &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inv_node(const MX &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::is_equal(const MX &x, const MX &y, casadi_int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::jacobian(const MX &f, const MX &x, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::kron(const MX &x, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mac(const MX &x, const MX &y, const MX &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const MX &e, const std::vector< MX > &boundary, const Dict &options) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const std::vector< MX > &e, const std::vector< MX > &boundary, const Dict &options) {
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
%exception  casadi::MX::mrdivide(const MX &a, const MX &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::mtimes(const MX &x, const MX &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::n_nodes(const MX &x) {
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
%exception  casadi::MX::pinv(const MX &A, const std::string &lsolver="qr", const Dict &dict=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::polyval(const MX &p, const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::print_operator(const MX &x, const std::vector< std::string > &args) {
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
%exception  casadi::MX::reset_input() const {
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
%exception  casadi::MX::set_temp(casadi_int t) const {
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
%exception  casadi::Map::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_default_in(casadi_int ind) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_function(const std::string &name) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::has_function(const std::string &fname) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::has_spfwd() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::has_sprev() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::info() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::is_a(const std::string &type, bool recursive) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::parallelization() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Map::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_default_in(casadi_int ind) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_function(const std::string &name) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::has_function(const std::string &fname) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::has_spfwd() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::has_sprev() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::parallelization() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MapSum::serialize_type(SerializingStream &s) const override {
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
%exception  casadi::Matrix< Scalar >::adj(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::all(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::any(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::binary(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::blockcat(const std::vector< std::vector< Matrix< Scalar > > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::chol(const Matrix< Scalar > &A) {
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
%exception  casadi::Matrix< Scalar >::forward(const std::vector< Matrix< Scalar > > &ex, const std::vector< Matrix< Scalar > > &arg, const std::vector< std::vector< Matrix< Scalar > > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order, const Matrix< Scalar > &w) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::gauss_quadrature(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Matrix< Scalar > &a, const Matrix< Scalar > &b, casadi_int order=5) {
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
%exception  casadi::Matrix< Scalar >::inv(const Matrix< Scalar > &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inv(const Matrix< Scalar > &A, const std::string &lsolver, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::inv_minor(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::is_equal(const Matrix< Scalar > &x, const Matrix< Scalar > &y, casadi_int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::jacobian(const Matrix< Scalar > &f, const Matrix< Scalar > &x, const Dict &opts=Dict()) {
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
%exception  casadi::Matrix< Scalar >::matrix_matrix(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::matrix_scalar(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
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
%exception  casadi::Matrix< Scalar >::print_operator(const Matrix< Scalar > &x, const std::vector< std::string > &args) {
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
%exception  casadi::Matrix< Scalar >::rand(const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::rectangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reshape(const Matrix< Scalar > &x, casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reshape(const Matrix< Scalar > &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::reverse(const std::vector< Matrix< Scalar > > &ex, const std::vector< Matrix< Scalar > > &arg, const std::vector< std::vector< Matrix< Scalar > > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::scalar() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::scalar_matrix(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::serialize(std::ostream &stream) const {
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
%exception  casadi::Matrix< Scalar >::trace(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::triangle(const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::unary(casadi_int op, const Matrix< Scalar > &x) {
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
%exception  casadi::Newton::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Newton::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::disp_more(std::ostream &stream) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::getReducedHessian() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_default_in(casadi_int ind) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::integer_support() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::is_a(const std::string &type, bool recursive) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::kkt() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::nlpsol_codegen_body(CodeGenerator &g) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::uses_output() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OmpMap::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OmpMap::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OmpMap::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OmpMap::is_a(const std::string &type, bool recursive) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OmpMap::parallelization() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::create_function(const std::string &fname, const std::vector< std::string > &s_in, const std::vector< std::string > &s_out, const Function::AuxOut &aux=Function::AuxOut()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::expand() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::finalize() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::generate_dependencies(const std::string &fname, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::get_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::get_function(const std::string &name) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::has_function(const std::string &fname) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::jit_dependencies(const std::string &fname) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::monitored(const std::string &name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::oracle() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::set_function(const Function &fcn) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::set_function(const Function &fcn, const std::string &fname, bool jit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Conic  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Dple  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Dple  >::serialize_type(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Expm  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Expm  >::serialize_type(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< ImporterInternal  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Integrator  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Interpolant  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< LinsolInternal  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Nlpsol  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Rootfinder  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Simulator  >::plugin_name() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Simulator  >::serialize_type(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Printable::repr(const Derived &obj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PrintableObject::repr(const PrintableObject< Derived > &obj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PrintableObject::str(const PrintableObject< Derived > &obj) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::alloc_mem() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::change_option(const std::string &option_name, const GenericType &option_value) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::checkout() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::clear_mem() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::construct(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::finalize() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::generate_options(bool is_temp=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::get_options() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::has_option(const std::string &option_name) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::memory(int ind) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::print_option(const std::string &name, std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::print_options(std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::print_time(const std::map< std::string, FStats > &fstats) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::release(int mem) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::serialize(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::serialize_base_function() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::QpToNlp::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrqp::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrqp::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrqp::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrqp::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrqp::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrqp::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrqp::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::getConic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::print_iteration() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf, double dx_norm, double rg, casadi_int ls_trials, bool ls_success) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::ad_forward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::ad_reverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool always_inline, bool never_inline) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::has_spfwd() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::has_sprev() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::serialize_base_function() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::uses_output() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RungeKutta::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RungeKutta::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RungeKutta::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RungeKutta::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::RungeKutta::setupFG() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::has_duplicates() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::reset_input() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::assignIfDuplicate(const SXElem &scalar, casadi_int depth=1) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::assignNoDelete(const SXElem &scalar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::dep(casadi_int ch=0) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::disp(std::ostream &stream, bool more=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::get() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::get_temp() const {
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
%exception  casadi::SXElem::is_regular() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_symbolic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::is_zero() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::mark() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SXElem::marked() const {
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
%exception  casadi::SXElem::set_temp(int t) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::getConic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::get() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::getCount() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::print_ptr(std::ostream &stream=casadi::uout()) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::swap(SharedObject &other) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObject::weak() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectInternal::class_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectInternal::disp(std::ostream &stream, bool more) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectInternal::getCount() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectInternal::weak() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ShellCompiler::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ShellCompiler::get_function(const std::string &symname) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ShellCompiler::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ShellCompiler::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ShellCompiler::library() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ShellCompiler::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::create_advanced(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::get_name_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::get_name_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::p() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::sp_jac_dae() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::t() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::u() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::x() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::y() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Simulator::z() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotC2d::alloc_memory() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotC2d::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotC2d::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotC2d::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotDple::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotDple::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotDple::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotDple::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotDple::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotExpm::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotExpm::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotExpm::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SlicotExpm::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::calc_fd() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::calc_stepsize(double abstol) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::get_abstol() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::has_err() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::n_pert() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::pert(casadi_int k, double h) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Smoothing::pert(const std::string &k) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::blockcat(const std::vector< std::vector< Sparsity > > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::colind() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::combine(const Sparsity &y, bool f0x_is_zero, bool function0_is_zero, std::vector< unsigned char > &mapping) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::diagcat(const std::vector< Sparsity > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::diagsplit(const Sparsity &x, const std::vector< casadi_int > &offset1, const std::vector< casadi_int > &offset2) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::file_format(const std::string &filename, const std::string &format_hint, const std::set< std::string > &file_formats) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::find(bool ind1=SWIG_IND1) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::find(std::vector< casadi_int > &loc, bool ind1=false) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::get() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::getCache() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::getEmpty() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::getScalar() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::getScalarSparse() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::horzcat(const std::vector< Sparsity > &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::horzsplit(const Sparsity &x, const std::vector< casadi_int > &offset) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::intersect(const Sparsity &y, std::vector< unsigned char > &mapping) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::kron(const Sparsity &a, const Sparsity &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::mac(const Sparsity &x, const Sparsity &y, const Sparsity &z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::mtimes(const Sparsity &x, const Sparsity &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::norm_0_mul(const Sparsity &x, const Sparsity &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::reshape(const Sparsity &x, casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::reshape(const Sparsity &x, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::row() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::serialize(std::ostream &stream) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::sprank(const Sparsity &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::sum1(const Sparsity &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::sum2(const Sparsity &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::tril(const Sparsity &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::triu(const Sparsity &x, bool includeDiagonal=true) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::unite(const Sparsity &y, std::vector< unsigned char > &mapping) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::vertcat(const std::vector< Sparsity > &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sparsity::vertsplit(const Sparsity &x, const std::vector< casadi_int > &offset) {
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
%exception  casadi::Sqpmethod::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_calc_gamma_1(CodeGenerator &cg) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_qp_ela_solve(CodeGenerator &cg, const std::string &H, const std::string &g, const std::string &lbdz, const std::string &ubdz, const std::string &A, const std::string &x_opt, const std::string &dlam) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_qp_solve(CodeGenerator &cg, const std::string &H, const std::string &g, const std::string &lbdz, const std::string &ubdz, const std::string &A, const std::string &x_opt, const std::string &dlam, int mode) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_solve_elastic_mode(CodeGenerator &cg, int mode) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::getConic() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::print_iteration() const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf, double dx_norm, double rg, casadi_int ls_trials, bool ls_success, bool so_succes, std::string info) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsSimulator::get_abstol() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsSimulator::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsSimulator::get_reltol() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SundialsSimulator::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::disp_more(std::ostream &stream) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::get_n_in() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::get_n_out() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::has_forward(casadi_int nfwd) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::info() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Switch::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SymbolicQr::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ThreadMap::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ThreadMap::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ThreadMap::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ThreadMap::is_a(const std::string &type, bool recursive) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ThreadMap::parallelization() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::XmlFile::parse(const std::string &filename) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::assign_vector(const std::vector< S > &s, std::vector< D > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::casadi_max(casadi_int x, casadi_int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::casadi_min(casadi_int x, casadi_int y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::check_exposed(T t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::constpow(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::copy_vector(const std::vector< S > &s, std::vector< D > &d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::copysign(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::cumsum0(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::diff(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::dot(const std::vector< T > &a, const std::vector< T > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::enum_names() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::erfinv(double x) {
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
%exception  casadi::has_enum(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::if_else(double x, double y, double z) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::if_else_zero(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::init_vector(std::vector< S > &d, const std::vector< D > &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::integrator(const std::string &name, const std::string &solver, const Function &dae, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::isUnique(const std::vector< T > &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_equal(double x, double y, casadi_int depth=0) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_regular(N_Vector v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::join(const std::vector< T > &a, const std::vector< T > &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::join(const std::vector< T > &a, const std::vector< T > &b, const std::vector< T > &c) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::linspace(std::vector< T > &v, const F &first, const L &last) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::nlpsol(const std::string &name, const std::string &solver, const Function &nlp, const Dict &opts=Dict()) {
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
%exception  casadi::operation_checker(casadi_int op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::operation_getter(casadi_int op) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::permute(const std::vector< T > &a, const std::vector< casadi_int > &order) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::printme(double x, double y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::product(const std::vector< T > &values) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::replace_mat(const M &arg, const Sparsity &inp, casadi_int npar) {
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
%exception  casadi::sign(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::simplify(double x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sort(const std::vector< T > &values, std::vector< T > &sorted_values, std::vector< casadi_int > &indices, bool invert_indices=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::sq(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const T &v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::str(const T &v, bool more) {
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
%exception  casadi::text2set(const std::string &text) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::text2type(const std::string &text) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::text2vector(const std::string &text) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::to_enum(const std::string &s, const std::string &s_def="") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::twice(const T &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_slice(const std::vector< T > &v, const std::vector< casadi_int > &i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::vector_static_cast(const std::vector< S > &rhs) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::zip(const std::vector< std::string > &id, const std::vector< T > &mat) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  std::mul_overflows(const T &a, const T &b) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::BSplineInterpolant::BSplineInterpolant(const std::string &name, const std::vector< double > &grid, const std::vector< casadi_int > &offset, const std::vector< double > &values, casadi_int m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::BackwardDiff::BackwardDiff(const std::string &name, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Blocksqp::Blocksqp(const std::string &name, const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CentralDiff::CentralDiff(const std::string &name, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ClangCompiler::ClangCompiler(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ClangXCompiler::ClangXCompiler(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Collocation::Collocation(const std::string &name, const Function &dae) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Conic::Conic(const std::string &name, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CustomNlpsol::CustomNlpsol(const std::string &name, const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::CvodesSimulator::CvodesSimulator(const std::string &name, const Function &dae, const std::vector< double > &grid) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DeserializerBase::DeserializerBase(std::unique_ptr< std::istream > stream) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::DllLibrary::DllLibrary(const std::string &bin_name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Dple::Dple(const std::string &name, const SpDict &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Expm::Expm(const std::string &name, const Sparsity &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::External::External(const std::string &name, const Importer &li) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FStats::FStats() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FastNewton::FastNewton(const std::string &name, const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FiniteDiff::FiniteDiff(const std::string &name, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::FixedStepIntegrator::FixedStepIntegrator(const std::string &name, const Function &dae) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ForwardDiff::ForwardDiff(const std::string &name, casadi_int n) {
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
%exception casadi::FunctionBuffer::FunctionBuffer(const FunctionBuffer &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericExternal::GenericExternal(DeserializingStream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericExternal::GenericExternal(const std::string &name, const Importer &li) {
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
%exception casadi::GenericType::GenericType(const std::vector< Function > &f) {
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
%exception casadi::GenericType::GenericType(const std::vector< std::vector< casadi_int > > &ivv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(const std::vector< std::vector< double > > &dv) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(double d) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::GenericType::GenericType(int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(const std::string &name, const Function &dae) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ImplicitToNlp::ImplicitToNlp(const std::string &name, const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Integrator::Integrator(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Interpolant::Interpolant(const std::string &name, const std::vector< double > &grid, const std::vector< casadi_int > &offset, const std::vector< double > &values, casadi_int m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Ipqp::Ipqp(const std::string &name, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::JitFunction::JitFunction(const std::string &name, const std::string &body, const std::vector< std::string > &name_in, const std::vector< std::string > &name_out, const std::vector< Sparsity > &sparsity_in, const std::vector< Sparsity > &sparsity_out) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LapackLu::LapackLu(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LapackQr::LapackQr(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinearInterpolant::LinearInterpolant(const std::string &name, const std::vector< double > &grid, const std::vector< casadi_int > &offset, const std::vector< double > &values, casadi_int m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinearInterpolantJac::LinearInterpolantJac(DeserializingStream &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinearInterpolantJac::LinearInterpolantJac(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolLdl::LinsolLdl(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolQr::LinsolQr(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolTridiag::LinsolTridiag(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Lsqr::Lsqr(const std::string &name, const Sparsity &sp) {
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
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const Scalar &val, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const Sparsity &sp, const std::vector< Scalar > &d, bool dummy) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(const std::pair< casadi_int, casadi_int > &rc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Matrix< Scalar >::Matrix(std::initializer_list< Scalar > x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Newton::Newton(const std::string &name, const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::NlImporter::NlImporter(NlpBuilder &nlp, const std::string &filename, const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Nlpsol::Nlpsol(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OmpMap::OmpMap(const std::string &name, const Function &f, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Opti::Opti(const Opti &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OracleFunction::OracleFunction(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToNlp::QpToNlp(const std::string &name, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Qrqp::Qrqp(const std::string &name, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Qrsqp::Qrsqp(const std::string &name, const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Rootfinder::Rootfinder(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RungeKutta::RungeKutta(const std::string &name, const Function &dae) {
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
%exception casadi::ScopedTiming::ScopedTiming(FStats &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Scpgen::Scpgen(const std::string &name, const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SerializerBase::SerializerBase(std::unique_ptr< std::ostream > stream, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SharedObject::SharedObject(const SharedObject &ref) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ShellCompiler::ShellCompiler(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Simulator::Simulator(const std::string &name, const Function &oracle, const std::vector< double > &grid) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SlicotC2d::SlicotC2d() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SlicotC2d::SlicotC2d(const std::string &name, const DMDict &AB) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SlicotDple::SlicotDple() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SlicotDple::SlicotDple(const SpDict &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SlicotDple::SlicotDple(const std::string &name, const SpDict &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SlicotExpm::SlicotExpm() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SlicotExpm::SlicotExpm(const std::string &name, const Sparsity &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Smoothing::Smoothing(const std::string &name, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Sqpmethod::Sqpmethod(const std::string &name, const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SundialsSimulator::SundialsSimulator(const std::string &name, const Function &dae, const std::vector< double > &grid) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Switch::Switch(const std::string &name, const std::vector< Function > &f, const Function &f_def) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicQr::SymbolicQr(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ThreadMap::ThreadMap(const std::string &name, const Function &f, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::null_ptr_on_copy< T >::null_ptr_on_copy() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::null_ptr_on_copy< T >::null_ptr_on_copy(const null_ptr_on_copy< T > &rhs) {
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