%exception  casadi::BSplineInterpolant::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::has_codegen() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::has_jacobian() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::plugin_name() const override {
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
%exception  casadi::Collocation::setupFG() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::generateNativeCode(std::ostream &file) const  {
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
%exception  casadi::Conic::integer_support() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConicActiveSet::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConicActiveSet::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConicActiveSet::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConicActiveSet::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ConicActiveSet::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::can_have_meta() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::DllLibrary::get_function(const std::string &symname) override {
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
%exception  casadi::Expm::getJacSparsity(casadi_int iind, casadi_int oind, bool symmetric) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const override {
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
%exception  casadi::Expm::has_reverse(casadi_int nadj) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Expm::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::External::factory(const std::string &name, const std::vector< std::string > &s_in, const std::vector< std::string > &s_out, const Function::AuxOut &aux, const Dict &opts) const override {
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
%exception  casadi::FStats::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::tic() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::toc() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_input(const std::string &s, const MatType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_output(const std::string &s, const MatType &e) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::get_input(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::get_output(const std::string &s) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::has_in(const std::string &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::has_out(const std::string &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::name_in() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::name_out() const  {
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
%exception  casadi::FixedStepIntegrator::getExplicit() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::getExplicitB() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FixedStepIntegrator::init(const Dict &opts) override {
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
%exception  casadi::FunctionInternal::ad_weight() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::adjViaJac(casadi_int nadj) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::all_scalar() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::alloc(const Function &f, bool persistent=false) {
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
%exception  casadi::FunctionInternal::call(const std::vector< M > &arg, std::vector< M > &res, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_forward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_forward(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &fsens, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_gen(const MXVector &arg, MXVector &res, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_gen(const std::vector< Matrix< D > > &arg, std::vector< Matrix< D > > &res, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_reverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_reverse(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_arg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_res(const std::vector< M > &res) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen(CodeGenerator &g, const std::string &fname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_body(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_declarations(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_decref(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_incref(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_meta(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_name(const CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_sparsities(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::definition() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::disp(std::ostream &stream, bool more) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::disp_more(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::eval_dm(const std::vector< DM > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::eval_mx(const MXVector &arg, MXVector &res, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::export_code(const std::string &lang, std::ostream &stream, const Dict &options) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::factory(const std::string &name, const std::vector< std::string > &s_in, const std::vector< std::string > &s_out, const Function::AuxOut &aux, const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::finalize(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::forward(casadi_int nfwd) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::free_mx() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::free_sx() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::fwdViaJac(casadi_int nfwd) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::fwd_seed(casadi_int nfwd) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generate_dependencies(const std::string &fname, const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::generate_lifted(Function &vdef_fcn, Function &vinit_fcn) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getAdaptorSolverName() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsity(casadi_int iind, casadi_int oind, bool symmetric) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityGen(casadi_int iind, casadi_int oind, bool symmetric, casadi_int gr_i=1, casadi_int gr_o=1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchical(casadi_int iind, casadi_int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::getJacSparsityHierarchicalSymm(casadi_int iind, casadi_int oind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_abstol() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_default_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_forward(casadi_int nfwd, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_free() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_function() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_function(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jac(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jacobian(const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_jacobian_sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_max_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_min_in(casadi_int ind) const  {
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
%exception  casadi::FunctionInternal::get_options() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_partition(casadi_int iind, casadi_int oind, Sparsity &D1, Sparsity &D2, bool compact, bool symmetric, bool allow_forward, bool allow_reverse) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_reltol() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_reverse(casadi_int nadj, const std::string &name, const std::vector< std::string > &inames, const std::vector< std::string > &onames, const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_sparsity_in(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::get_sparsity_out(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_codegen() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_derivative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_forward(casadi_int nfwd) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_free() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_function(const std::string &fname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_jac() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_jacobian() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_reverse(casadi_int nadj) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_spfwd() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::has_sprev() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::index_in(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::index_out(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::info() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_MX(casadi_int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_constant(casadi_int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_id(casadi_int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_input(casadi_int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::instruction_output(casadi_int k) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::is_a(const std::string &type, bool recursive) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jac() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacobian() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jit_dependencies(const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mapsum_mx(const std::vector< MX > &arg, const std::string &parallelization) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::matching_arg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::matching_res(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_in() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_out() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mx_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::n_instructions() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::n_nodes() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_in() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_out() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nnz_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_in() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_out() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::numel_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::oracle() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::print_dimensions(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::print_option(const std::string &name, std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::print_options(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_arg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_aseed(const std::vector< std::vector< M > > &aseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_fseed(const std::vector< std::vector< M > > &fseed) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_res(const std::vector< M > &res) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::reverse(casadi_int nadj) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::self() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::serialize(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::serialize_header(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::signature(const std::string &fname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size1_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size1_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size2_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size2_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::size_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::slice(const std::string &name, const std::vector< casadi_int > &order_in, const std::vector< casadi_int > &order_out, const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sp_weight() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sparsity_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sparsity_jac(casadi_int iind, casadi_int oind, bool compact, bool symmetric) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sparsity_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_in() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_out() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sx_out(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolicAdjSeed(casadi_int nadj, const std::vector< MatType > &v) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::symbolic_output(const std::vector< MX > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_arg() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_iw() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_res() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_w() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::sz_work(size_t &sz_arg, size_t &sz_res, size_t &sz_iw, size_t &sz_w) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::uses_output() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::which_depends(const std::string &s_in, const std::vector< std::string > &s_out, casadi_int order, bool tr=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::wrap() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::get_sparsity_in(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::get_sparsity_out(casadi_int i) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericTypeBase::getType() const  {
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
%exception  casadi::ImporterInternal::body(const std::string &symname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::can_have_meta() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::construct(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::disp(std::ostream &stream, bool more) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::get_meta(const std::string &cmd, casadi_int ind=-1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::get_options() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::has_function(const std::string &symname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::has_meta(const std::string &cmd, casadi_int ind=-1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::init(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::inlined(const std::string &symname) const  {
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
%exception  casadi::ImporterInternal::to_text(const std::string &cmd, casadi_int ind=-1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::aug_adj(casadi_int nadj) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::aug_fwd(casadi_int nfwd) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::getDerivativeOptions(bool fwd) const  {
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
%exception  casadi::Integrator::p() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::q() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rq() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rx() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::rz() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::sp_jac_dae() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::sp_jac_rdae() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::t() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::x() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::z() const  {
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
%exception  casadi::Interpolant::init(const Dict &opts) override {
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
%exception  casadi::LinearInterpolantJac::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::colind() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::disp(std::ostream &stream, bool more) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::disp_more(std::ostream &stream) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::generate(CodeGenerator &g, const std::string &A, const std::string &x, casadi_int nrhs, bool tr) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::ncol() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::nnz() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::nrow() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolInternal::row() const  {
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
%exception  casadi::LinsolLdl::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolLdl::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::generate(CodeGenerator &g, const std::string &A, const std::string &x, casadi_int nrhs, bool tr) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinsolQr::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Lsqr::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Lsqr::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Lsqr::plugin_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_temp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::has_duplicates() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::reset_input() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::set_temp(casadi_int t) const  {
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
%exception  casadi::Map::parallelization() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::binary(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::matrix_matrix(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::matrix_scalar(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::scalar_matrix(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< T >::unary(casadi_int op, const Matrix< Scalar > &x) {
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
%exception  casadi::Nlpsol::alloc_mem() const override {
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
%exception  casadi::Nlpsol::integer_support() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::kkt() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::setOptionsFromFile(const std::string &file) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::uses_output() const override {
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
%exception  casadi::OracleFunction::finalize(const Dict &opts) override {
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
%exception  casadi::OracleFunction::monitored(const std::string &name) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::oracle() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::set_function(const Function &fcn) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::OracleFunction::set_function(const Function &fcn, const std::string &fname, bool jit=false) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Conic  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Dple  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Expm  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Integrator  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Interpolant  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Nlpsol  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Rootfinder  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::alloc_mem() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::checkout() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::clear_mem() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::construct(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::finalize(const Dict &opts) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::get_options() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::memory(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::release(casadi_int mem) const  {
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
%exception  casadi::Rootfinder::ad_forward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &fsens, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Rootfinder::ad_reverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool always_inline, bool never_inline) const  {
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
%exception  casadi::RungeKutta::setupFG() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::dep(casadi_int ch=0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::element_hash() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::has_duplicates() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_commutative() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_leaf() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_op(casadi_int op) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_regular() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_smooth() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_symbolic() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::is_valid_input() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::n_dep() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::op() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::print_split(std::vector< std::string > &output_nz, std::vector< std::string > &output_inter) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SX::reset_input() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Scpgen::getConic() const  {
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
%exception  casadi::SharedObject::print_ptr(std::ostream &stream=casadi::uout()) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectInternal::class_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectInternal::disp(std::ostream &stream, bool more) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SharedObjectInternal::getCount() const  {
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
%exception  casadi::ShellCompiler::plugin_name() const override {
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
%exception  casadi::SparseStorage< DataType >::clear() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::elem(casadi_int rr, casadi_int cc) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::has_nz(casadi_int rr, casadi_int cc) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::nonzeros() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::nonzeros() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::reserve(casadi_int nnz) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::reserve(casadi_int nnz, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::resize(casadi_int nrow, casadi_int ncol) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::SparseStorage< DataType >::sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::alloc_mem() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::getConic() const  {
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
%exception  casadi::Sqpmethod::print_iteration() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf, double dx_norm, double reg, casadi_int ls_trials, bool ls_success) const  {
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
%exception  casadi::check_exposed(T t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::combine(const Dict &first, const Dict &second) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_regular(N_Vector v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::matrixName< SXElem >() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::replace_mat(const M &arg, const Sparsity &inp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::zip(const std::vector< std::string > &id, const std::vector< T > &mat) {
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
%exception casadi::Collocation::Collocation(const std::string &name, const Function &dae) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Conic::Conic(const std::string &name, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ConicActiveSet::ConicActiveSet(const std::string &name, const std::map< std::string, Sparsity > &st) {
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
%exception casadi::Factory< MatType >::Factory(const Function::AuxOut &aux) {
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
%exception casadi::GenericExternal::GenericExternal(const std::string &name, const Importer &li) {
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
%exception casadi::LinearInterpolantJac::LinearInterpolantJac(const std::string &name) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolLdl::LinsolLdl(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::LinsolQr::LinsolQr(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Lsqr::Lsqr(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Newton::Newton(const std::string &name, const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Nlpsol::Nlpsol(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OracleFunction::OracleFunction(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::QpToNlp::QpToNlp(const std::string &name, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Rootfinder::Rootfinder(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::RungeKutta::RungeKutta(const std::string &name, const Function &dae) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Scpgen::Scpgen(const std::string &name, const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::ShellCompiler::ShellCompiler(const std::string &name) {
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
%exception casadi::SparseStorage< DataType >::SparseStorage() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const SparseStorage< DataType > &m) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SparseStorage< DataType >::SparseStorage(const Sparsity &sparsity, const DataType &val=DataType(0)) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Sqpmethod::Sqpmethod(const std::string &name, const Function &nlp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Switch::Switch(const std::string &name, const std::vector< Function > &f, const Function &f_def) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::SymbolicQr::SymbolicQr(const std::string &name, const Sparsity &sp) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}