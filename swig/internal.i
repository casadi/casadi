%include "windows.i"
%exception  casadi::BSplineInterpolant::class_name() const override {

 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::BSplineInterpolant::codegen_declarations(CodeGenerator &g) const override {
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
%exception  casadi::BSplineInterpolant::is_diff_in(casadi_int i) override {
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
%exception  casadi::Conic::is_a(const std::string &type, bool recursive) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Conic::psd_support() const  {
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
%exception  casadi::External::any_symbol_found() const  {
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
%exception  casadi::FStats::reset() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::tic() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FStats::toc() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_input(const std::string &s, const MatType &e, bool is_diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::add_output(const std::string &s, const MatType &e, bool is_diff) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Factory< MatType >::calculate(const Dict &opts=Dict()) {
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
%exception  casadi::FunctionInternal::call_gen(const MXVector &arg, MXVector &res, casadi_int npar, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_gen(const std::vector< Matrix< D > > &arg, std::vector< Matrix< D > > &res, casadi_int npar, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_reverse(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &aseed, std::vector< std::vector< MX > > &asens, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::call_reverse(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &aseed, std::vector< std::vector< SX > > &asens, bool always_inline, bool never_inline) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_arg(const std::vector< M > &arg, casadi_int &npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_mat(const Sparsity &arg, const Sparsity &inp, casadi_int &npar) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::check_res(const std::vector< M > &res, casadi_int &npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen(CodeGenerator &g, const std::string &fname) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_alloc_mem(CodeGenerator &g) const  {
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
%exception  casadi::FunctionInternal::codegen_free_mem(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_incref(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_init_mem(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_mem(CodeGenerator &g, const std::string &index="mem") const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_mem_type() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_meta(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_name(const CodeGenerator &g, bool ns=true) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::codegen_sparsities(CodeGenerator &g) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_arg(const std::map< std::string, M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_arg(const std::vector< M > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_res(const std::map< std::string, M > &res) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::convert_res(const std::vector< M > &res) const  {
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
%exception  casadi::FunctionInternal::dm_in() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dm_in(casadi_int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dm_out() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::dm_out(casadi_int ind) const  {
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
%exception  casadi::FunctionInternal::finalize() override {
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
%exception  casadi::FunctionInternal::generate_options(bool is_temp=false) const override {
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
%exception  casadi::FunctionInternal::has_eval_dm() const  {
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
%exception  casadi::FunctionInternal::has_jacobian_sparsity() const  {
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
%exception  casadi::FunctionInternal::incache(const std::string &fname, Function &f, const std::string &suffix="") const  {
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
%exception  casadi::FunctionInternal::instructions_sx() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::is_a(const std::string &type, bool recursive) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::is_diff_in(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::is_diff_out(casadi_int i) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jac() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacobian() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacobian_sparsity() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jacobian_sparsity_filter(const Sparsity &sp) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::jit_dependencies(const std::string &fname) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::map(casadi_int n, const std::string &parallelization) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::mapsum_mx(const std::vector< MX > &arg, const std::string &parallelization) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::matching_arg(const std::vector< M > &arg, casadi_int &npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::matching_res(const std::vector< M > &arg, casadi_int &npar) const  {
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
%exception  casadi::FunctionInternal::nz_in(const std::vector< DM > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nz_in(const std::vector< double > &arg) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nz_out(const std::vector< DM > &res) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::nz_out(const std::vector< double > &res) const  {
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
%exception  casadi::FunctionInternal::project_arg(const std::vector< M > &arg, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::project_res(const std::vector< M > &arg, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_arg(const std::vector< M > &arg, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_aseed(const std::vector< std::vector< M > > &aseed, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_aseed(const std::vector< std::vector< M >> &aseed, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_fseed(const std::vector< std::vector< M > > &fseed, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_fseed(const std::vector< std::vector< M >> &fseed, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::replace_res(const std::vector< M > &res, casadi_int npar) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::reverse(casadi_int nadj) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::self() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::serialize_body(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::FunctionInternal::serialize_type(SerializingStream &s) const override {
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
%exception  casadi::FunctionInternal::tocache(const Function &f, const std::string &suffix="") const  {
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
%exception  casadi::FunctionInternal::wrap_as_needed(const Dict &opts) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::any_symbol_found() const override {
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
%exception  casadi::GenericExternal::init_external() override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericExternal::serialize_type(SerializingStream &s) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::conditional(const MatType &ind, const std::vector< MatType > &x, const MatType &x_default, bool short_circuit=false) {
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
%exception  casadi::GenericMatrix::inv(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::inv(const MatType &A, const std::string &lsolver, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::inv_minor(const MatType &A) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::jacobian(const MatType &ex, const MatType &arg, const Dict &opts=Dict()) {
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
%exception  casadi::GenericMatrix::reverse(const std::vector< MatType > &ex, const std::vector< MatType > &arg, const std::vector< std::vector< MatType > > &v, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::shared(const std::vector< MatType > &ex, std::vector< MatType > &ex_output, std::vector< MatType > &v, std::vector< MatType > &vdef, const std::string &v_prefix="v_", const std::string &v_suffix="") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::shared(std::vector< MatType > &ex, std::vector< MatType > &v, std::vector< MatType > &vdef, const std::string &v_prefix="v_", const std::string &v_suffix="") {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::simplify(const MatType &x) {
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
%exception  casadi::GenericMatrix::unite(const MatType &A, const MatType &B) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::GenericMatrix::which_depends(const MatType &expr, const MatType &var, casadi_int order, bool tr) {
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
%exception  casadi::ImporterInternal::finalize() {
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
%exception  casadi::ImporterInternal::library() const  {
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
%exception  casadi::ImporterInternal::serialize(SerializingStream &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::serialize_body(SerializingStream &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::serialize_type(SerializingStream &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ImporterInternal::to_text(const std::string &cmd, casadi_int ind=-1) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::algebraic_state_init(const MX &x0, const MX &z0) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::algebraic_state_output(const MX &Z) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::aug_adj(casadi_int nadj) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::aug_fwd(casadi_int nfwd) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::create_advanced(const Dict &opts) {
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
%exception  casadi::Integrator::t() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::x() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Integrator::z() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::arg_grid() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::arg_values() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::class_name() const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::coeff_size() const  {
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
%exception  casadi::Interpolant::has_parametric_grid() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::has_parametric_values() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::init(const Dict &opts) override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Interpolant::is_diff_in(casadi_int i) override {
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
%exception  casadi::LinearInterpolantJac::has_parametric_grid() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::LinearInterpolantJac::has_parametric_values() const  {
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
%exception  casadi::MX::bspline(const MX &x, const DM &coeffs, const std::vector< std::vector< double > > &knots, const std::vector< casadi_int > &degree, casadi_int m, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::bspline(const MX &x, const MX &coeffs, const std::vector< std::vector< double > > &knots, const std::vector< casadi_int > &degree, casadi_int m, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::convexify(const MX &H, const Dict &opts=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::evalf(const MX &expr) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::find(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::get_temp() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const MX &ex, const std::vector< MX > &v, const std::vector< MX > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::graph_substitute(const std::vector< MX > &ex, const std::vector< MX > &v, const std::vector< MX > &vdef) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::has_duplicates() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::inv_node(const MX &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::lift(const MX &x, const MX &x_guess) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::low(const MX &v, const MX &p, const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const MX &e, const std::vector< MX > &boundary=std::vector< MX >(), const Dict &options=Dict()) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::MX::matrix_expand(const std::vector< MX > &e, const std::vector< MX > &boundary=std::vector< MX >(), const Dict &options=Dict()) {
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
%exception  casadi::Map::parallelization() const  {
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
%exception  casadi::MapSum::parallelization() const  {
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
%exception  casadi::Matrix::taylor(const Matrix< Scalar > &ex, const Matrix< Scalar > &x) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::binary(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::matrix_matrix(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::matrix_scalar(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::scalar_matrix(casadi_int op, const Matrix< Scalar > &x, const Matrix< Scalar > &y) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Matrix< Scalar >::unary(casadi_int op, const Matrix< Scalar > &x) {
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
%exception  casadi::Nlpsol::integer_support() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::is_a(const std::string &type, bool recursive) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::kkt() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Nlpsol::nlpsol_codegen_body(CodeGenerator &g) const  {
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
%exception  casadi::OracleFunction::monitored(const std::string &name) const  {
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
%exception  casadi::PluginInterface< Conic  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Dple  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Dple  >::serialize_type(SerializingStream &s) const {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Expm  >::plugin_name() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::PluginInterface< Expm  >::serialize_type(SerializingStream &s) const {
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
%exception  casadi::ProtoFunction::finalize() {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::generate_options(bool is_temp=false) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::get_options() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::memory(int ind) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::print_time(const std::map< std::string, FStats > &fstats) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::release(int mem) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::serialize(SerializingStream &s) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::ProtoFunction::serialize_base_function() const  {
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
%exception  casadi::Qrsqp::getConic() const  {
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
%exception  casadi::Qrsqp::print_iteration() const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Qrsqp::print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf, double dx_norm, double rg, casadi_int ls_trials, bool ls_success) const  {
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
%exception  casadi::SX::has_duplicates() const  {
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
%exception  casadi::ShellCompiler::library() const override {
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
%exception  casadi::Sqpmethod::codegen_body(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_declarations(CodeGenerator &g) const override {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::codegen_qp_solve(CodeGenerator &cg, const std::string &H, const std::string &g, const std::string &lbdz, const std::string &ubdz, const std::string &A, const std::string &x_opt, const std::string &dlam) const  {
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
%exception  casadi::Sqpmethod::print_iteration(casadi_int iter, double obj, double pr_inf, double du_inf, double dx_norm, double rg, casadi_int ls_trials, bool ls_success) const  {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::Sqpmethod::serialize_body(SerializingStream &s) const override {
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
%exception  casadi::casadi_cvx_scalar(T1 epsilon, casadi_int reflect, T1 eig) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::check_exposed(T t) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::is_regular(N_Vector v) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception  casadi::replace_mat(const M &arg, const Sparsity &inp, casadi_int npar) {
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
%exception casadi::GenericExternal::GenericExternal(DeserializingStream &s) {
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
%exception casadi::Newton::Newton(const std::string &name, const Function &f) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::Nlpsol::Nlpsol(const std::string &name, const Function &oracle) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}
%exception casadi::OmpMap::OmpMap(const std::string &name, const Function &f, casadi_int n) {
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
%exception casadi::ScopedTiming::ScopedTiming(FStats &f) {
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
%exception casadi::ThreadMap::ThreadMap(const std::string &name, const Function &f, casadi_int n) {
 CATCH_OR_NOT(INTERNAL_MSG() $action) 
}