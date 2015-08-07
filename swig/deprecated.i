%exception  casadi::GenericMatrix< MX  >::size() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< MX  >::sparse(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< MX  >::sparse(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< MatType >::size() const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< MatType >::sparse(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< MatType >::sparse(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::size() const {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::sparse(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::sparse(int nrow=1, int ncol=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Sparsity::reserve(int nnz, int ncol) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Sparsity::sparse(const std::pair< int, int > &rc) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Sparsity::sparse(int nrow, int ncol=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::CleSolver::CleSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::ControlSimulator::ControlSimulator(const Function &dae, const Function &output_fcn, const Matrix< double > &grid) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::ControlSimulator::ControlSimulator(const Function &dae, const Matrix< double > &grid) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::CustomFunction::CustomFunction(const CustomEvaluate &c_fcn, const std::pair< SparsityDict, std::vector< std::string > > &inputscheme, const std::pair< SparsityDict, std::vector< std::string > > &outputscheme) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::CustomFunction::CustomFunction(const CustomEvaluate &c_fcn, const std::pair< SparsityDict, std::vector< std::string > > &inputscheme, const std::vector< Sparsity > &outputscheme) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::CustomFunction::CustomFunction(const CustomEvaluate &c_fcn, const std::vector< Sparsity > &inputscheme, const std::pair< SparsityDict, std::vector< std::string > > &outputscheme) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::CustomFunction::CustomFunction(const CustomEvaluate &c_fcn, const std::vector< Sparsity > &inputscheme, const std::vector< Sparsity > &outputscheme) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::DleSolver::DleSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::DpleSolver::DpleSolver(const std::string &solver, const std::map< std::string, std::vector< Sparsity > > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::ImplicitFunction::ImplicitFunction(const std::string &solver, const Function &f, const Function &jac=Function(), const LinearSolver &linsol=LinearSolver()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::Integrator::Integrator(const std::string &solver, const Function &f, const Function &g=Function()) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::LinearSolver::LinearSolver(const std::string &solver, const Sparsity &sp, int nrhs=1) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::LpSolver::LpSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::LrDleSolver::LrDleSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::LrDpleSolver::LrDpleSolver(const std::string &solver, const std::map< std::string, std::vector< Sparsity > > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::MXFunction::MXFunction(const std::pair< MXDict, std::vector< std::string > > &arg, const std::pair< MXDict, std::vector< std::string > > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::MXFunction::MXFunction(const std::pair< MXDict, std::vector< std::string > > &arg, const std::vector< MX > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::MXFunction::MXFunction(const std::vector< MX > &arg, const std::pair< MXDict, std::vector< std::string > > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::MXFunction::MXFunction(const std::vector< MX > &arg, const std::vector< MX > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::NlpSolver::NlpSolver(const std::string &solver, const Function &nlp) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::QcqpSolver::QcqpSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::QpSolver::QpSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::SXFunction::SXFunction(const std::pair< SXDict, std::vector< std::string > > &arg, const std::pair< SXDict, std::vector< std::string > > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::SXFunction::SXFunction(const std::pair< SXDict, std::vector< std::string > > &arg, const std::vector< SX > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::SXFunction::SXFunction(const std::vector< SX > &arg, const std::pair< SXDict, std::vector< std::string > > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::SXFunction::SXFunction(const std::vector< SX > &arg, const std::vector< SX > &res) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::SdpSolver::SdpSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::SdqpSolver::SdqpSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::Simulator::Simulator(const Integrator &integrator, const Function &output_fcn, const Matrix< double > &grid) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::Simulator::Simulator(const Integrator &integrator, const Matrix< double > &grid) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::SocpSolver::SocpSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception casadi::StabilizedQpSolver::StabilizedQpSolver(const std::string &solver, const std::map< std::string, Sparsity > &st) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}