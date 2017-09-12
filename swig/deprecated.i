%exception  casadi::Function::fullJacobian() const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::spCanEvaluate(bool fwd) const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Sparsity::print_compact(std::ostream &stream=casadi::userOut()) const  {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}