%exception  casadi::GenericMatrix< MX  >::size() const {
 CATCH_OR_NOT(DEPRECATED_MSG("")) 
}
%exception  casadi::GenericMatrix< MatType >::size() const  {
 CATCH_OR_NOT(DEPRECATED_MSG("")) 
}
%exception  casadi::GenericMatrix< Matrix< DataType >  >::size() const {
 CATCH_OR_NOT(DEPRECATED_MSG("")) 
}