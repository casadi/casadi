%exception  casadi::Function::forward(const std::vector< DM > &arg, const std::vector< DM > &res, const std::vector< std::vector< DM > > &fseed, std::vector< std::vector< DM > > &output_fsens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::forward(const std::vector< MX > &arg, const std::vector< MX > &res, const std::vector< std::vector< MX > > &fseed, std::vector< std::vector< MX > > &output_fsens, bool always_inline=false, bool never_inline=false) {
 CATCH_OR_NOT(DEPRECATED_MSG("") $action)
}
%exception  casadi::Function::forward(const std::vector< SX > &arg, const std::vector< SX > &res, const std::vector< std::vector< SX > > &fseed, std::vector< std::vector< SX > > &output_fsens, bool always_inline=false, bool never_inline=false) {
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