// Rust types accepted where casadi converts (see the traits)
%typemap(rsaccept) casadi::Matrix<double> "DM | Sparsity | f64 | i64 | bool | Vec<f64> | Vec<i64> | Vec<Vec<f64>>"
%typemap(rsaccept) casadi::Matrix<casadi::SXElem> "SX | DM | Sparsity | f64 | i64 | bool | Vec<f64> | Vec<i64> | Vec<Vec<f64>>"
%typemap(rsaccept) casadi::MX "MX | DM | Sparsity | f64 | i64 | bool | Vec<f64> | Vec<i64> | Vec<Vec<f64>>"
%typemap(rstype) casadi::Matrix<casadi_int> "DM"
%typemap(rsaccept) casadi::Matrix<casadi_int> "i64 | Vec<i64> | Vec<Vec<i64>> | Sparsity"
%typemap(rstype) casadi::SXElem "SX"
%typemap(rsaccept) casadi::SXElem "SX | f64 | i64"
%typemap(rsaccept) casadi::Slice "Slice | i64"
%typemap(rstype) casadi::GenericType "rt::Value"
%typemap(rsaccept) casadi::GenericType "Any"
// Numeric evaluation returns a Result, symbolic evaluation panics
%feature("rust:result") casadi::Function::call(const std::vector<DM>&, std::vector<DM>&, bool, bool) const;
%feature("rust:result") casadi::Function::call(const DMDict&, DMDict&, bool, bool) const;
%rename(call_symbolic) casadi::Function::call(const std::vector<SX>&, std::vector<SX>&, bool, bool) const;
%rename(call_symbolic) casadi::Function::call(const std::vector<MX>&, std::vector<MX>&, bool, bool) const;
%rename(call_symbolic) casadi::Function::call(const SXDict&, SXDict&, bool, bool) const;
%rename(call_symbolic) casadi::Function::call(const MXDict&, MXDict&, bool, bool) const;
// Overloads the plain/_with rule cannot name; (nrow, ncol) pairs duplicate two-integer overloads
%typemap(rsaccept) std::pair<casadi_int, casadi_int> ""
%ignore casadi::Slice::Slice(casadi_int, bool);
%ignore casadi::Function::Function(const std::string&);
%rename(named) casadi::Function::Function(const std::string&, const std::vector<SX>&, const std::vector<SX>&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&);
%rename(named) casadi::Function::Function(const std::string&, const std::vector<MX>&, const std::vector<MX>&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&);
%rename(from_dict) casadi::Function::Function(const std::string&, const std::map<std::string, SX>&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&);
%rename(from_dict) casadi::Function::Function(const std::string&, const std::map<std::string, MX>&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&);

%rename(from_sparsity) casadi::MX::MX(const Sparsity&);
%rename(from_ccs) casadi::Sparsity::Sparsity(casadi_int, casadi_int, const std::vector<casadi_int>&, const std::vector<casadi_int>&, bool);
%rename(sym_sparsity) sym(const std::string&, const Sparsity&);
%rename(sym_vec) sym(const std::string&, const Sparsity&, casadi_int);
%rename(sym_vec) sym(const std::string&, casadi_int, casadi_int, casadi_int);
%rename(sym_vecvec) sym(const std::string&, const Sparsity&, casadi_int, casadi_int);
%rename(sym_vecvec) sym(const std::string&, casadi_int, casadi_int, casadi_int, casadi_int);
%rename(map_threads) casadi::Function::map(casadi_int, const std::string&, casadi_int) const;
%rename(map_reduce) casadi::Function::map(casadi_int, const std::vector<bool>&, const std::vector<bool>&, const Dict&) const;
%rename(mapaccum_n) casadi::Function::mapaccum(const std::string&, casadi_int, casadi_int, const Dict&) const;
%rename(variable_sparsity) casadi::Opti::variable(const Sparsity&, const std::string&);
%rename(parameter_sparsity) casadi::Opti::parameter(const Sparsity&, const std::string&);
%rename(to_function_named) casadi::Opti::to_function(const std::string&, const std::vector<MX>&, const std::vector<MX>&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&);
%rename(to_function_dict) casadi::Opti::to_function(const std::string&, const std::map<std::string, MX>&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&);
%rename(add_var) casadi::DaeBuilder::add(const std::string&, const std::string&, const std::string&, const Dict&);
%rename(symvar_of) casadi::OptiAdvanced::symvar(const MX&) const;
%rename(variable_symbol) casadi::Opti::variable(const MX&, const std::string&);
%rename(parameter_symbol) casadi::Opti::parameter(const MX&, const std::string&);
%rename(clear_constraints) casadi::Opti::subject_to();
%rename(zeros_sparsity) zeros(const Sparsity&);
%rename(ones_sparsity) ones(const Sparsity&);
%rename(inf_sparsity) inf(const Sparsity&);
%rename(nan_sparsity) nan(const Sparsity&);
%rename(rand_sparsity) rand(const Sparsity&);
%rename(reshape_sparsity) casadi::casadi_reshape(const Sparsity&, const Sparsity&);
%rename(reshape_sparsity) casadi::casadi_reshape(const MX&, const Sparsity&);
%rename(reshape_sparsity) casadi::casadi_reshape(const Matrix<double>&, const Sparsity&);
%rename(reshape_sparsity) casadi::casadi_reshape(const Matrix<SXElem>&, const Sparsity&);
%rename(logic_and) casadi_and;
%rename(logic_or) casadi_or;
%rename(logic_not) casadi_not;
%rename(logic_all) casadi_all;
%rename(logic_any) casadi_any;
// Overloads that are not a truncation of the _with form
%rename(map_named) casadi::Function::map(const std::string&, const std::string&, casadi_int, const std::vector<casadi_int>&, const std::vector<casadi_int>&, const Dict&) const;
%rename(map_named) casadi::Function::map(const std::string&, const std::string&, casadi_int, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&) const;
%rename(mapaccum_named) casadi::Function::mapaccum(const std::string&, casadi_int, const Dict&) const;
%rename(export_code_file) casadi::Function::export_code(const std::string&, const std::string&, const Dict&) const;
%rename(integrator_opts) casadi::integrator(const std::string&, const std::string&, const SXDict&, const Dict&);
%rename(integrator_opts) casadi::integrator(const std::string&, const std::string&, const MXDict&, const Dict&);
%rename(integrator_opts) casadi::integrator(const std::string&, const std::string&, const Function&, const Dict&);
%rename(external_opts) casadi::external(const std::string&, const Dict&);
%rename(subject_to_scaled) casadi::Opti::subject_to(const MX&, const DM&, const Dict&);
%rename(subject_to_scaled) casadi::Opti::subject_to(const std::vector<MX>&, const DM&, const Dict&);
%rename(add_causality) casadi::DaeBuilder::add(const std::string&, const std::string&, const Dict&);
%rename(add_expr) casadi::DaeBuilder::add(const std::string&, const std::string&, const std::string&, const MX&, const Dict&);
%rename(create_default) casadi::DaeBuilder::create() const;
%rename(create_default) casadi::GraphBuilder::create() const;
%rename(create_io) casadi::DaeBuilder::create(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&) const;
%rename(create_io) casadi::GraphBuilder::create(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&, const Dict&) const;
%rename(erase_linear) erase(const std::vector<casadi_int>&, bool);
%ignore casadi::DaeBuilder::create(const std::string&, const std::vector<std::string>&, const std::vector<std::string>&, bool, bool) const;
// Overridden from Rust through CallbackImpl::eval; the wrapper calls the C++ implementation
%rename(eval_dm) casadi::Callback::eval;

%rename(uminus) operator-() const;

// Opti values return as DM
%typemap(out, doc="DM", tsstub_out="DM", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = swig::from< casadi::DM >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'DM'.");
}
