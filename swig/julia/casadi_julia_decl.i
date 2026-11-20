// casadi runs GC-safe: it may wait on threads that call back into Julia
%exception {
  try {
    SwigJlGcSafe gc_safe;
    $action
  } catch(const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

// Julia docstrings are Markdown: prototypes go in code fences
%feature("customdoc:main", "$brief\n\n```\n$overview```\n$main");
%feature("customdoc:group", "\n---\n\n```\n$group```\n$main\n");

// Dicts, nested vectors and GenericType return as the Julia value
%define %julia_value_return(xType...)
%typemap(out, ctype="jl_value_t*", jltype="Any", noblock=1, fragment="casadi_all") xType %{ if(!($result = swig::from< xType >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output."); %}
%typemap(out, ctype="jl_value_t*", jltype="Any", noblock=1, fragment="casadi_all") const xType& %{ if(!($result = swig::from< xType >(*$1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output."); %}
%typemap(jlout) xType, const xType& "_check($call)"
%enddef
%julia_value_return(casadi::GenericType)
%julia_value_return(std::vector<casadi::GenericType>)
%julia_value_return(std::map<std::string, casadi::Matrix<double> >)
%julia_value_return(std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%julia_value_return(std::map<std::string, casadi::MX>)
%julia_value_return(std::map<std::string, casadi::Sparsity>)
%julia_value_return(std::vector<std::vector<casadi::Matrix<double> > >)
%julia_value_return(std::vector<std::vector<casadi::Matrix<casadi::SXElem> > >)
%julia_value_return(std::vector<std::vector<casadi::MX> >)
%julia_value_return(std::vector<std::vector<casadi::Sparsity> >)
%julia_value_return(std::vector<std::vector<casadi_int> >)
%julia_value_return(std::vector<std::vector<double> >)
%julia_value_return(std::map<std::string, casadi::GenericType>)
%typemap(jlout) casadi::Dict, const casadi::Dict& "_check($call)::AbstractDict"
// std::size_t (hashes) like size_t
%typemap(out, ctype="long long", jltype="Clonglong") std::size_t %{ $result = (long long)$1; %}
%typemap(jlout) std::size_t "Int(_check($call))"
// Matrix<casadi_int> as a Julia Matrix{Int}
%typemap(out, ctype="void *", jltype="Ptr{Cvoid}") casadi::Matrix<casadi_int> %{ $result = (void*) new casadi::Matrix<double>($1); %}
%typemap(jlout) casadi::Matrix<casadi_int> "Base.round.(Int, Base.Matrix{Float64}(DM(_check($call))))"
// an option value is any Julia value GenericType converts from
%typemap(jlparam) const casadi::GenericType& option_value "Any"
// director returns convert as Base.convert does (true for 1, 2 for 2.0), not as the strict input typemaps
%typemap(directorout) bool, casadi_int, double %{ $result = Swig::swig_director_as< $1_ltype >($input); %}
// Julia parameter unions for the matrix types, nested like Python's precedence
%typemap(jlparam) casadi::Matrix<double>, const casadi::Matrix<double>& "Union{DM,Real,AbstractArray{<:Real},Vector{Any}}"
%typemap(jlparam) casadi::Matrix<casadi::SXElem>, const casadi::Matrix<casadi::SXElem>& "Union{SX,DM,Real,AbstractArray{<:Real},Vector{Any}}"
%typemap(jlparam) casadi::MX, const casadi::MX& "Union{MX,DM,Real,AbstractArray{<:Real},Vector{Any}}"
// a Callback is a Function, as in C++
%typemap(jlparam) casadi::Function, const casadi::Function& "Union{Function,CallbackFunction}"
// Containers of matrices: the common element types; other values (Vector{Any}) reach the runtime dispatcher
%define %julia_partial_containers(T, V, D)
%typemap(jlparam, partial="1") std::vector< T >, const std::vector< T >&
  "Union{" V "Vector{DM},Vector{Float64},Vector{Int64}}"
%typemap(jlparam, partial="1") std::map<std::string, T >, const std::map<std::string, T >&
  "Union{" D "Dict{String,DM},Dict{String,Float64},Dict{String,Int64}}"
%enddef
%julia_partial_containers(casadi::Matrix<double>, "", "")
%julia_partial_containers(casadi::Matrix<casadi::SXElem>, "Vector{SX},", "Dict{String,SX},")
%julia_partial_containers(casadi::MX, "Vector{MX},", "Dict{String,MX},")
%typemap(jlparam, partial="1") std::vector<casadi::Sparsity>, const std::vector<casadi::Sparsity>& "Vector{Sparsity}"

// Julia concatenates with [a b] / [a; b]; builders stay as _horzcat/_vertcat
%rename(_vertcat) casadi_vertcat;
%rename(_horzcat) casadi_horzcat;
// Python's names; Julia's own spelling is & | ! (and .& .| .!)
%rename(logic_and) casadi_and;
%rename(logic_or) casadi_or;
%rename(logic_not) casadi_not;
%rename(logic_all) casadi_all;
%rename(logic_any) casadi_any;

// inv is the matrix inverse in Julia; elementwise inverse is 1 ./ x
%ignore casadi::MX::inv() const;

// users subtype the abstract Callback, as Python subclasses Callback; the proxy is CallbackFunction
%rename(CallbackFunction) casadi::Callback;
%feature("julia:director", "Callback") casadi::Callback;
// a Callback cannot be copied (as in C++); CallbackFunction(obj) takes a Julia Callback
%ignore casadi::Callback::Callback(const Callback&);
// eval_buffer: arg and res reach Julia as vectors of Vector{Float64} views (nothing for a null)
%{
namespace casadi {
  jl_value_t* julia_buffers(double** p, const std::vector<casadi_int>& sz) {
    jl_value_t* vt = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
    jl_array_t* a = jl_alloc_array_1d(jl_apply_array_type((jl_value_t*)jl_any_type, 1), sz.size());
    JL_GC_PUSH2(&vt, &a);
    for (size_t i = 0; i < sz.size(); ++i) {
      jl_value_t* e = p[i] ? (jl_value_t*)jl_ptr_to_array_1d(vt, p[i], sz[i], 0) : jl_nothing;
      jl_array_ptr_set(a, i, e);
    }
    JL_GC_POP();
    return (jl_value_t*)a;
  }
} // namespace casadi
%}
%typemap(directorin, noblock=1) (const double** arg, const std::vector<casadi_int>& sizes_arg) {
  $input = casadi::julia_buffers(const_cast<double**>($1), $2);
}
%typemap(directorin, noblock=1) (double** res, const std::vector<casadi_int>& sizes_res) {
  $input = casadi::julia_buffers($1, $2);
}

// variable/parameter/subject_to record the call site; the Julia layer wraps these
%rename(_variable) casadi::Opti::variable;
%rename(_parameter) casadi::Opti::parameter;
%rename(_subject_to) casadi::Opti::subject_to;
// the Opti proxy keeps the OptiCallback it registers, as Python's self._cb
%feature("julia:keep") casadi::Opti::callback_class;
// copy(::Opti) is an independent problem; the copy constructor shares it
%feature("julia:copy", "copy") casadi::Opti;

// internal: opti.callback takes any Julia callable
%rename(_OptiCallback) casadi::OptiCallback;
%rename(_callback_class) casadi::Opti::callback_class;
%rename(_call) casadi::OptiCallback::call;

// keep the casadi_ prefix: solver(...) is the setter
%rename(casadi_solver) casadi::OptiAdvanced::casadi_solver;
make_property(casadi::OptiAdvanced, casadi_solver, Function);

%typemap(out) casadi::native_DM;
// Native Julia value: Float64, Vector or Matrix
%typemap(jlout) casadi::native_DM "(function (d::DM) n = nonzeros(densify(d)); r, c = size1(d), size2(d); r * c == 1 ? n[1] : (r == 1 || c == 1) ? n : Base.reshape(n, r, c) end)(DM(_check($call)))"
