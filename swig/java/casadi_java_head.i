// Data structure in the target language holding data
#define GUESTOBJECT _jobject
#define CASADI_NOT_NONE(p) p

%begin %{
#define SWIG_JAVA_GUEST_EXCEPTION "org/casadi/CasadiException"
%}
%include <javaguest.swg>
%include <javacleaner.swg>
%include <enums.swg>
// Java has single inheritance: the empty mixin bases become interfaces, so each
// class extends its first real base (DM extends GenDM, MX extends GenMX)
%include <swiginterface.i>
%interface_impl(casadi::PrintableCommon);
%interface_impl(casadi::GenericExpressionCommon);
%interface_impl(casadi::GenericMatrixCommon);
%interface_impl(casadi::SparsityInterfaceCommon);
%interface_impl(casadi::MatrixCommon);
// toString() is str(), as Python's __str__ (an interface cannot override Object methods)
%typemap(javacode) casadi::Matrix<double>, casadi::Matrix<casadi::SXElem>, casadi::Sparsity,
  casadi::Function, casadi::Slice, casadi::DaeBuilder, casadi::Importer, casadi::Linsol,
  casadi::NlpBuilder, casadi::Opti, casadi::OptiSol, casadi::XmlFile, casadi::Resource, casadi::GraphBuilder,
  casadi::ModelicaParser %{
  @Override
  public String toString() {
    return str();
  }
%}
// f(a, b, ...) and f(map) as Python's __call__; Kotlin spells it f(a, b)
%extend casadi::Function {
  %proxycode %{
  /** One output is returned as is, several as a List; numbers, numeric arrays and Lists of numbers become DM. */
  public Object invoke(Object... args) {
    java.util.List<GenericMatrix> in = new java.util.ArrayList<>();
    for (Object a : args) in.add(swigMatrix(a));
    java.util.List<GenericMatrix> r = call(in);
    return r.size() == 1 ? r.get(0) : r;
  }

  /** Call with inputs by name; returns the outputs by name. */
  public java.util.Map<String, GenericMatrix> invoke(java.util.Map<String, ?> args) {
    java.util.Map<String, GenericMatrix> in = new java.util.LinkedHashMap<>();
    for (java.util.Map.Entry<String, ?> e : args.entrySet()) in.put(e.getKey(), swigMatrix(e.getValue()));
    return call(in);
  }

  private static GenericMatrix swigMatrix(Object x) {
    if (x instanceof GenericMatrix) return (GenericMatrix) x;
    if (x instanceof Number) return new DM(((Number) x).doubleValue());
    if (x instanceof double[]) return new DM((double[]) x);
    if (x instanceof double[][]) return new DM((double[][]) x);
    if (x instanceof long[]) return new DM(java.util.Arrays.stream((long[]) x).asDoubleStream().toArray());
    if (x instanceof int[]) return new DM(java.util.Arrays.stream((int[]) x).asDoubleStream().toArray());
    if (x instanceof java.util.List && ((java.util.List<?>) x).stream().allMatch(e -> e instanceof Number))
      return new DM(((java.util.List<?>) x).stream().mapToDouble(e -> ((Number) e).doubleValue()).toArray());
    throw new IllegalArgumentException("Cannot convert " + (x == null ? "null" : x.getClass().getName()) + " to a matrix");
  }
  %}
}
%rename(hash) __hash__;
%rename(nonzero) __nonzero__;
// DM truth value as Python's DM.__bool__ (an empty scalar is false)
%ignore casadi::Matrix<double>::__nonzero__;
// C++ streams and raw pointers have no Java form
%ignore casadi::StreamStateGuard;
%ignore casadi::normalized_setup;
%ignore casadi::normalized_in;
%ignore casadi::conic_debug(const Function&, std::ostream&);
%ignore casadi::Options::best_matches;
%ignore casadi::to_slice2;
%ignore casadi::hash_sparsity(casadi_int, casadi_int, const casadi_int*, const casadi_int*);
%ignore deserialize(std::istream&);
%ignore deserialize(DeserializingStream&);
%ignore serialize(SerializingStream&) const;
%ignore casadi::SerializerBase::pack(const std::vector<Linsol>&);
%ignore casadi::DeserializerBase::blind_unpack_linsol_vector;
%ignore casadi::DeserializerBase::unpack_linsol_vector;
%ignore test_cast;
%ignore casadi::MX::simplify_ref_count(std::vector<MX>&, std::vector<MX>&, const Dict&);
%ignore casadi::MX::simplify_const_folding(std::vector<MX>&, std::vector<MX>&, const Dict&);
%ignore casadi::MX::simplify_combine_terms(std::vector<MX>&, std::vector<MX>&, const Dict&);
// MX extends GenMX, so SharedObject's members are spelled out
%warnfilter(813) casadi::MX;
%extend casadi::MX {
  std::string toString() const { return $self->get_str(); }
  std::string class_name() const { return $self->class_name(); }
  bool is_null() const { return $self->is_null(); }
  casadi_int hash() const { return $self->__hash__(); }
}
%extend casadi::GenericType { std::string toString() const { return $self->get_str(); } }
