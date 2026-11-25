// Kotlin operators (plus, minus, times, div, rem) and named relatives, through
// the module functions, as Python's operators go through _casadi.<op>
%define %java_binop(Self, Op, F)
  public Self Op(Self y) { return ca.F(this, y); }
  public Self Op(double y) { return ca.F(this, y); }
%enddef
%define %java_binop_dm(Op, F)
  %java_binop(DM, Op, F)
  public SX Op(SX y) { return ca.F(this, y); }
  public MX Op(MX y) { return ca.F(this, y); }
%enddef
%define %java_binop_sym(Self, Op, F)
  %java_binop(Self, Op, F)
  public Self Op(DM y) { return ca.F(this, y); }
%enddef
%define %java_binops(M, Self)
  M(Self, plus, plus) M(Self, minus, minus) M(Self, times, times) M(Self, div, rdivide) M(Self, rdivide, rdivide)
  M(Self, rem, mod) M(Self, mtimes, mtimes) M(Self, power, power)
  M(Self, lt, lt) M(Self, le, le) M(Self, gt, gt) M(Self, ge, ge) M(Self, eq, eq) M(Self, ne, ne)
%enddef
%define %java_binop_dm2(Self, Op, F) %java_binop_dm(Op, F) %enddef
// Math as methods, for Matrix<M>
%define %java_math(Self)
  public Self dot(Self y) { return ca.dot(this, y); }
  public Self fmin(Self y) { return ca.fmin(this, y); }
  public Self fmax(Self y) { return ca.fmax(this, y); }
  public Self atan2(Self x) { return ca.atan2(this, x); }
  public Self sin() { return ca.sin(this); }
  public Self cos() { return ca.cos(this); }
  public Self tan() { return ca.tan(this); }
  public Self asin() { return ca.asin(this); }
  public Self acos() { return ca.acos(this); }
  public Self atan() { return ca.atan(this); }
  public Self sinh() { return ca.sinh(this); }
  public Self cosh() { return ca.cosh(this); }
  public Self tanh() { return ca.tanh(this); }
  public Self exp() { return ca.exp(this); }
  public Self log() { return ca.log(this); }
  public Self sqrt() { return ca.sqrt(this); }
  public Self abs() { return ca.abs(this); }
  public Self floor() { return ca.floor(this); }
  public Self ceil() { return ca.ceil(this); }
  public Self sign() { return ca.sign(this); }
  public Self sum1() { return ca.sum1(this); }
  public Self sum2() { return ca.sum2(this); }
  public Self sumsqr() { return ca.sumsqr(this); }
  public Self norm_1() { return ca.norm_1(this); }
  public Self norm_2() { return ca.norm_2(this); }
  public Self norm_inf() { return ca.norm_inf(this); }
  public Self mmin() { return ca.mmin(this); }
  public Self mmax() { return ca.mmax(this); }
%enddef
// Operation codes as the Operation enum
%define %java_ops(Self)
  public boolean is_op(Operation op) { return is_op(op.swigValue()); }
  public static Self binary(Operation op, Self x, Self y) { return binary(op.swigValue(), x, y); }
  public static Self unary(Operation op, Self x) { return unary(op.swigValue(), x); }
%enddef
%extend casadi::Matrix<double> {
  %proxycode {
  %java_binops(%java_binop_dm2, DM)
  %java_ops(DM)
  %java_math(DM)
  }
}
%extend casadi::Matrix<casadi::SXElem> {
  %proxycode {
  %java_binops(%java_binop_sym, SX)
  %java_ops(SX)
  %java_math(SX)
  }
}
%extend casadi::MX {
  %proxycode {
  %java_binops(%java_binop_sym, MX)
  %java_ops(MX)
  %java_math(MX)
  }
}

%extend casadi::DeserializerBase {
  %proxycode %{
  /** Unpacks the next item, whatever its type */
  public Object unpack() {
    String name = SerializerBase.type_to_string(pop_type());
    try {
      return getClass().getMethod("blind_unpack_" + name).invoke(this);
    } catch (java.lang.reflect.InvocationTargetException e) {
      if (e.getCause() instanceof RuntimeException) throw (RuntimeException) e.getCause();
      throw new RuntimeException(e.getCause());
    } catch (ReflectiveOperationException e) {
      throw new IllegalStateException("cannot unpack " + name, e);
    }
  }
  %}
}

// A DM: a Java value has one static type
%typemap(jni) casadi::native_DM "jobject"
%typemap(jtype) casadi::native_DM "Object"
%typemap(jstype) casadi::native_DM "DM"
%typemap(javaout) casadi::native_DM { return $imclassname.swigCast($jnicall); }

%clear casadi::MX variable, casadi::MX parameter;
%clear const casadi::MX& g, const std::vector<casadi::MX>& g;

%extend casadi::Opti {
  %proxycode %{
  private OptiCallback swigCallback;

  /** Calls fh with the iteration number after each solver iteration; null removes it */
  public void callback(java.util.function.LongConsumer fh) {
    if (fh == null) {
      callback_class();
      swigCallback = null;
      return;
    }
    swigCallback = new OptiCallback() {
      @Override
      public void call(long i) {
        fh.accept(i);
      }
    };
    callback_class(swigCallback);
  }
  %}
}

namespace casadi{
  %extend Matrix<double> {
    // double[] as a column, double[][] one array per row
    Matrix(const std::vector<double>& x) { return new casadi::Matrix<double>(x); }
    Matrix(const std::vector< std::vector<double> >& x) { return new casadi::Matrix<double>(x); }
    Matrix(const casadi::Sparsity& sp, const std::vector<double>& nz) { return new casadi::Matrix<double>(sp, nz); }
    bool nonzero() const {
      casadi_assert($self->numel()==1, "Only a scalar can be cast to a float");
      return $self->nnz()!=0 && $self->scalar()!=0;
    }
  }
  %extend WeakRef {
    WeakRef(const casadi::MX& x) { return new casadi::WeakRef(x); }
  }
  %extend MX {
    MX(const std::vector<double>& x) { return new casadi::MX(casadi::DM(x)); }
    MX(const std::vector< std::vector<double> >& x) { return new casadi::MX(casadi::DM(x)); }
    // Parametric indexing (0-based)
    MX get(const casadi::MX& rr) const { casadi::MX m; $self->get(m, false, rr); return m; }
    MX get(const casadi::MX& rr, const casadi::Slice& cc) const { casadi::MX m; $self->get(m, false, rr, cc); return m; }
    MX get(const casadi::Slice& rr, const casadi::MX& cc) const { casadi::MX m; $self->get(m, false, rr, cc); return m; }
    MX get(const casadi::MX& rr, const casadi::MX& cc) const { casadi::MX m; $self->get(m, false, rr, cc); return m; }
  }
} // namespace casadi
