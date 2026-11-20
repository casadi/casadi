// Every wrapper turns C++ exceptions into CasadiException (javaguest.swg)
%exception {
  $action
}

%fragment(SWIG_From_frag(bool));
%fragment(SWIG_From_frag(long long));
%fragment(SWIG_From_frag(double));
%{
  namespace casadi {
    inline jclass java_class(const char* name) {
      JNIEnv* e = SWIG_JavaGuestEnvChecked();
      jclass c = e->FindClass(name);
      SWIG_JavaGuestCheck(e);
      jclass g = static_cast<jclass>(e->NewGlobalRef(c));
      e->DeleteLocalRef(c);
      return g;
    }
    // casadi printout through System.out and System.err
    static void javalogger(const char* s, std::streamsize num, bool error) {
      JNIEnv* e = SWIG_JavaGuestEnv();
      if (!e || e->ExceptionCheck()) return Logger::writeDefault(s, num, error);
      static jclass sys = java_class("java/lang/System");
      static jclass ps = java_class("java/io/PrintStream");
      static jfieldID out = e->GetStaticFieldID(sys, "out", "Ljava/io/PrintStream;");
      static jfieldID err = e->GetStaticFieldID(sys, "err", "Ljava/io/PrintStream;");
      static jmethodID print = e->GetMethodID(ps, "print", "(Ljava/lang/String;)V");
      jobject stream = e->GetStaticObjectField(sys, error ? err : out);
      jobject str = SWIG_JavaGuestNewString(std::string(s, num));
      e->CallVoidMethod(stream, print, str);
      e->DeleteLocalRef(str);
      e->DeleteLocalRef(stream);
      if (e->ExceptionCheck()) e->ExceptionClear();
    }
    static void javaflush(bool error) {
      JNIEnv* e = SWIG_JavaGuestEnv();
      if (!e || e->ExceptionCheck()) return Logger::flushDefault(error);
      static jclass sys = java_class("java/lang/System");
      static jclass ps = java_class("java/io/PrintStream");
      static jfieldID out = e->GetStaticFieldID(sys, "out", "Ljava/io/PrintStream;");
      static jfieldID err = e->GetStaticFieldID(sys, "err", "Ljava/io/PrintStream;");
      static jmethodID flush = e->GetMethodID(ps, "flush", "()V");
      jobject stream = e->GetStaticObjectField(sys, error ? err : out);
      e->CallVoidMethod(stream, flush);
      e->DeleteLocalRef(stream);
      if (e->ExceptionCheck()) e->ExceptionClear();
    }
    // Thread.interrupt() of the calling thread interrupts casadi
    static bool javacheckinterrupted() {
      JNIEnv* e = SWIG_JavaGuestEnv();
      if (!e || e->ExceptionCheck()) return false;
      static jclass th = java_class("java/lang/Thread");
      static jmethodID cur = e->GetStaticMethodID(th, "currentThread", "()Ljava/lang/Thread;");
      static jmethodID is = e->GetMethodID(th, "isInterrupted", "()Z");
      jobject t = e->CallStaticObjectMethod(th, cur);
      if (e->ExceptionCheck()) return false;
      bool r = e->CallBooleanMethod(t, is);
      e->DeleteLocalRef(t);
      if (e->ExceptionCheck()) e->ExceptionClear();
      return r;
    }
    static void javaclearinterrupted() {
      JNIEnv* e = SWIG_JavaGuestEnv();
      if (!e || e->ExceptionCheck()) return;
      static jclass th = java_class("java/lang/Thread");
      static jmethodID clear = e->GetStaticMethodID(th, "interrupted", "()Z");
      e->CallStaticBooleanMethod(th, clear);
      if (e->ExceptionCheck()) e->ExceptionClear();
    }
  } // namespace casadi
%}
%init %{
  casadi::Logger::writeFun = casadi::javalogger;
  casadi::Logger::flush = casadi::javaflush;
  casadi::InterruptHandler::checkInterrupted = casadi::javacheckinterrupted;
  casadi::InterruptHandler::clearInterrupted = casadi::javaclearinterrupted;
  // Plugins and headers next to the casadi_java library that NativeLoader loaded
  {
    JNIEnv *env = SWIG_JavaGuestEnvChecked();
    jclass c = env->FindClass("org/casadi/NativeLoader");
    jfieldID f = c ? env->GetStaticFieldID(c, "directory", "Ljava/lang/String;") : 0;
    jobject d = f ? env->GetStaticObjectField(c, f) : 0;
    jmethodID m = c ? env->GetStaticMethodID(c, "includeDirectory", "()Ljava/lang/String;") : 0;
    jobject inc = m ? env->CallStaticObjectMethod(c, m) : 0;
    env->ExceptionClear();
    if (d) casadi::GlobalOptions::setCasadiPath(SWIG_JavaGuestString(d));
    if (inc) casadi::GlobalOptions::setCasadiIncludePath(SWIG_JavaGuestString(inc));
  }
%}
%pragma(java) jniclasscode=%{
  static {
    NativeLoader.load();
  }
%}
// Free functions live in the module class, like java.lang.Math
%pragma(java) moduleclassmodifiers="public final class"
%pragma(java) modulecode=%{
  private ca() {}

  public static final double inf = Double.POSITIVE_INFINITY;
  public static final double nan = Double.NaN;
  public static final double pi = Math.PI;

  public static <T> T vertcat(java.util.List<T> x) { return x.isEmpty() ? caJNI.swigCast(new DM(0, 1)) : swigVertcat(x); }
  public static <T> T horzcat(java.util.List<T> x) { return x.isEmpty() ? caJNI.swigCast(new DM(1, 0)) : swigHorzcat(x); }
  public static <T> T diagcat(java.util.List<T> x) { return x.isEmpty() ? caJNI.swigCast(new DM(0, 0)) : swigDiagcat(x); }
  public static <T> T veccat(java.util.List<T> x) { return x.isEmpty() ? caJNI.swigCast(new DM(0, 1)) : swigVeccat(x); }

  public static DM vertcat() { return new DM(0, 1); }
  public static DM horzcat() { return new DM(1, 0); }
  public static DM diagcat() { return new DM(0, 0); }
  public static DM veccat() { return new DM(0, 1); }

  // Varargs forms of the List concatenations, as Python's vertcat(a, b)
  public static DM vertcat(DM... x) { return vertcat(java.util.Arrays.asList(x)); }
  public static SX vertcat(SX... x) { return vertcat(java.util.Arrays.asList(x)); }
  public static MX vertcat(MX... x) { return vertcat(java.util.Arrays.asList(x)); }
  public static Sparsity vertcat(Sparsity... x) { return vertcat(java.util.Arrays.asList(x)); }
  public static DM horzcat(DM... x) { return horzcat(java.util.Arrays.asList(x)); }
  public static SX horzcat(SX... x) { return horzcat(java.util.Arrays.asList(x)); }
  public static MX horzcat(MX... x) { return horzcat(java.util.Arrays.asList(x)); }
  public static Sparsity horzcat(Sparsity... x) { return horzcat(java.util.Arrays.asList(x)); }
  public static DM diagcat(DM... x) { return diagcat(java.util.Arrays.asList(x)); }
  public static SX diagcat(SX... x) { return diagcat(java.util.Arrays.asList(x)); }
  public static MX diagcat(MX... x) { return diagcat(java.util.Arrays.asList(x)); }
  public static Sparsity diagcat(Sparsity... x) { return diagcat(java.util.Arrays.asList(x)); }
  public static DM veccat(DM... x) { return veccat(java.util.Arrays.asList(x)); }
  public static SX veccat(SX... x) { return veccat(java.util.Arrays.asList(x)); }
  public static MX veccat(MX... x) { return veccat(java.util.Arrays.asList(x)); }
%}
// Merged overloads over the matrix types take <T extends GenericMatrix>
%feature("java:typebound", "GenericMatrix");
%typemap(javainterfaces) casadi::Matrix<double>, casadi::Matrix<casadi::SXElem>, casadi::MX "AutoCloseable, GenericMatrix"

%java_type("String", "Ljava/lang/String;", std::string)
%java_type("java.util.List<String>", "Ljava/util/List;", std::vector<std::string>)
%java_type("java.util.List<java.util.List<String>>", "Ljava/util/List;", std::vector<std::vector<std::string> >)
%java_type("boolean", "Z", bool)
%java_type("boolean[]", "[Z", std::vector<bool>)
%java_type("boolean[][]", "[[Z", std::vector<std::vector<bool> >)
%java_type("long", "J", casadi_int)
%java_type("long[]", "[J", std::pair<casadi_int,casadi_int>)
%java_type("long[]", "[J", std::vector<casadi_int>)
%java_type("long[][]", "[[J", std::vector<std::vector<casadi_int> >)
%java_type("double", "D", double)
%java_type("double[]", "[D", std::vector<double>)
%java_type("double[][]", "[[D", std::vector<std::vector<double> >)
%java_type("Sparsity", "L$packagepath/Sparsity;", casadi::Sparsity)
%java_type("java.util.List<Sparsity>", "Ljava/util/List;", std::vector<casadi::Sparsity>)
%java_type("java.util.List<java.util.List<Sparsity>>", "Ljava/util/List;", std::vector<std::vector<casadi::Sparsity> >)
%java_type_map("Sparsity", std::map<std::string, casadi::Sparsity>)
%java_type_map("java.util.List<Sparsity>", std::map<std::string, std::vector<casadi::Sparsity> >)
%java_type_entry("java.util.Map<String,Sparsity>", "java.util.List<String>", std::pair<std::map<std::string, casadi::Sparsity>, std::vector<std::string> >)
%java_type_widen("DM", "L$packagepath/DM;", "double", casadi::Matrix<double>)
%java_type("java.util.List<DM>", "Ljava/util/List;", std::vector<casadi::Matrix<double> >)
%java_type("java.util.List<java.util.List<DM>>", "Ljava/util/List;", std::vector<std::vector<casadi::Matrix<double> > >)
%java_type_map("DM", std::map<std::string, casadi::Matrix<double> >)
%java_type("SX", "L$packagepath/SX;", casadi::SXElem)
%java_type("java.util.List<SX>", "Ljava/util/List;", std::vector<casadi::SXElem>)
%java_type_widen("SX", "L$packagepath/SX;", "DM,double", casadi::Matrix<casadi::SXElem>)
%java_type("java.util.List<SX>", "Ljava/util/List;", std::vector<casadi::Matrix<casadi::SXElem> >)
%java_type("java.util.List<java.util.List<SX>>", "Ljava/util/List;", std::vector<std::vector<casadi::Matrix<casadi::SXElem> > >)
%java_type_map("SX", std::map<std::string, casadi::Matrix<casadi::SXElem> >)
%java_type_widen("MX", "L$packagepath/MX;", "DM,double", casadi::MX)
%java_type("java.util.List<MX>", "Ljava/util/List;", std::vector<casadi::MX>)
%java_type("java.util.List<java.util.List<MX>>", "Ljava/util/List;", std::vector<std::vector<casadi::MX> >)
%java_type_map("MX", std::map<std::string, casadi::MX>)
%java_type_entry("MX", "MX", std::pair<casadi::MX, casadi::MX>)
%java_type("DM", "L$packagepath/DM;", casadi::Matrix<casadi_int>)
%java_type("Object", "Ljava/lang/Object;", casadi::GenericType)
%java_type("java.util.List<Object>", "Ljava/util/List;", std::vector<casadi::GenericType>)
%java_type("java.util.List<java.util.List<Object>>", "Ljava/util/List;", std::vector<std::vector<casadi::GenericType> >)
%java_type_widen("Slice", "L$packagepath/Slice;", "long", casadi::Slice)
%java_type("Function", "L$packagepath/Function;", casadi::Function)
%java_type("java.util.List<Function>", "Ljava/util/List;", std::vector<casadi::Function>)
%java_type_entry("Function", "Function", std::pair<casadi::Function, casadi::Function>)
%java_type_map("Object", std::map<std::string, casadi::GenericType>)
// Options are read from any Map with String keys
%typemap(jstype) const std::map<std::string, casadi::GenericType>& "java.util.Map<String,?>"
%java_type_map("java.util.List<String>", std::map<std::string, std::vector<std::string> >)

// Concatenations of an empty List give a DM, as in Python (see the module code)
%rename(swigVertcat) casadi_vertcat;
%rename(swigHorzcat) casadi_horzcat;
%rename(swigDiagcat) casadi_diagcat;
%rename(swigVeccat) casadi_veccat;
%javamethodmodifiers casadi_vertcat "private";
%javamethodmodifiers casadi_horzcat "private";
%javamethodmodifiers casadi_diagcat "private";
%javamethodmodifiers casadi_veccat "private";
%rename(logic_and) casadi_and;
%rename(logic_or) casadi_or;
%rename(logic_not) casadi_not;
%rename(logic_all) casadi_all;
%rename(logic_any) casadi_any;

// Kotlin operator names; Sparsity and Slice keep is_equal / equality members unwrapped
%rename(unaryMinus) operator-() const;
%rename(unaryPlus) operator+() const;
%rename(plus) casadi::Sparsity::operator+;
%rename(times) casadi::Sparsity::operator*;
%rename(minus) casadi::Slice::operator-;
%rename(times) casadi::Slice::operator*;
%ignore operator==;
%ignore operator!=;
%ignore casadi::DaeBuilder::operator();
// Object.finalize() is called by the garbage collector
%rename(finalize_) casadi::Callback::finalize;
// Slice(long, boolean) would clash with the proxy constructor (pointer, ownership)
%ignore casadi::Slice::Slice(casadi_int i, bool ind1);
%extend casadi::Slice {
  Slice(casadi_int i) { return new casadi::Slice(i); }
}
// Interface-internal macros
%ignore CASADI_INT_TYPE_STR; %ignore MNAME; %ignore L_INT; %ignore L_BOOL; %ignore LL; %ignore LR;
%ignore L_DICT; %ignore L_DOUBLE; %ignore L_STR; %ignore LABEL; %ignore SWIG_IND1;

// Indexing (0-based) as Kotlin's get/set operators: x[i], x[i, j], x[i] = v
%define %java_index_pair(Type, R, C)
    Type get(const R& rr, const C& cc) const {
      Type m;
      $self->get(m, false, rr, cc);
      return m;
    }
    void set(const R& rr, const C& cc, const Type& m) { $self->set(m, false, rr, cc); }
%enddef

// eval_buffer: arg and res reach Java as Lists of direct DoubleBuffers over casadi's memory (null for a null)
%{
namespace casadi {
  jobject java_buffers(double** p, const std::vector<casadi_int>& sz) {
    JNIEnv* e = SWIG_JavaGuestEnvChecked();
    static jclass bb = java_class("java/nio/ByteBuffer"), bo = java_class("java/nio/ByteOrder");
    static jmethodID order = e->GetMethodID(bb, "order", "(Ljava/nio/ByteOrder;)Ljava/nio/ByteBuffer;");
    static jmethodID as_double = e->GetMethodID(bb, "asDoubleBuffer", "()Ljava/nio/DoubleBuffer;");
    static jmethodID native_order = e->GetStaticMethodID(bo, "nativeOrder", "()Ljava/nio/ByteOrder;");
    jobject no = e->CallStaticObjectMethod(bo, native_order);
    SWIG_JavaGuestCheck(e);
    jobject r = SWIG_JavaGuestListNew(static_cast<jint>(sz.size()));
    for (size_t i = 0; i < sz.size(); ++i) {
      jobject d = 0;
      if (p[i]) {
        jobject b = e->NewDirectByteBuffer(p[i], static_cast<jlong>(sz[i] * sizeof(double)));
        SWIG_JavaGuestCheck(e);
        jobject o = e->CallObjectMethod(b, order, no);
        SWIG_JavaGuestCheck(e);
        d = e->CallObjectMethod(o, as_double);
        SWIG_JavaGuestCheck(e);
        e->DeleteLocalRef(o);
        e->DeleteLocalRef(b);
      }
      SWIG_JavaGuestListAppend(r, d);
    }
    e->DeleteLocalRef(no);
    return r;
  }
  // Pointers and lengths of a List of direct buffers
  bool java_from_buffers(jobject list, std::vector<double*>& p, std::vector<casadi_int>& sz) {
    JNIEnv* e = SWIG_JavaGuestEnvChecked();
    if (!list || !SWIG_JavaGuestIsList(list)) return false;
    jint n = SWIG_JavaGuestListSize(list);
    p.resize(n);
    sz.resize(n);
    for (jint i = 0; i < n; ++i) {
      jobject b = SWIG_JavaGuestListGet(list, i);
      p[i] = b ? static_cast<double*>(e->GetDirectBufferAddress(b)) : 0;
      sz[i] = b ? static_cast<casadi_int>(e->GetDirectBufferCapacity(b)) : 0;
      if (b && !p[i]) {
        SWIG_JavaGuestRelease(b);
        return false;
      }
      SWIG_JavaGuestRelease(b);
    }
    return true;
  }
} // namespace casadi
%}
%define %java_buffers(PTR, SIZES)
%typemap(jni) (PTR, SIZES) "jobject"
%typemap(jtype) (PTR, SIZES) "Object"
%typemap(jstype) (PTR, SIZES) "java.util.List<java.nio.DoubleBuffer>"
%typemap(javain) (PTR, SIZES) "$javainput"
%typemap(javadirectorin) (PTR, SIZES) "$imclassname.swigCast($jniinput)"
%typemap(directorin, descriptor="Ljava/util/List;", noblock=1) (PTR, SIZES) {
  $input = casadi::java_buffers(const_cast<double**>($1), $2);
}
%typemap(in, noblock=1) (PTR, SIZES) (std::vector<double*> p, std::vector<casadi_int> sz) {
  if (!casadi::java_from_buffers($input, p, sz)) SWIG_exception_fail(SWIG_TypeError, "Expected a List of direct DoubleBuffers");
  $1 = casadi::get_ptr(p);
  $2 = &sz;
}
%enddef
%java_buffers(const double** arg, const std::vector<casadi_int>& sizes_arg)
%java_buffers(double** res, const std::vector<casadi_int>& sizes_res)

// FunctionBuffer reads and writes direct DoubleBuffers in place; the proxy keeps them alive
%ignore _function_buffer_eval;
%ignore casadi::FunctionBuffer::_self;
%nocopyctor casadi::FunctionBuffer;
%rename(eval) casadi::FunctionBuffer::_eval;
%rename(swigSetArg) casadi::FunctionBuffer::set_arg;
%rename(swigSetRes) casadi::FunctionBuffer::set_res;
%javamethodmodifiers casadi::FunctionBuffer::set_arg "private";
%javamethodmodifiers casadi::FunctionBuffer::set_res "private";
%define %java_direct_buffer(PTR)
%typemap(jni) (PTR, casadi_int size) "jobject"
%typemap(jtype) (PTR, casadi_int size) "Object"
%typemap(jstype) (PTR, casadi_int size) "java.nio.DoubleBuffer"
%typemap(javain) (PTR, casadi_int size) "$javainput"
%typemap(in, noblock=1) (PTR, casadi_int size) {
  $1 = $input ? static_cast<double*>(jenv->GetDirectBufferAddress($input)) : 0;
  if (!$1) SWIG_exception_fail(SWIG_TypeError, "Expected a direct DoubleBuffer");
  $2 = static_cast<casadi_int>(jenv->GetDirectBufferCapacity($input) * sizeof(double));
}
%enddef
%java_direct_buffer(const double* a)
%java_direct_buffer(double* a)
%typemap(javacode) casadi::FunctionBuffer %{
  private final java.util.Map<Long, java.nio.DoubleBuffer> swigBuffers = new java.util.HashMap<>();

  /** Binds a direct DoubleBuffer holding the nonzeros of input i */
  public void set_arg(long i, java.nio.DoubleBuffer a) {
    swigSetArg(i, a);
    swigBuffers.put(i, a);
  }

  /** Binds a direct DoubleBuffer that receives the nonzeros of output i */
  public void set_res(long i, java.nio.DoubleBuffer a) {
    swigSetRes(i, a);
    swigBuffers.put(-1 - i, a);
  }

  @Override
  public String toString() {
    return "FunctionBuffer";
  }
%}
%extend casadi::Function {
  %proxycode %{
  /** A FunctionBuffer for evaluating with minimal overhead */
  public FunctionBuffer buffer() {
    return new FunctionBuffer(this);
  }
  %}
}

// Opti values return as DM
%typemap(out, doc="DM", tsstub_out="DM", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = swig::from< casadi::DM >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'DM'.");
}

// Source location of the Java caller in the user dict of new variables, parameters and constraints
%typemap(javaout) casadi::MX casadi::Opti::variable, casadi::MX casadi::Opti::parameter {
  MX r = $imclassname.swigCast($jnicall);
  update_user_dict(r, swigStacktrace(false));
  return r;
}
%typemap(javaout) void casadi::Opti::subject_to {
    try {
      $jnicall;
    } finally {
      update_user_dict(g, swigStacktrace(true));
    }
  }
// subject_to() has no constraint to record
%ignore casadi::Opti::subject_to();
%javamethodmodifiers casadi::Opti::swig_subject_to "private";
%extend casadi::Opti {
  void swig_subject_to() { $self->subject_to(); }
  %proxycode %{
  /**
   * Clear all constraints.
   *
   * Removes every constraint added so far, including their dual variables; the objective,
   * variables, parameters, values and initial guesses are kept.
   */
  public void subject_to() {
    swig_subject_to();
  }
  %}
}
%extend casadi::Opti {
  %proxycode %{
private static java.util.Map<String, Object> swigStacktrace(boolean all) {
  java.util.List<java.util.Map<String, Object>> st = new java.util.ArrayList<>();
  for (StackTraceElement e : new Throwable().getStackTrace()) {
    if (e.getClassName().equals(Opti.class.getName())) continue;
    java.util.Map<String, Object> f = new java.util.LinkedHashMap<>();
    f.put("file", String.valueOf(e.getFileName()));
    f.put("line", (long) e.getLineNumber());
    f.put("name", e.getMethodName());
    st.add(f);
    if (!all) break;
  }
  java.util.Map<String, Object> meta = new java.util.LinkedHashMap<>();
  meta.put("stacktrace", st);
  return meta;
}
  %}
}
