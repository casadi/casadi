// Python's names; and/or/not would read as C# keywords in docs
%rename(logic_and) casadi_and;
%rename(logic_or) casadi_or;
%rename(logic_not) casadi_not;
%rename(logic_all) casadi_all;
%rename(logic_any) casadi_any;

// C# operators are defined on the proxies below; the C++ ones stay unwrapped
%rename(uminus) operator-() const;
%rename(uplus) operator+() const;
%rename(plus) casadi::Sparsity::operator+;
%rename(times) casadi::Sparsity::operator*;
%rename(minus) casadi::Slice::operator-;
%rename(times) casadi::Slice::operator*;
%ignore operator==;
%ignore operator!=;
%ignore casadi::DaeBuilder::operator();

// DM from a list of SX constants: rarely used, and it would make ca.DM([1, 2, 3]) ambiguous with DM(double[])
%ignore casadi::Matrix<double>::Matrix(const std::vector<casadi::SXElem>&);

// Java's names for the Python dunders
%rename(hash) __hash__;
%rename(nonzero) __nonzero__;

// In C#, MX extends GenMX, not SharedObject: a WeakRef of an MX upcasts in C++
%extend casadi::WeakRef {
  WeakRef(const casadi::MX& shared) { return new casadi::WeakRef(shared); }
}

// C++ streams and raw pointers have no C# form
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
%ignore test_cast;

// Interface-internal macros
%ignore CASADI_INT_TYPE_STR; %ignore MNAME; %ignore L_INT; %ignore L_BOOL; %ignore LL; %ignore LR;
%ignore L_DICT; %ignore L_DOUBLE; %ignore L_STR; %ignore LABEL; %ignore SWIG_IND1;

// casadi printout goes through Console.Out and Console.Error, so .NET redirections see it
%{
namespace casadi {
  typedef void (*dotnet_write_t)(const char* s, int n, int error);
  typedef void (*dotnet_flush_t)(int error);
  static dotnet_write_t dotnet_write = 0;
  static dotnet_flush_t dotnet_flush = 0;
  static void dotnetlogger(const char* s, std::streamsize num, bool error) {
    if (!dotnet_write) return Logger::writeDefault(s, num, error);
    dotnet_write(s, static_cast<int>(num), error ? 1 : 0);
  }
  static void dotnetflush(bool error) {
    if (!dotnet_flush) return Logger::flushDefault(error);
    dotnet_flush(error ? 1 : 0);
  }
} // namespace casadi
%}
%wrapper %{
extern "C" SWIGEXPORT void casadi_dotnet_logger(casadi::dotnet_write_t write, casadi::dotnet_flush_t flush) {
  casadi::dotnet_write = write;
  casadi::dotnet_flush = flush;
  casadi::Logger::writeFun = casadi::dotnetlogger;
  casadi::Logger::flush = casadi::dotnetflush;
}
%}
%insert("dotnetmodule") %{
  [System.Runtime.InteropServices.UnmanagedFunctionPointer(System.Runtime.InteropServices.CallingConvention.Cdecl)]
  internal delegate void SwigLogWrite(System.IntPtr s, int n, int error);
  [System.Runtime.InteropServices.UnmanagedFunctionPointer(System.Runtime.InteropServices.CallingConvention.Cdecl)]
  internal delegate void SwigLogFlush(int error);
  [System.Runtime.InteropServices.DllImport(SwigNative.Lib, CallingConvention = System.Runtime.InteropServices.CallingConvention.Cdecl)]
  static extern void casadi_dotnet_logger(SwigLogWrite write, SwigLogFlush flush);
  // Kept alive for the whole process: C++ holds pointers to them
  static readonly SwigLogWrite logWrite = (s, n, error) =>
    (error != 0 ? System.Console.Error : System.Console.Out).Write(System.Runtime.InteropServices.Marshal.PtrToStringUTF8(s, n));
  static readonly SwigLogFlush logFlush = error => (error != 0 ? System.Console.Error : System.Console.Out).Flush();
  internal static void InstallLogger() { casadi_dotnet_logger(logWrite, logFlush); }
%}
%insert("dotnetinit") %{
    ca.InstallLogger();
%}

// Opti values return as DM (implicitly convertible to double, double[], double[,])
%typemap(out, doc="DM", tsstub_out="DM", noblock=1, fragment="casadi_all") casadi::native_DM {
  if(!($result = swig::from< casadi::DM >($1))) SWIG_exception_fail(SWIG_TypeError,"Failed to convert output to type 'DM'.");
}

// The free functions live in the module class (static class ca): using static CasADi.ca;
%insert("dotnetmodule") %{
  public const double inf = double.PositiveInfinity;
  public const double nan = double.NaN;
  public const double pi = System.Math.PI;

  // pow as in System.Math (casadi's name is power)
  public static double pow(double x, double n) => System.Math.Pow(x, n);
  public static DM pow(DM x, DM n) => power(x, n);
  public static SX pow(SX x, SX n) => power(x, n);
  public static MX pow(MX x, MX n) => power(x, n);

  // Factories of the array constructors, as the generated ones of the C++ constructors:
  // DM([1, 2, 3]) is a column, DM([[1, 2], [3, 4]]) a matrix from its rows
  public static DM DM(double[] x) => new DM(x);
  public static DM DM(double[][] rows) => new DM(rows);
  public static DM DM(double[,] x) => new DM(x);
  public static SX SX(double[] x) => new SX(x);
  public static SX SX(double[][] rows) => new SX(rows);
  public static SX SX(double[,] x) => new SX(x);
  public static MX MX(double[] x) => new MX(x);
  public static MX MX(double[][] rows) => new MX(rows);
  public static MX MX(double[,] x) => new MX(x);

  // Python's names
  public static double fabs(double x) => System.Math.Abs(x);
  public static DM fabs(DM x) => abs(x);
  public static SX fabs(SX x) => abs(x);
  public static MX fabs(MX x) => abs(x);
  public static DM vcat(System.Collections.Generic.IEnumerable<DM> x) => vertcat(x);
  public static SX vcat(System.Collections.Generic.IEnumerable<SX> x) => vertcat(x);
  public static MX vcat(System.Collections.Generic.IEnumerable<MX> x) => vertcat(x);
  public static DM vcat(System.Collections.Generic.IEnumerable<double> x) => new DM(System.Linq.Enumerable.ToArray(x));
  public static DM hcat(System.Collections.Generic.IEnumerable<DM> x) => horzcat(x);
  public static SX hcat(System.Collections.Generic.IEnumerable<SX> x) => horzcat(x);
  public static MX hcat(System.Collections.Generic.IEnumerable<MX> x) => horzcat(x);
  public static DM hcat(System.Collections.Generic.IEnumerable<double> x) => new DM(System.Linq.Enumerable.ToArray(x)).T();

  /// <summary>Readable text of any value casadi returns, as Python's str(): dictionaries
  /// (call results, stats(), options), lists, arrays and casadi objects</summary>
  public static string str(object x) {
    switch (x) {
      case null: return "None";
      case string s: return "'" + s + "'";
      case bool b: return b ? "True" : "False";
      case double d: return d.ToString("R", System.Globalization.CultureInfo.InvariantCulture);
      case float f: return ((double)f).ToString("R", System.Globalization.CultureInfo.InvariantCulture);
      case System.IFormattable n when !(x is System.Enum): return n.ToString(null, System.Globalization.CultureInfo.InvariantCulture);
      case SwigObject o: return o.ToString();
      case double[,] m: {
        var rows = new System.Collections.Generic.List<string>();
        for (int i = 0; i < m.GetLength(0); ++i) {
          var r = new System.Collections.Generic.List<string>();
          for (int j = 0; j < m.GetLength(1); ++j) r.Add(str(m[i, j]));
          rows.Add("[" + string.Join(", ", r) + "]");
        }
        return "[" + string.Join(", ", rows) + "]";
      }
      case System.Collections.IDictionary d: {
        var items = new System.Collections.Generic.List<string>();
        foreach (System.Collections.DictionaryEntry e in d) items.Add(str(e.Key) + ": " + str(e.Value));
        return "{" + string.Join(", ", items) + "}";
      }
      case System.Runtime.CompilerServices.ITuple t: {
        var items = new System.Collections.Generic.List<string>();
        for (int i = 0; i < t.Length; ++i) items.Add(str(t[i]));
        return "(" + string.Join(", ", items) + ")";
      }
      case System.Collections.IEnumerable e: {
        var items = new System.Collections.Generic.List<string>();
        foreach (object v in e) items.Add(str(v));
        return "[" + string.Join(", ", items) + "]";
      }
    }
    return x.ToString();
  }

  // Concatenations of nothing are empty DMs, as in Python
  public static DM vertcat() => new DM(0, 1);
  public static DM horzcat() => new DM(1, 0);
  public static DM diagcat() => new DM(0, 0);
  public static DM veccat() => new DM(0, 1);
  // Concatenations with any number of arguments, as Python's vertcat(a, b)
  public static DM vertcat(params DM[] x) => x.Length == 0 ? new DM(0, 1) : vertcat((System.Collections.Generic.IEnumerable<DM>)x);
  public static SX vertcat(params SX[] x) => vertcat((System.Collections.Generic.IEnumerable<SX>)x);
  public static MX vertcat(params MX[] x) => vertcat((System.Collections.Generic.IEnumerable<MX>)x);
  public static DM horzcat(params DM[] x) => x.Length == 0 ? new DM(1, 0) : horzcat((System.Collections.Generic.IEnumerable<DM>)x);
  public static SX horzcat(params SX[] x) => horzcat((System.Collections.Generic.IEnumerable<SX>)x);
  public static MX horzcat(params MX[] x) => horzcat((System.Collections.Generic.IEnumerable<MX>)x);
  public static DM diagcat(params DM[] x) => x.Length == 0 ? new DM(0, 0) : diagcat((System.Collections.Generic.IEnumerable<DM>)x);
  public static SX diagcat(params SX[] x) => diagcat((System.Collections.Generic.IEnumerable<SX>)x);
  public static MX diagcat(params MX[] x) => diagcat((System.Collections.Generic.IEnumerable<MX>)x);
  public static DM veccat(params DM[] x) => x.Length == 0 ? new DM(0, 1) : veccat((System.Collections.Generic.IEnumerable<DM>)x);
  public static SX veccat(params SX[] x) => veccat((System.Collections.Generic.IEnumerable<SX>)x);
  public static MX veccat(params MX[] x) => veccat((System.Collections.Generic.IEnumerable<MX>)x);
%}

// var (a, b) = f.call(x, y): the outputs of a call (any List) unpack into variables, as in Python
%insert("dotnet") %{
public static class SwigListDeconstruct {
  public static void Deconstruct<T>(this System.Collections.Generic.List<T> l, out T a, out T b) {
    Check(l, 2); a = l[0]; b = l[1];
  }
  public static void Deconstruct<T>(this System.Collections.Generic.List<T> l, out T a, out T b, out T c) {
    Check(l, 3); a = l[0]; b = l[1]; c = l[2];
  }
  public static void Deconstruct<T>(this System.Collections.Generic.List<T> l, out T a, out T b, out T c, out T d) {
    Check(l, 4); a = l[0]; b = l[1]; c = l[2]; d = l[3];
  }
  public static void Deconstruct<T>(this System.Collections.Generic.List<T> l, out T a, out T b, out T c, out T d, out T e) {
    Check(l, 5); a = l[0]; b = l[1]; c = l[2]; d = l[3]; e = l[4];
  }
  public static void Deconstruct<T>(this System.Collections.Generic.List<T> l, out T a, out T b, out T c, out T d, out T e, out T f) {
    Check(l, 6); a = l[0]; b = l[1]; c = l[2]; d = l[3]; e = l[4]; f = l[5];
  }
  public static void Deconstruct<T>(this System.Collections.Generic.List<T> l, out T a, out T b, out T c, out T d, out T e, out T f, out T g) {
    Check(l, 7); a = l[0]; b = l[1]; c = l[2]; d = l[3]; e = l[4]; f = l[5]; g = l[6];
  }
  public static void Deconstruct<T>(this System.Collections.Generic.List<T> l, out T a, out T b, out T c, out T d, out T e, out T f, out T g, out T h) {
    Check(l, 8); a = l[0]; b = l[1]; c = l[2]; d = l[3]; e = l[4]; f = l[5]; g = l[6]; h = l[7];
  }
  static void Check<T>(System.Collections.Generic.List<T> l, int n) {
    if (l.Count != n) throw new System.ArgumentException("cannot unpack " + l.Count + " values into " + n + " variables");
  }
}
%}

// opti.callback(f): f(iteration) at every iteration, as Python's opti.callback(fh); f is any C# method
// taking a long (or a lambda). OptiCallback and callback_class are the plumbing behind it, hidden from
// users. As Python's self._cb, the Opti keeps the wrapper of f alive (C++ holds only a pointer).
%feature("dotnet:internal") casadi::OptiCallback;
%feature("dotnet:internal") casadi::Opti::callback_class;
%extend casadi::Opti {
  %proxycode %{
  OptiCallback swigCallback;
  /// <summary>Call f(iteration) at every solver iteration, as Python's opti.callback(fh):
  /// opti.callback(MyMonitor) with void MyMonitor(long iter), or opti.callback(i => ...)</summary>
  public void callback(System.Action<long> f) {
    var cb = new SwigCallback(f);
    callback_class(cb);
    swigCallback = cb;
  }
  /// <summary>Remove the iteration callback</summary>
  public void callback() {
    callback_class();
    swigCallback = null;
  }
  sealed class SwigCallback : OptiCallback {
    readonly System.Action<long> f;
    public SwigCallback(System.Action<long> f) { this.f = f; }
    internal override void call(long i) => f(i);
  }
  %}
}

// opti.bounded(lb, x, ub) on the instance, as in Python (C++ has it static)
%ignore casadi::Opti::bounded;
%rename(bounded) casadi::Opti::swig_bounded;
%extend casadi::Opti {
  casadi::MX swig_bounded(const casadi::MX& lb, const casadi::MX& expr, const casadi::MX& ub) const {
    return casadi::Opti::bounded(lb, expr, ub);
  }
}

// Arithmetic, comparisons and indexing of the matrix types; numbers and arrays convert implicitly
%define %dotnet_matrix_ops(T)
  public static T operator +(T x, T y) => ca.plus(x, y);
  public static T operator -(T x, T y) => ca.minus(x, y);
  /// <summary>Elementwise product, as Python's *; ca.mtimes is the matrix product</summary>
  public static T operator *(T x, T y) => ca.times(x, y);
  public static T operator /(T x, T y) => ca.rdivide(x, y);
  public static T operator -(T x) => ca.minus(new T(0.0), x);
  public static T operator +(T x) => x;
  /// <summary>Elementwise equality as an expression (x == y in constraints), as in Python; not a
  /// reference comparison: test for null with x is null</summary>
  public static T operator ==(T x, T y) => ca.eq(x, y);
  public static T operator !=(T x, T y) => ca.ne(x, y);
  public static T operator <(T x, T y) => ca.lt(x, y);
  public static T operator <=(T x, T y) => ca.le(x, y);
  public static T operator >(T x, T y) => ca.gt(x, y);
  public static T operator >=(T x, T y) => ca.ge(x, y);

  /// <summary>Entry k (0-based, column-major; negative counts from the end)</summary>
  public T this[long k] {
    get => get(new DM((double)k));
    set => set(new DM((double)k), value);
  }
  // An int literal matches this exactly (it would be ambiguous between long and Index)
  public T this[int k] {
    get => this[(long)k];
    set => this[(long)k] = value;
  }
  /// <summary>Entry k counted from the end with ^: x[^1] is the last entry</summary>
  public T this[System.Index k] {
    get => this[k.IsFromEnd ? -(long)k.Value : k.Value];
    set => this[k.IsFromEnd ? -(long)k.Value : k.Value] = value;
  }
  /// <summary>Entries in a range of linear indices, e.g. x[1..3], x[^2..]</summary>
  public T this[Slice k] {
    get => get(k);
    set => set(k, value);
  }
  /// <summary>Entries at the given linear indices</summary>
  public T this[long[] k] {
    get => get((DM)k);
    set => set((DM)k, value);
  }
  /// <summary>Entry (i, j), 0-based</summary>
  public T this[long i, long j] {
    get => get(new DM((double)i), new DM((double)j));
    set => set(new DM((double)i), new DM((double)j), value);
  }
  public T this[Slice i, Slice j] {
    get => get(i, j);
    set => set(i, j, value);
  }
  public T this[long i, Slice j] {
    get => get(new Slice(i), j);
    set => set(new Slice(i), j, value);
  }
  public T this[Slice i, long j] {
    get => get(i, new Slice(j));
    set => set(i, new Slice(j), value);
  }
  /// <summary>Rows/columns at index lists, as Python's x[1, [0, 10, 20]]</summary>
  public T this[long i, long[] j] {
    get => get(ca.DM((double)i), (DM)j);
    set => set(ca.DM((double)i), (DM)j, value);
  }
  public T this[long[] i, long j] {
    get => get((DM)i, ca.DM((double)j));
    set => set((DM)i, ca.DM((double)j), value);
  }
  public T this[long[] i, long[] j] {
    get => get((DM)i, (DM)j);
    set => set((DM)i, (DM)j, value);
  }
  public T this[Slice i, long[] j] {
    get => get(i, (DM)j);
    set => set(i, (DM)j, value);
  }
  public T this[long[] i, Slice j] {
    get => get((DM)i, j);
    set => set((DM)i, j, value);
  }

  public override string ToString() => str();
%enddef

// Host values convert to DM through DM's own conversion only: the constructor dispatcher would also
// try std::vector<std::vector<double>> and turn a 1xn array into an nx1 one
%extend casadi::Matrix<double> {
  static casadi::Matrix<double> swig_from_value(const casadi::Matrix<double>& x) { return x; }
}
%extend casadi::Matrix<double> {
  %proxycode %{
  // Arrays convert in C++, as the copy constructor does for any numeric value (shape kept)
  internal static DM FromValue(object x) => SwigRT.Call1<DM>(SwigNative.swigcs_DM_swig_from_value, 1, x);
  /// <summary>A column vector, as Python's DM([1, 2, 3])</summary>
  public DM(double[] x) : base(SwigRT.Owned.Pending) {
    // null reaches C++, which rejects it (SwigException), as for any other null argument
    if (x != null && x.Length == 0) SwigRT.Construct(this, SwigNative.swigcs_new_DM, 0L, 1L);
    else SwigRT.Construct(this, SwigNative.swigcs_DM_swig_from_value, x);
  }
  /// <summary>A dense matrix, row by row, as Python's DM([[1, 2], [3, 4]])</summary>
  public DM(double[,] x) : base(SwigRT.Owned.Pending) {
    SwigRT.Construct(this, SwigNative.swigcs_DM_swig_from_value, x);
  }
  /// <summary>A dense matrix from its rows, as Python's DM([[1, 2], [3, 4]]): DM([[1, 2], [3, 4]])</summary>
  public DM(double[][] rows) : this(Rows(rows)) {}
  internal static double[,] Rows(double[][] rows) {
    if (rows == null) return null;
    for (int i = 0; i < rows.Length; ++i)
      if (rows[i] == null) throw new System.ArgumentException("row " + i + " is null", nameof(rows));
    int nr = rows.Length, nc = nr == 0 ? 0 : rows[0].Length;
    var m = new double[nr, nc];
    for (int i = 0; i < nr; ++i) {
      if (rows[i].Length != nc) throw new System.ArgumentException("rows of different lengths: " + nc + " and " + rows[i].Length);
      for (int j = 0; j < nc; ++j) m[i, j] = rows[i][j];
    }
    return m;
  }
  public static implicit operator DM(double x) => new DM(x);
  /// <summary>A column vector</summary>
  public static implicit operator DM(double[] x) => x is null ? null : new DM(x);
  public static implicit operator DM(long[] x) => x is null ? null : x.Length == 0 ? new DM(0, 1) : FromValue(x);
  public static implicit operator DM(double[,] x) => x is null ? null : new DM(x);
  public static implicit operator DM(System.Collections.Generic.List<double> x) => x is null ? null : (DM)x.ToArray();
  public static explicit operator double(DM x) => x.to_double();
  /// <summary>Dense values, row by row (as Java's full())</summary>
  public double[,] full() => ToArray();
  /// <summary>Dense values, row by row</summary>
  public double[,] ToArray() {
    long nr = size1(), nc = size2();
    double[] d = ca.densify(this).nonzeros().ToArray();
    double[,] r = new double[nr, nc];
    for (long j = 0; j < nc; ++j) for (long i = 0; i < nr; ++i) r[i, j] = d[i + j * nr];
    return r;
  }
  %}
  %proxycode {
  %dotnet_matrix_ops(DM)
  }
}
%extend casadi::Matrix<casadi::SXElem> {
  %proxycode %{
  public static implicit operator SX(double x) => new SX(x);
  public static implicit operator SX(DM x) => new SX(x);
  public static implicit operator SX(double[] x) => new SX((DM)x);
  public static implicit operator SX(double[,] x) => new SX((DM)x);
  /// <summary>A constant column vector / dense matrix (row by row)</summary>
  public SX(double[] x) : base(SwigRT.Owned.Pending) { SwigRT.Construct(this, SwigNative.swigcs_new_SX, new DM(x)); }
  public SX(double[,] x) : base(SwigRT.Owned.Pending) { SwigRT.Construct(this, SwigNative.swigcs_new_SX, new DM(x)); }
  public SX(double[][] rows) : base(SwigRT.Owned.Pending) { SwigRT.Construct(this, SwigNative.swigcs_new_SX, new DM(rows)); }
  %}
  %proxycode {
  %dotnet_matrix_ops(SX)
  }
}
%extend casadi::MX {
  %proxycode %{
  public static implicit operator MX(double x) => new MX(x);
  public static implicit operator MX(DM x) => new MX(x);
  public static implicit operator MX(double[] x) => new MX((DM)x);
  public static implicit operator MX(double[,] x) => new MX((DM)x);
  /// <summary>A constant column vector / dense matrix (row by row)</summary>
  public MX(double[] x) : base(SwigRT.Owned.Pending) { SwigRT.Construct(this, SwigNative.swigcs_new_MX, new DM(x)); }
  public MX(double[,] x) : base(SwigRT.Owned.Pending) { SwigRT.Construct(this, SwigNative.swigcs_new_MX, new DM(x)); }
  public MX(double[][] rows) : base(SwigRT.Owned.Pending) { SwigRT.Construct(this, SwigNative.swigcs_new_MX, new DM(rows)); }
  %}
  %proxycode {
  %dotnet_matrix_ops(MX)
  }
}

// Ranges index like Python slices: x[1..3], x[^1], x[..]
%extend casadi::Slice {
  %proxycode %{
  // No implicit long -> Slice: an integer argument would be ambiguous between the Slice and DM index overloads
  /// <summary>One entry counted with an Index, e.g. x[.., ^1]</summary>
  public static implicit operator Slice(System.Index i) => new Slice(i.IsFromEnd ? -(long)i.Value : i.Value);
  public static implicit operator Slice(System.Range r) {
    long start = r.Start.IsFromEnd ? -r.Start.Value : r.Start.Value;
    long stop = r.End.IsFromEnd ? (r.End.Value == 0 ? long.MaxValue : -r.End.Value) : r.End.Value;
    return new Slice(start, stop);
  }
  public override string ToString() => str();
  %}
}

%extend casadi::Sparsity {
  %proxycode %{
  public override string ToString() => str();
  %}
}

// Printing shows casadi's description, as Python's str() and Java's toString()
// (SharedObject covers Function, Opti, DaeBuilder, Importer, Linsol, XmlFile, Resource, ...)
%extend casadi::SharedObject {
  %proxycode %{
  public override string ToString() => str();
  %}
}
// Classes with their own str() (C# str() is not virtual, SharedObject's would print nothing)
%extend casadi::Opti {
  %proxycode %{
  public override string ToString() => str();
  %}
}
%extend casadi::OptiSol {
  %proxycode %{
  public override string ToString() => str();
  %}
}
%extend casadi::NlpBuilder {
  %proxycode %{
  public override string ToString() => str();
  %}
}

// Opti records the C# caller of variable/parameter/subject_to in the user dict ("stacktrace"),
// as Python, MATLAB and Java do; error messages then point at the user's code
%rename(swig_variable) casadi::Opti::variable;
%rename(swig_parameter) casadi::Opti::parameter;
%rename(swig_subject_to) casadi::Opti::subject_to;
%extend casadi::Opti {
  %proxycode %{
  static System.Collections.Generic.Dictionary<string, object> SwigStacktrace(bool all) {
    var st = new System.Collections.Generic.List<object>();
    foreach (var fr in new System.Diagnostics.StackTrace(1, true).GetFrames()) {
      var m = fr.GetMethod();
      if (m == null || m.DeclaringType == null || typeof(Opti).IsAssignableFrom(m.DeclaringType)) continue;
      st.Add(new System.Collections.Generic.Dictionary<string, object> {
        ["file"] = System.IO.Path.GetFileName(fr.GetFileName() ?? ""),
        ["line"] = (long)fr.GetFileLineNumber(),
        ["name"] = m.Name });
      if (!all) break;
    }
    return new System.Collections.Generic.Dictionary<string, object> { ["stacktrace"] = st };
  }
  static System.Collections.Generic.IDictionary<string, object> SwigOpts(System.Collections.Generic.IDictionary<string, object> o) =>
    o ?? new System.Collections.Generic.Dictionary<string, object>();

  /// <summary>Create a decision variable (symbol): n x m, with an attribute ("full", "symmetric")</summary>
  public MX variable(long n = 1, long m = 1, string attribute = "full") {
    MX r = swig_variable(n, m, attribute);
    update_user_dict(r, SwigStacktrace(false));
    return r;
  }
  public MX variable(Sparsity sp, string attribute = "full") {
    MX r = swig_variable(sp, attribute);
    update_user_dict(r, SwigStacktrace(false));
    return r;
  }
  public MX variable(MX symbol, string attribute = "full") {
    MX r = swig_variable(symbol, attribute);
    update_user_dict(r, SwigStacktrace(false));
    return r;
  }
  /// <summary>Create a parameter (symbol): fixed during optimization, set with set_value</summary>
  public MX parameter(long n = 1, long m = 1, string attribute = "full") {
    MX r = swig_parameter(n, m, attribute);
    update_user_dict(r, SwigStacktrace(false));
    return r;
  }
  public MX parameter(Sparsity sp, string attribute = "full") {
    MX r = swig_parameter(sp, attribute);
    update_user_dict(r, SwigStacktrace(false));
    return r;
  }
  public MX parameter(MX symbol, string attribute = "full") {
    MX r = swig_parameter(symbol, attribute);
    update_user_dict(r, SwigStacktrace(false));
    return r;
  }
  /// <summary>Add constraints: g is a (chained) comparison such as x &gt;= 0 or 0 &lt;= x &lt;= 1</summary>
  public void subject_to(MX g, System.Collections.Generic.IDictionary<string, object> options = null) {
    try { swig_subject_to(g, SwigOpts(options)); } finally { update_user_dict(g, SwigStacktrace(true)); }
  }
  public void subject_to(System.Collections.Generic.IEnumerable<MX> g, System.Collections.Generic.IDictionary<string, object> options = null) {
    var l = new System.Collections.Generic.List<MX>(g);
    try { swig_subject_to(l, SwigOpts(options)); } finally { update_user_dict(l, SwigStacktrace(true)); }
  }
  public void subject_to(MX g, DM linear_scale, System.Collections.Generic.IDictionary<string, object> options = null) {
    try { swig_subject_to(g, linear_scale, SwigOpts(options)); } finally { update_user_dict(g, SwigStacktrace(true)); }
  }
  public void subject_to(System.Collections.Generic.IEnumerable<MX> g, DM linear_scale, System.Collections.Generic.IDictionary<string, object> options = null) {
    var l = new System.Collections.Generic.List<MX>(g);
    try { swig_subject_to(l, linear_scale, SwigOpts(options)); } finally { update_user_dict(l, SwigStacktrace(true)); }
  }
  /// <summary>Clear all constraints (the objective, variables, parameters and values are kept)</summary>
  public void subject_to() => swig_subject_to();
  %}
}

// Generic unpack, as Python's and Java's: the next item, whatever its type
%extend casadi::DeserializerBase {
  %proxycode %{
  /// <summary>Unpacks the next item, whatever its type</summary>
  public object unpack() {
    string name = SerializerBase.type_to_string(pop_type());
    var m = GetType().GetMethod("blind_unpack_" + name, System.Type.EmptyTypes);
    if (m == null) throw new System.InvalidOperationException("cannot unpack " + name);
    try {
      return m.Invoke(this, null);
    } catch (System.Reflection.TargetInvocationException e) when (e.InnerException != null) {
      System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(e.InnerException).Throw();
      throw;
    }
  }
  %}
}

// Plugins and headers next to the casadi_dotnet library (as Java's NativeLoader directory)
%{
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <sys/stat.h>
namespace casadi {
  static std::string dotnet_library_directory() {
    std::string path;
#ifdef _WIN32
    HMODULE h = 0;
    char buf[MAX_PATH];
    if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                           reinterpret_cast<LPCSTR>(&dotnet_library_directory), &h)
        && GetModuleFileNameA(h, buf, MAX_PATH)) path = buf;
#else
    Dl_info info;
    if (dladdr(reinterpret_cast<void*>(&dotnet_library_directory), &info) && info.dli_fname) path = info.dli_fname;
#endif
    size_t k = path.find_last_of("/\\");
    return k == std::string::npos ? std::string() : path.substr(0, k);
  }
  static bool dotnet_is_dir(const std::string& p) {
    struct stat s;
    return stat(p.c_str(), &s) == 0 && (s.st_mode & S_IFDIR);
  }
} // namespace casadi
%}
%init %{
  {
    std::string dir = casadi::dotnet_library_directory();
    if (!dir.empty()) {
      casadi::GlobalOptions::setCasadiPath(dir);
      for (const char* rel : {"/include", "/../include"}) {
        if (casadi::dotnet_is_dir(dir + rel)) {
          casadi::GlobalOptions::setCasadiIncludePath(dir + rel);
          break;
        }
      }
    }
  }
%}

// f.Call(x, y): the first output, as Python's f(x, y) for a single-output Function
%extend casadi::Function {
  %proxycode %{
  // f.call(a, b): all outputs (a collection expression f.call([a, b]) of DMs would be ambiguous with the SX/MX overloads)
  public System.Collections.Generic.List<DM> call(params DM[] args) => call((System.Collections.Generic.IEnumerable<DM>)args);
  public System.Collections.Generic.List<SX> call(params SX[] args) => call((System.Collections.Generic.IEnumerable<SX>)args);
  public System.Collections.Generic.List<MX> call(params MX[] args) => call((System.Collections.Generic.IEnumerable<MX>)args);
  public DM Call(params DM[] args) => call((System.Collections.Generic.IEnumerable<DM>)args)[0];
  public SX Call(params SX[] args) => call((System.Collections.Generic.IEnumerable<SX>)args)[0];
  public MX Call(params MX[] args) => call((System.Collections.Generic.IEnumerable<MX>)args)[0];
  public override string ToString() => str();
  %}
}
