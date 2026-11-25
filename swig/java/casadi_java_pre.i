// Every default argument as an overload
%feature("compactdefaultargs","0");

// Docstrings keep the \lang groups of the target language
%feature("docstring:lang", "java");

 // Director arguments with the JNI descriptor of their Java type
%define %java_type_director(xDesc, xType...)
%typemap(directorin, descriptor=xDesc, noblock=1, fragment="casadi_all") xType, const xType, const xType& { $input = swig::from< xType >($1); }
%enddef
 // Java type of a casadi-typemapped C++ type (parameters, returns, outputs, directors)
%define %java_type(xJava, xDesc, xType...)
%typemap(jstype) xType, const xType, const xType&, xType &OUTPUT, xType &INOUT xJava
// OUTPUT1..6 copied OUTPUT before Java types were known
%typemap(jstype) xType &OUTPUT1, xType &OUTPUT2, xType &OUTPUT3, xType &OUTPUT4, xType &OUTPUT5, xType &OUTPUT6 xJava
%java_type_director(xDesc, xType)
%enddef
 // ... whose parameters also take the Java types in xWiden
%define %java_type_widen(xJava, xDesc, xWiden, xType...)
%java_type(xJava, xDesc, xType)
%typemap(jstype, widen=xWiden) const xType&, xType xJava
%enddef
 // java.util.Map<String, xValue> and java.util.Map.Entry<xKey, xValue>
%define %java_type_map(xValue, xType...)
%java_type("java.util.Map<String," xValue ">", "Ljava/util/Map;", xType)
%enddef
%define %java_type_entry(xKey, xValue, xType...)
%java_type("java.util.Map.Entry<" xKey "," xValue ">", "Ljava/util/Map$Entry;", xType)
%enddef

%fragment("casadi_containers", "header", fragment="casadi_aux", fragment="StdVectorTraits",
          fragment="casadi_map_traits", fragment="StdPairTraits", fragment=SWIG_Traits_frag(bool),
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(double),
          fragment=SWIG_Traits_frag(std::string)) {
  namespace swig {
    // (casadi_int, casadi_int) is a long[2], or a Map.Entry
    template<> struct traits_asptr<std::pair<casadi_int, casadi_int> > {
      static int asptr(jobject p, std::pair<casadi_int, casadi_int>** m) {
        if (!p || !SWIG_JavaGuestIsLongs(p)) return traits_asptr_stdpair<casadi_int, casadi_int>::asptr(p, m);
        std::vector<long long> v = SWIG_JavaGuestLongs(p);
        if (v.size()!=2) return SWIG_ERROR;
        if (m) *m = new std::pair<casadi_int, casadi_int>(v[0], v[1]);
        return m ? SWIG_NEWOBJ : SWIG_OK;
      }
    };
    template<> struct traits_from<std::pair<casadi_int, casadi_int> > {
      static jobject from(const std::pair<casadi_int, casadi_int>& v) {
        std::vector<long long> r(2);
        r[0] = v.first;
        r[1] = v.second;
        return SWIG_JavaGuestNewLongs(r);
      }
    };
  } // namespace swig
}

// Host conversions: is_null, full/sparse, casadi_host_conv overloads
%fragment("casadi_host", "header", fragment="casadi_decl") {
  namespace casadi {
    bool is_null(GUESTOBJECT *p) {
      return !p;
    }

    GUESTOBJECT* full(const DM& m, bool simplify) {
      return 0;
    }

    // Convert to a sparse matrix
    GUESTOBJECT* sparse(const DM& m) {
      return 0;
    }
  } // namespace casadi

  namespace swig {
    inline int casadi_host_conv(GUESTOBJECT* p, casadi::Slice* m) {
      long long i;
      if (!SWIG_IsOK(swig::asval(p, m ? &i : 0))) return false;
      if (m) {
        m->start = i;
        m->stop = m->start+1;
        if (m->stop==0) m->stop = std::numeric_limits<casadi_int>::max();
      }
      return true;
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::DM* m) {
      // double[] or List of numbers: column; double[][] or nested: one entry per row
      {
        std::vector<double> t;
        if (swig::assign(p, m ? &t : 0)) {
          if (m) *m = t.empty() ? casadi::Matrix<double>(0, 0) : casadi::Matrix<double>(t);
          return true;
        }
      }
      {
        std::vector< std::vector<double> > t;
        if (swig::assign(p, m ? &t : 0)) {
          if (m) *m = casadi::Matrix<double>(t);
          return true;
        }
      }
      return false;
    }

    // An empty Java array of objects still has an element type
    inline bool java_empty_generic(jobject p, casadi::GenericType* m) {
      using casadi::GenericType;
      using casadi::java_class;
      JNIEnv* e = SWIG_JavaGuestEnvChecked();
      static jclass s = java_class("[Ljava/lang/String;"), f = java_class("[Lorg/casadi/Function;"),
        dd = java_class("[[D"), ll = java_class("[[J"), ss = java_class("[[Ljava/lang/String;");
      GenericType g;
      if (e->IsInstanceOf(p, s)) g = std::vector<std::string>();
      else if (e->IsInstanceOf(p, f)) g = std::vector<casadi::Function>();
      else if (e->IsInstanceOf(p, dd)) g = std::vector<std::vector<double> >();
      else if (e->IsInstanceOf(p, ll)) g = std::vector<std::vector<casadi_int> >();
      else if (e->IsInstanceOf(p, ss)) g = std::vector<std::vector<std::string> >();
      else return false;
      if (e->GetArrayLength(static_cast<jarray>(p)) != 0) return false;
      if (m) *m = g;
      return true;
    }
    inline int casadi_host_conv(GUESTOBJECT* p, casadi::GenericType* m) {
      return java_empty_generic(p, m) ? 1 : -1;
    }
  } // namespace swig
}

// Target typemaps of each casadi-typemapped type
%define CASADI_TARGET_INPUT_TYPEMAPS(xName, xStubIn, xTsStub, xType...)
 // Values cross as Object and convert in C++; %java_type sets the Java type
%typemap(jni) xType, const xType, const xType&, xType &INOUT "jobject"
%typemap(jtype) xType, const xType, const xType&, xType &INOUT "Object"
%typemap(jstype) xType, const xType, const xType&, xType &OUTPUT, xType &INOUT "Object"
%typemap(javain) xType, const xType&, xType &INOUT "$javainput"
%typemap(javaout) xType, const xType, const xType& { return $imclassname.swigCast($jnicall); }
%typemap(javadirectorin) xType, const xType, const xType& "$imclassname.swigCast($jniinput)"
%typemap(javadirectorout) xType, const xType "$javacall"
%enddef

%define CASADI_TARGET_OUTPUT_TYPEMAPS(xName, xStubIn, xStubOut, xTsStub, xType...)
 // null is a Java value
%typemap(out, noblock=1, fragment="casadi_all") xType, const xType { $result = swig::from< xType >($1); }
%typemap(out, noblock=1, fragment="casadi_all") const xType& { $result = swig::from< xType >(*$1); }
%java_type_director("Ljava/lang/Object;", xType)
%enddef

// Indexing helpers of the matrix classes
%define %matrix_helpers(Type, xIndex)
    Type get(const casadi::Slice& rr) const {
      Type m;
      $self->get(m, false, rr);
      return m;
    }
    Type get(const casadi::Matrix<casadi_int>& rr) const {
      Type m;
      $self->get(m, false, rr);
      return m;
    }
    Type get(const casadi::Sparsity& sp) const {
      Type m;
      $self->get(m, false, sp);
      return m;
    }
    void set(const casadi::Slice& rr, const Type& m) { $self->set(m, false, rr); }
    void set(const casadi::Matrix<casadi_int>& rr, const Type& m) { $self->set(m, false, rr); }
    void set(const casadi::Sparsity& sp, const Type& m) { $self->set(m, false, sp); }
    %java_index_pair(Type, casadi::Slice, casadi::Slice)
    %java_index_pair(Type, casadi::Slice, casadi::Matrix<casadi_int>)
    %java_index_pair(Type, casadi::Matrix<casadi_int>, casadi::Slice)
    %java_index_pair(Type, casadi::Matrix<casadi_int>, casadi::Matrix<casadi_int>)
    // An integer index is an integer, not a Slice, as in Python
    Type get(casadi_int rr) const {
      Type m;
      $self->get(m, false, casadi::Matrix<casadi_int>(rr));
      return m;
    }
    // long[] indices, as Python's x[[0, 5]]
    Type get(const std::vector<casadi_int>& rr) const {
      Type m;
      $self->get(m, false, casadi::Matrix<casadi_int>(rr));
      return m;
    }
    void set(const std::vector<casadi_int>& rr, const Type& m) { $self->set(m, false, casadi::Matrix<casadi_int>(rr)); }
    Type get(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc) const {
      Type m;
      $self->get(m, false, casadi::Matrix<casadi_int>(rr), casadi::Matrix<casadi_int>(cc));
      return m;
    }
    void set(const std::vector<casadi_int>& rr, const std::vector<casadi_int>& cc, const Type& m) {
      $self->set(m, false, casadi::Matrix<casadi_int>(rr), casadi::Matrix<casadi_int>(cc));
    }
    Type get(const std::vector<casadi_int>& rr, const casadi::Slice& cc) const {
      Type m;
      $self->get(m, false, casadi::Matrix<casadi_int>(rr), cc);
      return m;
    }
    Type get(const casadi::Slice& rr, const std::vector<casadi_int>& cc) const {
      Type m;
      $self->get(m, false, rr, casadi::Matrix<casadi_int>(cc));
      return m;
    }
    Type get_nz(const std::vector<casadi_int>& kk) const {
      Type m;
      $self->get_nz(m, false, casadi::Matrix<casadi_int>(kk));
      return m;
    }
    void set_nz(const std::vector<casadi_int>& kk, const Type& m) { $self->set_nz(m, false, casadi::Matrix<casadi_int>(kk)); }
    Type get_nz(const casadi::Slice& kk) const {
      Type m;
      $self->get_nz(m, false, kk);
      return m;
    }
    Type get_nz(const casadi::Matrix<casadi_int>& kk) const {
      Type m;
      $self->get_nz(m, false, kk);
      return m;
    }
    Type get_nz(casadi_int kk) const {
      Type m;
      $self->get_nz(m, false, casadi::Matrix<casadi_int>(kk));
      return m;
    }
    void set_nz(const casadi::Slice& kk, const Type& m) { $self->set_nz(m, false, kk); }
    void set_nz(const casadi::Matrix<casadi_int>& kk, const Type& m) { $self->set_nz(m, false, kk); }
    void set_nz(casadi_int kk, const Type& m) { $self->set_nz(m, false, casadi::Matrix<casadi_int>(kk)); }
%enddef
