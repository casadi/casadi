// Docstrings keep the \lang groups of the target language
%feature("docstring:lang", "matlab");

// Prototypes in MATLAB syntax
%feature("customdoc:proto:constructor", "new_obj = $NAME($in)");
%feature("customdoc:proto:single_out", "$out = $NAME($in)");
%feature("customdoc:proto:normal", "[$out] = $NAME($in)");
%feature("customdoc:main", "    $NAME $brief\n\n$overview\n$main");

// mxCreateString: '' is returned as 0x0
%fragment(SWIG_From_frag(std::string), "header") {
SWIGINTERNINLINE mxArray* SWIG_From_dec(std::string)(const std::string& s) {
  return mxCreateString(s.c_str());
}
}
// casadi_int is a double in MATLAB
%fragment(SWIG_Traits_frag(long long), "header", fragment=SWIG_AsVal_frag(long long), fragment="StdTraits") {
namespace swig {
  template <> struct traits< long long > {
    typedef value_category category;
    static const char* type_name() { return "long long"; }
  };
  template <> struct traits_asval< long long > {
    static int asval(mxArray* obj, long long* val) { return SWIG_AsVal(long long)(obj, val); }
  };
  template <> struct traits_from< long long > {
    static mxArray* from(const long long& val) { return mxCreateDoubleScalar(static_cast<double>(val)); }
  };
}
}

%{
  namespace casadi {
    // Redirect printout to mexPrintf
    static void mexlogger(const char* s, std::streamsize num, bool error) {
      if (!casadi::InterruptHandler::is_main_thread()) {
        casadi::Logger::writeDefault(s, num, error);
        return;
      }
      mexPrintf("%.*s", static_cast<int>(num), s);
    }

#ifdef HAVE_OCTAVE
    // Flush the command window buffer (needed in gui mode)
    static void mexflush(bool error) {
    }
    // Never for Octave
    static bool mexcheckinterrupted() {
      return false;
    }
    void mexclearinterrupted() {

    }
#else
    // Undocumented matlab feature
    extern "C" bool utIsInterruptPending(void);
    extern "C" void utSetInterruptPending(bool);

    static bool mexcheckinterrupted() {
      if (!casadi::InterruptHandler::is_main_thread()) return false;
      return utIsInterruptPending();
    }

    void mexclearinterrupted() {
      utSetInterruptPending(false);
    }

    // Flush the command window buffer (needed in gui mode)
    static void mexflush(bool error) {
      if (!casadi::InterruptHandler::is_main_thread()) {
        casadi::Logger::flushDefault(error);
        return;
      }
      if (!mexcheckinterrupted()) {
        if (mexEvalString("drawnow('update');pause(0.0001);")) {
          utSetInterruptPending(true);
        }
      }
    }
#endif

  }
%}
%init %{
  // Get full path
  mxArray *fullpath, *fullpath_cmd = mxCreateString("fullpath");
  mexCallMATLAB(1, &fullpath, 1, &fullpath_cmd, "mfilename");
  mxDestroyArray(fullpath_cmd);
  std::string path = mxArrayToString(fullpath);
  mxDestroyArray(fullpath);

  // Get file separator
  mxArray *filesep;
  mexCallMATLAB(1, &filesep, 0, 0, "filesep");
  std::string sep = mxArrayToString(filesep);
  mxDestroyArray(filesep);

  // Truncate at separator
  path = path.substr(0, path.rfind(sep));

  // Strip a superfluous +casadi (Octave on Windows)
  if (path.rfind(sep)!=std::string::npos && path.substr(path.rfind(sep)+1)=="+casadi")
    path = path.substr(0, path.rfind(sep));

  // Set library path
  casadi::GlobalOptions::setCasadiPath(path);
  casadi::GlobalOptions::setCasadiIncludePath(path+sep+"include");

  // Matlab is index-one based
  casadi::GlobalOptions::start_index = 1;

  // @jgillis: please document
  mxArray *warning_rhs[] = {mxCreateString("error"),

                            mxCreateString("SWIG:OverloadError")};
  mexCallMATLAB(0, 0, 2, warning_rhs, "warning");
  mxDestroyArray(warning_rhs[0]);
  mxDestroyArray(warning_rhs[1]);


  // Set logger functions
  casadi::Logger::writeFun = casadi::mexlogger;
  casadi::Logger::flush = casadi::mexflush;

  // @jgillis: please document
  casadi::InterruptHandler::checkInterrupted = casadi::mexcheckinterrupted;
  casadi::InterruptHandler::clearInterrupted = casadi::mexclearinterrupted;

  casadi::InterruptHandler::is_main_thread();

%}

%fragment("casadi_containers", "header", fragment="casadi_aux", fragment="StdVectorTraits",
          fragment="casadi_map_traits", fragment="StdPairTraits", fragment=SWIG_Traits_frag(bool),
          fragment=SWIG_Traits_frag(long long), fragment=SWIG_Traits_frag(double),
          fragment=SWIG_Traits_frag(std::string)) {
  namespace swig {
    // Containers without a %template'd traits fragment (e.g. std::vector<Dict>)
    template<class T, class A> struct traits<std::vector<T, A> > {
      typedef pointer_category category;
      static const char* type_name() { return "std::vector"; }
    };
    template<class K, class T, class C, class A> struct traits<std::map<K, T, C, A> > {
      typedef pointer_category category;
      static const char* type_name() { return "std::map"; }
    };
    template<class T, class U> struct traits<std::pair<T, U> > {
      typedef pointer_category category;
      static const char* type_name() { return "std::pair"; }
    };

    inline bool casadi_is_vector(mxArray* p) {
      return mxGetNumberOfDimensions(p)==2 && (mxGetM(p)<=1 || mxGetN(p)<=1);
    }

    // Cells only: native arrays are handled by the specializations below
    template<class Seq> inline int casadi_asptr_cell(mxArray* p, Seq** m) {
      if (!p || !mxIsCell(p)) return SWIG_ERROR;
      return traits_asptr_stdseq<Seq>::asptr(p, m);
    }

    // Double vector or cell
    template<> struct traits_asptr<std::vector<double> > {
      static int asptr(mxArray* p, std::vector<double>** m) {
        if (p && mxIsDouble(p) && casadi_is_vector(p)) {
          const double* d = static_cast<const double*>(mxGetData(p));
          if (m) *m = new std::vector<double>(d, d + mxGetM(p)*mxGetN(p));
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return casadi_asptr_cell(p, m);
      }
    };

    // Integral double vector, non-scalar logical vector or cell
    template<> struct traits_asptr<std::vector<casadi_int> > {
      static int asptr(mxArray* p, std::vector<casadi_int>** m) {
        if (p && mxIsDouble(p) && casadi_is_vector(p)) {
          const double* d = static_cast<const double*>(mxGetData(p));
          size_t n = mxGetM(p)*mxGetN(p);
          for (size_t i=0; i<n; ++i) {
            if (d[i]!=static_cast<casadi_int>(d[i])) return casadi_asptr_cell(p, m);
          }
          if (m) *m = new std::vector<casadi_int>(d, d + n);
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        if (p && mxIsLogical(p) && !mxIsLogicalScalar(p) && casadi_is_vector(p)) {
          const mxLogical* d = static_cast<const mxLogical*>(mxGetData(p));
          if (m) *m = new std::vector<casadi_int>(d, d + mxGetM(p)*mxGetN(p));
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return casadi_asptr_cell(p, m);
      }
    };

    // Double vector of zeros and ones, non-scalar logical vector or cell
    template<> struct traits_asptr<std::vector<bool> > {
      static int asptr(mxArray* p, std::vector<bool>** m) {
        if (p && mxIsDouble(p) && casadi_is_vector(p)) {
          const double* d = static_cast<const double*>(mxGetData(p));
          size_t n = mxGetM(p)*mxGetN(p);
          for (size_t i=0; i<n; ++i) {
            if (d[i]!=0 && d[i]!=1) return casadi_asptr_cell(p, m);
          }
          if (m) *m = new std::vector<bool>(d, d + n);
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        if (p && mxIsLogical(p) && !mxIsLogicalScalar(p) && casadi_is_vector(p)) {
          const mxLogical* d = static_cast<const mxLogical*>(mxGetData(p));
          if (m) *m = new std::vector<bool>(d, d + mxGetM(p)*mxGetN(p));
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return casadi_asptr_cell(p, m);
      }
    };

    // Returned as a char matrix, one string per row
    template<> struct traits_from<std::vector<std::string> > {
      static mxArray* from(const std::vector<std::string>& v) {
        std::vector<const char*> s(v.size());
        for (size_t i=0; i<s.size(); ++i) s[i] = v[i].c_str();
        return mxCreateCharMatrixFromStrings(s.size(), s.empty() ? 0 : &s[0]);
      }
    };

    // Cell of vectors, else any cell with one row per inner vector
    template<class T> struct traits_asptr<std::vector<std::vector<T> > > {
      typedef std::vector<std::vector<T> > VV;
      static int asptr(mxArray* p, VV** m) {
        int res = traits_asptr_stdseq<VV>::asptr(p, m);
        if (SWIG_IsOK(res) || !p || !mxIsCell(p)) return res;
        size_t nrow = mxGetM(p), ncol = mxGetN(p);
        VV* v = m ? new VV(nrow, std::vector<T>(ncol)) : 0;
        for (size_t i=0; i<nrow*ncol; ++i) {
          mxArray* e = mxGetCell(p, i);
          T t;
          if (!e || !SWIG_IsOK(swig::asval<T>(e, v ? &t : 0))) {
            delete v;
            return SWIG_ERROR;
          }
          if (v) (*v)[i % nrow][i / nrow] = t;
        }
        if (!m) return SWIG_OK;
        *m = v;
        return SWIG_NEWOBJ;
      }
    };

    // (casadi_int, casadi_int) is a 1-by-2 double matrix with integral entries
    template<> struct traits_asptr<std::pair<casadi_int, casadi_int> > {
      static int asptr(mxArray* p, std::pair<casadi_int, casadi_int>** m) {
        if (p && mxIsDouble(p) && mxGetNumberOfDimensions(p)==2 && !mxIsSparse(p)
            && mxGetM(p)==1 && mxGetN(p)==2) {
          const double* d = static_cast<const double*>(mxGetData(p));
          casadi_int first = static_cast<casadi_int>(d[0]);
          casadi_int second = static_cast<casadi_int>(d[1]);
          if (d[0]!=first || d[1]!=second) return SWIG_ERROR;
          if (m) *m = new std::pair<casadi_int, casadi_int>(first, second);
          return m ? SWIG_NEWOBJ : SWIG_OK;
        }
        return SWIG_ERROR;
      }
    };
    template<> struct traits_from<std::pair<casadi_int, casadi_int> > {
      static mxArray* from(const std::pair<casadi_int, casadi_int>& v) {
        mxArray* ret = mxCreateDoubleMatrix(1, 2, mxREAL);
        double* d = static_cast<double*>(mxGetData(ret));
        d[0] = v.first;
        d[1] = v.second;
        return ret;
      }
    };
  } // namespace swig
}

// Host conversions: is_null, full/sparse, casadi_host_conv overloads
%fragment("casadi_host", "header", fragment="casadi_decl") {
  namespace casadi {
    bool is_null(GUESTOBJECT *p) {
      if (p == 0) return true;
      return false;
    }

    GUESTOBJECT* full(const DM& m, bool simplify) {
      mxArray *p  = mxCreateDoubleMatrix(m.size1(), m.size2(), mxREAL);
      double* d = static_cast<double*>(mxGetData(p));
      casadi_densify(m.ptr(), m.sparsity(), d, false); // Column-major
      return p;
    }

    // Convert to a sparse matrix
    GUESTOBJECT* sparse(const DM& m) {
      mxArray *p  = mxCreateSparse(m.size1(), m.size2(), m.nnz(), mxREAL);
      casadi::casadi_copy(m.ptr(), m.nnz(), static_cast<double*>(mxGetData(p)));
      std::copy(m.colind(), m.colind()+m.size2()+1, mxGetJc(p));
      std::copy(m.row(), m.row()+m.nnz(), mxGetIr(p));
      return p;
    }

    Sparsity get_sparsity(const mxArray* p) {
      // Get sparsity pattern
      size_t nrow = mxGetM(p);
      size_t ncol = mxGetN(p);

      if (mxIsSparse(p)) {
        // Sparse storage in MATLAB
        mwIndex *Jc = mxGetJc(p);
        mwIndex *Ir = mxGetIr(p);

        // Store in vectors
        std::vector<casadi_int> colind(ncol+1);
        std::copy(Jc, Jc+colind.size(), colind.begin());
        std::vector<casadi_int> row(colind.back());
        std::copy(Ir, Ir+row.size(), row.begin());

        // Create pattern and return
        return Sparsity(nrow, ncol, colind, row);
      } else {
        return Sparsity::dense(nrow, ncol);
      }
    }

    size_t getNNZ(const mxArray* p) {
      // Dimensions
      size_t nrow = mxGetM(p);
      size_t ncol = mxGetN(p);
      if (mxIsSparse(p)) {
        // Sparse storage in MATLAB
        mwIndex *Jc = mxGetJc(p);
        return Jc[ncol];
      } else {
        return nrow*ncol;
      }
    }
  } // namespace casadi

  namespace swig {
    inline int casadi_host_conv(GUESTOBJECT* p, casadi::DM* m) {
      using namespace casadi;
      // MATLAB double matrix (sparse or dense)
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
        if (m) {
          *m = casadi::DM(get_sparsity(p));
          double* data = static_cast<double*>(mxGetData(p));
          casadi_copy(data, m->nnz(), m->ptr());
        }
        return true;
      }

      // No match
      return false;
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::IM* m) {
      using namespace casadi;
      // In MATLAB, it is common to use floating point values to represent integers
      if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
        double* data = static_cast<double*>(mxGetData(p));

        // Check if all integers
        bool all_integers=true;
        size_t sz = getNNZ(p);
        for (size_t i=0; i<sz; ++i) {
          if (data[i] != casadi_int(data[i])) {
            all_integers = false;
            break;
          }
        }

        // If successful
        if (all_integers) {
          if (m) {
            *m = casadi::IM(get_sparsity(p));
            for (size_t i=0; i<sz; ++i) {
              (*m)->at(i) = casadi_int(data[i]);
            }
          }
          return true;
        }
      }

      return -1;
    }

    inline GUESTOBJECT* casadi_host_from(const casadi::GenericType& a) {
      using namespace casadi;
      switch (a.getType()) {
      case OT_VOIDPTR:
        return mxCreateDoubleScalar(0);
      default: return 0;
      }
    }
  } // namespace swig
}

// Type names in error messages and docs
  #define L_INT "int"
  #define L_BOOL "bool"
  #define LPAIR(A,B) "(" A "," B ")"
  #define L_DOUBLE "double"
  #define L_DICT "struct"
  #define LDICT(ARG) L_DICT ":" ARG
  #define LL "{"
  #define LR "}"
  #define L_STR "char"
  #define LABEL "[int,int]"

// MATLAB two-index paren overloads, with ':' passed as char
%define %_paren2(Type, RT, R, CT, C)
    const Type paren(RT rr, CT cc) const { Type m; $self->get(m, true, R, C); return m;}
%enddef
%define %_paren_asgn2(Type, RT, R, CT, C)
    void paren_asgn(const Type& m, RT rr, CT cc) { $self->set(m, true, R, C);}
%enddef
%define %_paren_combos(Macro, Type, I)
    Macro(Type, char, casadi::char2Slice(rr), I, cc)
    Macro(Type, I, rr, char, casadi::char2Slice(cc))
    Macro(Type, I, rr, I, cc)
%enddef

// Indexing helpers of the matrix classes
%define %matrix_helpers(Type, xIndex)
    // Get a submatrix (index-1)
    const Type paren(char rr) const {
      casadi_assert_dev(rr==':');
      return vec(*$self);
    }
    const Type paren(const Matrix<casadi_int>& rr) const {
      Type m;
      $self->get(m, true, rr);
      return m;
    }
    const Type paren(const Sparsity& sp) const {
      Type m;
      $self->get(m, true, sp);
      return m;
    }
    %_paren2(Type, char, casadi::char2Slice(rr), char, casadi::char2Slice(cc))
    %_paren_combos(%_paren2, Type, const Matrix<casadi_int>&)

    // Set a submatrix (index-1)
    void paren_asgn(const Type& m, char rr) {
      casadi_assert_dev(rr==':');
      $self->set(m, false, casadi::IM(casadi::range($self->numel())));
    }
    void paren_asgn(const Type& m, const Matrix<casadi_int>& rr) { $self->set(m, true, rr);}
    void paren_asgn(const Type& m, const Sparsity& sp) { $self->set(m, true, sp);}
    %_paren_asgn2(Type, char, casadi::char2Slice(rr), char, casadi::char2Slice(cc))
    %_paren_combos(%_paren_asgn2, Type, const Matrix<casadi_int>&)

    // Get nonzeros (index-1)
    const Type brace(char rr) const { Type m; $self->get_nz(m, true, casadi::char2Slice(rr)); return m;}
    const Type brace(const Matrix<casadi_int>& rr) const { Type m; $self->get_nz(m, true, rr); return m;}

    // Set nonzeros (index-1)
    void setbrace(const Type& m, char rr) { $self->set_nz(m, true, casadi::char2Slice(rr));}
    void setbrace(const Type& m, const Matrix<casadi_int>& rr) { $self->set_nz(m, true, rr);}

    // 'end' function (needed for end syntax in MATLAB)
    inline casadi_int end(casadi_int i, casadi_int n) const {
      return n==1 ? $self->numel() : i==1 ? $self->size1() : $self->size2();
    }


    // Needed for brace syntax to access nonzeros
    casadi_int numel(casadi_int k) const {
      return 1;
    }

    // Needed for brace syntax to access nonzeros
    casadi_int numel(char rr) const {
      casadi_assert_dev(rr==':');
      return 1;
    }

    // Needed for brace syntax to access nonzeros
    casadi_int numel(const std::vector<casadi_int> &k) const {
      return 1;
    }

    // Needed because original numel call gets hidden by the above extend overloads
    casadi_int numel() const {
      return $self->numel();
    }


    // Transpose using the A' syntax in addition to A.'
    Type ctranspose() const { return $self->T();}

%enddef

// Sparsity interface functions in MATLAB form
%define CASADI_TARGET_SPARSITY_INTERFACE_FUN(DECL, FLAG, M)
#if FLAG & IS_GLOBAL
 // Package functions for calls without a casadi object first, e.g. offset({...}) or diagcat()
 DECL std::vector<casadi_int >
 casadi_offset(const std::vector< M > &v, bool vert=true) {
 return offset(v, vert);
 }
 DECL M casadi_blockcat(const M& A, const M& B, const M& C, const M& D) {
 return vertcat(horzcat(A, B), horzcat(C, D));
 }
 DECL M casadi_diagcat(const std::vector< M > &A) {
 return diagcat(A);
 }
 DECL M casadi_veccat(const std::vector< M >& x) {
 return veccat(x);
 }
#endif
#if FLAG & IS_MEMBER
 // [] for a dimension: deduce it, as reshape(x, [], 1) in MATLAB
 DECL M casadi_reshape(const M& a, const std::vector<casadi_int>& nrow, casadi_int ncol) {
 casadi_assert(nrow.size()<=1, "reshape: expected a scalar or [] for the number of rows");
 return reshape(a, nrow.empty() ? -1 : nrow[0], ncol);
 }
 DECL M casadi_reshape(const M& a, casadi_int nrow, const std::vector<casadi_int>& ncol) {
 casadi_assert(ncol.size()<=1, "reshape: expected a scalar or [] for the number of columns");
 return reshape(a, nrow, ncol.empty() ? -1 : ncol[0]);
 }
#endif
%enddef

// length and 1-based sum(A, dim)
%define CASADI_TARGET_SUM_FUN(DECL, M)
      DECL casadi_int casadi_length(const M &v) {
        if (v.is_empty()) return 0;
        return std::max(v.size1(), v.size2());
      }
      DECL M casadi_sum(const M& x, casadi_int dim) {
        if (dim==1) return sum1(x);
        if (dim==2) return sum2(x);
        casadi_error(
          "Expected sum(A,1), sum(A,2), sum(A,\"all\") got " + casadi::str(dim) + " instead.");
      }
      DECL M casadi_sum(const M& x, const std::string& dim) {
        casadi_assert(dim=="all",
          "Expected sum(...,'all'), got '" + dim + "' instead.");
        return sum(x);
      }
      DECL M casadi_sum(const M& x) {
        if (x.is_vector()) return sum(x);
        return sum1(x);
      }
%enddef

// Member functions as static functions of the mixin classes
%define CASADI_TARGET_MEMBER_FUNS
// Wrap (static) member functions
%feature("nonstatic");
namespace casadi {
  %extend SparsityInterfaceCommon {
    SPARSITY_INTERFACE_ALL(static inline, IS_MEMBER)
  }
  %extend GenericExpressionCommon {
    GENERIC_EXPRESSION_ALL(static inline, IS_MEMBER)
  }
  %extend GenericMatrixCommon {
    GENERIC_MATRIX_ALL(static inline, IS_MEMBER)
  }
  %extend MatrixCommon {
    MATRIX_ALL(static inline, IS_MEMBER)
  }
  %extend MX {
    MX_ALL(static inline, IS_MEMBER)
    const MX brace(const casadi::MX& rr) const { casadi::MX m; $self->get_nz(m, true, rr); return m;}
    void brace_asgn(const MX& m, const casadi::MX& rr) { $self->set_nz(m, true, rr); }
    const MX paren(const casadi::MX& rr) const {
      casadi::MX m;
      $self->get(m, true, rr);
      return m;
    }
    %_paren_combos(%_paren2, casadi::MX, const casadi::MX&)
    // Needed for brace syntax to access nonzeros
    casadi_int numel(const MX &k) const {
      return 1;
    }
  }
} // namespace casadi
%feature("nonstatic", "");
%enddef
