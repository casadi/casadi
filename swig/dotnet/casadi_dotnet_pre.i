// Host conversions: is_null, full/sparse, casadi_host_conv overloads
%fragment("casadi_host", "header", fragment="casadi_decl") {
  namespace casadi {
    bool is_null(GUESTOBJECT *p) {
      return !p || p->kind == SWIG_RV_NIL;
    }

    GUESTOBJECT* full(const DM& m, bool simplify) {
      SwigRV *p = SWIG_Dotnet_NewNum(SWIG_RV_F64, m.size1(), m.size2());
      if (m.numel()) casadi_densify(m.ptr(), m.sparsity(), &p->f64[0], false);
      return p;
    }

    // Convert to a sparse matrix
    GUESTOBJECT* sparse(const DM& m) {
      SwigRV *p = new SwigRV();
      p->kind = SWIG_RV_NUM;
      p->rows = m.size1();
      p->cols = m.size2();
      p->sparse = true;
      p->colind.assign(m.colind(), m.colind()+m.size2()+1);
      p->row.assign(m.row(), m.row()+m.nnz());
      p->f64 = m.nonzeros();
      return p;
    }

    Sparsity get_sparsity(const SwigRV* p) {
      if (!p->sparse) return Sparsity::dense(p->rows, p->cols);
      return Sparsity(p->rows, p->cols, std::vector<casadi_int>(p->colind.begin(), p->colind.end()),
                      std::vector<casadi_int>(p->row.begin(), p->row.end()));
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
      using namespace casadi;
      // Numeric matrix (dense or sparse)
      if (SWIG_Dotnet_IsNum(p)) {
        if (m) {
          *m = DM(get_sparsity(p));
          if (p->elem==SWIG_RV_F64) std::copy(p->f64.begin(), p->f64.end(), m->ptr());
          else std::copy(p->i64.begin(), p->i64.end(), m->ptr());
        }
        return true;
      }
      if (p->kind!=SWIG_RV_LIST) return false;

      // List of equally long numeric rows
      if (!p->items.empty() && SWIG_Dotnet_IsVector(p->items[0]) && p->items[0]->numel()!=1) {
        casadi_int ncol = p->items[0]->numel();
        std::vector<double> rows, r;
        for (SwigRV* e : p->items) {
          if (!SWIG_Dotnet_IsVector(e) || e->numel()!=ncol || !swig::assign(e, &r)) return false;
          rows.insert(rows.end(), r.begin(), r.end());
        }
        if (m) *m = DM::reshape(DM(rows), ncol, p->items.size()).T();
        return true;
      }

      // List of numbers: a column
      std::vector<double> col;
      if (!swig::assign(p, m ? &col : 0)) return false;
      if (m) *m = DM(col);
      return true;
    }

    inline int casadi_host_conv(GUESTOBJECT* p, casadi::IM* m) {
      using namespace casadi;
      // Integer or bool matrix (dense or sparse)
      if (SWIG_Dotnet_IsNum(p) && p->elem!=SWIG_RV_F64) {
        if (m) {
          *m = IM(get_sparsity(p));
          std::copy(p->i64.begin(), p->i64.end(), m->ptr());
        }
        return true;
      }

      return -1;
    }

    inline GUESTOBJECT* casadi_host_from(const casadi::GenericType& a) {
      using namespace casadi;
      switch (a.getType()) {
      case OT_NULL:
      case OT_VOIDPTR:
        return new SwigRV();
      default: return 0;
      }
    }
  } // namespace swig
}

// Docstrings keep the \lang groups of the target language
%feature("docstring:lang", "csharp");

// Indexing helpers of the matrix classes: get/set by Slice, IM, Sparsity and integers (0-based),
// behind the C# indexers of casadi_dotnet_decl.i
%define %dotnet_index_pair(Type, R, C)
    Type get(const R& rr, const C& cc) const {
      Type m;
      $self->get(m, false, rr, cc);
      return m;
    }
    void set(const R& rr, const C& cc, const Type& m) { $self->set(m, false, rr, cc); }
%enddef
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
    %dotnet_index_pair(Type, casadi::Slice, casadi::Slice)
    %dotnet_index_pair(Type, casadi::Slice, casadi::Matrix<casadi_int>)
    %dotnet_index_pair(Type, casadi::Matrix<casadi_int>, casadi::Slice)
    %dotnet_index_pair(Type, casadi::Matrix<casadi_int>, casadi::Matrix<casadi_int>)
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
    void set_nz(const casadi::Slice& kk, const Type& m) { $self->set_nz(m, false, kk); }
    void set_nz(const casadi::Matrix<casadi_int>& kk, const Type& m) { $self->set_nz(m, false, kk); }
%enddef

// Zero-argument getters as C# properties: opti.f, opti.x, sol.opti
%define make_property(class, name, type)
  %feature("dotnet:property") class ## :: ## name;
%enddef

// Type names in error messages and docs
  #define L_INT "long"
  #define L_BOOL "bool"
  #define LPAIR(A,B) "(" A "," B ")"
  #define L_DOUBLE "double"
  #define L_DICT "IDictionary<string,object>"
  #define LDICT(ARG) "IDictionary<string," ARG ">"
  #define LL "IEnumerable<"
  #define LR ">"
  #define L_STR "string"
  #define LABEL LPAIR(L_INT, L_INT)
