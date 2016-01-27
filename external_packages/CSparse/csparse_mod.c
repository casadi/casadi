#include "csparse_mod.h"

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }

/* C = alpha*A + beta*B */
void cs_add(cs *C, const cs *A, const cs *B, double alpha, double beta) {
  int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
  double *x, *Bx, *Cx ;
  m = A->m;
  anz = A->p[A->n] ;
  n = B->n;
  Bp = B->p;
  Bx = B->x;
  bnz = Bp[n] ;
  w = cs_calloc (m, sizeof (int)) ;                       /* get workspace */
  values = (A->x != NULL) && (Bx != NULL) ;
  x = values ? cs_malloc (m, sizeof (double)) : NULL ;    /* get workspace */
  cs_spalloc(C, m, n, anz + bnz, values);
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (j = 0 ; j < n ; j++) {
    Cp [j] = nz ;                   /* column j of C starts here */
    nz = cs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
    nz = cs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
    if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
  }
  Cp [n] = nz ;                       /* finalize the last column of C */
  cs_sprealloc (C, 0) ;               /* remove extra space from C */
  cs_free (w) ;                       /* free workspace */
  cs_free (x) ;
}

/* clear w */
static int cs_wclear (int mark, int lemax, int *w, int n) {
  int k ;
  if (mark < 2 || (mark + lemax < 0)) {
    for (k = 0 ; k < n ; k++) if (w [k] != 0) w [k] = 1 ;
    mark = 2 ;
  }
  return (mark) ;     /* at this point, w [0..n-1] < mark holds */
}

/* keep off-diagonal entries; drop diagonal entries */
static int cs_diag (int i, int j, double aij, void *other) { return (i != j) ; }

/* p = amd(A+A') if symmetric is true, or amd(A'A) otherwise */
/* order 0:natural, 1:Chol, 2:LU, 3:QR */
int *cs_amd (int order, const cs *A) {
  cs *C, *A2, *AT ;
  int *Cp, *Ci, *last, *W, *len, *nv, *next, *P, *head, *elen, *degree, *w,
    *hhead, *ATp, *ATi, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
    k2, k3, jlast, ln, dense, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
    ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, n, m, t ;
  int h ;
  /* --- Construct matrix C ----------------------------------------------- */
  AT = cs_calloc(1, sizeof (cs));
  cs_transpose (A, AT, 0) ;              /* compute A' */
  m = A->m ; n = A->n ;
  dense = CS_MAX (16, 10 * sqrt ((double) n)) ;   /* find dense threshold */
  dense = CS_MIN (n-2, dense) ;
  C = cs_calloc(1, sizeof (cs));
  if (order == 1 && n == m) {
    cs_add (C, A, AT, 0, 0);          /* C = A+A' */
  } else if (order == 2) {
    ATp = AT->p ;                       /* drop dense columns from AT */
    ATi = AT->i ;
    for (p2 = 0, j = 0 ; j < m ; j++) {
      p = ATp [j] ;                   /* column j of AT starts here */
      ATp [j] = p2 ;                  /* new column j starts here */
      if (ATp [j+1] - p > dense) continue ;   /* skip dense col j */
      for ( ; p < ATp [j+1] ; p++) ATi [p2++] = ATi [p] ;
    }
    ATp [m] = p2 ;                      /* finalize AT */
    A2 = cs_calloc(1, sizeof (cs));
    cs_transpose(AT, A2, 0) ;         /* A2 = AT' */
    cs_multiply (C, AT, A2); /* C=A'*A with no dense rows */
    cs_spfree (A2) ;
  } else {
    cs_multiply(C, AT, A) ;           /* C=A'*A */
  }
  cs_spfree (AT) ;
  cs_fkeep (C, &cs_diag, NULL) ;          /* drop diagonal entries */
  Cp = C->p ;
  cnz = Cp [n] ;
  P = cs_malloc (n+1, sizeof (int)) ;     /* allocate result */
  W = cs_malloc (8*(n+1), sizeof (int)) ; /* get workspace */
  t = cnz + cnz/5 + 2*n ;                 /* add elbow room to C */
  cs_sprealloc (C, t);
  len  = W           ; nv     = W +   (n+1) ; next   = W + 2*(n+1) ;
  head = W + 3*(n+1) ; elen   = W + 4*(n+1) ; degree = W + 5*(n+1) ;
  w    = W + 6*(n+1) ; hhead  = W + 7*(n+1) ;
  last = P ;                              /* use P as workspace for last */
  /* --- Initialize quotient graph ---------------------------------------- */
  for (k = 0 ; k < n ; k++) len [k] = Cp [k+1] - Cp [k] ;
  len [n] = 0 ;
  nzmax = C->nzmax ;
  Ci = C->i ;
  for (i = 0 ; i <= n ; i++) {
    head [i] = -1 ;                     /* degree list i is empty */
    last [i] = -1 ;
    next [i] = -1 ;
    hhead [i] = -1 ;                    /* hash list i is empty */
    nv [i] = 1 ;                        /* node i is just one node */
    w [i] = 1 ;                         /* node i is alive */
    elen [i] = 0 ;                      /* Ek of node i is empty */
    degree [i] = len [i] ;              /* degree of node i */
  }
  mark = cs_wclear (0, 0, w, n) ;         /* clear w */
  elen [n] = -2 ;                         /* n is a dead element */
  Cp [n] = -1 ;                           /* n is a root of assembly tree */
  w [n] = 0 ;                             /* n is a dead element */
  /* --- Initialize degree lists ------------------------------------------ */
  for (i = 0 ; i < n ; i++) {
    d = degree [i] ;
    if (d == 0) {
      /* node i is empty */
      elen [i] = -2 ;                 /* element i is dead */
      nel++ ;
      Cp [i] = -1 ;                   /* i is a root of assembly tree */
      w [i] = 0 ;
    } else if (d > dense) {
      /* node i is dense */
      nv [i] = 0 ;                    /* absorb i into element n */
      elen [i] = -1 ;                 /* node i is dead */
      nel++ ;
      Cp [i] = CS_FLIP (n) ;
      nv [n]++ ;
    } else {
      if (head [d] != -1) last [head [d]] = i ;
      next [i] = head [d] ;           /* put node i in degree list d */
      head [d] = i ;
    }
  }

  /* while (selecting pivots) do */
  while (nel < n) {
    /* --- Select node of minimum approximate degree -------------------- */
    for (k = -1 ; mindeg < n && (k = head [mindeg]) == -1 ; mindeg++) ;
    if (next [k] != -1) last [next [k]] = -1 ;
    head [mindeg] = next [k] ;          /* remove k from degree list */
    elenk = elen [k] ;                  /* elenk = |Ek| */
    nvk = nv [k] ;                      /* # of nodes k represents */
    nel += nvk ;                        /* nv[k] nodes of A eliminated */
    /* --- Garbage collection ------------------------------------------- */
    if (elenk > 0 && cnz + mindeg >= nzmax) {
      for (j = 0 ; j < n ; j++) {
        if ((p = Cp [j]) >= 0) {
          /* j is a live node or element */
          Cp [j] = Ci [p] ;       /* save first entry of object */
          Ci [p] = CS_FLIP (j) ;  /* first entry is now CS_FLIP(j) */
        }
      }

      /* scan all of memory */
      for (q = 0, p = 0 ; p < cnz ; ) {
        if ((j = CS_FLIP (Ci [p++])) >= 0) {
          /* found object j */
          Ci [q] = Cp [j] ;       /* restore first entry of object */
          Cp [j] = q++ ;          /* new pointer to object j */
          for (k3 = 0 ; k3 < len [j]-1 ; k3++) Ci [q++] = Ci [p++] ;
        }
      }
      cnz = q ;                       /* Ci [cnz...nzmax-1] now free */
    }
    /* --- Construct new element ---------------------------------------- */
    dk = 0 ;
    nv [k] = -nvk ;                     /* flag k as in Lk */
    p = Cp [k] ;
    pk1 = (elenk == 0) ? p : cnz ;      /* do in place if elen[k] == 0 */
    pk2 = pk1 ;
    for (k1 = 1 ; k1 <= elenk + 1 ; k1++) {
      if (k1 > elenk) {
        e = k ;                     /* search the nodes in k */
        pj = p ;                    /* list of nodes starts at Ci[pj]*/
        ln = len [k] - elenk ;      /* length of list of nodes in k */
      } else {
        e = Ci [p++] ;              /* search the nodes in e */
        pj = Cp [e] ;
        ln = len [e] ;              /* length of list of nodes in e */
      }
      for (k2 = 1 ; k2 <= ln ; k2++) {
        i = Ci [pj++] ;
        if ((nvi = nv [i]) <= 0) continue ; /* node i dead, or seen */
        dk += nvi ;                 /* degree[Lk] += size of node i */
        nv [i] = -nvi ;             /* negate nv[i] to denote i in Lk*/
        Ci [pk2++] = i ;            /* place i in Lk */
        if (next [i] != -1) last [next [i]] = last [i] ;
        if (last [i] != -1) {
          /* remove i from degree list */
          next [last [i]] = next [i] ;
        } else {
          head [degree [i]] = next [i] ;
        }
      }
      if (e != k) {
        Cp [e] = CS_FLIP (k) ;      /* absorb e into k */
        w [e] = 0 ;                 /* e is now a dead element */
      }
    }
    if (elenk != 0) cnz = pk2 ;         /* Ci [cnz...nzmax] is free */
    degree [k] = dk ;                   /* external degree of k - |Lk\i| */
    Cp [k] = pk1 ;                      /* element k is in Ci[pk1..pk2-1] */
    len [k] = pk2 - pk1 ;
    elen [k] = -2 ;                     /* k is now an element */
    /* --- Find set differences ----------------------------------------- */
    mark = cs_wclear (mark, lemax, w, n) ;  /* clear w if necessary */
    for (pk = pk1 ; pk < pk2 ; pk++) {
      /* scan 1: find |Le\Lk| */
      i = Ci [pk] ;
      if ((eln = elen [i]) <= 0) continue ;/* skip if elen[i] empty */
      nvi = -nv [i] ;                      /* nv [i] was negated */
      wnvi = mark - nvi ;
      /* scan Ei */
      for (p = Cp [i] ; p <= Cp [i] + eln - 1 ; p++) {
        e = Ci [p] ;
        if (w [e] >= mark) {
          w [e] -= nvi ;          /* decrement |Le\Lk| */
        } else if (w [e] != 0) {
          /* ensure e is a live element */
          w [e] = degree [e] + wnvi ; /* 1st time e seen in scan 1 */
        }
      }
    }
    /* --- Degree update ------------------------------------------------ */
    /* scan2: degree update */
    for (pk = pk1 ; pk < pk2 ; pk++) {
      i = Ci [pk] ;                   /* consider node i in Lk */
      p1 = Cp [i] ;
      p2 = p1 + elen [i] - 1 ;
      pn = p1 ;
      /* scan Ei */
      for (h = 0, d = 0, p = p1 ; p <= p2 ; p++) {
        e = Ci [p] ;
        /* e is an unabsorbed element */
        if (w [e] != 0) {
          dext = w [e] - mark ;   /* dext = |Le\Lk| */
          if (dext > 0)
            {
              d += dext ;         /* sum up the set differences */
              Ci [pn++] = e ;     /* keep e in Ei */
              h += e ;            /* compute the hash of node i */
            }  else {
            Cp [e] = CS_FLIP (k) ;  /* aggressive absorb. e->k */
            w [e] = 0 ;             /* e is a dead element */
          }
        }
      }
      elen [i] = pn - p1 + 1 ;        /* elen[i] = |Ei| */
      p3 = pn ;
      p4 = p1 + len [i] ;
      /* prune edges in Ai */
      for (p = p2 + 1 ; p < p4 ; p++) {
        j = Ci [p] ;
        if ((nvj = nv [j]) <= 0) continue ; /* node j dead or in Lk */
        d += nvj ;                  /* degree(i) += |j| */
        Ci [pn++] = j ;             /* place j in node list of i */
        h += j ;                    /* compute hash for node i */
      }

      /* check for mass elimination */
      if (d == 0) {
        Cp [i] = CS_FLIP (k) ;      /* absorb i into k */
        nvi = -nv [i] ;
        dk -= nvi ;                 /* |Lk| -= |i| */
        nvk += nvi ;                /* |k| += nv[i] */
        nel += nvi ;
        nv [i] = 0 ;
        elen [i] = -1 ;             /* node i is dead */
      } else {
        degree [i] = CS_MIN (degree [i], d) ;   /* update degree(i) */
        Ci [pn] = Ci [p3] ;         /* move first node to end */
        Ci [p3] = Ci [p1] ;         /* move 1st el. to end of Ei */
        Ci [p1] = k ;               /* add k as 1st element in of Ei */
        len [i] = pn - p1 + 1 ;     /* new len of adj. list of node i */
        h = ((h<0) ? (-h):h) % n ;  /* finalize hash of i */
        next [i] = hhead [h] ;      /* place i in hash bucket */
        hhead [h] = i ;
        last [i] = h ;              /* save hash of i in last[i] */
      }
    }                                   /* scan2 is done */
    degree [k] = dk ;                   /* finalize |Lk| */
    lemax = CS_MAX (lemax, dk) ;
    mark = cs_wclear (mark+lemax, lemax, w, n) ;    /* clear w */
    /* --- Supernode detection ------------------------------------------ */
    for (pk = pk1 ; pk < pk2 ; pk++) {
      i = Ci [pk] ;
      if (nv [i] >= 0) continue ;         /* skip if i is dead */
      h = last [i] ;                      /* scan hash bucket of node i */
      i = hhead [h] ;
      hhead [h] = -1 ;                    /* hash bucket will be empty */
      for ( ; i != -1 && next [i] != -1 ; i = next [i], mark++) {
        ln = len [i] ;
        eln = elen [i] ;
        for (p = Cp [i]+1 ; p <= Cp [i] + ln-1 ; p++) w [Ci [p]] = mark;
        jlast = i ;
        for (j = next [i] ; j != -1 ; ) /* compare i with all j */
          {
            ok = (len [j] == ln) && (elen [j] == eln) ;
            for (p = Cp [j] + 1 ; ok && p <= Cp [j] + ln - 1 ; p++)
              {
                if (w [Ci [p]] != mark) ok = 0 ;    /* compare i and j*/
              }
            if (ok)                     /* i and j are identical */
              {
                Cp [j] = CS_FLIP (i) ;  /* absorb j into i */
                nv [i] += nv [j] ;
                nv [j] = 0 ;
                elen [j] = -1 ;         /* node j is dead */
                j = next [j] ;          /* delete j from hash bucket */
                next [jlast] = j ;
              }
            else
              {
                jlast = j ;             /* j and i are different */
                j = next [j] ;
              }
          }
      }
    }
    /* --- Finalize new element------------------------------------------ */
    /* finalize Lk */
    for (p = pk1, pk = pk1 ; pk < pk2 ; pk++) {
      i = Ci [pk] ;
      if ((nvi = -nv [i]) <= 0) continue ;/* skip if i is dead */
      nv [i] = nvi ;                      /* restore nv[i] */
      d = degree [i] + dk - nvi ;         /* compute external degree(i) */
      d = CS_MIN (d, n - nel - nvi) ;
      if (head [d] != -1) last [head [d]] = i ;
      next [i] = head [d] ;               /* put i back in degree list */
      last [i] = -1 ;
      head [d] = i ;
      mindeg = CS_MIN (mindeg, d) ;       /* find new minimum degree */
      degree [i] = d ;
      Ci [p++] = i ;                      /* place i in Lk */
    }
    nv [k] = nvk ;                      /* # nodes absorbed into k */
    /* length of adj list of element k*/
    if ((len [k] = p-pk1) == 0) {
      Cp [k] = -1 ;                   /* k is a root of the tree */
      w [k] = 0 ;                     /* k is now a dead element */
    }
    if (elenk != 0) cnz = p ;           /* free unused space in Lk */
  }
  /* --- Postordering ----------------------------------------------------- */
  for (i = 0 ; i < n ; i++) Cp [i] = CS_FLIP (Cp [i]) ;/* fix assembly tree */
  for (j = 0 ; j <= n ; j++) head [j] = -1 ;
  /* place unordered nodes in lists */
  for (j = n ; j >= 0 ; j--) {
    if (nv [j] > 0) continue ;          /* skip if j is an element */
    next [j] = head [Cp [j]] ;          /* place j in list of its parent */
    head [Cp [j]] = j ;
  }
  /* place elements in lists */
  for (e = n ; e >= 0 ; e--) {
    if (nv [e] <= 0) continue ;         /* skip unless e is an element */
    if (Cp [e] != -1)
      {
        next [e] = head [Cp [e]] ;      /* place e in list of its parent */
        head [Cp [e]] = e ;
      }
  }
  /* postorder the assembly tree */
  for (k = 0, i = 0 ; i <= n ; i++) {
    if (Cp [i] == -1) k = cs_tdfs (i, k, head, next, P, w) ;
  }
  cs_spfree(C);
  cs_free(W);
  return P;
}

/* L = chol (A, [pinv parent cp]), pinv is optional */
csn *cs_chol (const cs *A, const css *S) {
  double d, lki, *Lx, *x, *Cx ;
  int top, i, p, k, n, *Li, *Lp, *cp, *pinv, *s, *c, *parent, *Cp, *Ci ;
  cs *L, *C, *E ;
  csn *N ;
  n = A->n ;
  N = cs_calloc (1, sizeof (csn)) ;       /* allocate result */
  c = cs_malloc (2*n, sizeof (int)) ;     /* get int workspace */
  x = cs_malloc (n, sizeof (double)) ;    /* get double workspace */
  cp = S->cp ; pinv = S->pinv ; parent = S->parent ;
  C = pinv ? cs_symperm (A, pinv, 1) : ((cs *) A) ;
  E = pinv ? C : NULL ;           /* E is alias for A, or a copy E=A(p,p) */
  s = c + n ;
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  N->L = L = cs_calloc(1, sizeof (cs));
  cs_spalloc(L, n, n, cp [n], 1) ;    /* allocate result */
  Lp = L->p ; Li = L->i ; Lx = L->x ;
  for (k = 0 ; k < n ; k++) Lp [k] = c [k] = cp [k] ;
  /* compute L(k,:) for L*L' = C */
  for (k = 0 ; k < n ; k++) {
    /* --- Nonzero pattern of L(k,:) ------------------------------------ */
    top = cs_ereach (C, k, parent, s, c) ;      /* find pattern of L(k,:) */
    x [k] = 0 ;                                 /* x (0:k) is now zero */
    for (p = Cp [k] ; p < Cp [k+1] ; p++)       /* x = full(triu(C(:,k))) */
      {
        if (Ci [p] <= k) x [Ci [p]] = Cx [p] ;
      }
    d = x [k] ;                     /* d = C(k,k) */
    x [k] = 0 ;                     /* clear x for k+1st iteration */
    /* --- Triangular solve --------------------------------------------- */
    /* solve L(0:k-1,0:k-1) * x = C(:,k) */
    for ( ; top < n ; top++) {
      i = s [top] ;               /* s [top..n-1] is pattern of L(k,:) */
      lki = x [i] / Lx [Lp [i]] ; /* L(k,i) = x (i) / L(i,i) */
      x [i] = 0 ;                 /* clear x for k+1st iteration */
      for (p = Lp [i] + 1 ; p < c [i] ; p++) {
        x [Li [p]] -= Lx [p] * lki ;
      }
      d -= lki * lki ;            /* d = d - L(k,i)*L(k,i) */
      p = c [i]++ ;
      Li [p] = k ;                /* store L(k,i) in column i */
      Lx [p] = lki ;
    }
    /* --- Compute L(k,k) ----------------------------------------------- */
    if (d <= 0) {
      cs_spfree(E);
      cs_free (c);
      cs_free (x);
      cs_nfree(N);
      return NULL; /* not pos def */
    }
    p = c [k]++ ;
    Li [p] = k ;                /* store L(k,k) = sqrt (d) in column k */
    Lx [p] = sqrt (d) ;
  }
  Lp [n] = cp [n] ;               /* finalize L */
  cs_spfree(E);
  cs_free(c);
  cs_free(x);
  return N;
}

/* x=A\b where A is symmetric positive definite; b overwritten with solution */
int cs_cholsol (int order, const cs *A, double *b) {
  double *x ;
  css *S ;
  csn *N ;
  int n, ok ;
  n = A->n ;
  S = cs_schol (order, A) ;               /* ordering and symbolic analysis */
  N = cs_chol (A, S) ;                    /* numeric Cholesky factorization */
  x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
  ok = (S && N && x) ;
  if (ok) {
    cs_ipvec (S->pinv, b, x, n) ;   /* x = P*b */
    cs_lsolve (N->L, x) ;           /* x = L\x */
    cs_ltsolve (N->L, x) ;          /* x = L'\x */
    cs_pvec (S->pinv, x, b, n) ;    /* b = P'*x */
  }
  cs_free (x) ;
  cs_sfree (S) ;
  cs_nfree (N) ;
  return (ok) ;
}

/* column counts of LL'=A or LL'=A'A, given parent & post ordering */
#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
static void init_ata (cs *AT, const int *post, int *w, int **head, int **next) {
  int i, k, p, m = AT->n, n = AT->m, *ATp = AT->p, *ATi = AT->i ;
  *head = w+4*n, *next = w+5*n+1 ;
  for (k = 0 ; k < n ; k++) w [post [k]] = k ;    /* invert post */
  for (i = 0 ; i < m ; i++) {
    for (k = n, p = ATp[i] ; p < ATp[i+1] ; p++) k = CS_MIN (k, w [ATi[p]]);
    (*next) [i] = (*head) [k] ;     /* place row i in linked list k */
    (*head) [k] = i ;
  }
}

int *cs_counts (const cs *A, const int *parent, const int *post, int ata) {
  int i, j, k, n, m, J, s, p, q, jleaf, *ATp, *ATi, *maxfirst, *prevleaf,
    *ancestor, *head = NULL, *next = NULL, *colcount, *w, *first, *delta ;
  cs *AT ;
  m = A->m ; n = A->n ;
  s = 4*n + (ata ? (n+m+1) : 0) ;
  delta = colcount = cs_malloc (n, sizeof (int)) ;    /* allocate result */
  w = cs_malloc (s, sizeof (int)) ;                   /* get workspace */
  AT = cs_calloc(1, sizeof (cs));
  cs_transpose(A, AT, 0) ;                          /* AT = A' */
  ancestor = w ; maxfirst = w+n ; prevleaf = w+2*n ; first = w+3*n ;
  for (k = 0 ; k < s ; k++) w [k] = -1 ;      /* clear workspace w [0..s-1] */
  /* find first [j] */
  for (k = 0 ; k < n ; k++) {
    j = post [k] ;
    delta [j] = (first [j] == -1) ? 1 : 0 ;  /* delta[j]=1 if j is a leaf */
    for ( ; j != -1 && first [j] == -1 ; j = parent [j]) first [j] = k ;
  }
  ATp = AT->p ; ATi = AT->i ;
  if (ata) init_ata (AT, post, w, &head, &next) ;
  for (i = 0 ; i < n ; i++) ancestor [i] = i ; /* each node in its own set */
  for (k = 0 ; k < n ; k++) {
    j = post [k] ;          /* j is the kth node in postordered etree */
    if (parent [j] != -1) delta [parent [j]]-- ;    /* j is not a root */
    /* J=j for LL'=A case */
    for (J = HEAD (k,j) ; J != -1 ; J = NEXT (J)) {
      for (p = ATp [J] ; p < ATp [J+1] ; p++) {
        i = ATi [p] ;
        q = cs_leaf (i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
        if (jleaf >= 1) delta [j]++ ;   /* A(i,j) is in skeleton */
        if (jleaf == 2) delta [q]-- ;   /* account for overlap in q */
      }
    }
    if (parent [j] != -1) ancestor [j] = parent [j] ;
  }

  /* sum up delta's of each child */
  for (j = 0 ; j < n ; j++) {
    if (parent [j] != -1) colcount [parent [j]] += colcount [j] ;
  }
  cs_spfree(AT);
  cs_free(w);
  return colcount;
} 

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
double cs_cumsum (int *p, int *c, int n) {
  int i, nz = 0 ;
  double nz2 = 0 ;
  if (!p || !c) return (-1) ;     /* check inputs */
  for (i = 0 ; i < n ; i++) {
    p [i] = nz ;
    nz += c [i] ;
    nz2 += c [i] ;              /* also in double to avoid int overflow */
    c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
  }
  p [n] = nz ;
  return (nz2) ;                  /* return sum (c [0..n-1]) */
}

/* depth-first-search of the graph of a matrix, starting at node j */
int cs_dfs (int j, cs *G, int top, int *xi, int *pstack, const int *pinv) {
  int i, p, p2, done, jnew, head = 0, *Gp, *Gi ;
  Gp = G->p ; Gi = G->i ;
  xi [0] = j ;                /* initialize the recursion stack */
  while (head >= 0) {
    j = xi [head] ;         /* get j from the top of the recursion stack */
    jnew = pinv ? (pinv [j]) : j ;
    if (!CS_MARKED (Gp, j)) {
      CS_MARK (Gp, j) ;       /* mark node j as visited */
      pstack [head] = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew]) ;
    }
    done = 1 ;                  /* node j done if no unvisited neighbors */
    p2 = (jnew < 0) ? 0 : CS_UNFLIP (Gp [jnew+1]) ;
    /* examine all neighbors of j */
    for (p = pstack [head] ; p < p2 ; p++) {
      i = Gi [p] ;            /* consider neighbor node i */
      if (CS_MARKED (Gp, i)) continue ;   /* skip visited node i */
      pstack [head] = p ;     /* pause depth-first search of node j */
      xi [++head] = i ;       /* start dfs at node i */
      done = 0 ;              /* node j is not done */
      break ;                 /* break, to start dfs (i) */
    }
    /* depth-first search at node j is done */
    if (done) {
      head-- ;            /* remove j from the recursion stack */
      xi [--top] = j ;    /* and place in the output stack */
    }
  }
  return (top) ;
}

/* breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3) */
static int cs_bfs (const cs *A, int n, int *wi, int *wj, int *queue,
                   const int *imatch, const int *jmatch, int mark) {
  int *Ap, *Ai, head = 0, tail = 0, j, i, p, j2 ;
  cs *C ;
  /* place all unmatched nodes in queue */
  for (j = 0 ; j < n ; j++) {
    if (imatch [j] >= 0) continue ; /* skip j if matched */
    wj [j] = 0 ;                    /* j in set C0 (R0 if transpose) */
    queue [tail++] = j ;            /* place unmatched col j in queue */
  }
  if (tail == 0) return (1) ;         /* quick return if no unmatched nodes */
  if (mark == 1) {
    C = (cs *)A;
  } else {
    C = cs_calloc(1, sizeof (cs));
    cs_transpose(A, C, 0);
  }
  Ap = C->p ; Ai = C->i ;
  /* while queue is not empty */
  while (head < tail) {
    j = queue [head++] ;            /* get the head of the queue */
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      i = Ai [p] ;
      if (wi [i] >= 0) continue ; /* skip if i is marked */
      wi [i] = mark ;             /* i in set R1 (C3 if transpose) */
      j2 = jmatch [i] ;           /* traverse alternating path to j2 */
      if (wj [j2] >= 0) continue ;/* skip j2 if it is marked */
      wj [j2] = mark ;            /* j2 in set C1 (R3 if transpose) */
      queue [tail++] = j2 ;       /* add j2 to queue */
    }
  }
  if (mark != 1) cs_spfree (C) ;      /* free A' if it was created */
  return (1) ;
}

/* collect matched rows and columns into p and q */
static void cs_matched (int n, const int *wj, const int *imatch, int *p, int *q,
                        int *cc, int *rr, int set, int mark) {
  int kc = cc [set], j ;
  int kr = rr [set-1] ;
  for (j = 0 ; j < n ; j++)  {
    if (wj [j] != mark) continue ;      /* skip if j is not in C set */
    p [kr++] = imatch [j] ;
    q [kc++] = j ;
  }
  cc [set+1] = kc ;
  rr [set] = kr ;
}

/* collect unmatched rows into the permutation vector p */
static void cs_unmatched (int m, const int *wi, int *p, int *rr, int set) {
  int i, kr = rr [set] ;
  for (i = 0 ; i < m ; i++) if (wi [i] == 0) p [kr++] = i ;
  rr [set+1] = kr ;
}

/* return 1 if row i is in R2 */
static int cs_rprune (int i, int j, double aij, void *other) {
  int *rr = (int *) other ;
  return (i >= rr [1] && i < rr [2]) ;
}

/* Given A, compute coarse and then fine dmperm */
csd *cs_dmperm (const cs *A, int seed) {
  int m, n, i, j, k, cnz, nc, *jmatch, *imatch, *wi, *wj, *pinv, *Cp, *Ci,
    *ps, *rs, nb1, nb2, *p, *q, *cc, *rr, *r, *s, ok ;
  cs *C ;
  csd *D, *scc ;
  /* --- Maximum matching ------------------------------------------------- */
  m = A->m ; n = A->n ;
  D = cs_calloc (1, sizeof (csd));
  cs_dalloc(D, m, n) ;                      /* allocate result */
  p = D->p ; q = D->q ; r = D->r ; s = D->s ; cc = D->cc ; rr = D->rr ;
  jmatch = cs_maxtrans (A, seed) ;            /* max transversal */
  imatch = jmatch + m ;                       /* imatch = inverse of jmatch */
  /* --- Coarse decomposition --------------------------------------------- */
  wi = r ; wj = s ;                           /* use r and s as workspace */
  for (j = 0 ; j < n ; j++) wj [j] = -1 ;     /* unmark all cols for bfs */
  for (i = 0 ; i < m ; i++) wi [i] = -1 ;     /* unmark all rows for bfs */
  cs_bfs (A, n, wi, wj, q, imatch, jmatch, 1) ;       /* find C1, R1 from C0*/
  ok = cs_bfs (A, m, wj, wi, p, jmatch, imatch, 3) ;  /* find R3, C3 from R0*/
  cs_unmatched (n, wj, q, cc, 0) ;                    /* unmatched set C0 */
  cs_matched (n, wj, imatch, p, q, cc, rr, 1, 1) ;    /* set R1 and C1 */
  cs_matched (n, wj, imatch, p, q, cc, rr, 2, -1) ;   /* set R2 and C2 */
  cs_matched (n, wj, imatch, p, q, cc, rr, 3, 3) ;    /* set R3 and C3 */
  cs_unmatched (m, wi, p, rr, 3) ;                    /* unmatched set R0 */
  cs_free (jmatch) ;
  /* --- Fine decomposition ----------------------------------------------- */
  pinv = cs_pinv (p, m) ;         /* pinv=p' */
  C = cs_permute (A, pinv, q, 0) ;/* C=A(p,q) (it will hold A(R2,C2)) */
  cs_free (pinv) ;
  Cp = C->p ;
  nc = cc [3] - cc [2] ;          /* delete cols C0, C1, and C3 from C */
  if (cc [2] > 0) for (j = cc [2] ; j <= cc [3] ; j++) Cp [j-cc[2]] = Cp [j] ;
  C->n = nc ;
  /* delete rows R0, R1, and R3 from C */
  if (rr [2] - rr [1] < m) {
    cs_fkeep (C, cs_rprune, rr) ;
    cnz = Cp [nc] ;
    Ci = C->i ;
    if (rr [1] > 0) for (k = 0 ; k < cnz ; k++) Ci [k] -= rr [1] ;
  }
  C->m = nc ;
  scc = cs_scc (C) ;              /* find strongly connected components of C*/
  /* --- Combine coarse and fine decompositions --------------------------- */
  ps = scc->p ;                   /* C(ps,ps) is the permuted matrix */
  rs = scc->r ;                   /* kth block is rs[k]..rs[k+1]-1 */
  nb1 = scc->nb  ;                /* # of blocks of A(R2,C2) */
  for (k = 0 ; k < nc ; k++) wj [k] = q [ps [k] + cc [2]] ;
  for (k = 0 ; k < nc ; k++) q [k + cc [2]] = wj [k] ;
  for (k = 0 ; k < nc ; k++) wi [k] = p [ps [k] + rr [1]] ;
  for (k = 0 ; k < nc ; k++) p [k + rr [1]] = wi [k] ;
  nb2 = 0 ;                       /* create the fine block partitions */
  r [0] = s [0] = 0 ;
  if (cc [2] > 0) nb2++ ;         /* leading coarse block A (R1, [C0 C1]) */
  /* coarse block A (R2,C2) */
  for (k = 0 ; k < nb1 ; k++) {
    r [nb2] = rs [k] + rr [1] ; /* A (R2,C2) splits into nb1 fine blocks */
    s [nb2] = rs [k] + cc [2] ;
    nb2++ ;
  }
  if (rr [2] < m) {
    r [nb2] = rr [2] ;          /* trailing coarse block A ([R3 R0], C3) */
    s [nb2] = cc [3] ;
    nb2++ ;
  }
  r [nb2] = m ;
  s [nb2] = n ;
  D->nb = nb2 ;
  cs_dfree (scc) ;
  cs_spfree (C) ;                     /* free temporary matrix */
  return D;
}

static int cs_tol (int i, int j, double aij, void *tol) {
  return (fabs (aij) > *((double *) tol)) ;
}

int cs_droptol (cs *A, double tol) {
  return (cs_fkeep (A, &cs_tol, &tol)) ;    /* keep all large entries */
}

static int cs_nonzero (int i, int j, double aij, void *other) {
  return (aij != 0) ;
}
int cs_dropzeros (cs *A) {
  return (cs_fkeep (A, &cs_nonzero, NULL)) ;  /* keep all nonzero entries */
} 

/* remove duplicate entries from A */
int cs_dupl (cs *A) {
  int i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w ;
  double *Ax ;
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  w = cs_malloc (m, sizeof (int)) ;           /* get workspace */
  for (i = 0 ; i < m ; i++) w [i] = -1 ;      /* row i not yet seen */
  for (j = 0 ; j < n ; j++) {
    q = nz ;                                /* column j will start at q */
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      i = Ai [p] ;                        /* A(i,j) is nonzero */
      if (w [i] >= q) {
        Ax [w [i]] += Ax [p] ;          /* A(i,j) is a duplicate */
      } else {
        w [i] = nz ;                    /* record where row i occurs */
        Ai [nz] = i ;                   /* keep A(i,j) */
        Ax [nz++] = Ax [p] ;
      }
    }
    Ap [j] = q ;                            /* record start of column j */
  }
  Ap [n] = nz ;                               /* finalize A */
  cs_free (w) ;                               /* free workspace */
  cs_sprealloc(A, 0);              /* remove extra space from A */
  return 1;
}

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
int cs_ereach (const cs *A, int k, const int *parent, int *s, int *w) {
  int i, p, n, len, top, *Ap, *Ai ;
  top = n = A->n ; Ap = A->p ; Ai = A->i ;
  CS_MARK (w, k) ;                /* mark node k as visited */
  for (p = Ap [k] ; p < Ap [k+1] ; p++) {
    i = Ai [p] ;                /* A(i,k) is nonzero */
    if (i > k) continue ;       /* only use upper triangular part of A */
    /* traverse up etree*/
    for (len = 0 ; !CS_MARKED (w,i) ; i = parent [i]) {
      s [len++] = i ;         /* L(k,i) is nonzero */
      CS_MARK (w, i) ;        /* mark i as visited */
    }
    while (len > 0) s [--top] = s [--len] ; /* push path onto stack */
  }
  for (p = top ; p < n ; p++) CS_MARK (w, s [p]) ;    /* unmark all nodes */
  CS_MARK (w, k) ;                /* unmark node k */
  return (top) ;                  /* s [top..n-1] contains pattern of L(k,:)*/
}

/* compute the etree of A (using triu(A), or A'A without forming A'A */
int *cs_etree (const cs *A, int ata) {
  int i, k, p, m, n, inext, *Ap, *Ai, *w, *parent, *ancestor, *prev ;
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ;
  parent = cs_malloc (n, sizeof (int)) ;              /* allocate result */
  w = cs_malloc (n + (ata ? m : 0), sizeof (int)) ;   /* get workspace */
  ancestor = w ; prev = w + n ;
  if (ata) for (i = 0 ; i < m ; i++) prev [i] = -1 ;
  for (k = 0 ; k < n ; k++) {
    parent [k] = -1 ;                   /* node k has no parent yet */
    ancestor [k] = -1 ;                 /* nor does k have an ancestor */
    for (p = Ap [k] ; p < Ap [k+1] ; p++) {
      i = ata ? (prev [Ai [p]]) : (Ai [p]) ;
      /* traverse from i to k */
      for ( ; i != -1 && i < k ; i = inext) {
        inext = ancestor [i] ;              /* inext = ancestor of i */
        ancestor [i] = k ;                  /* path compression */
        if (inext == -1) parent [i] = k ;   /* no anc., parent is k */
      }
      if (ata) prev [Ai [p]] = k ;
    }
  }
  cs_free(w);
  return parent;
}

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
int cs_fkeep (cs *A, int (*fkeep) (int, int, double, void *), void *other) {
  int j, p, nz = 0, n, *Ap, *Ai ;
  double *Ax ;
  n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++) {
    p = Ap [j] ;                        /* get current location of col j */
    Ap [j] = nz ;                       /* record new location of col j */
    for ( ; p < Ap [j+1] ; p++) {
      if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other)) {
        if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
        Ai [nz++] = Ai [p] ;
      }
    }
  }
  Ap [n] = nz ;                           /* finalize A */
  cs_sprealloc (A, 0) ;                   /* remove extra space from A */
  return (nz) ;
}

/* y = A*x+y */
int cs_gaxpy (const cs *A, const double *x, double *y) {
  int p, j, n, *Ap, *Ai ;
  double *Ax ;
  n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++) {
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      y [Ai [p]] += Ax [p] * x [j] ;
    }
  }
  return (1) ;
}

/* apply the ith Householder vector to x */
int cs_happly (const cs *V, int i, double beta, double *x) {
  int p, *Vp, *Vi ;
  double *Vx, tau = 0 ;
  Vp = V->p ; Vi = V->i ; Vx = V->x ;
  for (p = Vp [i] ; p < Vp [i+1] ; p++) {
    /* tau = v'*x */
    tau += Vx [p] * x [Vi [p]] ;
  }
  tau *= beta ;                           /* tau = beta*(v'*x) */
  /* x = x - v*tau */
  for (p = Vp [i] ; p < Vp [i+1] ; p++) {
    x [Vi [p]] -= Vx [p] * tau ;
  }
  return (1) ;
}

/* create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
 * where (I-beta*v*v')*x = s*e1.  See Algo 5.1.1, Golub & Van Loan, 3rd ed. */
double cs_house (double *x, double *beta, int n) {
  double s, sigma = 0 ;
  int i ;
  if (!x || !beta) return (-1) ;          /* check inputs */
  for (i = 1 ; i < n ; i++) sigma += x [i] * x [i] ;
  if (sigma == 0) {
    s = fabs (x [0]) ;                  /* s = |x(0)| */
    (*beta) = (x [0] <= 0) ? 2 : 0 ;
    x [0] = 1 ;
  } else {
    s = sqrt (x [0] * x [0] + sigma) ;  /* s = norm (x) */
    x [0] = (x [0] <= 0) ? (x [0] - s) : (-sigma / (x [0] + s)) ;
    (*beta) = -1. / (s * x [0]) ;
  }
  return (s) ;
}

/* x(p) = b, for dense vectors x and b; p=NULL denotes identity */
int cs_ipvec (const int *p, const double *b, double *x, int n) {
  int k ;
  if (!x || !b) return (0) ;                              /* check inputs */
  for (k = 0 ; k < n ; k++) x [p ? p [k] : k] = b [k] ;
  return (1) ;
}

/* consider A(i,j), node j in ith row subtree and return lca(jprev,j) */
int cs_leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
             int *ancestor, int *jleaf) {
  int q, s, sparent, jprev ;
  if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1) ;
  *jleaf = 0 ;
  if (i <= j || first [j] <= maxfirst [i]) return (-1) ;  /* j not a leaf */
  maxfirst [i] = first [j] ;      /* update max first[j] seen so far */
  jprev = prevleaf [i] ;          /* jprev = previous leaf of ith subtree */
  prevleaf [i] = j ;
  *jleaf = (jprev == -1) ? 1: 2 ; /* j is first or subsequent leaf */
  if (*jleaf == 1) return (i) ;   /* if 1st leaf, q = root of ith subtree */
  for (q = jprev ; q != ancestor [q] ; q = ancestor [q]) ;
  for (s = jprev ; s != q ; s = sparent) {
    sparent = ancestor [s] ;    /* path compression */
    ancestor [s] = q ;
  }
  return (q) ;                    /* q = least common ancester (jprev,j) */
}

/* solve Lx=b where x and b are dense.  x=b on input, solution on output. */
int cs_lsolve (const cs *L, double *x) {
  int p, j, n, *Lp, *Li ;
  double *Lx ;
  n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
  for (j = 0 ; j < n ; j++) {
    x [j] /= Lx [Lp [j]] ;
    for (p = Lp [j]+1 ; p < Lp [j+1] ; p++) {
      x [Li [p]] -= Lx [p] * x [j] ;
    }
  }
  return (1) ;
}

/* solve L'x=b where x and b are dense.  x=b on input, solution on output. */
int cs_ltsolve (const cs *L, double *x) {
  int p, j, n, *Lp, *Li ;
  double *Lx ;
  n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
  for (j = n-1 ; j >= 0 ; j--) {
    for (p = Lp [j]+1 ; p < Lp [j+1] ; p++) {
      x [j] -= Lx [p] * x [Li [p]] ;
    }
    x [j] /= Lx [Lp [j]] ;
  }
  return (1) ;
}

/* [L,U,pinv]=lu(A, [q lnz unz]). lnz and unz can be guess */
csn *cs_lu (const cs *A, const css *S, double tol) {
  cs *L, *U ;
  csn *N ;
  double pivot, *Lx, *Ux, *x,  a, t ;
  int *Lp, *Li, *Up, *Ui, *pinv, *xi, *q, n, ipiv, k, top, p, i, col, lnz,unz;
  n = A->n ;
  q = S->q ; lnz = S->lnz ; unz = S->unz ;
  x = cs_malloc (n, sizeof (double)) ;            /* get double workspace */
  xi = cs_malloc (2*n, sizeof (int)) ;            /* get int workspace */
  N = cs_calloc (1, sizeof (csn)) ;               /* allocate result */
  N->L = L = cs_calloc(1, sizeof (cs));
  cs_spalloc (L, n, n, lnz, 1) ;       /* allocate result L */
  N->U = U = cs_calloc(1, sizeof (cs));
  cs_spalloc (U, n, n, unz, 1) ;       /* allocate result U */
  N->pinv = pinv = cs_malloc (n, sizeof (int)) ;  /* allocate result pinv */
  Lp = L->p ; Up = U->p ;
  for (i = 0 ; i < n ; i++) x [i] = 0 ;           /* clear workspace */
  for (i = 0 ; i < n ; i++) pinv [i] = -1 ;       /* no rows pivotal yet */
  for (k = 0 ; k <= n ; k++) Lp [k] = 0 ;         /* no cols of L yet */
  lnz = unz = 0 ;
  /* compute L(:,k) and U(:,k) */
  for (k = 0 ; k < n ; k++) {
    /* --- Triangular solve --------------------------------------------- */
    Lp [k] = lnz ;              /* L(:,k) starts here */
    Up [k] = unz ;              /* U(:,k) starts here */
    if (lnz + n > L->nzmax) {
      cs_sprealloc(L, 2*L->nzmax + n);
    } else if (unz + n > U->nzmax) {
      cs_sprealloc(U, 2*U->nzmax + n);
    }
    Li = L->i ; Lx = L->x ; Ui = U->i ; Ux = U->x ;
    col = q ? (q [k]) : k ;
    top = cs_spsolve (L, A, col, xi, x, pinv, 1) ;  /* x = L\A(:,col) */
    /* --- Find pivot --------------------------------------------------- */
    ipiv = -1 ;
    a = -1 ;
    for (p = top ; p < n ; p++) {
      i = xi [p] ;            /* x(i) is nonzero */
      if (pinv [i] < 0)       /* row i is not yet pivotal */
        {
          if ((t = fabs (x [i])) > a)
            {
              a = t ;         /* largest pivot candidate so far */
              ipiv = i ;
            }
        } else {
        /* x(i) is the entry U(pinv[i],k) */
        Ui [unz] = pinv [i] ;
        Ux [unz++] = x [i] ;
      }
    }
    if (ipiv == -1 || a <= 0) {
      cs_free(xi);
      cs_free(x);
      cs_nfree(N);
      return NULL;
    }
    if (pinv [col] < 0 && fabs (x [col]) >= a*tol) ipiv = col ;
    /* --- Divide by pivot ---------------------------------------------- */
    pivot = x [ipiv] ;          /* the chosen pivot */
    Ui [unz] = k ;              /* last entry in U(:,k) is U(k,k) */
    Ux [unz++] = pivot ;
    pinv [ipiv] = k ;           /* ipiv is the kth pivot row */
    Li [lnz] = ipiv ;           /* first entry in L(:,k) is L(k,k) = 1 */
    Lx [lnz++] = 1 ;
    for (p = top ; p < n ; p++) /* L(k+1:n,k) = x / pivot */
      {
        i = xi [p] ;
        if (pinv [i] < 0)       /* x(i) is an entry in L(:,k) */
          {
            Li [lnz] = i ;      /* save unpermuted row in L */
            Lx [lnz++] = x [i] / pivot ;    /* scale pivot column */
          }
        x [i] = 0 ;             /* x [0..n-1] = 0 for next k */
      }
  }
  /* --- Finalize L and U ------------------------------------------------- */
  Lp [n] = lnz ;
  Up [n] = unz ;
  Li = L->i ;                     /* fix row indices of L for final pinv */
  for (p = 0 ; p < lnz ; p++) Li [p] = pinv [Li [p]] ;
  cs_sprealloc (L, 0) ;           /* remove extra space from L and U */
  cs_sprealloc (U, 0) ;
  cs_free(xi) ;                       /* free workspace */
  cs_free(x) ;
  return N;
}

/* x=A\b where A is unsymmetric; b overwritten with solution */
int cs_lusol (int order, const cs *A, double *b, double tol) {
  double *x ;
  css *S ;
  csn *N ;
  int n, ok ;
  n = A->n ;
  S = cs_sqr (order, A, 0) ;              /* ordering and symbolic analysis */
  N = cs_lu (A, S, tol) ;                 /* numeric LU factorization */
  x = cs_malloc (n, sizeof (double)) ;    /* get workspace */
  ok = (S && N && x) ;
  if (ok) {
    cs_ipvec (N->pinv, b, x, n) ;       /* x = b(p) */
    cs_lsolve (N->L, x) ;               /* x = L\x */
    cs_usolve (N->U, x) ;               /* x = U\x */
    cs_ipvec (S->q, x, b, n) ;          /* b(q) = x */
  }
  cs_free (x) ;
  cs_sfree (S) ;
  cs_nfree (N) ;
  return (ok) ;
}

/* wrapper for malloc */
void *cs_malloc (int n, size_t size) {
  return (malloc (CS_MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *cs_calloc (int n, size_t size) {
  return (calloc (CS_MAX (n,1), size)) ;
}

/* wrapper for free */
void *cs_free (void *p) {
  if (p) free (p) ;       /* free p if it is not already NULL */
  return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void *cs_realloc (void *p, int n, size_t size) {
  return realloc(p, CS_MAX (n,1) * size);
}

/* find an augmenting path starting at column k and extend the match if found */
static void cs_augment (int k, const cs *A, int *jmatch, int *cheap, int *w,
                        int *js, int *is, int *ps) {
  int found = 0, p, i = -1, *Ap = A->p, *Ai = A->i, head = 0, j ;
  js [0] = k ;                        /* start with just node k in jstack */
  while (head >= 0) {
    /* --- Start (or continue) depth-first-search at node j ------------- */
    j = js [head] ;                 /* get j from top of jstack */
    if (w [j] != k)                 /* 1st time j visited for kth path */
      {
        w [j] = k ;                 /* mark j as visited for kth path */
        for (p = cheap [j] ; p < Ap [j+1] && !found ; p++)
          {
            i = Ai [p] ;            /* try a cheap assignment (i,j) */
            found = (jmatch [i] == -1) ;
          }
        cheap [j] = p ;             /* start here next time j is traversed*/
        if (found)
          {
            is [head] = i ;         /* column j matched with row i */
            break ;                 /* end of augmenting path */
          }
        ps [head] = Ap [j] ;        /* no cheap match: start dfs for j */
      }
    /* --- Depth-first-search of neighbors of j ------------------------- */
    for (p = ps [head] ; p < Ap [j+1] ; p++)
      {
        i = Ai [p] ;                /* consider row i */
        if (w [jmatch [i]] == k) continue ; /* skip jmatch [i] if marked */
        ps [head] = p + 1 ;         /* pause dfs of node j */
        is [head] = i ;             /* i will be matched with j if found */
        js [++head] = jmatch [i] ;  /* start dfs at column jmatch [i] */
        break ;
      }
    if (p == Ap [j+1]) head-- ;     /* node j is done; pop from stack */
  }                                   /* augment the match if path found: */
  if (found) for (p = head ; p >= 0 ; p--) jmatch [is [p]] = js [p] ;
}

/* find a maximum transveral */
int *cs_maxtrans (const cs *A, int seed) {
  /*[jmatch [0..m-1]; imatch [0..n-1]]*/
  int i, j, k, n, m, p, n2 = 0, m2 = 0, *Ap, *jimatch, *w, *cheap, *js, *is,
    *ps, *Ai, *Cp, *jmatch, *imatch, *q ;
  cs *C ;
  n = A->n ; m = A->m ; Ap = A->p ; Ai = A->i ;
  w = jimatch = cs_calloc (m+n, sizeof (int)) ;   /* allocate result */
  /* count nonempty rows and columns */
  for (k = 0, j = 0 ; j < n ; j++) {
    n2 += (Ap [j] < Ap [j+1]) ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      w [Ai [p]] = 1 ;
      k += (j == Ai [p]) ;        /* count entries already on diagonal */
    }
  }

  /* quick return if diagonal zero-free */
  if (k == CS_MIN (m,n)) {
    jmatch = jimatch ; imatch = jimatch + m ;
    for (i = 0 ; i < k ; i++) jmatch [i] = i ;
    for (      ; i < m ; i++) jmatch [i] = -1 ;
    for (j = 0 ; j < k ; j++) imatch [j] = j ;
    for (      ; j < n ; j++) imatch [j] = -1 ;
    return jimatch;
  }
  for (i = 0 ; i < m ; i++) m2 += w [i] ;
  if (m2 < n2) {
    C = cs_calloc(1, sizeof (cs));
    cs_transpose (A, C, 0);
  } else {
    C = (cs *)A;
  }
  n = C->n ; m = C->m ; Cp = C->p ;
  jmatch = (m2 < n2) ? jimatch + n : jimatch ;
  imatch = (m2 < n2) ? jimatch : jimatch + m ;
  w = cs_malloc (5*n, sizeof (int)) ;             /* get workspace */
  cheap = w + n ; js = w + 2*n ; is = w + 3*n ; ps = w + 4*n ;
  for (j = 0 ; j < n ; j++) cheap [j] = Cp [j] ;  /* for cheap assignment */
  for (j = 0 ; j < n ; j++) w [j] = -1 ;          /* all columns unflagged */
  for (i = 0 ; i < m ; i++) jmatch [i] = -1 ;     /* nothing matched yet */
  q = cs_randperm (n, seed) ;                     /* q = random permutation */
  /* augment, starting at column q[k] */
  for (k = 0 ; k < n ; k++) {
    cs_augment (q ? q [k]: k, C, jmatch, cheap, w, js, is, ps) ;
  }
  cs_free (q) ;
  for (j = 0 ; j < n ; j++) imatch [j] = -1 ;     /* find row match */
  for (i = 0 ; i < m ; i++) if (jmatch [i] >= 0) imatch [jmatch [i]] = i ;
  if (m2 < n2) cs_spfree(C);
  cs_free(w);
  return jimatch;
}

/* C = A*B */
void cs_multiply (cs *C, const cs *A, const cs *B) {
  int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
  double *x, *Bx, *Cx ;
   ;
  m = A->m ; anz = A->p [A->n] ;
  n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
  w = cs_calloc (m, sizeof (int)) ;                    /* get workspace */
  values = (A->x != NULL) && (Bx != NULL) ;
  x = values ? cs_malloc (m, sizeof (double)) : NULL ; /* get workspace */
  cs_spalloc (C, m, n, anz + bnz, values) ;        /* allocate result */
  Cp = C->p ;
  for (j = 0 ; j < n ; j++) {
    if (nz + m > C->nzmax) {
      cs_sprealloc (C, 2*(C->nzmax)+m);
    }
    Ci = C->i ; Cx = C->x ;         /* C->i and C->x may be reallocated */
    Cp [j] = nz ;                   /* column j of C starts here */
    for (p = Bp [j] ; p < Bp [j+1] ; p++) {
      nz = cs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
    }
    if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
  }
  Cp [n] = nz ;                       /* finalize the last column of C */
  cs_sprealloc (C, 0) ;               /* remove extra space from C */
  cs_free(w);
  cs_free(x);
}

/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
double cs_norm (const cs *A) {
  int p, j, n, *Ap ;
  double *Ax,  norm = 0, s ;
  n = A->n ; Ap = A->p ; Ax = A->x ;
  for (j = 0 ; j < n ; j++) {
    for (s = 0, p = Ap [j] ; p < Ap [j+1] ; p++) s += fabs (Ax [p]) ;
    norm = CS_MAX (norm, s) ;
  }
  return (norm) ;
}

/* C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1. */
cs *cs_permute (const cs *A, const int *pinv, const int *q, int values) {
  int t, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci ;
  double *Cx, *Ax ;
  cs *C ;
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  C = cs_calloc(1, sizeof (cs));
  cs_spalloc (C, m, n, Ap [n], values && Ax != NULL) ;  /* alloc result */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (k = 0 ; k < n ; k++) {
    Cp [k] = nz ;                   /* column k of C is column q[k] of A */
    j = q ? (q [k]) : k ;
    for (t = Ap [j] ; t < Ap [j+1] ; t++) {
      if (Cx) Cx [nz] = Ax [t] ;  /* row i of A is row pinv[i] of C */
      Ci [nz++] = pinv ? (pinv [Ai [t]]) : Ai [t] ;
    }
  }
  Cp [n] = nz ;                       /* finalize the last column of C */
  return C;
}

/* pinv = p', or p = pinv' */
int *cs_pinv (int const *p, int n) {
  int k, *pinv ;
  if (!p) return (NULL) ;                     /* p = NULL denotes identity */
  pinv = cs_malloc (n, sizeof (int)) ;        /* allocate result */
  for (k = 0 ; k < n ; k++) pinv [p [k]] = k ;/* invert the permutation */
  return (pinv) ;                             /* return result */
}

/* post order a forest */
int *cs_post (const int *parent, int n) {
  int j, k = 0, *post, *w, *head, *next, *stack ;
  post = cs_malloc (n, sizeof (int)) ;                /* allocate result */
  w = cs_malloc (3*n, sizeof (int)) ;                 /* get workspace */
  head = w ; next = w + n ; stack = w + 2*n ;
  for (j = 0 ; j < n ; j++) head [j] = -1 ;           /* empty linked lists */
  for (j = n-1 ; j >= 0 ; j--)            /* traverse nodes in reverse order*/
    {
      if (parent [j] == -1) continue ;    /* j is a root */
      next [j] = head [parent [j]] ;      /* add j to list of its parent */
      head [parent [j]] = j ;
    }
  for (j = 0 ; j < n ; j++)
    {
      if (parent [j] != -1) continue ;    /* skip j if it is not a root */
      k = cs_tdfs (j, k, head, next, post, stack) ;
    }
  cs_free(w);
  return post;
}

/* x = b(p), for dense vectors x and b; p=NULL denotes identity */
int cs_pvec (const int *p, const double *b, double *x, int n) {
  int k ;
  if (!x || !b) return (0) ;                              /* check inputs */
  for (k = 0 ; k < n ; k++) x [k] = b [p ? p [k] : k] ;
  return (1) ;
}

/* sparse QR factorization [V,beta,pinv,R] = qr (A) */
csn *cs_qr (const cs *A, const css *S) {
  double *Rx, *Vx, *Ax, *x,  *Beta ;
  int i, k, p, m, n, vnz, p1, top, m2, len, col, rnz, *s, *leftmost, *Ap, *Ai,
    *parent, *Rp, *Ri, *Vp, *Vi, *w, *pinv, *q ;
  cs *R, *V ;
  csn *N ;
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  q = S->q ; parent = S->parent ; pinv = S->pinv ; m2 = S->m2 ;
  vnz = S->lnz ; rnz = S->unz ; leftmost = S->leftmost ;
  w = cs_malloc (m2+n, sizeof (int)) ;            /* get int workspace */
  x = cs_malloc (m2, sizeof (double)) ;           /* get double workspace */
  N = cs_calloc (1, sizeof (csn)) ;               /* allocate result */
  s = w + m2 ;                                    /* s is size n */
  for (k = 0 ; k < m2 ; k++) x [k] = 0 ;          /* clear workspace x */
  N->L = V = cs_calloc(1, sizeof (cs));
  cs_spalloc (V, m2, n, vnz, 1) ;      /* allocate result V */
  N->U = R = cs_calloc(1, sizeof (cs));
  cs_spalloc (R, m2, n, rnz, 1) ;      /* allocate result R */
  N->B = Beta = cs_malloc (n, sizeof (double)) ;  /* allocate result Beta */
  Rp = R->p ; Ri = R->i ; Rx = R->x ;
  Vp = V->p ; Vi = V->i ; Vx = V->x ;
  for (i = 0 ; i < m2 ; i++) w [i] = -1 ; /* clear w, to mark nodes */
  rnz = 0 ; vnz = 0 ;
  /* compute V and R */
  for (k = 0 ; k < n ; k++) {
    Rp [k] = rnz ;                      /* R(:,k) starts here */
    Vp [k] = p1 = vnz ;                 /* V(:,k) starts here */
    w [k] = k ;                         /* add V(k,k) to pattern of V */
    Vi [vnz++] = k ;
    top = n ;
    col = q ? q [k] : k ;
    /* find R(:,k) pattern */
    for (p = Ap [col] ; p < Ap [col+1] ; p++) {
      i = leftmost [Ai [p]] ;         /* i = min(find(A(i,q))) */
      /* traverse up to k */
      for (len = 0 ; w [i] != k ; i = parent [i]) {
        s [len++] = i ;
        w [i] = k ;
      }
      /* push path on stack */
      while (len > 0) s [--top] = s [--len] ; 
      i = pinv [Ai [p]] ;             /* i = permuted row of A(:,col) */
      x [i] = Ax [p] ;                /* x (i) = A(:,col) */
      if (i > k && w [i] < k)         /* pattern of V(:,k) = x (k+1:m) */
        {
          Vi [vnz++] = i ;            /* add i to pattern of V(:,k) */
          w [i] = k ;
        }
    }
    /* for each i in pattern of R(:,k) */
    for (p = top ; p < n ; p++) {
      i = s [p] ;                     /* R(i,k) is nonzero */
      cs_happly (V, i, Beta [i], x) ; /* apply (V(i),Beta(i)) to x */
      Ri [rnz] = i ;                  /* R(i,k) = x(i) */
      Rx [rnz++] = x [i] ;
      x [i] = 0 ;
      if (parent [i] == k) vnz = cs_scatter (V, i, 0, w, NULL, k, V, vnz);
    }

    /* gather V(:,k) = x */
    for (p = p1 ; p < vnz ; p++) {
      Vx [p] = x [Vi [p]] ;
      x [Vi [p]] = 0 ;
    }
    Ri [rnz] = k ;                     /* R(k,k) = norm (x) */
    Rx [rnz++] = cs_house (Vx+p1, Beta+k, vnz-p1) ; /* [v,beta]=house(x) */
  }
  Rp [n] = rnz ;                          /* finalize R */
  Vp [n] = vnz ;                          /* finalize V */
  cs_free(w);
  cs_free(x);
  return N;
}

/* x=A\b where A can be rectangular; b overwritten with solution */
int cs_qrsol (int order, const cs *A, double *b) {
  double *x ;
  css *S ;
  csn *N ;
  cs *AT = NULL ;
  int k, m, n, ok ;
  n = A->n ;
  m = A->m ;
  if (m >= n) {
    S = cs_sqr (order, A, 1) ;          /* ordering and symbolic analysis */
    N = cs_qr (A, S) ;                  /* numeric QR factorization */
    x = cs_calloc (S ? S->m2 : 1, sizeof (double)) ;    /* get workspace */
    ok = (S && N && x) ;
    if (ok)
      {
        cs_ipvec (S->pinv, b, x, m) ;   /* x(0:m-1) = b(p(0:m-1) */
        for (k = 0 ; k < n ; k++)       /* apply Householder refl. to x */
          {
            cs_happly (N->L, k, N->B [k], x) ;
          }
        cs_usolve (N->U, x) ;           /* x = R\x */
        cs_ipvec (S->q, x, b, n) ;      /* b(q(0:n-1)) = x(0:n-1) */
      }
  } else {
    AT = cs_calloc(1, sizeof (cs));
    cs_transpose (A, AT, 1) ;          /* Ax=b is underdetermined */
    S = cs_sqr (order, AT, 1) ;         /* ordering and symbolic analysis */
    N = cs_qr (AT, S) ;                 /* numeric QR factorization of A' */
    x = cs_calloc (S ? S->m2 : 1, sizeof (double)) ;    /* get workspace */
    ok = (AT && S && N && x) ;
    if (ok) {
      cs_pvec (S->q, b, x, m) ;       /* x(q(0:m-1)) = b(0:m-1) */
      cs_utsolve (N->U, x) ;          /* x = R'\x */
      for (k = m-1 ; k >= 0 ; k--)    /* apply Householder refl. to x */
        {
          cs_happly (N->L, k, N->B [k], x) ;
        }
      cs_pvec (S->pinv, x, b, n) ;    /* b(0:n-1) = x(p(0:n-1)) */
    }
  }
  cs_free (x) ;
  cs_sfree (S) ;
  cs_nfree (N) ;
  cs_spfree (AT) ;
  return (ok) ;
}

/* return a random permutation vector, the identity perm, or p = n-1:-1:0.
 * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
 * p = random permutation.  */
int *cs_randperm (int n, int seed) {
  int *p, k, j, t ;
  if (seed == 0) return (NULL) ;      /* return p = NULL (identity) */
  p = cs_malloc (n, sizeof (int)) ;   /* allocate result */
  for (k = 0 ; k < n ; k++) p [k] = n-k-1 ;
  if (seed == -1) return (p) ;        /* return reverse permutation */
  srand (seed) ;                      /* get new random number seed */
  for (k = 0 ; k < n ; k++) {
    j = k + (rand ( ) % (n-k)) ;    /* j = rand integer in range k to n-1 */
    t = p [j] ;                     /* swap p[k] and p[j] */
    p [j] = p [k] ;
    p [k] = t ;
  }
  return (p) ;
}

/* xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
 * xi [n...2n-1] used as workspace */
int cs_reach (cs *G, const cs *B, int k, int *xi, const int *pinv) {
  int p, n, top, *Bp, *Bi, *Gp ;
  n = G->n ; Bp = B->p ; Bi = B->i ; Gp = G->p ;
  top = n ;
  for (p = Bp [k] ; p < Bp [k+1] ; p++) {
    /* start a dfs at unmarked node i */
    if (!CS_MARKED (Gp, Bi [p])) {
      top = cs_dfs (Bi [p], G, top, xi, xi+n, pinv) ;
    }
  }
  for (p = top ; p < n ; p++) CS_MARK (Gp, xi [p]) ;  /* restore G */
  return (top) ;
}

/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
int cs_scatter (const cs *A, int j, double beta, int *w, double *x, int mark,
                cs *C, int nz) {
  int i, p, *Ap, *Ai, *Ci ;
  double *Ax ;
  Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
  for (p = Ap [j] ; p < Ap [j+1] ; p++) {
    i = Ai [p] ;                            /* A(i,j) is nonzero */
    if (w [i] < mark) {
      w [i] = mark ;                      /* i is new entry in column j */
      Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
      if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
    }
    else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
  }
  return (nz) ;
}

/* find the strongly connected components of a square matrix */
csd *cs_scc (cs *A) {
  /* matrix A temporarily modified, then restored */
  int n, i, k, b, nb = 0, top, *xi, *pstack, *p, *r, *Ap, *ATp, *rcopy, *Blk ;
  cs *AT ;
  csd *D ;
  n = A->n ; Ap = A->p ;
  D = cs_calloc (1, sizeof (csd));
  cs_dalloc(D, n, 0) ;                          /* allocate result */
  AT = cs_calloc(1, sizeof (cs));
  cs_transpose (A, AT, 0) ;                      /* AT = A' */
  xi = cs_malloc (2*n+1, sizeof (int)) ;          /* get workspace */
  Blk = xi ; rcopy = pstack = xi + n ;
  p = D->p ; r = D->r ; ATp = AT->p ;
  top = n ;
  /* first dfs(A) to find finish times (xi) */
  for (i = 0 ; i < n ; i++) {
    if (!CS_MARKED (Ap, i)) top = cs_dfs (i, A, top, xi, pstack, NULL) ;
  }
  for (i = 0 ; i < n ; i++) CS_MARK (Ap, i) ; /* restore A; unmark all nodes*/
  top = n ;
  nb = n ;
  /* dfs(A') to find strongly connnected comp */
  for (k = 0 ; k < n ; k++) {
    i = xi [k] ;            /* get i in reverse order of finish times */
    if (CS_MARKED (ATp, i)) continue ;  /* skip node i if already ordered */
    r [nb--] = top ;        /* node i is the start of a component in p */
    top = cs_dfs (i, AT, top, p, pstack, NULL) ;
  }
  r [nb] = 0 ;                /* first block starts at zero; shift r up */
  for (k = nb ; k <= n ; k++) r [k-nb] = r [k] ;
  D->nb = nb = n-nb ;         /* nb = # of strongly connected components */
  /* sort each block in natural order */
  for (b = 0 ; b < nb ; b++) {
    for (k = r [b] ; k < r [b+1] ; k++) Blk [p [k]] = b ;
  }
  for (b = 0 ; b <= nb ; b++) rcopy [b] = r [b] ;
  for (i = 0 ; i < n ; i++) p [rcopy [Blk [i]]++] = i ;
  cs_spfree (AT) ;                     /* free temporary matrix */
  cs_free (xi) ;                       /* free workspace */
  return D;
}

/* ordering and symbolic analysis for a Cholesky factorization */
css *cs_schol (int order, const cs *A) {
  int n, *c, *post, *P ;
  cs *C ;
  css *S ;
  n = A->n ;
  S = cs_calloc (1, sizeof (css)) ;       /* allocate result S */
  P = cs_amd (order, A) ;                 /* P = amd(A+A'), or natural */
  S->pinv = cs_pinv (P, n) ;              /* find inverse permutation */
  cs_free (P) ;
  if (order && !S->pinv) {
    cs_sfree(S);
    return NULL;
  }
  C = cs_symperm (A, S->pinv, 0) ;        /* C = spones(triu(A(P,P))) */
  S->parent = cs_etree (C, 0) ;           /* find etree of C */
  post = cs_post (S->parent, n) ;         /* postorder the etree */
  c = cs_counts (C, S->parent, post, 0) ; /* find column counts of chol(C) */
  cs_free (post) ;
  cs_spfree (C) ;
  S->cp = cs_malloc (n+1, sizeof (int)) ; /* allocate result S->cp */
  S->unz = S->lnz = cs_cumsum (S->cp, c, n) ; /* find column pointers for L */
  cs_free (c) ;
  if (S->lnz >= 0) {
    return S;
  } else {
    cs_sfree(S);
    return NULL;
  }
}

/* solve Gx=b(:,k), where G is either upper (lo=0) or lower (lo=1) triangular */
int cs_spsolve (cs *G, const cs *B, int k, int *xi, double *x, const int *pinv,
                int lo) {
  int j, J, p, q, px, top, n, *Gp, *Gi, *Bp, *Bi ;
  double *Gx, *Bx ;
  Gp = G->p ; Gi = G->i ; Gx = G->x ; n = G->n ;
  Bp = B->p ; Bi = B->i ; Bx = B->x ;
  top = cs_reach (G, B, k, xi, pinv) ;        /* xi[top..n-1]=Reach(B(:,k)) */
  for (p = top ; p < n ; p++) x [xi [p]] = 0 ;    /* clear x */
  for (p = Bp [k] ; p < Bp [k+1] ; p++) x [Bi [p]] = Bx [p] ; /* scatter B */
  for (px = top ; px < n ; px++) {
    j = xi [px] ;                               /* x(j) is nonzero */
    J = pinv ? (pinv [j]) : j ;                 /* j maps to col J of G */
    if (J < 0) continue ;                       /* column J is empty */
    x [j] /= Gx [lo ? (Gp [J]) : (Gp [J+1]-1)] ;/* x(j) /= G(j,j) */
    p = lo ? (Gp [J]+1) : (Gp [J]) ;            /* lo: L(j,j) 1st entry */
    q = lo ? (Gp [J+1]) : (Gp [J+1]-1) ;        /* up: U(j,j) last entry */
    for ( ; p < q ; p++)
      {
        x [Gi [p]] -= Gx [p] * x [j] ;          /* x(i) -= G(i,j) * x(j) */
      }
  }
  return (top) ;                                  /* return top of stack */
}

/* compute nnz(V) = S->lnz, S->pinv, S->leftmost, S->m2 from A and S->parent */
static int cs_vcount (const cs *A, css *S) {
  int i, k, p, pa, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i, *next, *head,
    *tail, *nque, *pinv, *leftmost, *w, *parent = S->parent ;
  S->pinv = pinv = cs_malloc (m+n, sizeof (int)) ;        /* allocate pinv, */
  S->leftmost = leftmost = cs_malloc (m, sizeof (int)) ;  /* and leftmost */
  w = cs_malloc (m+3*n, sizeof (int)) ;   /* get workspace */
  next = w ; head = w + m ; tail = w + m + n ; nque = w + m + 2*n ;
  for (k = 0 ; k < n ; k++) head [k] = -1 ;   /* queue k is empty */
  for (k = 0 ; k < n ; k++) tail [k] = -1 ;
  for (k = 0 ; k < n ; k++) nque [k] = 0 ;
  for (i = 0 ; i < m ; i++) leftmost [i] = -1 ;
  for (k = n-1 ; k >= 0 ; k--) {
    for (p = Ap [k] ; p < Ap [k+1] ; p++) {
      /* leftmost[i] = min(find(A(i,:))) */
      leftmost [Ai [p]] = k;
    }
  }

  /* scan rows in reverse order */
  for (i = m-1 ; i >= 0 ; i--) {
    pinv [i] = -1 ;                     /* row i is not yet ordered */
    k = leftmost [i] ;
    if (k == -1) continue ;             /* row i is empty */
    if (nque [k]++ == 0) tail [k] = i ; /* first row in queue k */
    next [i] = head [k] ;               /* put i at head of queue k */
    head [k] = i ;
  }
  S->lnz = 0 ;
  S->m2 = m ;
  /* find row permutation and nnz(V)*/
  for (k = 0 ; k < n ; k++) {
    i = head [k] ;                      /* remove row i from queue k */
    S->lnz++ ;                          /* count V(k,k) as nonzero */
    if (i < 0) i = S->m2++ ;            /* add a fictitious row */
    pinv [i] = k ;                      /* associate row i with V(:,k) */
    if (--nque [k] <= 0) continue ;     /* skip if V(k+1:m,k) is empty */
    S->lnz += nque [k] ;                /* nque [k] is nnz (V(k+1:m,k)) */
    /* move all rows to parent of k */
    if ((pa = parent [k]) != -1) {
      if (nque [pa] == 0) tail [pa] = tail [k] ;
      next [tail [k]] = head [pa] ;
      head [pa] = next [i] ;
      nque [pa] += nque [k] ;
    }
  }
  for (i = 0 ; i < m ; i++) if (pinv [i] < 0) pinv [i] = k++ ;
  cs_free (w) ;
  return (1) ;
}

/* symbolic ordering and analysis for QR or LU */
css *cs_sqr (int order, const cs *A, int qr) {
  int n, k, ok = 1, *post ;
  css *S ;
  n = A->n ;
  S = cs_calloc (1, sizeof (css)) ;       /* allocate result S */
  S->q = cs_amd (order, A) ;              /* fill-reducing ordering */
  if (order && !S->q) {
    cs_sfree(S);
    return NULL;
  }
  /* QR symbolic analysis */
  if (qr) {
    cs *C = order ? cs_permute (A, NULL, S->q, 0) : ((cs *) A) ;
    S->parent = cs_etree (C, 1) ;       /* etree of C'*C, where C=A(:,q) */
    post = cs_post (S->parent, n) ;
    S->cp = cs_counts (C, S->parent, post, 1) ;  /* col counts chol(C'*C) */
    cs_free (post) ;
    ok = C && S->parent && S->cp && cs_vcount (C, S) ;
    if (ok) for (S->unz = 0, k = 0 ; k < n ; k++) S->unz += S->cp [k] ;
    if (order) cs_spfree (C) ;
  } else {
    S->unz = 4*(A->p [n]) + n ;         /* for LU factorization only, */
    S->lnz = S->unz ;                   /* guess nnz(L) and nnz(U) */
  }
  if (ok) {
    return S;
  } else {
    cs_sfree(S);
    return NULL;
  }
}

/* C = A(p,p) where A and C are symmetric the upper part stored; pinv not p */
cs *cs_symperm (const cs *A, const int *pinv, int values) {
  int i, j, p, q, i2, j2, n, *Ap, *Ai, *Cp, *Ci, *w ;
  double *Cx, *Ax ;
  cs *C ;
  n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  C = cs_calloc(1, sizeof (cs));
  cs_spalloc (C, n, n, Ap [n], values && (Ax != NULL)) ; /* alloc result*/
  w = cs_calloc (n, sizeof (int)) ;                   /* get workspace */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  /* count entries in each column of C */
  for (j = 0 ; j < n ; j++) {
    j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      i = Ai [p] ;
      if (i > j) continue ;       /* skip lower triangular part of A */
      i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
      w [CS_MAX (i2, j2)]++ ;     /* column count of C */
    }
  }
  cs_cumsum (Cp, w, n) ;              /* compute column pointers of C */
  for (j = 0 ; j < n ; j++) {
    j2 = pinv ? pinv [j] : j ;      /* column j of A is column j2 of C */
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      i = Ai [p] ;
      if (i > j) continue ;       /* skip lower triangular part of A*/
      i2 = pinv ? pinv [i] : i ;  /* row i of A is row i2 of C */
      Ci [q = w [CS_MAX (i2, j2)]++] = CS_MIN (i2, j2) ;
      if (Cx) Cx [q] = Ax [p] ;
    }
  }
  cs_free (w);
  return C;
}

/* depth-first search and postorder of a tree rooted at node j */
int cs_tdfs (int j, int k, int *head, const int *next, int *post, int *stack) {
  int i, p, top = 0 ;
  if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
  stack [0] = j ;                 /* place j on the stack */
  /* while (stack is not empty) */
  while (top >= 0) {
    p = stack [top] ;           /* p = top of stack */
    i = head [p] ;              /* i = youngest child of p */
    if (i == -1) {
      top-- ;                 /* p has no unordered children left */
      post [k++] = p ;        /* node p is the kth postordered node */
    } else {
      head [p] = next [i] ;   /* remove i from children of p */
      stack [++top] = i ;     /* start dfs on child node i */
    }
  }
  return (k) ;
}

/* C = A' */
void cs_transpose (const cs *A, cs *C, int values) {
  int p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
  double *Cx, *Ax ;
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  cs_spalloc(C, n, m, Ap [n], values && Ax) ;       /* allocate result */
  w = cs_calloc (m, sizeof (int)) ;                      /* get workspace */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
  cs_cumsum (Cp, w, m) ;                                 /* row pointers */
  for (j = 0 ; j < n ; j++) {
    for (p = Ap [j] ; p < Ap [j+1] ; p++) {
      Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
      if (Cx) Cx [q] = Ax [p] ;
    }
  }
  cs_free(w);
}

/* sparse Cholesky update/downdate, L*L' + sigma*w*w' (sigma = +1 or -1) */
int cs_updown (cs *L, int sigma, const cs *C, const int *parent) {
  int n, p, f, j, *Lp, *Li, *Cp, *Ci ;
  double *Lx, *Cx, alpha, beta = 1, delta, gamma, w1, w2, *w, beta2 = 1 ;
  Lp = L->p ; Li = L->i ; Lx = L->x ; n = L->n ;
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  if ((p = Cp [0]) >= Cp [1]) return (1) ;        /* return if C empty */
  w = cs_malloc (n, sizeof (double)) ;            /* get workspace */
  f = Ci [p] ;
  for ( ; p < Cp [1] ; p++) f = CS_MIN (f, Ci [p]) ;  /* f = min (find (C)) */
  for (j = f ; j != -1 ; j = parent [j]) w [j] = 0 ;  /* clear workspace w */
  for (p = Cp [0] ; p < Cp [1] ; p++) w [Ci [p]] = Cx [p] ; /* w = C */
  /* walk path f up to root */
  for (j = f ; j != -1 ; j = parent [j]) {
    p = Lp [j] ;
    alpha = w [j] / Lx [p] ;                    /* alpha = w(j) / L(j,j) */
    beta2 = beta*beta + sigma*alpha*alpha ;
    if (beta2 <= 0) break ;                     /* not positive definite */
    beta2 = sqrt (beta2) ;
    delta = (sigma > 0) ? (beta / beta2) : (beta2 / beta) ;
    gamma = sigma * alpha / (beta2 * beta) ;
    Lx [p] = delta * Lx [p] + ((sigma > 0) ? (gamma * w [j]) : 0) ;
    beta = beta2 ;
    for (p++ ; p < Lp [j+1] ; p++) {
      w1 = w [Li [p]] ;
      w [Li [p]] = w2 = w1 - alpha * Lx [p] ;
      Lx [p] = delta * Lx [p] + gamma * ((sigma > 0) ? w1 : w2) ;
    }
  }
  cs_free (w) ;
  return (beta2 > 0) ;
}

/* solve Ux=b where x and b are dense.  x=b on input, solution on output. */
int cs_usolve (const cs *U, double *x) {
  int p, j, n, *Up, *Ui ;
  double *Ux ;
  n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
  for (j = n-1 ; j >= 0 ; j--) {
    x [j] /= Ux [Up [j+1]-1] ;
    for (p = Up [j] ; p < Up [j+1]-1 ; p++) {
      x [Ui [p]] -= Ux [p] * x [j] ;
    }
  }
  return (1) ;
}

/* allocate a sparse matrix */
void cs_spalloc(cs *A, int m, int n, int nzmax, int values) {
  A->m = m;
  A->n = n;
  A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
  A->p = cs_malloc(n+1, sizeof (int));
  A->i = cs_malloc(nzmax, sizeof (int)) ;
  A->x = values ? cs_malloc (nzmax, sizeof (double)) : NULL ;
}

/* change the max # of entries sparse matrix */
void cs_sprealloc (cs *A, int nzmax) {
  if (nzmax <= 0) nzmax = A->p[A->n];
  A->i = cs_realloc(A->i, nzmax, sizeof (int));
  if (A->x) A->x = cs_realloc (A->x, nzmax, sizeof (double));
  A->nzmax = nzmax ;
}

/* free a sparse matrix */
void cs_spfree(cs *A) {
  if (!A) return;
  cs_free(A->p);
  cs_free(A->i);
  cs_free(A->x);
  cs_free(A);
}

/* free a numeric factorization */
void cs_nfree (csn *N) {
  if (!N) return;
  cs_spfree(N->L) ;
  cs_spfree(N->U) ;
  cs_free(N->pinv) ;
  cs_free(N->B) ;
  cs_free(N);
}

/* free a symbolic factorization */
void cs_sfree (css *S) {
  if (!S) return;
  cs_free (S->pinv) ;
  cs_free (S->q) ;
  cs_free (S->parent) ;
  cs_free (S->cp) ;
  cs_free (S->leftmost) ;
  cs_free (S);
}

/* allocate a cs_dmperm or cs_scc result */
void cs_dalloc (csd *D, int m, int n) {
  D->p = cs_malloc (m, sizeof (int)) ;
  D->r = cs_malloc (m+6, sizeof (int)) ;
  D->q = cs_malloc (n, sizeof (int)) ;
  D->s = cs_malloc (n+6, sizeof (int)) ;
}

/* free a cs_dmperm or cs_scc result */
void cs_dfree (csd *D) {
  if (!D) return;
  cs_free (D->p) ;
  cs_free (D->q) ;
  cs_free (D->r) ;
  cs_free (D->s) ;
  cs_free (D);
}

/* solve U'x=b where x and b are dense.  x=b on input, solution on output. */
int cs_utsolve (const cs *U, double *x) {
  int p, j, n, *Up, *Ui ;
  double *Ux ;
  n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
  for (j = 0 ; j < n ; j++) {
    for (p = Up [j] ; p < Up [j+1]-1 ; p++) {
      x [j] -= Ux [p] * x [Ui [p]] ;
    }
    x [j] /= Ux [Up [j+1]-1] ;
  }
  return (1) ;
}
