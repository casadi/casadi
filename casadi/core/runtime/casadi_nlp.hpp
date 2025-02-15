//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//


// SYMBOL "nlpsol_detect_bounds_prob"
template<typename T1>
struct casadi_nlpsol_detect_bounds_prob {
  casadi_int sz_arg;
  casadi_int sz_res;
  casadi_int sz_iw;
  casadi_int sz_w;
  // Original number og constraints
  casadi_int ng;
  // Number of bounds
  casadi_int nb;
  const casadi_int *target_x;
  const casadi_int *target_g;
  const char *is_simple;

  int (*callback)(const T1** arg, T1** res, casadi_int* iw, T1* w, void* callback_data);
  void* callback_data;
};
// C-REPLACE "casadi_nlpsol_detect_bounds_prob<T1>" "struct casadi_nlpsol_detect_bounds_prob"

// SYMBOL "nlpsol_prob"
template<typename T1>
struct casadi_nlpsol_prob {
  casadi_int nx, ng, np;

  casadi_nlpsol_detect_bounds_prob<T1> detect_bounds;
};
// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"

// SYMBOL "nlpsol_detect_bounds_data"
template<typename T1>
struct casadi_nlpsol_detect_bounds_data {
  // Work vectors
  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1* w;

  // Simple bounds g(x)=A(b)x+b(p);
  // a[i]*x[target_x[i]]+b[i]
  T1* a;
  T1* b;
  casadi_int* target_l;
  casadi_int* target_u;
  T1* lam_xl;
  T1* lam_xu;
};
// C-REPLACE "casadi_nlpsol_detect_bounds_data<T1>" "struct casadi_nlpsol_detect_bounds_data"
// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"

// SYMBOL "nlpsol_data"
template<typename T1>
struct casadi_nlpsol_data {
  // Problem structure
  const casadi_nlpsol_prob<T1>* prob;
  casadi_oracle_data<T1>* oracle;
  // Variable bounds
  T1 *lbz, *ubz;
  // Current primal solution
  T1 *z;
  // Current dual solution
  T1 *lam;
  // Storage for f when null pointer passed
  T1 objective;

  // NLP data, pointers to arg (no allocations needed)
  const T1 *p, *lbx, *ubx, *lbg, *ubg, *x0, *lam_x0, *lam_g0;
  // NLP results, pointers to res (no allocations needed)
  T1 *f, *x, *g, *lam_x, *lam_g, *lam_p;

  casadi_nlpsol_detect_bounds_data<T1> detect_bounds;
};
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"

// SYMBOL "nlpsol_work"
template<typename T1>
void casadi_nlpsol_work(const casadi_nlpsol_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res,
    casadi_int* sz_iw, casadi_int* sz_w) {
  // Reset sz_arg, sz_res
  *sz_arg = *sz_res = 0;
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  *sz_w += p->nx + p->ng; // z
  *sz_w += p->nx + p->ng; // lbz
  *sz_w += p->nx + p->ng; // ubz
  *sz_w += p->nx + p->ng; // lam

  if (p->detect_bounds.ng) {
    *sz_arg += p->detect_bounds.sz_arg;
    *sz_res += p->detect_bounds.sz_res;
    *sz_iw += p->detect_bounds.sz_iw;
    *sz_w += p->detect_bounds.sz_w;

    *sz_w += p->detect_bounds.nb; // a;
    *sz_w += p->detect_bounds.nb; // b;
    *sz_iw += p->nx; // target_l
    *sz_iw += p->nx; // target_u
    *sz_w += p->nx; // lam_xl;
    *sz_w += p->nx; // lam_xu;
  }

}


// SYMBOL "nlpsol_init"
template<typename T1>
void casadi_nlpsol_init(casadi_nlpsol_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  // Local variables
  casadi_int nx, ng;
  const casadi_nlpsol_prob<T1>* p = d->prob;
  nx = p->nx;
  ng = p->ng;
  // Get matrix number of nonzeros
  d->z = *w; *w += nx + ng;
  d->lbz = *w; *w += nx + ng;
  d->ubz = *w; *w += nx + ng;
  d->lam = *w; *w += nx + ng;

  if (p->detect_bounds.ng) {
    d->detect_bounds.arg = *arg; *arg += p->detect_bounds.sz_arg;
    d->detect_bounds.res = *res; *res += p->detect_bounds.sz_res;
    d->detect_bounds.iw = *iw; *iw += p->detect_bounds.sz_iw;
    d->detect_bounds.w = *w; *w += p->detect_bounds.sz_w;

    d->detect_bounds.a = *w; *w += p->detect_bounds.nb;
    d->detect_bounds.b = *w; *w += p->detect_bounds.nb;
    d->detect_bounds.target_l = *iw; *iw += p->nx;
    d->detect_bounds.target_u = *iw; *iw += p->nx;
    d->detect_bounds.lam_xl = *w; *w += nx;
    d->detect_bounds.lam_xu = *w; *w += nx;
  }

}

// SYMBOL "nlpsol_detect_bounds_before"
template<typename T1>
int casadi_detect_bounds_before(casadi_nlpsol_data<T1>* d_nlp) {
  const casadi_nlpsol_prob<T1>* p_nlp = d_nlp->prob;
  casadi_nlpsol_detect_bounds_data<T1>* d_bounds = &d_nlp->detect_bounds;
  const casadi_nlpsol_detect_bounds_prob<T1>* p_bounds = &p_nlp->detect_bounds;

  casadi_int nx = p_nlp->nx;
  d_bounds->arg[0] = d_nlp->z;
  d_bounds->arg[1] = d_nlp->p;
  d_bounds->res[0] = d_bounds->a;
  d_bounds->res[1] = d_bounds->b;
  p_bounds->callback(d_bounds->arg, d_bounds->res,
    d_bounds->iw, d_bounds->w, p_bounds->callback_data);

  for (casadi_int i=0;i<p_bounds->nb;++i) {
    if (d_bounds->a[i]==0) {
      casadi_int k = p_bounds->target_g[i];
      if (d_nlp->lbg[k]>d_bounds->b[i]) return 1;
      if (d_nlp->ubg[k]<d_bounds->b[i]) return 1;
    }
  }

  T1* lbz = d_nlp->lbz+nx;
  T1* ubz = d_nlp->ubz+nx;
  T1* lam = d_nlp->lam+nx;

  for (casadi_int i=0;i<nx;++i) {
    d_bounds->lam_xl[i] = d_nlp->lam_x0 ? (d_nlp->lam_x0[i]<0)*d_nlp->lam_x0[i] : 0.;
    d_bounds->lam_xu[i] = d_nlp->lam_x0 ? (d_nlp->lam_x0[i]>0)*d_nlp->lam_x0[i] : 0.;
  }

  for (casadi_int i=0;i<nx;++i) {
    d_bounds->target_l[i] = i;
    d_bounds->target_u[i] = i;
  }

  // Update lbz/ubz
  casadi_int k=0;
  for (casadi_int i=0;i<p_bounds->ng;++i) {
    if (p_bounds->is_simple[i]) {
      // Update lbz/ubz

      T1 lb = (d_nlp->lbg[i]-d_bounds->b[k])/d_bounds->a[k];
      T1 ub = (d_nlp->ubg[i]-d_bounds->b[k])/d_bounds->a[k];
      if (d_bounds->a[k]<0) {
        T1 tmp = lb;
        lb = ub;
        ub = tmp;
      }
      casadi_int j = p_bounds->target_x[k];
      lb += d_nlp->z[j];
      ub += d_nlp->z[j];

      if (lb==d_nlp->lbz[j]) {
        if (d_nlp->lam_g0) d_bounds->lam_xl[j] += (d_nlp->lam_g0[i]<0)*d_nlp->lam_g0[i];
      } else if (lb>d_nlp->lbz[j]) {
        d_nlp->lbz[j] = lb;
        d_bounds->target_l[j] = nx+i;
        if (d_nlp->lam_g0) d_bounds->lam_xl[j] = (d_nlp->lam_g0[i]<0)*d_nlp->lam_g0[i];
      }

      if (ub==d_nlp->ubz[j]) {
        if (d_nlp->lam_g0) d_bounds->lam_xu[j] += (d_nlp->lam_g0[i]>0)*d_nlp->lam_g0[i];
      } else if (ub<d_nlp->ubz[j]) {
        d_nlp->ubz[j] = ub;
        d_bounds->target_u[j] = nx+i;
        if (d_nlp->lam_g0) d_bounds->lam_xu[j] = (d_nlp->lam_g0[i]>0)*d_nlp->lam_g0[i];
      }

      k++;
    } else {

      // Update lbz/ubz
      *lbz++ = d_nlp->lbg[i];
      *ubz++ = d_nlp->ubg[i];

      if (d_nlp->lam_g0) *lam++ = d_nlp->lam_g0[i];
    }
  }

  for (casadi_int i=0;i<nx;++i) {
    d_nlp->lam[i] = d_bounds->lam_xl[i]+d_bounds->lam_xu[i];
  }
  return 0;
}

// SYMBOL "nlpsol_detect_bounds_after"
template<typename T1>
int casadi_detect_bounds_after(casadi_nlpsol_data<T1>* d_nlp) {
  const casadi_nlpsol_prob<T1>* p_nlp = d_nlp->prob;
  casadi_nlpsol_detect_bounds_data<T1>* d_bounds = &d_nlp->detect_bounds;
  const casadi_nlpsol_detect_bounds_prob<T1>* p_bounds = &p_nlp->detect_bounds;
  casadi_int nx = p_nlp->nx;

  casadi_fill(d_nlp->lam_x, nx, 0.);
  casadi_fill(d_nlp->lam_g, p_bounds->ng, 0.);

  casadi_int k = 0;
  casadi_int k_normal = 0;
  for (casadi_int i=0;i<p_bounds->ng;++i) {
    if (p_bounds->is_simple[i]) {
      casadi_int j = p_bounds->target_x[k];
      if (d_nlp->g) {
        d_nlp->g[i] = d_bounds->a[k]*d_nlp->z[j]-d_bounds->b[k];
        if (d_nlp->x0) d_nlp->g[i] += d_nlp->x0[j];
      }
      k++;
    } else {
      if (d_nlp->g) d_nlp->g[i] = d_nlp->z[nx+k_normal];
      if (d_nlp->lam_g) d_nlp->lam_g[i] = d_nlp->lam[nx+k_normal];
      k_normal++;
    }
  }

  for (casadi_int i=0;i<nx;++i) {
    if (d_bounds->target_l[i]<nx) {
      if (d_nlp->lam_x) d_nlp->lam_x[i] += (d_nlp->lam[i]<0)*d_nlp->lam[i];
    } else {
      if (d_nlp->lam_g)
        d_nlp->lam_g[d_bounds->target_l[i]-nx] += (d_nlp->lam[i]<0)*d_nlp->lam[i];
    }
    if (d_bounds->target_u[i]<nx) {
      if (d_nlp->lam_x) d_nlp->lam_x[i] += (d_nlp->lam[i]>0)*d_nlp->lam[i];
    } else {
      if (d_nlp->lam_g)
        d_nlp->lam_g[d_bounds->target_u[i]-nx] += (d_nlp->lam[i]>0)*d_nlp->lam[i];
    }
  }

  k = 0;
  for (casadi_int i=0;i<p_bounds->ng;++i) {
    if (p_bounds->is_simple[i]) {
      if (d_nlp->lam_g) d_nlp->lam_g[i] /= d_bounds->a[k];
      k++;
    }
  }

  return 0;
}
