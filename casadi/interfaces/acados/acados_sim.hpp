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

// acados sim glue shared by the VM and generated code. The acados objects are laid out in
// casadi work memory every call; acados evaluates the model through casadi_acados_ext_eval,
// which dispatches to CasADi functions via calc_function (OracleCallback / oracle callback).
// Requires acados_c/sim_interface.h, acados/sim/sim_{erk,irk}_integrator.h, blasfeo_d_aux.h

// C-REPLACE "casadi_acados_sim_prob<T1>" "struct casadi_acados_sim_prob"
// C-REPLACE "casadi_acados_sim_data<T1>" "struct casadi_acados_sim_data"
// C-REPLACE "casadi_acados_ext_fun<T1>" "struct casadi_acados_ext_fun"
// C-REPLACE "casadi_oracle_data<T1>" "struct casadi_oracle_data"
// C-REPLACE "OracleCallback" "struct casadi_oracle_callback"
// C-REPLACE "calc_function" "casadi_oracle_call"
// C-REPLACE "static_cast<casadi_acados_ext_fun<T1>*>" "(struct casadi_acados_ext_fun*) "
// C-REPLACE "static_cast<casadi_acados_sim_data<T1>*>" "(struct casadi_acados_sim_data*) "
// C-REPLACE "static_cast<sim_collocation_type>" "(sim_collocation_type) "
// C-REPLACE "reinterpret_cast<char*>" "(char*) "
// C-REPLACE "reinterpret_cast<sim_out*>" "(sim_out*) "
// C-REPLACE "reinterpret_cast<sim_info*>" "(sim_info*) "
// C-REPLACE "reinterpret_cast<T1*>" "(casadi_real*) "
// C-REPLACE "reinterpret_cast<size_t>" "(size_t) "
// C-REPLACE "static_cast<size_t>" "(size_t) "
// C-REPLACE "const_cast<T1*>" "(casadi_real*) "
// C-REPLACE "static_cast<const T1*>" "(const casadi_real*) "
// C-REPLACE "static_cast<T1*>" "(casadi_real*) "
// C-REPLACE "static_cast<struct colmaj_args*>" "(struct colmaj_args*) "
// C-REPLACE "static_cast<struct blasfeo_dmat*>" "(struct blasfeo_dmat*) "
// C-REPLACE "static_cast<struct blasfeo_dmat_args*>" "(struct blasfeo_dmat_args*) "
// C-REPLACE "static_cast<struct blasfeo_dvec*>" "(struct blasfeo_dvec*) "
// C-REPLACE "static_cast<struct blasfeo_dvec_args*>" "(struct blasfeo_dvec_args*) "
// C-REPLACE "casadi_acados_ext_eval<T1>" "casadi_acados_ext_eval"
// C-REPLACE "casadi_acados_ext_ws<T1>" "casadi_acados_ext_ws"
// C-REPLACE "casadi_acados_ext_set_ws<T1>" "casadi_acados_ext_set_ws"
// C-REPLACE "casadi_acados_take<T1>" "casadi_acados_take"
// C-REPLACE "casadi_acados_sim_nraw<T1>" "casadi_acados_sim_nraw"

template<typename T1>
struct casadi_acados_sim_prob {
  // Scheme: irk (else erk), stages, steps, Newton iterations, sim_collocation_type
  int irk, num_stages, num_steps, newton_iter, collocation_type;
  int nx, nu, nz;
  // Sensitivities acados memory is sized for
  int sens_forw, sens_adj, sens_hess;
  // Bytes of the acados objects, measured once on the host
  casadi_int sim_bytes;
  // Model functions (the hessian, k = 3, only with sens_hess):
  // callback, acados model field, n_in, n_out, (nrow, ncol) of inputs/outputs
  OracleCallback cb[4];
  const char* field[4];
  const casadi_int* dims[4];
  // Setup-derived: number of model functions, dense conversion buffer length of each
  casadi_int nfun, sz_buf[4];
};

template<typename T1> struct casadi_acados_sim_data;

// acados external function for model function k
template<typename T1>
struct casadi_acados_ext_fun {
  external_function_generic base;  // must be first: acados casts self
  casadi_acados_sim_data<T1>* d;
  int k;
};

template<typename T1>
struct casadi_acados_sim_data {
  const casadi_acados_sim_prob<T1>* prob;
  // Model function evaluation: oracle work, parameters, failure flag
  casadi_oracle_data<T1>* oracle;
  const T1* p;
  // Guess for the algebraic states
  const T1* z0;
  int fun_flag;
  // Dense conversion buffers, one per model function
  T1* buf[4];
  casadi_acados_ext_fun<T1> fun[4];
  // Memory for the acados objects; if measure, only the bytes needed are computed
  char* raw;
  casadi_int nraw;
  int measure;
  sim_config* config;
  void* dims;
  void* opts;
  sim_in* in;
  sim_out* out;
  sim_solver* solver;
};

// SYMBOL "acados_sim_setup"
template<typename T1>
void casadi_acados_sim_setup(casadi_acados_sim_prob<T1>* p) {
  casadi_int k, i, n_in, n_out;
  p->nfun = p->sens_hess ? 4 : 3;
  for (k = 0; k < 4; ++k) p->sz_buf[k] = 0;
  for (k = 0; k < p->nfun; ++k) {
    n_in = p->dims[k][0];
    n_out = p->dims[k][1];
    for (i = 0; i < n_in + n_out; ++i) p->sz_buf[k] += p->dims[k][2 + 2*i] * p->dims[k][3 + 2*i];
  }
}

// SYMBOL "acados_sim_nraw"
// Work vector length of the acados objects (with alignment slack)
template<typename T1>
casadi_int casadi_acados_sim_nraw(const casadi_acados_sim_prob<T1>* p) {
  return (p->sim_bytes + 63) / 8 + 8;
}

// SYMBOL "acados_sim_work"
template<typename T1>
void casadi_acados_sim_work(const casadi_acados_sim_prob<T1>* p, casadi_int* sz_w) {
  casadi_int k;
  for (k = 0; k < 4; ++k) *sz_w += p->sz_buf[k];
  *sz_w += casadi_acados_sim_nraw<T1>(p);
}

// SYMBOL "acados_sim_set_work"
template<typename T1>
void casadi_acados_sim_set_work(casadi_acados_sim_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  const casadi_acados_sim_prob<T1>* p = d->prob;
  casadi_int k;
  for (k = 0; k < 4; ++k) {
    d->buf[k] = *w; *w += p->sz_buf[k];
  }
  d->measure = 0;
  d->nraw = 8 * casadi_acados_sim_nraw<T1>(p);
  d->raw = reinterpret_cast<char*>(*w); *w += casadi_acados_sim_nraw<T1>(p);
}

// Convert an acados external function argument to a dense column-major nr-by-nc buffer
// SYMBOL "acados_unpack"
template<typename T1>
void casadi_acados_unpack(ext_fun_arg_t t, void* in, int nr, int nc, T1* b) {
  int c;
  switch (t) {
    case COLMAJ:
      casadi_copy(static_cast<const T1*>(in), nr * nc, b);
      break;
    case COLMAJ_ARGS:
    {
      struct colmaj_args* a = static_cast<struct colmaj_args*>(in);
      for (c = 0; c < nc; ++c) casadi_copy(a->A + c * a->lda, nr, b + c * nr);
      break;
    }
    case BLASFEO_DMAT:
      blasfeo_unpack_dmat(nr, nc, static_cast<struct blasfeo_dmat*>(in), 0, 0, b, nr);
      break;
    case BLASFEO_DMAT_ARGS:
    {
      struct blasfeo_dmat_args* a = static_cast<struct blasfeo_dmat_args*>(in);
      blasfeo_unpack_dmat(nr, nc, a->A, a->ai, a->aj, b, nr);
      break;
    }
    case BLASFEO_DVEC:
      blasfeo_unpack_dvec(nr, static_cast<struct blasfeo_dvec*>(in), 0, b, 1);
      break;
    case BLASFEO_DVEC_ARGS:
    {
      struct blasfeo_dvec_args* a = static_cast<struct blasfeo_dvec_args*>(in);
      blasfeo_unpack_dvec(nr, a->x, a->xi, b, 1);
      break;
    }
    default:
      break;
  }
}

// Convert a dense column-major nr-by-nc buffer to an acados external function argument
// SYMBOL "acados_pack"
template<typename T1>
void casadi_acados_pack(ext_fun_arg_t t, const T1* b, int nr, int nc, void* out) {
  int c;
  switch (t) {
    case COLMAJ:
      casadi_copy(b, nr * nc, static_cast<T1*>(out));
      break;
    case COLMAJ_ARGS:
    {
      struct colmaj_args* a = static_cast<struct colmaj_args*>(out);
      for (c = 0; c < nc; ++c) casadi_copy(b + c * nr, nr, a->A + c * a->lda);
      break;
    }
    case BLASFEO_DMAT:
      blasfeo_pack_dmat(nr, nc, const_cast<T1*>(b), nr,
        static_cast<struct blasfeo_dmat*>(out), 0, 0);
      break;
    case BLASFEO_DMAT_ARGS:
    {
      struct blasfeo_dmat_args* a = static_cast<struct blasfeo_dmat_args*>(out);
      blasfeo_pack_dmat(nr, nc, const_cast<T1*>(b), nr, a->A, a->ai, a->aj);
      break;
    }
    case BLASFEO_DVEC:
      blasfeo_pack_dvec(nr, const_cast<T1*>(b), 1, static_cast<struct blasfeo_dvec*>(out), 0);
      break;
    case BLASFEO_DVEC_ARGS:
    {
      struct blasfeo_dvec_args* a = static_cast<struct blasfeo_dvec_args*>(out);
      blasfeo_pack_dvec(nr, const_cast<T1*>(b), 1, a->x, a->xi);
      break;
    }
    default:
      break;
  }
}

// Copy a block of a column-major matrix with leading dimension ld, optionally transposed
// SYMBOL "acados_block"
template<typename T1>
void casadi_acados_block(const T1* src, casadi_int ld, casadi_int r0, casadi_int c0,
    casadi_int nr, casadi_int nc, int tr, T1* dst) {
  casadi_int c, k;
  for (c = 0; c < nc; ++c) {
    for (k = 0; k < nr; ++k) {
      dst[k + c * nr] = tr ? src[c0 + c + (r0 + k) * ld] : src[r0 + k + (c0 + c) * ld];
    }
  }
}

// SYMBOL "acados_ext_eval"
// Evaluate model function k for acados: convert arguments, calc_function, convert results
template<typename T1>
void casadi_acados_ext_eval(void* self, ext_fun_arg_t* type_in, void** in,
    ext_fun_arg_t* type_out, void** out) {
  casadi_acados_ext_fun<T1>* e = static_cast<casadi_acados_ext_fun<T1>*>(self);
  casadi_acados_sim_data<T1>* d = e->d;
  const casadi_int* dims = d->prob->dims[e->k];
  casadi_int n_in = dims[0], n_out = dims[1], i, nr, nc;
  casadi_oracle_data<T1>* o = d->oracle;
  T1* b = d->buf[e->k];
  // Inputs, CasADi p last
  for (i = 0; i < n_in - 1; ++i) {
    nr = dims[2 + 2*i];
    nc = dims[3 + 2*i];
    if (type_in[i] == IGNORE_ARGUMENT) {
      o->arg[i] = 0;
    } else {
      casadi_acados_unpack(type_in[i], in[i], nr, nc, b);
      o->arg[i] = b;
    }
    b += nr * nc;
  }
  b += dims[2 + 2*(n_in - 1)] * dims[3 + 2*(n_in - 1)];
  o->arg[n_in - 1] = d->p;
  for (i = 0; i < n_out; ++i) {
    o->res[i] = type_out[i] == IGNORE_ARGUMENT ? 0 : b;
    b += dims[2 + 2*(n_in + i)] * dims[3 + 2*(n_in + i)];
  }
  if (calc_function(&d->prob->cb[e->k], o)) d->fun_flag = 1;
  for (i = 0; i < n_out; ++i) {
    if (o->res[i]) casadi_acados_pack(type_out[i], o->res[i],
      dims[2 + 2*(n_in + i)], dims[3 + 2*(n_in + i)], out[i]);
  }
}

// SYMBOL "acados_ext_ws"
template<typename T1>
size_t casadi_acados_ext_ws(void* self) { return 0; }

// SYMBOL "acados_ext_set_ws"
template<typename T1>
void casadi_acados_ext_set_ws(void* self, void* work) {}

// SYMBOL "acados_take"
// Chunk of n bytes (rounded up to 64) from [*c, end), null if it does not fit; *c 64-aligned.
// Zeroed if clear: acados relies on zeroed memory (calloc)
template<typename T1>
char* casadi_acados_take(char** c, char* end, size_t n, int clear) {
  char* r = *c;
  n = (n + 63) / 64 * 64;
  if (static_cast<size_t>(end - r) < n) return 0;
  *c = r + n;
  if (clear) memset(r, 0, n);
  return r;
}

// SYMBOL "acados_sim_init"
// Create the acados objects in d->raw; returns the bytes needed (independent of the
// alignment of d->raw, including worst-case alignment slack), -1 if d->nraw is too small.
// If d->measure, only config, dims, opts and sim_in are created, sim_out and the solver sized
template<typename T1>
casadi_int casadi_acados_sim_init(casadi_acados_sim_data<T1>* d) {
  const casadi_acados_sim_prob<T1>* p = d->prob;
  char *base = reinterpret_cast<char*>((reinterpret_cast<size_t>(d->raw) + 63) / 64 * 64);
  char *c = base, *end = d->raw + d->nraw, *ptr;
  int nx = p->nx, nz = p->nz, nf = p->nx + p->nu;
  T1* x;
  size_t o, n_out, n_solver;
  if (base > end) return -1;
  int k, tmp;
  bool b;
  double T = 1;
  sim_collocation_type ct = static_cast<sim_collocation_type>(p->collocation_type);
  d->fun_flag = 0;
  // Same steps as sim_config_create, sim_dims_create, ..., sim_solver_create, without calloc
  if (!(ptr = casadi_acados_take<T1>(&c, end, sim_config_calculate_size(), 1))) return -1;
  d->config = sim_config_assign(ptr);
  if (p->irk) {
    sim_irk_config_initialize_default(d->config);
  } else {
    sim_erk_config_initialize_default(d->config);
  }
  if (!(ptr = casadi_acados_take<T1>(&c, end, d->config->dims_calculate_size(), 1))) return -1;
  d->dims = d->config->dims_assign(d->config, ptr);
  sim_dims_set(d->config, d->dims, "nx", &p->nx);
  sim_dims_set(d->config, d->dims, "nu", &p->nu);
  sim_dims_set(d->config, d->dims, "nz", &p->nz);
  ptr = casadi_acados_take<T1>(&c, end, d->config->opts_calculate_size(d->config, d->dims), 1);
  if (!ptr) return -1;
  d->opts = d->config->opts_assign(d->config, d->dims, ptr);
  d->config->opts_initialize_default(d->config, d->dims, d->opts);
  tmp = p->newton_iter;
  sim_opts_set(d->config, d->opts, "newton_iter", &tmp);
  sim_opts_set(d->config, d->opts, "collocation_type", &ct);
  tmp = p->num_stages;
  sim_opts_set(d->config, d->opts, "num_stages", &tmp);
  tmp = p->num_steps;
  sim_opts_set(d->config, d->opts, "num_steps", &tmp);
  // Memory is sized for the sensitivities enabled at creation
  b = p->sens_forw;
  sim_opts_set(d->config, d->opts, "sens_forw", &b);
  b = p->sens_adj;
  sim_opts_set(d->config, d->opts, "sens_adj", &b);
  b = p->sens_hess;
  sim_opts_set(d->config, d->opts, "sens_hess", &b);
  b = false;
  sim_opts_set(d->config, d->opts, "sens_algebraic", &b);
  sim_opts_set(d->config, d->opts, "output_z", &b);
  ptr = casadi_acados_take<T1>(&c, end, sim_in_calculate_size(d->config, d->dims), 1);
  if (!ptr) return -1;
  d->in = sim_in_assign(d->config, d->dims, ptr);
  // sim_out as in sim_out_assign, but S_hess, (nx+nu)^2, only with sens_hess
  o = (sizeof(sim_out) + sizeof(sim_info) + 7) / 8 * 8;
  n_out = o + sizeof(T1) * (nx + nx * nf + nf + (p->sens_hess ? nf * nf : 0) + nf + nz
    + nz * nf);
  if (d->measure) {
    c += (n_out + 63) / 64 * 64;
  } else {
    if (!(ptr = casadi_acados_take<T1>(&c, end, n_out, 1))) return -1;
    d->out = reinterpret_cast<sim_out*>(ptr);
    d->out->info = reinterpret_cast<sim_info*>(ptr + sizeof(sim_out));
    x = reinterpret_cast<T1*>(ptr + o);
    d->out->xn = x; x += nx;
    d->out->S_forw = x; x += nx * nf;
    d->out->S_adj = x; x += nf;
    d->out->S_hess = 0;
    if (p->sens_hess) {
      d->out->S_hess = x; x += nf * nf;
    }
    d->out->grad = x; x += nf;
    d->out->zn = x; x += nz;
    d->out->S_algebraic = x;
  }
  sim_in_set(d->config, d->dims, d->in, "T", &T);
  for (k = 0; k < p->nfun; ++k) {
    d->fun[k].base.evaluate = casadi_acados_ext_eval<T1>;
    d->fun[k].base.get_external_workspace_requirement = casadi_acados_ext_ws<T1>;
    d->fun[k].base.set_external_workspace = casadi_acados_ext_set_ws<T1>;
    d->fun[k].d = d;
    d->fun[k].k = k;
    if (d->config->model_set(d->in->model, p->field[k], &d->fun[k])) return -1;
  }
  d->config->opts_update(d->config, d->dims, d->opts);
  n_solver = sim_calculate_size(d->config, d->dims, d->opts, d->in);
  if (d->measure) return (c - base) + (n_solver + 63) / 64 * 64 + 63;
  if (!(ptr = casadi_acados_take<T1>(&c, end, n_solver, 1))) return -1;
  d->solver = sim_assign(d->config, d->dims, d->opts, d->in, ptr);
  // Seed for forward sensitivities: identity (memory is zeroed)
  for (k = 0; k < p->nx; ++k) d->in->S_forw[k + k * p->nx] = 1;
  if (sim_precompute(d->solver, d->in, d->out)) return -1;
  // Guess for the algebraic states, kept by acados from one solve to the next
  if (nz && d->z0) sim_solver_set(d->solver, "z", const_cast<T1*>(d->z0));
  return (c - base) + 63;
}

// SYMBOL "acados_sim_solve"
// One solve over an interval of length T; signature of casadi_acados_chain_data::solve
template<typename T1>
int casadi_acados_sim_solve(void* ctx, T1 T, const T1* x0, const T1* p, const T1* u, T1* xf,
    T1* S_forw, const T1* seed_adj, T1* S_adj, T1* S_hess) {
  casadi_acados_sim_data<T1>* d = static_cast<casadi_acados_sim_data<T1>*>(ctx);
  bool sens_hess = S_hess != 0;
  bool sens_adj = S_adj != 0 || sens_hess;
  bool sens_forw = S_forw != 0 || sens_hess;
  int flag;
  sim_opts_set(d->config, d->opts, "sens_forw", &sens_forw);
  sim_opts_set(d->config, d->opts, "sens_adj", &sens_adj);
  sim_opts_set(d->config, d->opts, "sens_hess", &sens_hess);
  sim_in_set(d->config, d->dims, d->in, "T", &T);
  sim_in_set(d->config, d->dims, d->in, "x", const_cast<T1*>(x0));
  sim_in_set(d->config, d->dims, d->in, "u", const_cast<T1*>(u));
  if (sens_adj) sim_in_set(d->config, d->dims, d->in, "seed_adj", const_cast<T1*>(seed_adj));
  d->p = p;
  flag = sim_solve(d->solver, d->in, d->out);
  // Failure of acados or of a model function evaluation
  if (flag || d->fun_flag) return 1;
  if (xf) sim_out_get(d->config, d->dims, d->out, "xn", xf);
  if (S_forw) sim_out_get(d->config, d->dims, d->out, "S_forw", S_forw);
  if (S_adj) sim_out_get(d->config, d->dims, d->out, "S_adj", S_adj);
  if (S_hess) sim_out_get(d->config, d->dims, d->out, "S_hess", S_hess);
  return 0;
}
