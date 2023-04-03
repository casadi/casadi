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


// SYMBOL "jac_prob"
template<typename T1>
struct casadi_jac_prob {
  // Number of outputs, i.e. rows of the Jacobian
  casadi_int n_out;
  // Number of inputs, i.e. columns of the Jacobian
  casadi_int n_in;
  // Number of colors
  casadi_int n_color;
  // Extended Jacobian sparsity
  const casadi_int* sp_ext;
  // Jacobian coloring
  const casadi_int* coloring;
  // Nominal values for inputs, if any
  const T1* nom_in;
  // Index mapping for outputs (i.e. Jacobian rows), if any
  const size_t* map_out;
  // Index mapping for inputs (i.e.  Jacobian columns), if any
  const size_t* map_in;
};
// C-REPLACE "casadi_jac_prob<T1>" "struct casadi_jac_prob"

// SYMBOL "jac_setup"
template<typename T1>
void casadi_jac_setup(casadi_jac_prob<T1>* p, const casadi_int* sp_ext,
    const casadi_int* coloring) {
  // Set pointers
  p->sp_ext = sp_ext;
  p->coloring = coloring;
  // Dimensions are given by the sparsity patterns
  p->n_out = sp_ext[0];
  p->n_in = sp_ext[1];
  p->n_color = coloring[1];
  // The following defaults to null
  p->nom_in = 0;
  p->map_out = 0;
  p->map_in = 0;
}

// SYMBOL "jac_work"
template<typename T1>
void casadi_jac_work(const casadi_jac_prob<T1>* p, casadi_int* sz_iw, casadi_int* sz_w) {
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  // Work vectors in data struct
  *sz_iw += p->n_in;  // iseed
  *sz_w += p->n_in;  // seed
  *sz_iw += p->n_out;  // isens
  *sz_w += p->n_out;  // sens
  *sz_w += p->n_out;  // scal
  *sz_iw += p->n_out;  // wrt
  *sz_iw += p->n_out;  // nzind
}

// SYMBOL "jac_data"
template<typename T1>
struct casadi_jac_data {
  // Number of seeds, sensitivities for the current color
  casadi_int nseed, nsens;
  // Inputs that are being seeded
  casadi_int *iseed;
  // Set of seeds for the seeded inputs
  T1 *seed;
  // Set of outputs for which sensitivities are calculated
  casadi_int *isens;
  // Set of values for the calculated sensitivities
  T1 *sens;
  // Scaling factors for calculated sensitivities
  T1 *scal;
  // Input corresponding to calculated sensitivities
  casadi_int *wrt;
  // Jacobian nonzero corresponding to calculated sensitivities
  casadi_int *nzind;
};
// C-REPLACE "casadi_jac_data<T1>" "struct casadi_jac_data"

// SYMBOL "jac_init"
template<typename T1>
void casadi_jac_init(const casadi_jac_prob<T1>* p, casadi_jac_data<T1>* d,
    casadi_int** iw, T1** w) {
  // Set work vectors
  d->iseed = *iw; *iw += p->n_in;
  d->seed = *w; *w += p->n_in;
  d->isens = *iw; *iw += p->n_out;
  d->sens = *w; *w += p->n_out;
  d->scal = *w; *w += p->n_out;
  d->wrt = *iw; *iw += p->n_out;
  d->nzind = *iw; *iw += p->n_out;
}

// SYMBOL "jac_pre"
template<typename T1>
void casadi_jac_pre(const casadi_jac_prob<T1>* p, casadi_jac_data<T1>* d, casadi_int c) {
  // Local variables
  casadi_int i, kc, vin, vout, Jk;
  double nom, inv_nom;
  const casadi_int *color_colind, *color_row, *jac_colind, *jac_row;
  // Extract sparsities
  color_colind = p->coloring + 2;
  color_row = color_colind + p->n_color + 1;
  jac_colind = p->sp_ext + 2;
  jac_row = jac_colind + p->n_in + 1;
  // Loop over input indices for color
  d->nseed = d->nsens = 0;
  for (kc = color_colind[c]; kc < color_colind[c + 1]; ++kc) {
    vin = color_row[kc];
    // Nominal value, used as a seed for the column
    nom = p->nom_in ? p->nom_in[vin] : 1;
    inv_nom = 1. / nom;
    // Collect seeds for column
    d->seed[d->nseed] = nom;
    d->iseed[d->nseed] = vin;
    d->nseed++;
    // Request corresponding outputs
    for (Jk = jac_colind[vin]; Jk < jac_colind[vin + 1]; ++Jk) {
      vout = jac_row[Jk];
      d->scal[d->nsens] = inv_nom;
      d->isens[d->nsens] = vout;
      d->wrt[d->nsens] = vin;
      d->nzind[d->nsens] = Jk;
      d->nsens++;
    }
  }
  // Map indices
  if (p->map_in) {
    for (i = 0; i < d->nseed; ++i) d->iseed[i] = p->map_in[d->iseed[i]];
    for (i = 0; i < d->nsens; ++i) d->wrt[i] = p->map_in[d->wrt[i]];
  }
  if (p->map_out) {
    for (i = 0; i < d->nsens; ++i) d->isens[i] = p->map_out[d->isens[i]];
  }
}

// SYMBOL "jac_scale"
template<typename T1>
void casadi_jac_scale(const casadi_jac_prob<T1>* p, casadi_jac_data<T1>* d) {
  // Local variables
  casadi_int i;
  // Scale derivatives
  for (i = 0; i < d->nsens; ++i) d->sens[i] *= d->scal[i];
}

// SYMBOL "get_sub"
template<typename T1>
void casadi_get_sub(T1* sub, const casadi_int* sp_a, const T1* nz_a,
    casadi_int rbegin, casadi_int rend, casadi_int cbegin, casadi_int cend) {
  // Local variables
  casadi_int nc, r, c, k;
  const casadi_int *colind, *row;
  // Quick return if null
  if (sub == 0) return;
  // Extract sparsity
  nc = sp_a[1];
  colind = sp_a + 2;
  row = colind + nc + 1;
  // Loop over columns
  for (c = cbegin; c < cend; ++c) {
    // Loop over nonzeros
    for (k = colind[c]; k < colind[c + 1]; ++k) {
      // Get row
      r = row[k];
      // Save, if in range
      if (r >= rbegin && r < rend) *sub++ = nz_a[k];
    }
  }
}
