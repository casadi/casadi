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



// SYMBOL "convexify_strategy_t"
typedef enum {
  CVX_REGULARIZE,
  CVX_EIGEN_CLIP,
  CVX_EIGEN_REFLECT
} casadi_convexify_strategy_t;

// SYMBOL "convexify_type_in_t"
typedef enum {
  CVX_SYMM, CVX_TRIL, CVX_TRIU
} casadi_convexify_type_in_t;


// SYMBOL "convexify_config"
template<typename T1>
struct casadi_convexify_config {
  casadi_convexify_strategy_t strategy;
  casadi_convexify_type_in_t type_in;
  const casadi_int * Hsp;
  const casadi_int * Hrsp;
  T1 margin;
  // Projection of Hessian sparsity needed? (cache)
  int Hsp_project;
  // Reordering of Hessian needed for scc? (cache)
  int scc_transform;
  /// Block structure of Hessian for certain convexification methods
  const casadi_int *scc_offset;
  const casadi_int *scc_mapping;
  casadi_int scc_offset_size;
  /// For eigen-* convexification strategies: maximum iterations for symmetric Schur decomposition
  // Needs to be "big enough"
  casadi_int max_iter_eig;
  int verbose;
};
// C-REPLACE "casadi_convexify_config<T1>" "struct casadi_convexify_config"


// SYMBOL "convexify_eval"
template<typename T1>
int convexify_eval(const casadi_convexify_config<T1>* c, const T1* Hin, T1* Hout, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
    casadi_int i, j, k, kk, block_size, offset;
    int ret;
    T1 reg, e;
    T1 *H_block, *w_cvx;

    casadi_int Hrsp_nnz = c->Hrsp[2+c->Hrsp[1]];
    casadi_int nnz = c->Hsp[2+c->Hsp[1]];

    if (c->Hsp_project) {
      if (Hin==Hout) {
        casadi_copy(Hin, Hrsp_nnz, w);
        casadi_project(w, c->Hrsp, Hout, c->Hsp, w+Hrsp_nnz);
      } else {
        casadi_project(Hin, c->Hrsp, Hout, c->Hsp, w);
      }
    } else {
      if (Hin!=Hout) casadi_copy(Hin, nnz, Hout);
    }

    if (c->strategy==CVX_REGULARIZE) {
      // Determine regularization parameter with Gershgorin theorem
      reg = c->margin-casadi_lb_eig(c->Hsp, Hout);
      if (reg > 0) casadi_regularize(c->Hsp, Hout, reg);
    } else if (c->strategy==CVX_EIGEN_REFLECT || c->strategy==CVX_EIGEN_CLIP) {
      offset = 0;

      // Loop over Hessian blocks
      for (k=0;k<c->scc_offset_size-1;++k) {
        block_size = c->scc_offset[k+1]-c->scc_offset[k];

        H_block = w;
        w_cvx = w;

        // Set w_cvx to dense Hessian block from Hout
        if (c->scc_transform) {
          kk=0;
          if (c->type_in==CVX_SYMM) {
            // Loop over columns of block
            for (i=0;i<block_size;++i) {
              // Loop over elements in column
              for (j=0;j<block_size;++j, ++kk) {
                H_block[kk] = Hout[c->scc_mapping[offset+kk]];
              }
            }
          } else if (c->type_in==CVX_TRIU) {
            // Loop over columns of block
            for (i=0;i<block_size;++i) {
              // Loop over elements in column
              for (j=0;j<i+1;++j, ++kk) {
                  e = Hout[c->scc_mapping[offset+kk]];
                  H_block[i*block_size+j] = e;
                  H_block[i+block_size*j] = e;
              }
            }
          } else {
            // Loop over columns of block
            for (i=0;i<block_size;++i) {
              for (j=i;j<block_size;++j, ++kk) {
                e = Hout[c->scc_mapping[offset+kk]];
                H_block[i*block_size+j] = e;
                H_block[i+block_size*j] = e;
              }
            }
          }
          w_cvx += block_size*block_size;
        } else {
          H_block = Hout+offset;
        }

        // Perform convexification
        ret = casadi_cvx(block_size, H_block, c->margin, 1e-10,
          c->strategy==CVX_EIGEN_REFLECT, c->max_iter_eig, w_cvx, iw);
        if (ret) return ret;

        // Fill in upper-rectangular part
        for (i=0;i<block_size;++i) {
          for (j=0;j<i+1;++j) {
            H_block[block_size*i+j] = H_block[block_size*j+i];
          }
        }

        // Put results back in Hout
        if (c->scc_transform) {
          kk=0;
          if (c->type_in==CVX_SYMM) {
            // Loop over columns of block
            for (i=0;i<block_size;++i) {
              // Loop over elements in column
              for (j=0;j<block_size;++j, ++kk) {
                Hout[c->scc_mapping[offset+kk]] = H_block[kk];
              }
            }
          } else if (c->type_in==CVX_TRIU) {
            // Loop over columns of block
            for (i=0;i<block_size;++i) {
              // Loop over elements in column
              for (j=0;j<i+1;++j, ++kk) {
                Hout[c->scc_mapping[offset+kk]] = H_block[block_size*i+j];
              }
            }
          } else {
            // Loop over columns of block
            for (i=0;i<block_size;++i) {
              // Loop over elements in column
              for (j=i;j<block_size;++j, ++kk) {
                Hout[c->scc_mapping[offset+kk]] = H_block[block_size*i+j];
              }
            }
          }
        }

        if (c->type_in==CVX_SYMM) {
          offset += block_size*block_size;
        } else {
          offset += block_size*(block_size+1)/2;
        }
      }
    }

    return 0;
  }
