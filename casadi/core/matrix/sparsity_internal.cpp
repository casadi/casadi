/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "sparsity_internal.hpp"
#include "sparsity_tools.hpp"
#include "../std_vector_tools.hpp"
#include <climits>
#include <cstdlib>
#include <cmath>
#include "matrix.hpp"

using namespace std;

namespace casadi {

  int SparsityInternal::size() const {
    return row_.size();
  }

  int SparsityInternal::numel() const {
    return nrow_*ncol_;
  }

  void SparsityInternal::repr(ostream &stream) const {
    stream << "Compressed Column Storage: ";
    printCompact(stream);
  }

  void SparsityInternal::printCompact(std::ostream &stream) const {
    // Print dimensions
    stream << nrow_ << "x" << ncol_;

    // Print shape
    if (isEmpty()) {
      // Print nothing, shape clear anyway
    } else if (isDense()) {
      stream << ", dense";
    } else if (isDiagonal()) {
      stream << ", diagonal";
    } else {
      stream << ", " << row_.size() << " nnz";
    }
  }

  void SparsityInternal::sanityCheck(bool complete) const {
    casadi_assert_message(nrow_ >=0,
                          "SparsityInternal: number of rows must be positive, but got "
                          << nrow_ << ".");
    casadi_assert_message(ncol_>=0 ,
                          "SparsityInternal: number of columns must be positive, but got "
                          << ncol_ << ".");
    if (colind_.size() != ncol_+1) {
      std::stringstream s;
      s << "SparsityInternal:Compressed Column Storage is not sane. The following must hold:"
        << std::endl;
      s << "  colind.size() = ncol + 1, but got   colind.size() = "
        << colind_.size() << "   and   ncol = "  << ncol_ << std::endl;
      s << "  Note that the signature is as follows: Sparsity (nrow, ncol, colind, row)."
        << std::endl;
      casadi_error(s.str());
    }
    if (complete) {

      if (colind_.size()>0) {
        for (int k=1;k<colind_.size();k++) {
          casadi_assert_message(colind_[k]>=colind_[k-1],
                                "SparsityInternal:Compressed Column Storage is not sane. "
                                "colind must be monotone. Note that the signature is as "
                                "follows: Sparsity (nrow, ncol, colind, row).");
        }

        casadi_assert_message(colind_[0]==0,
                              "SparsityInternal:Compressed Column Storage is not sane. "
                              "First element of colind must be zero. Note that the signature "
                              "is as follows: Sparsity (nrow, ncol, colind, row).");
        if (colind_[(colind_.size()-1)]!=row_.size()) {
          std::stringstream s;
          s << "SparsityInternal:Compressed Column Storage is not sane. The following must hold:"
            << std::endl;
          s << "  colind[lastElement] = row.size(), but got   colind[lastElement] = "
            << colind_[(colind_.size()-1)] << "   and   row.size() = "  << row_.size() << std::endl;
          s << "  Note that the signature is as follows: Sparsity (nrow, ncol, colind, row)."
            << std::endl;
          casadi_error(s.str());
        }
        if (row_.size()>ncol_*nrow_) {
          std::stringstream s;
          s << "SparsityInternal:Compressed Column Storage is not sane. The following must hold:"
            << std::endl;
          s << "  row.size() <= nrow * ncol, but got   row.size()  = " << row_.size()
            << "   and   ncol * nrow = "  << nrow_*ncol_ << std::endl;
          s << "  Note that the signature is as follows: Sparsity (nrow, ncol, colind, row)."
            << std::endl;
          casadi_error(s.str());
        }
      }
      for (int k=0;k<row_.size();k++) {
        if (row_[k]>=nrow_ || row_[k] < 0) {
          std::stringstream s;
          s << "SparsityInternal:Compressed Column Storage is not sane. The following must hold:"
            << std::endl;
          s << "  0 <= row[i] < nrow for each i, but got   row[i] = " << row_[k]
            << "   and   nrow = "  << nrow_ << std::endl;
          s << "  Note that the signature is as follows: Sparsity (nrow, ncol, colind, row)."
            << std::endl;
          casadi_error(s.str());
        }
      }

    }
  }


  void SparsityInternal::print(ostream &stream) const {
    repr(stream);
    stream << endl;
    stream << "colind: " << colind_ << endl;
    stream << "row:    " << row_ << endl;
  }

  vector<int> SparsityInternal::getCol() const {
    vector<int> col(size());
    for (int r=0; r<ncol_; ++r) {
      for (int el = colind_[r]; el < colind_[r+1]; ++el) {
        col[el] = r;
      }
    }
    return col;
  }

  Sparsity SparsityInternal::transpose() const {
    // Dummy mapping
    vector<int> mapping;

    return transpose(mapping);
  }

  Sparsity SparsityInternal::transpose(vector<int>& mapping, bool invert_mapping) const {
    // Get the sparsity of the transpose in sparse triplet form
    const vector<int>& trans_col = row_;
    vector<int> trans_row = getCol();

    // Create the sparsity pattern
    return Sparsity::triplet(ncol_, nrow_, trans_row, trans_col, mapping, invert_mapping);
  }

  std::vector<int> SparsityInternal::eliminationTree(bool ata) const {
    // Allocate result
    vector<int> parent(ncol_);

    // Allocate workspace
    vector<int> ancestor(ncol_);
    vector<int> prev(ata ? nrow_ : 0, -1);

    // Loop over cols
    for (int k=0; k<ncol_; ++k) {
      // Start with no parent or ancestor
      parent[k] = -1;
      ancestor[k] = -1;

      // Loop over nonzeros
      for (int p=colind_[k]; p<colind_[k+1]; ++p) {

        // What is this?
        int i=ata ? (prev[row_[p]]) : (row_[p]);

        // Transverse from i to k
        while (i!=-1 && i<k) {

          // Next i is the ancestor of i
          int inext = ancestor[i];

          // Path compression
          ancestor[i] = k;

          // No ancestor, parent is k
          if (inext==-1)
            parent[i] = k;

          // Update i
          i=inext;
        }

        // What is this?
        if (ata) {
          prev[row_[p]] = k;
        }
      }
    }

    return parent;

  }

  int SparsityInternal::depthFirstSearch(int j, int top, std::vector<int>& xi,
                                         std::vector<int>& pstack, const std::vector<int>& pinv,
                                         std::vector<bool>& marked) const {
    int head = 0;

    // initialize the recursion stack
    xi[0] = j;
    while (head >= 0) {

      // get j from the top of the recursion stack
      j = xi[head];
      int jnew = !pinv.empty() ? (pinv[j]) : j;
      if (!marked[j]) {

        // mark node j as visited
        marked[j]=true;
        pstack[head] = (jnew < 0) ? 0 : colind_[jnew];
      }

      // node j done if no unvisited neighbors
      int done = 1;
      int p2 = (jnew < 0) ? 0 : colind_[jnew+1];

      // examine all neighbors of j
      for (int p = pstack[head]; p< p2; ++p) {

        // consider neighbor node i
        int i = row_[p];

        // skip visited node i
        if (marked[i]) continue ;

        // pause depth-first search of node j
        pstack[head] = p;

        // start dfs at node i
        xi[++head] = i;

        // node j is not done
        done = 0;

        // break, to start dfs (i)
        break;
      }

      //depth-first search at node j is done
      if (done) {
        // remove j from the recursion stack
        head--;

        // and place in the output stack
        xi[--top] = j ;
      }
    }
    return (top) ;
  }

  int SparsityInternal::stronglyConnectedComponents(std::vector<int>& p,
                                                    std::vector<int>& r) const {
    // NOTE: This implementation has been copied from CSparse and then modified,
    // it needs cleaning up to be proper C++
    vector<int> tmp;

    Sparsity AT = transpose();

    vector<int> xi(2*ncol_+1);
    vector<int>& Blk = xi;

    vector<int> pstack(ncol_+1);

    p.resize(ncol_);
    r.resize(ncol_+6);

    vector<bool> marked(ncol_, false);

    int top = ncol_;

    //first dfs(A) to find finish times (xi)
    for (int i = 0; i<ncol_; ++i) {
      if (!marked[i]) {
        top = depthFirstSearch(i, top, xi, pstack, tmp, marked);
      }
    }

    //restore A; unmark all nodes
    fill(marked.begin(), marked.end(), false);

    top = ncol_;
    int nb = ncol_;

    // dfs(A') to find strongly connnected comp
    for (int k=0 ; k < ncol_ ; ++k) {
      // get i in reverse order of finish times
      int i = xi[k];

      // skip node i if already ordered
      if (marked[i]) continue;

      // node i is the start of a component in p
      r[nb--] = top;
      top = AT.depthFirstSearch(i, top, p, pstack, tmp, marked);
    }

    // first block starts at zero; shift r up
    r[nb] = 0;
    for (int k = nb ; k <= ncol_ ; ++k)
      r[k-nb] = r[k] ;

    // nb = # of strongly connected components
    nb = ncol_-nb;

    // sort each block in natural order
    for (int b = 0 ; b < nb ; b++) {
      for (int k = r[b]; k<r[b+1] ; ++k)
        Blk[p[k]] = b ;
    }

    // Get p; shift r down (side effect)
    for (int i=0; i<ncol_; ++i) {
      p[r[Blk[i]]++] = i;
    }

    // Shift up r
    r.resize(nb+1);
    for (int i=nb; i>0; --i) {
      r[i]=r[i-1];
    }
    r[0]=0;

    return nb;
  }

  void SparsityInternal::breadthFirstSearch(int n, std::vector<int>& wi, std::vector<int>& wj,
                                            std::vector<int>& queue, const std::vector<int>& imatch,
                                            const std::vector<int>& jmatch, int mark) const {
    // NOTE: This implementation has been copied from CSparse and then modified,
    // it needs cleaning up to be proper C++
    int head = 0, tail = 0, j, i, p, j2 ;

    // place all unmatched nodes in queue
    for (j=0; j<n; ++j) {
      // skip j if matched
      if (imatch[j] >= 0) continue;

      // j in set C0 (R0 if transpose)
      wj[j] = 0;

      // place unmatched row j in queue
      queue[tail++] = j;
    }

    // quick return if no unmatched nodes
    if (tail == 0) return;

    Sparsity trans;
    const SparsityInternal *C;
    if (mark == 1) {
      C = this;
    } else {
      trans = transpose();
      C = static_cast<const SparsityInternal *>(trans.get());
    }

    // while queue is not empty
    while (head < tail) {

      // get the head of the queue
      j = queue[head++];
      for (p = C->colind_[j] ; p < C->colind_[j+1] ; p++) {
        i = C->row_[p] ;

        // skip if i is marked
        if (wi[i] >= 0) continue;

        // i in set R1 (C3 if transpose)
        wi[i] = mark;

        // traverse alternating path to j2
        j2 = jmatch[i];

        // skip j2 if it is marked
        if (wj[j2] >= 0) continue;

        // j2 in set C1 (R3 if transpose)
        wj[j2] = mark;

        // add j2 to queue
        queue[tail++] = j2;
      }
    }
  }

  void SparsityInternal::matched(int n, const std::vector<int>& wj,
                                 const std::vector<int>& imatch, std::vector<int>& p,
                                 std::vector<int>& q, std::vector<int>& cc, std::vector<int>& rr,
                                 int set, int mark) {
    // NOTE: This implementation has been copied from CSparse and then modified,
    // it needs cleaning up to be proper C++
    int kc = cc[set];
    int kr = rr[set-1] ;
    for (int j=0; j<n; ++j) {
      // skip if j is not in C set
      if (wj[j] != mark) continue;

      p[kr++] = imatch[j] ;
      q[kc++] = j ;
    }

    cc[set+1] = kc ;
    rr[set] = kr ;
  }

  void SparsityInternal::unmatched(int m, const std::vector<int>& wi, std::vector<int>& p,
                                   std::vector<int>& rr, int set) {
    // NOTE: This implementation has been copied from CSparse and then modified,
    // it needs cleaning up to be proper C++
    int i, kr = rr[set] ;
    for (i=0; i<m; i++)
      if (wi[i] == 0)
        p[kr++] = i;

    rr[set+1] = kr;
  }

  int SparsityInternal::rprune(int i, int j, double aij, void *other) {
    // NOTE: This implementation has been copied from CSparse and then modified,
    // it needs cleaning up to be proper C++
    vector<int> &rr = *static_cast<vector<int> *>(other);
    return (i >= rr[1] && i < rr[2]) ;
  }

  void SparsityInternal::augmentingPath(int k, std::vector<int>& jmatch, int *cheap,
                                        std::vector<int>& w, int *js, int *is, int *ps) const {
    // NOTE: This implementation has been copied from CSparse and then modified,
    // it needs cleaning up to be proper C++

    int found = 0, p, i = -1, head = 0, j ;

    // start with just node k in jstack
    js[0] = k ;

    while (head >= 0) {
      // --- Start (or continue) depth-first-search at node j -------------

      // get j from top of jstack
      j = js[head];

      // 1st time j visited for kth path
      if (w[j] != k) {

        // mark j as visited for kth path
        w[j] = k;
        for (p = cheap[j] ; p < colind_[j+1] && !found; ++p) {
          i = row_[p] ;            /* try a cheap assignment (i, j) */
          found = (jmatch[i] == -1) ;
        }

        // start here next time j is traversed
        cheap[j] = p;
        if (found) {
          // row j matched with col i
          is[head] = i;

          // end of augmenting path
          break;
        }

        // no cheap match: start dfs for j
        ps[head] = colind_[j];
      }

      // --- Depth-first-search of neighbors of j -------------------------
      for (p = ps[head]; p<colind_[j+1]; ++p) {

        // consider col i
        i = row_[p];

        // skip jmatch[i] if marked
        if (w[jmatch[i]] == k) continue;

        // pause dfs of node j
        ps[head] = p + 1;

        // i will be matched with j if found
        is[head] = i;

        // start dfs at row jmatch[i]
        js[++head] = jmatch[i];
        break ;
      }

      // node j is done; pop from stack
      if (p == colind_[j+1]) head--;
    } // augment the match if path found:

    if (found)
      for (p = head; p>=0; --p)
        jmatch[is[p]] = js[p];
  }

  void SparsityInternal::maxTransversal(std::vector<int>& imatch, std::vector<int>& jmatch,
                                        Sparsity& trans, int seed) const {
    // NOTE: This implementation has been copied from CSparse and then modified,
    // it needs cleaning up to be proper C++

    int n2 = 0, m2 = 0;

    // allocate result
    jmatch.resize(nrow_);
    imatch.resize(ncol_);
    vector<int> w(nrow_+ncol_);

    // count nonempty cols and rows
    int k=0;
    for (int j=0; j<ncol_; ++j) {
      n2 += (colind_[j] < colind_[j+1]);
      for (int p=colind_[j]; p < colind_[j+1]; ++p) {
        w[row_[p]] = 1;

        // count entries already on diagonal
        k += (j == row_[p]);
      }
    }

    // quick return if diagonal zero-free
    if (k == std::min(nrow_, ncol_)) {
      int i;
      for (i=0; i<k; ++i) jmatch[i] = i;
      for (;    i<nrow_; ++i) jmatch[i] = -1;

      int j;
      for (j=0; j<k; ++j) imatch[j] = j;
      for (;    j<ncol_; ++j) imatch[j] = -1;
    }

    for (int i=0; i<nrow_; ++i) m2 += w[i];

    // transpose if needed
    if (m2 < n2 && trans.isNull())
      trans = transpose();

    // Get pointer to sparsity
    const SparsityInternal* C = m2 < n2 ? static_cast<const SparsityInternal*>(trans.get()) : this;

    std::vector<int>& Cjmatch = m2 < n2 ? imatch : jmatch;
    std::vector<int>& Cimatch = m2 < n2 ? jmatch : imatch;

    // get workspace
    w.resize(5 * C->ncol_);

    int *cheap = &w.front() + C->ncol_;
    int *js = &w.front() + 2*C->ncol_;
    int *is = &w.front() + 3*C->ncol_;
    int *ps = &w.front() + 4*C->ncol_;

    // for cheap assignment
    for (int j=0; j<C->ncol_; ++j)
      cheap[j] = C->colind_[j];

    // all rows unflagged
    for (int j=0; j<C->ncol_; ++j)
      w[j] = -1;

    // nothing matched yet
    for (int i=0; i<C->nrow_; ++i)
      Cjmatch[i] = -1;

    // q = random permutation
    std::vector<int> q = randomPermutation(C->ncol_, seed);

    // augment, starting at row q[k]
    for (k=0; k<C->ncol_; ++k) {
      C->augmentingPath(!q.empty() ? q[k]: k, Cjmatch, cheap, w, js, is, ps);
    }

    // find col match
    for (int j=0; j<C->ncol_; ++j)
      Cimatch[j] = -1;

    for (int i = 0; i<C->nrow_; ++i)
      if (Cjmatch[i] >= 0)
        Cimatch[Cjmatch[i]] = i;
  }

  int SparsityInternal::dulmageMendelsohnUpper(std::vector<int>& rowperm,
                                               std::vector<int>& colperm,
                                               std::vector<int>& rowblock,
                                               std::vector<int>& colblock,
                                               std::vector<int>& coarse_rowblock,
                                               std::vector<int>& coarse_colblock, int seed) const {
    // The transpose of the expression
    Sparsity trans;

    // Part 1: Maximum matching

    // col permutation
    rowperm.resize(nrow_);

    // row permutation
    colperm.resize(ncol_);

    // size nb+1, block k is cols r[k] to r[k+1]-1 in A(p, q)
    rowblock.resize(nrow_+6);

    // size nb+1, block k is rows s[k] to s[k+1]-1 in A(p, q)
    colblock.resize(ncol_+6);

    // coarse col decomposition
    coarse_rowblock.resize(5);
    fill(coarse_rowblock.begin(), coarse_rowblock.end(), 0);

    // coarse row decomposition
    coarse_colblock.resize(5);
    fill(coarse_colblock.begin(), coarse_colblock.end(), 0);

    // max transversal
    vector<int> imatch, jmatch;
    maxTransversal(imatch, jmatch, trans, seed);

    // Coarse decomposition

    // use rowblock and colblock as workspace
    vector<int>& wi = rowblock;
    vector<int>& wj = colblock;

    // unmark all rows for bfs
    for (int j=0; j<ncol_; ++j)
      wj[j] = -1;

    // unmark all cols for bfs
    for (int i=0; i<nrow_; ++i)
      wi[i] = -1 ;

    // find C1, R1 from C0
    breadthFirstSearch(ncol_, wi, wj, colperm, imatch, jmatch, 1);

    // find R3, C3 from R0
    breadthFirstSearch(nrow_, wj, wi, rowperm, jmatch, imatch, 3);

    // unmatched set C0
    unmatched(ncol_, wj, colperm, coarse_colblock, 0);

    // set R1 and C1
    matched(ncol_, wj, imatch, rowperm, colperm, coarse_colblock, coarse_rowblock, 1, 1);

    // set R2 and C2
    matched(ncol_, wj, imatch, rowperm, colperm, coarse_colblock, coarse_rowblock, 2, -1);

    // set R3 and C3
    matched(ncol_, wj, imatch, rowperm, colperm, coarse_colblock, coarse_rowblock, 3, 3);

    // unmatched set R0
    unmatched(nrow_, wi, rowperm, coarse_rowblock, 3);

    // --- Fine decomposition -----------------------------------------------
    // pinv=p'
    vector<int> pinv = invertPermutation(rowperm);

    // C=A(p, q) (it will hold A(R2, C2))
    Sparsity C = permute(pinv, colperm, 0);

    vector<int>& colind_C = C.colindRef();

    // delete rows C0, C1, and C3 from C
    int nc = coarse_colblock[3] - coarse_colblock[2];
    if (coarse_colblock[2] > 0) {
      for (int j = coarse_colblock[2]; j <= coarse_colblock[3]; ++j)
        colind_C[j-coarse_colblock[2]] = colind_C[j];
    }
    C->ncol_ = nc;

    C->colind_.resize(nc+1);
    // delete cols R0, R1, and R3 from C
    if (coarse_rowblock[2] - coarse_rowblock[1] < nrow_) {
      C->drop(rprune, &coarse_rowblock);
      int cnz = colind_C[nc];
      vector<int>& row_C = C->row_;
      if (coarse_rowblock[1] > 0)
        for (int k=0; k<cnz; ++k)
          row_C[k] -= coarse_rowblock[1];
    }
    C->row_.resize(C->colind_.back());
    C->nrow_ = nc ;

    // find strongly connected components of C
    vector<int> scc_p, scc_r;
    int scc_nb = C->stronglyConnectedComponents(scc_p, scc_r);

    // --- Combine coarse and fine decompositions ---------------------------

    // C(ps, ps) is the permuted matrix
    vector<int> ps = scc_p;

    // kth block is rs[k]..rs[k+1]-1
    vector<int> rs = scc_r;

    // # of blocks of A(R2, C2)
    int nb1 = scc_nb;

    for (int k=0; k<nc; ++k)
      wj[k] = colperm[ps[k] + coarse_colblock[2]];

    for (int k=0; k<nc; ++k)
      colperm[k + coarse_colblock[2]] = wj[k];

    for (int k=0; k<nc; ++k)
      wi[k] = rowperm[ps[k] + coarse_rowblock[1]];

    for (int k=0; k<nc; ++k)
      rowperm[k + coarse_rowblock[1]] = wi[k];

    // create the fine block partitions
    int nb2 = 0;
    rowblock[0] = colblock[0] = 0;

    // leading coarse block A (R1, [C0 C1])
    if (coarse_colblock[2] > 0)
      nb2++ ;

    // coarse block A (R2, C2)
    for (int k=0; k<nb1; ++k) {
      // A (R2, C2) splits into nb1 fine blocks
      rowblock[nb2] = rs[k] + coarse_rowblock[1];
      colblock[nb2] = rs[k] + coarse_colblock[2] ;
      nb2++ ;
    }

    if (coarse_rowblock[2] < nrow_) {
      // trailing coarse block A ([R3 R0], C3)
      rowblock[nb2] = coarse_rowblock[2];
      colblock[nb2] = coarse_colblock[3];
      nb2++ ;
    }

    rowblock[nb2] = nrow_;
    colblock[nb2] = ncol_ ;

    // Shrink rowblock and colblock
    rowblock.resize(nb2+1);
    colblock.resize(nb2+1);
    return nb2;
  }

  std::vector<int> SparsityInternal::randomPermutation(int n, int seed) {
    // Return object
    std::vector<int> p;

    // return p = empty (identity)
    if (seed==0) return p;

    // allocate result
    p.resize(n);

    for (int k=0; k<n; ++k)
      p[k] = n-k-1;

    // return reverse permutation
    if (seed==-1) return p;
    #if defined(_WIN32)
    srand(seed);
    #else
    unsigned int seedu = seed;
    #endif

    for (int k=0; k<n; ++k) {
      // j = rand int in range k to n-1
      #if defined(_WIN32)
      int j = k + (rand() % (n-k)); // NOLINT(runtime/threadsafe_fn)
      #else
      int j = k + (rand_r(&seedu) % (n-k));
      #endif
      // swap p[k] and p[j]
      int t = p[j];
      p[j] = p[k];
      p[k] = t;
    }

    return p;
  }

  std::vector<int> SparsityInternal::invertPermutation(const std::vector<int>& p) {
    // pinv = p', or p = pinv'

    // allocate result
    vector<int> pinv(p.size());

    // invert the permutation
    for (int k=0; k<p.size(); ++k)
      pinv[p[k]] = k;

    // return result
    return pinv;
  }

  Sparsity SparsityInternal::permute(const std::vector<int>& pinv,
                                     const std::vector<int>& q, int values) const {
    // alloc result
    Sparsity C = Sparsity::sparse(nrow_, ncol_);

    // Col offset
    vector<int>& colind_C = C.colindRef();

    // Row for each nonzero
    vector<int>& row_C = C.rowRef();
    row_C.resize(size());

    int nz = 0;
    for (int k = 0; k<ncol_; ++k) {
      // row k of C is row q[k] of A
      colind_C[k] = nz;

      int j = !q.empty() ? (q[k]) : k;

      for (int t = colind_[j]; t<colind_[j+1]; ++t) {
        row_C[nz++] = !pinv.empty() ? (pinv[row_[t]]) : row_[t] ;
      }
    }

    // finalize the last row of C
    colind_C[ncol_] = nz;
    return C;
  }

  int SparsityInternal::drop(int (*fkeep)(int, int, double, void *), void *other) {
    int nz = 0;

    for (int j = 0; j<ncol_; ++j) {
      // get current location of row j
      int p = colind_[j];

      // record new location of row j
      colind_[j] = nz;
      for ( ; p < colind_[j+1] ; ++p) {
        if (fkeep(row_[p], j, 1, other)) {
          // keep A(i, j)
          row_[nz++] = row_[p] ;
        }
      }
    }

    // finalize A
    colind_[ncol_] = nz;
    return nz ;
  }

  int SparsityInternal::leaf(int i, int j, const int *first, int *maxfirst, int *prevleaf,
                             int *ancestor, int *jleaf) {
    int q, s, sparent, jprev ;
    if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1) ;
    *jleaf = 0 ;
    if (i <= j || first[j] <= maxfirst[i]) return (-1) ;  /* j not a leaf */
    maxfirst[i] = first[j] ;      /* update max first[j] seen so far */
    jprev = prevleaf[i] ;          /* jprev = previous leaf of ith subtree */
    prevleaf[i] = j ;
    *jleaf = (jprev == -1) ? 1: 2 ; /* j is first or subsequent leaf */
    if (*jleaf == 1) return (i) ;   /* if 1st leaf, q = root of ith subtree */
    for (q = jprev ; q != ancestor[q] ; q = ancestor[q]) ;
    for (s = jprev ; s != q ; s = sparent) {
        sparent = ancestor[s] ;    /* path compression */
        ancestor[s] = q ;
    }
    return (q) ;                    /* q = least common ancester (jprev, j) */
  }

  int SparsityInternal::vcount(std::vector<int>& pinv,
                               std::vector<int>& parent,
                               std::vector<int>& leftmost,
                               int& S_m2, double& S_lnz) const {
    int i, k, p, pa;
    int n = ncol_, m = nrow_;
    const int* Ap = &colind_.front();
    const int* Ai = &row_.front();

    // allocate pinv
    pinv.resize(m+n);
    fill(pinv.begin(), pinv.end(), 0);

    // and leftmost
    leftmost.resize(m);

    // get workspace
    vector<int> w(m+3*n);

    int *next = &w.front();
    int *head = &w.front() + m;
    int *tail = &w.front() + m + n;
    int *nque = &w.front() + m + 2*n;

    // queue k is empty
    for (k = 0 ; k < n ; k++)
      head[k] = -1;

    for (k=0; k<n; ++k)
      tail[k] = -1;

    for (k=0; k<n; ++k)
      nque[k] = 0;

    for (i=0; i<m; ++i)
      leftmost[i] = -1;

    for (k=n-1; k>=0; --k) {
      for (p=Ap[k]; p<Ap[k+1]; ++p) {
        // leftmost[i] = min(find(A(i, :)))
        leftmost[Ai[p]] = k;
      }
    }

    // scan cols in reverse order
    for (i = m-1; i >= 0; i--) {
      // col i is not yet ordered
      pinv[i] = -1;
      k = leftmost[i] ;

      // col i is empty
      if (k == -1) continue;

      // first col in queue k
      if (nque[k]++ == 0)
        tail[k] = i;

      // put i at head of queue k
      next[i] = head[k];
      head[k] = i;
    }
    S_lnz = 0;
    S_m2 = m;

    // find col permutation and nnz(V)
    for (k=0; k<n; ++k) {
      // remove col i from queue k
      i = head[k];

      // count V(k, k) as nonzero
      S_lnz++;

      // add a fictitious col
      if (i < 0)
        i = S_m2++;

      // associate col i with V(:, k)
      pinv[i] = k;

      // skip if V(k+1:m, k) is empty
      if (--nque[k] <= 0) continue;

      // nque [k] is nnz (V(k+1:m, k))
      S_lnz += nque[k];

      // move all cols to parent of k
      if ((pa = parent[k]) != -1) {
        if (nque[pa] == 0)
          tail[pa] = tail[k];

        next[tail[k]] = head[pa] ;
        head[pa] = next[i] ;
        nque[pa] += nque[k] ;
      }
    }
    for (i=0; i<m ; ++i)
      if (pinv[i] < 0)
        pinv[i] = k++;

    pinv.resize(m);
    return 1;
  }

  std::vector<int> SparsityInternal::postorder(const std::vector<int>& parent, int n) {
    int j, k = 0, *head, *next, *stack ;

    // allocate result
    vector<int> post(n);

    // get workspace
    vector<int> w(3*n);

    head = &w.front() ;
    next = &w.front() + n ;
    stack = &w.front() + 2*n;

    // empty linked lists
    for (j=0; j<n; ++j)
      head[j] = -1;

    // traverse nodes in reverse order
    for (j=n-1; j>=0; --j) {
      // j is a root
      if (parent[j] == -1) continue;

      // add j to list of its parent
      next[j] = head[parent[j]];
      head[parent[j]] = j ;
    }

    for (j=0; j<n; ++j) {
      // skip j if it is not a root
      if (parent[j] != -1) continue;

      k = depthFirstSearchAndPostorder(j, k, head, next, &post.front(), stack);
    }

    // success; return post
    return post;
  }

  int SparsityInternal::depthFirstSearchAndPostorder(int j, int k, int *head,
                                                     const int *next, int *post, int *stack) {
    int i, p, top = 0;

    // place j on the stack
    stack[0] = j;

    // while (stack is not empty)
    while (top >= 0) {
      // p = top of stack
      p = stack[top];

      // i = youngest child of p
      i = head[p];
      if (i == -1) {
        // p has no unordered children left
        top--;

        // node p is the kth postordered node
        post[k++] = p;
      } else {
        // remove i from children of p
        head[p] = next[i];

        // start dfs on child node i
        stack[++top] = i;
      }
    }

    return k;
  }

  void SparsityInternal::init_ata(const int *post, int *w, int **head, int **next) const {
    int i, k, p, m = ncol_, n = nrow_;
    const int *ATp = &colind_.front();
    const int *ATi = &row_.front();
    *head = w+4*n, *next = w+5*n+1;

    // invert post
    for (k=0; k<n; ++k)
      w[post[k]] = k;

    for (i=0; i<m; ++i) {
      for (k=n, p=ATp[i]; p<ATp[i+1]; ++p)
        k = std::min(k, w[ATi[p]]);

      // place col i in linked list k
      (*next)[i] = (*head)[k];
      (*head)[k] = i ;
    }
  }

#define HEAD(k, j) (ata ? head[k] : j)
#define NEXT(J)   (ata ? next[J] : -1)
  std::vector<int> SparsityInternal::counts(const int *parent, const int *post, int ata) const {
    int i, j, k, n, m, J, s, p, q, jleaf=0, *maxfirst, *prevleaf, *ancestor,
        *head = 0, *next = 0, *first;

    m = nrow_;
    n = ncol_;
    s = 4*n + (ata ? (n+m+1) : 0);

    // allocate result
    vector<int> rowcount(n);
    vector<int>& delta = rowcount;

    // get workspace
    vector<int> w(s);

    // AT = A'
    Sparsity AT = transpose();

    ancestor = &w.front();
    maxfirst = &w.front()+n;
    prevleaf = &w.front()+2*n;
    first = &w.front()+3*n;

    // clear workspace w[0..s-1]
    for (k=0; k<s; ++k)
      w[k] = -1;

    // find first[j]
    for (k=0; k<n; ++k) {
      j = post[k];

      // delta[j]=1 if j is a leaf
      delta[j] = (first[j] == -1) ? 1 : 0;

      for (; j!=-1 && first[j] == -1; j=parent[j])
        first[j] = k;
    }

    const int* ATp = &AT.colind().front();
    const int* ATi = &AT.row().front();
    if (ata) AT->init_ata(post, &w.front(), &head, &next);

    // each node in its own set
    for (i=0; i<n; ++i)
      ancestor[i] = i;

    for (k=0; k<n; ++k) {
      // j is the kth node in postordered etree
      j = post[k];

      // j is not a root
      if (parent[j] != -1)
        delta[parent[j]]--;

      // J=j for LL'=A case
      for (J=HEAD(k, j); J != -1; J=NEXT(J)) {
        for (p = ATp[J]; p<ATp[J+1]; ++p) {
          i = ATi[p] ;
          q = leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);

          // A(i, j) is in skeleton
          if (jleaf >= 1)
            delta[j]++ ;

          // account for overlap in q
          if (jleaf == 2)
            delta[q]-- ;
        }
      }
      if (parent[j] != -1)
        ancestor[j] = parent[j] ;
    }

    // sum up delta's of each child
    for (j = 0 ; j < n ; ++j) {
      if (parent[j] != -1)
        rowcount[parent[j]] += rowcount[j] ;
    }

    // success
    return rowcount;
  }
#undef HEAD
#undef NEXT

  int SparsityInternal::wclear(int mark, int lemax, int *w, int n) {
    int k ;
    if (mark < 2 || (mark + lemax < 0)) {
        for (k = 0 ; k < n ; k++) if (w[k] != 0) w[k] = 1 ;
        mark = 2 ;
      }
    return (mark) ;     /* at this point, w[0..n-1] < mark holds */
  }

  int SparsityInternal::diag(int i, int j, double aij, void *other) {
    return (i != j) ;
  }

#define CS_FLIP(i) (-(i)-2)

  std::vector<int> SparsityInternal::approximateMinimumDegree(int order) const {

    int *Cp, *Ci, *last, *len, *nv, *next, *head, *elen, *degree, *w;
    int *hhead, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1;
    int k2, k3, jlast, ln, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi;
    int ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q;

    unsigned int h;

    //-- Construct matrix C -----------------------------------------------
    Sparsity AT = transpose() ;              // compute A'

    int m = nrow_;
    int n = ncol_;
    // find dense threshold:
    int dense = std::max(16, 10 * static_cast<int>(sqrt(static_cast<double>(n))));
    dense = std::min(n-2, dense);
    Sparsity C;
    if (order == 1 && n == m) {
      // C = A+A
      C = patternCombine(AT, false, false);
    } else if (order==2) {

      // drop dense rows from AT
      vector<int>& AT_colind = AT.colindRef();
      vector<int>& AT_row = AT.rowRef();
      for (p2=0, j=0; j<m; ++j) {

        // row j of AT starts here
        p = AT_colind[j];

        // new row j starts here
        AT_colind[j] = p2;

        // skip dense row j
        if (AT_colind[j+1] - p > dense)
          continue ;

        for ( ; p < AT_colind[j+1] ; p++)
          AT_row[p2++] = AT_row[p];
      }

      // finalize AT
      AT_colind[m] = p2;

      // Resize row vector
      AT_row.resize(p2);

      // A2 = AT'
      Sparsity A2 = AT->transpose();

      // C=A'*A with no dense cols
      C = AT->multiply(A2);
    } else {
      // C=A'*A
      C = AT->multiply(shared_from_this<Sparsity>());
    }

    // Free memory
    AT = Sparsity();
    // drop diagonal entries
    C->drop(diag, 0);

    Cp = &C.colindRef().front();
    cnz = Cp[n] ;

    // allocate result
    vector<int> P(n+1);

    // get workspace
    vector<int> W(8*(n+1));

    len = &W.front();
    nv = &W.front()  + (n+1);
    next = &W.front() + 2*(n+1);
    head = &W.front() + 3*(n+1);
    elen = &W.front() + 4*(n+1);
    degree = &W.front() + 5*(n+1);
    w= &W.front() + 6*(n+1);
    hhead = &W.front() + 7*(n+1);

    // use P as workspace for last
    last = &P.front();

    // --- Initialize quotient graph ----------------------------------------
    for (k=0; k<n; ++k)
      len[k] = Cp[k+1]-Cp[k];

    len[n] = 0;
    nzmax = C.size();
    Ci = &C.rowRef().front() ;
    for (i=0; i<=n; ++i) {
      // degree list i is empty
      head[i] = -1;
      last[i] = -1;
      next[i] = -1;

      // hash list i is empty
      hhead[i] = -1;

      // node i is just one node
      nv[i] = 1;

      // node i is alive
      w[i] = 1;

      // Ek of node i is empty
      elen[i] = 0;

      // degree of node i
      degree[i] = len[i];
    }

    // clear w
    mark = wclear(0, 0, w, n);

    // n is a dead element
    elen[n] = -2;

    // n is a root of assembly tree
    Cp[n] = -1;

    // n is a dead element
    w[n] = 0;

    // --- Initialize degree lists ------------------------------------------
    for (i = 0; i<n; ++i) {
      d = degree[i];

      // node i is empty
      if (d == 0) {
        // element i is dead
        elen[i] = -2;
        nel++;

        // i is a root of assembly tree
        Cp[i] = -1;
        w[i] = 0;
      } else if (d > dense) { // node i is dense
        // absorb i into element n
        nv[i] = 0;

        // node i is dead
        elen[i] = -1;
        nel++;
        Cp[i] = CS_FLIP(n) ;
        nv[n]++ ;
      } else {
        if (head[d] != -1)
          last[head[d]] = i;

        // put node i in degree list d
        next[i] = head[d];
        head[d] = i;
      }
    }

    // while (selecting pivots) do
    while (nel < n) {
      // --- Select node of minimum approximate degree --------------------
      for (k = -1 ; mindeg < n && (k = head[mindeg]) == -1; mindeg++) {}

      if (next[k] != -1) last[next[k]] = -1 ;

      // remove k from degree list
      head[mindeg] = next[k];

      // elenk = |Ek|
      elenk = elen[k];

      // # of nodes k represents
      nvk = nv[k];

      // nv[k] nodes of A eliminated
      nel += nvk;

      // --- Garbage collection -------------------------------------------
      if (elenk > 0 && cnz + mindeg >= nzmax) {
        for (j=0; j<n; ++j) {
          // j is a live node or element
          if ((p = Cp[j]) >= 0) {
            // save first entry of object
            Cp[j] = Ci[p];

            // first entry is now CS_FLIP(j)
            Ci[p] = CS_FLIP(j);
          }
        }

        // scan all of memory
        for (q = 0, p = 0 ; p < cnz ;) {
          // found object j
          if ((j = CS_FLIP(Ci[p++])) >= 0) {
            // restore first entry of object
            Ci[q] = Cp[j];

            // new pointer to object j
            Cp[j] = q++;
            for (k3 = 0 ; k3 < len[j]-1 ; k3++)
              Ci[q++] = Ci[p++] ;
          }
        }

        // Ci[cnz...nzmax-1] now free
        cnz = q;
      }

      // --- Construct new element ----------------------------------------
      dk = 0 ;

      // flag k as in Lk
      nv[k] = -nvk;
      p = Cp[k] ;

      // do in place if elen[k] == 0
      pk1 = (elenk == 0) ? p : cnz;
      pk2 = pk1 ;
      for (k1 = 1 ; k1 <= elenk + 1 ; k1++) {
        if (k1 > elenk) {
          // search the nodes in k
          e = k;

          // list of nodes starts at Ci[pj]
          pj = p;

          // length of list of nodes in k
          ln = len[k] - elenk;
        } else {
          // search the nodes in e
          e = Ci[p++];
          pj = Cp[e] ;

          // length of list of nodes in e
          ln = len[e];
        }

        for (k2=1; k2<=ln ; ++k2) {
          i = Ci[pj++] ;

          // node i dead, or seen
          if ((nvi = nv[i]) <= 0) continue;

          // degree[Lk] += size of node i
          dk += nvi;

          // negate nv[i] to denote i in Lk
          nv[i] = -nvi;

          // place i in Lk
          Ci[pk2++] = i;

          if (next[i] != -1)
            last[next[i]] = last[i];

          // remove i from degree list
          if (last[i] != -1) {
            next[last[i]] = next[i] ;
          } else {
            head[degree[i]] = next[i] ;
          }
        }

        if (e != k) {
          // absorb e into k
          Cp[e] = CS_FLIP(k);

          // e is now a dead element
          w[e] = 0;
        }
      }

      // Ci[cnz...nzmax] is free
      if (elenk != 0)
        cnz = pk2;

      // external degree of k - |Lk\i|
      degree[k] = dk;

      // element k is in Ci[pk1..pk2-1]
      Cp[k] = pk1;
      len[k] = pk2 - pk1 ;

      // k is now an element
      elen[k] = -2;

      // --- Find set differences -----------------------------------------

      // clear w if necessary
      mark = wclear(mark, lemax, w, n);

      // scan 1: find |Le\Lk|
      for (pk = pk1 ; pk < pk2 ; ++pk) {
        i = Ci[pk] ;

        // skip if elen[i] empty
        if ((eln = elen[i]) <= 0)
          continue;

        // nv[i] was negated
        nvi = -nv[i];

        wnvi = mark - nvi ;

        // scan Ei
        for (p = Cp[i] ; p <= Cp[i] + eln - 1 ; ++p) {
          e = Ci[p];
          if (w[e] >= mark) {
            // decrement |Le\Lk|
            w[e] -= nvi;
          } else if (w[e] != 0) {        /* ensure e is a live element */
            w[e] = degree[e] + wnvi ; /* 1st time e seen in scan 1 */
          }
        }
      }

      // --- Degree update ------------------------------------------------
      // scan2: degree update
      for (pk = pk1 ; pk < pk2 ; ++pk) {
        // consider node i in Lk
        i = Ci[pk];
        p1 = Cp[i] ;
        p2 = p1 + elen[i] - 1 ;
        pn = p1 ;

        // scan Ei
        for (h = 0, d = 0, p = p1 ; p <= p2 ; p++) {
          e = Ci[p] ;

          // e is an unabsorbed element
          if (w[e] != 0) {
            // dext = |Le\Lk|
            dext = w[e] - mark;
            if (dext > 0) {
              // sum up the set differences
              d += dext;

              // keep e in Ei
              Ci[pn++] = e;

              // compute the hash of node i
              h += e;

            } else {
              // aggressive absorb. e->k
              Cp[e] = CS_FLIP(k);

              // e is a dead element
              w[e] = 0;
            }
          }
        }

        // elen[i] = |Ei|
        elen[i] = pn - p1 + 1;
        p3 = pn ;
        p4 = p1 + len[i] ;

        // prune edges in Ai
        for (p = p2 + 1 ; p < p4 ; p++) {
          j = Ci[p] ;

          // node j dead or in Lk
          if ((nvj = nv[j]) <= 0)
            continue;

          // degree(i) += |j|
          d += nvj;

          // place j in node list of i
          Ci[pn++] = j;

          // compute hash for node i
          h += j;
        }

        // check for mass elimination
        if (d == 0) {
          // absorb i into k
          Cp[i] = CS_FLIP(k);
          nvi = -nv[i];

          // |Lk| -= |i|
          dk -= nvi;

          // |k| += nv[i]
          nvk += nvi;
          nel += nvi;
          nv[i] = 0 ;

          // node i is dead
          elen[i] = -1;
        } else {
          // update degree(i)
          degree[i] = std::min(degree[i], d);

          // move first node to end
          Ci[pn] = Ci[p3];

          // move 1st el. to end of Ei
          Ci[p3] = Ci[p1];

          // add k as 1st element in of Ei
          Ci[p1] = k;

          // new len of adj. list of node i
          len[i] = pn - p1 + 1;

          // finalize hash of i
          h %= n;

          // place i in hash bucket
          next[i] = hhead[h];
          hhead[h] = i ;

          // save hash of i in last[i]
          last[i] = h;
        }
      } // scan2 is done

      // finalize |Lk|
      degree[k] = dk;
      lemax = std::max(lemax, dk);

      // clear w
      mark = wclear(mark+lemax, lemax, w, n);

      // --- Supernode detection ------------------------------------------
      for (pk = pk1 ; pk < pk2 ; pk++) {
        i = Ci[pk] ;

        // skip if i is dead
        if (nv[i] >= 0)
          continue;

        // scan hash bucket of node i
        h = last[i];
        i = hhead[h];

        // hash bucket will be empty
        hhead[h] = -1;
        for ( ; i != -1 && next[i] != -1 ; i = next[i], mark++) {
          ln = len[i] ;
          eln = elen[i] ;
          for (p = Cp[i]+1 ; p <= Cp[i] + ln-1 ; p++)
            w[Ci[p]] = mark;

          jlast = i;

          // compare i with all j
          for (j = next[i] ; j != -1 ;) {
            ok = (len[j] == ln) && (elen[j] == eln) ;
            for (p = Cp[j] + 1 ; ok && p <= Cp[j] + ln - 1 ; p++) {
              if (w[Ci[p]] != mark) ok = 0 ;    /* compare i and j*/
            }

            // i and j are identical
            if (ok) {
              // absorb j into i
              Cp[j] = CS_FLIP(i);
              nv[i] += nv[j] ;
              nv[j] = 0;

              // node j is dead
              elen[j] = -1;

              // delete j from hash bucket
              j = next[j];
              next[jlast] = j ;
            } else {
              // j and i are different
              jlast = j;
              j = next[j] ;
            }
          }
        }
      }

      //  --- Finalize new element------------------------------------------
      // finalize Lk
      for (p = pk1, pk = pk1 ; pk < pk2 ; pk++) {
        i = Ci[pk] ;

        // skip if i is dead
        if ((nvi = -nv[i]) <= 0)
          continue;

        // restore nv[i]
        nv[i] = nvi;

        // compute external degree(i)
        d = degree[i] + dk - nvi ;
        d = std::min(d, n - nel - nvi);
        if (head[d] != -1)
          last[head[d]] = i;

        // put i back in degree list
        next[i] = head[d];
        last[i] = -1 ;
        head[d] = i ;

        // find new minimum degree
        mindeg = std::min(mindeg, d);
        degree[i] = d ;

        // place i in Lk
        Ci[p++] = i;
      }

      // # nodes absorbed into k
      nv[k] = nvk;

      // length of adj list of element k
      if ((len[k] = p-pk1) == 0) {
        // k is a root of the tree
        Cp[k] = -1;

        // k is now a dead element
        w[k] = 0;
      }

      // free unused space in Lk
      if (elenk != 0)
        cnz = p;
    }

    // --- Postordering -----------------------------------------------------

    // fix assembly tree
    for (i=0; i<n; ++i)
      Cp[i] = CS_FLIP(Cp[i]);

    for (j = 0 ; j <= n ; j++)
      head[j] = -1 ;

    // place unordered nodes in lists
    for (j = n ; j >= 0 ; j--) {
      // skip if j is an element
      if (nv[j] > 0) continue;

      // place j in list of its parent
      next[j] = head[Cp[j]];
      head[Cp[j]] = j;
    }

    // place elements in lists
    for (e = n ; e >= 0 ; e--) {
      // skip unless e is an element
      if (nv[e] <= 0) continue;
      if (Cp[e] != -1) {
        // place e in list of its parent
        next[e] = head[Cp[e]];
        head[Cp[e]] = e ;
      }
    }

    // postorder the assembly tree
    for (k = 0, i = 0 ; i <= n ; i++) {
      if (Cp[i] == -1)
        k = depthFirstSearchAndPostorder(i, k, head, next, &P.front(), w) ;
    }

    return P;
  }

  int SparsityInternal::scatter(int j, std::vector<int>& w, int mark, Sparsity& C, int nz) const {
    int i, p;
    const int *Ap = &colind_.front();
    const int *Ai = &row_.front();
    int *Ci = &C.rowRef().front();

    for (p = Ap[j]; p<Ap[j+1]; ++p) {
      // A(i, j) is nonzero
      i = Ai[p];

      if (w[i] < mark) {
        // i is new entry in row j
        w[i] = mark;

        // add i to pattern of C(:, j)
        Ci[nz++] = i;
      }
    }
    return nz;
  }

  Sparsity SparsityInternal::multiply(const Sparsity& B) const {
    int nz = 0;
    casadi_assert_message(ncol_ == B.size1(), "Dimension mismatch.");
    int m = nrow_;
    int anz = colind_[ncol_];
    int n = B.size2();
    const int* Bp = &B.colind().front();
    const int* Bi = &B.row().front();
    int bnz = Bp[n];

    // get workspace
    vector<int> w(m);

    // allocate result
    Sparsity C = Sparsity::sparse(m, n);
    C.colindRef().resize(anz + bnz);

    int* Cp = &C.colindRef().front();
    for (int j=0; j<n; ++j) {
      if (nz+m > C.size()) {
        C.rowRef().resize(2*(C.size())+m);
      }

      // row j of C starts here
      Cp[j] = nz;
      for (int p = Bp[j] ; p<Bp[j+1] ; ++p) {
        nz = scatter(Bi[p], w, j+1, C, nz);
      }
    }

    // finalize the last row of C
    Cp[n] = nz;
    C.rowRef().resize(nz);

    // Success
    return C;
  }

  void SparsityInternal::prefactorize(int order, int qr, std::vector<int>& S_pinv,
                                      std::vector<int>& S_q, std::vector<int>& S_parent,
                                      std::vector<int>& S_cp, std::vector<int>& S_leftmost,
                                      int& S_m2, double& S_lnz, double& S_unz) const {
    int k;
    int n = ncol_;
    vector<int> post;

    // fill-reducing ordering
    if (order!=0) {
      S_q = approximateMinimumDegree(order);
    }

    // QR symbolic analysis
    if (qr) {
      Sparsity C;
      if (order!=0) {
        std::vector<int> pinv_tmp;
        C = permute(pinv_tmp, S_q, 0);
      } else {
        C = shared_from_this<Sparsity>();
      }

      // etree of C'*C, where C=A(:, q)
      S_parent = C->eliminationTree(1);

      post = postorder(S_parent, n);

      // row counts chol(C'*C)
      S_cp = C->counts(&S_parent.front(), &post.front(), 1);
      post.clear();

      C->vcount(S_pinv, S_parent, S_leftmost, S_m2, S_lnz);
      for (S_unz = 0, k = 0; k<n; k++)
        S_unz += S_cp[k];

      // int overflow guard
      casadi_assert(S_lnz >= 0);
      casadi_assert(S_unz >= 0);
    } else {
      // for LU factorization only
      S_unz = 4*(colind_[n]) + n ;

      // guess nnz(L) and nnz(U)
      S_lnz = S_unz;
    }
  }

  Sparsity SparsityInternal::getDiag(std::vector<int>& mapping) const {
    if (ncol_==nrow_) {
      // Return object
      Sparsity ret = Sparsity::sparse(0, 1);
      ret.reserve(std::min(size(), ncol_), ncol_);

      // Mapping
      mapping.clear();

      // Loop over nonzero
      for (int i=0; i<ncol_; ++i) {

        // Enlarge the return matrix
        ret.resize(i+1, 1);

        // Get to the right nonzero of the col
        int el = colind_[i];
        while (el<colind_[i+1] && row_[el]<i) {
          el++;
        }

        if (el>=size()) return ret;

        // Add element if nonzero on diagonal
        if (row_[el]==i) {
          ret.getNZ(i, 0);
          mapping.push_back(el);
        }
      }

      return ret;

    } else if (ncol_==1 || nrow_==1) {
      Sparsity trans;
      const SparsityInternal *sp;

      // Have a col vector
      if (ncol_ == 1) {
        sp = this;
      } else {
        trans = transpose();
        sp = static_cast<const SparsityInternal *>(trans.get());
      }

      // Return object
      mapping.clear();
      mapping.resize(size());

      std::vector<int> colind(sp->nrow_+1, 0);
      std::vector<int> row(sp->size());

      int i_prev = 0;

      // Loop over nonzero
      for (int k=0;k<size();k++) {
        mapping[k]=k; // mapping will just be a range(size())

        int i = sp->row_[k];
        std::fill(colind.begin()+i_prev+1, colind.begin()+i+1, k);
        row[k]=i;
        i_prev = i;
      }
      std::fill(colind.begin()+i_prev+1, colind.end(), size());

      return Sparsity(sp->nrow_, sp->nrow_, colind, row);
    } else {
      casadi_error("diag: wrong argument shape. Expecting square matrix or vector-like, but got "
                   << dimString() << " instead.");
    }
  }

  std::string SparsityInternal::dimString() const {
    std::stringstream ss;
    if (numel()==size()) {
      ss << nrow_ << "-by-" << ncol_ << " (dense)";
    } else {
      ss << nrow_ << "-by-" << ncol_ << " (" << size() << "/" << numel() << " nz)";
    }
    return ss.str();
  }

  Sparsity SparsityInternal::patternProductNew(const Sparsity& y) const {
    // Dimensions of the result
    int d1 = nrow_;
    int d2 = y.size2();

    // Quick return if both are dense
    if (isDense() && y.isDense()) {
      return !isEmpty() && !y.isEmpty() ? Sparsity::dense(d1, d2) :
        Sparsity::sparse(d1, d2);
    }

    // Quick return if first factor is diagonal
    if (isDiagonal()) return y;

    // Quick return if second factor is diagonal
    if (y.isDiagonal()) return shared_from_this<Sparsity>();

    // Direct access to the vectors
    const vector<int> &x_row = row_;
    const vector<int> &x_colind = colind_;
    const vector<int> &y_row = y.row();
    const vector<int> &y_colind = y.colind();

    // Sparsity pattern of the result
    vector<int> row, col;

    // Temporary vector for avoiding duplicate nonzeros
    vector<int> tmp(d1, -1);

    // Loop over the nonzeros of y
    for (int cc=0; cc<d2; ++cc) {
      for (int kk=y_colind[cc]; kk<y_colind[cc+1]; ++kk) {
        int rr = y_row[kk];

        // Loop over corresponding columns of x
        for (int kk1=x_colind[rr]; kk1<x_colind[rr+1]; ++kk1) {
          int rr1 = x_row[kk1];

          // Add to pattern if not already encountered
          if (tmp[rr1]!=cc) {
            tmp[rr1] = cc;
            row.push_back(rr1);
            col.push_back(cc);
          }
        }
      }
    }

    // Assemble sparsity pattern and return
    return Sparsity::triplet(d1, d2, row, col);
  }

  Sparsity SparsityInternal::patternProduct(const Sparsity& x_trans,
                                            vector< vector< pair<int, int> > >& mapping) const {
    // return object
    Sparsity ret = shared_from_this<Sparsity>().patternProduct(x_trans);

    // Get the vectors for the return pattern
    const vector<int>& ret_row = ret.row();
    const vector<int>& ret_colind = ret.colind();

    // Direct access to the arrays
    const vector<int> &x_col = x_trans.row();
    const vector<int> &y_colind = colind_;
    const vector<int> &y_row = row_;
    const vector<int> &x_rowind = x_trans.colind();

    // Clear the mapping
    mapping.resize(ret.size());

    // the entry of the matrix to be calculated
    vector< pair<int, int> > d;

    // loop over the col of the resulting matrix)
    for (int i=0; i<ncol_; ++i) {
      // Loop over nonzeros
      for (int el=ret_colind[i]; el<ret_colind[i+1]; ++el) {
        int j = ret_row[el];
        int el1 = y_colind[i];
        int el2 = x_rowind[j];
        d.clear();
        while (el1 < y_colind[i+1] && el2 < x_rowind[j+1]) { // loop over non-zero elements
          int j1 = y_row[el1];
          int i2 = x_col[el2];
          if (j1==i2) {
            d.push_back(pair<int, int>(el1++, el2++));
          } else if (j1<i2) {
            el1++;
          } else {
            el2++;
          }
        }
        mapping[el] = d;
      }
    }

    return ret;
  }

  bool SparsityInternal::isScalar(bool scalar_and_dense) const {
    return ncol_==1 && nrow_==1 && (!scalar_and_dense || size()==1);
  }

  bool SparsityInternal::isDense() const {
    return size() == numel();
  }

  bool SparsityInternal::isEmpty(bool both) const {
    return both ? ncol_==0 && nrow_==0 : ncol_==0 || nrow_==0;
  }

  bool SparsityInternal::isDiagonal() const {
    // Check if matrix is square
    if (ncol_ != nrow_) return false;

    // Check if correct number of non-zeros (one per col)
    if (size() != ncol_) return false;

    // Check that the row indices are correct
    for (int i=0; i<size(); ++i) {
      if (row_[i]!=i)
        return false;
    }

    // Make sure that the col indices are correct
    for (int i=0; i<ncol_; ++i) {
      if (colind_[i]!=i)
        return false;
    }

    // Diagonal if reached this point
    return true;
  }

  bool SparsityInternal::isSquare() const {
    return ncol_ == nrow_;
  }

  bool SparsityInternal::isSymmetric() const {
    return isTranspose(*this);
  }

  int SparsityInternal::sizeL() const {
    int nnz = 0;
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el = colind_[cc+1]-1; el>=colind_[cc] && row_[el]>=cc; --el) nnz++;
    }
    return nnz;
  }

  int SparsityInternal::sizeD() const {
    int nnz = 0;
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el = colind_[cc]; el < colind_[cc+1]; ++el) {
        nnz += row_[el]==cc;
      }
    }
    return nnz;
  }

  int SparsityInternal::sizeU() const {
    int nnz = 0;
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el = colind_[cc]; el < colind_[cc+1] && row_[el]<=cc; ++el) nnz ++;
    }
    return nnz;
  }

  std::pair<int, int> SparsityInternal::shape() const {
    return std::pair<int, int>(nrow_, ncol_);
  }

  vector<int> SparsityInternal::erase(const vector<int>& jj, const vector<int>& ii) {
    if (!inBounds(jj, nrow_)) {
      casadi_error("Slicing [jj, ii] out of bounds. Your jj contains " <<
                   *std::min_element(jj.begin(), jj.end()) << " up to " <<
                   *std::max_element(jj.begin(), jj.end()) <<
                   ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(ii, ncol_)) {
      casadi_error("Slicing [jj, ii] out of bounds. Your ii contains "
                   << *std::min_element(ii.begin(), ii.end()) << " up to "
                   << *std::max_element(ii.begin(), ii.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    // Mapping
    vector<int> mapping;

    // Quick return if no elements
    if (numel()==0)
      return mapping;

    // Reserve memory
    mapping.reserve(size());

    // Number of non-zeros
    int nz=0;

    // Cols to be erased
    vector<int>::const_iterator ie = ii.begin();

    // First and last index for the col
    int el_first=0, el_last=0;

    // Loop over cols
    for (int i=0; i<ncol_; ++i) {
      // Update beginning and end of non-zero indices
      el_first = el_last;
      el_last = colind_[i+1];

      // Is it a col that can be deleted
      bool deletable_col = ie!=ii.end() && *ie==i;
      if (deletable_col) {
        ie++;

        // Rows to be erased
        vector<int>::const_iterator je = jj.begin();

        // Loop over nonzero elements of the col
        for (int el=el_first; el<el_last; ++el) {
          // Row
          int j=row_[el];

          // Continue to the next row to skip
          for (; je!=jj.end() && *je<j; ++je) {}

          // Remove row if necessary
          if (je!=jj.end() && *je==j) {
            je++;
            continue;
          }

          // Save old nonzero for each new nonzero
          mapping.push_back(el);

          // Update row and increase nonzero counter
          row_[nz++] = j;
        }
      } else {
        // Loop over nonzero elements of the col
        for (int el=el_first; el<el_last; ++el) {
          // Row
          int j=row_[el];

          // Save old nonzero for each new nonzero
          mapping.push_back(el);

          // Update row and increase nonzero counter
          row_[nz++] = j;
        }
      }

      // Register last nonzero of the col
      colind_[i+1]=nz;
    }

    // Truncate row matrix
    row_.resize(nz);

    return mapping;
  }

  vector<int> SparsityInternal::getNZ(const vector<int>& jj, const vector<int>& ii) const {
    if (!inBounds(jj, nrow_)) {
      casadi_error("Slicing [jj, ii] out of bounds. Your jj contains "
                   << *std::min_element(jj.begin(), jj.end()) << " up to "
                   << *std::max_element(jj.begin(), jj.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(ii, ncol_)) {
      casadi_error("Slicing [jj, ii] out of bounds. Your ii contains "
                   << *std::min_element(ii.begin(), ii.end()) << " up to "
                   << *std::max_element(ii.begin(), ii.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    std::vector<int> jj_sorted;
    std::vector<int> jj_sorted_index;

    sort(jj, jj_sorted, jj_sorted_index);

    vector<int> ret(ii.size()*jj.size());

    int stride = jj.size();

    for (int i=0;i<ii.size();++i) {
      int it = ii[i];
      int el=colind_[it];
      for (int j=0;j<jj_sorted.size();++j) {
        int jt=jj_sorted[j];
        // Continue to the non-zero element
        for (; el<colind_[it+1] && row_[el]<jt; ++el) {}
        // Add the non-zero element, if there was an element in the location exists
        if (el<colind_[it+1] && row_[el]== jt) {
          ret[i*stride+jj_sorted_index[j]] = el;
        } else {
          ret[i*stride+jj_sorted_index[j]] = -1;
        }
      }
    }
    return ret;
  }

  Sparsity SparsityInternal::sub(const vector<int>& jj, const vector<int>& ii,
                                 vector<int>& mapping) const {
    if (!inBounds(jj, nrow_)) {
      casadi_error("Slicing [jj, ii] out of bounds. Your jj contains "
                   << *std::min_element(jj.begin(), jj.end()) << " up to "
                   << *std::max_element(jj.begin(), jj.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(ii, ncol_)) {
      casadi_error("Slicing [jj, ii] out of bounds. Your ii contains "
                   << *std::min_element(ii.begin(), ii.end()) << " up to "
                   << *std::max_element(ii.begin(), ii.end())
                   << ", which is outside of the matrix shape " << dimString() << ".");
    }

    if (static_cast<double>(ii.size())*static_cast<double>(jj.size()) > size()) {
      // Typical use case:
      // a = SX::sym("a", sp_diag(50000))
      // a[:, :]
      return sub2(jj, ii, mapping);
    } else {
      // Typical use case:
      // a = DMatrix.ones(1000, 1000)
      // a[[0, 1],[0, 1]]
      return sub1(jj, ii, mapping);
    }

  }

  Sparsity SparsityInternal::sub2(const vector<int>& jj, const vector<int>& ii,
                                  vector<int>& mapping) const {
    std::vector<int> jj_sorted;
    std::vector<int> jj_sorted_index;

    std::vector<int> ii_sorted;
    std::vector<int> ii_sorted_index;

    sort(jj, jj_sorted, jj_sorted_index, false);
    sort(ii, ii_sorted, ii_sorted_index, false);

    std::vector<int> jjlookup = lookupvector(jj_sorted, nrow_);

    // count the number of non-zeros
    int nnz = 0;

    // loop over the cols of the slice
    for (int i=0;i<ii.size();++i) {
      int it = ii_sorted[i];
      for (int el=colind_[it]; el<colind_[it+1]; ++el) { // loop over the non-zeros of the matrix
        int j = row_[el];
        int ji = jjlookup[j];
        if (ji!=-1) {
          int jv = jj_sorted[ji];
          while (ji>=0 && jv == jj_sorted[ji--]) nnz++;
        }
      }
    }

    mapping.resize(nnz);

    vector<int> cols(nnz);
    vector<int> rows(nnz);

    int k = 0;
    // loop over the col of the slice
    for (int i=0;i<ii.size();++i) {
      int it = ii_sorted[i];
      for (int el=colind_[it]; el<colind_[it+1]; ++el) { // loop over the non-zeros of the matrix
        int jt = row_[el];
        int ji = jjlookup[jt];
        if (ji!=-1) {
          int jv = jj_sorted[ji];
          while (ji>=0 && jv == jj_sorted[ji]) {
            rows[k] = jj_sorted_index[ji];
            cols[k] = ii_sorted_index[i];
            mapping[k] = el;
            k++;
            ji--;
          }
        }
      }
    }

    std::vector<int> sp_mapping;
    std::vector<int> mapping_ = mapping;
    Sparsity ret = Sparsity::triplet(jj.size(), ii.size(), rows, cols, sp_mapping);

    for (int i=0;i<mapping.size();++i)
      mapping[i] = mapping_[sp_mapping[i]];

    // Create sparsity pattern
    return ret;
  }

  Sparsity SparsityInternal::sub1(const vector<int>& jj, const vector<int>& ii,
                                  vector<int>& mapping) const {

    std::vector<int> jj_sorted;
    std::vector<int> jj_sorted_index;

    std::vector<int> ii_sorted;
    std::vector<int> ii_sorted_index;

    sort(jj, jj_sorted, jj_sorted_index, false);
    sort(ii, ii_sorted, ii_sorted_index, false);

    // count the number of non-zeros
    int nnz = 0;

    for (int i=0;i<ii.size();++i) {
      int it = ii_sorted[i];
      int el=colind_[it];
      for (int j=0;j<jj_sorted.size();++j) {
        int jt=jj_sorted[j];
        // Continue to the non-zero element
        for (; el<colind_[it+1] && row_[el]<jt; ++el) {}
        // Add the non-zero element, if there was an element in the location exists
        if (el<colind_[it+1] && row_[el]== jt) {
          nnz++;
        }
      }
    }

    mapping.resize(nnz);

    vector<int> cols(nnz);
    vector<int> rows(nnz);


    int k=0;
    for (int i=0;i<ii.size();++i) {
      int it = ii_sorted[i];
      int K = colind_[it];
      for (int j=0;j<jj_sorted.size();++j) {
        int jt=jj_sorted[j];
        // Continue to the non-zero element
        for (; K<colind_[it+1] && row_[K]<jt; ++K) {}
        // Add the non-zero element, if there was an element in the location exists
        if (K<colind_[it+1] && row_[K]== jt) {
          rows[k] = jj_sorted_index[j];
          cols[k] = ii_sorted_index[i];
          mapping[k] = K;
          k++;
        }
      }
    }

    std::vector<int> sp_mapping;
    std::vector<int> mapping_ = mapping;
    Sparsity ret = Sparsity::triplet(jj.size(), ii.size(), rows, cols, sp_mapping);

    for (int i=0;i<mapping.size();++i)
      mapping[i] = mapping_[sp_mapping[i]];

    // Create sparsity pattern
    return ret;
  }

  Sparsity SparsityInternal::patternCombine(const Sparsity& y, bool f0x_is_zero,
                                            bool function0_is_zero) const {
    static vector<unsigned char> mapping;
    return patternCombineGen1<false>(y, f0x_is_zero, function0_is_zero, mapping);
  }

  Sparsity SparsityInternal::patternCombine(const Sparsity& y, bool f0x_is_zero,
                                            bool function0_is_zero,
                                            vector<unsigned char>& mapping) const {
    return patternCombineGen1<true>(y, f0x_is_zero, function0_is_zero, mapping);
  }

  template<bool with_mapping>
  Sparsity SparsityInternal::patternCombineGen1(const Sparsity& y, bool f0x_is_zero,
                                                bool function0_is_zero,
                                                std::vector<unsigned char>& mapping) const {

    // Quick return if identical
    if (isEqual(y)) {
      if (with_mapping) {
        mapping.resize(y.size());
        fill(mapping.begin(), mapping.end(), 1 | 2);
      }
      return y;
    }

    if (f0x_is_zero) {
      if (function0_is_zero) {
        return patternCombineGen<with_mapping, true, true>(y, mapping);
      } else {
        return patternCombineGen<with_mapping, true, false>(y, mapping);
      }
    } else if (function0_is_zero) {
      return patternCombineGen<with_mapping, false, true>(y, mapping);
    } else {
      return patternCombineGen<with_mapping, false, false>(y, mapping);
    }
  }

  template<bool with_mapping, bool f0x_is_zero, bool function0_is_zero>
  Sparsity SparsityInternal::patternCombineGen(const Sparsity& y,
                                               vector<unsigned char>& mapping) const {

    // Assert dimensions
    casadi_assert_message(ncol_==y.size2() && nrow_==y.size1(), "Dimension mismatch");

    // Sparsity pattern of the argument
    const vector<int>& y_colind = y.colind();
    const vector<int>& y_row = y.row();

    // Sparsity pattern of the result
    vector<int> ret_colind(ncol_+1, 0);
    vector<int> ret_row;

    // Clear the mapping
    if (with_mapping) mapping.clear();

    // Loop over cols of both patterns
    for (int i=0; i<ncol_; ++i) {
      // Non-zero element of the two matrices
      int el1 = colind_[i];
      int el2 = y_colind[i];

      // End of the non-zero elements of the col for the two matrices
      int el1_last = colind_[i+1];
      int el2_last = y_colind[i+1];

      // Loop over the non-zeros of both matrices
      while (el1<el1_last || el2<el2_last) {
        // Get the rows
        int row1 = el1<el1_last ? row_[el1] : nrow_;
        int row2 = el2<el2_last ? y_row[el2] : nrow_;

        // Add to the return matrix
        if (row1==row2) { //  both nonzero
          ret_row.push_back(row1);
          if (with_mapping) mapping.push_back( 1 | 2);
          el1++; el2++;
        } else if (row1<row2) { //  only first argument is nonzero
          if (!function0_is_zero) {
            ret_row.push_back(row1);
            if (with_mapping) mapping.push_back(1);
          } else {
            if (with_mapping) mapping.push_back(1 | 4);
          }
          el1++;
        } else { //  only second argument is nonzero
          if (!f0x_is_zero) {
            ret_row.push_back(row2);
            if (with_mapping) mapping.push_back(2);
          } else {
            if (with_mapping) mapping.push_back(2 | 4);
          }
          el2++;
        }
      }

      // Save the index of the last nonzero on the col
      ret_colind[i+1] = ret_row.size();
    }

    // Return cached object
    return Sparsity(nrow_, ncol_, ret_colind, ret_row);
  }

  bool SparsityInternal::isEqual(const Sparsity& y) const {
    // Quick true if the objects are the same
    if (this == y.get()) return true;

    // Otherwise, compare the patterns
    return isEqual(y.size1(), y.size2(), y.colind(), y.row());
  }

  Sparsity SparsityInternal::patternInverse() const {
    // Quick return clauses
    if (isEmpty()) return Sparsity::dense(nrow_, ncol_);
    if (isDense()) return Sparsity::sparse(nrow_, ncol_);

    // Sparsity of the result
    std::vector<int> row_ret;
    std::vector<int> colind_ret=colind_;

    // Loop over cols
    for (int i=0;i<ncol_;++i) {
      // Update colind vector of the result
      colind_ret[i+1]=colind_ret[i]+nrow_-(colind_[i+1]-colind_[i]);

      // Counter of new row indices
      int j=0;

      // Loop over all nonzeros
      for (int k=colind_[i];k<colind_[i+1];++k) {

        // Try to reach current nonzero
        while (j<row_[k])  {
          // And meanwhile, add nonzeros to the result
          row_ret.push_back(j);
          j++;
        }
        j++;
      }
      // Process the remainder up to the row size
      while (j < nrow_)  {
        row_ret.push_back(j);
        j++;
      }
    }

    // Return result
    return Sparsity(nrow_, ncol_, colind_ret, row_ret);
  }


  bool SparsityInternal::isEqual(int nrow, int ncol, const std::vector<int>& colind,
                                 const std::vector<int>& row) const {
    // First check dimensions and number of non-zeros
    if (size()!=row.size() || ncol_!=ncol || nrow_!=nrow)
      return false;

    // Check if dense
    if (size()==numel())
      return true;

    // Check the number of non-zeros per col
    if (!equal(colind_.begin(), colind_.end(), colind.begin()))
      return false;

    // Finally check the row indices
    if (!equal(row_.begin(), row_.end(), row.begin()))
      return false;

    // Equal if reached this point
    return true;
  }

  void SparsityInternal::reserve(int nnz, int ncol) {
    row_.reserve(nnz);
    colind_.reserve(ncol+1);
  }

  void SparsityInternal::append(const SparsityInternal& sp) {
    // NOTE: Works also if this == &sp

    casadi_assert_message(ncol_ == sp.ncol_,
      "SparsityInternal::append(sp): column sizes must match but got " << ncol_ <<
      " for lhs, and " << sp.ncol_ << " for rhs.");

    // Get current number of non-zeros
    int sz = row_.size();

    // Add row indices
    row_.resize(sz + sp.row_.size());
    for (int i=sz; i<row_.size(); ++i)
      row_[i] = sp.row_[i-sz] + nrow_;

    // Update column indices
    colind_[1] = row_.size();

    // Update dimensions
    nrow_ += sp.nrow_;
  }

  void SparsityInternal::appendColumns(const SparsityInternal& sp) {
    // NOTE: Works also if this == &sp

    casadi_assert_message(nrow_== sp.nrow_,
      "SparsityInternal::appendColumns(sp): row sizes must match but got " << nrow_ <<
      " for lhs, and " << sp.nrow_ << " for rhs.");

    // Get current number of non-zeros
    int sz = row_.size();

    // Add row indices
    row_.resize(sz + sp.row_.size());
    for (int i=sz; i<row_.size(); ++i) row_[i] = sp.row_[i-sz];

    // Add column indices
    colind_.resize(ncol_ + sp.ncol_ + 1);
    for (int i = ncol_+1; i<colind_.size(); ++i)
      colind_[i] = sp.colind_[i-ncol_] + sz;

    // Update dimensions
    ncol_ += sp.ncol_;
  }

  void SparsityInternal::enlargeColumns(int ncol, const std::vector<int>& ii) {
    // Assert dimensions
    casadi_assert(ii.size() == ncol_);

    // Update dimensions
    ncol_ = ncol;

    // Sparsify the cols
    colind_.resize(ncol+1, size());

    // Quick return if matrix had no cols before
    if (ii.empty()) return;

    int ik=ii.back(); // need only to update from the last new index
    int nz=size(); // number of nonzeros up till this col
    for (int i=ii.size()-1; i>=0; --i) {
      // Update colindex for new cols
      for (; ik>ii[i]; --ik) {
        colind_[ik] = nz;
      }

      // Update non-zero counter
      nz = colind_[i];

      // Update colindex for old cols
      colind_[ii[i]] = nz;
    }

    // Append zeros to the beginning
    for (; ik>=0; --ik) {
      colind_[ik] = 0;
    }
  }

  void SparsityInternal::enlargeRows(int nrow, const std::vector<int>& jj) {
    // Assert dimensions
    casadi_assert(jj.size() == nrow_);

    // Update dimensions
    nrow_ = nrow;

    // Begin by sparsify the rows
    for (int k=0; k<row_.size(); ++k) {
      row_[k] = jj[row_[k]];
    }
  }

  Sparsity SparsityInternal::makeDense(std::vector<int>& mapping) const {
    mapping.resize(size());
    for (int i=0; i<ncol_; ++i) {
      for (int el=colind_[i]; el<colind_[i+1]; ++el) {
        int j = row_[el];
        mapping[el] = j + i*nrow_;
      }
    }

    return Sparsity::dense(nrow_, ncol_);
  }

  int SparsityInternal::getNZ(int rr, int cc) const {
    // If negative index, count from the back
    if (rr<0) rr += nrow_;
    if (cc<0) cc += ncol_;

    // Check consistency
    casadi_assert_message(rr>=0 && rr<nrow_, "Row index " << rr
                          << " out of bounds [0, " << nrow_ << ")");
    casadi_assert_message(cc>=0 && cc<ncol_, "Column index " << cc
                          << " out of bounds [0, " << ncol_ << ")");

    // Quick return if matrix is dense
    if (isDense()) return rr+cc*nrow_;

    // Quick return if past the end
    if (colind_[cc]==size() || (colind_[cc+1]==size() && row_.back()<rr)) return -1;

    // Find sparse element
    for (int ind=colind_[cc]; ind<colind_[cc+1]; ++ind) {
      if (row_[ind] == rr) {
        return ind;     // element exists
      } else if (row_[ind] > rr) {
        break;          // break at the place where the element should be added
      }
    }
    return -1;
  }

  Sparsity SparsityInternal::reshape(int nrow, int ncol) const {
    casadi_assert_message(numel() == nrow*ncol,
                          "reshape: number of elements must remain the same. Old shape is "
                          << dimString() << ". New shape is " << nrow << "x" << ncol
                          << "=" << nrow*ncol << ".");
    Sparsity ret = Sparsity::sparse(nrow, ncol);
    ret.reserve(size(), ncol);

    std::vector<int> col(size());
    std::vector<int> row(size());
    for (int i=0; i<ncol_; ++i) {
      for (int el=colind_[i]; el<colind_[i+1]; ++el) {
        int j = row_[el];

        // Element number
        int k_ret = j+i*nrow_;

        // Col and row in the new matrix
        int i_ret = k_ret/nrow;
        int j_ret = k_ret%nrow;
        col[el] = i_ret;
        row[el] = j_ret;
      }
    }
    return Sparsity::triplet(nrow, ncol, row, col);
  }

  void SparsityInternal::resize(int nrow, int ncol) {
    if (ncol != ncol_ || nrow != nrow_) {
      if (ncol < ncol_ || nrow < nrow_) {
        // Col and row index of the new
        vector<int> row_new, colind_new(ncol+1, 0);

        // Loop over the cols which may contain nonzeros
        int i;
        for (i=0; i<ncol_ && i<ncol; ++i) {
          // First nonzero element of the col
          colind_new[i] = row_new.size();

          // Record rows of the nonzeros
          for (int el=colind_[i]; el<colind_[i+1] && row_[el]<nrow; ++el) {
            row_new.push_back(row_[el]);
          }
        }

        // Save col-indices for the rest of the cols
        for (; i<ncol+1; ++i) {
          colind_new[i] = row_new.size();
        }

        // Save the sparsity
        ncol_ = ncol;
        nrow_ = nrow;
        row_.swap(row_new);
        colind_.swap(colind_new);

      } else {
        // Make larger: Very cheap operation
        ncol_ = ncol;
        nrow_ = nrow;
        colind_.resize(ncol_+1, size());
      }
    }
  }

  bool SparsityInternal::rowsSequential(bool strictly) const {
    for (int i=0; i<ncol_; ++i) {
      int lastrow = -1;
      for (int k=colind_[i]; k<colind_[i+1]; ++k) {

        // check if not in sequence
        if (row_[k] < lastrow)
          return false;

        // Check if duplicate
        if (strictly && row_[k] == lastrow)
          return false;

        // update last row of the col
        lastrow = row_[k];
      }
    }

    // sequential if reached this point
    return true;
  }

  void SparsityInternal::removeDuplicates(std::vector<int>& mapping) {
    casadi_assert(mapping.size()==size());

    // Nonzero counter without duplicates
    int k_strict=0;

    // Loop over cols
    for (int i=0; i<ncol_; ++i) {

      // Last row encountered on the col so far
      int lastrow = -1;

      // Save new col offset (cannot set it yet, since we will need the old value below)
      int new_colind = k_strict;

      // Loop over nonzeros (including duplicates)
      for (int k=colind_[i]; k<colind_[i+1]; ++k) {

        // Make sure that the rows appear sequentially
        casadi_assert_message(row_[k] >= lastrow, "rows are not sequential");

        // Skip if duplicate
        if (row_[k] == lastrow)
          continue;

        // update last row encounterd on the col
        lastrow = row_[k];

        // Update mapping
        mapping[k_strict] = mapping[k];

        // Update row index
        row_[k_strict] = row_[k];

        // Increase the strict nonzero counter
        k_strict++;
      }

      // Update col offset
      colind_[i] = new_colind;
    }

    // Finalize the sparsity pattern
    colind_[ncol_] = k_strict;
    row_.resize(k_strict);
    mapping.resize(k_strict);
  }

  void SparsityInternal::getElements(std::vector<int>& loc, bool col_major) const {

    // Element for each nonzero
    loc.resize(size());

    // Loop over cols
    for (int i=0; i<ncol_; ++i) {

      // Loop over the nonzeros
      for (int el=colind_[i]; el<colind_[i+1]; ++el) {

        // Get row
        int j = row_[el];

        // Get the element
        if (col_major) {
          loc[el] = j+i*nrow_;
        } else {
          loc[el] = i+j*ncol_;
        }
      }
    }
  }

  void SparsityInternal::getNZInplace(std::vector<int>& indices) const {
    // Quick return if no elements
    if (indices.empty()) return;

    // Make a sanity check
    int last=-1;
    for (vector<int>::iterator it=indices.begin(); it!=indices.end(); ++it) {
      if (*it>=0) {
        int el_col = *it % ncol_;
        int el_row = *it / ncol_;
        int el = nrow_*el_col + el_row;
        casadi_assert_message(el>=last, "Elements must be sorted col-wise in non-decreasing order");
        last = el;
      }
    }

    // Quick return if no elements
    if (last<0) return;

    // Iterator to input/output
    vector<int>::iterator it=indices.begin();
    while (*it<0) it++; // first non-ignored

    // Current element sought
    int el_col = *it % ncol_;
    int el_row = *it / ncol_;

    // Loop over cols
    for (int i=0; i<ncol_; ++i) {

      // Loop over the nonzeros
      for (int el=colind_[i]; el<colind_[i+1] && el_col<=i; ++el) {

        // Get row
        int j = row_[el];

        // Add leading elements not in pattern
        while (i>el_col || (i==el_col && j>el_row)) {
          // Mark as not found
          *it = -1;

          // Increase index and terminate if end of vector reached
          if (++it==indices.end()) return;

          // Next element sought
          el_col = *it % ncol_;
          el_row = *it / ncol_;
        }

        // Add elements in pattern
        while (i==el_col && j==el_row) {
          // Save element index
          *it = el;

          // Increase index and terminate if end of vector reached
          do {
            if (++it==indices.end()) return;
          } while (*it<0);

          // Next element sought
          el_col = *it % ncol_;
          el_row = *it / ncol_;

        }
      }
    }

    // Add trailing elements not in pattern
    fill(it, indices.end(), -1);
  }

  Sparsity SparsityInternal::unidirectionalColoring(const Sparsity& AT, int cutoff) const {

    // Allocate temporary vectors
    vector<int> forbiddenColors;
    forbiddenColors.reserve(ncol_);
    vector<int> color(ncol_, 0);

    // Access the sparsity of the transpose
    const vector<int>& AT_colind = AT.colind();
    const vector<int>& AT_row = AT.row();

    // Loop over cols
    for (int i=0; i<ncol_; ++i) {

      // Loop over nonzero elements
      for (int el=colind_[i]; el<colind_[i+1]; ++el) {

        // Get row
        int c = row_[el];

        // Loop over previous cols that have an element in row c
        for (int el_prev=AT_colind[c]; el_prev<AT_colind[c+1]; ++el_prev) {

          // Get the col
          int i_prev = AT_row[el_prev];

          // Escape loop if we have arrived at the current col
          if (i_prev>=i)
            break;

          // Get the color of the col
          int color_prev = color[i_prev];

          // Mark the color as forbidden for the current col
          forbiddenColors[color_prev] = i;
        }
      }

      // Get the first nonforbidden color
      int color_i;
      for (color_i=0; color_i<forbiddenColors.size(); ++color_i) {
        // Break if color is ok
        if (forbiddenColors[color_i]!=i) break;
      }
      color[i] = color_i;

      // Add color if reached end
      if (color_i==forbiddenColors.size()) {
        forbiddenColors.push_back(0);

        // Cutoff if too many colors
        if (forbiddenColors.size()>cutoff) {
          return Sparsity();
        }
      }
    }

    // Create return sparsity containing the coloring
    Sparsity ret = Sparsity::sparse(ncol_, forbiddenColors.size());
    vector<int>& colind = ret.colindRef();
    vector<int>& row = ret.rowRef();

    // Get the number of rows for each col
    for (int i=0; i<color.size(); ++i) {
      colind[color[i]+1]++;
    }

    // Cumsum
    for (int j=0; j<forbiddenColors.size(); ++j) {
      colind[j+1] += colind[j];
    }

    // Get row for each col
    row.resize(color.size());
    for (int j=0; j<row.size(); ++j) {
      row[colind[color[j]]++] = j;
    }

    // Swap index back one step
    for (int j=colind.size()-2; j>=0; --j) {
      colind[j+1] = colind[j];
    }
    colind[0] = 0;

    // Return the coloring
    return ret;
  }

  Sparsity SparsityInternal::starColoring2(int ordering, int cutoff) const {
    casadi_assert_warning(ncol_==nrow_,
                          "StarColoring requires a square matrix, but got "
                          << dimString() << ".");

    // TODO(Joel): What we need here, is a distance-2 smallest last ordering
    // Reorder, if necessary
    if (ordering!=0) {
      casadi_assert(ordering==1);

      // Ordering
      vector<int> ord = largestFirstOrdering();

      // Create a new sparsity pattern
      Sparsity sp_permuted = pmult(ord, true, true, true);

      // Star coloring for the permuted matrix
      Sparsity ret_permuted = sp_permuted.starColoring2(0);

      // Permute result back
      return ret_permuted.pmult(ord, true, false, false);
    }

    // Allocate temporary vectors
    vector<int> forbiddenColors;
    forbiddenColors.reserve(ncol_);
    vector<int> color(ncol_, -1);

    vector<int> firstNeighborP(ncol_, -1);
    vector<int> firstNeighborQ(ncol_, -1);
    vector<int> firstNeighborQ_el(ncol_, -1);

    vector<int> treated(ncol_, -1);
    vector<int> hub(sizeU(), -1);

    vector<int> Tmapping;
    transpose(Tmapping);

    vector<int> star(size());
    int k = 0;
    for (int i=0; i<ncol_; ++i) {
      for (int j_el=colind_[i]; j_el<colind_[i+1]; ++j_el) {
        int j = row_[j_el];
        if (i<j) {
          star[j_el] = k;
          star[Tmapping[j]] = k;
          k++;
        }
      }
    }



    int starID = 0;

    // 3: for each v \in V do
    for (int v=0; v<ncol_; ++v) {

      // 4: for each colored w \in N1(v) do
      for (int w_el=colind_[v]; w_el<colind_[v+1]; ++w_el) {
          int w = row_[w_el];
          int colorW = color[w];
          if (colorW==-1) continue;

          // 5: forbiddenColors[color[w]] <- v
          forbiddenColors[colorW] = v;

          // 6: (p, q) <- firstNeighbor[color[w]]
          int p = firstNeighborP[colorW];
          int q = firstNeighborQ[colorW];

          // 7: if p = v then    <   Case 1
          if (v==p) {

            // 8: if treated[q] != v then
            if (treated[q]!=v) {

              // 9: treat(v, q)  < forbid colors of neighbors of q

                // treat@2: for each colored x \in N1 (q) do
                for (int x_el=colind_[q]; x_el<colind_[q+1]; ++x_el) {
                  int x = row_[x_el];
                  if (color[x]==-1) continue;

                  // treat@3: forbiddenColors[color[x]] <- v
                  forbiddenColors[color[x]] = v;
                }

                // treat@4: treated[q] <- v
                treated[q] = v;

            }
            // 10: treat(v, w) < forbid colors of neighbors of w

              // treat@2: for each colored x \in N1 (w) do
              for (int x_el=colind_[w]; x_el<colind_[w+1]; ++x_el) {
                int x = row_[x_el];
                if (color[x]==-1) continue;

                // treat@3: forbiddenColors[color[x]] <- v
                forbiddenColors[color[x]] = v;
              }

              // treat@4: treated[w] <- v
              treated[w] = v;

          // 11: else
          } else {

            // 12: firstNeighbor[color[w]] <- (v, w)
            firstNeighborP[colorW] = v;
            firstNeighborQ[colorW] = w;
            firstNeighborQ_el[colorW] = w_el;

            // 13: for each colored vertex x \in N1 (w) do
            int x_el_end = colind_[w+1];
            int x, colorx;
            for (int x_el=colind_[w]; x_el < x_el_end; ++x_el) {
              x = row_[x_el];
              colorx = color[x];
              if (colorx==-1 || x==v) continue;

              // 14: if x = hub[star[wx]] then potential Case 2
              if (hub[star[x_el]]==x) {

                // 15: forbiddenColors[color[x]] <- v
                forbiddenColors[colorx] = v;

              }
            }
          }

      }

      // 16: color[v] <- min {c > 0 : forbiddenColors[c] != v}
      bool new_color = true;
      for (int color_i=0; color_i<forbiddenColors.size(); ++color_i) {
        // Break if color is ok
        if (forbiddenColors[color_i]!=v) {
          color[v] = color_i;
          new_color = false;
          break;
        }
      }

      // New color if reached end
      if (new_color) {
        color[v] = forbiddenColors.size();
        forbiddenColors.push_back(-1);

        // Cutoff if too many colors
        if (forbiddenColors.size()>cutoff) {
          return Sparsity();
        }
      }

      // 17: updateStars(v)

        // updateStars@2: for each colored w \in N1 (v) do
        for (int w_el=colind_[v]; w_el<colind_[v+1]; ++w_el) {
            int w = row_[w_el];
            int colorW = color[w];
            if (colorW==-1) continue;

            // updateStars@3: if exits x \in N1 (w) where x = v and color[x] = color[v] then
            bool check = false;
            int x;
            int x_el;
            for (x_el=colind_[w]; x_el<colind_[w+1]; ++x_el) {
              x = row_[x_el];
              if (x==v || color[x]!=color[v]) continue;
              check = true;
              break;
            }
            if (check) {

              // updateStars@4: hub[star[wx]] <- w
              int starwx = star[x_el];
              hub[starwx] = w;

              // updateStars@5: star[vw] <- star[wx]
              star[w_el]  = starwx;
              star[Tmapping[w_el]] = starwx;

            // updateStars@6: else
            } else {

              // updateStars@7: (p, q) <- firstNeighbor[color[w]]
              int p = firstNeighborP[colorW];
              int q = firstNeighborQ[colorW];
              int q_el = firstNeighborQ_el[colorW];

              // updateStars@8: if (p = v) and (q = w) then
              if (p==v && q!=w) {

                // updateStars@9: hub[star[vq]] <- v
                int starvq = star[q_el];
                hub[starvq] = v;

                // updateStars@10: star[vw] <- star[vq]
                star[w_el]  = starvq;
                star[Tmapping[w_el]] = starvq;

              // updateStars@11: else
              } else {

                // updateStars@12: starID <- starID + 1
                starID+= 1;

                // updateStars@13: star[vw] <- starID
                star[w_el] = starID;
                star[Tmapping[w_el]]= starID;

              }

            }

         }

    }

    // Create return sparsity containing the coloring
    Sparsity ret = Sparsity::sparse(ncol_, forbiddenColors.size());
    vector<int>& colind = ret.colindRef();
    vector<int>& row = ret.rowRef();

    // Get the number of rows for each col
    for (int i=0; i<color.size(); ++i) {
      colind[color[i]+1]++;
    }

    // Cumsum
    for (int j=0; j<forbiddenColors.size(); ++j) {
      colind[j+1] += colind[j];
    }

    // Get row for each col
    row.resize(color.size());
    for (int j=0; j<row.size(); ++j) {
      row[colind[color[j]]++] = j;
    }

    // Swap index back one step
    for (int j=colind.size()-2; j>=0; --j) {
      colind[j+1] = colind[j];
    }
    colind[0] = 0;

    // Return the coloring
    return ret;

  }


  Sparsity SparsityInternal::starColoring(int ordering, int cutoff) const {
    casadi_assert_warning(ncol_==nrow_, "StarColoring requires a square matrix, but got "
                          << dimString() << ".");
    // Reorder, if necessary
    if (ordering!=0) {
      casadi_assert(ordering==1);

      // Ordering
      vector<int> ord = largestFirstOrdering();

      // Create a new sparsity pattern
      Sparsity sp_permuted = pmult(ord, true, true, true);

      // Star coloring for the permuted matrix
      Sparsity ret_permuted = sp_permuted.starColoring(0);

      // Permute result back
      return ret_permuted.pmult(ord, true, false, false);
    }

    // Allocate temporary vectors
    vector<int> forbiddenColors;
    forbiddenColors.reserve(ncol_);
    vector<int> color(ncol_, -1);

    // 4: for i <- 1 to |V | do
    for (int i=0; i<ncol_; ++i) {

      // 5: for each w \in N1 (vi) do
      for (int w_el=colind_[i]; w_el<colind_[i+1]; ++w_el) {
        int w = row_[w_el];

        // 6: if w is colored then
        if (color[w]!=-1) {

          // 7: forbiddenColors[color[w]] <- v
          forbiddenColors[color[w]] = i;

        } // 8: end if

        // 9: for each colored vertex x \in N1 (w) do
        for (int x_el=colind_[w]; x_el<colind_[w+1]; ++x_el) {
          int x = row_[x_el];
          if (color[x]==-1) continue;

          // 10: if w is not colored then
          if (color[w]==-1) {

            //11: forbiddenColors[color[x]] <- vi
            forbiddenColors[color[x]] = i;

          } else { // 12: else

            // 13: for each colored vertex y \in N1 (x), y != w do
            for (int y_el=colind_[x]; y_el<colind_[x+1]; ++y_el) {
              int y = row_[y_el];
              if (color[y]==-1 || y==w) continue;

              // 14: if color[y] = color[w] then
              if (color[y]==color[w]) {

                // 15: forbiddenColors[color[x]] <- vi
                forbiddenColors[color[x]] = i;

                // 16: break
                break;

              } // 17: end if

            } // 18: end for

          } // 19: end if

        } // 20 end for

      } // 21 end for

      // 22: color[v] <- min {c > 0 : forbiddenColors[c] = v}
      bool new_color = true;
      for (int color_i=0; color_i<forbiddenColors.size(); ++color_i) {
        // Break if color is ok
        if (forbiddenColors[color_i]!=i) {
          color[i] = color_i;
          new_color = false;
          break;
        }
      }

      // New color if reached end
      if (new_color) {
        color[i] = forbiddenColors.size();
        forbiddenColors.push_back(-1);

        // Cutoff if too many colors
        if (forbiddenColors.size()>cutoff) {
          return Sparsity();
        }
      }

    } // 23 end for

    // Number of colors used
    int num_colors = forbiddenColors.size();

    // Return sparsity in sparse triplet format
    return Sparsity::triplet(ncol_, num_colors, range(color.size()), color);
  }

  std::vector<int> SparsityInternal::largestFirstOrdering() const {
    vector<int> degree = colind_;
    int max_degree = 0;
    for (int k=0; k<ncol_; ++k) {
      degree[k] = degree[k+1]-degree[k];
      max_degree = max(max_degree, 1+degree[k]);
    }
    degree.resize(ncol_);

    // Vector for binary sort
    vector<int> degree_count(max_degree+1, 0);
    for (vector<int>::const_iterator it=degree.begin(); it!=degree.end(); ++it) {
      degree_count.at(*it+1)++;
    }

    // Cumsum to get the offset for each degree
    for (int d=0; d<max_degree; ++d) {
      degree_count[d+1] += degree_count[d];
    }

    // Now a bucket sort
    vector<int> ordering(ncol_);
    for (int k=ncol_-1; k>=0; --k) {
      ordering[degree_count[degree[k]]++] = k;
    }

    // Invert the ordering
    vector<int>& reverse_ordering = degree_count; // reuse memory
    reverse_ordering.resize(ordering.size());
    copy(ordering.begin(), ordering.end(), reverse_ordering.rbegin());

    // Return the ordering
    return reverse_ordering;
  }

  Sparsity SparsityInternal::pmult(const std::vector<int>& p, bool permute_rows,
                                   bool permute_cols, bool invert_permutation) const {
    // Invert p, possibly
    vector<int> p_inv;
    if (invert_permutation) {
      p_inv.resize(p.size());
      for (int k=0; k<p.size(); ++k) {
        p_inv[p[k]] = k;
      }
    }
    const vector<int>& pp = invert_permutation ? p_inv : p;

    // Get cols
    vector<int> col = getCol();

    // Sparsity of the return matrix
    vector<int> new_row(row_.size()), new_col(col.size());

    // Possibly permute cols
    if (permute_cols) {
      // Assert dimensions
      casadi_assert(p.size()==ncol_);

      // Permute
      for (int k=0; k<col.size(); ++k) {
        new_col[k] = pp[col[k]];
      }

    } else {
      // No permutation of cols
      copy(col.begin(), col.end(), new_col.begin());
    }

    // Possibly permute rows
    if (permute_rows) {
      // Assert dimensions
      casadi_assert(p.size()==nrow_);

      // Permute
      for (int k=0; k<row_.size(); ++k) {
        new_row[k] = pp[row_[k]];
      }

    } else {
      // No permutation of rows
      copy(row_.begin(), row_.end(), new_row.begin());
    }

    // Return permuted matrix
    return Sparsity::triplet(nrow_, ncol_, new_row, new_col);
  }

  bool SparsityInternal::isTranspose(const SparsityInternal& y) const {
    // Assert dimensions and number of nonzeros
    if (ncol_!=y.nrow_ || nrow_!=y.ncol_ || size()!=y.size())
      return false;

    // Quick return if empty interior or dense
    if (size()==0 || isDense())
      return true;

    // Run algorithm on the pattern with the least number of rows
    if (nrow_>ncol_) return y.isTranspose(*this);

    // Index counter for columns of the possible transpose
    vector<int> y_col_count(y.ncol_, 0);

    // Loop over the columns
    for (int i=0; i<ncol_; ++i) {

      // Loop over the nonzeros
      for (int el=colind_[i]; el<colind_[i+1]; ++el) {

        // Get the row
        int j=row_[el];

        // Get the element of the possible transpose
        int el_y = y.colind_[j] + y_col_count[j]++;

        // Quick return if element doesn't exist
        if (el_y>=y.colind_[j+1]) return false;

        // Get the row of the possible transpose
        int j_y = y.row_[el_y];

        // Quick return if mismatch
        if (j_y != i) return false;
      }
    }

    // Transpose if reached this point
    return true;
  }

  bool SparsityInternal::isReshape(const SparsityInternal& y) const {
    // Quick true if the objects are the same or are empty
    if (this==&y || isEmpty()) return true;

    // If same number of rows, check if patterns are equal
    if (nrow_==y.nrow_)
      return isEqual(y.nrow_, y.ncol_, y.colind_, y.row_);

    // Assert dimensions and number of nonzeros
    if (ncol_*nrow_!=y.nrow_*y.ncol_ || size()!=y.size())
      return false;

    // Quick return if empty interior or dense
    if (size()==0 || isDense()) return true;

    // Loop over the elements
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el=colind_[cc]; el<colind_[cc+1]; ++el) {
        int rr=row_[el];

        // Get row and column of y
        int loc = rr+nrow_*cc;
        int rr_y = loc % y.nrow_;
        int cc_y = loc / y.nrow_;

        // Make sure matching
        if (rr_y != y.row_[el] || el<y.colind_[cc_y] || el>=y.colind_[cc_y+1])
          return false;
      }
    }

    // Reshape if reached this point
    return true;
  }

  void SparsityInternal::spy(std::ostream &stream) const {

    // Index counter for each column
    std::vector<int> cind = colind_;

    // Loop over rows
    for (int rr=0; rr<nrow_; ++rr) {

      // Loop over columns
      for (int cc=0; cc<ncol_; ++cc) {
        // Check if nonzero
        if (cind[cc]<colind_[cc+1] && row_[cind[cc]]==rr) {
          stream << "*";
          cind[cc]++;
        } else {
          stream << ".";
        }
      }

      // End of row
      stream << endl;
    }
  }

  void SparsityInternal::spyMatlab(const std::string& mfile_name) const {
    // Create the .m file
    ofstream mfile;
    mfile.open(mfile_name.c_str());

    // Header
    mfile << "% This function was automatically generated by CasADi" << endl;

    // Print dimensions
    mfile << "n = " << nrow_ << ";" << endl;
    mfile << "m = " << ncol_ << ";" << endl;

    // Matlab indices are one-based
    const int index_offset = 1;

    // Print cols
    mfile << "i = [";
    bool first = true;
    for (int i=0; i<ncol_; ++i) {
      for (int el=colind_[i]; el<colind_[i+1]; ++el) {
        if (!first) mfile << ", ";
        mfile << (i+index_offset);
        first = false;
      }
    }
    mfile << "];" << endl;

    // Print rows
    mfile << "j = [";
    first = true;
    for (vector<int>::const_iterator j=row_.begin(); j!=row_.end(); ++j) {
      if (!first) mfile << ", ";
      mfile << (*j+index_offset);
      first = false;
    }
    mfile << "];" << endl;

    // Generate matrix
    mfile << "A = sparse(i, j, ones(size(i)), m, n)';" << endl;

    // Issue spy command
    mfile << "spy(A);" << endl;

    mfile.close();
  }

  std::size_t SparsityInternal::hash() const {
    return hash_sparsity(nrow_, ncol_, colind_, row_);
  }

  bool SparsityInternal::isTril() const {
    // loop over columns
    for (int i=0; i<ncol_; ++i) {
      if (colind_[i] != colind_[i+1]) { // if there are any elements of the column
        // check row of the top-most element of the column
        int row = row_[colind_[i]];

        // not lower triangular if row>i
        if (row<i) return false;
      }
    }
    // all columns ok
    return true;
  }

  bool SparsityInternal::isTriu() const {
    // loop over columns
    for (int i=0; i<ncol_; ++i) {
      if (colind_[i] != colind_[i+1]) { // if there are any elements of the column
        // check row of the bottom-most element of the column
        int row = row_[colind_[i+1]-1];

        // not upper triangular if row>i
        if (row>i) return false;
      }
    }
    // all columns ok
    return true;
  }

  Sparsity SparsityInternal::getTril(bool includeDiagonal) const {
    vector<int> colind, row;
    colind.reserve(ncol_+1);
    colind.push_back(0);
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el=colind_[cc]; el<colind_[cc+1]; ++el) {
        int rr=row_[el];
        if (rr>cc || (includeDiagonal && rr==cc)) {
          row.push_back(rr);
        }
      }
      colind.push_back(row.size());
    }
    return Sparsity(nrow_, ncol_, colind, row);
  }

  Sparsity SparsityInternal::getTriu(bool includeDiagonal) const {
    vector<int> colind, row;
    colind.reserve(ncol_+1);
    colind.push_back(0);
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el=colind_[cc]; el<colind_[cc+1]; ++el) {
        int rr=row_[el];
        if (rr<cc || (includeDiagonal && rr==cc)) {
          row.push_back(rr);
        }
      }
      colind.push_back(row.size());
    }
    return Sparsity(nrow_, ncol_, colind, row);
  }

  std::vector<int> SparsityInternal::getLowerNZ() const {
    vector<int> ret;
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el = colind_[cc]; el<colind_[cc+1]; ++el) {
        if (row_[el]>=cc) {
          ret.push_back(el);
        }
      }
    }
    return ret;
  }

  std::vector<int> SparsityInternal::getUpperNZ() const {
    vector<int> ret;
    for (int cc=0; cc<ncol_; ++cc) {
      for (int el = colind_[cc]; el<colind_[cc+1] && row_[el]<=cc; ++el) {
        ret.push_back(el);
      }
    }
    return ret;
  }

  int SparsityInternal::bandwidthU() const {
    int bw = 0;
    for (int cc=0; cc<ncol_; ++cc) {
      if (colind_[cc] != colind_[cc+1]) { // if there are any elements of the column
        int rr = row_[colind_[cc]];
        bw = std::max(bw, cc-rr);
      }
    }
    return bw;
  }

  int SparsityInternal::bandwidthL() const {
    int bw = 0;
    for (int cc=0; cc<ncol_; ++cc) {
      if (colind_[cc] != colind_[cc+1]) { // if there are any elements of the column
        int rr = row_[colind_[cc+1]-1];
        bw = std::max(bw, rr-cc);
      }
    }
    return bw;
  }

} // namespace casadi


