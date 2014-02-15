/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "crs_sparsity_internal.hpp"
#include "sparsity_tools.hpp"
#include "../stl_vector_tools.hpp"
#include <climits>
#include <cstdlib>
#include <cmath>
#include "matrix.hpp"

//#include "../external_packages/ColPack/ReducedHeader.h"

using namespace std;

namespace CasADi{

  int CRSSparsityInternal::size() const{
    return col_.size();
  }
    
  int CRSSparsityInternal::numel() const{
    return nrow_*ncol_;
  }
    
  void CRSSparsityInternal::repr(ostream &stream) const{
    stream << "Compressed Row Storage: " << nrow_ << "-by-" << ncol_ << " matrix, " << col_.size() << " structural non-zeros";
  }

  void CRSSparsityInternal::sanityCheck(bool complete) const{
    casadi_assert_message(nrow_>=0 ,"CRSSparsityInternal: number of rows must be positive, but got " << nrow_ << ".");
    casadi_assert_message(ncol_ >=0,"CRSSparsityInternal: number of cols must be positive, but got " << ncol_ << ".");
    if (rowind_.size() != nrow_+1) {
      std::stringstream s;
      s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
      s << "  rowind.size() = nrow + 1, but got   rowind.size() = " << rowind_.size() << "   and   nrow = "  << nrow_ << std::endl;
      s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
      casadi_error(s.str());
    }
    if (complete) {
  
      if (rowind_.size()>0) {
        for (int k=1;k<rowind_.size();k++) {
          casadi_assert_message(rowind_[k]>=rowind_[k-1], "CRSSparsityInternal:Compressed Row Storage is not sane. rowind must be monotone. Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind).");
        }
      
        casadi_assert_message(rowind_[0]==0, "CRSSparsityInternal:Compressed Row Storage is not sane. First element of rowind must be zero. Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind).");
        if (rowind_[(rowind_.size()-1)]!=col_.size()) {
          std::stringstream s;
          s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
          s << "  rowind[lastElement] = col.size(), but got   rowind[lastElement] = " << rowind_[(rowind_.size()-1)] << "   and   col.size() = "  << col_.size() << std::endl;
          s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
          casadi_error(s.str());
        }
        if (col_.size()>nrow_*ncol_) {
          std::stringstream s;
          s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
          s << "  col.size() <= nrow * ncol, but got   col.size()  = " << col_.size() << "   and   nrow * ncol = "  << nrow_*ncol_ << std::endl;
          s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
          casadi_error(s.str());
        }
      }
      for (int k=0;k<col_.size();k++) {
        if (col_[k]>=ncol_ || col_[k] < 0) {
          std::stringstream s;
          s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
          s << "  0 <= col[i] < ncol for each i, but got   col[i] = " << col_[k] << "   and   ncol = "  << ncol_ << std::endl;
          s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
          casadi_error(s.str());
        }
      }
  
    }
  }


  void CRSSparsityInternal::print(ostream &stream) const{
    repr(stream);
    stream << endl;
    stream << "col:    " << col_ << endl;
    stream << "rowind: " << rowind_ << endl;
  }

  vector<int> CRSSparsityInternal::getRow() const{
    vector<int> row(size());
    for(int r=0; r<nrow_; ++r){
      for(int el = rowind_[r]; el < rowind_[r+1]; ++el){
        row[el] = r;
      }
    }
    return row;
  }

  CRSSparsity CRSSparsityInternal::transpose() const{
    // Dummy mapping
    vector<int> mapping;

    return transpose(mapping);
  }

  CRSSparsity CRSSparsityInternal::transpose(vector<int>& mapping, bool invert_mapping) const{
    // Get the sparsity of the transpose in sparse triplet form
    const vector<int>& trans_row = col_;
    vector<int> trans_col = getRow();

    // Create the sparsity pattern
    return sp_triplet(ncol_,nrow_,trans_row,trans_col,mapping,invert_mapping);

  }

  std::vector<int> CRSSparsityInternal::eliminationTree(bool ata) const{
    // Allocate result
    vector<int> parent(nrow_);
  
    // Allocate workspace 
    vector<int> ancestor(nrow_);
    vector<int> prev(ata ? ncol_ : 0, -1);
  
    // Loop over rows
    for(int k=0; k<nrow_; ++k){
      // Start with no parent or ancestor
      parent[k] = -1;
      ancestor[k] = -1;
    
      // Loop over nonzeros
      for(int p=rowind_[k]; p<rowind_[k+1]; ++p){
      
        // What is this?
        int i=ata ? (prev[col_[p]]) : (col_[p]);
      
        // Transverse from i to k
        while(i!=-1 && i<k){
        
          // Next i is the ancestor of i
          int inext = ancestor[i];

          // Path compression
          ancestor[i] = k;
        
          // No ancestor, parent is k
          if(inext==-1) 
            parent[i] = k;
        
          // Update i
          i=inext;
        }
      
        // What is this?
        if(ata){
          prev[col_[p]] = k;
        }
      }
    }
  
    return parent;
  
  }

  int CRSSparsityInternal::depthFirstSearch(int j, int top, std::vector<int>& xi, std::vector<int>& pstack, const std::vector<int>& pinv, std::vector<bool>& marked) const{
    int head = 0;
  
    // initialize the recursion stack
    xi[0] = j;
    while (head >= 0){
    
      // get j from the top of the recursion stack 
      j = xi[head];
      int jnew = !pinv.empty() ? (pinv[j]) : j;
      if (!marked[j]){
      
        // mark node j as visited
        marked[j]=true;
        pstack[head] = (jnew < 0) ? 0 : rowind_[jnew];
      }
    
      // node j done if no unvisited neighbors
      int done = 1;
      int p2 = (jnew < 0) ? 0 : rowind_[jnew+1];
    
      // examine all neighbors of j
      for(int p = pstack[head]; p< p2; ++p){

        // consider neighbor node i
        int i = col_[p];
      
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
      if(done){
        // remove j from the recursion stack
        head--;
      
        // and place in the output stack
        xi[--top] = j ;
      }
    }
    return (top) ;
  }

  int CRSSparsityInternal::stronglyConnectedComponents(std::vector<int>& p, std::vector<int>& r) const{
    // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
    vector<int> tmp;

    CRSSparsity AT = transpose();
  
    vector<int> xi(2*nrow_+1);
    vector<int>& Blk = xi;
  
    vector<int> pstack(nrow_+1);
  
    p.resize(nrow_);
    r.resize(nrow_+6);
  
    vector<bool> marked(nrow_,false);
  
    int top = nrow_;
  
    //first dfs(A) to find finish times (xi)
    for(int i = 0; i<nrow_; ++i){
      if(!marked[i]){
        top = depthFirstSearch(i, top, xi, pstack, tmp, marked);
      }
    }

    //restore A; unmark all nodes
    fill(marked.begin(),marked.end(),false);
  
    top = nrow_;
    int nb = nrow_;

    // dfs(A') to find strongly connnected comp 
    for(int k=0 ; k < nrow_ ; ++k){
      // get i in reverse order of finish times
      int i = xi[k];
    
      // skip node i if already ordered
      if(marked[i]) continue;
    
      // node i is the start of a component in p
      r[nb--] = top;
      top = AT.depthFirstSearch(i, top, p, pstack, tmp, marked);
    }
  
    // first block starts at zero; shift r up
    r[nb] = 0;
    for (int k = nb ; k <= nrow_ ; ++k) 
      r[k-nb] = r[k] ;
  
    // nb = # of strongly connected components
    nb = nrow_-nb;
  
    // sort each block in natural order
    for(int b = 0 ; b < nb ; b++){
      for (int k = r[b]; k<r[b+1] ; ++k) 
        Blk[p[k]] = b ;
    }
  
    // Get p; shift r down (side effect)
    for(int i=0; i<nrow_; ++i){
      p[r[Blk[i]]++] = i;
    }
  
    // Shift up r
    r.resize(nb+1);
    for(int i=nb; i>0; --i){
      r[i]=r[i-1];
    }
    r[0]=0;
  
    return nb;
  }

  void CRSSparsityInternal::breadthFirstSearch(int n, std::vector<int>& wi, std::vector<int>& wj, std::vector<int>& queue, const std::vector<int>& imatch, const std::vector<int>& jmatch, int mark) const{
    // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
    int head = 0, tail = 0, j, i, p, j2 ;
  
    // place all unmatched nodes in queue
    for(j=0; j<n; ++j){
      // skip j if matched
      if(imatch[j] >= 0) continue;
    
      // j in set C0 (R0 if transpose)
      wj[j] = 0;
    
      // place unmatched col j in queue
      queue[tail++] = j;
    }
  
    // quick return if no unmatched nodes
    if(tail == 0) return;
  
    CRSSparsity trans;
    const CRSSparsityInternal *C;
    if(mark == 1){
      C = this;
    } else {
      trans = transpose();
      C = static_cast<const CRSSparsityInternal *>(trans.get());
    }
  
    // while queue is not empty
    while (head < tail){
    
      // get the head of the queue
      j = queue[head++];
      for(p = C->rowind_[j] ; p < C->rowind_[j+1] ; p++){
        i = C->col_[p] ;
      
        // skip if i is marked
        if (wi [i] >= 0) continue;
      
        // i in set R1 (C3 if transpose)
        wi [i] = mark;
      
        // traverse alternating path to j2
        j2 = jmatch [i];
      
        // skip j2 if it is marked
        if(wj [j2] >= 0) continue;
      
        // j2 in set C1 (R3 if transpose)
        wj[j2] = mark;
      
        // add j2 to queue
        queue [tail++] = j2;
      }
    }
  }

  void CRSSparsityInternal::matched(int n, const std::vector<int>& wj, const std::vector<int>& imatch, std::vector<int>& p, std::vector<int>& q, std::vector<int>& cc, std::vector<int>& rr, int set, int mark){
    // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
    int kc = cc[set];
    int kr = rr[set-1] ;
    for(int j=0; j<n; ++j){
      // skip if j is not in C set 
      if (wj[j] != mark) continue;
    
      p[kr++] = imatch[j] ;
      q[kc++] = j ;
    }
  
    cc[set+1] = kc ;
    rr[set] = kr ;
  }

  void CRSSparsityInternal::unmatched(int m, const std::vector<int>& wi, std::vector<int>& p, std::vector<int>& rr, int set){
    // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
    int i, kr = rr[set] ;
    for (i=0; i<m; i++) 
      if (wi[i] == 0) 
        p[kr++] = i;
    
    rr[set+1] = kr;
  }

  int CRSSparsityInternal::rprune(int i, int j, double aij, void *other){
    // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
    vector<int> &rr = *static_cast<vector<int> *>(other);
    return (i >= rr[1] && i < rr[2]) ;
  }

  void CRSSparsityInternal::augmentingPath(int k, std::vector<int>& jmatch, int *cheap, std::vector<int>& w, int *js, int *is, int *ps) const{
    // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

    int found = 0, p, i = -1, head = 0, j ;
  
    // start with just node k in jstack
    js[0] = k ;
  
    while (head >= 0){
      // --- Start (or continue) depth-first-search at node j -------------
    
      // get j from top of jstack
      j = js[head];
    
      // 1st time j visited for kth path
      if (w [j] != k){
      
        // mark j as visited for kth path 
        w[j] = k;
        for(p = cheap [j] ; p < rowind_[j+1] && !found; ++p){
          i = col_[p] ;            /* try a cheap assignment (i,j) */
          found = (jmatch [i] == -1) ;
        }
      
        // start here next time j is traversed
        cheap[j] = p;
        if(found){
          // column j matched with row i
          is[head] = i;
        
          // end of augmenting path
          break;
        }
      
        // no cheap match: start dfs for j
        ps[head] = rowind_[j];
      }
    
      // --- Depth-first-search of neighbors of j -------------------------
      for(p = ps[head]; p<rowind_[j+1]; ++p){
      
        // consider row i
        i = col_[p];
      
        // skip jmatch [i] if marked
        if(w[jmatch[i]] == k) continue;
      
        // pause dfs of node j
        ps[head] = p + 1;
      
        // i will be matched with j if found
        is[head] = i;
      
        // start dfs at column jmatch [i]
        js[++head] = jmatch [i];
        break ;
      }
    
      // node j is done; pop from stack
      if(p == rowind_[j+1]) head--;
    } // augment the match if path found:
  
    if(found)
      for(p = head; p>=0; --p)
        jmatch[is[p]] = js[p];
  }

  void CRSSparsityInternal::maxTransversal(std::vector<int>& imatch, std::vector<int>& jmatch, CRSSparsity& trans, int seed) const{
    // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

    int n2 = 0, m2 = 0;
  
    // allocate result
    jmatch.resize(ncol_);
    imatch.resize(nrow_);
    vector<int> w(ncol_+nrow_);
  
    // count nonempty rows and columns
    int k=0;
    for(int j=0; j<nrow_; ++j){
      n2 += (rowind_[j] < rowind_[j+1]);
      for(int p=rowind_[j]; p < rowind_[j+1]; ++p){
        w[col_[p]] = 1;
      
        // count entries already on diagonal
        k += (j == col_[p]);
      }
    }
  
    // quick return if diagonal zero-free
    if(k == std::min(ncol_,nrow_)){
      int i;
      for(i=0; i<k; ++i) jmatch[i] = i;
      for(;    i<ncol_; ++i) jmatch[i] = -1;

      int j;
      for(j=0; j<k; ++j) imatch[j] = j;
      for(;    j<nrow_; ++j) imatch[j] = -1;
    }

    for(int i=0; i<ncol_; ++i) m2 += w[i];
  
    // transpose if needed
    if(m2 < n2 && trans.isNull())
      trans = transpose();
  
    // Get pointer to sparsity
    const CRSSparsityInternal* C = m2 < n2 ? static_cast<const CRSSparsityInternal*>(trans.get()) : this;
  
    std::vector<int>& Cjmatch = m2 < n2 ? imatch : jmatch;
    std::vector<int>& Cimatch = m2 < n2 ? jmatch : imatch;
  
    // get workspace 
    w.resize( 5 * C->nrow_);

    int *cheap = &w.front() + C->nrow_;
    int *js = &w.front() + 2*C->nrow_;
    int *is = &w.front() + 3*C->nrow_; 
    int *ps = &w.front() + 4*C->nrow_;

    // for cheap assignment
    for(int j=0; j<C->nrow_; ++j) 
      cheap[j] = C->rowind_[j];
  
    // all columns unflagged 
    for(int j=0; j<C->nrow_; ++j)
      w[j] = -1;
  
    // nothing matched yet
    for(int i=0; i<C->ncol_; ++i)
      Cjmatch[i] = -1;

    // q = random permutation 
    std::vector<int> q = randomPermutation(C->nrow_,seed);

    // augment, starting at column q[k]
    for(k=0; k<C->nrow_; ++k){
      C->augmentingPath(!q.empty() ? q[k]: k, Cjmatch, cheap, w, js, is, ps);
    }

    // find row match
    for(int j=0; j<C->nrow_; ++j)
      Cimatch[j] = -1;
  
    for(int i = 0; i<C->ncol_; ++i)
      if(Cjmatch [i] >= 0)
        Cimatch[Cjmatch[i]] = i;
  }

  int CRSSparsityInternal::dulmageMendelsohn(std::vector<int>& rowperm, std::vector<int>& colperm, std::vector<int>& rowblock, std::vector<int>& colblock, std::vector<int>& coarse_rowblock, std::vector<int>& coarse_colblock, int seed) const{
    // The transpose of the expression
    CRSSparsity trans;
  
    // Part 1: Maximum matching

    // row permutation 
    colperm.resize(ncol_);
  
    // column permutation 
    rowperm.resize(nrow_);
  
    // size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
    colblock.resize(ncol_+6);
  
    // size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
    rowblock.resize(nrow_+6);

    // coarse row decomposition
    coarse_colblock.resize(5);
    fill(coarse_colblock.begin(),coarse_colblock.end(),0);
  
    // coarse column decomposition
    coarse_rowblock.resize(5);
    fill(coarse_rowblock.begin(),coarse_rowblock.end(),0);

    // max transversal
    vector<int> imatch, jmatch;
    maxTransversal(imatch,jmatch,trans,seed);
  
    // Coarse decomposition
  
    // use colblock and rowblock as workspace
    vector<int>& wi = colblock;
    vector<int>& wj = rowblock;
  
    // unmark all cols for bfs
    for(int j=0; j<nrow_; ++j)
      wj[j] = -1;
  
    // unmark all rows for bfs
    for(int i=0; i<ncol_; ++i)
      wi[i] = -1 ;
  
    // find C1, R1 from C0
    breadthFirstSearch(nrow_, wi, wj, rowperm, imatch, jmatch, 1);

    // find R3, C3 from R0
    breadthFirstSearch(ncol_, wj, wi, colperm, jmatch, imatch, 3);

    // unmatched set C0
    unmatched(nrow_, wj, rowperm, coarse_rowblock, 0);

    // set R1 and C1
    matched(nrow_, wj, imatch, colperm, rowperm, coarse_rowblock, coarse_colblock, 1, 1);

    // set R2 and C2
    matched(nrow_, wj, imatch, colperm, rowperm, coarse_rowblock, coarse_colblock, 2, -1);

    // set R3 and C3
    matched(nrow_, wj, imatch, colperm, rowperm, coarse_rowblock, coarse_colblock, 3, 3);

    // unmatched set R0
    unmatched(ncol_, wi, colperm, coarse_colblock, 3);
  
    // --- Fine decomposition -----------------------------------------------
    // pinv=p'
    vector<int> pinv = invertPermutation(colperm);

    // C=A(p,q) (it will hold A(R2,C2))
    CRSSparsity C = permute(pinv, rowperm, 0);

    vector<int>& rowind_C = C.rowindRef();

    // delete cols C0, C1, and C3 from C 
    int nc = coarse_rowblock[3] - coarse_rowblock[2];
    if(coarse_rowblock[2] > 0){
      for(int j = coarse_rowblock[2]; j <= coarse_rowblock[3]; ++j)
        rowind_C[j-coarse_rowblock[2]] = rowind_C[j];
    }
    C->nrow_ = nc;

    C->rowind_.resize(nc+1);
    // delete rows R0, R1, and R3 from C
    if(coarse_colblock[2] - coarse_colblock[1] < ncol_){
      C->drop(rprune, &coarse_colblock);
      int cnz = rowind_C[nc];
      vector<int>& col_C = C->col_;
      if(coarse_colblock[1] > 0)
        for(int k=0; k<cnz; ++k)
          col_C[k] -= coarse_colblock[1];
    }
    C->col_.resize(C->rowind_.back());
    C->ncol_ = nc ;

    // find strongly connected components of C
    vector<int> scc_p, scc_r;
    int scc_nb = C->stronglyConnectedComponents(scc_p, scc_r);
  
    // --- Combine coarse and fine decompositions ---------------------------
  
    // C(ps,ps) is the permuted matrix
    vector<int> ps = scc_p;
  
    // kth block is rs[k]..rs[k+1]-1
    vector<int> rs = scc_r;
  
    // # of blocks of A(R2,C2)
    int nb1 = scc_nb;

    for(int k=0; k<nc; ++k)
      wj[k] = rowperm[ps[k] + coarse_rowblock[2]];
  
    for(int k=0; k<nc; ++k)
      rowperm[k + coarse_rowblock[2]] = wj[k];
  
    for(int k=0; k<nc; ++k)
      wi[k] = colperm[ps[k] + coarse_colblock[1]];
  
    for(int k=0; k<nc; ++k)
      colperm[k + coarse_colblock[1]] = wi[k];
  
    // create the fine block partitions
    int nb2 = 0;
    colblock[0] = rowblock[0] = 0;

    // leading coarse block A (R1, [C0 C1])
    if(coarse_rowblock[2] > 0)
      nb2++ ;
  
    // coarse block A (R2,C2)
    for(int k=0; k<nb1; ++k){
      // A (R2,C2) splits into nb1 fine blocks 
      colblock[nb2] = rs[k] + coarse_colblock[1];
      rowblock[nb2] = rs[k] + coarse_rowblock[2] ;
      nb2++ ;
    }
  
    if(coarse_colblock[2] < ncol_){
      // trailing coarse block A ([R3 R0], C3)
      colblock[nb2] = coarse_colblock[2];
      rowblock[nb2] = coarse_rowblock[3];
      nb2++ ;
    }
  
    colblock[nb2] = ncol_;
    rowblock[nb2] = nrow_ ;
  
    // Shrink colblock and rowblock
    colblock.resize(nb2+1);
    rowblock.resize(nb2+1);
    return nb2;
  }

  std::vector<int> CRSSparsityInternal::randomPermutation(int n, int seed){
    // Return object
    std::vector<int> p;
  
    // return p = empty (identity)
    if(seed==0) return p;
  
    // allocate result
    p.resize(n);
  
    for(int k=0; k<n; ++k) 
      p[k] = n-k-1;
  
    // return reverse permutation
    if(seed==-1) return p;
  
    // get new random number seed
    srand(seed);
  
    for(int k=0; k<n; ++k){
      // j = rand int in range k to n-1
      int j = k + (rand ( ) % (n-k));
    
      // swap p[k] and p[j]
      int t = p[j];
      p[j] = p[k];
      p[k] = t;
    }
  
    return p;
  }

  std::vector<int> CRSSparsityInternal::invertPermutation(const std::vector<int>& p){
    // pinv = p', or p = pinv'

    // allocate result
    vector<int> pinv(p.size());
  
    // invert the permutation
    for(int k=0; k<p.size(); ++k)
      pinv[p[k]] = k;
  
    // return result
    return pinv;
  }

  CRSSparsity CRSSparsityInternal::permute(const std::vector<int>& pinv, const std::vector<int>& q, int values) const{
    // alloc result
    CRSSparsity C = CRSSparsity(nrow_,ncol_);
  
    // Row offset
    vector<int>& rowind_C = C.rowindRef();
  
    // Column for each nonzero
    vector<int>& col_C = C.colRef();
    col_C.resize(size());

    int nz = 0;
    for(int k = 0; k<nrow_; ++k){
      // column k of C is column q[k] of A
      rowind_C[k] = nz;
    
      int j = !q.empty() ? (q[k]) : k;
    
      for(int t = rowind_[j]; t<rowind_[j+1]; ++t){
        col_C[nz++] = !pinv.empty() ? (pinv[col_[t]]) : col_[t] ;
      }
    }
  
    // finalize the last column of C
    rowind_C[nrow_] = nz;
    return C;
  }

  int CRSSparsityInternal::drop(int (*fkeep) (int, int, double, void *), void *other){
    int nz = 0;
  
    for(int j = 0; j<nrow_; ++j){
      // get current location of col j
      int p = rowind_[j];
    
      // record new location of col j
      rowind_[j] = nz;
      for ( ; p < rowind_[j+1] ; ++p){
        if (fkeep(col_[p], j, 1, other)){
          // keep A(i,j)
          col_[nz++] = col_[p] ;
        }
      }
    }
  
    // finalize A
    rowind_[nrow_] = nz;
    return nz ;
  }

  int CRSSparsityInternal::leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf){
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
    for (s = jprev ; s != q ; s = sparent)
      {
        sparent = ancestor [s] ;    /* path compression */
        ancestor [s] = q ;
      }
    return (q) ;                    /* q = least common ancester (jprev,j) */
  }

  int CRSSparsityInternal::vcount(std::vector<int>& pinv, std::vector<int>& parent, std::vector<int>& leftmost, int& S_m2, double& S_lnz) const{
    int i, k, p, pa;
    int n = nrow_, m = ncol_;
    const int* Ap = &rowind_.front();
    const int* Ai = &col_.front();

    // allocate pinv
    pinv.resize(m+n);
    fill(pinv.begin(),pinv.end(),0);
  
    // and leftmost
    leftmost.resize(m);
  
    // get workspace
    vector<int> w(m+3*n);

    int *next = &w.front();
    int *head = &w.front() + m;
    int *tail = &w.front() + m + n;
    int *nque = &w.front() + m + 2*n;

    // queue k is empty
    for(k = 0 ; k < n ; k++)
      head[k] = -1;
  
    for(k=0; k<n; ++k)
      tail[k] = -1;
  
    for(k=0; k<n; ++k)
      nque[k] = 0;
  
    for(i=0; i<m; ++i)
      leftmost[i] = -1;
  
    for(k=n-1; k>=0; --k){
      for(p=Ap[k]; p<Ap[k+1]; ++p){
        // leftmost[i] = min(find(A(i,:)))
        leftmost[Ai[p]] = k;
      }
    }
  
    // scan rows in reverse order
    for (i = m-1; i >= 0; i--){
      // row i is not yet ordered
      pinv[i] = -1;
      k = leftmost [i] ;
    
      // row i is empty
      if (k == -1) continue;
    
      // first row in queue k
      if(nque[k]++ == 0)
        tail[k] = i;
    
      // put i at head of queue k 
      next[i] = head[k];
      head[k] = i;
    }
    S_lnz = 0;
    S_m2 = m;
  
    // find row permutation and nnz(V)
    for(k=0; k<n; ++k){
      // remove row i from queue k
      i = head[k];
    
      // count V(k,k) as nonzero 
      S_lnz++;
    
      // add a fictitious row
      if(i < 0)
        i = S_m2++;
    
      // associate row i with V(:,k)
      pinv[i] = k;
    
      // skip if V(k+1:m,k) is empty
      if(--nque[k] <= 0) continue;
    
      // nque [k] is nnz (V(k+1:m,k))
      S_lnz += nque[k];
    
      // move all rows to parent of k
      if((pa = parent[k]) != -1){
        if(nque[pa] == 0)
          tail[pa] = tail[k];
      
        next[tail[k]] = head[pa] ;
        head[pa] = next[i] ;
        nque[pa] += nque[k] ;
      }
    }
    for(i=0; i<m ; ++i)
      if(pinv[i] < 0)
        pinv[i] = k++;
    
    pinv.resize(m);
    return 1;
  }

  std::vector<int> CRSSparsityInternal::postorder(const std::vector<int>& parent, int n){
    int j, k = 0, *head, *next, *stack ;
  
    // allocate result
    vector<int> post(n);
  
    // get workspace
    vector<int> w(3*n);
  
    head = &w.front() ;
    next = &w.front() + n ;
    stack = &w.front() + 2*n;
  
    // empty linked lists
    for(j=0; j<n; ++j)
      head[j] = -1;
  
    // traverse nodes in reverse order
    for (j=n-1; j>=0; --j){
      // j is a root
      if (parent [j] == -1) continue;
    
      // add j to list of its parent
      next[j] = head[parent[j]];
      head[parent[j]] = j ;
    }
  
    for(j=0; j<n; ++j){
      // skip j if it is not a root
      if (parent [j] != -1) continue;

      k = depthFirstSearchAndPostorder(j, k, head, next, &post.front(), stack);
    }
  
    // success; return post
    return post;
  }

  int CRSSparsityInternal::depthFirstSearchAndPostorder(int j, int k, int *head, const int *next, int *post, int *stack){
    int i, p, top = 0;
  
    // place j on the stack
    stack[0] = j;
  
    // while (stack is not empty)
    while(top >= 0){
      // p = top of stack
      p = stack[top];
    
      // i = youngest child of p 
      i = head[p];
      if (i == -1){
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

  void CRSSparsityInternal::init_ata(const int *post, int *w, int **head, int **next) const{
    int i, k, p, m = nrow_, n = ncol_;
    const int *ATp = &rowind_.front();
    const int *ATi = &col_.front();
    *head = w+4*n, *next = w+5*n+1;
  
    // invert post
    for(k=0; k<n; ++k)
      w[post[k]] = k;
  
    for(i=0; i<m; ++i){
      for(k=n, p=ATp[i]; p<ATp[i+1]; ++p)
        k = std::min(k, w[ATi[p]]);
    
      // place row i in linked list k
      (*next)[i] = (*head)[k];
      (*head)[k] = i ;
    }
  }

#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
  std::vector<int> CRSSparsityInternal::counts(const int *parent, const int *post, int ata) const{
    int i, j, k, n, m, J, s, p, q, jleaf, *maxfirst, *prevleaf, *ancestor, *head = NULL, *next = NULL, *first;

    m = ncol_;
    n = nrow_;
    s = 4*n + (ata ? (n+m+1) : 0);

    // allocate result
    vector<int> colcount(n);
    vector<int>& delta = colcount;
  
    // get workspace
    vector<int> w(s);
  
    // AT = A'
    CRSSparsity AT = transpose();

    ancestor = &w.front();
    maxfirst = &w.front()+n;
    prevleaf = &w.front()+2*n;
    first = &w.front()+3*n;
  
    // clear workspace w [0..s-1]
    for(k=0; k<s; ++k)
      w[k] = -1;
  
    // find first [j]
    for(k=0; k<n; ++k){
      j = post[k];

      // delta[j]=1 if j is a leaf
      delta[j] = (first[j] == -1) ? 1 : 0;

      for (; j!=-1 && first [j] == -1; j=parent[j])
        first[j] = k;
    }

    const int* ATp = &AT.rowind().front();
    const int* ATi = &AT.col().front();
    if (ata) AT->init_ata(post, &w.front(), &head, &next);
  
    // each node in its own set
    for(i=0; i<n; ++i)
      ancestor[i] = i;
  
    for(k=0; k<n; ++k){
      // j is the kth node in postordered etree
      j = post[k];
    
      // j is not a root
      if (parent [j] != -1)
        delta[parent [j]]--;
      
      // J=j for LL'=A case
      for(J=HEAD(k,j); J != -1; J=NEXT(J)){
        for(p = ATp[J]; p<ATp[J+1]; ++p){
          i = ATi [p] ;
          q = leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);

          // A(i,j) is in skeleton
          if(jleaf >= 1)
            delta[j]++ ;
        
          // account for overlap in q
          if(jleaf == 2)
            delta [q]-- ;
        }
      }
      if(parent[j] != -1)
        ancestor[j] = parent[j] ;
    }
  
    // sum up delta's of each child
    for(j = 0 ; j < n ; ++j){
      if (parent[j] != -1)
        colcount[parent [j]] += colcount[j] ;
    }
  
    // success
    return colcount;
  }
#undef HEAD
#undef NEXT

  int CRSSparsityInternal::wclear (int mark, int lemax, int *w, int n){
    int k ;
    if (mark < 2 || (mark + lemax < 0))
      {
        for (k = 0 ; k < n ; k++) if (w [k] != 0) w [k] = 1 ;
        mark = 2 ;
      }
    return (mark) ;     /* at this point, w [0..n-1] < mark holds */
  }

  int CRSSparsityInternal::diag (int i, int j, double aij, void *other){
    return (i != j) ;
  }

#define CS_FLIP(i) (-(i)-2)

  std::vector<int> CRSSparsityInternal::approximateMinimumDegree(int order) const{
  
    int *Cp, *Ci, *last, *len, *nv, *next, *head, *elen, *degree, *w;
    int *hhead, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1;
    int k2, k3, jlast, ln, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi;
    int ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, t;

    unsigned int h;

    //-- Construct matrix C -----------------------------------------------
    CRSSparsity AT = transpose() ;              // compute A'
  
    int m = ncol_;
    int n = nrow_;
    int dense = std::max(16, 10 * int(sqrt(double(n)))) ;   // find dense threshold
    dense = std::min(n-2, dense);
    CRSSparsity C;
    if(order == 1 && n == m){
      // C = A+A
      C = patternCombine(AT,false,false);
    } else if(order==2){
    
      // drop dense columns from AT
      vector<int>& AT_rowind = AT.rowindRef();
      vector<int>& AT_col = AT.colRef();
      for(p2=0, j=0; j<m; ++j){
      
        // column j of AT starts here
        p = AT_rowind[j];
      
        // new column j starts here
        AT_rowind[j] = p2;
      
        // skip dense col j
        if(AT_rowind[j+1] - p > dense)
          continue ;
      
        for ( ; p < AT_rowind[j+1] ; p++)
          AT_col[p2++] = AT_col[p];
      }
    
      // finalize AT
      AT_rowind[m] = p2;
    
      // Resize column vector
      AT_col.resize(p2);
    
      // A2 = AT'
      CRSSparsity A2 = AT->transpose();
    
      // C=A'*A with no dense rows
      C = AT->multiply(A2);
    } else {
      // C=A'*A
      C = AT->multiply(shared_from_this<CRSSparsity>());
    }
  
    // Free memory
    AT = CRSSparsity();
    // drop diagonal entries
    C->drop(diag, NULL);
  
    Cp = &C.rowindRef().front();
    cnz = Cp[n] ;

    // allocate result
    vector<int> P(n+1);

    // get workspace 
    vector<int> W(8*(n+1));

    // add elbow room to C
    t = cnz + cnz/5 + 2*n ;

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
    for(k=0; k<n; ++k)
      len[k] = Cp[k+1]-Cp[k];
  
    len[n] = 0;
    nzmax = C.size();
    Ci = &C.colRef().front() ;
    for(i=0; i<=n; ++i){
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
    mark = wclear (0, 0, w, n);
  
    // n is a dead element
    elen[n] = -2;
  
    // n is a root of assembly tree
    Cp[n] = -1;
  
    // n is a dead element
    w[n] = 0;
  
    // --- Initialize degree lists ------------------------------------------
    for(i = 0; i<n; ++i){
      d = degree[i];
    
      // node i is empty
      if(d == 0){
        // element i is dead
        elen [i] = -2;
        nel++;
      
        // i is a root of assembly tree 
        Cp[i] = -1;
        w[i] = 0;
      } else if (d > dense) { // node i is dense
        // absorb i into element n
        nv[i] = 0;
      
        // node i is dead
        elen [i] = -1;
        nel++;
        Cp[i] = CS_FLIP (n) ;
        nv[n]++ ;
      } else {
        if(head[d] != -1)
          last[head[d]] = i;
      
        // put node i in degree list d
        next[i] = head[d];
        head[d] = i;
      }
    }
  
    // while (selecting pivots) do
    while(nel < n){
      // --- Select node of minimum approximate degree -------------------- 
      for(k = -1 ; mindeg < n && (k = head [mindeg]) == -1; mindeg++);
    
      if(next [k] != -1) last [next [k]] = -1 ;
    
      // remove k from degree list
      head[mindeg] = next[k];
    
      // elenk = |Ek|
      elenk = elen[k];
    
      // # of nodes k represents
      nvk = nv[k];
    
      // nv[k] nodes of A eliminated 
      nel += nvk;
    
      // --- Garbage collection -------------------------------------------
      if(elenk > 0 && cnz + mindeg >= nzmax){
        for(j=0; j<n; ++j){
          // j is a live node or element
          if((p = Cp [j]) >= 0){
            // save first entry of object
            Cp[j] = Ci[p];
          
            // first entry is now CS_FLIP(j)
            Ci[p] = CS_FLIP(j);
          }
        }
      
        // scan all of memory
        for (q = 0, p = 0 ; p < cnz ; ){
          // found object j
          if ((j = CS_FLIP (Ci [p++])) >= 0){
            // restore first entry of object 
            Ci [q] = Cp[j];
          
            // new pointer to object j
            Cp [j] = q++;
            for (k3 = 0 ; k3 < len [j]-1 ; k3++)
              Ci [q++] = Ci [p++] ;
          }
        }
      
        // Ci [cnz...nzmax-1] now free
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
      for (k1 = 1 ; k1 <= elenk + 1 ; k1++){
        if (k1 > elenk){
          // search the nodes in k
          e = k;
        
          // list of nodes starts at Ci[pj]
          pj = p;
        
          // length of list of nodes in k
          ln = len [k] - elenk;
        } else {
          // search the nodes in e
          e = Ci [p++];
          pj = Cp [e] ;
        
          // length of list of nodes in e
          ln = len [e];
        }
      
        for (k2=1; k2<=ln ; ++k2){
          i = Ci [pj++] ;
        
          // node i dead, or seen
          if ((nvi = nv [i]) <= 0) continue;
        
          // degree[Lk] += size of node i
          dk += nvi;
        
          // negate nv[i] to denote i in Lk
          nv [i] = -nvi;
        
          // place i in Lk
          Ci [pk2++] = i;
        
          if(next[i] != -1)
            last[next[i]] = last[i];
        
          // remove i from degree list 
          if (last[i] != -1){
            next[last[i]] = next[i] ;
          } else {
            head[degree[i]] = next[i] ;
          }
        }
      
        if(e != k){
          // absorb e into k
          Cp [e] = CS_FLIP (k);
        
          // e is now a dead element
          w [e] = 0;
        }
      }
    
      // Ci [cnz...nzmax] is free
      if(elenk != 0)
        cnz = pk2;
    
      // external degree of k - |Lk\i|
      degree [k] = dk;
    
      // element k is in Ci[pk1..pk2-1]
      Cp[k] = pk1;
      len [k] = pk2 - pk1 ;
    
      // k is now an element
      elen [k] = -2;
    
      // --- Find set differences ----------------------------------------- 
    
      // clear w if necessary
      mark = wclear(mark, lemax, w, n);
    
      // scan 1: find |Le\Lk|
      for(pk = pk1 ; pk < pk2 ; ++pk){
        i = Ci[pk] ;
      
        // skip if elen[i] empty 
        if ((eln = elen [i]) <= 0)
          continue;
      
        // nv [i] was negated
        nvi = -nv [i];
      
        wnvi = mark - nvi ;
      
        // scan Ei
        for(p = Cp [i] ; p <= Cp [i] + eln - 1 ; ++p){
          e = Ci[p];
          if (w[e] >= mark){
            // decrement |Le\Lk|
            w [e] -= nvi;
          } else if(w [e] != 0){        /* ensure e is a live element */
            w [e] = degree [e] + wnvi ; /* 1st time e seen in scan 1 */
          }
        }
      }
    
      // --- Degree update ------------------------------------------------
      // scan2: degree update
      for(pk = pk1 ; pk < pk2 ; ++pk){
        // consider node i in Lk
        i = Ci [pk];
        p1 = Cp [i] ;
        p2 = p1 + elen [i] - 1 ;
        pn = p1 ;
      
        // scan Ei
        for (h = 0, d = 0, p = p1 ; p <= p2 ; p++){
          e = Ci [p] ;
        
          // e is an unabsorbed element
          if (w [e] != 0){
            // dext = |Le\Lk|
            dext = w [e] - mark;
            if (dext > 0){
              // sum up the set differences
              d += dext;
            
              // keep e in Ei 
              Ci[pn++] = e;
            
              // compute the hash of node i
              h += e;
            
            } else {
              // aggressive absorb. e->k
              Cp [e] = CS_FLIP (k);
            
              // e is a dead element
              w [e] = 0;
            }
          }
        }
      
        // elen[i] = |Ei|
        elen [i] = pn - p1 + 1;
        p3 = pn ;
        p4 = p1 + len [i] ;
      
        // prune edges in Ai
        for (p = p2 + 1 ; p < p4 ; p++){
          j = Ci [p] ;
        
          // node j dead or in Lk
          if ((nvj = nv [j]) <= 0)
            continue;
        
          // degree(i) += |j|
          d += nvj;
        
          // place j in node list of i
          Ci[pn++] = j;
        
          // compute hash for node i
          h += j;
        }
      
        // check for mass elimination
        if (d == 0){
          // absorb i into k 
          Cp[i] = CS_FLIP(k);
          nvi = -nv [i];
        
          // |Lk| -= |i|
          dk -= nvi;
        
          // |k| += nv[i]
          nvk += nvi;
          nel += nvi;
          nv[i] = 0 ;
        
          // node i is dead
          elen [i] = -1;
        } else {
          // update degree(i)
          degree [i] = std::min(degree [i], d);
        
          // move first node to end
          Ci[pn] = Ci[p3];
        
          // move 1st el. to end of Ei
          Ci[p3] = Ci [p1];
        
          // add k as 1st element in of Ei
          Ci[p1] = k;
        
          // new len of adj. list of node i
          len[i] = pn - p1 + 1;
        
          // finalize hash of i
          h %= n;
        
          // place i in hash bucket
          next[i] = hhead [h];
          hhead [h] = i ;
        
          // save hash of i in last[i]
          last[i] = h;
        }
      } // scan2 is done
    
      // finalize |Lk|
      degree [k] = dk;
      lemax = std::max(lemax, dk);
    
      // clear w
      mark = wclear(mark+lemax, lemax, w, n);
    
      // --- Supernode detection ------------------------------------------
      for(pk = pk1 ; pk < pk2 ; pk++){
        i = Ci[pk] ;
      
        // skip if i is dead 
        if (nv [i] >= 0)
          continue;
      
        // scan hash bucket of node i
        h = last [i];
        i = hhead [h];
      
        // hash bucket will be empty
        hhead [h] = -1;
        for ( ; i != -1 && next [i] != -1 ; i = next [i], mark++){
          ln = len [i] ;
          eln = elen [i] ;
          for (p = Cp [i]+1 ; p <= Cp [i] + ln-1 ; p++)
            w [Ci [p]] = mark;
        
          jlast = i;
        
          // compare i with all j
          for (j = next [i] ; j != -1 ; ){
            ok = (len [j] == ln) && (elen [j] == eln) ;
            for (p = Cp [j] + 1 ; ok && p <= Cp [j] + ln - 1 ; p++){
              if (w [Ci [p]] != mark) ok = 0 ;    /* compare i and j*/
            }
          
            // i and j are identical
            if (ok){
              // absorb j into i
              Cp [j] = CS_FLIP (i);
              nv [i] += nv [j] ;
              nv [j] = 0;
            
              // node j is dead
              elen [j] = -1;
            
              // delete j from hash bucket
              j = next [j];
              next [jlast] = j ;
            } else {
              // j and i are different
              jlast = j;
              j = next [j] ;
            }
          }
        }
      }
    
      //  --- Finalize new element------------------------------------------
      // finalize Lk
      for (p = pk1, pk = pk1 ; pk < pk2 ; pk++){
        i = Ci [pk] ;
      
        // skip if i is dead
        if ((nvi = -nv [i]) <= 0)
          continue;
      
        // restore nv[i]
        nv [i] = nvi;
      
        // compute external degree(i)
        d = degree [i] + dk - nvi ;
        d = std::min(d, n - nel - nvi);
        if(head [d] != -1)
          last[head[d]] = i;
      
        // put i back in degree list
        next [i] = head [d];
        last [i] = -1 ;
        head [d] = i ;
      
        // find new minimum degree 
        mindeg = std::min(mindeg, d);
        degree [i] = d ;
      
        // place i in Lk
        Ci[p++] = i;
      }
    
      // # nodes absorbed into k
      nv [k] = nvk;
    
      // length of adj list of element k
      if ((len [k] = p-pk1) == 0){
        // k is a root of the tree
        Cp [k] = -1;
      
        // k is now a dead element
        w[k] = 0;
      }
    
      // free unused space in Lk
      if (elenk != 0)
        cnz = p;
    }
  
    // --- Postordering -----------------------------------------------------
  
    // fix assembly tree
    for(i=0; i<n; ++i)
      Cp[i] = CS_FLIP(Cp[i]);
  
    for (j = 0 ; j <= n ; j++)
      head [j] = -1 ;
  
    // place unordered nodes in lists
    for (j = n ; j >= 0 ; j--){
      // skip if j is an element
      if (nv [j] > 0) continue;
    
      // place j in list of its parent
      next[j] = head[Cp [j]];
      head[Cp[j]] = j;
    }
  
    // place elements in lists
    for (e = n ; e >= 0 ; e--){
      // skip unless e is an element
      if (nv[e] <= 0) continue;
      if (Cp[e] != -1){
        // place e in list of its parent
        next[e] = head[Cp [e]];
        head[Cp[e]] = e ;
      }
    }
  
    // postorder the assembly tree
    for(k = 0, i = 0 ; i <= n ; i++){
      if (Cp[i] == -1)
        k = depthFirstSearchAndPostorder(i, k, head, next, &P.front(), w) ;
    }
  
    return P;
  }

  int CRSSparsityInternal::scatter(int j, std::vector<int>& w, int mark, CRSSparsity& C, int nz) const{
    int i, p;
    const int *Ap = &rowind_.front();
    const int *Ai = &col_.front();
    int *Ci = &C.colRef().front();
  
    for(p = Ap[j]; p<Ap[j+1]; ++p){
      // A(i,j) is nonzero
      i = Ai [p];
    
      if (w[i] < mark){
        // i is new entry in column j
        w[i] = mark;
      
        // add i to pattern of C(:,j)
        Ci[nz++] = i;
      }
    }
    return nz;
  }

  CRSSparsity CRSSparsityInternal::multiply(const CRSSparsity& B) const{
    int nz = 0;
    casadi_assert_message(nrow_ == B.size2(), "Dimension mismatch.");
    int m = ncol_;
    int anz = rowind_[nrow_];
    int n = B.size1();
    const int* Bp = &B.rowind().front();
    const int* Bi = &B.col().front();
    int bnz = Bp[n];

    // get workspace
    vector<int> w(m);

    // allocate result
    CRSSparsity C(n,m);
    C.rowindRef().resize(anz + bnz);
  
    int* Cp = &C.rowindRef().front();
    for(int j=0; j<n; ++j){
      if(nz+m > C.size()){
        C.colRef().resize(2*(C.size())+m);
      }
    
      // column j of C starts here
      Cp[j] = nz;
      for(int p = Bp[j] ; p<Bp[j+1] ; ++p){
        nz = scatter(Bi[p], w, j+1, C, nz);
      }
    }
  
    // finalize the last column of C
    Cp[n] = nz;
    C.colRef().resize(nz);

    // Success
    return C;
  }

  void CRSSparsityInternal::prefactorize(int order, int qr, std::vector<int>& S_pinv, std::vector<int>& S_q, std::vector<int>& S_parent, std::vector<int>& S_cp, std::vector<int>& S_leftmost, int& S_m2, double& S_lnz, double& S_unz) const{
    int k;
    int n = nrow_;
    vector<int> post;
  
    // fill-reducing ordering
    if(order!=0){
      S_q = approximateMinimumDegree(order);
    }

    // QR symbolic analysis
    if (qr){
      CRSSparsity C;
      if(order!=0){
        std::vector<int> pinv_tmp;
        C = permute(pinv_tmp, S_q, 0);
      } else {
        C = shared_from_this<CRSSparsity>();
      }
    
      // etree of C'*C, where C=A(:,q)
      S_parent = C->eliminationTree(1);
    
      post = postorder(S_parent, n);

      // col counts chol(C'*C)
      S_cp = C->counts(&S_parent.front(), &post.front(), 1);
      post.clear();
    
      C->vcount(S_pinv, S_parent, S_leftmost, S_m2, S_lnz);
      for(S_unz = 0, k = 0; k<n; k++)
        S_unz += S_cp[k];
      
      // int overflow guard
      casadi_assert(S_lnz >= 0);
      casadi_assert(S_unz >= 0);
    } else {
      // for LU factorization only
      S_unz = 4*(rowind_[n]) + n ;
    
      // guess nnz(L) and nnz(U)
      S_lnz = S_unz;
    }
  }

  CRSSparsity CRSSparsityInternal::diag(std::vector<int>& mapping) const{
    if (nrow_==ncol_) {
      // Return object
      CRSSparsity ret(0,1);
      ret.reserve(std::min(size(),nrow_),nrow_);
    
      // Mapping
      mapping.clear();
    
      // Loop over nonzero
      for(int i=0; i<nrow_; ++i){
      
        // Enlarge the return matrix
        ret.resize(i+1,1);
    
        // Get to the right nonzero of the row
        int el = rowind_[i];
        while(el<rowind_[i+1] && col_[el]<i){
          el++;
        }
      
        if (el>=size()) return ret;
      
        // Add element if nonzero on diagonal
        if(col_[el]==i){
          ret.getNZ(i,0);
          mapping.push_back(el);
        }
      }
    
      return ret;
    
    } else if (nrow_==1 || ncol_==1) {
      CRSSparsity trans;
      const CRSSparsityInternal *sp;
    
      // Have a row vector
      if(nrow_ == 1){
        sp = this;
      } else {
        trans = transpose();
        sp = static_cast<const CRSSparsityInternal *>(trans.get());
      }
    
      // Return object
      mapping.clear();
      mapping.resize(size());
    
      std::vector<int> rowind(sp->ncol_+1,0);
      std::vector<int> col(sp->size());
    
      int i_prev = 0;
        
      // Loop over nonzero
      for(int k=0;k<size();k++) {
        mapping[k]=k; // mapping will just be a range(size())
     
        int i = sp->col_[k];
        std::fill(rowind.begin()+i_prev+1,rowind.begin()+i+1,k);
        col[k]=i;
        i_prev = i;
      }
      std::fill(rowind.begin()+i_prev+1,rowind.end(),size());
    
      return CRSSparsity(sp->ncol_,sp->ncol_,col,rowind);
    } else {
      casadi_error("diag: wrong argument shape. Expecting square matrix or vector-like, but got " << dimString() << " instead.");
    }
  }

  std::string CRSSparsityInternal::dimString() const { 
    std::stringstream ss;
    if (numel()==size()) {
      ss << nrow_ << "-by-" << ncol_ << " (dense)";
    } else {
      ss << nrow_ << "-by-" << ncol_ << " (" << size() << "/" << numel() << " nz)";
    }
    return ss.str();
  }

  CRSSparsity CRSSparsityInternal::patternProduct(const CRSSparsity& y_trans) const{
    // Dimensions
    int x_nrow = nrow_;
    int y_ncol = y_trans.size1();

    // Quick return if both are dense
    if(dense() && y_trans.dense()){
      return CRSSparsity(x_nrow,y_ncol,!empty() && !y_trans.empty());
    }
  
    // return object
    CRSSparsity ret(x_nrow,y_ncol);
  
    // Get the vectors for the return pattern
    vector<int>& c = ret.colRef();
    vector<int>& r = ret.rowindRef();
  
    // Direct access to the arrays
    const vector<int> &x_col = col_;
    const vector<int> &y_row = y_trans.col();
    const vector<int> &x_rowind = rowind_;
    const vector<int> &y_colind = y_trans.rowind();

    // If the compiler supports C99, we shall use the long long datatype, which is 64 bit, otherwise long
#if __STDC_VERSION__ >= 199901L
    typedef unsigned long long int_t;
#else
    typedef unsigned long int_t;
#endif
  
    // Number of directions we can deal with at a time
    int nr = CHAR_BIT*sizeof(int_t); // the size of int_t in bits (CHAR_BIT is the number of bits per byte, usually 8)

    // Number of such groups needed
    int ng = x_nrow/nr;
    if(ng*nr != x_nrow) ng++;
  
    // Which columns exist in a row of the first factor
    vector<int_t> in_x_row(ncol_);
    vector<int_t> in_res_col(y_ncol);

    // Loop over the rows of the resulting matrix, nr rows at a time
    for(int rr=0; rr<ng; ++rr){

      // Mark the elements in the x row
      fill(in_x_row.begin(),in_x_row.end(),0); // NOTE: expensive?
      int_t b=1;
      for(int i=rr*nr; i<rr*nr+nr && i<x_nrow; ++i){
        for(int el1=x_rowind[i]; el1<x_rowind[i+1]; ++el1){
          in_x_row[x_col[el1]] |= b;
        }
        b <<= 1;
      }

      // Get the sparsity pattern for the set of rows
      fill(in_res_col.begin(),in_res_col.end(),0); // NOTE: expensive?
      for(int j=0; j<y_ncol; ++j){
      
        // Loop over the nonzeros of the column of the second factor
        for(int el2=y_colind[j]; el2<y_colind[j+1]; ++el2){
        
          // Get the row
          int i_y = y_row[el2];
        
          // Add nonzero if the element matches an element in the x row
          in_res_col[j] |= in_x_row[i_y];
        }
      }

      b = 1;
      for(int i=rr*nr; i<rr*nr+nr && i<x_nrow; ++i){
      
        // loop over the columns of the resulting matrix
        for(int j=0; j<y_ncol; ++j){
        
          // Save nonzero, if any
          if(in_res_col[j] & b){
            c.push_back(j);
          }
        }
        r[i+1] = c.size();
        b <<=1;
      }
    }
    return ret;
  }

  CRSSparsity CRSSparsityInternal::patternProduct(const CRSSparsity& y_trans, vector< vector< pair<int,int> > >& mapping) const{
    // return object
    CRSSparsity ret = patternProduct(y_trans);
  
    // Get the vectors for the return pattern
    const vector<int>& c = ret.col();
    const vector<int>& r = ret.rowind();
  
    // Direct access to the arrays
    const vector<int> &x_col = col_;
    const vector<int> &y_row = y_trans.col();
    const vector<int> &x_rowind = rowind_;
    const vector<int> &y_colind = y_trans.rowind();

    // Clear the mapping
    mapping.resize(ret.size());

    // the entry of the matrix to be calculated
    vector< pair<int,int> > d;

    // loop over the row of the resulting matrix)
    for(int i=0; i<nrow_; ++i){
      // Loop over nonzeros
      for(int el=r[i]; el<r[i+1]; ++el){
        int j = c[el];
        int el1 = x_rowind[i];
        int el2 = y_colind[j];
        d.clear();
        while(el1 < x_rowind[i+1] && el2 < y_colind[j+1]){ // loop over non-zero elements
          int j1 = x_col[el1];
          int i2 = y_row[el2];      
          if(j1==i2){
            d.push_back(pair<int,int>(el1++,el2++));
          } else if(j1<i2) {
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

  bool CRSSparsityInternal::scalar(bool scalar_and_dense) const{
    return nrow_==1 && ncol_==1 && (!scalar_and_dense || size()==1);
  }

  bool CRSSparsityInternal::dense() const{
    return size() == numel();
  }
  
  bool CRSSparsityInternal::empty() const{
    return numel()==0;
  }
  
  bool CRSSparsityInternal::null() const{
    return nrow_==0 && ncol_==0;
  }

  bool CRSSparsityInternal::diagonal() const{
    // Check if matrix is square
    if(nrow_ != ncol_) return false;
    
    // Check if correct number of non-zeros (one per row)
    if(size() != nrow_) return false;

    // Check that the column indices are correct
    for(int i=0; i<size(); ++i){
      if(col_[i]!=i)
        return false;
    }
   
    // Make sure that the row indices are correct
    for(int i=0; i<nrow_; ++i){
      if(rowind_[i]!=i)
        return false;
    }
  
    // Diagonal if reached this point
    return true;
  }

  bool CRSSparsityInternal::square() const{
    return nrow_ == ncol_;
  }

  int CRSSparsityInternal::sizeU() const{
    int nnz = 0;
    for(int r=0; r<nrow_; ++r){
      for(int el = rowind_[r]; el < rowind_[r+1]; ++el){
        nnz += col_[el]>=r;
      }
    }
    return nnz;
  }
  
  int CRSSparsityInternal::sizeD() const{
    int nnz = 0;
    for(int r=0; r<nrow_; ++r){
      for(int el = rowind_[r]; el < rowind_[r+1]; ++el){
        nnz += col_[el]==r;
      }
    }
    return nnz;
  }

  int CRSSparsityInternal::sizeL() const{
    int nnz = 0;
    for(int r=0; r<nrow_; ++r){
      for(int el = rowind_[r]; el < rowind_[r+1] && col_[el]<=r; ++el){
        nnz ++;
      }
    }
    return nnz;
  }

  std::pair<int,int> CRSSparsityInternal::shape() const{
    return std::pair<int,int>(nrow_,ncol_);
  }

  vector<int> CRSSparsityInternal::erase(const vector<int>& ii, const vector<int>& jj){
    if (!inBounds(ii,nrow_)) {
      casadi_error("Slicing [ii,jj] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(jj,ncol_)) {
      casadi_error("Slicing [ii,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    // Mapping
    vector<int> mapping;
  
    // Quick return if no elements
    if(numel()==0)
      return mapping;
  
    // Reserve memory
    mapping.reserve(size());
  
    // Number of non-zeros
    int nz=0;
  
    // Rows to be erased
    vector<int>::const_iterator ie = ii.begin();
  
    // First and last index for the row
    int el_first=0, el_last=0;
  
    // Loop over rows
    for(int i=0; i<nrow_; ++i){
      // Update beginning and end of non-zero indices
      el_first = el_last;
      el_last = rowind_[i+1];
    
      // Is it a row that can be deleted
      bool deletable_row = ie!=ii.end() && *ie==i;
      if(deletable_row){
        ie++;
      
        // Columns to be erased
        vector<int>::const_iterator je = jj.begin();

        // Loop over nonzero elements of the row
        for(int el=el_first; el<el_last; ++el){
          // Column
          int j=col_[el];
        
          // Continue to the next column to skip
          for(; je!=jj.end() && *je<j; ++je);
        
          // Remove column if necessary
          if(je!=jj.end() && *je==j){
            je++;
            continue;
          }
        
          // Save old nonzero for each new nonzero
          mapping.push_back(el);
        
          // Update column and increase nonzero counter
          col_[nz++] = j;
        }
      } else {
        // Loop over nonzero elements of the row
        for(int el=el_first; el<el_last; ++el){
          // Column
          int j=col_[el];
      
          // Save old nonzero for each new nonzero
          mapping.push_back(el);
      
          // Update column and increase nonzero counter
          col_[nz++] = j;
        }
      }
    
      // Register last nonzero of the row
      rowind_[i+1]=nz;
    }
  
    // Truncate column matrix
    col_.resize(nz);
  
    return mapping;
  }

  vector<int> CRSSparsityInternal::getNZ(const vector<int>& ii, const vector<int>& jj) const{
    if (!inBounds(ii,nrow_)) {
      casadi_error("Slicing [ii,jj] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(jj,ncol_)) {
      casadi_error("Slicing [ii,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    std::vector<int> jj_sorted;
    std::vector<int> jj_sorted_index;
  
    sort(jj, jj_sorted, jj_sorted_index);
  
    vector<int> ret(ii.size()*jj.size());
  
    int stride = jj.size();
  
    for(int i=0;i<ii.size();++i){
      int it = ii[i];
      int el=rowind_[it];
      for(int j=0;j<jj_sorted.size();++j){
        int jt=jj_sorted[j];
        // Continue to the non-zero element
        for(; el<rowind_[it+1] && col_[el]<jt; ++el){}
        // Add the non-zero element, if there was an element in the location exists
        if(el<rowind_[it+1] && col_[el]== jt) {
          ret[i*stride+jj_sorted_index[j]] = el;
        }
        else
          ret[i*stride+jj_sorted_index[j]] = -1;
      }
    }
    return ret;
  }

  CRSSparsity CRSSparsityInternal::sub(const vector<int>& ii, const vector<int>& jj, vector<int>& mapping) const{
    if (!inBounds(ii,nrow_)) {
      casadi_error("Slicing [ii,jj] out of bounds. Your ii contains " << *std::min_element(ii.begin(),ii.end()) << " up to " << *std::max_element(ii.begin(),ii.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
    if (!inBounds(jj,ncol_)) {
      casadi_error("Slicing [ii,jj] out of bounds. Your jj contains " << *std::min_element(jj.begin(),jj.end()) << " up to " << *std::max_element(jj.begin(),jj.end()) << ", which is outside of the matrix shape " << dimString() << ".");
    }
  
    if (double(ii.size())*double(jj.size()) > size()) {
      // Typical use case:
      // a = ssym("a",sp_diag(50000))
      // a[:,:]
      return sub2(ii,jj,mapping);
    } else {
      // Typical use case:
      // a = DMatrix.ones(1000,1000)
      // a[[0,1],[0,1]]
      return sub1(ii,jj,mapping);
    }

  }

  CRSSparsity CRSSparsityInternal::sub2(const vector<int>& ii, const vector<int>& jj, vector<int>& mapping) const{
    std::vector<int> jj_sorted;
    std::vector<int> jj_sorted_index;

    std::vector<int> ii_sorted;
    std::vector<int> ii_sorted_index;

    sort(jj, jj_sorted, jj_sorted_index, false);
    sort(ii, ii_sorted, ii_sorted_index, false);

    std::vector<int> jjlookup = lookupvector(jj_sorted,ncol_);
  
    // count the number of non-zeros
    int nnz = 0;

    // loop over the rows of the slice
    for(int i=0;i<ii.size();++i){
      int it = ii_sorted[i];
      for(int el=rowind_[it]; el<rowind_[it+1]; ++el){ // loop over the non-zeros of the matrix
        int j = col_[el];
        int ji = jjlookup[j];
        if (ji!=-1) {
          int jv = jj_sorted[ji];
          while (ji>=0 && jv == jj_sorted[ji--]) nnz++; 
        }
      }
    }
  
    mapping.resize(nnz);
  
    vector<int> rows(nnz);
    vector<int> cols(nnz);

    int k = 0;
    // loop over the row of the slice
    for(int i=0;i<ii.size();++i){
      int it = ii_sorted[i];
      for(int el=rowind_[it]; el<rowind_[it+1]; ++el){ // loop over the non-zeros of the matrix
        int jt = col_[el];
        int ji = jjlookup[jt];
        if (ji!=-1) {
          int jv = jj_sorted[ji];
          while (ji>=0 && jv == jj_sorted[ji]) {
            cols[k] = jj_sorted_index[ji];
            rows[k] = ii_sorted_index[i];
            mapping[k] = el;
            k++;
            ji--;
          }
        }
      }
    }
  
    std::vector<int> sp_mapping;
    std::vector<int> mapping_ = mapping;
    CRSSparsity ret = sp_triplet(ii.size(),jj.size(),rows,cols,sp_mapping);
  
    for (int i=0;i<mapping.size();++i)
      mapping[i] = mapping_[sp_mapping[i]];
  
    // Create sparsity pattern
    return ret;
  }

  CRSSparsity CRSSparsityInternal::sub1(const vector<int>& ii, const vector<int>& jj, vector<int>& mapping) const{

    std::vector<int> jj_sorted;
    std::vector<int> jj_sorted_index;

    std::vector<int> ii_sorted;
    std::vector<int> ii_sorted_index;

    sort(jj, jj_sorted, jj_sorted_index, false);
    sort(ii, ii_sorted, ii_sorted_index, false);
  
    // count the number of non-zeros
    int nnz = 0;
  
    for(int i=0;i<ii.size();++i){
      int it = ii_sorted[i];
      int el=rowind_[it];
      for(int j=0;j<jj_sorted.size();++j){
        int jt=jj_sorted[j];
        // Continue to the non-zero element
        for(; el<rowind_[it+1] && col_[el]<jt; ++el){}
        // Add the non-zero element, if there was an element in the location exists
        if(el<rowind_[it+1] && col_[el]== jt) {
          nnz++;
        }
      }
    }
  
    mapping.resize(nnz);
  
    vector<int> rows(nnz);
    vector<int> cols(nnz);
  

    int k=0;
    for(int i=0;i<ii.size();++i){
      int it = ii_sorted[i];
      int K = rowind_[it];
      for(int j=0;j<jj_sorted.size();++j){
        int jt=jj_sorted[j];
        // Continue to the non-zero element
        for(; K<rowind_[it+1] && col_[K]<jt; ++K){}
        // Add the non-zero element, if there was an element in the location exists
        if(K<rowind_[it+1] && col_[K]== jt) {
          cols[k] = jj_sorted_index[j];
          rows[k] = ii_sorted_index[i];
          mapping[k] = K;
          k++;
        }
      }
    }
  
    std::vector<int> sp_mapping;
    std::vector<int> mapping_ = mapping;
    CRSSparsity ret = sp_triplet(ii.size(),jj.size(),rows,cols,sp_mapping);
  
    for (int i=0;i<mapping.size();++i)
      mapping[i] = mapping_[sp_mapping[i]];
  
    // Create sparsity pattern
    return ret;
  }

  CRSSparsity CRSSparsityInternal::patternCombine(const CRSSparsity& y, bool f0x_is_zero, bool fx0_is_zero) const{
    static vector<unsigned char> mapping;
    return patternCombineGen1<false>(y, f0x_is_zero, fx0_is_zero, mapping);
  }

  CRSSparsity CRSSparsityInternal::patternCombine(const CRSSparsity& y, bool f0x_is_zero, bool fx0_is_zero, vector<unsigned char>& mapping) const{
    return patternCombineGen1<true>(y, f0x_is_zero, fx0_is_zero, mapping);    
  }
  
  template<bool with_mapping>
  CRSSparsity CRSSparsityInternal::patternCombineGen1(const CRSSparsity& y, bool f0x_is_zero, bool fx0_is_zero, std::vector<unsigned char>& mapping) const{

    // Quick return if identical
    if(isEqual(y)){
      if(with_mapping){
        mapping.resize(y.size());
        fill(mapping.begin(),mapping.end(), 1 | 2);
      }
      return y;
    }

    if(f0x_is_zero){
      if(fx0_is_zero){
        return patternCombineGen<with_mapping,true,true>(y,mapping);
      } else {
        return patternCombineGen<with_mapping,true,false>(y,mapping);
      }      
    } else if(fx0_is_zero){
      return patternCombineGen<with_mapping,false,true>(y,mapping);
    } else {
      return patternCombineGen<with_mapping,false,false>(y,mapping);
    }
  }
  
  template<bool with_mapping, bool f0x_is_zero, bool fx0_is_zero>
  CRSSparsity CRSSparsityInternal::patternCombineGen(const CRSSparsity& y, vector<unsigned char>& mapping) const{

    // Assert dimensions
    casadi_assert_message(nrow_==y.size1() && ncol_==y.size2(), "Dimension mismatch");
    
    // Sparsity pattern of the argument
    const vector<int>& y_rowind = y.rowind();
    const vector<int>& y_col = y.col();
    
    // Sparsity pattern of the result
    vector<int> ret_rowind(nrow_+1,0);
    vector<int> ret_col;
    
    // Clear the mapping
    if(with_mapping) mapping.clear();

    // Loop over rows of both patterns
    for(int i=0; i<nrow_; ++i){
      // Non-zero element of the two matrices
      int el1 = rowind_[i];
      int el2 = y_rowind[i];
      
      // End of the non-zero elements of the row for the two matrices
      int el1_last = rowind_[i+1];
      int el2_last = y_rowind[i+1];
      
      // Loop over the non-zeros of both matrices
      while(el1<el1_last || el2<el2_last){
        // Get the columns
        int col1 = el1<el1_last ? col_[el1] : ncol_;
        int col2 = el2<el2_last ? y_col[el2] : ncol_;

        // Add to the return matrix
        if(col1==col2){ //  both nonzero
          ret_col.push_back(col1);
          if(with_mapping) mapping.push_back( 1 | 2);
          el1++; el2++;
        } else if(col1<col2){ //  only first argument is nonzero
          if(!fx0_is_zero){
            ret_col.push_back(col1);
            if(with_mapping) mapping.push_back(1);
          } else {
            if(with_mapping) mapping.push_back(1 | 4);
          }
          el1++;
        } else { //  only second argument is nonzero
          if(!f0x_is_zero){
            ret_col.push_back(col2);
            if(with_mapping) mapping.push_back(2);
          } else {
            if(with_mapping) mapping.push_back(2 | 4);
          }
          el2++;
        }
      }
      
      // Save the index of the last nonzero on the row
      ret_rowind[i+1] = ret_col.size();
    }
    
    // Return cached object
    return CRSSparsity(nrow_, ncol_, ret_col, ret_rowind);
  }

  bool CRSSparsityInternal::isEqual(const CRSSparsity& y) const{
    // Quick true if the objects are the same
    if(this == y.get()) return true;  
  
    // Otherwise, compare the patterns
    return isEqual(y.size1(),y.size2(),y.col(),y.rowind());
  }
  
  CRSSparsity CRSSparsityInternal::patternInverse() const {
    // Quick return clauses
    if (empty()) return CRSSparsity(nrow_,ncol_,true);
    if (dense()) return CRSSparsity(nrow_,ncol_,false);
    
    // Sparsity of the result
    std::vector<int> col_ret;
    std::vector<int> rowind_ret=rowind_;
    
    // Loop over rows
    for (int i=0;i<nrow_;++i) {
      // Update rowind vector of the result
      rowind_ret[i+1]=rowind_ret[i]+ncol_-(rowind_[i+1]-rowind_[i]);
      
      // Counter of new column indices
      int j=0;
      
      // Loop over all nonzeros
      for (int k=rowind_[i];k<rowind_[i+1];++k) {
      
        // Try to reach current nonzero
        while(j<col_[k])  {
          // And meanwhile, add nonzeros to the result
          col_ret.push_back(j);
          j++;
        }
        j++;
      } 
      // Process the remainder up to the column size
      while(j < ncol_)  {
        col_ret.push_back(j);
        j++;
      }
    }
    
    // Return result
    return CRSSparsity(nrow_,ncol_,col_ret,rowind_ret);    
  }


  bool CRSSparsityInternal::isEqual(int nrow, int ncol, const std::vector<int>& col, const std::vector<int>& rowind) const{
    // First check dimensions and number of non-zeros
    if(size()!=col.size() || nrow_!=nrow || ncol_!=ncol)
      return false;

    // Check if dense
    if(size()==numel())
      return true;
  
    // Check the number of non-zeros per row
    if(!equal(rowind_.begin(),rowind_.end(),rowind.begin()))
      return false;
  
    // Finally check the column indices
    if(!equal(col_.begin(),col_.end(),col.begin()))
      return false;
  
    // Equal if reached this point
    return true;
  }

  void CRSSparsityInternal::reserve(int nnz, int nrow){
    col_.reserve(nnz);
    rowind_.reserve(nrow+1);
  }

  void CRSSparsityInternal::append(const CRSSparsity& sp){
    // Assert dimensions
    casadi_assert_message(ncol_==sp.size2(),"CRSSparsityInternal::append: Dimension mismatch. You attempt to append a shape " << sp.dimString() << " to a shape " << dimString() << ". The number of columns must match.");
  
    // Get current number of non-zeros
    int sz = size();
  
    // Add column indices
    col_.insert(col_.end(),sp.col().begin(),sp.col().end());
  
    // Add row indices
    rowind_.pop_back();
    rowind_.insert(rowind_.end(),sp.rowind().begin(),sp.rowind().end());
    for(int i = nrow_; i<rowind_.size(); ++i)
      rowind_[i] += sz;
  
    // Update dimensions
    nrow_ += sp.size1();
  }

  void CRSSparsityInternal::enlargeRows(int nrow, const std::vector<int>& ii){
    // Assert dimensions
    casadi_assert(ii.size() == nrow_);

    // Update dimensions
    nrow_ = nrow;

    // Sparsify the rows
    rowind_.resize(nrow+1,size());
  
    // Quick return if matrix had no rows before
    if(ii.empty()) return;
  
    int ik=ii.back(); // need only to update from the last new index
    int nz=size(); // number of nonzeros up till this row
    for(int i=ii.size()-1; i>=0; --i){
      // Update rowindex for new rows
      for(; ik>ii[i]; --ik){
        rowind_[ik] = nz;
      }
    
      // Update non-zero counter
      nz = rowind_[i];
    
      // Update rowindex for old rows
      rowind_[ii[i]] = nz;
    }
  
    // Append zeros to the beginning
    for(; ik>=0; --ik){
      rowind_[ik] = 0;
    }
  }

  void CRSSparsityInternal::enlargeColumns(int ncol, const std::vector<int>& jj){
    // Assert dimensions
    casadi_assert(jj.size() == ncol_);
  
    // Update dimensions
    ncol_ = ncol;

    // Begin by sparsify the columns
    for(int k=0; k<col_.size(); ++k){
      col_[k] = jj[col_[k]];
    }
  }

  CRSSparsity CRSSparsityInternal::makeDense(std::vector<int>& mapping) const{
    mapping.resize(size());
    for(int i=0; i<nrow_; ++i){
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){
        int j = col_[el];
        mapping[el] = j + i*ncol_;
      }
    }
  
    return CRSSparsity(nrow_,ncol_,true);
  }

  int CRSSparsityInternal::getNZ(int i, int j) const{
    casadi_assert_message(i<nrow_,"First index (" << i  << ") out of bounds. Attempting to slice [" << i << "," << j << " ] out of shape " << dimString() << ".");
    casadi_assert_message(j<ncol_,"Second index (" << j  << ") out of bounds.  Attempting to slice [" << i << "," << j << " ] out of shape " << dimString() << ".");
  
    if (i<0) i += nrow_;
    if (j<0) j += ncol_;
  
    // Quick return if matrix is dense
    if(numel()==size())
      return j+i*ncol_;
  
    // Quick return if past the end
    if(rowind_[i]==size() || (rowind_[i+1]==size() && col_.back()<j)){
      return -1;
    }

    // Find sparse element
    for(int ind=rowind_[i]; ind<rowind_[i+1]; ++ind){
      if(col_[ind] == j){
        return ind;     // element exists
      }
      else if(col_[ind] > j)
        break;                // break at the place where the element should be added
    }
    return -1;
  }

  CRSSparsity CRSSparsityInternal::reshape(int n, int m) const{
    casadi_assert_message(numel() == n*m, "reshape: number of elements must remain the same. Old shape is " << dimString() << ". New shape is " << n << "x" << m << "=" << n*m << ".");
    CRSSparsity ret(n,m);
    ret.reserve(size(), n);
  
    std::vector<int> row(size());
    std::vector<int> col(size());
    for(int i=0; i<nrow_; ++i){
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){
        int j = col_[el];
      
        // Element number
        int k_ret = j+i*ncol_;
      
        // Row and column in the new matrix
        int i_ret = k_ret/m;
        int j_ret = k_ret%m;
        row[el] = i_ret;
        col[el] = j_ret;
      }
    }
    return sp_triplet(n,m,row,col);
  }

  void CRSSparsityInternal::resize(int nrow, int ncol){
    if(nrow != nrow_ || ncol != ncol_){
      if(nrow < nrow_ || ncol < ncol_){
        // Row and column index of the new
        vector<int> col_new, rowind_new(nrow+1,0);

        // Loop over the rows which may contain nonzeros
        int i;
        for(i=0; i<nrow_ && i<nrow; ++i){
          // First nonzero element of the row
          rowind_new[i] = col_new.size();
        
          // Record columns of the nonzeros
          for(int el=rowind_[i]; el<rowind_[i+1] && col_[el]<ncol; ++el){
            col_new.push_back(col_[el]);
          }
        }
      
        // Save row-indices for the rest of the rows
        for(; i<nrow+1; ++i){
          rowind_new[i] = col_new.size();
        }
        
        // Save the sparsity
        nrow_ = nrow;
        ncol_ = ncol;
        col_.swap(col_new);
        rowind_.swap(rowind_new);
      
      } else {
        // Make larger: Very cheap operation
        nrow_ = nrow;
        ncol_ = ncol;
        rowind_.resize(nrow_+1,size());
      }
    }
  }

  bool CRSSparsityInternal::columnsSequential(bool strictly) const{
    for(int i=0; i<nrow_; ++i){
      int lastcol = -1;
      for(int k=rowind_[i]; k<rowind_[i+1]; ++k){
      
        // check if not in sequence
        if(col_[k] < lastcol)
          return false;

        // Check if duplicate
        if(strictly && col_[k] == lastcol)
          return false;

        // update last column of the row
        lastcol = col_[k]; 
      }
    }
  
    // sequential if reached this point
    return true;
  }

  void CRSSparsityInternal::removeDuplicates(std::vector<int>& mapping){
    casadi_assert(mapping.size()==size());
  
    // Nonzero counter without duplicates
    int k_strict=0;
  
    // Loop over rows
    for(int i=0; i<nrow_; ++i){
    
      // Last column encountered on the row so far
      int lastcol = -1;
    
      // Save new row offset (cannot set it yet, since we will need the old value below)
      int new_rowind = k_strict;
    
      // Loop over nonzeros (including duplicates)
      for(int k=rowind_[i]; k<rowind_[i+1]; ++k){
      
        // Make sure that the columns appear sequentially
        casadi_assert_message(col_[k] >= lastcol, "columns are not sequential");

        // Skip if duplicate
        if(col_[k] == lastcol)
          continue;

        // update last column encounterd on the row
        lastcol = col_[k]; 

        // Update mapping
        mapping[k_strict] = mapping[k];
      
        // Update column index
        col_[k_strict] = col_[k];
      
        // Increase the strict nonzero counter
        k_strict++;
      }
    
      // Update row offset
      rowind_[i] = new_rowind;
    }
  
    // Finalize the sparsity pattern
    rowind_[nrow_] = k_strict;
    col_.resize(k_strict);
    mapping.resize(k_strict);
  }

  void CRSSparsityInternal::getElements(std::vector<int>& loc, bool row_major) const{

    // Element for each nonzero
    loc.resize(size());
    
    // Loop over rows
    for(int i=0; i<nrow_; ++i){
      
      // Loop over the nonzeros
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){
        
        // Get column
        int j = col_[el];
        
        // Get the element
        if(row_major){
          loc[el] = j+i*ncol_;
        } else {
          loc[el] = i+j*nrow_;
        }
      }
    }
  }

  void CRSSparsityInternal::getNZInplace(std::vector<int>& indices) const{
    // Quick return if no elements
    if(indices.empty()) return;

    // Make a sanity check
    int last=-1;
    for(vector<int>::iterator it=indices.begin(); it!=indices.end(); ++it){
      if(*it>=0){
        int el_row = *it % nrow_;
        int el_col = *it / nrow_;
        int el = ncol_*el_row + el_col;
        casadi_assert_message(el>=last,"Elements must be sorted row-wise in non-decreasing order");
        last = el;
      }
    }

    // Quick return if no elements
    if(last<0) return;

    // Iterator to input/output
    vector<int>::iterator it=indices.begin();
    while(*it<0) it++; // first non-ignored

    // Current element sought
    int el_row = *it % nrow_;
    int el_col = *it / nrow_;
  
    // Loop over rows
    for(int i=0; i<nrow_; ++i){
    
      // Loop over the nonzeros
      for(int el=rowind_[i]; el<rowind_[i+1] && el_row<=i; ++el){
        
        // Get column
        int j = col_[el];
      
        // Add leading elements not in pattern
        while(i>el_row || (i==el_row && j>el_col)){
          // Mark as not found
          *it = -1;
        
          // Increase index and terminate if end of vector reached
          if(++it==indices.end()) return;
        
          // Next element sought
          el_row = *it % nrow_;
          el_col = *it / nrow_;
        }

        // Add elements in pattern
        while(i==el_row && j==el_col){
          // Save element index
          *it = el;

          // Increase index and terminate if end of vector reached
          do{
            if(++it==indices.end()) return;
          } while(*it<0);

          // Next element sought
          el_row = *it % nrow_;
          el_col = *it / nrow_;
        
        }
      }
    }
  
    // Add trailing elements not in pattern
    fill(it,indices.end(),-1);
  }

  CRSSparsity CRSSparsityInternal::unidirectionalColoring(const CRSSparsity& AT, int cutoff) const{
  
    // Allocate temporary vectors
    vector<int> forbiddenColors;
    forbiddenColors.reserve(nrow_);
    vector<int> color(nrow_,0);
  
    // Access the sparsity of the transpose
    const vector<int>& AT_rowind = AT.rowind();
    const vector<int>& AT_col = AT.col();
  
    // Loop over rows
    for(int i=0; i<nrow_; ++i){
    
      // Loop over nonzero elements
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){
      
        // Get column
        int c = col_[el];
        
        // Loop over previous rows that have an element in column c
        for(int el_prev=AT_rowind[c]; el_prev<AT_rowind[c+1]; ++el_prev){
        
          // Get the row
          int i_prev = AT_col[el_prev];
        
          // Escape loop if we have arrived at the current row
          if(i_prev>=i)
            break;
        
          // Get the color of the row
          int color_prev = color[i_prev];
        
          // Mark the color as forbidden for the current row
          forbiddenColors[color_prev] = i;
        }
      }
    
      // Get the first nonforbidden color
      int color_i;
      for(color_i=0; color_i<forbiddenColors.size(); ++color_i){
        // Break if color is ok
        if(forbiddenColors[color_i]!=i) break;
      }
      color[i] = color_i;
    
      // Add color if reached end
      if(color_i==forbiddenColors.size()){
        forbiddenColors.push_back(0);
        
        // Cutoff if too many colors
        if(forbiddenColors.size()>cutoff){
          return CRSSparsity();
        }
      }
    }
  
    // Create return sparsity containing the coloring
    CRSSparsity ret(forbiddenColors.size(),nrow_);
    vector<int>& rowind = ret.rowindRef();
    vector<int>& col = ret.colRef();
  
    // Get the number of columns for each row
    for(int i=0; i<color.size(); ++i){
      rowind[color[i]+1]++;
    }
  
    // Cumsum
    for(int j=0; j<forbiddenColors.size(); ++j){
      rowind[j+1] += rowind[j];
    }
  
    // Get column for each row
    col.resize(color.size());
    for(int j=0; j<col.size(); ++j){
      col[rowind[color[j]]++] = j;
    }
  
    // Swap index back one step
    for(int j=rowind.size()-2; j>=0; --j){
      rowind[j+1] = rowind[j];
    }
    rowind[0] = 0;
  
    // Return the coloring
    return ret;
  }
  
  CRSSparsity CRSSparsityInternal::starColoring2(int ordering, int cutoff) const{
    casadi_assert_warning(ncol_==nrow_,"StarColoring requires a square matrix, but got " << dimString() << ".");

    // TODO What we need here, is a distance-2 smallest last ordering
    // Reorder, if necessary
    if(ordering!=0){
      casadi_assert(ordering==1);
    
      // Ordering
      vector<int> ord = largestFirstOrdering();

      // Create a new sparsity pattern 
      CRSSparsity sp_permuted = pmult(ord,true,true,true);
    
      // Star coloring for the permuted matrix
      CRSSparsity ret_permuted = sp_permuted.starColoring2(0);
        
      // Permute result back
      return ret_permuted.pmult(ord,false,true,false);
    }
    
    // Allocate temporary vectors
    vector<int> forbiddenColors;
    forbiddenColors.reserve(nrow_);
    vector<int> color(nrow_,-1);
    
    vector<int> firstNeighborP(nrow_,-1);
    vector<int> firstNeighborQ(nrow_,-1);
    vector<int> firstNeighborQ_el(nrow_,-1);
    
    vector<int> treated(nrow_,-1);
    vector<int> hub(sizeL(),-1);

    vector<int> Tmapping;
    transpose(Tmapping);
    
    vector<int> star(size());
    int k = 0;
    for(int i=0; i<nrow_; ++i){ 
      for(int j_el=rowind_[i]; j_el<rowind_[i+1]; ++j_el){ 
        int j = col_[j_el];
        if (i<j) {
          star[j_el] = k;
          star[Tmapping[j]] = k;         
          k++;
        }
      }
    }
    

    
    int starID = 0;

    // 3: for each v \in V do
    for(int v=0; v<nrow_; ++v){ 
      
      // 4: for each colored w \in N1(v) do
      for(int w_el=rowind_[v]; w_el<rowind_[v+1]; ++w_el){ 
          int w = col_[w_el];
          int colorW = color[w];
          if(colorW==-1) continue;
          
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
                for(int x_el=rowind_[q]; x_el<rowind_[q+1]; ++x_el){
                  int x = col_[x_el];
                  if(color[x]==-1) continue;
                  
                  // treat@3: forbiddenColors[color[x]] <- v
                  forbiddenColors[color[x]] = v;
                }
                
                // treat@4: treated[q] <- v
                treated[q] = v;

            }
            // 10: treat(v, w) < forbid colors of neighbors of w
            
              // treat@2: for each colored x \in N1 (w) do
              for(int x_el=rowind_[w]; x_el<rowind_[w+1]; ++x_el){
                int x = col_[x_el];
                if(color[x]==-1) continue;
                
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
            int x_el_end = rowind_[w+1]; 
            int x, colorx;
            for(int x_el=rowind_[w]; x_el < x_el_end; ++x_el){
              x = col_[x_el];
              colorx = color[x];
              if(colorx==-1 || x==v) continue;
              
              // 14: if x = hub[star[wx]] then potential Case 2
              if (hub[star[x_el]]==x) {

                // 15: forbiddenColors[color[x]] <- v
                forbiddenColors[colorx] = v;
          
              }
            }
          }
          
      }
      
      // 16: color[v] <- min{c > 0 : forbiddenColors[c] != v}
      bool new_color = true;
      for(int color_i=0; color_i<forbiddenColors.size(); ++color_i){
        // Break if color is ok
        if(forbiddenColors[color_i]!=v){
          color[v] = color_i;
          new_color = false;
          break;
        }
      }
      
      // New color if reached end
      if(new_color){
        color[v] = forbiddenColors.size();
        forbiddenColors.push_back(-1);

        // Cutoff if too many colors
        if(forbiddenColors.size()>cutoff){
          return CRSSparsity();
        }
      }
      
      // 17: updateStars(v)
      
        // updateStars@2: for each colored w \in N1 (v) do
        for(int w_el=rowind_[v]; w_el<rowind_[v+1]; ++w_el){ 
            int w = col_[w_el];
            int colorW = color[w];
            if(colorW==-1) continue;
            
            // updateStars@3: if exits x \in N1 (w) where x = v and color[x] = color[v] then
            bool check = false;
            int x;
            int x_el;
            for(x_el=rowind_[w]; x_el<rowind_[w+1]; ++x_el){
              x = col_[x_el];
              if(x==v || color[x]!=color[v]) continue;
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
    CRSSparsity ret(forbiddenColors.size(),nrow_);
    vector<int>& rowind = ret.rowindRef();
    vector<int>& col = ret.colRef();
  
    // Get the number of columns for each row
    for(int i=0; i<color.size(); ++i){
      rowind[color[i]+1]++;
    }
  
    // Cumsum
    for(int j=0; j<forbiddenColors.size(); ++j){
      rowind[j+1] += rowind[j];
    }
  
    // Get column for each row
    col.resize(color.size());
    for(int j=0; j<col.size(); ++j){
      col[rowind[color[j]]++] = j;
    }
  
    // Swap index back one step
    for(int j=rowind.size()-2; j>=0; --j){
      rowind[j+1] = rowind[j];
    }
    rowind[0] = 0;
  
    // Return the coloring
    return ret;
    
  }
  

  CRSSparsity CRSSparsityInternal::starColoring(int ordering, int cutoff) const{
    casadi_assert_warning(ncol_==nrow_,"StarColoring requires a square matrix, but got " << dimString() << ".");
    // Reorder, if necessary
    if(ordering!=0){
      casadi_assert(ordering==1);
    
      // Ordering
      vector<int> ord = largestFirstOrdering();

      // Create a new sparsity pattern 
      CRSSparsity sp_permuted = pmult(ord,true,true,true);
    
      // Star coloring for the permuted matrix
      CRSSparsity ret_permuted = sp_permuted.starColoring(0);
        
      // Permute result back
      return ret_permuted.pmult(ord,false,true,false);
    }
  
    // Allocate temporary vectors
    vector<int> forbiddenColors;
    forbiddenColors.reserve(nrow_);
    vector<int> color(nrow_,-1);
    
    // 4: for i <- 1 to |V | do
    for(int i=0; i<nrow_; ++i){
        
      // 5: for each w \in N1 (vi ) do
      for(int w_el=rowind_[i]; w_el<rowind_[i+1]; ++w_el){
        int w = col_[w_el];
              
        // 6: if w is colored then
        if(color[w]!=-1){
        
          // 7: forbiddenColors[color[w]] <- v
          forbiddenColors[color[w]] = i;
        
        } // 8: end if

        // 9: for each colored vertex x \in N1 (w) do
        for(int x_el=rowind_[w]; x_el<rowind_[w+1]; ++x_el){
          int x = col_[x_el];
          if(color[x]==-1) continue;
        
          // 10: if w is not colored then
          if(color[w]==-1){
          
            //11: forbiddenColors[color[x]] <- vi
            forbiddenColors[color[x]] = i;
          
          } else { // 12: else
          
            // 13: for each colored vertex y \in N1 (x), y != w do
            for(int y_el=rowind_[x]; y_el<rowind_[x+1]; ++y_el){
              int y = col_[y_el];
              if(color[y]==-1 || y==w) continue;
            
              // 14: if color[y] = color[w] then
              if(color[y]==color[w]){
              
                // 15: forbiddenColors[color[x]] <- vi
                forbiddenColors[color[x]] = i;
              
                // 16: break
                break;
              
              } // 17: end if
            
            } // 18: end for
          
          } // 19: end if

        } // 20 end for
      
      } // 21 end for
    
      // 22: color[v] <- min{c > 0 : forbiddenColors[c] = v}
      bool new_color = true;
      for(int color_i=0; color_i<forbiddenColors.size(); ++color_i){
        // Break if color is ok
        if(forbiddenColors[color_i]!=i){
          color[i] = color_i;
          new_color = false;
          break;
        }
      }
    
      // New color if reached end
      if(new_color){
        color[i] = forbiddenColors.size();
        forbiddenColors.push_back(-1);

        // Cutoff if too many colors
        if(forbiddenColors.size()>cutoff){
          return CRSSparsity();
        }
      }
  
    } // 23 end for

    // Number of colors used
    int num_colors = forbiddenColors.size();

    // Return sparsity in sparse triplet format
    return sp_triplet(num_colors,nrow_,color,range(color.size()));
  }

  std::vector<int> CRSSparsityInternal::largestFirstOrdering() const{
    vector<int> degree = rowind_;
    int max_degree = 0;
    for(int k=0; k<nrow_; ++k){
      degree[k] = degree[k+1]-degree[k];
      max_degree = max(max_degree,1+degree[k]);
    }
    degree.resize(nrow_);

    // Vector for binary sort
    vector<int> degree_count(max_degree+1,0);
    for(vector<int>::const_iterator it=degree.begin(); it!=degree.end(); ++it){
      degree_count.at(*it+1)++;
    }

    // Cumsum to get the offset for each degree
    for(int d=0; d<max_degree; ++d){
      degree_count[d+1] += degree_count[d];
    }
  
    // Now a bucket sort
    vector<int> ordering(nrow_);
    for(int k=nrow_-1; k>=0; --k){
      ordering[degree_count[degree[k]]++] = k;
    }
  
    // Invert the ordering
    vector<int>& reverse_ordering = degree_count; // reuse memory
    reverse_ordering.resize(ordering.size());
    copy(ordering.begin(),ordering.end(),reverse_ordering.rbegin());
  
    // Return the ordering
    return reverse_ordering;
  }

  CRSSparsity CRSSparsityInternal::pmult(const std::vector<int>& p, bool permute_rows, bool permute_columns, bool invert_permutation) const{
    // Invert p, possibly
    vector<int> p_inv;
    if(invert_permutation){
      p_inv.resize(p.size());
      for(int k=0; k<p.size(); ++k){
        p_inv[p[k]] = k;
      }
    }
    const vector<int>& pp = invert_permutation ? p_inv : p;

    // Get rows
    vector<int> row = getRow();
  
    // Sparsity of the return matrix
    vector<int> new_col(col_.size()), new_row(row.size());

    // Possibly permute rows
    if(permute_rows){
      // Assert dimensions
      casadi_assert(p.size()==nrow_);
    
      // Permute
      for(int k=0; k<row.size(); ++k){
        new_row[k] = pp[row[k]];
      }
    
    } else {
      // No permutation of rows
      copy(row.begin(),row.end(),new_row.begin());
    }
  
    // Possibly permute columns
    if(permute_columns){
      // Assert dimensions
      casadi_assert(p.size()==ncol_);
    
      // Permute
      for(int k=0; k<col_.size(); ++k){
        new_col[k] = pp[col_[k]];
      }
    
    } else {
      // No permutation of columns
      copy(col_.begin(),col_.end(),new_col.begin());
    }
  
    // Return permuted matrix
    return sp_triplet(nrow_,ncol_,new_row,new_col);
  }

  bool CRSSparsityInternal::isTranspose(const CRSSparsityInternal& y) const{
    // Assert dimensions and number of nonzeros
    if(nrow_!=y.ncol_ || ncol_!=y.nrow_ || size()!=y.size())
      return false;
  
    // Quick return if empty or dense
    if(size()==0 || dense())
      return true;
    
    // Run algorithm on the pattern with the least number of columns
    if(ncol_>nrow_) return y.isTranspose(*this);

    // Index counter for row of the possible transpose
    vector<int> y_row_count(y.nrow_,0);
  
    // Loop over the rows
    for(int i=0; i<nrow_; ++i){
    
      // Loop over the nonzeros
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){
     
        // Get the column
        int j=col_[el];
      
        // Get the element of the possible transpose
        int el_y = y.rowind_[j] + y_row_count[j]++;
      
        // Quick return if element doesn't exist
        if(el_y>=y.rowind_[j+1]) return false;
      
        // Get the column of the possible transpose
        int j_y = y.col_[el_y];
      
        // Quick return if mismatch
        if(j_y != i) return false;
      }
    }
  
    // Transpose if reached this point
    return true;
  }

  void CRSSparsityInternal::spyMatlab(const std::string& mfile_name) const{
    // Create the .m file
    ofstream mfile;
    mfile.open(mfile_name.c_str());

    // Header
    mfile << "% This function was automatically generated by CasADi" << endl;
 
    // Print dimensions
    mfile << "m = " << nrow_ << ";" << endl;
    mfile << "n = " << ncol_ << ";" << endl;

    // Matlab indices are one-based
    const int index_offset = 1;
  
    // Print rows
    mfile << "i = [";
    bool first = true;
    for(int i=0; i<nrow_; ++i){
      for(int el=rowind_[i]; el<rowind_[i+1]; ++el){
        if(!first) mfile << ",";
        mfile << (i+index_offset);
        first = false;
      }
    }
    mfile << "];" << endl;
    
    // Print columns
    mfile << "j = [";
    first = true;
    for(vector<int>::const_iterator j=col_.begin(); j!=col_.end(); ++j){
      if(!first) mfile << ",";
      mfile << (*j+index_offset);
      first = false;
    }
    mfile << "];" << endl;
  
    // Generate matrix
    mfile << "A = sparse(i,j,ones(size(i)),m,n);" << endl; 
  
    // Issue spy command
    mfile << "spy(A);" << endl; 

    mfile.close();
  }

  std::size_t CRSSparsityInternal::hash() const{
    return hash_sparsity(nrow_,ncol_,col_,rowind_);
  }

  
} // namespace CasADi


