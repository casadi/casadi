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

using namespace std;

namespace CasADi{

int CRSSparsityInternal::size() const{
  return col_.size();
}
    
void CRSSparsityInternal::repr(ostream &stream) const{
  stream << "Compressed Row Storage: " << nrow_ << "-by-" << ncol_ << " matrix, " << col_.size() << " structural non-zeros";
}

void CRSSparsityInternal::sanityCheck(bool complete) const{
  if (rowind_.size() != nrow_+1) {
    std::stringstream s;
    s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
    s << "  rowind.size() = nrow + 1, but got   rowind.size() = " << rowind_.size() << "   and   nrow = "  << nrow_ << std::endl;
    s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
    throw CasadiException(s.str());
  }
  if (complete) {
  
    if (rowind_.size()>0) {
      for (int k=1;k<rowind_.size();k++) {
        if (rowind_[k]<rowind_[k-1]) {
          throw CasadiException("CRSSparsityInternal:Compressed Row Storage is not sane. rowind must be monotone. Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind).");
        }
      }
      
      if (rowind_[0]!=0) {
        throw CasadiException("CRSSparsityInternal:Compressed Row Storage is not sane. First element of rowind must be zero. Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind).");
      }
      if (rowind_[(rowind_.size()-1)]!=col_.size()) {
        std::stringstream s;
        s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
        s << "  rowind[lastElement] = col.size(), but got   rowind[lastElement] = " << rowind_[(rowind_.size()-1)] << "   and   col.size() = "  << col_.size() << std::endl;
        s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
        throw CasadiException(s.str());
      }
      if (col_.size()>nrow_*ncol_) {
        std::stringstream s;
        s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
        s << "  col.size() <= nrow * ncol, but got   col.size()  = " << col_.size() << "   and   nrow * ncol = "  << nrow_*ncol_ << std::endl;
        s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
        throw CasadiException(s.str());
      }
    }
    for (int k=0;k<col_.size();k++) {
      if (col_[k]>=ncol_ || col_[k] < 0) {
        std::stringstream s;
        s << "CRSSparsityInternal:Compressed Row Storage is not sane. The following must hold:" << std::endl;
        s << "  0 <= col[i] < ncol for each i, but got   col[i] = " << col_[k] << "   and   ncol = "  << ncol_ << std::endl;
        s << "  Note that the signature is as follows: CRSSparsity (nrow, ncol, col, rowind)." << std::endl;
        throw CasadiException(s.str());
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

CRSSparsity CRSSparsityInternal::transpose(vector<int>& mapping) const{
  // Get the sparsity of the transpose in sparse triplet form
  const vector<int>& trans_row = col_;
  vector<int> trans_col = getRow();

  // Create the sparsity pattern
  return sp_triplet(ncol_,nrow_,trans_row,trans_col,mapping);

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
  
  vector<int> rcopy(nrow_+1);
  vector<int>& pstack = rcopy;
  
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
  
  for(int b=0; b <= nb; ++b){
    rcopy[b] = r[b] ;
  }
  
  for(int i=0; i<nrow_; ++i){
    p[rcopy[Blk[i]]++] = i;
  }
  
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
  vector<int> jimatch(ncol_+nrow_);
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
    jmatch[i] = -1;

  // q = random permutation 
  std::vector<int> q = randomPermutation(C->nrow_,seed);

  // augment, starting at column q[k]
  for(k=0; k<C->nrow_; ++k){
    C->augmentingPath(!q.empty() ? q[k]: k, jmatch, cheap, w, js, is, ps);
  }

  // find row match
  for(int j=0; j<C->nrow_; ++j)
    imatch[j] = -1;
  
  for(int i = 0; i<C->ncol_; ++i)
    if(jmatch [i] >= 0)
      imatch[jmatch[i]] = i;
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
  if(coarse_rowblock[2] > 0)
    for(int j = coarse_rowblock[2]; j <= coarse_rowblock[3]; ++j)
      rowind_C[j-coarse_rowblock[2]] = rowind_C[j];
  
  C->nrow_ = nc;
  C->rowind_.resize(nc+1);
  C->col_.resize(C->rowind_.back());

  // delete rows R0, R1, and R3 from C
  if(coarse_colblock[2] - coarse_colblock[1] < ncol_){
    C->drop(rprune, &coarse_colblock);
    int cnz = rowind_C[nc];
    vector<int>& col_C = C->col_;
    if(coarse_colblock[1] > 0)
      for(int k=0; k<cnz; ++k)
        col_C[k] -= coarse_colblock[1];
  }
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
    pinv [i] = k;
    
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

#if 0


/* symbolic ordering and analysis for QR or LU */
css *cs_sqr (int order, const cs *A, int qr)
{
    int n, k, ok = 1, *post ;
    css *S ;
    if (!CS_CSC (A)) return (NULL) ;        /* check inputs */
    n = A->n ;
    S = cs_calloc (1, sizeof (css)) ;       /* allocate result S */
    if (!S) return (NULL) ;                 /* out of memory */
    S->q = cs_amd (order, A) ;              /* fill-reducing ordering */
    if (order && !S->q) return (cs_sfree (S)) ;
    if (qr)                                 /* QR symbolic analysis */
    {
        cs *C = order ? cs_permute (A, NULL, S->q, 0) : ((cs *) A) ;
        S->parent = cs_etree (C, 1) ;       /* etree of C'*C, where C=A(:,q) */
        post = cs_post (S->parent, n) ;
        S->cp = cs_counts (C, S->parent, post, 1) ;  /* col counts chol(C'*C) */
        cs_free (post) ;
        ok = C && S->parent && S->cp && cs_vcount (C, S) ;
        if (ok) for (S->unz = 0, k = 0 ; k < n ; k++) S->unz += S->cp [k] ;
        ok = ok && S->lnz >= 0 && S->unz >= 0 ;     /* int overflow guard */
        if (order) cs_spfree (C) ;
    }
    else
    {
        S->unz = 4*(A->p [n]) + n ;         /* for LU factorization only, */
        S->lnz = S->unz ;                   /* guess nnz(L) and nnz(U) */
    }
    return (ok ? S : cs_sfree (S)) ;        /* return result S */
}
#endif





















} // namespace CasADi


