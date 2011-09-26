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
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
  int i, p, p2, done, jnew, head = 0;
  
  // initialize the recursion stack
  xi [0] = j;
  while (head >= 0){
    
    // get j from the top of the recursion stack 
    j = xi[head];
    jnew = !pinv.empty() ? (pinv[j]) : j;
    if (!marked[j]){
      
      // mark node j as visited
      marked[j]=true;
      pstack[head] = (jnew < 0) ? 0 : rowind_[jnew];
    }
    
    // node j done if no unvisited neighbors
    done = 1;
    p2 = (jnew < 0) ? 0 : rowind_[jnew+1];
    
    // examine all neighbors of j
    for (p = pstack[head] ; p < p2 ; p++){

      // consider neighbor node i
      i = col_[p];
      
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

void CRSSparsityInternal::breadthFirstSearch(int n, int *wi, int *wj, int *queue, const int *imatch, const int *jmatch, int mark) const{
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
  int head = 0, tail = 0, j, i, p, j2 ;
  
  // place all unmatched nodes in queue
  for (j=0; j<n; ++j){
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
  
  const int *Ap = &C->rowind_.front();
  const int *Ai = &C->col_.front();
  
  // while queue is not empty
  while (head < tail){
    
    // get the head of the queue
    j = queue[head++];
    for(p = Ap [j] ; p < Ap [j+1] ; p++){
      i = Ai [p] ;
      
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

void CRSSparsityInternal::matched(int n, const int *wj, const int *imatch, int *p, int *q, int *cc, int *rr, int set, int mark){
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

void CRSSparsityInternal::unmatched(int m, const int *wi, int *p, int *rr, int set){
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
  int i, kr = rr[set] ;
  for (i=0; i<m; i++) 
    if (wi[i] == 0) 
      p[kr++] = i;
    
  rr[set+1] = kr;
}

int CRSSparsityInternal::rprune(int i, int j, double aij, void *other){
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
  int *rr = (int *) other;
  return (i >= rr[1] && i < rr[2]) ;
}

void CRSSparsityInternal::augmentingPath(int k, std::vector<int>& jmatch, int *cheap, std::vector<int>& w, int *js, int *is, int *ps) const{
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

  int found = 0, p, i = -1, head = 0, j ;
  
  const int *Ap = &rowind_.front();
  const int *Ai = &col_.front();
  
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
      for(p = cheap [j] ; p < Ap[j+1] && !found; ++p){
        i = Ai [p] ;            /* try a cheap assignment (i,j) */
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
      ps[head] = Ap[j];
    }
    
    // --- Depth-first-search of neighbors of j -------------------------
    for(p = ps[head]; p<Ap[j+1]; ++p){
      
      // consider row i
      i = Ai[p];
      
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
    if(p == Ap[j+1]) head--;
  } // augment the match if path found:
  
  if(found)
    for(p = head; p>=0; --p)
      jmatch[is[p]] = js[p];
}

void CRSSparsityInternal::maxTransversal(std::vector<int>& imatch, std::vector<int>& jmatch, CRSSparsity& trans, int seed) const{
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

  int n2 = 0, m2 = 0;
  
  //cs *C ;
  int n = nrow_;
  int m = ncol_;
  const int *Ap = &rowind_.front();
  const int *Ai = &col_.front();
  
  // allocate result
  vector<int> jimatch(m+n);
  jmatch.resize(m);
  imatch.resize(n);
  vector<int> w(m+n);
  
  // count nonempty rows and columns
  int k=0;
  for(int j=0; j<n; ++j){
    n2 += (Ap[j] < Ap[j+1]);
    for(int p=Ap[j]; p < Ap[j+1]; ++p){
      w[Ai[p]] = 1;
      
      // count entries already on diagonal
      k += (j == Ai [p]);
    }
  }
  
  // quick return if diagonal zero-free
  if(k == std::min(m,n)){
    int i;
    for(i=0; i<k; ++i) jmatch[i] = i;
    for(;    i<m; ++i) jmatch[i] = -1;

    int j;
    for(j=0; j<k; ++j) imatch[j] = j;
    for(;    j<n; ++j) imatch[j] = -1;
  }

  for(int i=0; i<m; ++i) m2 += w[i];
  
  // transpose if needed
  if(m2 < n2 && trans.isNull())
    trans = transpose();
  
  // Get pointer to sparsity
  const CRSSparsityInternal* C = m2 < n2 ? static_cast<const CRSSparsityInternal*>(trans.get()) : this;
  
  n = C->nrow_;
  m = C->ncol_;
  const int* Cp = &C->rowind_.front();

  std::vector<int>& Cjmatch = m2 < n2 ? imatch : jmatch;
  std::vector<int>& Cimatch = m2 < n2 ? jmatch : imatch;
  
  // get workspace 
  w.resize(5*n);

  int *cheap = &w.front() + n;
  int *js = &w.front() + 2*n;
  int *is = &w.front() + 3*n; 
  int *ps = &w.front() + 4*n;

  // for cheap assignment
  for(int j=0; j<n; ++j) 
    cheap[j] = Cp[j];
  
  // all columns unflagged 
  for(int j=0; j<n; ++j)
    w[j] = -1;
  
  // nothing matched yet
  for(int i=0; i<m; ++i)
    jmatch[i] = -1;

  // q = random permutation 
  std::vector<int> q = randomPermutation(n,seed);

  // augment, starting at column q[k]
  for(k=0; k<n; ++k){
    C->augmentingPath(!q.empty() ? q[k]: k, jmatch, cheap, w, js, is, ps);
  }

  // find row match
  for(int j=0; j<n; ++j)
    imatch[j] = -1;
  
  for(int i = 0; i<m; ++i)
    if(jmatch [i] >= 0)
      imatch[jmatch[i]] = i;
}

int CRSSparsityInternal::dulmageMendelsohn(std::vector<int>& rowperm, std::vector<int>& colperm, std::vector<int>& rowblock, std::vector<int>& colblock, std::vector<int>& coarse_rowblock, std::vector<int>& coarse_colblock, int seed) const{
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

  int k, cnz, nc, nb2, ok ;

  // The transpose of the expression
  CRSSparsity trans;
  
  
  //  csd *D, *scc ;

  // Part 1: Maximum matching
  int m = ncol_;
  int n = nrow_;

  // row permutation 
  vector<int> p(m);
  
  // column permutation 
  vector<int> q(n);
  
  // size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q)
  vector<int> r(m+6);
  
  // size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q)
  vector<int> s(n+6);

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
  
  // use r and s as workspace
  vector<int>& wi = r;
  vector<int>& wj = s;
  
  // unmark all cols for bfs
  for(int j=0; j<n; ++j)
    wj[j] = -1;
  
  // unmark all rows for bfs
  for(int i=0; i<m; ++i)
    wi[i] = -1 ;
  
  // find C1, R1 from C0
  breadthFirstSearch(n, &wi.front(), &wj.front(), &q.front(), &imatch.front(), &jmatch.front(), 1);

  // find R3, C3 from R0
  breadthFirstSearch(m, &wj.front(), &wi.front(), &p.front(), &jmatch.front(), &imatch.front(), 3);

  // unmatched set C0
  unmatched(n, &wj.front(), &q.front(), &coarse_rowblock.front(), 0);

  // set R1 and C1
  matched(n, &wj.front(), &imatch.front(), &p.front(), &q.front(), &coarse_rowblock.front(), &coarse_colblock.front(), 1, 1);

  // set R2 and C2
  matched(n, &wj.front(), &imatch.front(), &p.front(), &q.front(), &coarse_rowblock.front(), &coarse_colblock.front(), 2, -1);

  // set R3 and C3
  matched(n, &wj.front(), &imatch.front(), &p.front(), &q.front(), &coarse_rowblock.front(), &coarse_colblock.front(), 3, 3);

  // unmatched set R0
  unmatched(m, &wi.front(), &p.front(), &coarse_colblock.front(), 3);
  
  // --- Fine decomposition -----------------------------------------------
  // pinv=p'
  vector<int> pinv = invertPermutation(&p.front(), m);

  // C=A(p,q) (it will hold A(R2,C2))
  CRSSparsity C = permute(&pinv.front(), &q.front(), 0);

  vector<int>& Cp = C.rowindRef();

  // delete cols C0, C1, and C3 from C 
  nc = coarse_rowblock[3] - coarse_rowblock[2];
  if(coarse_rowblock[2] > 0)
    for(int j = coarse_rowblock[2]; j <= coarse_rowblock[3]; ++j)
      Cp[j-coarse_rowblock[2]] = Cp[j];
  
  C->nrow_ = nc;
  C->rowind_.resize(nc+1);
  C->col_.resize(C->rowind_.back());

  // delete rows R0, R1, and R3 from C
  if(coarse_colblock[2] - coarse_colblock[1] < m){
    C->drop(rprune, &coarse_colblock.front()) ;
    cnz = Cp[nc];
    vector<int>& Ci = C->col_;
    if(coarse_colblock[1] > 0)
      for(k=0; k<cnz; ++k)
        Ci[k] -= coarse_colblock[1];
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

  for(k=0; k<nc; ++k)
    wj [k] = q[ps[k] + coarse_rowblock[2]];
  
  for(k=0; k<nc; ++k)
    q[k + coarse_rowblock[2]] = wj[k];
  
  for(k=0; k<nc; ++k)
    wi[k] = p[ps[k] + coarse_colblock[1]];
  
  for(k=0; k<nc; ++k)
    p[k + coarse_colblock[1]] = wi[k];
  
  // create the fine block partitions
  nb2 = 0;
  r[0] = s[0] = 0;

  // leading coarse block A (R1, [C0 C1])
  if(coarse_rowblock[2] > 0)
    nb2++ ;
  
  // coarse block A (R2,C2)
  for (k=0; k<nb1; ++k){
    // A (R2,C2) splits into nb1 fine blocks 
    r [nb2] = rs [k] + coarse_colblock[1];
    s [nb2] = rs [k] + coarse_rowblock[2] ;
    nb2++ ;
  }
  
  if(coarse_colblock[2] < m){
    // trailing coarse block A ([R3 R0], C3)
    r[nb2] = coarse_colblock[2];
    s[nb2] = coarse_rowblock[3];
    nb2++ ;
  }
  
  r[nb2] = m ;
  s[nb2] = n ;
  
  
  // Copy to output
  rowperm = vector<int>(q.begin(), q.end());
  colperm = vector<int>(p.begin(), p.end());
  rowblock = vector<int>(s.begin(), s.end());
  colblock = vector<int>(r.begin(), r.end());
  return nb2;
}

std::vector<int> CRSSparsityInternal::randomPermutation(int n, int seed){
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

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

std::vector<int> CRSSparsityInternal::invertPermutation(const int *p, int n){
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

// pinv = p', or p = pinv'

  // allocate result
  vector<int> pinv(n);
  
  // invert the permutation
  for(int k=0; k<n; ++k)
    pinv[p[k]] = k;
  
  // return result
  return pinv;
}

CRSSparsity CRSSparsityInternal::permute(const int *pinv, const int *q, int values) const{
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++

int m = ncol_;
  int n = nrow_;
  const int *Ap = &rowind_.front();
  const int *Ai = &col_.front();
  
  // alloc result
  CRSSparsity C = CRSSparsity(nrow_,ncol_);
  
  // Row offset
  vector<int>& Cp = C.rowindRef();
  
  // Column for each nonzero
  vector<int>& Ci = C.colRef();
  Ci.resize(size());

  int nz = 0;
  for(int k = 0; k<n ; ++k){
    // column k of C is column q[k] of A
    Cp[k] = nz;
    
    int j = q ? (q[k]) : k;
    
    for(int t = Ap[j]; t<Ap[j+1]; ++t){
      Ci[nz++] = pinv ? (pinv[Ai[t]]) : Ai[t] ;
    }
  }
  
  // finalize the last column of C
  Cp[n] = nz;
  return C;
}

int CRSSparsityInternal::drop(int (*fkeep) (int, int, double, void *), void *other){
  // NOTE: This implementation has been copied from CSparse and then modified, it needs cleaning up to be proper C++
  
  int j, p, nz = 0, n, *Ap, *Ai;
  n = nrow_;
  Ap = &rowind_.front();
  Ai = &col_.front();
  
  for (j = 0 ; j < n ; j++){
    // get current location of col j
    p = Ap[j];
    
    // record new location of col j
    Ap[j] = nz;
    for ( ; p < Ap [j+1] ; ++p){
      if (fkeep (Ai [p], j, 1, other)){
        // keep A(i,j)
        Ai[nz++] = Ai[p] ;
      }
    }
  }
  
  // finalize A
  Ap[n] = nz;
  return nz ;
}


} // namespace CasADi


