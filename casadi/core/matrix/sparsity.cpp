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
#include "../matrix/matrix.hpp"
#include "../std_vector_tools.hpp"
#include <climits>

using namespace std;

namespace casadi {

  /// \cond INTERNAL
  // Singletons
  class EmptySparsity : public Sparsity {
  public:
    EmptySparsity() {
      std::vector<int> colind(1, 0), row;
      assignNode(new SparsityInternal(0, 0, colind, row));
    }
  };

  class ScalarSparsity : public Sparsity {
  public:
    ScalarSparsity() {
      std::vector<int> colind(2), row(1, 0);
      colind[0] = 0;
      colind[1] = 1;
      assignNode(new SparsityInternal(1, 1, colind, row));
    }
  };

  class ScalarSparseSparsity : public Sparsity {
  public:
    ScalarSparseSparsity() {
      std::vector<int> colind(2, 0), row;
      assignNode(new SparsityInternal(1, 1, colind, row));
    }
  };
  /// \endcond

  Sparsity::Sparsity(int dummy) {
    casadi_assert(dummy==0);
  }

  Sparsity Sparsity::create(SparsityInternal *node) {
    Sparsity ret;
    ret.assignNode(node);
    return ret;
  }

  Sparsity::Sparsity(int nrow, int ncol, const std::vector<int>& colind,
                     const std::vector<int>& row) {
    assignCached(nrow, ncol, colind, row);
    sanityCheck(true);
  }

  void Sparsity::reCache() {
    assignCached(size1(), size2(), colind(), row());
  }

  SparsityInternal* Sparsity::operator->() {
    makeUnique();
    return static_cast<SparsityInternal*>(SharedObject::operator->());
  }

  const SparsityInternal* Sparsity::operator->() const {
    return static_cast<const SparsityInternal*>(SharedObject::operator->());
  }

  SparsityInternal& Sparsity::operator*() {
    makeUnique();
    return *static_cast<SparsityInternal*>(get());
  }

  const SparsityInternal& Sparsity::operator*() const {
    return *static_cast<const SparsityInternal*>(get());
  }

  bool Sparsity::testCast(const SharedObjectNode* ptr) {
    return dynamic_cast<const SparsityInternal*>(ptr)!=0;
  }

  int Sparsity::size1() const {
    return (*this)->nrow_;
  }

  int Sparsity::size2() const {
    return (*this)->ncol_;
  }

  int Sparsity::numel() const {
    return (*this)->numel();
  }

  bool Sparsity::isEmpty(bool both) const {
    return (*this)->isEmpty(both);
  }

  int Sparsity::size() const {
    return (*this)->size();
  }

  std::pair<int, int> Sparsity::shape() const {
    return (*this)->shape();
  }

  const std::vector<int>& Sparsity::row() const {
    return (*this)->row_;
  }

  const std::vector<int>& Sparsity::colind() const {
    return (*this)->colind_;
  }

  std::vector<int>& Sparsity::rowRef() {
    makeUnique();
    return (*this)->row_;
  }

  std::vector<int>& Sparsity::colindRef() {
    makeUnique();
    return (*this)->colind_;
  }

  int Sparsity::row(int el) const {
    return row().at(el);
  }

  int Sparsity::colind(int col) const {
    return colind().at(col);
  }

  void Sparsity::sanityCheck(bool complete) const {
    (*this)->sanityCheck(complete);
  }

  void Sparsity::resize(int nrow, int ncol) {
    makeUnique();
    (*this)->resize(nrow, ncol);
  }

  int Sparsity::getNZ(int rr, int cc) {
    // If negative index, count from the back
    if (rr<0) rr += size1();
    if (cc<0) cc += size2();

    // Check consistency
    casadi_assert_message(rr>=0 && rr<size1(), "Row index out of bounds");
    casadi_assert_message(cc>=0 && cc<size2(), "Column index out of bounds");

    // Quick return if matrix is dense
    if (isDense()) return rr+cc*size1();

    // Quick return if we are adding an element to the end
    if (colind(cc)==size() || (colind(cc+1)==size() && row().back()<rr)) {
      std::vector<int>& rowv = rowRef();
      std::vector<int>& colindv = colindRef();
      rowv.push_back(rr);
      for (int c=cc; c<size2(); ++c) {
        colindv[c+1]++;
      }
      return rowv.size()-1;
    }

    // go to the place where the element should be
    int ind;
    for (ind=colind(cc); ind<colind(cc+1); ++ind) { // better: loop from the back to the front
      if (row(ind) == rr) {
        return ind; // element exists
      } else if (row(ind) > rr) {
        break;                // break at the place where the element should be added
      }
    }

    // insert the element
    std::vector<int>& rowv = rowRef();
    std::vector<int>& colindv = colindRef();
    rowv.insert(rowv.begin()+ind, rr);
    for (int c=cc+1; c<size2()+1; ++c)
      colindv[c]++;

    // Return the location of the new element
    return ind;
  }

  bool Sparsity::hasNZ(int rr, int cc) const {
    return (*this)->getNZ(rr, cc)!=-1;
  }


  int Sparsity::getNZ(int rr, int cc) const {
    return (*this)->getNZ(rr, cc);
  }

  Sparsity Sparsity::reshape(int nrow, int ncol) const {
    return (*this)->reshape(nrow, ncol);
  }

  std::vector<int> Sparsity::getNZ(const std::vector<int>& rr, const std::vector<int>& cc) const {
    return (*this)->getNZ(rr, cc);
  }

  bool Sparsity::isScalar(bool scalar_and_dense) const {
    return (*this)->isScalar(scalar_and_dense);
  }

  bool Sparsity::isDense() const {
    return (*this)->isDense();
  }

  bool Sparsity::isDiagonal() const {
    return (*this)->isDiagonal();
  }

  bool Sparsity::isSquare() const {
    return (*this)->isSquare();
  }

  bool Sparsity::isSymmetric() const {
    return (*this)->isSymmetric();
  }

  bool Sparsity::isTril() const {
    return (*this)->isTril();
  }

  bool Sparsity::isTriu() const {
    return (*this)->isTriu();
  }

  Sparsity Sparsity::sub(const std::vector<int>& jj, const std::vector<int>& ii,
                         std::vector<int>& mapping) const {
    return (*this)->sub(jj, ii, mapping);
  }

  std::vector<int> Sparsity::erase(const std::vector<int>& jj, const std::vector<int>& ii) {
    makeUnique();
    return (*this)->erase(jj, ii);
  }

  int Sparsity::sizeL() const {
    return (*this)->sizeL();
  }

  int Sparsity::sizeU() const {
    return (*this)->sizeU();
  }

  int Sparsity::sizeD() const {
    return (*this)->sizeD();
  }

  std::vector<int> Sparsity::getCol() const {
    return (*this)->getCol();
  }

  void Sparsity::getCCS(std::vector<int>& colind, std::vector<int>& row) const {
    colind = this->colind();
    row = this->row();
  }

  void Sparsity::getCRS(std::vector<int>& rowind, std::vector<int>& col) const {
    transpose().getCCS(rowind, col);
  }


  void Sparsity::getTriplet(std::vector<int>& row, std::vector<int>& col) const {
    row = this->row();
    col = this->getCol();
  }

  Sparsity Sparsity::transpose(std::vector<int>& mapping, bool invert_mapping) const {
    return (*this)->transpose(mapping, invert_mapping);
  }

  Sparsity Sparsity::transpose() const {
    return (*this)->transpose();
  }

  Sparsity Sparsity::patternCombine(const Sparsity& y, bool f0x_is_zero,
                                    bool function0_is_zero,
                                    std::vector<unsigned char>& mapping) const {
    return (*this)->patternCombine(y, f0x_is_zero, function0_is_zero, mapping);
  }

  Sparsity Sparsity::patternCombine(const Sparsity& y, bool f0x_is_zero,
                                    bool function0_is_zero) const {
    return (*this)->patternCombine(y, f0x_is_zero, function0_is_zero);
  }

  Sparsity Sparsity::patternUnion(const Sparsity& y, std::vector<unsigned char>& mapping) const {
    return (*this)->patternCombine(y, false, false, mapping);
  }

  Sparsity Sparsity::patternUnion(const Sparsity& y) const {
    return (*this)->patternCombine(y, false, false);
  }

  Sparsity Sparsity::patternIntersection(const Sparsity& y,
                                         std::vector<unsigned char>& mapping) const {
    return (*this)->patternCombine(y, true, true, mapping);
  }

  Sparsity Sparsity::patternIntersection(const Sparsity& y) const {
    return (*this)->patternCombine(y, true, true);
  }

  Sparsity Sparsity::patternProductNew(const Sparsity& y) const {
    return (*this)->patternProductNew(y);
  }

  Sparsity Sparsity::patternProduct(const Sparsity& x_trans) const {
    return x_trans.T().patternProductNew(*this);
  }

  Sparsity Sparsity::patternProduct(const Sparsity& x_trans,
                                    std::vector< std::vector< pair<int, int> > >& mapping) const {
    return (*this)->patternProduct(x_trans, mapping);
  }

  bool Sparsity::isEqual(const Sparsity& y) const {
    return (*this)->isEqual(y);
  }

  bool Sparsity::isEqual(int nrow, int ncol, const std::vector<int>& colind,
                         const std::vector<int>& row) const {
    return (*this)->isEqual(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::operator+(const Sparsity& b) const {
    return patternUnion(b);
  }

  Sparsity Sparsity::operator*(const Sparsity& b) const {
    std::vector< unsigned char > mapping;
    return patternIntersection(b, mapping);
  }

  Sparsity Sparsity::patternInverse() const {
    return (*this)->patternInverse();
  }

  void Sparsity::reserve(int nnz, int ncol) {
    makeUnique();
    (*this)->reserve(nnz, ncol);
  }

  void Sparsity::append(const Sparsity& sp) {
    if (sp.size1()==0 && sp.size2()==0) {
      // Appending pattern is empty
      return;
    } else if (size1()==0 && size2()==0) {
      // This is empty
      *this = sp;
    } else {
      casadi_assert_message(size2()==sp.size2(),
                            "Sparsity::append: Dimension mismatch. "
                            "You attempt to append a shape " << sp.dimString()
                            << " to a shape " << dimString()
                            << ". The number of columns must match.");
      if (sp.size1()==0) {
        // No rows to add
        return;
      } else if (size1()==0) {
        // No rows before
        *this = sp;
      } else if (isVector()) {
        // Append to vector (efficient)
        makeUnique();
        (*this)->append(*sp);
      } else {
        // Append to matrix (inefficient)
        *this = vertcat(*this, sp);
      }
    }
  }

  void Sparsity::appendColumns(const Sparsity& sp) {
    if (sp.size1()==0 && sp.size2()==0) {
      // Appending pattern is empty
      return;
    } else if (size1()==0 && size2()==0) {
      // This is empty
      *this = sp;
    } else {
      casadi_assert_message(size1()==sp.size1(),
                            "Sparsity::appendColumns: Dimension mismatch. You attempt to "
                            "append a shape " << sp.dimString() << " to a shape "
                            << dimString() << ". The number of rows must match.");
      if (sp.size2()==0) {
        // No columns to add
        return;
      } else if (size2()==0) {
        // No columns before
        *this = sp;
      } else {
        // Append to matrix (efficient)
        makeUnique();
        (*this)->appendColumns(*sp);
      }
    }
  }

  Sparsity::CachingMap& Sparsity::getCache() {
    static CachingMap ret;
    return ret;
  }

  const Sparsity& Sparsity::getScalar() {
    static ScalarSparsity ret;
    return ret;
  }

  const Sparsity& Sparsity::getScalarSparse() {
    static ScalarSparseSparsity ret;
    return ret;
  }

  const Sparsity& Sparsity::getEmpty() {
    static EmptySparsity ret;
    return ret;
  }

  void Sparsity::enlarge(int nrow, int ncol, const std::vector<int>& jj,
                         const std::vector<int>& ii) {
    enlargeColumns(ncol, ii);
    enlargeRows(nrow, jj);
  }

  void Sparsity::enlargeColumns(int ncol, const std::vector<int>& ii) {
    makeUnique();
    (*this)->enlargeColumns(ncol, ii);
  }

  void Sparsity::enlargeRows(int nrow, const std::vector<int>& jj) {
    makeUnique();
    (*this)->enlargeRows(nrow, jj);
  }

  Sparsity Sparsity::diag(int nrow, int ncol) {
    // Smallest dimension
    int n = min(nrow, ncol);

    // Column offset
    vector<int> colind(ncol+1, n);
    for (int cc=0; cc<n; ++cc) colind[cc] = cc;

    // Row
    vector<int> row = range(n);

    // Create pattern from vectors
    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::makeDense(std::vector<int>& mapping) const {
    return (*this)->makeDense(mapping);
  }

  std::string Sparsity::dimString() const {
    return (*this)->dimString();
  }

  Sparsity Sparsity::getDiag(std::vector<int>& mapping) const {
    return (*this)->getDiag(mapping);
  }

  std::vector<int> Sparsity::eliminationTree(bool ata) const {
    return (*this)->eliminationTree(ata);
  }

  int Sparsity::depthFirstSearch(int j, int top, std::vector<int>& xi,
                                 std::vector<int>& pstack, const std::vector<int>& pinv,
                                 std::vector<bool>& marked) const {
    return (*this)->depthFirstSearch(j, top, xi, pstack, pinv, marked);
  }

  int Sparsity::stronglyConnectedComponents(std::vector<int>& p, std::vector<int>& r) const {
    return (*this)->stronglyConnectedComponents(p, r);
  }

  int Sparsity::dulmageMendelsohn(std::vector<int>& rowperm, std::vector<int>& colperm,
                                  std::vector<int>& rowblock, std::vector<int>& colblock,
                                  std::vector<int>& coarse_rowblock,
                                  std::vector<int>& coarse_colblock, int seed) const {
    return (*this)->dulmageMendelsohn(rowperm, colperm, rowblock, colblock,
                                      coarse_rowblock, coarse_colblock, seed);
  }

  bool Sparsity::rowsSequential(bool strictly) const {
    return (*this)->rowsSequential(strictly);
  }

  void Sparsity::removeDuplicates(std::vector<int>& mapping) {
    makeUnique();
    (*this)->removeDuplicates(mapping);
  }

  std::vector<int> Sparsity::getElements(bool col_major) const {
    std::vector<int> loc;
    getElements(loc, col_major);
    return loc;
  }

  void Sparsity::getElements(std::vector<int>& loc, bool col_major) const {
    (*this)->getElements(loc, col_major);
  }

  void Sparsity::getNZInplace(std::vector<int>& indices) const {
    (*this)->getNZInplace(indices);
  }

  Sparsity Sparsity::unidirectionalColoring(const Sparsity& AT, int cutoff) const {
    if (AT.isNull()) {
      return (*this)->unidirectionalColoring(transpose(), cutoff);
    } else {
      return (*this)->unidirectionalColoring(AT, cutoff);
    }
  }

  Sparsity Sparsity::starColoring(int ordering, int cutoff) const {
    return (*this)->starColoring(ordering, cutoff);
  }

  Sparsity Sparsity::starColoring2(int ordering, int cutoff) const {
    return (*this)->starColoring2(ordering, cutoff);
  }

  std::vector<int> Sparsity::largestFirstOrdering() const {
    return (*this)->largestFirstOrdering();
  }

  Sparsity Sparsity::pmult(const std::vector<int>& p, bool permute_rows, bool permute_cols,
                           bool invert_permutation) const {
    return (*this)->pmult(p, permute_rows, permute_cols, invert_permutation);
  }

  void Sparsity::spyMatlab(const std::string& mfile) const {
    (*this)->spyMatlab(mfile);
  }

  void Sparsity::spy(std::ostream &stream) const {
    (*this)->spy(stream);
  }

  bool Sparsity::isTranspose(const Sparsity& y) const {
    return (*this)->isTranspose(*y);
  }

  bool Sparsity::isReshape(const Sparsity& y) const {
    return (*this)->isReshape(*y);
  }

  std::size_t Sparsity::hash() const {
    return (*this)->hash();
  }

  void Sparsity::assignCached(int nrow, int ncol, const std::vector<int>& colind,
                              const std::vector<int>& row) {

    // Scalars and empty patterns are handled separately
    if (ncol==0 && nrow==0) {
      // If empty
      *this = getEmpty();
      return;
    } else if (ncol==1 && nrow==1) {
      if (row.empty()) {
        // If sparse scalar
        *this = getScalarSparse();
        return;
      } else {
        // If dense scalar
        *this = getScalar();
        return;
      }
    }

    // Hash the pattern
    std::size_t h = hash_sparsity(nrow, ncol, colind, row);

    // Get a reference to the cache
    CachingMap& cache = getCache();

    // Record the current number of buckets (for garbage collection below)
#ifdef USE_CXX11
    int bucket_count_before = cache.bucket_count();
#endif // USE_CXX11

    // WORKAROUND, functions do not appear to work when bucket_count==0
#ifdef USE_CXX11
    if (bucket_count_before>0) {
#endif // USE_CXX11

      // Find the range of patterns equal to the key (normally only zero or one)
      pair<CachingMap::iterator, CachingMap::iterator> eq = cache.equal_range(h);

      // Loop over maching patterns
      for (CachingMap::iterator i=eq.first; i!=eq.second; ++i) {

        // Get a weak reference to the cached sparsity pattern
        WeakRef& wref = i->second;

        // Check if the pattern still exists
        if (wref.alive()) {

          // Get an owning reference to the cached pattern
          Sparsity ref = shared_cast<Sparsity>(wref.shared());

          // Check if the pattern matches
          if (ref.isEqual(nrow, ncol, colind, row)) {

            // Found match!
            assignNode(ref.get());
            return;

          } else {
            // There are two options, either the pattern has changed or there is
            // a hash collision, so let's rehash the pattern
            std::size_t h_ref = ref.hash();

            if (h_ref!=h) { // The sparsity pattern has changed (the most likely event)

              // Create a new pattern
              assignNode(new SparsityInternal(nrow, ncol, colind, row));

              // Cache this pattern instead of the old one
              wref = *this;

              // Recache the old sparsity pattern
              // TODO(Joel): recache "ref"
              return;

            } else { // There is a hash rowision (unlikely, but possible)
              // Leave the pattern alone, continue to the next matching pattern
              continue;
            }
          }
        } else {
          // Check if one of the other cache entries indeed has a matching sparsity
          CachingMap::iterator j=i;
          j++; // Start at the next matching key
          for (; j!=eq.second; ++j) {
            if (j->second.alive()) {

              // Recover cached sparsity
              Sparsity ref = shared_cast<Sparsity>(j->second.shared());

              // Match found if sparsity matches
              if (ref.isEqual(nrow, ncol, colind, row)) {
                assignNode(ref.get());
                return;
              }
            }
          }

          // The cached entry has been deleted, create a new one
          assignNode(new SparsityInternal(nrow, ncol, colind, row));

          // Cache this pattern
          wref = *this;

          // Return
          return;
        }
      }

      // END WORKAROUND
#ifdef USE_CXX11
    }
#endif // USE_CXX11

    // No matching sparsity pattern could be found, create a new one
    assignNode(new SparsityInternal(nrow, ncol, colind, row));

    // Cache this pattern
    //cache.insert(eq.second, std::pair<std::size_t, WeakRef>(h, ret));
    cache.insert(std::pair<std::size_t, WeakRef>(h, *this));

    // Garbage collection (currently only supported for unordered_multimap)
#ifdef USE_CXX11
    int bucket_count_after = cache.bucket_count();

    // We we increased the number of buckets, take time to garbage-collect deleted references
    if (bucket_count_before!=bucket_count_after) {
      CachingMap::const_iterator i=cache.begin();
      while (i!=cache.end()) {
        if (!i->second.alive()) {
          i = cache.erase(i);
        } else {
          i++;
        }
      }
    }
#endif // USE_CXX11
  }

  void Sparsity::clearCache() {
    getCache().clear();
  }

  Sparsity Sparsity::getTril(bool includeDiagonal) const {
    return (*this)->getTril(includeDiagonal);
  }

  Sparsity Sparsity::getTriu(bool includeDiagonal) const {
    return (*this)->getTriu(includeDiagonal);
  }

  std::vector<int> Sparsity::getLowerNZ() const {
    return (*this)->getLowerNZ();
  }

  std::vector<int> Sparsity::getUpperNZ() const {
    return (*this)->getUpperNZ();
  }

  std::size_t hash_sparsity(int nrow, int ncol, const std::vector<int>& colind,
                            const std::vector<int>& row) {
    // Condense the sparsity pattern to a single, deterministric number
    std::size_t ret=0;
    hash_combine(ret, nrow);
    hash_combine(ret, ncol);
    hash_combine(ret, colind);
    hash_combine(ret, row);
    return ret;
  }

  Sparsity Sparsity::dense(int nrow, int ncol) {
    // Column offset
    std::vector<int> colind(ncol+1);
    for (int cc=0; cc<ncol+1; ++cc) colind[cc] = cc*nrow;

    // Row
    std::vector<int> row(ncol*nrow);
    for (int cc=0; cc<ncol; ++cc)
      for (int rr=0; rr<nrow; ++rr)
        row[rr+cc*nrow] = rr;

    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::sparse(int nrow, int ncol) {
    std::vector<int> row, colind(ncol+1, 0);
    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::triu(int n) {
    casadi_assert_message(n>=0, "Sparsity::triu expects a positive integer as argument");
    int nrow=n, ncol=n;
    std::vector<int> colind, row;
    colind.reserve(ncol+1);
    row.reserve((n*(n+1))/2);

    // Loop over columns
    colind.push_back(0);
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over rows for the upper triangular half
      for (int rr=0; rr<=cc; ++rr) {
        row.push_back(rr);
      }
      colind.push_back(row.size());
    }

    // Return the pattern
    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::tril(int n) {
    casadi_assert_message(n>=0, "Sparsity::tril expects a positive integer as argument");
    int nrow=n, ncol=n;
    std::vector<int> colind, row;
    colind.reserve(ncol+1);
    row.reserve((n*(n+1))/2);

    // Loop over columns
    colind.push_back(0);
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over rows for the lower triangular half
      for (int rr=cc; rr<nrow; ++rr) {
        row.push_back(rr);
      }
      colind.push_back(row.size());
    }

    // Return the pattern
    return Sparsity(nrow, ncol, colind, row);
  }

  Sparsity Sparsity::band(int n, int p) {
    casadi_assert_message(n>=0, "Sparsity::band expects a positive integer as argument");
    casadi_assert_message((p<0? -p : p)<n,
                          "Sparsity::band: position of band schould be smaller then size argument");

    int nc = n-(p<0? -p : p);

    std::vector< int >          row(nc);

    int offset = max(p, 0);
    for (int i=0;i<nc;i++) {
      row[i]=i+offset;
    }

    std::vector< int >          colind(n+1);

    offset = min(p, 0);
    for (int i=0;i<n+1;i++) {
      colind[i]=max(min(i+offset, nc), 0);
    }

    return Sparsity(n, n, colind, row);

  }

  Sparsity Sparsity::banded(int n, int p) {
    // This is not an efficient implementation
    Sparsity ret = Sparsity::sparse(n, n);
    for (int i=-p;i<=p;++i) {
      ret = ret + Sparsity::band(n, i);
    }
    return ret;
  }

  Sparsity Sparsity::unit(int n, int el) {
    Sparsity ret = Sparsity::sparse(n);
    ret.getNZ(el, 0);
    return ret;
  }

  Sparsity Sparsity::rowcol(const std::vector<int>& row, const std::vector<int>& col,
                            int nrow, int ncol) {
    std::vector<int> all_rows, all_cols;
    all_rows.reserve(row.size()*col.size());
    all_cols.reserve(row.size()*col.size());
    for (std::vector<int>::const_iterator c_it=col.begin(); c_it!=col.end(); ++c_it) {
      casadi_assert_message(*c_it>=0 && *c_it<ncol, "Sparsity::rowcol: Column index out of bounds");
      for (std::vector<int>::const_iterator r_it=row.begin(); r_it!=row.end(); ++r_it) {
        casadi_assert_message(*r_it>=0 && *r_it<nrow, "Sparsity::rowcol: Row index out of bounds");
        all_rows.push_back(*r_it);
        all_cols.push_back(*c_it);
      }
    }
    return Sparsity::triplet(nrow, ncol, all_rows, all_cols);
  }

  Sparsity Sparsity::triplet(int nrow, int ncol, const std::vector<int>& row,
                             const std::vector<int>& col, std::vector<int>& mapping,
                             bool invert_mapping) {
    // Assert dimensions
    casadi_assert_message(col.size()==row.size(), "inconsistent lengths");

    // Create the return sparsity pattern and access vectors
    Sparsity ret = Sparsity::sparse(nrow, ncol);
    std::vector<int> &r_colind = ret.colindRef();
    std::vector<int> &r_row = ret.rowRef();
    r_row.reserve(row.size());

    // Consistency check and check if elements are already perfectly ordered with no duplicates
    int last_col=-1, last_row=-1;
    bool perfectly_ordered=true;
    for (int k=0; k<col.size(); ++k) {
      // Consistency check
      casadi_assert_message(col[k]>=0 && col[k]<ncol, "Col index out of bounds");
      casadi_assert_message(row[k]>=0 && row[k]<nrow, "Row index out of bounds");

      // Check if ordering is already perfect
      perfectly_ordered = perfectly_ordered && (col[k]<last_col ||
                                                (col[k]==last_col && row[k]<=last_row));
      last_col = col[k];
      last_row = row[k];
    }

    // Quick return if perfectly ordered
    if (perfectly_ordered) {
      // Save rows
      r_row.resize(row.size());
      copy(row.begin(), row.end(), r_row.begin());

      // Find offset index
      int el=0;
      for (int i=0; i<ncol; ++i) {
        while (el<col.size() && col[el]==i) el++;
        r_colind[i+1] = el;
      }

      // Identity mapping
      mapping.resize(row.size());
      for (int k=0; k<row.size(); ++k)
        mapping[k] = k;

      // Quick return
      return ret;
    }

    // Reuse data
    std::vector<int>& mapping1 = invert_mapping ? r_row : mapping;
    std::vector<int>& mapping2 = invert_mapping ? mapping : r_row;

    // Make sure that enough memory is allocated to use as a work vector
    mapping1.reserve(std::max(nrow+1, static_cast<int>(col.size())));

    // Number of elements in each row
    std::vector<int>& rowcount = mapping1; // reuse memory
    rowcount.resize(nrow+1);
    fill(rowcount.begin(), rowcount.end(), 0);
    for (std::vector<int>::const_iterator it=row.begin(); it!=row.end(); ++it) {
      rowcount[*it+1]++;
    }

    // Cumsum to get index offset for each row
    for (int i=0; i<nrow; ++i) {
      rowcount[i+1] += rowcount[i];
    }

    // New row for each old row
    mapping2.resize(row.size());
    for (int k=0; k<row.size(); ++k) {
      mapping2[rowcount[row[k]]++] = k;
    }

    // Number of elements in each col
    std::vector<int>& colcount = r_colind; // reuse memory, r_colind is already the right size
                                           // and is filled with zeros
    for (std::vector<int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it) {
      colcount[col[*it]+1]++;
    }

    // Cumsum to get index offset for each col
    for (int i=0; i<ncol; ++i) {
      colcount[i+1] += colcount[i];
    }

    // New col for each old col
    mapping1.resize(col.size());
    for (std::vector<int>::const_iterator it=mapping2.begin(); it!=mapping2.end(); ++it) {
      mapping1[colcount[col[*it]]++] = *it;
    }

    // Current element in the return matrix
    int r_el = 0;
    r_row.resize(col.size());

    // Current nonzero
    std::vector<int>::const_iterator it=mapping1.begin();

    // Loop over cols
    r_colind[0] = 0;
    for (int i=0; i<ncol; ++i) {

      // Previous row (to detect duplicates)
      int j_prev = -1;

      // Loop over nonzero elements of the col
      while (it!=mapping1.end() && col[*it]==i) {

        // Get the element
        int el = *it;
        it++;

        // Get the row
        int j = row[el];

        // If not a duplicate, save to return matrix
        if (j!=j_prev)
          r_row[r_el++] = j;

        if (invert_mapping) {
          // Save to the inverse mapping
          mapping2[el] = r_el-1;
        } else {
          // If not a duplicate, save to the mapping vector
          if (j!=j_prev)
            mapping1[r_el-1] = el;
        }

        // Save row
        j_prev = j;
      }

      // Update col offset
      r_colind[i+1] = r_el;
    }

    // Resize the row vector
    r_row.resize(r_el);

    // Resize mapping matrix
    if (!invert_mapping) {
      mapping1.resize(r_el);
    }

    return ret;
  }

  Sparsity Sparsity::triplet(int nrow, int ncol, const std::vector<int>& row,
                             const std::vector<int>& col) {
    std::vector<int> mapping;
    return Sparsity::triplet(nrow, ncol, row, col, mapping, false);
  }

  bool Sparsity::isSingular() const {
    casadi_assert_message(isSquare(), "isSingular: only defined for square matrices, but got "
                          << dimString());
    return rank(*this)!=size2();
  }

  std::vector<int> Sparsity::compress() const {
    // Get the sparsity pattern
    int nrow = this->size1();
    int ncol = this->size2();
    const vector<int>& colind = this->colind();
    const vector<int>& row = this->row();

    // Create compressed pattern
    vector<int> ret;
    ret.reserve(1 + 1 + colind.size() + row.size());
    ret.push_back(nrow);
    ret.push_back(ncol);
    ret.insert(ret.end(), colind.begin(), colind.end());
    ret.insert(ret.end(), row.begin(), row.end());
    return ret;
  }

  Sparsity Sparsity::compressed(const std::vector<int>& v) {
    // Check consistency
    casadi_assert(v.size() >= 2);
    //int nrow = v[0];
    int ncol = v[1];
    casadi_assert(v.size() >= 2 + ncol+1);
    int nnz = v[2 + ncol];
    casadi_assert(v.size() == 2 + ncol+1 + nnz);

    // Call array version
    return compressed(&v.front());
  }

  Sparsity Sparsity::compressed(const int* v) {
    // Get sparsity pattern
    int nrow = v[0];
    int ncol = v[1];
    const int *colind = v+2;
    int nnz = colind[ncol];
    const int *row = v + 2 + ncol+1;

    // Construct sparsity pattern
    return Sparsity(nrow, ncol, vector<int>(colind, colind+ncol+1), vector<int>(row, row+nnz));
  }

  void Sparsity::printCompact(std::ostream &stream) const {
    (*this)->printCompact(stream);
  }

  int Sparsity::bandwidthU() const {
    return (*this)->bandwidthU();
  }

  int Sparsity::bandwidthL() const {
    return (*this)->bandwidthL();
  }

} // namespace casadi
