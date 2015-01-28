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

  Sparsity::Sparsity(int nrow, int ncol) {
    std::vector<int> row, colind(ncol+1, 0);
    assignCached(nrow, ncol, colind, row);
    sanityCheck(true);
  }

  Sparsity::Sparsity(const std::pair<int, int>& rc) {
    std::vector<int> row, colind(rc.second+1, 0);
    assignCached(rc.first, rc.second, colind, row);
    sanityCheck(true);
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
    return (*this)->size1();
  }

  int Sparsity::size2() const {
    return (*this)->size2();
  }

  int Sparsity::numel() const {
    return (*this)->numel();
  }

  bool Sparsity::isEmpty(bool both) const {
    return (*this)->isEmpty(both);
  }

  int Sparsity::nnz() const {
    return (*this)->nnz();
  }

  std::pair<int, int> Sparsity::shape() const {
    return (*this)->shape();
  }

  const int* Sparsity::row() const {
    return (*this)->row();
  }

  const int* Sparsity::colind() const {
    return (*this)->colind();
  }

  int Sparsity::row(int el) const {
    if (el<0 || el>=nnz()) {
      std::stringstream ss;
      ss <<  "Sparsity::row: Index " << el << " out of range [0," << nnz() << ")";
      throw std::out_of_range(ss.str());
    }
    return row()[el];
  }

  int Sparsity::colind(int cc) const {
    if (cc<0 || cc>size2()) {
      std::stringstream ss;
      ss << "Sparsity::colind: Index " << cc << " out of range [0," << size2() << "]";
      throw std::out_of_range(ss.str());
    }
    return colind()[cc];
  }

  void Sparsity::sanityCheck(bool complete) const {
    (*this)->sanityCheck(complete);
  }

  void Sparsity::resize(int nrow, int ncol) {
    if (size1()!=nrow || size2() != ncol) {
      *this = (*this)->zz_resize(nrow, ncol);
    }
  }

  int Sparsity::addNZ(int rr, int cc) {
    // If negative index, count from the back
    if (rr<0) rr += size1();
    if (cc<0) cc += size2();

    // Check consistency
    casadi_assert_message(rr>=0 && rr<size1(), "Row index out of bounds");
    casadi_assert_message(cc>=0 && cc<size2(), "Column index out of bounds");

    // Quick return if matrix is dense
    if (isDense()) return rr+cc*size1();

    // Quick return if we are adding an element to the end
    if (colind(cc)==nnz() || (colind(cc+1)==nnz() && row(nnz()-1)<rr)) {
      std::vector<int> rowv = getRow();
      std::vector<int> colindv = getColind();
      rowv.push_back(rr);
      for (int c=cc; c<size2(); ++c) {
        colindv[c+1]++;
      }
      assignCached(size1(), size2(), colindv, rowv);
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
    std::vector<int> rowv = getRow();
    std::vector<int> colindv = getColind();
    rowv.insert(rowv.begin()+ind, rr);
    for (int c=cc+1; c<size2()+1; ++c)
      colindv[c]++;

    // Return the location of the new element
    assignCached(size1(), size2(), colindv, rowv);
    return ind;
  }

  bool Sparsity::hasNZ(int rr, int cc) const {
    return getNZ(rr, cc)!=-1;
  }


  int Sparsity::getNZ(int rr, int cc) const {
    return (*this)->getNZ(rr, cc);
  }

  Sparsity Sparsity::zz_reshape(const Sparsity& sp) const {
    casadi_assert(isReshape(sp));
    return sp;
  }

  Sparsity Sparsity::zz_reshape(int nrow, int ncol) const {
    return (*this)->zz_reshape(nrow, ncol);
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

  Sparsity Sparsity::sub(const std::vector<int>& rr, const Sparsity& sp,
                         std::vector<int>& mapping, bool ind1) const {
    return (*this)->sub(rr, *sp, mapping, ind1);
  }

  Sparsity Sparsity::sub(const std::vector<int>& rr, const std::vector<int>& cc,
                         std::vector<int>& mapping, bool ind1) const {
    return (*this)->sub(rr, cc, mapping, ind1);
  }

  std::vector<int> Sparsity::erase(const std::vector<int>& rr, const std::vector<int>& cc,
                                   bool ind1) {
    vector<int> mapping;
    *this = (*this)->zz_erase(rr, cc, ind1, mapping);
    return mapping;
  }

  std::vector<int> Sparsity::erase(const std::vector<int>& rr, bool ind1) {
    vector<int> mapping;
    *this = (*this)->zz_erase(rr, ind1, mapping);
    return mapping;
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

  std::vector<int> Sparsity::getColind() const {
    return (*this)->getColind();
  }

  std::vector<int> Sparsity::getCol() const {
    return (*this)->getCol();
  }

  std::vector<int> Sparsity::getRow() const {
    return (*this)->getRow();
  }

  void Sparsity::getCCS(std::vector<int>& cind, std::vector<int>& r) const {
    cind.resize(size2()+1);
    const int* colind = this->colind();
    copy(colind, colind+cind.size(), cind.begin());
    r.resize(nnz());
    const int* row = this->row();
    copy(row, row+r.size(), r.begin());
  }

  void Sparsity::getCRS(std::vector<int>& rind, std::vector<int>& c) const {
    T().getCCS(rind, c);
  }


  void Sparsity::getTriplet(std::vector<int>& r, std::vector<int>& c) const {
    r.resize(nnz());
    const int* row = this->row();
    copy(row, row+r.size(), r.begin());
    c = this->getCol();
  }

  Sparsity Sparsity::transpose(std::vector<int>& mapping, bool invert_mapping) const {
    return (*this)->transpose(mapping, invert_mapping);
  }

  Sparsity Sparsity::T() const {
    return (*this)->T();
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

  Sparsity Sparsity::patternProduct(const Sparsity& y) const {
    return (*this)->patternProduct(y);
  }

  bool Sparsity::isEqual(const Sparsity& y) const {
    return (*this)->isEqual(y);
  }

  bool Sparsity::isEqual(int nrow, int ncol, const std::vector<int>& colind,
                         const std::vector<int>& row) const {
    return (*this)->isEqual(nrow, ncol, colind, row);
  }

  bool Sparsity::isEqual(int nrow, int ncol, const int* colind, const int* row) const {
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
        // Append to vector (inefficient)
        *this = (*this)->zz_appendVector(*sp);
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
        // Append to matrix (expensive)
        *this = (*this)->zz_appendColumns(*sp);
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

  void Sparsity::enlarge(int nrow, int ncol, const std::vector<int>& rr,
                         const std::vector<int>& cc, bool ind1) {
    enlargeColumns(ncol, cc, ind1);
    enlargeRows(nrow, rr, ind1);
  }

  void Sparsity::enlargeColumns(int ncol, const std::vector<int>& cc, bool ind1) {
    casadi_assert(cc.size() == size2());
    if (cc.empty()) {
      *this = Sparsity(size1(), ncol);
    } else {
      *this = (*this)->zz_enlargeColumns(ncol, cc, ind1);
    }
  }

  void Sparsity::enlargeRows(int nrow, const std::vector<int>& rr, bool ind1) {
    casadi_assert(rr.size() == size1());
    if (rr.empty()) {
      *this = Sparsity(nrow, size2());
    } else {
      *this = (*this)->zz_enlargeRows(nrow, rr, ind1);
    }
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

  std::vector<int> Sparsity::find(bool ind1) const {
    std::vector<int> loc;
    find(loc, ind1);
    return loc;
  }

  void Sparsity::find(std::vector<int>& loc, bool ind1) const {
    (*this)->find(loc, ind1);
  }

  void Sparsity::getNZ(std::vector<int>& indices) const {
    (*this)->getNZ(indices);
  }

  Sparsity Sparsity::unidirectionalColoring(const Sparsity& AT, int cutoff) const {
    if (AT.isNull()) {
      return (*this)->unidirectionalColoring(T(), cutoff);
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
    casadi_assert(colind.size()==ncol+1);
    casadi_assert(row.size()==colind.back());
    assignCached(nrow, ncol, getPtr(colind), getPtr(row));
  }

  void Sparsity::assignCached(int nrow, int ncol, const int* colind, const int* row) {
    // Scalars and empty patterns are handled separately
    if (ncol==0 && nrow==0) {
      // If empty
      *this = getEmpty();
      return;
    } else if (ncol==1 && nrow==1) {
      if (colind[ncol]==0) {
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

  Sparsity Sparsity::zz_tril(bool includeDiagonal) const {
    return (*this)->zz_tril(includeDiagonal);
  }

  Sparsity Sparsity::zz_triu(bool includeDiagonal) const {
    return (*this)->zz_triu(includeDiagonal);
  }

  std::vector<int> Sparsity::getLowerNZ() const {
    return (*this)->getLowerNZ();
  }

  std::vector<int> Sparsity::getUpperNZ() const {
    return (*this)->getUpperNZ();
  }


  std::size_t hash_sparsity(int nrow, int ncol, const std::vector<int>& colind,
                            const std::vector<int>& row) {
    return hash_sparsity(nrow, ncol, getPtr(colind), getPtr(row));
  }

  std::size_t hash_sparsity(int nrow, int ncol, const int* colind, const int* row) {
    // Condense the sparsity pattern to a single, deterministric number
    std::size_t ret=0;
    hash_combine(ret, nrow);
    hash_combine(ret, ncol);
    hash_combine(ret, colind, ncol+1);
    hash_combine(ret, row, colind[ncol]);
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
    return Sparsity(nrow, ncol);
  }

  Sparsity Sparsity::upper(int n) {
    casadi_assert_message(n>=0, "Sparsity::upper expects a positive integer as argument");
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

  Sparsity Sparsity::lower(int n) {
    casadi_assert_message(n>=0, "Sparsity::lower expects a positive integer as argument");
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
    Sparsity ret = Sparsity(n, n);
    for (int i=-p;i<=p;++i) {
      ret = ret + Sparsity::band(n, i);
    }
    return ret;
  }

  Sparsity Sparsity::unit(int n, int el) {
    std::vector<int> row(1, el), colind(2);
    colind[0] = 0;
    colind[1] = 1;
    return Sparsity(n, 1, colind, row);
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
    std::vector<int> r_colind(ncol+1, 0);
    std::vector<int> r_row;
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
      return Sparsity(nrow, ncol, r_colind, r_row);
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

    return Sparsity(nrow, ncol, r_colind, r_row);
  }

  Sparsity Sparsity::triplet(int nrow, int ncol, const std::vector<int>& row,
                             const std::vector<int>& col) {
    std::vector<int> mapping;
    return Sparsity::triplet(nrow, ncol, row, col, mapping, false);
  }

  bool Sparsity::isSingular() const {
    casadi_assert_message(isSquare(), "isSingular: only defined for square matrices, but got "
                          << dimString());
    return sprank(*this)!=size2();
  }

  std::vector<int> Sparsity::compress() const {
    // Get the sparsity pattern
    int nrow = this->size1();
    int ncol = this->size2();
    int sz = this->nnz();
    const int* colind = this->colind();
    const int* row = this->row();

    // Create compressed pattern
    vector<int> ret;
    ret.reserve(2 + ncol + 1 + sz);
    ret.push_back(nrow);
    ret.push_back(ncol);
    ret.insert(ret.end(), colind, colind+ncol+1);
    ret.insert(ret.end(), row, row + sz);
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

  Sparsity Sparsity::zz_horzcat(const std::vector<Sparsity> & sp) {
    if (sp.empty()) {
      return Sparsity();
    } else {
      Sparsity ret = sp[0];
      for (int i=1; i<sp.size(); ++i) {
        ret.appendColumns(sp[i]);
      }
      return ret;
    }
  }

  Sparsity Sparsity::zz_vertcat(const std::vector<Sparsity> & sp) {
    if (sp.empty()) {
      return Sparsity();
    } else if (sp[0].isVector()) {
      Sparsity ret = sp[0];
      for (int i=1; i<sp.size(); ++i) {
        ret.append(sp[i]);
      }
      return ret;
    } else {
      Sparsity ret = sp[0].T();
      for (int i=1; i<sp.size(); ++i) {
        ret.appendColumns(sp[i].T());
      }
      return ret.T();
    }
  }

  Sparsity Sparsity::zz_diagcat(const std::vector< Sparsity > &v) {
    int n = 0;
    int m = 0;

    std::vector<int> colind(1, 0);
    std::vector<int> row;

    int nz = 0;
    for (int i=0;i<v.size();++i) {
      const int* colind_ = v[i].colind();
      int ncol = v[i].size2();
      const int* row_ = v[i].row();
      int sz = v[i].nnz();
      for (int k=1; k<ncol+1; ++k) {
        colind.push_back(colind_[k]+nz);
      }
      for (int k=0; k<sz; ++k) {
        row.push_back(row_[k]+m);
      }
      n+= v[i].size2();
      m+= v[i].size1();
      nz+= v[i].nnz();
    }

    return Sparsity(m, n, colind, row);
  }

  std::vector<Sparsity> Sparsity::zz_horzsplit(const std::vector<int>& offset) const {
    // Consistency check
    casadi_assert(offset.size()>=1);
    casadi_assert(offset.front()==0);
    casadi_assert_message(offset.back()==size2(),
                          "horzsplit(Sparsity, std::vector<int>): Last elements of offset "
                          "(" << offset.back() << ") must equal the number of columns "
                          "(" << size2() << ")");
    casadi_assert(isMonotone(offset));

    // Number of outputs
    int n = offset.size()-1;

    // Get the sparsity of the input
    const int* colind_x = colind();
    const int* row_x = row();

    // Allocate result
    std::vector<Sparsity> ret;
    ret.reserve(n);

    // Sparsity pattern as CCS vectors
    vector<int> colind, row;
    int ncol, nrow = size1();

    // Get the sparsity patterns of the outputs
    for (int i=0; i<n; ++i) {
      int first_col = offset[i];
      int last_col = offset[i+1];
      ncol = last_col - first_col;

      // Construct the sparsity pattern
      colind.resize(ncol+1);
      copy(colind_x+first_col, colind_x+last_col+1, colind.begin());
      for (vector<int>::iterator it=colind.begin()+1; it!=colind.end(); ++it) *it -= colind[0];
      colind[0] = 0;
      row.resize(colind.back());
      copy(row_x+colind_x[first_col], row_x+colind_x[last_col], row.begin());

      // Append to the list
      ret.push_back(Sparsity(nrow, ncol, colind, row));
    }

    // Return (RVO)
    return ret;
  }

  std::vector<Sparsity> Sparsity::zz_vertsplit(const std::vector<int>& offset) const {
    std::vector<Sparsity> ret = horzsplit(T(), offset);
    for (std::vector<Sparsity>::iterator it=ret.begin(); it!=ret.end(); ++it) {
      *it = it->T();
    }
    return ret;
  }

  Sparsity Sparsity::zz_blockcat(const std::vector< std::vector< Sparsity > > &v) {
    std::vector< Sparsity > ret;
    for (int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }

  Sparsity Sparsity::zz_vecNZ() const {
    return Sparsity::dense(nnz());
  }

  std::vector<Sparsity> Sparsity::zz_diagsplit(const std::vector<int>& offset1,
                                               const std::vector<int>& offset2) const {
    // Consistency check
    casadi_assert(offset1.size()>=1);
    casadi_assert(offset1.front()==0);
    casadi_assert_message(offset1.back()==size1(),
                          "diagsplit(Sparsity, offset1, offset2): Last elements of offset1 "
                          "(" << offset1.back() << ") must equal the number of rows "
                          "(" << size1() << ")");
    casadi_assert_message(offset2.back()==size2(),
                          "diagsplit(Sparsity, offset1, offset2): Last elements of offset2 "
                          "(" << offset2.back() << ") must equal the number of rows "
                          "(" << size2() << ")");
    casadi_assert(isMonotone(offset1));
    casadi_assert(isMonotone(offset2));
    casadi_assert(offset1.size()==offset2.size());

    // Number of outputs
    int n = offset1.size()-1;

    // Return value
    std::vector<Sparsity> ret;

    // Caveat: this is a very silly implementation
    IMatrix x = IMatrix::zeros(*this);

    for (int i=0;i<n;++i) {
      ret.push_back(x(Slice(offset1[i], offset1[i+1]),
                      Slice(offset2[i], offset2[i+1])).sparsity());
    }

    return ret;
  }

  int Sparsity::zz_sprank() const {
    std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
    return coarse_colblock.at(3);
  }

} // namespace casadi
