/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#ifndef CASADI_STRUCT_HPP
#define CASADI_STRUCT_HPP

#include "im.hpp"
#include "dm.hpp"
#include "sx.hpp"
#include "mx.hpp"
#include "generic_type.hpp"
#include <functional>
#include <map>
#include <memory>

namespace casadi {

  /** \brief Layout of a flat vector as a tree of named, possibly repeated, matrix entries

      A path addresses a part of the layout.
      Its elements are, in order of descent:
      - an entry name, a list of names or ":" (all names) at a record
      - an index, a list of indices or a slice ":", "a:b", "a:b:c" at a repeated entry
      - one (linear) or two (row, column) matrix indices at a leaf

      The result of a path is an index matrix into the flat vector:
      - a leaf yields its (sparse) matrix shape
      - a record yields a column
      - several repetitions are concatenated horizontally
      - several names are concatenated vertically, as columns
      Omitted repetition indices select all repetitions.
      "@horzcat", "@vertcat", "@veccat" or "@blockcat" before a selection overrides how it is
      concatenated; "@blockcat" stacks it vertically and the next one horizontally.

      \author Joris Gillis
      \date 2026
  */
  class CASADI_EXPORT Struct
    : public SWIG_IF_ELSE(PrintableCommon, Printable<Struct>) {
  public:
    /// Empty structure
    Struct();

    /// Structure of scalar entries
    explicit Struct(const std::vector<std::string>& names);

    /// Matrix, optionally storing only the upper triangle of a symmetric pattern
    static Struct leaf(const Sparsity& sp, bool symmetric=false);

    ///@{
    /// Dense matrix with rows and columns labeled by (record) structures
    static Struct matrix(const Struct& rows, const Struct& cols, bool symmetric=false);
    static Struct matrix(const Struct& rows, casadi_int ncol=1);
    ///@}

    ///@{
    /// Add an entry, optionally repeated along the dimensions in repeat
    void add(const std::string& name, const Struct& s,
      const std::vector<casadi_int>& repeat=std::vector<casadi_int>());
    void add(const std::string& name, const Sparsity& sp,
      const std::vector<casadi_int>& repeat=std::vector<casadi_int>());
    void add(const std::string& name, casadi_int nrow=1, casadi_int ncol=1,
      const std::vector<casadi_int>& repeat=std::vector<casadi_int>());
    ///@}

    /// Interleave the repetitions of entries (along their first dimension) in the flat vector
    void interleave(const std::vector<std::string>& names);

    /// Length of the flat vector
    casadi_int nnz() const { return nnz_;}

    /// Indices into the flat vector addressed by a path
    Matrix<casadi_int> index(const std::vector<GenericType>& path, bool ind1=SWIG_IND1) const;

    /// Path of a single element of the flat vector
    std::vector<GenericType> path(casadi_int i, bool ind1=SWIG_IND1) const;

    /// Human-readable path of each element of the flat vector, 0-based like symbol names
    std::vector<std::string> labels() const;

    /// Is the structure a matrix, rather than a record?
    bool is_leaf() const { return is_leaf_;}

    /// Entry names, in declaration order
    std::vector<std::string> names() const;

    /// Repetition of an entry
    std::vector<casadi_int> repeat(const std::string& name) const;

    /// Structure of an entry
    Struct child(const std::string& name) const;

    /// Sparsity of a matrix
    const Sparsity& sparsity() const;

    /// Readable name of the class
    static std::string type_name() {return "Struct";}

    /// Print a description of the object
    void disp(std::ostream& stream, bool more=false) const;

    /// Get string representation
    std::string get_str(bool more=false) const {
      std::stringstream ss;
      disp(ss, more);
      return ss.str();
    }

#ifndef SWIG
    /// Matrix at an offset into the flat vector, with a stored sparsity
    template<class T> using Leaf = std::function<T(casadi_int, const Sparsity&)>;

    /// Part addressed by a path, assembled from matrices
    template<class T>
    T get(const std::vector<GenericType>& path, bool ind1, const Leaf<T>& leaf) const;

    /// Visit each matrix in flat order: offset, stored sparsity and path
    typedef std::function<void(casadi_int, const Sparsity&, const std::vector<GenericType>&)>
      LeafVisitor;
    void for_each_leaf(const LeafVisitor& f) const;

    /// Path element for a slice
    static std::string slice_str(const Slice& s);

  private:
    struct Entry {
      std::string name;
      std::vector<casadi_int> repeat;
      std::shared_ptr<const Struct> s;
      // Repetitions along the first dimension: count, size and offset of each
      casadi_int n, size;
      std::vector<casadi_int> start;
    };
    // Matrix
    bool is_leaf_, symmetric_;
    Sparsity sp_;
    std::shared_ptr<const Struct> rows_, cols_;
    // Record
    std::vector<Entry> entries_;
    std::vector<std::vector<casadi_int> > groups_;
    // Flat order: start, entry and first repetition index of each block
    std::vector<casadi_int> start_, block_entry_, block_rep_;
    casadi_int nnz_;

    void update();
    casadi_int find(const std::string& name) const;
    template<class T>
    T get(casadi_int offset, const std::vector<GenericType>& p, casadi_int k, bool ind1,
      const Leaf<T>& leaf, int cat) const;
    template<class T>
    T get(casadi_int e, casadi_int offset, const std::vector<GenericType>& p, casadi_int k,
      std::vector<casadi_int>& rep, bool ind1, const Leaf<T>& leaf, int cat) const;
    template<class T>
    T get_leaf(casadi_int offset, const std::vector<GenericType>& p, casadi_int k, bool ind1,
      const Leaf<T>& leaf) const;
    Sparsity stored() const;
    // Indices along a matrix dimension, resolving names against labels
    static IM dim_index(const GenericType& key, const std::shared_ptr<const Struct>& labels,
      casadi_int n, bool ind1);
    void path(casadi_int i, std::vector<GenericType>& p) const;
    void for_each_leaf(casadi_int offset, std::vector<GenericType>& p, const LeafVisitor& f) const;
#endif // SWIG
  };

  /** \brief Flat vector with a Struct layout

      \author Joris Gillis
      \date 2026
  */
  template<class M>
  class CASADI_EXPORT StructValue
    : public SWIG_IF_ELSE(PrintableCommon, Printable<StructValue<M> >) {
  public:
    /// Wrap a flat vector, or fill with a scalar
    StructValue(const Struct& s, const M& data);

    /// Symbolic, one primitive per matrix, named by prefix and path; only accepts symbols
    static StructValue sym(const Struct& s, const std::string& prefix="");

#if !defined(SWIG) || defined(SWIGPYTHON)
    /// Read and write a dense matrix in place, which must outlive the result
    static StructValue view(const Struct& s, M* target);
#endif

    /// Layout
    Struct structure() const { return s_;}

    /// Flat vector
    const M& cat() const;

    /// Get the part addressed by a path
    M get(const std::vector<GenericType>& path, bool ind1=SWIG_IND1) const;

    /** \brief Set the part addressed by a path

        Scalars are repeated, rows repeated vertically, matrices repeated horizontally
        and vectors transposed as needed.
    */
    void set(const std::vector<GenericType>& path, const M& value, bool ind1=SWIG_IND1);

    /// Get flat vector entries, shaped like ind
    M get_nz(const Matrix<casadi_int>& ind, bool ind1=SWIG_IND1) const;

    /// Set flat vector entries, as for set
    void set_nz(const Matrix<casadi_int>& ind, const M& value, bool ind1=SWIG_IND1);

    /// Readable name of the class
    static std::string type_name();

    /// Print a description of the object
    void disp(std::ostream& stream, bool more=false) const;

    /// Get string representation
    std::string get_str(bool more=false) const {
      std::stringstream ss;
      disp(ss, more);
      return ss.str();
    }

#ifndef SWIG
    /// Part addressed by a path, assignable
    class Ref : public M {
    public:
      Ref(StructValue& v, const std::vector<GenericType>& p) : M(v.get(p)), v_(v), p_(p) {}
      const M& operator=(const M& y) { v_.set(p_, y); return y;}
      const M& operator=(const Ref& y) { return operator=(static_cast<const M&>(y));}
    private:
      StructValue& v_;
      std::vector<GenericType> p_;
    };

    template<typename... A> Ref operator()(const A&... a) { return Ref(*this, keys(a...));}
    template<typename... A> M operator()(const A&... a) const { return get(keys(a...));}

    /// Conversion to the flat vector
    explicit operator const M&() const { return cat();}

  private:
    template<typename T> static GenericType key(const T& a) { return a;}
    static GenericType key(const Slice& a) { return Struct::slice_str(a);}
    // Not {key(a)...}: a single key may bind via GenericType's conversion to a vector
    template<typename... A> static std::vector<GenericType> keys(const A&... a) {
      return std::vector<GenericType>(std::initializer_list<GenericType>{key(a)...});
    }
    Struct s_;
    // Stored nonzeros of each matrix, and their offset in the flat vector
    std::vector<M> leaves_;
    std::vector<casadi_int> offset_;
    mutable M cat_;
    mutable bool cached_;
    // Viewed matrix, if any
    M* target_;
    // Only symbols may be assigned
    bool symbolic_;
    // Values assigned to whole matrices at once, by first flat index
    std::map<casadi_int, std::pair<IM, M> > assigned_;
    casadi_int leaf(casadi_int offset) const;
    M leaf_value(casadi_int offset, const Sparsity& sp) const;
    void assign(const IM& ind, const M& value, const std::string& what);
#endif // SWIG
  };

  typedef StructValue<DM> StructDM;
  typedef StructValue<SX> StructSX;
  typedef StructValue<MX> StructMX;

#ifndef CASADI_STRUCT_CPP
  extern template class StructValue<DM>;
  extern template class StructValue<SX>;
  extern template class StructValue<MX>;
#endif // CASADI_STRUCT_CPP

} // namespace casadi

#endif // CASADI_STRUCT_HPP
