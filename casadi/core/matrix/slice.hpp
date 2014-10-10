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


#ifndef CASADI_SLICE_HPP
#define CASADI_SLICE_HPP

#include <vector>
#include "../casadi_exception.hpp"
#include "../printable_object.hpp"
#include <limits>
#include <iostream>

namespace casadi {

  /** \brief Class representing a Slice
   *
   * Note that Python or Octave do not need to use this class.
   * They can just use slicing utility from the host language ( M[0:6]  in Python, M(1:7) )
   */
  class CASADI_CORE_EXPORT Slice : public PrintableObject<Slice> {
  public:
    /// Default constructor - all elements
    Slice();

    /// A single element
    Slice(int i);

    /// A slice
    Slice(int start, int stop, int step=1);

    /// Construct from an index vector (requires isSlice(v) to be true)
    explicit Slice(const std::vector<int>& v);

    /// Construct nested slices from an index vector (requires isSlice2(v) to be true)
    explicit Slice(const std::vector<int>& v, Slice& outer);

    /// Check if an index vector can be represented more efficiently as a slice
    static bool isSlice(const std::vector<int>& v);

    /// Check if an index vector can be represented more efficiently as two nested slices
    static bool isSlice2(const std::vector<int>& v);

    /// Get a vector of indices
    std::vector<int> getAll(int len) const;

    /// Get a vector of indices (nested slice)
    std::vector<int> getAll(const Slice& outer, int len) const;

    /// Check equality
    bool operator==(const Slice& other) const
    { return start_==other.start_ && stop_==other.stop_ && step_==other.step_;}

    /// Check inequality
    bool operator!=(const Slice& other) const { return !(*this == other);}

    /// Print a representation of the object
    void repr(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print a description of the object
    void print(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// start value: negative values will get added to length
    int start_;
    /// stop value: use std::numeric_limits<int>::max() to indicate unboundedness
    int stop_;
    int step_;
  };

#ifndef SWIG
  static Slice ALL;
#endif // SWIG

  /// \cond INTERNAL
  /**  Class representing a non-regular (and thus non-slice) index list
   */
  class CASADI_CORE_EXPORT IndexList {
  private:

  public:
    enum Type {NILL, INT, SLICE, IVECTOR};
    /// Constructor
    IndexList();
    explicit IndexList(int i);
    explicit IndexList(const std::vector<int> &i);
    explicit IndexList(const Slice &i);

    /// Get a vector of indices
    std::vector<int> getAll(int len) const;

    /// Data members (all public)
    Slice slice;
    int i;
    std::vector<int> iv;
    Type type;
  };
  /// \endcond

} // namespace casadi


#ifdef SWIG
%template(Pair_Slice_Int) std::pair<casadi::Slice, int>;
%template(Pair_Int_Slice) std::pair<int, casadi::Slice>;
%template(Pair_Slice_Slice) std::pair<casadi::Slice, casadi::Slice>;
#endif // SWIG

#endif // CASADI_SLICE_HPP
