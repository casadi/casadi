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


#ifndef CASADI_MISC_HPP
#define CASADI_MISC_HPP

#include "exception.hpp"
#include "casadi_common.hpp"

/** \brief Convenience tools for C++ Standard Library vectors

    \author Joel Andersson
    \date 2010-2011

    \identifier{1l5} */

namespace casadi {

#ifndef SWIG

template<typename T>
class scoped_checkout {
public:
  scoped_checkout(const T& proto) : proto_(proto) {
    mem = proto_.checkout();
  }

  scoped_checkout(scoped_checkout&& that) : mem(that.mem), proto_(that.proto_) {
    that.mem = -1;
  }

  scoped_checkout(const scoped_checkout& that) = delete;

  ~scoped_checkout() {
    if (mem!=-1) proto_.release(mem);
  }

  operator casadi_int() const {
    return mem;
  }

private:
  int mem;
  const T& proto_;
};

  /**  \brief Range function

  * \param start
  * \param stop
  * \param step
  * \param len
  *
  * Consider a infinitely long list [start, start+step, start+2*step, ...]
  * Elements larger than or equal to stop are chopped off.
  *
       \identifier{1l6} */
  CASADI_EXPORT std::vector<casadi_int> range(casadi_int start, casadi_int stop, casadi_int step=1,
                                            casadi_int len=std::numeric_limits<casadi_int>::max());

  /** \brief Check if a vector matches a range
   *
      \identifier{1l7} */
  CASADI_EXPORT bool is_range(const std::vector<casadi_int>& v,
    casadi_int start, casadi_int stop, casadi_int step=1);

  CASADI_EXPORT std::string join(const std::vector<std::string>& l, const std::string& delim=",");

  /// Checks if s starts with p
  CASADI_EXPORT bool startswith(const std::string& s, const std::string& p);

  /// Replace all occurences of p with r in s
  CASADI_EXPORT std::string replace(const std::string& s,
    const std::string& p, const std::string& r);

  /**  \brief Range function

  * \param stop
  *
  * \return list [0, 1, 2...stop-1]

       \identifier{1l8} */
  CASADI_EXPORT std::vector<casadi_int> range(casadi_int stop);

  /// Check if all arguments are true
  CASADI_EXPORT bool all(const std::vector<bool> &v);
  /// Check if any arguments are true
  CASADI_EXPORT bool any(const std::vector<bool> &v);
  /// Invert all entries
  CASADI_EXPORT std::vector<bool> boolvec_not(const std::vector<bool> &v);
  /// And operation on boolean vector
  CASADI_EXPORT std::vector<bool> boolvec_and(const std::vector<bool> &lhs,
      const std::vector<bool> &rhs);
  /// Or operation on boolean vector
  CASADI_EXPORT std::vector<bool> boolvec_or(const std::vector<bool> &lhs,
      const std::vector<bool> &rhs);

  CASADI_EXPORT std::vector<casadi_int> boolvec_to_index(const std::vector<bool> &v);

  CASADI_EXPORT bool is_equally_spaced(const std::vector<double> &v);

  /// Computes a mapping for a (dense) tensor permutation
  CASADI_EXPORT std::vector<casadi_int> tensor_permute_mapping(const std::vector<casadi_int>& dims,
      const std::vector<casadi_int>& order);

  CASADI_EXPORT int to_int(casadi_int rhs);
  CASADI_EXPORT std::vector<int> to_int(const std::vector<casadi_int>& rhs);
  CASADI_EXPORT std::vector< std::vector<int> > to_int(
    const std::vector< std::vector<casadi_int> >& rhs);

  template<typename T, typename S>
  std::vector<T> vector_static_cast(const std::vector<S>& rhs);

  CASADI_EXPORT std::string str_bvec(bvec_t v);

  /**  \brief Slicing vector

  *  \param v Vector to slice
  *  \param i List of indices

       \identifier{1l9} */
  template<typename T>
  std::vector<T> vector_slice(const std::vector<T> &v, const std::vector<casadi_int> &i);

  /** \brief Reverse a list

      \identifier{1la} */
  template<typename T>
  std::vector<T> reverse(const std::vector<T> &v);

  /** \brief Join two lists

      \identifier{1lb} */
  template<typename T>
  std::vector<T> join(const std::vector<T> &a, const std::vector<T> &b);

  /** \brief Join three lists

      \identifier{1lc} */
  template<typename T>
  std::vector<T> join(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c);

  /** \brief Does the list represent a permutation?

      \identifier{250} */
  bool is_permutation(const std::vector<casadi_int> &order);

  /** \brief permute a list

      \identifier{1ld} */
  template<typename T>
  std::vector<T> permute(const std::vector<T> &a, const std::vector<casadi_int> &order);

  /** \brief inverse a permutation vector

      \identifier{251} */
  std::vector<casadi_int> invert_permutation(const std::vector<casadi_int> &a);

  /** \brief find nonzeros

      \identifier{1le} */
  template<typename T>
  std::vector<casadi_int> find(const std::vector<T> &v);

  #endif // SWIG

  /// Check if for each element of v holds: v_i < upper
  template<typename T>
  bool in_range(const std::vector<T> &v, casadi_int upper);

  /// Check if for each element of v holds: lower <= v_i < upper
  template<typename T>
  bool in_range(const std::vector<T> &v, casadi_int lower, casadi_int upper);

  // Assert that a indices are in a range
  #define casadi_assert_in_range(v, lower, upper) \
   casadi_assert(in_range(v, lower, upper), \
    "Out of bounds error. Got elements in range [" \
    + str(*std::min_element(v.begin(), v.end())) + ","\
    + str(*std::max_element(v.begin(), v.end())) + "], which is outside the range ["\
    + str(lower) + "," + str(upper) + ").")

    // Assert that a indices are bounded
    #define casadi_assert_bounded(v, upper) \
     casadi_assert(in_range(v, upper), \
      "Out of bounds error. Got elements in range [" \
      + str(*std::min_element(v.begin(), v.end())) + ","\
      + str(*std::max_element(v.begin(), v.end())) + "], which exceeds the upper bound "\
      + str(upper) + ".")

  /** \brief Returns the list of all i in [0, size[ not found in supplied list
  *
  * The supplied vector may contain duplicates and may be non-monotonous
  * The supplied vector will be checked for bounds
  * The result vector is guaranteed to be monotonously increasing

      \identifier{1lf} */
  CASADI_EXPORT std::vector<casadi_int> complement(const std::vector<casadi_int> &v,
                                                    casadi_int size);

  /** \brief Returns a vector for quickly looking up entries of supplied list
  *
  *  lookupvector[i]!=-1     <=>  v contains i
  *  v[lookupvector[i]] == i <=>  v contains i
  *
  *  Duplicates are treated by looking up last occurrence

      \identifier{1lg} */
  CASADI_EXPORT std::vector<casadi_int> lookupvector(const std::vector<casadi_int> &v,
                                                     casadi_int size);
  CASADI_EXPORT std::vector<casadi_int> lookupvector(const std::vector<casadi_int> &v);

  /** \brief Flatten a nested std::vector tot a single flattened vector
   * 
   * Contents of nested[i] ends up in flat[indices[i]]..flat[indices[i+1]-1]

      \identifier{1lh} */
  template<class T, class S>
  void flatten_nested_vector(const std::vector< std::vector<T> >& nested,
                            std::vector<S>& flat);

  /** \brief Flatten a nested std::vector tot a single flattened vector
   * 
   * Contents of nested[i] ends up in flat[indices[i]]..flat[indices[i+1]-1]

      \identifier{1li} */
  template<class T, class S, class I>
  void flatten_nested_vector(const std::vector< std::vector<T> >& nested,
                            std::vector<S>& flat,
                            std::vector<I>& indices);

  /// \cond INTERNAL
#ifndef SWIG
  /**
  Apply a function f to each element in a vector
  */
  template<class T>
  std::vector<T> applymap(T (*f)(const T&), const std::vector<T>& comp);

  /**
  Apply a function f to each element in a vector
  */
  template<class T>
  void applymap(void (*f)(T&), std::vector<T>& comp);
#endif // SWIG
  /// \endcond


  /// Check if the vector is strictly increasing
  template<typename T>
  bool is_increasing(const std::vector<T> &v);

  /// Check if the vector is strictly decreasing
  template<typename T>
  bool is_decreasing(const std::vector<T> &v);

  /// Check if the vector is non-increasing
  template<typename T>
  bool is_nonincreasing(const std::vector<T> &v);

  /// Check if the vector is non-decreasing
  template<typename T>
  bool is_nondecreasing(const std::vector<T> &v);

  /// Check if the vector is monotone
  template<typename T>
  bool is_monotone(const std::vector<T> &v);

  /// Check if the vector is strictly monotone
  template<typename T>
  bool is_strictly_monotone(const std::vector<T> &v);

  /// Check if the vector has negative entries
  template<typename T>
  bool has_negative(const std::vector<T> &v);

  /// Print vector, matlab style
  template<typename T>
  void write_matlab(std::ostream &stream, const std::vector<T> &v);

  /// Print matrix, matlab style
  template<typename T>
  void write_matlab(std::ostream &stream, const std::vector<std::vector<T> > &v);

  /// Read vector, matlab style
  template<typename T>
  void read_matlab(std::istream &stream, std::vector<T> &v);

  /// Read matrix, matlab style
  template<typename T>
  void read_matlab(std::ifstream &file, std::vector<std::vector<T> > &v);

#ifndef SWIG
  /// Matlab's linspace
  template<typename T, typename F, typename L>
  void linspace(std::vector<T> &v, const F& first, const L& last);

  /// \cond INTERNAL
  /// Get an pointer of sets of booleans from a double vector
  CASADI_EXPORT bvec_t* get_bvec_t(std::vector<double>& v);

  /// Get an pointer of sets of booleans from a double vector
  CASADI_EXPORT const bvec_t* get_bvec_t(const std::vector<double>& v);

  /// Bit-wise or operation on bvec_t array
  CASADI_EXPORT bvec_t bvec_or(const bvec_t* arg, casadi_int n);

  /// Get an pointer of sets of booleans from a double vector
  template<typename T>
  bvec_t* get_bvec_t(std::vector<T>& v);

  /// Get an pointer of sets of booleans from a double vector
  template<typename T>
  const bvec_t* get_bvec_t(const std::vector<T>& v);

  /// Get a pointer to the data contained in the vector
  template<typename T>
  T* get_ptr(std::vector<T> &v);

  /// Get a pointer to the data contained in the vector
  template<typename T>
  const T* get_ptr(const std::vector<T> &v);

  /// \endcond

  /** \brief Sort the data in a vector
  *
  * \param[in]  values the vector that needs sorting
  * \param[out] sorted_values the sorted vector
  * \param[out] indices The indices such that 'sorted_values= values[indices]'
  * \param[in] invert_indices Output indices such that 'sorted_values[indices=values'
  *
      \identifier{1lj} */
  template<typename T>
  void sort(const std::vector<T> &values, std::vector<T> &sorted_values,
            std::vector<casadi_int> &indices, bool invert_indices =false);

  /** \brief product
  *
      \identifier{1lk} */
  template<typename T>
  T product(const std::vector<T> &values);

  /** \brief sum
  *
      \identifier{1ll} */
  template<typename T>
  T sum(const std::vector<T> &values);

  /** \brief cumulative sum
  *
      \identifier{1lm} */
  template<typename T>
  std::vector<T> cumsum(const std::vector<T> &values);

  /** \brief diff
  *
      \identifier{1ln} */
  template<typename T>
  std::vector<T> diff(const std::vector<T> &values);

  /** \brief cumulative sum, starting with zero
  *
      \identifier{1lo} */
  template<typename T>
  std::vector<T> cumsum0(const std::vector<T> &values);
#endif //SWIG

  /// Checks if array does not contain NaN or Inf
  template<typename T>
  bool is_regular(const std::vector<T> &v) {
    for (auto&& vk : v) {
      if (vk!=vk || vk==std::numeric_limits<T>::infinity() ||
          vk==-std::numeric_limits<T>::infinity()) return false;
    }
    return true;
  }

  // Create a temporary file
  CASADI_EXPORT std::string temporary_file(const std::string& prefix, const std::string& suffix);

  CASADI_EXPORT void normalized_setup(std::istream& stream);
  CASADI_EXPORT void normalized_setup(std::ostream& stream);

  inline void normalized_out(std::ostream& stream, double val) {
    if (val==std::numeric_limits<double>::infinity()) {
      stream << "inf";
    } else if (val==-std::numeric_limits<double>::infinity()) {
      stream << "-inf";
    } else if (val!=val) {
      stream << "nan";
    } else {
      stream << val;
    }
  }
  inline int normalized_in(std::istream& stream, double& ret) {
    std::streampos start = stream.tellg();
    stream >> ret;
    // Failed to interpret as double?
    if (stream.fail()) {
      // Clear error flag
      stream.clear();
      // Reset stream position
      // Need to parse e.g "-inf"
      stream.seekg(start);
      // Might be a inf/nan
      std::string non_reg;
      stream >> non_reg;
      // Break on trailing whitespace
      if (stream.fail()) {
        if (stream.eof()) {
          ret = std::numeric_limits<double>::quiet_NaN();
          return -1; // End of stream
        } else {
          ret = std::numeric_limits<double>::quiet_NaN();
          return 1; // Failed to parse to string
        }
      }
      if (non_reg=="inf") {
        ret = std::numeric_limits<double>::infinity();
      } else if (non_reg=="-inf") {
        ret = -std::numeric_limits<double>::infinity();
      } else if (non_reg=="nan") {
        ret = std::numeric_limits<double>::quiet_NaN();
      } else {
        ret = std::numeric_limits<double>::quiet_NaN();
        return 2; // Failed to interpret as number
      }
    }
    return 0;
  }

} // namespace casadi

#ifndef SWIG
// In std namespace
namespace std {

  /// Enables flushing an std::vector to a stream (prints representation)
  template<typename T>
  ostream& operator<<(ostream& stream, const vector<T>& v) {
    stream << casadi::str(v);
    return stream;
  }

  /// Enables flushing an std::set to a stream (prints representation)
  template<typename T>
  ostream& operator<<(ostream& stream, const set<T>& v) {
    stream << casadi::str(v);
    return stream;
  }

  template<typename T1, typename T2>
  ostream& operator<<(ostream& stream, const pair<T1, T2>& p) {
    stream << casadi::str(p);
    return stream;
  }

  template<typename T1, typename T2>
  ostream& operator<<(ostream& stream, const std::map<T1, T2>& p) {
    stream << casadi::str(p);
    return stream;
  }

  template<typename T2>
  ostream& operator<<(ostream& stream, const std::map<std::string, T2>& p) {
    stream << casadi::str(p);
    return stream;
  }

  template<typename T>
  bool mul_overflows(const T& a, const T& b) {
    if (a==0 || b==0) return false;
    return abs(std::numeric_limits<T>::max()/a) < abs(b);
  }

} // namespace std

// Implementations
namespace casadi {

  template<typename T, typename S>
  std::vector<T> vector_static_cast(const std::vector<S>& rhs) {
    std::vector<T> ret;
    ret.reserve(rhs.size());
    for (auto e : rhs) ret.push_back(static_cast<T>(e));
    return ret;
  }

  template<typename T>
  std::vector<T> vector_slice(const std::vector<T> &v, const std::vector<casadi_int> &i) {
    std::vector<T> ret;
    ret.reserve(i.size());
    for (casadi_int k=0;k<i.size();++k) {
       casadi_int j = i[k];
       casadi_assert(j>=0,
         "vector_slice: Indices should be larger than zero."
         "You have " + str(j) + " at location " + str(k) + ".");
       casadi_assert(j<v.size(),
         "vector_slice: Indices should be larger than zero."
         "You have " + str(j) + " at location " + str(k) + ".");
       ret.push_back(v[j]);
    }
    return ret;
  }

  template<typename T>
  std::vector<T> reverse(const std::vector<T> &v) {
    std::vector<T> ret(v.size());
    std::reverse_copy(v.begin(), v.end(), ret.begin());
    return ret;
  }

  template<typename T>
  std::vector<T> join(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<T> ret = a;
    ret.insert(ret.end(), b.begin(), b.end());
    return ret;
  }

  template<typename T>
  std::vector<T> join(const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c) {
    std::vector<T> ret = a;
    ret.insert(ret.end(), b.begin(), b.end());
    ret.insert(ret.end(), c.begin(), c.end());
    return ret;
  }

  template<typename T>
  std::vector<T> permute(const std::vector<T> &a, const std::vector<casadi_int> &order) {
    casadi_assert_dev(order.size()==a.size());
    casadi_assert_dev(is_permutation(order));
    return vector_slice(a, order);
  }

  template<typename T>
  std::vector<casadi_int> find(const std::vector<T> &v) {
    std::vector<casadi_int> ret;
    for (casadi_int i=0;i<v.size();++i) {
      if (v[i]) ret.push_back(i);
    }
    return ret;
  }

#ifndef SWIG
  template<class T>
  std::vector<T> applymap(T (*f)(const T&) , const std::vector<T>& comp) {
    std::vector<T> ret(comp.size());
    std::transform(comp.begin(), comp.end(), ret.begin(), f);
    return ret;
  }

  template<class T>
  void applymap(void (*f)(T &), std::vector<T>& comp) {
    std::for_each(comp.begin(), comp.end(), f);
  }

  template<class S, class D>
  void copy_vector(const std::vector<S>& s, std::vector<D>& d) {
    casadi_assert(s.size()==d.size(), "Dimension mismatch.");
    std::copy(s.begin(), s.end(), d.begin());
  }

  template<class S, class D>
  void assign_vector(const std::vector<S>& s, std::vector<D>& d) {
    casadi_assert(d.empty(), "Receiving vector must be empty");
    d.resize(s.size());
    std::copy(s.begin(), s.end(), d.begin());
  }

  template<class S, class D>
  void copy_vector(const S* s, std::vector<D>& d) {
    for (casadi_int i=0;i<d.size();++i) {
      d[i] = static_cast<D>(s[i]);
    }
  }

  template<class S, class D>
  void init_vector(std::vector<S>& d, const std::vector<D>& s) {
    d.resize(s.size());
    std::copy(s.begin(), s.end(), d.begin());
  }
#endif //SWIG

  template<typename T>
  bool in_range(const std::vector<T> &v, casadi_int upper) {
    return in_range(v, 0, upper);
  }

  template<typename T>
  bool in_range(const std::vector<T> &v, casadi_int lower, casadi_int upper) {
    if (v.empty()) return true;
    casadi_int max = *std::max_element(v.begin(), v.end());
    if (max >= upper) return false;
    casadi_int min = *std::min_element(v.begin(), v.end());
    return (min >= lower);
  }

  template<class T, class S>
  void flatten_nested_vector(const std::vector< std::vector<T> >& nested,
                             std::vector<S>& flat) {
    // Count total elements in nested
    casadi_int N = 0;
    for (const auto& e : nested) {
      N += e.size();
    }

    // Populate flat, one nested section at a time
    flat.clear();
    flat.reserve(N);
    for (const auto& e : nested) {
      flat.insert(flat.end(), e.begin(), e.end());
    }
  }

  template<class T, class S, class I>
  void flatten_nested_vector(const std::vector< std::vector<T> >& nested,
                             std::vector<S>& flat,
                             std::vector<I>& indices) {
    // Delegate
    flatten_nested_vector(nested, flat);

    // Build up indices
    casadi_int N = nested.size();
    indices.resize(1, 0);
    indices.reserve(N+1);
    casadi_int offset = 0;
    for (const auto& e : nested) {
      offset += e.size();
      indices.push_back(offset);
    }
  }

  template<typename T>
  bool isUnique(const std::vector<T> &v) {
    std::set<T> s(v.begin(), v.end());
    return v.size()==s.size();
  }

  template<typename T>
  bool is_increasing(const std::vector<T> &v) {
    if (v.empty()) return true;
    T el = v[0];
    for (casadi_int i=1;i<v.size();++i) {
      if (!(v[i] > el)) return false;
      el = v[i];
    }
    return el==el; // nan -> false
  }

  template<typename T>
  bool is_decreasing(const std::vector<T> &v) {
    if (v.empty()) return true;
    T el = v[0];
    for (casadi_int i=1;i<v.size();++i) {
      if (!(v[i] < el)) return false;
      el = v[i];
    }
    return el==el; // nan -> false
  }

  template<typename T>
  bool is_nonincreasing(const std::vector<T> &v) {
    if (v.empty()) return true;
    T el = v[0];
    for (casadi_int i=1;i<v.size();++i) {
      if (!(v[i] <= el)) return false;
      el = v[i];
    }
    return el==el; // nan -> false
  }

  template<typename T>
  bool is_nondecreasing(const std::vector<T> &v) {
    if (v.empty()) return true;
    T el = v[0];
    for (casadi_int i=1;i<v.size();++i) {
      if (!(v[i] >= el)) return false;
      el = v[i];
    }
    return el==el; // nan -> false
  }

  template<typename T>
  bool is_monotone(const std::vector<T> &v) {
    return is_nondecreasing(v) || is_nonincreasing(v);
  }

  template<typename T>
  bool is_strictly_monotone(const std::vector<T> &v) {
    return is_decreasing(v) || is_increasing(v);
  }

  template<typename T>
  bool has_negative(const std::vector<T> &v) {
    for (std::size_t i=0; i<v.size(); ++i) {
      if (v[i]<0) return true;
    }
    return false;
  }

  template<typename T>
  void write_matlab(std::ostream &stream, const std::vector<T> &v) {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(stream, " "));
  }

  template<typename T>
  void write_matlab(std::ostream &stream, const std::vector<std::vector<T> > &v) {
    for (casadi_uint i=0; i<v.size(); ++i) {
      std::copy(v[i].begin(), v[i].end(), std::ostream_iterator<T>(stream, " "));
      stream << std::endl;
    }
  }

  template<typename T>
  void read_matlab(std::istream &stream, std::vector<T> &v) {
    v.clear();

    while (!stream.eof()) {
      T val;
      stream >> val;
      if (stream.fail()) {
        stream.clear();
        std::string s;
        stream >> s;
        if (s=="inf")
          val = std::numeric_limits<T>::infinity();
        else
          break;
      }
      v.push_back(val);
    }
  }

  template<typename T>
  void read_matlab(std::ifstream &file, std::vector<std::vector<T> > &v) {
    v.clear();
    std::string line;
    while (!getline(file, line, '\n').eof()) {
      std::istringstream reader(line);
      std::vector<T> lineData;

      while (!reader.eof()) {
        T val;
        reader >> val;
        if (reader.fail()) {
          reader.clear();
          std::string s;
          reader >> s;
          if (s=="inf")
            val = std::numeric_limits<T>::infinity();
          else
            break;
        }
        lineData.push_back(val);
      }
      v.push_back(lineData);
    }
  }

  template<typename T, typename F, typename L>
  void linspace(std::vector<T> &v, const F& first, const L& last) {
    if (v.size()<2)
        throw CasadiException("std::linspace: vector must contain at least two elements");

    // Increment
    T increment = (last-first)/T(v.size()-1);

    v[0] = first;
    for (unsigned i=1; i<v.size()-1; ++i)
      v[i] = v[i-1] + increment;
    v[v.size()-1] = last;
  }

  template<typename T>
  T* get_ptr(std::vector<T> &v) {
    if (v.empty())
      return nullptr;
    else
      return &v.front();
  }

  template<typename T>
  const T* get_ptr(const std::vector<T> &v) {
    if (v.empty())
      return nullptr;
    else
      return &v.front();
  }

  // Helper class
  template<typename T>
  struct sortCompare {
    const std::vector<T> &v_;
    sortCompare(const std::vector<T> &v) : v_(v) {}
    bool operator() (casadi_int i, casadi_int j) const { return v_[i]<v_[j];}
  };

  template<typename T>
  void sort(const std::vector<T> &values, std::vector<T> &sorted_values,
            std::vector<casadi_int> &indices, bool invert_indices) {
    // Call recursively if indices need to be inverted
    if (invert_indices) {
      std::vector<casadi_int> inverted;
      sort(values, sorted_values, inverted, false);
      indices.resize(inverted.size());
      for (size_t i=0; i<inverted.size(); ++i) {
        indices[inverted[i]] = i;
      }
      return;
    }

    // Create list of indices
    indices.resize(values.size());
    for (size_t i=0; i<indices.size(); ++i) indices[i] = i;

    // Sort this list by the values
    std::sort(indices.begin(), indices.end(), sortCompare<T>(values));

    // Sort the values accordingly
    sorted_values.resize(values.size());
    for (size_t i=0; i<values.size(); ++i) {
      sorted_values[i] = values[indices[i]];
    }
  }

  template<typename T>
  T product(const std::vector<T> &values) {
    T r = 1;
    for (casadi_int i=0;i<values.size();++i) r*=values[i];
    return r;
  }

  template<typename T>
  T sum(const std::vector<T> &values) {
    T r = 0;
    for (casadi_int i=0;i<values.size();++i) r+=values[i];
    return r;
  }

  template<typename T>
  std::vector<T> cumsum(const std::vector<T> &values) {
    std::vector<T> ret(values.size());
    T acc = 0;
    for (casadi_int i=0;i<values.size();++i) {
      acc+= values[i];
      ret[i] = acc;
    }
    return ret;
  }

  template<typename T>
  std::vector<T> cumsum0(const std::vector<T> &values) {
    std::vector<T> ret(values.size()+1, 0);
    T acc = 0;
    for (casadi_int i=0;i<values.size();++i) {
      acc+= values[i];
      ret[i+1] = acc;
    }
    return ret;
  }

  template<typename T>
  std::vector<T> diff(const std::vector<T> &values) {
    casadi_assert(!values.empty(), "Array must be non-empty");
    std::vector<T> ret(values.size()-1);
    for (casadi_int i=0;i<values.size()-1;++i) {
      ret[i] = values[i+1]-values[i];
    }
    return ret;
  }

  template<typename T>
  T dot(const std::vector<T>& a, const std::vector<T>& b) {
    T ret = 0;
    for (casadi_int k=0; k<a.size(); ++k) {
      ret += a[k]*b[k];
    }
    return ret;
  }

  template<typename T>
  T norm_inf(const std::vector<T>& x) {
    T ret = 0;
    for (casadi_int k=0; k<x.size(); ++k) {
      ret = fmax(ret, fabs(x[k]));
    }
    return ret;
  }

  template<typename T>
  T norm_1(const std::vector<T>& x) {
    T ret = 0;
    for (casadi_int k=0; k<x.size(); ++k) {
      ret += fabs(x[k]);
    }
    return ret;
  }

  template<typename T>
  T norm_2(const std::vector<T>& x) {
    T ret = 0;
    for (casadi_int k=0; k<x.size(); ++k) {
      ret += x[k]*x[k];
    }
    return sqrt(ret);
  }

  template<typename T>
  bvec_t* get_bvec_t(std::vector<T>& v) {
    casadi_assert(0, "get_bvec_t only supported for double");
  }

  template<typename T>
  const bvec_t* get_bvec_t(const std::vector<T>& v) {
    casadi_assert(0, "get_bvec_t only supported for double");
  }

  ///@{
  /// Readability typedefs
  typedef std::vector<std::string> StringVector;
  ///@}

} // namespace casadi
#endif // SWIG

#endif // CASADI_MISC_HPP
