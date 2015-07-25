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


#ifndef CASADI_SPARSITY_INTERFACE_HPP
#define CASADI_SPARSITY_INTERFACE_HPP

#include "../std_vector_tools.hpp"

namespace casadi {
  /** \brief Sparsity interface class

      This is a common base class for GenericMatrix (i.e. MX and Matrix<>) and Sparsity, introducing a
      uniform syntax and implementing common functionality using the curiously recurring template pattern
      (CRTP) idiom.\n

      \author Joel Andersson
      \date 2014
  */
  template<typename MatType>
  class CASADI_EXPORT SparsityInterface {
#ifndef SWIG
  protected:
    // Helper functions
    inline const MatType& self() const { return static_cast<const MatType&>(*this); }
    inline MatType& self() { return static_cast<MatType&>(*this); }
  public:
    // Create vector with 1 element
    inline friend std::vector<MatType> make_vector(const MatType& x0) {
      return std::vector<MatType>(1, x0);
    }

    // Create vector with 2 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1) {
      MatType x[] = {x0, x1};
      return std::vector<MatType>(x, x+2);
    }

    // Create vector with 3 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2) {
      MatType x[] = {x0, x1, x2};
      return std::vector<MatType>(x, x+3);
    }

    // Create vector with 4 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2, const MatType& x3) {
      MatType x[] = {x0, x1, x2, x3};
      return std::vector<MatType>(x, x+4);
    }

    // Create vector with 5 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2, const MatType& x3,
                                                   const MatType& x4) {
      MatType x[] = {x0, x1, x2, x3, x4};
      return std::vector<MatType>(x, x+5);
    }

    // Create vector with 6 elements
    inline friend std::vector<MatType> make_vector(const MatType& x0, const MatType& x1,
                                                   const MatType& x2, const MatType& x3,
                                                   const MatType& x4, const MatType& x5) {
      MatType x[] = {x0, x1, x2, x3, x4, x5};
      return std::vector<MatType>(x, x+6);
    }

    // Create vector from map and vector with key order
    inline friend std::vector<MatType> make_vector(const std::map<std::string, MatType>& m,
                                                   const std::vector<std::string>& s) {
      std::vector<MatType> ret(s.size());
      for (size_t i=0; i!=s.size(); ++i) {
        typename std::map<std::string, MatType>::const_iterator it=m.find(s[i]);
        if (it!=m.end()) {
          ret[i]=it->second;
        }
      }
      return ret;
    }

    // Create vector from map and vector with key order
    inline friend std::vector<std::vector<MatType> >
      make_vector(const std::map<std::string, std::vector<MatType> >& m,
                  const std::vector<std::string>& s) {
      std::vector<std::vector<MatType> > ret(s.size());
      for (size_t i=0; i!=s.size(); ++i) {
        typename std::map<std::string, std::vector<MatType> >::const_iterator it=m.find(s[i]);
        if (it!=m.end()) {
          ret[i]=it->second;
        }
      }
      return ret;
    }

    // Create vector from map and vector with key order
    inline friend std::vector<MatType>
      make_vector(const std::pair<std::map<std::string, MatType>, std::vector<std::string> >& ms) {
      return make_vector(ms.first, ms.second);
    }

    // Create vector from map and vector with key order
    inline friend std::vector<std::vector<MatType> >
      make_vector(const std::pair<std::map<std::string, std::vector<MatType> >,
                  std::vector<std::string> >& ms) {
      return make_vector(ms.first, ms.second);
    }

    // Assign 1 element from a vector
    template<typename T0>
      inline friend void assign_vector(T0& x0,
                                       const std::vector<MatType>& x) {
      x.at(0).get(x0);
    }

    // Assign 2 elements from a vector
    template<typename T0, typename T1>
      inline friend void assign_vector(T0& x0, T1& x1,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x);
      x.at(1).get(x1);
    }

    // Assign 3 elements from a vector
    template<typename T0, typename T1, typename T2>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x);
      x.at(2).get(x2);
    }

    // Assign 4 elements from a vector
    template<typename T0, typename T1, typename T2, typename T3>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2, T3& x3,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x2, x);
      x.at(3).get(x3);
    }

    // Assign 5 elements from a vector
    template<typename T0, typename T1, typename T2, typename T3, typename T4>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2, T3& x3, T4& x4,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x2, x3, x);
      x.at(4).get(x4);
    }

    // Assign 6 elements from a vector
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
      inline friend void assign_vector(T0& x0, T1& x1, T2& x2, T3& x3, T4& x4, T5& x5,
                                       const std::vector<MatType>& x) {
      assign_vector(x0, x1, x2, x3, x4, x);
      x.at(5).get(x5);
    }

    // Create map with 1 element
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0) {
      std::map<std::string, MatType> ret;
      ret[n0]=x0;
      return ret;
    }

    // Create map with 2 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1) {
      std::map<std::string, MatType> ret=make_map(n0, x0);
      ret[n1]=x1;
      return ret;
    }

    // Create map with 3 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1);
      ret[n2]=x2;
      return ret;
    }

    // Create map with 4 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2,
               const std::string& n3, const MatType& x3) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1, n2, x2);
      ret[n3]=x3;
      return ret;
    }

    // Create map with 5 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2,
               const std::string& n3, const MatType& x3,
               const std::string& n4, const MatType& x4) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1, n2, x2, n3, x3);
      ret[n4]=x4;
      return ret;
    }

    // Create map with 6 elements
    inline friend std::map<std::string, MatType>
      make_map(const std::string& n0, const MatType& x0,
               const std::string& n1, const MatType& x1,
               const std::string& n2, const MatType& x2,
               const std::string& n3, const MatType& x3,
               const std::string& n4, const MatType& x4,
               const std::string& n5, const MatType& x5) {
      std::map<std::string, MatType> ret=make_map(n0, x0, n1, x1, n2, x2, n3, x3, n4, x4);
      ret[n5]=x5;
      return ret;
    }

    // Assign 1 element from a map
    template<typename T0>
    inline friend void assign_map(const std::string& n0, T0& x0,
                                  const std::map<std::string, MatType>& x) {
      x.at(n0).get(x0);
    }

    // Assign 2 elements from a map
    template<typename T0, typename T1>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, x);
      x.at(n1).get(x1);
    }

    // Assign 3 elements from a map
    template<typename T0, typename T1, typename T2>
    inline friend void assign_map(const std::string& n0, T0& x0,
                                  const std::string& n1, T0& x1,
                                  const std::string& n2, T0& x2,
                                  const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, x);
      x.at(n2).get(x2);
    }

    // Assign 4 elements from a map
    template<typename T0, typename T1, typename T2, typename T3>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::string& n2, T0& x2,
                                    const std::string& n3, T0& x3,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, n2, x2, x);
      x.at(n3).get(x3);
    }

    // Assign 5 elements from a map
    template<typename T0, typename T1, typename T2, typename T3, typename T4>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::string& n2, T0& x2,
                                    const std::string& n3, T0& x3,
                                    const std::string& n4, T0& x4,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, n2, x2, n3, x3, x);
      x.at(n4).get(x4);
    }

    // Assign 6 elements from a map
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
      inline friend void assign_map(const std::string& n0, T0& x0,
                                    const std::string& n1, T0& x1,
                                    const std::string& n2, T0& x2,
                                    const std::string& n3, T0& x3,
                                    const std::string& n4, T0& x4,
                                    const std::string& n5, T0& x5,
                                    const std::vector<MatType>& x) {
      assign_map(n0, x0, n1, x1, n2, x2, n3, x3, n4, x4, x);
      x.at(n5).get(x5);
    }
#endif // SWIG

  public:

    /// \cond CLUTTER
    std::vector< std::vector< MatType > >
      zz_blocksplit(const std::vector<int>& vert_offset, const std::vector<int>& horz_offset) const;
    static MatType zz_veccat(const std::vector< MatType >& x);
    MatType zz_vec() const;
    MatType zz_repmat(int n, int m=1) const;
    static std::vector<int> zz_offset(const std::vector< MatType > &v, bool vert=true);
    /// \endcond

#ifndef SWIG
#include "sparsity_interface_friends.hpp"
#endif // SWIG
  };

#ifndef SWIG
  template<typename MatType>
  MatType SparsityInterface<MatType>::zz_vec() const {
    return reshape(self(), self().numel(), 1);
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::zz_repmat(int n, int m) const {
    MatType allrows = vertcat(std::vector<MatType>(n, self()));
    return horzcat(std::vector<MatType>(m, allrows));
  }

  template<typename MatType>
  std::vector< std::vector< MatType > >
  SparsityInterface<MatType>::zz_blocksplit(const std::vector<int>& vert_offset,
                                            const std::vector<int>& horz_offset) const {
    std::vector< MatType > rows = vertsplit(self(), vert_offset);
    std::vector< std::vector< MatType > > ret;
    for (int i=0;i<rows.size();++i) {
      ret.push_back(horzsplit(rows[i], horz_offset));
    }
    return ret;
  }

  template<typename MatType>
  std::vector<int>
  SparsityInterface<MatType>::zz_offset(const std::vector< MatType > &v, bool vert) {
    std::vector<int> ret(v.size()+1);
    ret[0]=0;
    for (int i=0; i<v.size(); ++i) {
      ret[i+1] = ret[i] + (vert ? v[i].size1() : v[i].size2());
    }
    return ret;
  }

  template<typename MatType>
  MatType SparsityInterface<MatType>::zz_veccat(const std::vector< MatType >& x) {
    std::vector< MatType > x_vec = x;
    for (typename std::vector< MatType >::iterator it=x_vec.begin();
         it!=x_vec.end(); ++it) {
      *it = vec(*it);
    }
    return vertcat(x_vec);
  }
#endif // SWIG

} // namespace casadi

#endif // CASADI_SPARSITY_INTERFACE_HPP
