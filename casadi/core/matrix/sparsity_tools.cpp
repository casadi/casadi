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


#include "sparsity_tools.hpp"
#include "../casadi_exception.hpp"
#include "../std_vector_tools.hpp"
#include "../matrix/matrix.hpp"
#include "../matrix/matrix_tools.hpp"

using namespace std;

namespace casadi {

  Sparsity reshape(const Sparsity& a, int nrow, int ncol) {
    return a.reshape(nrow, ncol);
  }

  Sparsity vec(const Sparsity& a) {
    return reshape(a, a.numel(), 1);
  }

  Sparsity mul(const Sparsity& a, const Sparsity &b) {
    return a.patternProductNew(b);
  }

  Sparsity mul(const std::vector<Sparsity>& s) {
    if (s.size()==0) return Sparsity();
    Sparsity ret = s[0];
    for (int i=1;i<s.size();++i) {
      ret = mul(ret, s[i]);
    }
    return ret;
  }

  int rank(const Sparsity& a) {
    std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    a.dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
    return coarse_colblock.at(3);
  }

  Sparsity horzcat(const std::vector<Sparsity> & sp) {
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

  Sparsity horzcat(const Sparsity & a, const Sparsity & b) {
    Sparsity ret = a;
    ret.appendColumns(b);
    return ret;
  }

  Sparsity vertcat(const std::vector<Sparsity> & sp) {
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

  Sparsity blockcat(const std::vector< std::vector< Sparsity > > &v) {
    std::vector< Sparsity > ret;
    for (int i=0; i<v.size(); ++i)
      ret.push_back(horzcat(v[i]));
    return vertcat(ret);
  }

  Sparsity vertcat(const Sparsity & a, const Sparsity & b) {
    if (a.isVector()) {
      Sparsity ret = a;
      ret.append(b);
      return ret;
    } else {
      Sparsity ret = a.T();
      ret.appendColumns(b.T());
      return ret.T();
    }
  }

  Sparsity blkdiag(const std::vector< Sparsity > &v) {
    int n = 0;
    int m = 0;

    std::vector<int> colind(1, 0);
    std::vector<int> row;

    int nz = 0;
    for (int i=0;i<v.size();++i) {
      const std::vector<int> &colind_ = v[i].colind();
      const std::vector<int> &row_ = v[i].row();
      for (int k=1;k<colind_.size();++k) {
        colind.push_back(colind_[k]+nz);
      }
      for (int k=0;k<row_.size();++k) {
        row.push_back(row_[k]+m);
      }
      n+= v[i].size2();
      m+= v[i].size1();
      nz+= v[i].size();
    }

    return Sparsity(m, n, colind, row);
  }

  Sparsity blkdiag(const Sparsity &a, const Sparsity &b) {

    std::vector<Sparsity> v;
    v.push_back(a);
    v.push_back(b);

    return blkdiag(v);
  }

  std::vector<Sparsity> diagsplit(const Sparsity& sp,
     const std::vector<int>& offset1,
     const std::vector<int>& offset2) {
    // Consistency check
    casadi_assert(offset1.size()>=1);
    casadi_assert(offset1.front()==0);
    casadi_assert_message(offset1.back()==sp.size1(),
                          "diagsplit(Sparsity, offset1, offset2): Last elements of offset1 "
                          "(" << offset1.back() << ") must equal the number of rows "
                          "(" << sp.size1() << ")");
    casadi_assert_message(offset2.back()==sp.size2(),
                          "diagsplit(Sparsity, offset1, offset2): Last elements of offset2 "
                          "(" << offset2.back() << ") must equal the number of rows "
                          "(" << sp.size2() << ")");
    casadi_assert(isMonotone(offset1));
    casadi_assert(isMonotone(offset2));
    casadi_assert(offset1.size()==offset2.size());

    // Number of outputs
    int n = offset1.size()-1;

    // Return value
    std::vector<Sparsity> ret;

    // Caveat: this is a very silly implementation
    IMatrix x = IMatrix::zeros(sp);

    for (int i=0;i<n;++i) {
      ret.push_back(x(
          Slice(offset1[i], offset1[i+1]),
          Slice(offset2[i], offset2[i+1])).sparsity());
    }

    return ret;
  }

  std::vector<Sparsity> horzsplit(const Sparsity& sp, const std::vector<int>& offset) {
    // Consistency check
    casadi_assert(offset.size()>=1);
    casadi_assert(offset.front()==0);
    casadi_assert_message(offset.back()==sp.size2(),
                          "horzsplit(Sparsity, std::vector<int>): Last elements of offset "
                          "(" << offset.back() << ") must equal the number of columns "
                          "(" << sp.size2() << ")");
    casadi_assert(isMonotone(offset));

    // Number of outputs
    int n = offset.size()-1;

    // Get the sparsity of the input
    const vector<int>& colind_x = sp.colind();
    const vector<int>& row_x = sp.row();

    // Allocate result
    std::vector<Sparsity> ret;
    ret.reserve(n);

    // Sparsity pattern as CCS vectors
    vector<int> colind, row;
    int ncol, nrow = sp.size1();

    // Get the sparsity patterns of the outputs
    for (int i=0; i<n; ++i) {
      int first_col = offset[i];
      int last_col = offset[i+1];
      ncol = last_col - first_col;

      // Construct the sparsity pattern
      colind.resize(ncol+1);
      copy(colind_x.begin()+first_col, colind_x.begin()+last_col+1, colind.begin());
      for (vector<int>::iterator it=colind.begin()+1; it!=colind.end(); ++it) *it -= colind[0];
      colind[0] = 0;
      row.resize(colind.back());
      copy(row_x.begin()+colind_x[first_col], row_x.begin()+colind_x[last_col], row.begin());

      // Append to the list
      ret.push_back(Sparsity(nrow, ncol, colind, row));
    }

    // Return (RVO)
    return ret;
  }

  std::vector<Sparsity> vertsplit(const Sparsity& sp, const std::vector<int>& offset) {
    std::vector<Sparsity> ret = horzsplit(sp.T(), offset);
    for (std::vector<Sparsity>::iterator it=ret.begin(); it!=ret.end(); ++it) {
      *it = it->T();
    }
    return ret;
  }

  std::vector<Sparsity> diagsplit(const Sparsity& x, int incr) {
    casadi_assert(incr>=1);
    casadi_assert_message(x.isSquare(),
      "diagsplit(x,incr)::input must be square but got " << x.dimString()  << ".");
    vector<int> offset2 = range(0, x.size2(), incr);
    offset2.push_back(x.size2());
    return diagsplit(x, offset2);
  }

  std::vector<Sparsity> diagsplit(const Sparsity& x, int incr1, int incr2) {
    casadi_assert(incr1>=1);
    casadi_assert(incr2>=1);
    vector<int> offset1 = range(0, x.size1(), incr1);
    offset1.push_back(x.size1());
    vector<int> offset2 = range(0, x.size2(), incr2);
    offset2.push_back(x.size2());
    return diagsplit(x, offset1, offset2);
  }

  std::vector<Sparsity> diagsplit(const Sparsity& x, const std::vector<int>& output_offset) {
    casadi_assert_message(x.isSquare(),
      "diagsplit(x,incr)::input must be square but got " << x.dimString()  << ".");
    return diagsplit(x, output_offset, output_offset);
  }

} // namespace casadi
