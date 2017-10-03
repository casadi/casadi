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


#include "mx.hpp"

#include "casadi_misc.hpp"

namespace casadi {
  std::vector<int> range(int start, int stop, int step, int len) {
    start = std::min(start, len);
    stop = std::min(stop, len);
    int nret = (stop-start)/step + ((stop-start)%step!=0);
    std::vector<int> ret(nret);
    int ind = start;
    for (std::vector<int>::iterator it=ret.begin(); it!=ret.end(); ++it) {
      *it = ind;
      ind += step;
    }
    return ret;
  }

  bool is_equally_spaced(const std::vector<double> &v) {
    if (v.size()<=1) return true;

    double margin = (v[v.size()-1]-v[0])*1e-14;

    for (int i=2;i<v.size();++i) {
      double ref = v[0]+(i*(v[v.size()-1]-v[0]))/(v.size()-1);
      if (abs(ref-v[i])>margin) return false;
    }
    return true;
  }

  std::vector<int> range(int stop) {
    return range(0, stop);
  }

  std::vector<int> complement(const std::vector<int> &v, int size) {
    casadi_assert(in_range(v, size),
                          "complement: out of bounds. Some elements in v fall out of [0, size[");
    std::vector<int> lookup(size, 0);
    std::vector<int> ret;

    for (int i=0;i<v.size();i++) {
      lookup[v[i]] = 1;
    }

    for (int i=0;i<size;i++) {
      if (lookup[i]==0) ret.push_back(i);
    }

    return ret;

  }

  std::vector<int> lookupvector(const std::vector<int> &v, int size) {
    casadi_assert(in_range(v, size),
                          "lookupvector: out of bounds. Some elements in v fall out of [0, size[");
    std::vector<int> lookup(size, -1);

    for (int i=0;i<v.size();i++) {
      lookup[v[i]] = i;
    }
    return lookup;
  }

  std::vector<int> lookupvector(const std::vector<int> &v) {
    casadi_assert_dev(!has_negative(v));
    return lookupvector(v, (*std::max_element(v.begin(), v.end()))+1);
  }

  bvec_t* get_bvec_t(std::vector<double>& v) {
    if (v.empty()) {
      return 0;
    } else {
      return reinterpret_cast<bvec_t*>(&v.front());
    }
  }

  /// Get an pointer of sets of booleans from a double vector
  const bvec_t* get_bvec_t(const std::vector<double>& v) {
    if (v.empty()) {
      return 0;
    } else {
      return reinterpret_cast<const bvec_t*>(&v.front());
    }
  }

  std::string join(const std::vector<std::string>& l, const std::string& delim) {
    std::stringstream ss;
    for (int i=0;i<l.size();++i) {
      if (i>0) ss << delim;
      ss << l[i];
    }
    return ss.str();
  }

} // namespace casadi
