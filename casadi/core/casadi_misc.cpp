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
#ifdef HAVE_MKSTEMPS
#include <unistd.h>
#endif
using namespace std;

namespace casadi {
  bool all(const std::vector<bool>& v) {
    for (auto && e : v) {
      if (!e) return false;
    }
    return true;
  }

  bool any(const std::vector<bool>& v) {
    for (auto && e : v) {
      if (e) return true;
    }
    return false;
  }

  std::vector<casadi_int> range(casadi_int start, casadi_int stop,
      casadi_int step, casadi_int len) {
    start = std::min(start, len);
    stop = std::min(stop, len);
    casadi_int nret = (stop-start)/step + ((stop-start)%step!=0);
    std::vector<casadi_int> ret(nret);
    casadi_int ind = start;
    for (std::vector<casadi_int>::iterator it=ret.begin(); it!=ret.end(); ++it) {
      *it = ind;
      ind += step;
    }
    return ret;
  }

  bool is_equally_spaced(const std::vector<double> &v) {
    if (v.size()<=1) return true;

    double margin = (v[v.size()-1]-v[0])*1e-14;

    for (casadi_int i=2;i<v.size();++i) {
      double ref = v[0]+(static_cast<double>(i)*(v[v.size()-1]-v[0]))/
        static_cast<double>(v.size()-1);
      if (abs(ref-v[i])>margin) return false;
    }
    return true;
  }

  std::vector<casadi_int> range(casadi_int stop) {
    return range(0, stop);
  }

  std::vector<casadi_int> complement(const std::vector<casadi_int> &v, casadi_int size) {
    casadi_assert(in_range(v, size),
                          "complement: out of bounds. Some elements in v fall out of [0, size[");
    std::vector<casadi_int> lookup(size, 0);
    std::vector<casadi_int> ret;

    for (casadi_int i=0;i<v.size();i++) {
      lookup[v[i]] = 1;
    }

    for (casadi_int i=0;i<size;i++) {
      if (lookup[i]==0) ret.push_back(i);
    }

    return ret;

  }

  std::vector<casadi_int> lookupvector(const std::vector<casadi_int> &v, casadi_int size) {
    casadi_assert(in_range(v, size),
                          "lookupvector: out of bounds. Some elements in v fall out of [0, size[");
    std::vector<casadi_int> lookup(size, -1);

    for (casadi_int i=0;i<v.size();i++) {
      lookup[v[i]] = i;
    }
    return lookup;
  }

  std::vector<casadi_int> lookupvector(const std::vector<casadi_int> &v) {
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
    for (casadi_int i=0;i<l.size();++i) {
      if (i>0) ss << delim;
      ss << l[i];
    }
    return ss.str();
  }

  std::string temporary_file(const std::string& prefix, const std::string& suffix) {
    #ifdef HAVE_MKSTEMPS
    // Preferred solution
    string ret = prefix + "XXXXXX" + suffix;
    if (mkstemps(&ret[0], static_cast<int>(suffix.size())) == -1) {
      casadi_error("Failed to create temporary file: '" + ret + "'");
    }
    return ret;
    #else // HAVE_MKSTEMPS
    // Fallback, may result in deprecation warnings
    return prefix + string(tmpnam(nullptr)) + suffix;
    #endif // HAVE_MKSTEMPS
  }


} // namespace casadi
