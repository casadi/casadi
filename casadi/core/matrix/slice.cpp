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


#include "slice.hpp"
#include "matrix_tools.hpp"

using namespace std;
namespace casadi {

  Slice::Slice() : start_(0), stop_(std::numeric_limits<int>::max()), step_(1) {
  }

  Slice::Slice(int i) : start_(i), stop_(i+1), step_(1) {
    if (i==-1) stop_ = std::numeric_limits<int>::max();
  }

  Slice::Slice(int start, int stop, int step) : start_(start), stop_(stop), step_(step) {
  }

  std::vector<int> Slice::getAll(int len, bool ind1) const {
    int start;
    int stop;
    if (start_==std::numeric_limits<int>::min()) {
      start = (step_ < 0) ? len - 1 : 0;
    } else {
      start = start_;
      if (start<0) start+=len;
    }
    if (stop_==std::numeric_limits<int>::max()) {
      stop = (step_ < 0) ? -1 : len;
    } else {
      stop = stop_;
      if (stop<0) stop+=len;
    }

    casadi_assert_message(stop<=len,
                          "Slice (start=" << start << ", stop=" << stop << ", step=" << step_
                          << ") out of bounds with supplied length of " << len);
    casadi_assert_message(start>=0,
                          "Slice (start=" << start << ", stop=" << stop << ", step=" << step_
                          << ") out of bounds with start<0.");
    if ((stop>=start && step_<0) || (stop<=start && step_>0)) return std::vector<int>();

    return range(start+ind1, stop+ind1, step_, len+ind1);
  }

  void Slice::repr(std::ostream& stream, bool trailing_newline) const {
    print(stream, trailing_newline);
  }

  void Slice::print(std::ostream& stream, bool trailing_newline) const {
    bool from_beginning = start_ == 0;
    bool till_end = stop_ == std::numeric_limits<int>::max();
    bool skip_none = step_==1;
    if (stop_==start_+1) {
      stream << start_;
    } else {
      if (!from_beginning) stream << start_;
      stream << ":";
      if (!till_end) stream << stop_;
      if (!skip_none) stream << ":" << step_;
    }
    if (trailing_newline) stream << std::endl;
  }

  Slice::Slice(const std::vector<int>& v, bool ind1) {
    casadi_assert_message(isSlice(v, ind1), "Cannot be represented as a Slice");
    if (v.size()==0) {
      start_=stop_=0;
      step_ = 1;
    } else if (v.size()==1) {
      start_ = v.front()-ind1;
      stop_ = start_ + 1;
      step_ = 1;
    } else {
      start_ = v[0]-ind1;
      step_ = v[1]-v[0];
      stop_ = start_ + step_*v.size();
    }
  }

  bool Slice::isSlice(const std::vector<int>& v, bool ind1) {
    // Always false if negative numbers or non-increasing
    int last_v = -1;
    for (int i=0; i<v.size(); ++i) {
      if (v[i]-ind1<=last_v) return false;
      last_v = v[i]-ind1;
    }

    // Always true if less than 2 elements
    if (v.size()<2) return true;

    // If two elements, true if they are different
    if (v.size()==2) return v[0]!=v[1];

    // We can now get the beginning, end and step
    int start = v[0]-ind1;
    int step = v[1]-v[0];
    //int stop = start + step*v.size();

    // Consistency check
    for (int i=2; i<v.size(); ++i) {
      if (v[i]-ind1!=start+i*step) return false;
    }

    // True if reached this point
    return true;
  }

  bool Slice::isSlice2(const std::vector<int>& v) {
    // Always true if 1D slice
    if (isSlice(v)) return true;

    // Always false if negative numbers or non-increasing
    int last_v = -1;
    for (int i=0; i<v.size(); ++i) {
      if (v[i]<=last_v) return false;
      last_v = v[i];
    }

    // Get the slices
    int start_outer = 0;
    int step_outer = -1;
    int start_inner = v.front();
    int step_inner = v[1]-v[0];
    int stop_inner = -1;
    for (int i=2; i<v.size(); ++i) {
      int predicted_v = start_inner+i*step_inner;
      if (v[i]!=predicted_v) {
        stop_inner = predicted_v;
        step_outer = v[i] - start_inner;
        break;
      }
    }
    casadi_assert(stop_inner>=0);

    // Get the end of the outer slice
    int stop_outer = v.back();
    do {
      if (step_outer>0) stop_outer++;
      else             stop_outer--;
    } while (stop_outer % step_outer!=0);

    // Check consistency
    std::vector<int>::const_iterator it=v.begin();
    for (int i=start_outer; i!=stop_outer; i+=step_outer) {
      for (int j=i+start_inner; j!=i+stop_inner; j+=step_inner) {
        // False if we've reached the end
        if (it==v.end()) return false;

        // Check if value matches
        if (*it++ != j) return false;
      }
    }

    // False if there are still elements not accounted for
    if (it!=v.end()) return false;

    // True if reached this point
    return true;
  }

  Slice::Slice(const std::vector<int>& v, Slice& outer) {
    casadi_assert_message(isSlice2(v), "Cannot be represented as a nested Slice");

    // If simple slice
    if (isSlice(v)) {
      *this = Slice(v);
      outer.start_ = 0;
      outer.step_ = outer.stop_ = stop_;
      return;
    }

    // Get the slices
    outer.start_ = 0;
    outer.step_ = -1;
    start_ = v.front();
    step_ = v[1]-v[0];
    stop_ = -1;
    for (int i=2; i<v.size(); ++i) {
      int predicted_v = start_+i*step_;
      if (v[i]!=predicted_v) {
        stop_ = predicted_v;
        outer.step_ = v[i] - start_;
        break;
      }
    }

    // Get the end of the outer slice
    outer.stop_ = v.back();
    do {
      if (outer.step_>0) outer.stop_++;
      else              outer.stop_--;
    } while (outer.stop_ % outer.step_!=0);
  }

  std::vector<int> Slice::getAll(const Slice& outer, int len) const {
    std::vector<int> ret;
    for (int i=outer.start_; i!=outer.stop_; i+=outer.step_) {
      for (int j=i+start_; j!=i+stop_; j+=step_) {
        ret.push_back(j);
      }
    }
    return ret;
  }

  bool Slice::isScalar(int len) const {
    int start = std::min(start_, len);
    int stop = std::min(stop_, len);
    int nret = (stop-start)/step_ + ((stop-start)%step_!=0);
    return nret==1;
  }

  int Slice::toScalar(int len) const {
    casadi_assert(isScalar(len));
    casadi_assert_message(start_ >= -len && start_ < len, "Slice::getScalar: out of bounds");
    return start_ >= 0 ? start_ : start_+len;
  }

} // namespace casadi

