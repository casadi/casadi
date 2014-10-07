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


#include "polynomial.hpp"
#include "std_vector_tools.hpp"
#include <sstream>

using namespace std;
namespace casadi {

  Polynomial::Polynomial(real_t scalar) : p_(1, scalar) {
  }

  Polynomial::Polynomial(real_t p0, real_t p1) {
    p_.resize(2);
    p_[0] = p0;
    p_[1] = p1;
  }

  Polynomial::Polynomial(real_t p0, real_t p1, real_t p2) {
    p_.resize(3);
    p_[0] = p0;
    p_[1] = p1;
    p_[2] = p2;
  }

  Polynomial::Polynomial(real_t p0, real_t p1, real_t p2, real_t p3) {
    p_.resize(2);
    p_[0] = p0;
    p_[1] = p1;
    p_[2] = p2;
    p_[3] = p3;
  }

  void Polynomial::repr(std::ostream &stream, bool trailing_newline) const {
    cout << "poly(" << p_ << ")" << endl;
    if (trailing_newline) stream << std::endl;
  }

  void Polynomial::print(std::ostream &stream, bool trailing_newline) const {
    for (int d=0; d<p_.size(); ++d) {
      if (d==0) {
        stream << p_[d];
      } else if (d==1) {
        stream << " + " << p_[d] << "*x";
      } else {
        stream << " + " << p_[d] << "*x^" << d;
      }
    }
    stream << endl;
    if (trailing_newline) stream << std::endl;
  }

  int Polynomial::degree() const {
    return p_.size()-1;
  }

  Polynomial::real_t Polynomial::toScalar() const {
    casadi_assert(degree()==0);
    return p_.front();
  }

  Polynomial Polynomial::operator*(const Polynomial& a) const {
    vector<real_t> p_ret(p_.size() + a.p_.size() - 1, 0);
    for (int d=0; d<p_.size(); ++d) {
      for (int d_a=0; d_a<a.p_.size(); ++d_a) {
        p_ret[d+d_a] += p_[d] * a.p_[d_a];
      }
    }
    return Polynomial(p_ret);
  }

  Polynomial& Polynomial::operator*=(const Polynomial& d) {
    return *this = *this*d;
  }

  Polynomial Polynomial::operator/(real_t d) const {
    Polynomial ret = *this;
    ret/=d;
    return ret;
  }

  Polynomial& Polynomial::operator/=(real_t d) {
    for (vector<real_t>::iterator it=p_.begin(); it!=p_.end(); ++it) {
      *it /= d;
    }
    return *this;
  }

  Polynomial Polynomial::operator+(const Polynomial& b) const {
    Polynomial ret = *this;
    return ret+=b;
  }

  Polynomial& Polynomial::operator+=(const Polynomial& b) {
    p_.resize(max(p_.size(), b.p_.size()), 0);
    transform(b.p_.begin(), b.p_.end(), p_.begin(), p_.begin(), std::plus<real_t>());
    trim();
    return *this;
  }

  Polynomial Polynomial::operator-(const Polynomial& b) const {
    Polynomial ret = *this;
    return ret-=b;
  }

  Polynomial& Polynomial::operator-=(const Polynomial& b) {
    p_.resize(max(p_.size(), b.p_.size()), 0);
    transform(p_.begin(), p_.begin()+b.p_.size(), b.p_.begin(), p_.begin(), std::minus<real_t>());
    trim();
    return *this;
  }

  void Polynomial::trim() {
    // Remove trailing zeros
    size_t new_size = p_.size();
    vector<real_t>::const_reverse_iterator it=p_.rbegin();
    while (it!=p_.rend() && 0==*it++) new_size--;
    p_.resize(new_size);
  }

  Polynomial Polynomial::derivative() const {
    vector<real_t> ret_p(p_.size()-1);
    for (int k=1; k<p_.size(); ++k) {
      ret_p[k-1] = k*p_[k];
    }
    return Polynomial(ret_p);
  }

  Polynomial Polynomial::anti_derivative() const {
    vector<real_t> ret_p(p_.size()+1);
    ret_p[0] = 0;
    for (int k=0; k<p_.size(); ++k) {
      ret_p[k+1] = p_[k]/(k+1);
    }
    return Polynomial(ret_p);
  }


} // namespace casadi



