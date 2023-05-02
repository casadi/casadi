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


#include "polynomial.hpp"
#include "casadi_misc.hpp"
#include <sstream>

namespace casadi {

  Polynomial::Polynomial(double scalar) : p_(1, scalar) {
  }

  Polynomial::Polynomial(double p0, double p1) {
    p_.resize(2);
    p_[0] = p0;
    p_[1] = p1;
  }

  Polynomial::Polynomial(double p0, double p1, double p2) {
    p_.resize(3);
    p_[0] = p0;
    p_[1] = p1;
    p_[2] = p2;
  }

  Polynomial::Polynomial(double p0, double p1, double p2, double p3) {
    p_.resize(2);
    p_[0] = p0;
    p_[1] = p1;
    p_[2] = p2;
    p_[3] = p3;
  }

  void Polynomial::disp(std::ostream& stream, bool more) const {
    if (more) {
      for (casadi_int d=0; d<p_.size(); ++d) {
        if (d==0) {
          stream << p_[d];
        } else if (d==1) {
          stream << " + " << p_[d] << "*x";
        } else {
          stream << " + " << p_[d] << "*x^" << d;
        }
      }
    } else {
      stream << p_;
    }
  }

  casadi_int Polynomial::degree() const {
    return p_.size()-1;
  }

  double Polynomial::scalar() const {
    casadi_assert_dev(degree()==0);
    return p_.front();
  }

  Polynomial Polynomial::operator*(const Polynomial& a) const {
    std::vector<double> p_ret(p_.size() + a.p_.size() - 1, 0);
    for (casadi_int d=0; d<p_.size(); ++d) {
      for (casadi_int d_a=0; d_a<a.p_.size(); ++d_a) {
        p_ret[d+d_a] += p_[d] * a.p_[d_a];
      }
    }
    return Polynomial(p_ret);
  }

  Polynomial& Polynomial::operator*=(const Polynomial& d) {
    return *this = *this*d;
  }

  Polynomial Polynomial::operator/(double d) const {
    Polynomial ret = *this;
    ret/=d;
    return ret;
  }

  Polynomial& Polynomial::operator/=(double d) {
    for (std::vector<double>::iterator it=p_.begin(); it!=p_.end(); ++it) {
      *it /= d;
    }
    return *this;
  }

  Polynomial Polynomial::operator+(const Polynomial& b) const {
    Polynomial ret = *this;
    return ret+=b;
  }

  Polynomial& Polynomial::operator+=(const Polynomial& b) {
    p_.resize(std::max(p_.size(), b.p_.size()), 0);
    transform(b.p_.begin(), b.p_.end(), p_.begin(), p_.begin(), std::plus<double>());
    trim();
    return *this;
  }

  Polynomial Polynomial::operator-(const Polynomial& b) const {
    Polynomial ret = *this;
    return ret-=b;
  }

  Polynomial& Polynomial::operator-=(const Polynomial& b) {
    p_.resize(std::max(p_.size(), b.p_.size()), 0);
    transform(p_.begin(), p_.begin()+b.p_.size(),
              b.p_.begin(), p_.begin(), std::minus<double>());
    trim();
    return *this;
  }

  void Polynomial::trim() {
    // Remove trailing zeros
    size_t new_size = p_.size();
    std::vector<double>::const_reverse_iterator it=p_.rbegin();
    while (it!=p_.rend() && 0==*it++) new_size--;
    p_.resize(new_size);
  }

  Polynomial Polynomial::derivative() const {
    std::vector<double> ret_p(p_.size()-1);
    for (casadi_int k=1; k<p_.size(); ++k) {
      ret_p[k-1] = static_cast<double>(k)*p_[k];
    }
    return Polynomial(ret_p);
  }

  Polynomial Polynomial::anti_derivative() const {
    std::vector<double> ret_p(p_.size()+1);
    ret_p[0] = 0;
    for (casadi_int k=0; k<p_.size(); ++k) {
      ret_p[k+1] = p_[k]/static_cast<double>(k+1);
    }
    return Polynomial(ret_p);
  }


} // namespace casadi
