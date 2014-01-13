/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

# include "printable_object.hpp"

namespace CasADi{

  /** \brief Helper class for differentiating and integrating polynomials
      \author Joel Andersson 
      \date 2014
  */
  class Polynomial : public PrintableObject{
  public:
    
    /// Construct a constant polynomial
    Polynomial(double scalar=1);

    /// Construct a linear polynomial
    Polynomial(double p0, double p1);

    /// Construct a quadratic polynomial
    Polynomial(double p0, double p1, double p2);

    /// Construct a cubic polynomial
    Polynomial(double p0, double p1, double p2, double p3);

    /// Construct from a vector of polynomial coefficients
    Polynomial(const std::vector<double>& coeff);

    /// Evaluate numerically
    template<typename T>
    T operator()(const T& x) const{
      std::vector<double>::const_reverse_iterator it = p_.rbegin();
      T ret = *it++;
      while(it!=p_.rend()){
        ret *= x;
        ret += *it++;
      }
      return ret;
    }
    
    /// Degree of the polynomial
    int degree() const;

    /// Get scalar value (error if degree()!=0)
    double toScalar() const;

    /// Create a new polynomial for the derivative
    Polynomial derivative() const;

    /// Create a new polynomial for the anti-derivative (primitive function)
    Polynomial anti_derivative() const;
    
    /// Remove excess zeros
    void trim();

#ifndef SWIG
    /// Print a destription of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;
#endif // SWIG

    // Add
    Polynomial operator+(const Polynomial& b) const;

    // Add (in-place)
    Polynomial& operator+=(const Polynomial& b);

    // Subtract
    Polynomial operator-(const Polynomial& b) const;

    // Subtract (in-place)
    Polynomial& operator-=(const Polynomial& b);

    // Multiply
    Polynomial operator*(const Polynomial& b) const;

    // Multiply (in-place)
    Polynomial& operator*=(const Polynomial& b);

    // Divide by constant
    Polynomial operator/(double b) const;

    // Divide by constant (in-place)
    Polynomial& operator/=(double b);


  protected:
    std::vector<double> p_;
  };

} // namespace CasADi


#endif // POLYNOMIAL_HPP

