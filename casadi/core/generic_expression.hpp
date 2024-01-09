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


#ifndef CASADI_GENERIC_EXPRESSION_HPP
#define CASADI_GENERIC_EXPRESSION_HPP

#include "calculus.hpp"

namespace casadi {

  /** \brief Empty Base

      This class is extended in SWIG.

      \identifier{ol} */
  struct CASADI_EXPORT GenericExpressionCommon {};

#ifndef SWIG
  /** \brief Expression interface
  *
  This is a common base class for SX, MX and Matrix<>, introducing a uniform syntax and implementing
  common functionality using the curiously recurring template pattern (CRTP) idiom.\n

  \author Joel Andersson
  \date 2012

      \identifier{om} */
template<typename ExType>
class GenericExpression : public GenericExpressionCommon {
  protected:
    // Helper functions
    inline const ExType& self() const { return static_cast<const ExType&>(*this); }
    inline ExType& self() { return static_cast<ExType&>(*this); }
  public:

/**
\addtogroup expression_tools
@{
*/

  ///@{
  /** \brief Addition: (x,y) -> x + y

      \identifier{on} */
  static ExType plus(const ExType &x, const ExType &y) {
    return ExType::binary(OP_ADD, x, y);
  }
  friend inline ExType plus(const ExType &x, const ExType &y) {
    return ExType::plus(x, y);
  }
  friend inline ExType operator+(const ExType &x, const ExType &y) {
    return plus(x, y);
  }
  inline ExType& operator+=(const ExType &y) { return self() = self() + y; }
  ///@}

  ///@{
  /** \brief Subtraction: (x,y) -> x - y

      \identifier{oo} */
  static ExType minus(const ExType &x, const ExType &y) {
    return ExType::binary(OP_SUB, x, y);
  }
  friend inline ExType minus(const ExType &x, const ExType &y) {
    return ExType::minus(x, y);
  }
  friend inline ExType operator-(const ExType &x, const ExType &y) {
    return minus(x, y);
  }
  inline ExType& operator-=(const ExType &y) { return self() = self() - y; }
  ///@}

  ///@{
  /** \brief Elementwise multiplication: (x,y) -> x .* y

      \identifier{op} */
  static ExType times(const ExType &x, const ExType &y) {
    return ExType::binary(OP_MUL, x, y);
  }
  friend inline ExType times(const ExType &x, const ExType &y) {
    return ExType::times(x, y);
  }
  friend inline ExType operator*(const ExType &x, const ExType &y) {
    return times(x, y);
  }
  inline ExType& operator*=(const ExType &y) {return self() = self() * y;}
  ///@}

  ///@{
  /** \brief Elementwise division: (x,y) -> x ./ y

      \identifier{oq} */
  static ExType rdivide(const ExType &x, const ExType &y) {
    return ExType::binary(OP_DIV, x, y);
  }
  friend inline ExType rdivide(const ExType &x, const ExType &y) {
    return ExType::rdivide(x, y);
  }
  friend inline ExType operator/(const ExType &x, const ExType &y) {
    return rdivide(x, y);
  }
  inline ExType& operator/=(const ExType &y) {return self() = self() / y;}
  ///@}

  ///@{
  /** \brief Logical less than: (x,y) -> x < y

      \identifier{or} */
  static ExType lt(const ExType &x, const ExType &y) {
    return ExType::binary(OP_LT, x, y);
  }
  friend inline ExType lt(const ExType &x, const ExType &y) {
    return ExType::lt(x, y);
  }
  friend inline ExType operator<(const ExType &x, const ExType &y) {
    return lt(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical less or equal to: (x,y) -> x <= y

      \identifier{os} */
  static ExType le(const ExType &x, const ExType &y) {
    return ExType::binary(OP_LE, x, y);
  }
  friend inline ExType le(const ExType &x, const ExType &y) {
    return ExType::le(x, y);
  }
  friend inline ExType operator<=(const ExType &x, const ExType &y) {
    return le(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical greater than: (x,y) -> x > y

      \identifier{ot} */
  static ExType gt(const ExType &x, const ExType &y) {
    return ExType::lt(y, x);
  }
  friend inline ExType gt(const ExType &x, const ExType &y) {
    return ExType::gt(x, y);
  }
  friend inline ExType operator>(const ExType &x, const ExType &y) {
    return gt(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical greater or equal to: (x,y) -> x >= y

      \identifier{ou} */
  static ExType ge(const ExType &x, const ExType &y) {
    return ExType::le(y, x);
  }
  friend inline ExType ge(const ExType &x, const ExType &y) {
    return ExType::ge(x, y);
  }
  friend inline ExType operator>=(const ExType &x, const ExType &y) {
    return ge(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical equal to: (x,y) -> x == y

      \identifier{ov} */
  static ExType eq(const ExType &x, const ExType &y) {
    return ExType::binary(OP_EQ, x, y);
  }
  friend inline ExType eq(const ExType &x, const ExType &y) {
    return ExType::eq(x, y);
  }
  friend inline ExType operator==(const ExType &x, const ExType &y) {
    return eq(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical not equal to: (x,y) -> x != y

      \identifier{ow} */
  static ExType ne(const ExType &x, const ExType &y) {
    return ExType::binary(OP_NE, x, y);
  }
  friend inline ExType ne(const ExType &x, const ExType &y) {
    return ExType::ne(x, y);
  }
  friend inline ExType operator!=(const ExType &x, const ExType &y) {
    return ne(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical `and`

   * Returns (an expression evaluating to) 1 if both
   * expressions are nonzero and 0 otherwise

      \identifier{ox} */
   static ExType logic_and(const ExType &x, const ExType &y) {
     return ExType::binary(OP_AND, x, y);
   }
   friend inline ExType logic_and(const ExType &x, const ExType &y) {
     return ExType::logic_and(x, y);
   }
   friend inline ExType operator&&(const ExType &x, const ExType &y) {
     return logic_and(x, y);
   }
  ///@}

  ///@{
  /** \brief  Logical `or`

   * returns (an expression evaluating to) 1 if at
   * least one expression is nonzero and 0 otherwise

      \identifier{oy} */
   static ExType logic_or(const ExType &x, const ExType &y) {
     return ExType::binary(OP_OR, x, y);
   }
   friend inline ExType logic_or(const ExType &x, const ExType &y) {
     return ExType::logic_or(x, y);
   }
   friend inline ExType operator||(const ExType &x, const ExType &y) {
     return logic_or(x, y);
   }
   ///@}

   ///@{
   /** \brief  Logical `not` x -> !x

    * Returns (an expression evaluating to) 1 if
    * expression is zero and 0 otherwise

       \identifier{oz} */
    static ExType logic_not(const ExType& x) {
      return ExType::unary(OP_NOT, x);
    }
    friend inline ExType logic_not(const ExType& x) {
      return ExType::logic_not(x);
    }
    inline ExType operator!() const {
      return logic_not(self());
    }
    ///@}

    ///@{
    /** \brief Absolute value: x -> abs(x)

        \identifier{p0} */
    static ExType abs(const ExType& x) {
      return ExType::unary(OP_FABS, x);
    }
    friend inline ExType abs(const ExType& x) {
      return ExType::abs(x);
    }
    friend inline ExType fabs(const ExType& x) {
      return abs(x);
    }
    ///@}

    ///@{
    /** \brief Square root: x -> sqrt(x)

        \identifier{p1} */
    static ExType sqrt(const ExType& x) {
      return ExType::unary(OP_SQRT, x);
    }
    friend inline ExType sqrt(const ExType& x) {
      return ExType::sqrt(x);
    }
    ///@}

    ///@{
    /** \brief Square: x -> x^2

        \identifier{p2} */
    static ExType sq(const ExType& x) {
      return ExType::unary(OP_SQ, x);
    }
    friend inline ExType sq(const ExType& x) {
      return ExType::sq(x);
    }
    ///@}

    ///@{
    /** \brief Sine: x -> sin(x)

        \identifier{p3} */
    static ExType sin(const ExType& x) {
      return ExType::unary(OP_SIN, x);
    }
    friend inline ExType sin(const ExType& x) {
      return ExType::sin(x);
    }
    ///@}

    ///@{
    /** \brief Cosine: x -> cos(x)

        \identifier{p4} */
    static ExType cos(const ExType& x) {
      return ExType::unary(OP_COS, x);
    }
    friend inline ExType cos(const ExType& x) {
      return ExType::cos(x);
    }
    ///@}

    ///@{
    /** \brief Tangent: x -> tan(x)

        \identifier{p5} */
    static ExType tan(const ExType& x) {
      return ExType::unary(OP_TAN, x);
    }
    friend inline ExType tan(const ExType& x) {
      return ExType::tan(x);
    }
    ///@}

    ///@{
    /** \brief Arc tangent: x -> atan(x)

        \identifier{p6} */
    static ExType atan(const ExType& x) {
      return ExType::unary(OP_ATAN, x);
    }
    friend inline ExType atan(const ExType& x) {
      return ExType::atan(x);
    }
    ///@}

    ///@{
    /** \brief Arc sine: x -> asin(x)

        \identifier{p7} */
    static ExType asin(const ExType& x) {
      return ExType::unary(OP_ASIN, x);
    }
    friend inline ExType asin(const ExType& x) {
      return ExType::asin(x);
    }
    ///@}

    ///@{
    /** \brief Arc cosine: x -> acos(x)

        \identifier{p8} */
    static ExType acos(const ExType& x) {
      return ExType::unary(OP_ACOS, x);
    }
    friend inline ExType acos(const ExType& x) {
      return ExType::acos(x);
    }
    ///@}

    ///@{
    /** \brief Hyperbolic tangent: x -> tanh(x)

        \identifier{p9} */
    static ExType tanh(const ExType& x) {
      return ExType::unary(OP_TANH, x);
    }
    friend inline ExType tanh(const ExType& x) {
      return ExType::tanh(x);
    }
    ///@}

    ///@{
    /** \brief Hyperbolic sin: x -> sinh(x)

        \identifier{pa} */
    static ExType sinh(const ExType& x) {
      return ExType::unary(OP_SINH, x);
    }
    friend inline ExType sinh(const ExType& x) {
      return ExType::sinh(x);
    }
    ///@}

    ///@{
    /** \brief Hyperbolic cosine: x -> cosh(x)

        \identifier{pb} */
    static ExType cosh(const ExType& x) {
      return ExType::unary(OP_COSH, x);
    }
    friend inline ExType cosh(const ExType& x) {
      return ExType::cosh(x);
    }
    ///@}

    ///@{
    /** \brief Inverse hyperbolic tangent: x -> atanh(x)

        \identifier{pc} */
    static ExType atanh(const ExType& x) {
      return ExType::unary(OP_ATANH, x);
    }
    friend inline ExType atanh(const ExType& x) {
      return ExType::atanh(x);
    }
    ///@}

    ///@{
    /** \brief Inverse hyperbolic sin: x -> asinh(x)

        \identifier{pd} */
    static ExType asinh(const ExType& x) {
      return ExType::unary(OP_ASINH, x);
    }
    friend inline ExType asinh(const ExType& x) {
      return ExType::asinh(x);
    }
    ///@}

    ///@{
    /** \brief Inverse hyperbolic cosine: x -> acosh(x)

        \identifier{pe} */
    static ExType acosh(const ExType& x) {
      return ExType::unary(OP_ACOSH, x);
    }
    friend inline ExType acosh(const ExType& x) {
      return ExType::acosh(x);
    }
    ///@}

    ///@{
    /** \brief Elementwise exponential: x -> exp(x)

        \identifier{pf} */
    static ExType exp(const ExType& x) {
      return ExType::unary(OP_EXP, x);
    }
    friend inline ExType exp(const ExType& x) {
      return ExType::exp(x);
    }
    ///@}

    ///@{
    /** \brief Natural logarithm: x -> log(x)

        \identifier{pg} */
    static ExType log(const ExType& x) {
      return ExType::unary(OP_LOG, x);
    }
    friend inline ExType log(const ExType& x) {
      return ExType::log(x);
    }
    ///@}

    ///@{
    /** \brief Base-10 logarithm: x -> log10(x)

        \identifier{ph} */
    static ExType log10(const ExType& x) {
      return log(x)*(1/std::log(10.));
    }
    friend inline ExType log10(const ExType& x) {
      return ExType::log10(x);
    }
    ///@}

    ///@{
    /** \brief Precision variant for natural logarithm: x -> log(x+1)

        \identifier{pi} */
    static ExType log1p(const ExType& x) {
      return ExType::unary(OP_LOG1P, x);
    }
    friend inline ExType log1p(const ExType& x) {
      return ExType::log1p(x);
    }
    ///@}

    ///@{
    /** \brief Precision variant for elementwise exponential: x -> exp(x)-1

        \identifier{pj} */
    static ExType expm1(const ExType& x) {
      return ExType::unary(OP_EXPM1, x);
    }
    friend inline ExType expm1(const ExType& x) {
      return ExType::expm1(x);
    }
    ///@}

    ///@{
    /** \brief Round down to nearest integer: x -> floor(x)

        \identifier{pk} */
    static ExType floor(const ExType& x) {
      return ExType::unary(OP_FLOOR, x);
    }
    friend inline ExType floor(const ExType& x) {
      return ExType::floor(x);
    }
    ///@}

    ///@{
    /** \brief Round up to nearest integer: x -> ceil(x)

        \identifier{pl} */
    static ExType ceil(const ExType& x) {
      return ExType::unary(OP_CEIL, x);
    }
    friend inline ExType ceil(const ExType& x) {
      return ExType::ceil(x);
    }
    ///@}

    ///@{
    /** \brief Error function: x -> erf(x)

        \identifier{pm} */
    static ExType erf(const ExType& x) {
      return ExType::unary(OP_ERF, x);
    }
    friend inline ExType erf(const ExType& x) {
      return ExType::erf(x);
    }
    ///@}

    ///@{
    /** \brief Inverse error function: x -> erfinv(x)

        \identifier{pn} */
    static ExType erfinv(const ExType& x) {
      return ExType::unary(OP_ERFINV, x);
    }
    friend inline ExType erfinv(const ExType& x) {
      return ExType::erfinv(x);
    }
    ///@}

    ///@{
    /** \brief Sign function:

        sign(x)   := -1 for x<0
        sign(x)   :=  1 for x>0,
        sign(0)   :=  0
        sign(NaN) :=  NaN

        \identifier{po} */
    static ExType sign(const ExType& x) {
      return ExType::unary(OP_SIGN, x);
    }
    friend inline ExType sign(const ExType& x) {
      return ExType::sign(x);
    }
    ///@}

    ///@{
    /** \brief Elementwise power: (x,y) -> x.^y

        \identifier{pp} */
    static ExType pow(const ExType& x, const ExType& y) {
      return ExType::binary(OP_POW, x, y);
    }
    friend inline ExType pow(const ExType& x, const ExType& y) {
      return ExType::pow(x, y);
    }
    ///@}

    ///@{
    /** \brief Remainder after division: (x,y) -> fmod(x,y)

    This Function follows the convention of https://en.cppreference.com/w/c/numeric/math/fmod

    Notably:
      - fmod(5,3)   -> 2
      - fmod(5,-3)  -> 2
      - fmod(-5,3)  -> -2
      - fmod(-5,-3) -> -2

    This is equivalent to Python's numpy.fmod and Matlab's rem.

    \seealso remainder

        \identifier{pq} */
    static ExType mod(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMOD, x, y);
    }
    friend inline ExType mod(const ExType& x, const ExType& y) {
      return ExType::mod(x, y);
    }
    friend inline ExType fmod(const ExType& x, const ExType& y) {
      return mod(x, y);
    }
    ///@}

    ///@{
    /** \brief Remainder after division: (x,y) -> remainder(x,y)

    This Function follows the convention of https://en.cppreference.com/w/c/numeric/math/remainder

    Notably:
      - remainder(5,3)   -> -1
      - remainder(5,-3)  -> -1
      - remainder(-5,3)  -> 1
      - remainder(-5,-3) -> 1

    This is equivalent to Python's math.remainder. There is no equivalence in Matlab.

    \seealso fmod

        \identifier{24x} */
    static ExType remainder(const ExType& x, const ExType& y) {
      return ExType::binary(OP_REMAINDER, x, y);
    }
    friend inline ExType remainder(const ExType& x, const ExType& y) {
      return ExType::remainder(x, y);
    }
    ///@}

    ///@{
    /** \brief Two argument arc tangent: (y,x) -> atan2(y,x)
     *
     * theta = atan2(y,x) corresponds to x = r cos(theta), y = r sin(theta)

        \identifier{pr} */
    static ExType atan2(const ExType& y, const ExType& x) {
      return ExType::binary(OP_ATAN2, y, x);
    }
    friend inline ExType atan2(const ExType& y, const ExType& x) {
      return ExType::atan2(y, x);
    }
    ///@}

    ///@{
    /** \brief Conditional assignment: (x,y) -> x ? y : 0

        \identifier{ps} */
    static ExType if_else_zero(const ExType& x, const ExType& y) {
      return ExType::binary(OP_IF_ELSE_ZERO, x, y);
    }
    friend inline ExType if_else_zero(const ExType& x, const ExType& y) {
      return ExType::if_else_zero(x, y);
    }
    ///@}

    ///@{
    /** \brief Smallest of two values: (x,y) -> min(x,y)

        \identifier{pt} */
    static ExType fmin(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMIN, x, y);
    }
    friend inline ExType fmin(const ExType& x, const ExType& y) {
      return ExType::fmin(x, y);
    }
    ///@}

    ///@{
    /** \brief Largest of two values: (x,y) -> max(x,y)

        \identifier{pu} */
    static ExType fmax(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMAX, x, y);
    }
    friend inline ExType fmax(const ExType& x, const ExType& y) {
      return ExType::fmax(x, y);
    }
    ///@}

    ///@{
    /** \brief Check if two nodes are equivalent up to a given depth.

     * Depth=0 checks if the expressions are identical, i.e. points to the same node.
     *
     * a = x*x
     * b = x*x
     *
     *  is_equal(a,b,0)  will return false, but a.is_equal(a,b,1) will return true

        \identifier{pv} */
     friend inline bool is_equal(const ExType& x, const ExType& y, casadi_int depth=0) {
       return ExType::is_equal(x, y, depth);
     }
     ///@}

     ///@{
     /// Copy sign
     static ExType copysign(const ExType& x, const ExType& y) {
       return ExType::binary(OP_COPYSIGN, x, y);
     }
     friend inline ExType copysign(const ExType& x, const ExType& y) {
       return ExType::copysign(x, y);
     }
     ///@}

     ///@{
     /// Elementwise power with const power
     static ExType constpow(const ExType& x, const ExType& y) {
       return ExType::binary(OP_CONSTPOW, x, y);
     }
     friend inline ExType constpow(const ExType& x, const ExType& y) {
       return ExType::constpow(x, y);
     }
     ///@}

     ///@{
     /// Debug printing
     static ExType printme(const ExType& x, const ExType& y) {
       return ExType::binary(OP_PRINTME, x, y);
     }
     friend inline ExType printme(const ExType& x, const ExType& y) {
       return ExType::binary(OP_PRINTME, x, y);
     }
     ///@}

    ///@{
    /** \brief Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

        \identifier{pw} */
    static ExType hypot(const ExType& x, const ExType& y) {
      return ExType::binary(OP_HYPOT, x, y);
    }
    friend inline ExType hypot(const ExType& x, const ExType& y) {
      return ExType::hypot(x, y);
    }
    ///@}


/** @} */



};
#endif // SWIG

} // namespace casadi

#endif // CASADI_GENERIC_EXPRESSION_HPP
