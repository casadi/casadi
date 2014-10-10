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


#ifndef CASADI_CONSTANT_SX_HPP
#define CASADI_CONSTANT_SX_HPP

#include "sx_node.hpp"
#include <cassert>

/// \cond INTERNAL

// Cashing of constants requires a map (preferably a hash map)
#ifdef USE_CXX11
// Using C++11 unordered_map (hash map)
#ifdef USE_TR1_HASHMAP
#include <tr1/unordered_map>
#define CACHING_MAP std::tr1::unordered_map
#else // USE_TR1_HASHMAP
#include <unordered_map>
#define CACHING_MAP std::unordered_map
#endif // USE_TR1_HASHMAP
#else // USE_CXX11
// Falling back to std::map (binary search tree)
#include <map>
#define CACHING_MAP std::map
#endif // USE_CXX11

namespace casadi {

/** \brief Represents a constant SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT ConstantSX : public SXNode {
public:

// Destructor
virtual ~ConstantSX() {}

/** \brief  Get the value must be defined */
virtual double getValue() const = 0;

/** \brief  Properties */
virtual bool isConstant() const { return true; }

/** \brief  Get the operation */
virtual int getOp() const { return OP_CONST;}

/** \brief Check if two nodes are equivalent up to a given depth */
virtual bool isEqual(const SXNode* node, int depth) const {
  const ConstantSX* n = dynamic_cast<const ConstantSX*>(node);
  return n && n->getValue()==getValue();
}

protected:

/** \brief  print */
virtual void print(std::ostream &stream, long& remaining_calls) const {
  stream << getValue();
}

};

/** \brief  DERIVED CLASSES */

/** \brief  Represents a constant real SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT RealtypeSX : public ConstantSX {
  private:
    /// Constructor is private, use "create" below
    explicit RealtypeSX(double value) : value(value) {}

  public:

    /// Destructor
    virtual ~RealtypeSX() {
      size_t num_erased = cached_constants_.erase(value);
      assert(num_erased==1);
      (void)num_erased;
    }

    /// Static creator function (use instead of constructor)
    inline static RealtypeSX* create(double value) {
      // Try to find the constant
      CACHING_MAP<double, RealtypeSX*>::iterator it = cached_constants_.find(value);

      // If not found, add it,
      if (it==cached_constants_.end()) {
        // Allocate a new object
        RealtypeSX* n = new RealtypeSX(value);

        // Add to hash_table
        cached_constants_.insert(it, std::make_pair(value, n));

        // Return it to caller
        return n;
      } else { // Else, returned the object
        return it->second;
      }
    }

    ///@{
    /** \brief  Get the value */
    virtual double getValue() const { return value;}
    virtual int getIntValue() const { return static_cast<int>(value);}
    ///@}

    virtual bool isAlmostZero(double tol) const { return fabs(value)<=tol; }

  protected:
    /** \brief Hash map of all constants currently allocated
     * (storage is allocated for it in sx_element.cpp) */
    static CACHING_MAP<double, RealtypeSX*> cached_constants_;

    /** \brief  Data members */
    double value;
};


/** \brief  Represents a constant integer SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT IntegerSX : public ConstantSX {
  private:
    /// Constructor is private, use "create" below
    explicit IntegerSX(int value) : value(value) {}

  public:

    /// Destructor
    virtual ~IntegerSX() {
      size_t num_erased = cached_constants_.erase(value);
      assert(num_erased==1);
      (void)num_erased;
    }

    /// Static creator function (use instead of constructor)
    inline static IntegerSX* create(int value) {
      // Try to find the constant
      CACHING_MAP<int, IntegerSX*>::iterator it = cached_constants_.find(value);

      // If not found, add it,
      if (it==cached_constants_.end()) {
        // Allocate a new object
        IntegerSX* n = new IntegerSX(value);

        // Add to hash_table
        cached_constants_.insert(it, std::make_pair(value, n));

        // Return it to caller
        return n;
      } else { // Else, returned the object
        return it->second;
      }
    }

    ///@{
    /** \brief  evaluate function */
    virtual double getValue() const {  return value; }
    virtual int getIntValue() const {  return value; }
    ///@}

    /** \brief  Properties */
    virtual bool isInteger() const { return true; }

  protected:

    /** \brief Hash map of all constants currently allocated
     * (storage is allocated for it in sx_element.cpp) */
    static CACHING_MAP<int, IntegerSX*> cached_constants_;

    /** \brief  Data members */
    int value;
};

/** \brief  Represents a zero SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT ZeroSX : public ConstantSX {
public:

  virtual ~ZeroSX() {}
  explicit ZeroSX() {}

  ///@{
  /** \brief  Get the value */
  virtual double getValue() const { return 0;}
  virtual int getIntValue() const { return 0;}
  ///@}

  ///@{
  /** \brief  Properties */
  virtual bool isInteger() const { return true; }
  virtual bool isZero() const { return true; }
  virtual bool isAlmostZero(double tol) const { return true; }
  ///@}
};


/** \brief  Represents a one SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT OneSX : public ConstantSX {
public:

  explicit OneSX() {}
  virtual ~OneSX() {}

  /** \brief  Get the value */
  virtual double getValue() const { return 1;}
  virtual int getIntValue() const { return 1;}

  /** \brief  Properties */
  virtual bool isInteger() const { return true; }
  virtual bool isOne() const { return true; }

};


/** \brief  Represents a minus one SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT MinusOneSX : public ConstantSX {
public:

  explicit MinusOneSX() {}
  virtual ~MinusOneSX() {}

  ///@{
  /** \brief  Get the value */
  virtual double getValue() const { return -1;}
  virtual int getIntValue() const { return -1;}
  ///@}

  ///@{
  /** \brief  Properties */
  virtual bool isInteger() const { return true; }
  virtual bool isMinusOne() const { return true; }
  ///@}

};


/** \brief  Represents an infinity SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT InfSX : public ConstantSX {
public:

  explicit InfSX() {}
  virtual ~InfSX() {}

  /** \brief  Get the value */
  virtual double getValue() const { return std::numeric_limits<double>::infinity();}

  /** \brief  Properties */
  virtual bool isInf() const { return true; }

};


/** \brief  Represents a minus infinity SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT MinusInfSX : public ConstantSX {
public:

  explicit MinusInfSX() {}
  virtual ~MinusInfSX() {}

  /** \brief  Get the value */
  virtual double getValue() const { return -std::numeric_limits<double>::infinity();}

  /** \brief  Properties */
  virtual bool isMinusInf() const { return true; }

};


/** \brief  Represents a not-a-number SX
  \author Joel Andersson
  \date 2010
*/
class CASADI_CORE_EXPORT NanSX : public ConstantSX {
public:

  explicit NanSX() {this->count++;}
  virtual ~NanSX() {this->count--;}

  /** \brief  Get the value */
  virtual double getValue() const { return std::numeric_limits<double>::quiet_NaN();}

  /** \brief  Properties */
  virtual bool isNan() const { return true; }

};

} // namespace casadi
/// \endcond

#endif // CASADI_CONSTANT_SX_HPP
