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

#ifndef CONSTANT_NODE_HPP
#define CONSTANT_NODE_HPP

#include "sx_node.hpp"

namespace CasADi{

/** \brief  constants */
const double double_inf = std::numeric_limits<double>::infinity();
const double double_nan = std::numeric_limits<double>::quiet_NaN();

/** \brief Reprsents a constant SX
  \author Joel Andersson 
  \date 2010
*/
class ConstantSXNode : public SXNode{
public:

// Destructor
virtual ~ConstantSXNode(){};
  
/** \brief  Get the value must be defined */
virtual double getValue() const = 0;

/** \brief  Properties */
virtual bool isConstant() const{ return true; }

protected:

/** \brief  print */
virtual void print(std::ostream &stream) const{
  stream << getValue();
}

};

/** \brief  DERIVED CLASSES */

/** \brief  Represents a constant real SX
  \author Joel Andersson 
  \date 2010
*/
class RealtypeSXNode : public ConstantSXNode{
public:

// Destructor
virtual ~RealtypeSXNode(){}

explicit RealtypeSXNode(double value) : value(value){} 

//@{
/** \brief  Get the value */
virtual double getValue() const{ return value;}
virtual int getIntValue() const{ return int(value);}
//@}

protected:

/** \brief  Data members */
double value;

};


/** \brief  Represents a constant integer SX
  \author Joel Andersson 
  \date 2010
*/
class IntegerSXNode : public ConstantSXNode{
  public:

    virtual ~IntegerSXNode(){}
    explicit IntegerSXNode(int value) : value(value){}

    //@{
    /** \brief  evaluate function */
    virtual double getValue() const{  return value; }
    virtual int getIntValue() const{  return value; }
    //@}

    /** \brief  Properties */
    virtual bool isInteger() const{ return true; }
  
  protected:

    /** \brief  Data members */
    int value;
};


/** \brief  Represents a zero SX
  \author Joel Andersson 
  \date 2010
*/
class ZeroSXNode : public ConstantSXNode{
public:

  virtual ~ZeroSXNode(){}
  explicit ZeroSXNode(){}

  //@{
  /** \brief  Get the value */
  virtual double getValue() const{ return 0;}
  virtual int getIntValue() const{ return 0;}
  //@}

  //@{
  /** \brief  Properties */
  virtual bool isInteger() const{ return true; }
  virtual bool isZero() const{ return true; }
  //@}
};


/** \brief  Represents a one SX
  \author Joel Andersson 
  \date 2010
*/
class OneSXNode : public ConstantSXNode{
public:

  explicit OneSXNode(){}
  virtual ~OneSXNode(){}

  /** \brief  Get the value */
  virtual double getValue() const{ return 1;}
  virtual int getIntValue() const{ return 1;}

  /** \brief  Properties */
  virtual bool isInteger() const{ return true; }
  virtual bool isOne() const{ return true; }

};


/** \brief  Represents a minus one SX
  \author Joel Andersson 
  \date 2010
*/
class MinusOneSXNode : public ConstantSXNode{
public:

  explicit MinusOneSXNode(){}
  virtual ~MinusOneSXNode(){}

  //@{
  /** \brief  Get the value */
  virtual double getValue() const{ return -1;}
  virtual int getIntValue() const{ return -1;}
  //@}

  //@{
  /** \brief  Properties */
  virtual bool isInteger() const{ return true; }
  virtual bool isMinusOne() const{ return true; }
  //@}

};


/** \brief  Represents an infinity SX
  \author Joel Andersson 
  \date 2010
*/
class InfSXNode : public ConstantSXNode{
public:

  explicit InfSXNode(){}
  virtual ~InfSXNode(){}

  /** \brief  Get the value */
  virtual double getValue() const{ return double_inf;}

  /** \brief  Properties */
  virtual bool isInf() const{ return true; }

};


/** \brief  Represents a minus infinity SX
  \author Joel Andersson 
  \date 2010
*/
class MinusInfSXNode : public ConstantSXNode{
public:
  
  explicit MinusInfSXNode(){}
  virtual ~MinusInfSXNode(){}

  /** \brief  Get the value */
  virtual double getValue() const{ return -double_inf;}

  /** \brief  Properties */
  virtual bool isMinusInf() const{ return true; }

};


/** \brief  Represents a not-a-number SX
  \author Joel Andersson 
  \date 2010
*/
class NanSXNode : public ConstantSXNode{
public:
  
  explicit NanSXNode(){}
  virtual ~NanSXNode(){}

  /** \brief  Get the value */
  virtual double getValue() const{ return double_nan;}

  /** \brief  Properties */
  virtual bool isNan() const{ return true; }

};

} // namespace CasADi


#endif // CONSTANT_SCALAR_HPP
