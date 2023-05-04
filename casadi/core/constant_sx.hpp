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


#ifndef CASADI_CONSTANT_SX_HPP
#define CASADI_CONSTANT_SX_HPP

#include "sx_node.hpp"
#include "serializing_stream.hpp"
#include <cassert>

/// \cond INTERNAL

// Cashing of constants requires a map
#include <unordered_map>
#define CACHING_MAP std::unordered_map

namespace casadi {

/** \brief Represents a constant SX

  \author Joel Andersson
  \date 2010

    \identifier{1jj} */
class ConstantSX : public SXNode {
public:

// Destructor
~ConstantSX() override {}

// Class name
std::string class_name() const override {return "ConstantSX";}

/** \brief  Get the value must be defined

    \identifier{1jk} */
double to_double() const override = 0;

/** \brief  Properties

    \identifier{1jl} */
bool is_constant() const override { return true; }

/** \brief  Get the operation

    \identifier{1jm} */
casadi_int op() const override { return OP_CONST;}

/** \brief Check if two nodes are equivalent up to a given depth

    \identifier{1jn} */
bool is_equal(const SXNode* node, casadi_int depth) const override {
  const ConstantSX* n = dynamic_cast<const ConstantSX*>(node);
  return n && n->to_double()==to_double();
}

protected:

/** \brief  Print expression

    \identifier{1jo} */
std::string print(const std::string& arg1, const std::string& arg2) const override {
   std::stringstream ss;
   ss << to_double();
   return ss.str();
 }

};

/** \brief  DERIVED CLASSES

    \identifier{1jp} */

/** \brief  Represents a constant real SX

  \author Joel Andersson
  \date 2010

    \identifier{1jq} */
class RealtypeSX : public ConstantSX {
  private:
    /// Constructor is private, use "create" below
    explicit RealtypeSX(double value) : value(value) {}

  public:

    /// Destructor
    ~RealtypeSX() override {
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
    /** \brief  Get the value

        \identifier{1jr} */
    double to_double() const override { return value;}
    casadi_int to_int() const override { return static_cast<casadi_int>(value);}
    ///@}

    bool is_almost_zero(double tol) const override { return fabs(value)<=tol; }

    void serialize_node(SerializingStream& s) const override {
      s.pack("ConstantSX::type", 'r');
      s.pack("ConstantSX::value", value);
    }

  protected:
    /** \brief Hash map of all constants currently allocated

     * (storage is allocated for it in sx_element.cpp)

        \identifier{1js} */
    static CACHING_MAP<double, RealtypeSX*> cached_constants_;

    /** \brief  Data members

        \identifier{1jt} */
    double value;
};


/** \brief  Represents a constant integer SX

  \author Joel Andersson
  \date 2010

    \identifier{1ju} */
class IntegerSX : public ConstantSX {
  private:
    /// Constructor is private, use "create" below
    explicit IntegerSX(casadi_int value) : value(static_cast<int>(value)) {
      casadi_assert(value<=std::numeric_limits<int>::max() &&
                    value>=std::numeric_limits<int>::min(), "Integer overflow");
    }

  public:

    /// Destructor
    ~IntegerSX() override {
      size_t num_erased = cached_constants_.erase(value);
      assert(num_erased==1);
      (void)num_erased;
    }

    /// Static creator function (use instead of constructor)
    inline static IntegerSX* create(casadi_int value) {
      // Try to find the constant
      CACHING_MAP<casadi_int, IntegerSX*>::iterator it = cached_constants_.find(value);

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
    /** \brief  evaluate function

        \identifier{1jv} */
    double to_double() const override {  return static_cast<double>(value); }
    casadi_int to_int() const override {  return static_cast<casadi_int>(value); }
    ///@}

    /** \brief  Properties

        \identifier{1jw} */
    bool is_integer() const override { return true; }

    void serialize_node(SerializingStream& s) const override {
      s.pack("ConstantSX::type", 'i');
      s.pack("ConstantSX::value", value);
    }

  protected:

    /** \brief Hash map of all constants currently allocated

     * (storage is allocated for it in sx_element.cpp)

        \identifier{1jx} */
    static CACHING_MAP<casadi_int, IntegerSX*> cached_constants_;

    /** \brief  Data members

        \identifier{1jy} */
    int value;
};

/** \brief  Represents a zero SX

  \author Joel Andersson
  \date 2010

    \identifier{1jz} */
class ZeroSX : public ConstantSX {
private:
  /* Private constructor (singleton class) */
  explicit ZeroSX() {this->count++;}
public:
  /* Get singleton instance */
  static ZeroSX* singleton() {
    static ZeroSX instance;
    return &instance;
  }
  /* Destructor */
  ~ZeroSX() override {this->count--;}
  ///@{
  /** \brief  Get the value

      \identifier{1k0} */
  double to_double() const override { return 0;}
  casadi_int to_int() const override { return 0;}
  ///@}

  ///@{
  /** \brief  Properties

      \identifier{1k1} */
  bool is_integer() const override { return true; }
  bool is_zero() const override { return true; }
  bool is_almost_zero(double tol) const override { return true; }
  ///@}

  void serialize_node(SerializingStream& s) const override {
    s.pack("ConstantSX::type", '0');
  }
};


/** \brief  Represents a one SX

  \author Joel Andersson
  \date 2010

    \identifier{1k2} */
class OneSX : public ConstantSX {
private:
  /* Private constructor (singleton class) */
  explicit OneSX() {this->count++;}
public:
  /* Get singleton instance */
  static OneSX* singleton() {
    static OneSX instance;
    return &instance;
  }
  /* Destructor */
  ~OneSX() override {this->count--;}
  /** \brief  Get the value

      \identifier{1k3} */
  double to_double() const override { return 1;}
  casadi_int to_int() const override { return 1;}

  /** \brief  Properties

      \identifier{1k4} */
  bool is_integer() const override { return true; }
  bool is_one() const override { return true; }

  void serialize_node(SerializingStream& s) const override {
    s.pack("ConstantSX::type", '1');
  }

};


/** \brief  Represents a minus one SX

  \author Joel Andersson
  \date 2010

    \identifier{1k5} */
class MinusOneSX : public ConstantSX {
private:
  /* Private constructor (singleton class) */
  explicit MinusOneSX() {this->count++;}
public:
  /* Get singleton instance */
  static MinusOneSX* singleton() {
    static MinusOneSX instance;
    return &instance;
  }
  /* Destructor */
  ~MinusOneSX() override {this->count--;}

  ///@{
  /** \brief  Get the value

      \identifier{1k6} */
  double to_double() const override { return -1;}
  casadi_int to_int() const override { return -1;}
  ///@}

  ///@{
  /** \brief  Properties

      \identifier{1k7} */
  bool is_integer() const override { return true; }
  bool is_minus_one() const override { return true; }
  ///@}

  void serialize_node(SerializingStream& s) const override {
    s.pack("ConstantSX::type", 'm');
  }

};


/** \brief  Represents an infinity SX

  \author Joel Andersson
  \date 2010

    \identifier{1k8} */
class InfSX : public ConstantSX {
private:
  /* Private constructor (singleton class) */
  explicit InfSX() {this->count++;}
public:
  /* Get singleton instance */
  static InfSX* singleton() {
    static InfSX instance;
    return &instance;
  }
  /* Destructor */
  ~InfSX() override {this->count--;}
  /** \brief  Get the value

      \identifier{1k9} */
  double to_double() const override { return std::numeric_limits<double>::infinity();}

  /** \brief  Properties

      \identifier{1ka} */
  bool is_inf() const override { return true; }

  void serialize_node(SerializingStream& s) const override {
    s.pack("ConstantSX::type", 'F');
  }

};


/** \brief  Represents a minus infinity SX

  \author Joel Andersson
  \date 2010

    \identifier{1kb} */
class MinusInfSX : public ConstantSX {
private:
  /* Private constructor (singleton class) */
  explicit MinusInfSX() {this->count++;}
public:
  /* Get singleton instance */
  static MinusInfSX* singleton() {
    static MinusInfSX instance;
    return &instance;
  }
  /* Destructor */
  ~MinusInfSX() override {this->count--;}

  /** \brief  Get the value

      \identifier{1kc} */
  double to_double() const override { return -std::numeric_limits<double>::infinity();}

  /** \brief  Properties

      \identifier{1kd} */
  bool is_minus_inf() const override { return true; }

  void serialize_node(SerializingStream& s) const override {
    s.pack("ConstantSX::type", 'f');
  }
};

/** \brief  Represents a not-a-number SX

  \author Joel Andersson
  \date 2010

    \identifier{1ke} */
class NanSX : public ConstantSX {
private:
  /* Private constructor (singleton class) */
  explicit NanSX() {this->count++;}
public:
  /* Get singleton instance */
  static NanSX* singleton() {
    static NanSX instance;
    return &instance;
  }
  /* Destructor */
  ~NanSX() override {this->count--;}
  /** \brief  Get the value

      \identifier{1kf} */
  double to_double() const override { return std::numeric_limits<double>::quiet_NaN();}

  /** \brief  Properties

      \identifier{1kg} */
  bool is_nan() const override { return true; }

  void serialize_node(SerializingStream& s) const override {
    s.pack("ConstantSX::type", 'n');
  }

};

inline SXNode* ConstantSX_deserialize(DeserializingStream& s) {
  char type;
  s.unpack("ConstantSX::type", type);
  switch (type) {
    case '1': return casadi_limits<SXElem>::one.get();
    case '0': return casadi_limits<SXElem>::zero.get();
    case 'r': {
      double value;
      s.unpack("ConstantSX::value", value);
      return RealtypeSX::create(value);
    }
    case 'i': {
      int value;
      s.unpack("ConstantSX::value", value);
      if (value==2) return casadi_limits<SXElem>::two.get();
      return IntegerSX::create(value);
    }
    case 'n': return casadi_limits<SXElem>::nan.get();
    case 'f': return casadi_limits<SXElem>::minus_inf.get();
    case 'F': return casadi_limits<SXElem>::inf.get();
    case 'm': return casadi_limits<SXElem>::minus_one.get();
    default: casadi_error("ConstantSX::deserialize error");
  }
}

} // namespace casadi
/// \endcond

#endif // CASADI_CONSTANT_SX_HPP
