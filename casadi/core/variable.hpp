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


#ifndef CASADI_VARIABLE_HPP
#define CASADI_VARIABLE_HPP

#include "mx.hpp"
#include <iostream>

namespace casadi {
  /// Helper class: Specify number of entries in an enum
  template<typename T>
  struct enum_traits {};

  // Helper function: Convert string to enum
  template<typename T>
  T to_enum(const std::string& s) {
    // Linear search over permitted values
    for (size_t i = 0; i < enum_traits<T>::n_enum; ++i) {
      if (s == to_string(static_cast<T>(i))) return static_cast<T>(i);
    }
    // Informative error message
    std::stringstream ss;
    ss << "No such enum: '" << s << "'. Permitted values: ";
    for (size_t i = 0; i < enum_traits<T>::n_enum; ++i) {
      // Separate strings
      if (i > 0) ss << ", ";
      // Print enum name
      ss << "'" << to_string(static_cast<T>(i)) << "'";
    }
    casadi_error(ss.str());
    return enum_traits<T>::n_enum;  // never reached
  }

  /// Causality: FMI 2.0 specification, section 2.2.7
  enum Causality {PARAMETER, CALCULATED_PARAMETER, INPUT, OUTPUT, LOCAL, INDEPENDENT, N_CAUSALITY};

  /// Number of entries
  template<> struct enum_traits<Causality> {
    static const Causality n_enum = N_CAUSALITY;
  };

  /// Convert to string
  CASADI_EXPORT std::string to_string(Causality v);

  /// Variability: FMI 2.0 specification, section 2.2.7
  enum Variability {CONSTANT, FIXED, TUNABLE, DISCRETE, CONTINUOUS, N_VARIABILITY};

  /// Number of entries
  template<> struct enum_traits<Variability> {
    static const Variability n_enum = N_VARIABILITY;
  };

  /// Convert to string
  CASADI_EXPORT std::string to_string(Variability v);

  /// Initial: FMI 2.0 specification, section 2.2.7
  enum Initial {EXACT, APPROX, CALCULATED, INITIAL_NA, N_INITIAL};

  /// Number of entries
  template<> struct enum_traits<Initial> {
    static const Initial n_enum = N_INITIAL;
  };

  /// Convert to string
  CASADI_EXPORT std::string to_string(Initial v);

  /// Dynamics of the variable
  enum Dynamics {ALGEBRAIC, DIFFERENTIAL};

  /// Variable category
  enum Category {
    /** Unknown, not set */
    CAT_UNKNOWN,
    /** A state derivative */
    CAT_DERIVATIVE,
    /** A differential state, i.e. a variable that appears differentiated in the model */
    CAT_STATE,
    /** An independent constant: <tt>constant Real c1 = 3</tt> */
    CAT_DEPENDENT_CONSTANT,
    /** A dependent constant <tt>constant Real c2=c1*3</tt>. */
    CAT_INDEPENDENT_CONSTANT,
    /** A dependent parameter <tt>parameter Real p1=p2</tt>*/
    CAT_DEPENDENT_PARAMETER,
    /** An independent parameter <tt>parameter Real p2=3</tt>*/
    CAT_INDEPENDENT_PARAMETER,
    /** An algebraic variable or input */
    CAT_ALGEBRAIC
  };

  /** \brief Holds expressions and meta-data corresponding to a physical quantity evolving in time
      \date 2012-2014
      \author Joel Andersson
   */
  struct CASADI_EXPORT Variable :
    public SWIG_IF_ELSE(PrintableCommon, Printable<Variable>) {
    /// Default constructor
    Variable() {}

    /// Constructor
    explicit Variable(const std::string& name,
      const Sparsity& sp = Sparsity::scalar(),
      const MX& v = MX(),
      const MX& d = MX());

    /// Variable name
    std::string name() const;

    /// Variable expression
    MX v;

    /// Derivative expression
    MX d;

    /// Nominal value
    double nominal;

    /// Value at time 0
    double start;

    /// Lower bound
    double min;

    /// Upper bound
    double max;

    /// Initial guess
    double guess;

    /// Derivative at time 0
    double derivative_start;

    /// Variable reference (XML)
    casadi_int valueReference;

    /// Description
    std::string description;

    /// Causality
    Causality causality;

    /// Variability
    Variability variability;

    // Initial value
    Initial initial;

    /// Variable category
    Category category;

    /// Unit
    std::string unit;

    /// Display unit
    std::string display_unit;

    /// Free attribute
    bool free;

    /// Readable name of the class
    std::string type_name() const {return "Variable";}

    /// Print a description of the object
    void disp(std::ostream& stream, bool more=false) const;

    /// Get string representation
    std::string get_str(bool more=false) const {
      std::stringstream ss;
      disp(ss, more);
      return ss.str();
    }

    // Default initial attribute, per specification
    static Initial default_initial(Variability variability, Causality causality);
  };
} // namespace casadi

#endif // CASADI_VARIABLE_HPP
