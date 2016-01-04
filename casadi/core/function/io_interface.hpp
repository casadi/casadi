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


#ifndef CASADI_IO_INTERFACE_HPP
#define CASADI_IO_INTERFACE_HPP

#include "../sx/sx_elem.hpp"
#include "../mx/mx.hpp"
#include "../options_functionality.hpp"

#include <exception>

namespace casadi {

  /** \brief Interface for accessing input and output data structures
      \author Joel Andersson
      \date 2013
  */
  template<class Derived>
  class CASADI_EXPORT IOInterface {
#ifndef SWIG
  protected:
    // Helper functions
    inline const Derived& self() const { return static_cast<const Derived&>(*this); }
    inline Derived& self() { return static_cast<Derived&>(*this); }
#endif // SWIG
  public:
    /**
    * \defgroup  iname
    *   \param[in] iname input name. Only allowed when an input scheme is set.
    */

    /**
    * \defgroup  oname
    *   \param[in] oname output name. Only allowed when an output scheme is set.
    */

    /**
    * \defgroup  iind
    *   \param[in] iind index within the range [0..n_in()-1]
    */

    /**
    * \defgroup  oind
    *   \param[in] oind index within the range [0..n_out()-1]
    */

    /**
    * \defgroup  Tvalset
    *   \param[in] val can be double, const std::vector<double>&, const Matrix<double>&, double *
    */

    /**
    * \defgroup  Tvalget
    *   \param[in] val can be double&, std::vector<double>&, Matrix<double>&, double *
    */

    /// \name Simple Getters & Setters
    ///
    /// @{
    /** \brief Get an input by index
    *
    *  @copydoc iind
    */
    Matrix<double> getInput(int iind=0) const { return self().input(iind); }
    /** \brief Get an input by name
    *
    *  @copydoc iname
    *
     */
    Matrix<double> getInput(const std::string &iname) const  { return self().input(iname); }

    /** \brief Get an output by index
    *
    *  @copydoc oind
    */
    Matrix<double> getOutput(int oind=0) const  { return self().output(oind); }
    /** \brief Get an output by name
    *
    *  @copydoc oname
    *
     */
    Matrix<double> getOutput(const std::string &oname) const  { return self().output(oname); }

#ifndef DOXYGENPROC
    void setInputNZ(const Matrix<double>& val, int iind=0) {
      try {
        self().input(iind).setNZ(val);
      } catch(std::exception& e) {
        casadi_error(e.what() << "Occurred at iind = " << iind << ".");
      }
    }
    void setOutputNZ(const Matrix<double>& val, int oind=0) {
      self().output(oind).setNZ(val);
    }
    void setInputNZ(const Matrix<double>& val, const std::string &iname) {
      setInputNZ(val, self().index_in(iname));
    }
    void setOutputNZ(const Matrix<double>& val, const std::string &oname) {
      setOutputNZ(val, self().index_out(oname));
    }

    void setInput(const Matrix<double>& val, int iind=0) {
      try {
        self().input(iind).set(val);
      } catch(std::exception& e) {
        casadi_error(e.what() << "Occurred at iind = " << iind << ".");
      }
    }
    void setOutput(const Matrix<double>& val, int oind=0) {
      self().output(oind).set(val);
    }
    void setInput(const Matrix<double>& val, const std::string &iname) {
      setInput(val, self().index_in(iname));
    }
    void setOutput(const Matrix<double>& val, const std::string &oname) {
      setOutput(val, self().index_out(oname));
    }
#endif // DOXYGENPROC
  };


} // namespace casadi


#endif // CASADI_IO_INTERFACE_HPP
