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

#include "../sx/sx_element.hpp"
#include "../mx/mx.hpp"
#include "../options_functionality.hpp"
#include "../function/schemes_helpers.hpp"
#include "../function/io_scheme.hpp"

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
    *   \param[in] iind index within the range [0..nIn()-1]
    */

    /**
    * \defgroup  oind
    *   \param[in] oind index within the range [0..nOut()-1]
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

#ifdef DOXYGENPROC
    /** \brief Set an input by index
    *
    *  @copydoc Tvalset
    *  @copydoc iind
    */

    void setInput(T val, int iind=0);

     /** \brief Set an output by index
     *
     * @copydoc Tvalset
     * @copydoc oind
     */
    void setOutput(T val, int oind=0);

    /** \brief Set an input by name
    *
    *  @copydoc Tvalset
    *  @copydoc iname
    *
     */
    void setInput(T val, const std::string &iname);

    /** \brief Set an output by name
    *
    *  @copydoc Tvalset
    *  @copydoc oname
    *
     */
    void setOutput(T val, const std::string &oname);
#endif

    /// @}

/// \cond INTERNAL
#ifdef DOXYGENPROC
    /// \name Advanced Getters
    ///
    /// @{

    /** \brief Get an input by index
    *
    *  @copydoc Tvalget
    *  @copydoc iind
    */

    void getInput(T val, int iind=0);

     /** \brief Get an output by index
     *
     * @copydoc Tvalget
     * @copydoc oind
     */
    void getOutput(T val, int oind=0);

    /** \brief Get an input by name
    *
    *  @copydoc Tvalget
    *  @copydoc iname
    *
     */
    void getInput(T val, const std::string &iname);

    /** \brief Get an output by name
    *
    *  @copydoc Tvalget
    *  @copydoc oname
    *
     */
    void getOutput(T val, const std::string &oname);
    /// @}
#endif
/// \endcond


#define SETTERS_SUB_DUMMY(T) \
    void setInput(T val, int iind=0);                                   \
    void setOutput(T val, int oind=0);                                  \
    void setInput(T val, const std::string &iname);                     \
    void setOutput(T val, const std::string &oname);

#define SETTERS_NZ(T) \
    void setInputNZ(T val, int iind=0) {                                \
      self().assertInit();                                              \
      try {                                                             \
        self().input(iind).setNZ(val);                                  \
      } catch(std::exception& e) {                                      \
        casadi_error(e.what() << "Occurred at iind = " << iind << "."); \
      }                                                                 \
    }                                                                   \
    void setOutputNZ(T val, int oind=0) {                               \
      self().assertInit(); self().output(oind).setNZ(val);              \
    }                                                                   \
    void setInputNZ(T val, const std::string &iname) {                  \
      setInputNZ(val, self().inputIndex(iname));                        \
    }                                                                   \
    void setOutputNZ(T val, const std::string &oname) {                 \
      setOutputNZ(val, self().outputIndex(oname));                      \
    }

#define SETTERS_SUB(T)                                                  \
    void setInput(T val, int iind=0) {                                  \
      self().assertInit();                                              \
      try {                                                             \
        self().input(iind).set(val);                                    \
      } catch(std::exception& e) {                                      \
        casadi_error(e.what() << "Occurred at iind = " << iind << "."); \
      }                                                                 \
    }                                                                   \
    void setOutput(T val, int oind=0) {                                 \
      self().assertInit(); self().output(oind).set(val);                \
    }                                                                   \
    void setInput(T val, const std::string &iname) {                    \
      setInput(val, self().inputIndex(iname));                          \
    }                                                                   \
    void setOutput(T val, const std::string &oname) {                   \
      setOutput(val, self().outputIndex(oname));                        \
    }

#ifndef DOXYGENPROC
#ifndef SWIG
    SETTERS_SUB_DUMMY(const double*)
    SETTERS_SUB_DUMMY(const std::vector<double>&)
    SETTERS_NZ(const double*) // NOLINT(readability/casting) - false positive
#endif // SWIG
    SETTERS_SUB(const Matrix<double>&)
    SETTERS_NZ(const Matrix<double>&)
#endif // DOXYGENPROC

#undef SETTERS_NZ
#undef SETTERS_SUB
#undef SETTERS_SUB_DUMMY

#define GETTERS_NZ(T)                                           \
      void getInputNZ(T val, int iind=0) const {                \
        self().assertInit(); self().input(iind).getNZ(val);     \
      }                                                         \
      void getOutputNZ(T val, int oind=0) const {               \
        self().assertInit(); self().output(oind).getNZ(val);    \
      }                                                         \
      void getInputNZ(T val, const std::string &iname) const {  \
        getInputNZ(val, self().inputIndex(iname));              \
      }                                                         \
      void getOutputNZ(T val, const std::string &oname) const { \
        getOutputNZ(val, self().outputIndex(oname));            \
      }

#define GETTERS_SUB(T)                                                  \
      void getInput(T val, int iind=0) const {                          \
        self().assertInit(); self().input(iind).get(val);               \
      }                                                                 \
      void getOutput(T val, int oind=0) const {                         \
        self().assertInit(); self().output(oind).get(val);              \
      }                                                                 \
      void getInput(T val, const std::string &iname) const {            \
        getInput(val, self().inputIndex(iname));                        \
      }                                                                 \
      void getOutput(T val, const std::string &oname) const {           \
        getOutput(val, self().outputIndex(oname));                      \
      }

#ifndef DOXYGENPROC
#ifndef SWIG
GETTERS_SUB(double&)
#ifndef SWIG
GETTERS_NZ(double*) // NOLINT(readability/casting) - false positive
GETTERS_NZ(std::vector<double>&)
#endif // SWIG
GETTERS_SUB(Matrix<double>&)
#endif // SWIG
#endif // DOXYGENPROC
#undef GETTERS_NZ
#undef GETTERS_SUB

  };


} // namespace casadi


#endif // CASADI_IO_INTERFACE_HPP
