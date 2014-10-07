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
  class CASADI_CORE_EXPORT IOInterface {
  public:

    /// \cond UNSAFE
    /** \brief [UNSAFE] Obtain reference to inputs
     * \seealso getInput, setInput
     */
    ///@{
    /// Access input argument
    inline const Matrix<double>& input(int iind=0) const { return inputS<true>(iind);}
    inline const Matrix<double>& input(const std::string &iname) const
    { return input(inputSchemeEntry(iname)); }
#ifdef SWIG
    %rename(inputRef) input;
#endif
    inline Matrix<double>& input(int iind=0) { return inputS<true>(iind);}
    inline Matrix<double>& input(const std::string &iname) { return input(inputSchemeEntry(iname));}
    ///@}

    /** \brief [UNSAFE] Obtain reference to outputs
     * \seealso getOutput, getOutput
     */
    ///@{
    /// Access output argument
    inline const Matrix<double>& output(int oind=0) const { return outputS<true>(oind);}
    inline const Matrix<double>& output(const std::string &oname) const
    { return output(outputSchemeEntry(oname));}
#ifdef SWIG
    %rename(outputRef) output;
#endif
    inline Matrix<double>& output(int oind=0) { return outputS<true>(oind);}
    inline Matrix<double>& output(const std::string &oname)
    { return output(outputSchemeEntry(oname));}
    ///@}
    /// \endcond

    /// Get the number of function inputs
    inline int getNumInputs() const
    { return static_cast<const Derived*>(this)->input_struct().data.size();}

    /// Get the number of function outputs
    inline int getNumOutputs() const
    { return static_cast<const Derived*>(this)->output_struct().data.size();}

    /// Set the number of function inputs
    inline void setNumInputs(int num_in)
    { static_cast<Derived*>(this)->input_struct().data.resize(num_in); }

    /// Set the number of function outputs
    inline void setNumOutputs(int num_out)
    { static_cast<Derived*>(this)->output_struct().data.resize(num_out); }

    /// \cond INTERNAL

    /** \brief  Access an input */
    template<bool check>
    DMatrix& inputS(int i) {
      if (check) {
        try {
          return static_cast<Derived*>(this)->input_struct().data.at(i);
        } catch(std::out_of_range&) {
          std::stringstream ss;
          ss <<  "In function " << static_cast<const Derived*>(this)->getOption("name")
             << ": input " << i << " not in interval [0, " << getNumInputs() << ")";
          if (!static_cast<const Derived*>(this)->isInit()) ss << std::endl
                                                               << "Did you forget to initialize?";
          throw CasadiException(ss.str());
        }
      } else {
        return static_cast<Derived*>(this)->input_struct().data[i];
      }
    }

    /** \brief  Const access an input */
    template<bool check>
    inline const DMatrix& inputS(int i) const {
      return const_cast<IOInterface<Derived>*>(this)->inputS<check>(i);
    }

    /** \brief  Access an output*/
    template<bool check>
    DMatrix& outputS(int i) {
      if (check) {
        try {
          return static_cast<Derived*>(this)->output_struct().data.at(i);
        } catch(std::out_of_range&) {
          std::stringstream ss;
          ss <<  "In function " << static_cast<const Derived*>(this)->getOption("name")
             << ": output " << i << " not in interval [0, " << getNumOutputs() << ")";
          if (!static_cast<const Derived*>(this)->isInit()) ss << std::endl
                                                               << "Did you forget to initialize?";
          throw CasadiException(ss.str());
        }
      } else {
        return static_cast<Derived*>(this)->output_struct().data[i];
      }
    }

    /** \brief  Const access an output*/
    template<bool check>
    inline const DMatrix& outputS(int i) const {
      return const_cast<IOInterface<Derived>*>(this)->outputS<check>(i);
    }

    /// \endcond

    /** \brief Set input scheme */
    void setInputScheme(const casadi::IOScheme &scheme) {
      casadi_assert(scheme.compatibleSize(getNumInputs()));
      static_cast<Derived*>(this)->inputScheme() = scheme;
    }

    /** \brief Set output scheme */
    void setOutputScheme(const casadi::IOScheme &scheme) {
      casadi_assert(scheme.compatibleSize(getNumOutputs()));
      static_cast<Derived*>(this)->outputScheme() = scheme;
    }

    /** \brief Get input scheme */
    casadi::IOScheme getInputScheme() const
    { return static_cast<const Derived*>(this)->inputScheme(); }

    /** \brief Get output scheme */
    casadi::IOScheme getOutputScheme() const
    { return static_cast<const Derived*>(this)->outputScheme(); }

    /// \cond INTERNAL
    /** \brief Find the index for a string describing a particular entry of an input scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FunctionInternal adheres to
     * SCHEME_NLPINput
     */
    int inputSchemeEntry(const std::string &name) const
    { return schemeEntry(static_cast<const Derived*>(this)->inputScheme(), name, true);}

    /** \brief Find the index for a string describing a particular entry of an output scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FunctionInternal adheres to
     * SCHEME_NLPINput
     */
    int outputSchemeEntry(const std::string &name) const
    { return schemeEntry(static_cast<const Derived*>(this)->outputScheme(), name, false);}

    /** \brief Find the index for a string describing a particular entry of a scheme
     *
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FunctionInternal adheres to
     * SCHEME_NLPINput
     */
    int schemeEntry(const casadi::IOScheme &scheme, const std::string &name, bool input) const {
      if (scheme.isNull()) casadi_error("Unable to look up '" <<  name<< "' in "
                                        << (input? "input": "output") << "scheme, as the "
                                        <<  (input? "input": "output")
                                        << " scheme of this function is unknown. "
                                        "You can only index with integers.");
      if (name=="") casadi_error("FunctionInternal::schemeEntry: you supplied an empty "
                                 "string as the name of a entry in "
                                 << scheme.name() << ". Available names are: "
                                 << scheme.entryNames() << ".");
      int n = scheme.index(name);
      if (n==-1) casadi_error("FunctionInternal::schemeEntry: could not find entry '"
                              << name << "' in " << scheme.name() << ". Available names are: "
                              << scheme.entryNames() << ".");
      return n;
    }
    /// \endcond

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
    *   \param[in] iind index within the range [0..getNumInputs()-1]
    */

    /**
    * \defgroup  oind
    *   \param[in] oind index within the range [0..getNumOutputs()-1]
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
    Matrix<double> getInput(int iind=0) const { return input(iind); }
    /** \brief Get an input by name
    *
    *  @copydoc iname
    *
     */
    Matrix<double> getInput(const std::string &iname) const  { return input(iname); }

    /** \brief Get an output by index
    *
    *  @copydoc oind
    */
    Matrix<double> getOutput(int oind=0) const  { return output(oind); }
    /** \brief Get an output by name
    *
    *  @copydoc oname
    *
     */
    Matrix<double> getOutput(const std::string &oname) const  { return output(oname); }

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


#define SETTERS(T)                                                              \
    void setInput(T val, int iind=0)                                            \
    { static_cast<const Derived*>(this)->assertInit();                          \
      try { input(iind).set(val); }                                             \
      catch(std::exception& e) {                                                \
        casadi_error(e.what() << "Occurred at iind = " << iind << ".");         \
      }                                                                         \
    }                                                                           \
    void setOutput(T val, int oind=0)                                           \
    { static_cast<const Derived*>(this)->assertInit(); output(oind).set(val); } \
    void setInput(T val, const std::string &iname)                              \
    { setInput(val, inputSchemeEntry(iname));  }                                 \
    void setOutput(T val, const std::string &oname)                             \
    { setOutput(val, outputSchemeEntry(oname)); }                                \

#ifndef DOXYGENPROC
    SETTERS(double) // NOLINT(readability/casting) - false positive
#ifndef SWIG
    SETTERS(const double*)
#endif // SWIG
    SETTERS(const std::vector<double>&)
    SETTERS(const Matrix<double>&)
#endif // DOXYGENPROC

#undef SETTERS

#define GETTERS(T)                                                             \
    void getInput(T val, int iind=0) const                                     \
    { static_cast<const Derived*>(this)->assertInit(); input(iind).get(val);}  \
    void getOutput(T val, int oind=0) const                                    \
    { static_cast<const Derived*>(this)->assertInit(); output(oind).get(val);} \
    void getInput(T val, const std::string &iname) const                       \
    { getInput(val, inputSchemeEntry(iname)); }                                 \
    void getOutput(T val, const std::string &oname) const                      \
    { getOutput(val, outputSchemeEntry(oname)); }                               \

#ifndef DOXYGENPROC
#ifndef SWIG
GETTERS(double&)
GETTERS(double*) // NOLINT(readability/casting) - false positive
GETTERS(std::vector<double>&)
GETTERS(Matrix<double>&)
#endif // SWIG
#endif // DOXYGENPROC
#undef GETTERS

  };


} // namespace casadi


#endif // CASADI_IO_INTERFACE_HPP
