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

#ifndef IO_INTERFACE_HPP
#define IO_INTERFACE_HPP

#include "../sx/sx_element.hpp"
#include "../mx/mx.hpp"
#include "../options_functionality.hpp"
#include "../fx/schemes_helpers.hpp"
#include "../fx/io_scheme.hpp"

namespace CasADi{
  
  /** \brief Interface for accessing input and output data structures
      \author Joel Andersson
      \date 2013
  */
  template<class Derived>
  class IOInterface {
  public:
    
    //@{
    /// Access input argument
    inline const Matrix<double>& input(int iind=0) const{ return inputS<true>(iind);}
    inline const Matrix<double>& input(const std::string &iname) const{  return input(inputSchemeEntry(iname)); }
#ifdef SWIG
    %rename(inputRef) input;
#endif
    inline Matrix<double>& input(int iind=0){ return inputS<true>(iind);}
    inline Matrix<double>& input(const std::string &iname){ return input(inputSchemeEntry(iname));}
    //@}
    
    //@{
    /// Access output argument
    inline const Matrix<double>& output(int oind=0) const{ return outputS<true>(oind);}
    inline const Matrix<double>& output(const std::string &oname) const{ return output(outputSchemeEntry(oname));}
#ifdef SWIG
    %rename(outputRef) output;
#endif
    inline Matrix<double>& output(int oind=0){ return outputS<true>(oind);}
    inline Matrix<double>& output(const std::string &oname){ return output(outputSchemeEntry(oname));}
    //@}

    /// Get the number of function inputs
    inline int getNumInputs() const{ return static_cast<const Derived*>(this)->input_struct().data.size();}

    /// Get the number of function outputs
    inline int getNumOutputs() const{ return static_cast<const Derived*>(this)->output_struct().data.size();}

    /// Set the number of function inputs
    inline void setNumInputs(int num_in){ static_cast<Derived*>(this)->input_struct().data.resize(num_in); }
    
    /// Set the number of function outputs
    inline void setNumOutputs(int num_out){ static_cast<Derived*>(this)->output_struct().data.resize(num_out); }

    /** \brief  Access an input */
    template<bool check>
    DMatrix& inputS(int i){
      if(check){
        try{
          return static_cast<Derived*>(this)->input_struct().data.at(i);
        } catch(std::out_of_range&){
          std::stringstream ss;
          ss <<  "In function " << static_cast<const Derived*>(this)->getOption("name") << ": input " << i << " not in interval [0," << getNumInputs() << ")";
          if (!static_cast<const Derived*>(this)->isInit()) ss << std::endl << "Did you forget to initialize?";
          throw CasadiException(ss.str());
        }
      } else {
        return static_cast<Derived*>(this)->input_struct().data[i];        
      }
    }

    /** \brief  Const access an input */
    template<bool check>
    inline const DMatrix& inputS(int i) const{
      return const_cast<IOInterface<Derived>*>(this)->inputS<check>(i);
    }
    
    /** \brief  Access an output*/
    template<bool check>
    DMatrix& outputS(int i){
      if(check){
        try{
          return static_cast<Derived*>(this)->output_struct().data.at(i);
        } catch(std::out_of_range&){
          std::stringstream ss;
          ss <<  "In function " << static_cast<const Derived*>(this)->getOption("name") << ": output " << i << " not in interval [0," << getNumOutputs() << ")";
          if (!static_cast<const Derived*>(this)->isInit()) ss << std::endl << "Did you forget to initialize?";
          throw CasadiException(ss.str());
        }
      } else {
        return static_cast<Derived*>(this)->output_struct().data[i];        
      }
    }

    /** \brief  Const access an output*/
    template<bool check>
    inline const DMatrix& outputS(int i) const{
      return const_cast<IOInterface<Derived>*>(this)->outputS<check>(i);
    }

    /** \brief Set input scheme */
    void setInputScheme(const CasADi::IOScheme &scheme){
      casadi_assert(scheme.compatibleSize(getNumInputs()));
      static_cast<Derived*>(this)->inputScheme() = scheme;
    }

    /** \brief Set output scheme */
    void setOutputScheme(const CasADi::IOScheme &scheme){ 
      casadi_assert(scheme.compatibleSize(getNumOutputs()));
      static_cast<Derived*>(this)->outputScheme() = scheme;
    }

    /** \brief Get input scheme */
    CasADi::IOScheme getInputScheme() const{ return static_cast<const Derived*>(this)->inputScheme(); }

    /** \brief Get output scheme */
    CasADi::IOScheme getOutputScheme() const{ return static_cast<const Derived*>(this)->outputScheme(); }

    /** \brief Find the index for a string describing a particular entry of an input scheme
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FXInternal adheres to SCHEME_NLPINput 
     */
    int inputSchemeEntry(const std::string &name) const{ return schemeEntry(static_cast<const Derived*>(this)->inputScheme(),name,true);}

    /** \brief Find the index for a string describing a particular entry of an output scheme
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FXInternal adheres to SCHEME_NLPINput 
     */
    int outputSchemeEntry(const std::string &name) const{ return schemeEntry(static_cast<const Derived*>(this)->outputScheme(),name,false);}

    /** \brief Find the index for a string describing a particular entry of a scheme
     * example:  schemeEntry("x_opt")  -> returns  NLP_SOLVER_X if FXInternal adheres to SCHEME_NLPINput 
     */
    int schemeEntry(const CasADi::IOScheme &scheme,const std::string &name,bool input) const{
      if (scheme.isNull()) casadi_error("Unable to look up '" <<  name<< "' in " << (input? "input": "output") << "scheme, as the " <<  (input? "input": "output") << " scheme of this function is unknown. You can only index with integers.");
      if (name=="") casadi_error("FXInternal::schemeEntry: you supplied an empty string as the name of a entry in " << scheme.name() << ". Available names are: " << scheme.entryNames() << ".");
      int n = scheme.index(name);
      if (n==-1) casadi_error("FXInternal::schemeEntry: could not find entry '" << name << "' in " << scheme.name() << ". Available names are: " << scheme.entryNames() << ".");
      return n;
    }

#ifdef DOXYGENPROC
    /// \name Setters
    /// Set/get an input, output, forward seed/sensitivity or adjoint seed/sensitivity\n
    /// T can be double&, double*, std::vector<double>&, Matrix<double> &\n
    /// Assumes a properly allocated val.\n
    /// 
    /// @{
    /** 
        \brief Reads in the input argument from val.
        T can be double&, double*, std::vector<double>&, Matrix<double> &\n
        Assumes a properly allocated val.\n
    */
    void setInput(T val, int iind=0) const;
    /** 
        \brief Reads in the output argument from val.
        T can be double&, double*, std::vector<double>&, Matrix<double> &\n
        Assumes a properly allocated val.\n
    */
    void setOutput(T val, int oind=0) const;
    /// @}

#endif

#define SETTERS(T)                                                        \
    void setInput(T val, int iind=0)             { static_cast<const Derived*>(this)->assertInit(); input(iind).set(val);  } \
    void setOutput(T val, int oind=0)            { static_cast<const Derived*>(this)->assertInit(); output(oind).set(val); } \
    void setInput(T val, const std::string &iname)             { setInput(val,inputSchemeEntry(iname));  } \
    void setOutput(T val, const std::string &oname)            { setOutput(val,outputSchemeEntry(oname)); } \

#ifndef DOXYGENPROC
    SETTERS(double);
#ifndef SWIG
    SETTERS(const double*);
#endif // SWIG
    SETTERS(const std::vector<double>&);
    SETTERS(const Matrix<double>&);
#endif // DOXYGENPROC

#undef SETTERS
    
#define GETTERS(T)                                                        \
    void getInput(T val, int iind=0) const             { static_cast<const Derived*>(this)->assertInit(); input(iind).get(val);} \
    void getOutput(T val, int oind=0) const            { static_cast<const Derived*>(this)->assertInit(); output(oind).get(val);} \
    void getInput(T val, const std::string &iname) const             { getInput(val,inputSchemeEntry(iname)); } \
    void getOutput(T val, const std::string &oname) const            { getOutput(val,outputSchemeEntry(oname)); } \
    
#ifndef DOXYGENPROC
#ifndef SWIG
GETTERS(double&);
GETTERS(double*);
GETTERS(std::vector<double>&);
GETTERS(Matrix<double>&);
#endif // SWIG
#endif // DOXYGENPROC
#undef GETTERS

#ifdef DOXYGENPROC
    /// \name Getters
    /// A group of accessor for numerical data that operate on preallocated data.\n
    /// get an input, output, forward seed/sensitivity or adjoint seed/sensitivity\n
    //    T can be double&, double*, std::vector<double>&, Matrix<double> &\n
    //    Assumes a properly allocated val.\n
    /// @{

    /** \brief Writes out the input argument into val. \noswig
        T can be double&, double*, std::vector<double>&, Matrix<double> &\n
        Assumes a properly allocated val.\n
    */
    void getInput(T val, int iind=0) const;
 
    /** 
        \brief Writes out the output argument into val. \noswig
        T can be double&, double*, std::vector<double>&, Matrix<double> &\n
        Assumes a properly allocated val.\n
    */
    void getOutput(T val, int oind=0) const;
    /// @}
    
 
  ///@{
  /// Get the input as a new Matrix \nocpp
  Matrix<double> getInput(int iind=0) const;
  /// Get the input as a new Matrix \nocpp
  Matrix<double> getInput(const std::string &iname) const ;
  /// Get the output as a new Matrix \nocpp
  Matrix<double> getOutput(int oind=0) const ;
  /// Get the output as a new Matrix \nocpp
  Matrix<double> getOutput(const std::string &oname) const ;
  ///@}
#endif

  };


} // namespace CasADi


#endif // IO_INTERFACE_HPP
