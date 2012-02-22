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

#ifndef X_FUNCTION_HPP
#define X_FUNCTION_HPP

#include "fx.hpp"
#include "../sx/sx.hpp"
#include <vector>

namespace CasADi{

/// Forward declaration of internal class
class XFunctionInternal;

/**   \brief Dynamically created function that can be expanded into a series of scalar operations. Base class for XFunction and MXFunction.
\author Joel Andersson 
\date 2011
*/

class XFunction : public FX{

public:
  /// Default constructor
  XFunction();

  /// Access functions of the node 
  XFunctionInternal* operator->();

  /// Const access functions of the node 
  const XFunctionInternal* operator->() const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// evaluate symbolically, SX type (overloaded)
  std::vector<SXMatrix> eval(const std::vector<SXMatrix>& arg){ return evalSX(arg);}

  /// evaluate symbolically, MX type (overloaded)
  MXVector eval(const MXVector& arg){return evalMX(arg);}
  
  /// evaluate symbolically, MX type (unambiguous)
  MXVector evalMX(const MXVector& arg);

  /// evaluate symbolically, SX type (unambiguous)
  std::vector<SXMatrix> evalSX(const std::vector<SXMatrix>& arg);
  
  /** \brief Evaluate symbolically with with directional derivatives, SX type
   * The first two arguments are the nondifferentiated inputs and results of the evaluation,
   * the next two arguments are a set of forward directional seeds and the resulting forward directional derivatives,
   * the length of the vector being the number of forward directions.
   * The next two arguments are a set of adjoint directional seeds and the resulting adjoint directional derivatives,
   * the length of the vector being the number of adjoint directions.
   * The first boolean argument allows the second argument to the functions to be used as an input instead of output,
   * assuming it is already known and the second boolean arguments allows constants to be eliminated during the 
   * evaluations (as the treatment of constants in CasADi will get more efficient, this will become unnecessary).
   */
  void evalSX(const SXMatrixVector& input, SXMatrixVector& output, 
	      const SXMatrixVectorVector& fwdSeed, SXMatrixVectorVector& fwdSens, 
	      const SXMatrixVectorVector& adjSeed, SXMatrixVectorVector& adjSens,
	      bool output_given=false, bool eliminate_constants=false);

  /** \brief Evaluate symbolically with with directional derivatives, MX type
   * The first two arguments are the nondifferentiated inputs and results of the evaluation,
   * the next two arguments are a set of forward directional seeds and the resulting forward directional derivatives,
   * the length of the vector being the number of forward directions.
   * The next two arguments are a set of adjoint directional seeds and the resulting adjoint directional derivatives,
   * the length of the vector being the number of adjoint directions.
   * The first boolean argument allows the second argument to the functions to be used as an input instead of output,
   * assuming it is already known and the second boolean arguments allows constants to be eliminated during the 
   * evaluations (as the treatment of constants in CasADi will get more efficient, this will become unnecessary).
   */
  void evalMX(const MXVector& input, MXVector& output, 
	      const MXVectorVector& fwdSeed, MXVectorVector& fwdSens, 
	      const MXVectorVector& adjSeed, MXVectorVector& adjSens,
	      bool output_given=false, bool eliminate_constants=false);
                        
  /** \brief Evaluate symbolically with with directional derivatives, SX type, overloaded
   * The first two arguments are the nondifferentiated inputs and results of the evaluation,
   * the next two arguments are a set of forward directional seeds and the resulting forward directional derivatives,
   * the length of the vector being the number of forward directions.
   * The next two arguments are a set of adjoint directional seeds and the resulting adjoint directional derivatives,
   * the length of the vector being the number of adjoint directions.
   * The first boolean argument allows the second argument to the functions to be used as an input instead of output,
   * assuming it is already known and the second boolean arguments allows constants to be eliminated during the 
   * evaluations (as the treatment of constants in CasADi will get more efficient, this will become unnecessary).
   */
  void eval(const SXMatrixVector& input, std::vector<SXMatrix>& output, 
	    const SXMatrixVectorVector& fwdSeed, SXMatrixVectorVector& fwdSens, 
	    const SXMatrixVectorVector& adjSeed, SXMatrixVectorVector& adjSens,
	    bool output_given=false, bool eliminate_constants=false);

  /** \brief Evaluate symbolically with with directional derivatives, MX type, overloaded
   * The first two arguments are the nondifferentiated inputs and results of the evaluation,
   * the next two arguments are a set of forward directional seeds and the resulting forward directional derivatives,
   * the length of the vector being the number of forward directions.
   * The next two arguments are a set of adjoint directional seeds and the resulting adjoint directional derivatives,
   * the length of the vector being the number of adjoint directions.
   * The first boolean argument allows the second argument to the functions to be used as an input instead of output,
   * assuming it is already known and the second boolean arguments allows constants to be eliminated during the 
   * evaluations (as the treatment of constants in CasADi will get more efficient, this will become unnecessary).
   */
  void eval(const MXVector& input, MXVector& output, 
	    const MXVectorVector& fwdSeed, MXVectorVector& fwdSens, 
	    const MXVectorVector& adjSeed, MXVectorVector& adjSens,
	    bool output_given=false, bool eliminate_constants=false);
  
#ifndef SWIG
  /// evaluate symbolically (pass and get non-zero entries) LEGACY - REMOVE
  std::vector< std::vector<SX> > eval(const std::vector< std::vector<SX> >& arg);

  /// evaluate symbolically, single input, single output 
  SXMatrix eval(const SXMatrix& arg);

  /// evaluate symbolically, single input, single output (pass and get non-zero entries) LEGACY - REMOVE
  std::vector<SX> eval(const std::vector<SX>& arg);
#endif // SWIG
  
};

} // namespace CasADi

#endif // X_FUNCTION_HPP
