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


#ifndef IP_METHOD_HPP
#define IP_METHOD_HPP

#include "core/function/nlp_solver.hpp"

namespace casadi{
  
class IPInternal;
  
/**
  \brief Interior point method
  This method is experimental only. Do not attempt to use if you do not intend to dive into the source code.
  The current purpose of the class is to show how an IP method can be implemeted in CasADi.
  If someone wants to take responsibility for this class and make it work, then please contact the CasADi developers.
  
  \author Joel Andersson
  \date 2012
*/
class IPMethod : public NlpSolver {
  public:
    /// Default constructor
    IPMethod();

    /// \brief Constuct an NLP with non-linear constraints and provided hessian approximation
    explicit IPMethod(const Function& F,         /**< F objective function: \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}]\f$*/
                      const Function& G          /**< constraint function (default only bound constraints): \f$ [\mathbf{R}^n] \mapsto [\mathbf{R}^m]\f$ */
                      );

    /// Access functions of the node
    IPInternal* operator->();
    const IPInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function 
    #ifdef SWIG
    %callback("%s_cb");
    #endif
    static NlpSolver creator(const Function& F, const Function& G, int dummy){ return IPMethod(F,G);}
    #ifdef SWIG
    %nocallback;
    #endif
};

} // namespace casadi

#endif //IP_METHOD_HPP
