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

#ifndef FX_TOOLS_HPP
#define FX_TOOLS_HPP

#include "fx.hpp"
#include "parallelizer.hpp"
#include "c_function.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

namespace CasADi{
  
class MultipleShooting{
  public:
    // Constructor
    MultipleShooting(const FX& fcn, const FX& mfcn, int ns, int nx, int nu);
    
    // Initialize
    void init();
    
#ifndef SWIG
    // Jacobian callback function
    static void jacobian_wrapper(CFunction &f, int fsens_order, int asens_order, void* user_data);
    
    // Jacobian of the NLP
    void jacobian(CFunction &f, int fsens_order, int asens_order);
#endif // SWIG

    // Discrete time dynamics
    FX fcn_;
    
    // Mayer term
    FX mfcn_;
    
    //Numboer of shooting nodes
    int ns_;
    
    // Number of controls
    int nu_;
    
    // Number of differential states
    int nx_;

    // Variable bound and initial guess
    std::vector<double> V_min_, V_max_, V_init_;
    
    // Constraint bounds
    std::vector<double> G_min_, G_max_;

    // NLP objective function
    MXFunction F_;
    
    // NLP constraint function
    MXFunction G_;

    // Jacobian of the NLP constraints
    CFunction J_;

    // Parallel evaluation of the Jacobian blocks
    Parallelizer JX_,JP_;
    
    // Mapping from the Jacobian blocks to the sparse Jacobian
    SXFunction J_mapping_;
    
    // Control bounds and initial guess
    std::vector<double> u_min_, u_max_, u_init_;
    
    // State bounds and initial guess
    std::vector<double> x_min_, x_max_, x_init_;

    //State bounds at the initial time
    std::vector<double> x0_min_, x0_max_;

    //State bounds at the final time
    std::vector<double> xf_min_, xf_max_;

    //Final time
    double tf_;
};
                        
                        
} // namespace CasADi


#endif // FX_TOOLS_HPP
