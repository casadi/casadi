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

#ifndef ACADO_INTEGRATOR_BACKEND_HPP
#define ACADO_INTEGRATOR_BACKEND_HPP

#include <acado/integrator/integrator.hpp>
#include "acado_internal.hpp"

namespace CasADi{

class AcadoIntegratorBackend : public ACADO::Integrator{
  public:
    /// Constructor
    explicit AcadoIntegratorBackend(void *user_data);
    
    /// Static create function
    static Integrator* create(void *user_data);
    
    /// Destructor
    virtual ~AcadoIntegratorBackend( );
  
  protected:
    
    // The calling acado interface
    CasADi::AcadoInternal* ocp_solver_;
    
    // Current derivative
    int ider_;

    /// Initialize the solver given a differential equation
    virtual ACADO::returnValue init( const ACADO::DifferentialEquation &rhs);

    // Deep copy of the integrator
    CasADi::Integrator integrator_;
    
public:
    
        
    /** Default constructor. */
    AcadoIntegratorBackend( const ACADO::DifferentialEquation &rhs_ );

    /** Copy constructor (deep copy). */
    AcadoIntegratorBackend( const AcadoIntegratorBackend& arg );

    /** Destructor. */

    /** Assignment operator (deep copy). */
    virtual AcadoIntegratorBackend& operator=( const AcadoIntegratorBackend& arg );

    /** The (virtual) copy constructor */
    virtual Integrator* clone() const;
    
    virtual ACADO::returnValue freezeMesh();
    virtual ACADO::returnValue freezeAll();
    virtual ACADO::returnValue unfreeze();
    virtual ACADO::returnValue step( int number  /**< the step number */ );
    virtual ACADO::returnValue stop();
    virtual ACADO::returnValue setDxInitialization( double *dx0 );
    virtual int getNumberOfSteps() const;
    virtual int getNumberOfRejectedSteps() const;
    virtual double getStepSize() const;

  protected:
    virtual ACADO::returnValue evaluate( const ACADO::Vector &x0, const ACADO::Vector &xa, const ACADO::Vector &p, const ACADO::Vector &u, const ACADO::Vector &w,  const ACADO::Grid   &t_);
    virtual ACADO::returnValue evaluateSensitivities();
    virtual ACADO::returnValue setProtectedForwardSeed( const ACADO::Vector &xSeed, const ACADO::Vector &pSeed, const ACADO::Vector &uSeed, const ACADO::Vector &wSeed, const int &order);
    virtual ACADO::returnValue setProtectedBackwardSeed( const ACADO::Vector &seed, const int &order );
    virtual ACADO::returnValue getProtectedX( ACADO::Vector *xEnd ) const;
    virtual ACADO::returnValue getProtectedForwardSensitivities( ACADO::Matrix *Dx, int order) const;
    virtual ACADO::returnValue getProtectedBackwardSensitivities( ACADO::Vector &Dx_x0, ACADO::Vector &Dx_p, ACADO::Vector &Dx_u, ACADO::Vector &Dx_w, int order) const;
    virtual int getDim() const;
    virtual int getDimX() const;
    virtual ACADO::returnValue setBackwardSeed2( const ACADO::Vector &seed);

};


} // namespace CasADi



// #include <acado/integrator/integrator_bdf.ipp>


#endif  // ACADO_INTEGRATOR_BACKEND_HPP

// end of file.
