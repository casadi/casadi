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

#ifndef COMPILED_FUNCTION_INTERNAL_HPP
#define COMPILED_FUNCTION_INTERNAL_HPP

#include "compiled_function.hpp"
#include "fx_internal.hpp"

#ifdef WITH_DL 
#ifdef _WIN32 // also for 64-bit
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL 

namespace CasADi{
  
class CompiledFunctionInternal : public FXInternal{
  friend class CompiledFunction;
  public:

    /** \brief  constructor */
    explicit CompiledFunctionInternal(const std::string& bin_name);

    /** \brief  clone function */
    virtual CompiledFunctionInternal* clone() const;

    /** \brief  Destructor */
    virtual ~CompiledFunctionInternal();

    /** \brief  Evaluate */
    virtual void evaluate(int nfdir, int nadir);
  
    /** \brief  Initialize */
    virtual void init();

  protected:

//@{
/** \brief  Function pointer types */
  typedef int (*evaluatePtr)(const double** x, double** r);
  typedef int (*initPtr)(int *n_in_, int *n_out_);
  typedef int (*getSparsityPtr)(int n_in, int *n_row, int *n_col, int **rowind, int **col);
//@}

  /** \brief  Name of binary */
  std::string bin_name_;

  /** \brief  Function pointers */
  evaluatePtr evaluate_;
    
#if defined(WITH_DL) && defined(_WIN32) // also for 64-bit
  typedef HINSTANCE handle_t;
#else
  typedef void* handle_t;
#endif

  /** \brief  handle to the dll */
  handle_t handle_;
  
  /** \brief  Array of pointers to the input */
  std::vector<const double*> input_array_;
  
  /** \brief  Array of pointers to the output */
  std::vector<double*> output_array_;
  
};

} // namespace CasADi


#endif // COMPILED_FUNCTION_INTERNAL_HPP
