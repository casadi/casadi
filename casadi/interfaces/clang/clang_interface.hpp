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


#ifndef CASADI_CLANG_INTERFACE_HPP
#define CASADI_CLANG_INTERFACE_HPP

#include "casadi/core/function/jit_function_internal.hpp"
#include <casadi/interfaces/clang/casadi_jitfunction_clang_export.h>

/** \defgroup plugin_JitFunction_clang
      Interface to the JIT compiler CLANG
*/

/** \pluginsection{JitFunction,clang} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{JitFunction,clang}


   \author Joris Gillis
   \date 2015
   *
   @copydoc JitFunction_doc
   @copydoc plugin_JitFunction_clang
   * */
  class CASADI_JITFUNCTION_CLANG_EXPORT ClangJitFunctionInterface : public JitFunctionInternal {
  public:

    /** \brief Constructor */
    explicit ClangJitFunctionInterface(const Function &f);

    /** \brief Clone */
    virtual ClangJitFunctionInterface* clone() const;

    /** \brief  Create a new JIT function */
    static JitFunctionInternal* creator(const Function &f)
    { return new ClangJitFunctionInterface(f);}

    /** \brief  Evaluate numerically, work vectors given */
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief Destructor */
    virtual ~ClangJitFunctionInterface();

    /** \brief Initialize */
    virtual void init();

    /// A documentation string
    static const std::string meta_doc;

  protected:
    typedef int (*sparsityPtr)(int i, int *n_row, int *n_col,
                               const int **colind, const int **row);
    typedef int (*workPtr)(int *n_iw, int *n_w);
    typedef int (*evalPtr)(const double** arg, double** res, int* iw, double* w);

    /** \brief Function pointer for initialization of external */
    typedef int (*initPtr)(int *f_type, int *n_in, int *n_out, int *sz_arg, int* sz_res);

    sparsityPtr sparsity_fun_;
    evalPtr eval_fun_;
    workPtr work_fun_;
    initPtr init_fun_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_CLANG_INTERFACE_HPP
