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


#ifndef CASADI_SLICOT_EXPM_HPP
#define CASADI_SLICOT_EXPM_HPP

#include "../../core/expm_impl.hpp"
#include "../../core/linsol.hpp"
#include <casadi/interfaces/slicot/casadi_expm_slicot_export.h>

/** \defgroup plugin_Expm_slicot
 *
*/

/** \pluginsection{Expm,slicot} */

/// \cond INTERNAL
namespace casadi {


  // Forward declaration
  class SlicotExpm;

  struct CASADI_EXPM_SLICOT_EXPORT SlicotExpmMemory {

    double *A, *H;
    double *dwork;
    int* iwork;

    /// Constructor
    SlicotExpmMemory() {}

    /// Destructor
    ~SlicotExpmMemory() {}
  };

  /** \brief \pluginbrief{Expm,slicot}
   *
   * An efficient solver for Discrete Periodic Lyapunov Equations using SLICOT
   *
   * @copydoc Expm_doc
   * @copydoc plugin_Expm_slicot

       \author Joris Gillis
      \date 2014

  */
  class CASADI_EXPM_SLICOT_EXPORT SlicotExpm : public Expm {
  public:
    /** \brief  Constructor */
    explicit SlicotExpm();

    /** \brief  Constructor
     * \param st \structargument{Expm}
     */
    SlicotExpm(const std::string& name, const Sparsity& A);

    /** \brief  Create a new QP Solver */
    static Expm* creator(const std::string& name,
                          const Sparsity& A) {
      return new SlicotExpm(name, A);
    }

    /** \brief  Destructor */
    ~SlicotExpm() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "slicot";}

    /** \brief  Initialize */
    void init(const Dict& opts) override;

    /** \brief Create memory block */
    void* alloc_memory() const override { return new SlicotExpmMemory();}

    /** \brief Free memory block */
    void free_memory(void *mem) const override { delete static_cast<SlicotExpmMemory*>(mem);}

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const override;

    /** \brief Initalize memory block */
    void init_memory(void* mem) const override;

    /** \brief  Evaluate numerically */
    void eval(const double** arg, double** res, int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;


  private:

    int n_;


  };

} // namespace casadi

/// \endcond
#endif // CASADI_SLICOT_EXPM_HPP
