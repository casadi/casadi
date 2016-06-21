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


#ifndef CASADI_MA27_INTERFACE_HPP
#define CASADI_MA27_INTERFACE_HPP

#include "casadi/core/function/linsol_internal.hpp"
#include <casadi/interfaces/hsl/casadi_linsol_ma27_export.h>

extern "C" {
  void ma27id_(int* ICNTL, double* CNTL);
  void ma27ad_(int *N, int *NZ, const int *IRN, const int* ICN,
               int *IW, int* LIW, int* IKEEP, int *IW1,
               int* NSTEPS, int* IFLAG, int* ICNTL,
               double* CNTL, int *INFO, double* OPS);
  void ma27bd_(int *N, int *NZ, const int *IRN, const int* ICN,
               double* A, int* LA, int* IW, int* LIW,
               int* IKEEP, int* NSTEPS, int* MAXFRT,
               int* IW1, int* ICNTL, double* CNTL,
               int* INFO);
  void ma27cd_(int *N, double* A, int* LA, int* IW,
               int* LIW, double* W, int* MAXFRT,
               double* RHS, int* IW1, int* NSTEPS,
               int* ICNTL, double* CNTL);
}

/** \defgroup plugin_Linsol_ma27
 * Interface to the sparse direct linear solver MA27
 * Works for symmetric indefinite systems
 * Partly adopted from qpOASES 3.2
 * \author Joel Andersson
 * \date 2016
 */

/** \pluginsection{Linsol,ma27} */
/// \cond INTERNAL
namespace casadi {
  struct CASADI_LINSOL_MA27_EXPORT Ma27Memory : public LinsolMemory {
    // Constructor
    Ma27Memory();

    // Destructor
    ~Ma27Memory();

    /* Work vector for MA27AD */
    std::vector<int> iw1;

    /* Number of nonzeros in the current linear system. */
    int nnz;

    /* matrix/factor for MA27 (A in MA27) */
    std::vector<double> nz;

    /* Row entries of matrix (IRN in MA27) */
    std::vector<int> irn;

    /* Column entries of matrix (JCN in MA27) */
    std::vector<int> jcn;

    /* integer control values (ICNRL in MA27) */
    int icntl[30];

    /* real control values (CNRL in MA27) */
    double cntl[5];

    /* integer work space (IW in MA27) */
    std::vector<int> iw;

    /* Real work space (W in MA27CD) */
    std::vector<double> w;

    /* IKEEP in MA27 */
    std::vector<int> ikeep;

    /* NSTEPS in MA27 */
    int nsteps;

    /* MAXFRT in MA27 */
    int maxfrt;

    /* number of negative eigenvalues */
    int neig;

    // Rank of matrix
    int rank;
  };

  /** \brief \pluginbrief{Linsol,ma27}
   * @copydoc Linsol_doc
   * @copydoc plugin_Linsol_ma27
   */
  class CASADI_LINSOL_MA27_EXPORT Ma27Interface : public LinsolInternal {
  public:

    // Create a linear solver given a sparsity pattern and a number of right hand sides
    Ma27Interface(const std::string& name);

    /** \brief  Create a new Linsol */
    static LinsolInternal* creator(const std::string& name) {
      return new Ma27Interface(name);
    }

    // Destructor
    virtual ~Ma27Interface();

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new Ma27Memory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<Ma27Memory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    // Set sparsity pattern
    virtual void reset(void* mem, const int* sp) const;

    // Factorize the linear system
    virtual void factorize(void* mem, const double* A) const;

    /// Number of negative eigenvalues
    virtual int neig(void* mem) const;

    /// Matrix rank
    virtual int rank(void* mem) const;

    // Solve the linear system
    virtual void solve(void* mem, double* x, int nrhs, bool tr) const;

    /// A documentation string
    static const std::string meta_doc;

    // Get name of the plugin
    virtual const char* plugin_name() const { return "ma27";}

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;
  };

} // namespace casadi

/// \endcond

#endif // CASADI_MA27_INTERFACE_HPP
