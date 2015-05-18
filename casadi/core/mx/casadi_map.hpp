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


#ifndef CASADI_MAP_HPP
#define CASADI_MAP_HPP

#include "multiple_output.hpp"
#include "../function/function.hpp"

/// \cond INTERNAL

namespace casadi {

  /**
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT Map : public MultipleOutput {
  public:
    /** \brief  Create map node */
    static std::vector<std::vector<MX> >
      create(const Function& fcn, const std::vector<std::vector<MX> > &arg,
             const std::string& parallelization);

    /** \brief  Destructor */
    virtual ~Map() {}

    /** \brief  Clone function */
    virtual Map* clone() const;

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(const SXElement** arg, SXElement** res,
                        int* iw, SXElement* w);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Get function reference */
    virtual const Function& getFunction() const;

    /** \brief  Get function input */
    virtual int getFunctionInput() const { return -1;}

    /** \brief  Get function output */
    virtual int getFunctionOutput() const { return -1;}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Number of outputs */
    virtual int nout() const;

    /** \brief  Get the sparsity of output oind */
    virtual const Sparsity& sparsity(int oind) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_MAP;}

    /** \brief Get number of temporary variables needed */
    virtual void nwork(size_t& n_arg, size_t& n_res, size_t& n_iw, size_t& n_w) const;

    /// Get parallelization
    virtual std::string parallelization() const { return "serial";}

  protected:
    /** \brief  Constructor */
    explicit Map(const Function& fcn, const std::vector<std::vector<MX> >& arg);

    // Function to be evaluated
    Function fcn_;

    /// Number of times called
    int n_;
  };

  /**
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT OmpMap : public Map {
  public:

    /** \brief  Constructor */
    explicit OmpMap(const Function& fcn, const std::vector<std::vector<MX> >& arg);

    /** \brief  Destructor */
    virtual ~OmpMap() {}

    /** \brief  Clone function */
    virtual OmpMap* clone() const;

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /** \brief Get number of temporary variables needed */
    virtual void nwork(size_t& n_arg, size_t& n_res, size_t& n_iw, size_t& n_w) const;

    /// Get parallelization
    virtual std::string parallelization() const { return "openmp";}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MAP_HPP
