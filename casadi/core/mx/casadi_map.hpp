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

    /** \brief  Constructor */
    explicit Map(const Function& fcn, const std::vector<std::vector<MX> >& arg);

    /** \brief  Destructor */
    virtual ~Map() {}

    /** \brief  Clone function */
    virtual Map* clone() const;

    /** \brief  Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /// Evaluate the function numerically
    virtual void evalD(cp_double* arg, p_double* res, int* itmp, double* rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(cp_SXElement* arg, p_SXElement* res,
                        int* itmp, SXElement* rtmp);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(cp_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(p_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Get function reference */
    virtual Function& getFunction();

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

    /// Get number of temporary variables needed
    virtual void nTmp(size_t& ni, size_t& nr);

    /// Get parallelization
    virtual std::string parallelization() const { return "serial";}

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
    virtual void evalD(cp_double* arg, p_double* res, int* itmp, double* rtmp);

    /// Get number of temporary variables needed
    virtual void nTmp(size_t& ni, size_t& nr);

    /// Get parallelization
    virtual std::string parallelization() const { return "openmp";}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_MAP_HPP
