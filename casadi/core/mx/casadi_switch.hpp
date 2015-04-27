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


#ifndef CASADI_SWITCH_HPP
#define CASADI_SWITCH_HPP

#include "casadi_call.hpp"

/// \cond INTERNAL

namespace casadi {

  /** Embeds a function call in an expression graph
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT Switch : public GenericCall {
  public:
    /** \brief  Create function call node */
    static std::vector<MX> create(const MX& ind, const std::vector<MX>& arg,
                                  const std::vector<Function>& f,
                                  const Function& f_def);

    /** \brief  Destructor */
    virtual ~Switch() {}

    /** \brief  Clone function */
    virtual Switch* clone() const;

    /** \brief  Print argument vector */
    static std::string printArg(const std::vector<std::string>& arg);

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Generate code for the operation */
    virtual void generate(const std::vector<int>& arg, const std::vector<int>& res,
                          CodeGenerator& g) const;

    /// Evaluate the function numerically
    virtual void evalD(cp_double* arg, p_double* res, int* itmp, double* rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(cp_SXElement* arg, p_SXElement* res, int* itmp, SXElement* rtmp);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(cp_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(p_bvec_t* arg, p_bvec_t* res, int* itmp, bvec_t* rtmp);

    /** \brief  Number of functions */
    virtual int numFunctions() const {return f_.size() + 1;}

    /** \brief  Get function reference */
    virtual Function& getFunction(int i);

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
    virtual int getOp() const { return OP_SWITCH;}

    /// Get number of temporary variables needed
    virtual void nTmp(size_t& ni, size_t& nr);

  private:
    /** \brief  Constructor */
    explicit Switch(const MX& ind, const std::vector<MX>& arg,
                    const std::vector<Function>& f, const Function& f_def);

    // Function to be evaluated
    std::vector<Function> f_;
    Function f_def_;

    // Input and output sparsities
    std::vector<Sparsity> sp_in_, sp_out_;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SWITCH_HPP
