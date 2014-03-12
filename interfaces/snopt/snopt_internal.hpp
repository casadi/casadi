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

#ifndef SNOPT_INTERNAL_HPP
#define SNOPT_INTERNAL_HPP

#include "symbolic/fx/nlp_solver_internal.hpp"

namespace CasADi{

  /**
     @copydoc NLPSolver_doc
  */
  class SnoptInternal : public NLPSolverInternal{

  public:
    // Constructor
    explicit SnoptInternal(const FX& nlp);

    // Destructor
    virtual ~SnoptInternal();

    // Clone function
    virtual SnoptInternal* clone() const;

    // Reset solver
    void reset();

    // (Re)initialize
    virtual void init();
    
    // Solve the NLP
    virtual void evaluate();

    
    virtual void setQPOptions();

    /// Read options from worhp parameter xml
    void setOptionsFromFile(const std::string & file);
  
    /// Exact Hessian?
    bool exact_hessian_;
  
    std::map<int,std::string> status_;
    std::map<std::string,opt_type> ops_;

    // Accummulated time since last reset:
    double t_eval_f_; // time spent in eval_f
    double t_eval_grad_f_; // time spent in eval_grad_f
    double t_eval_g_; // time spent in eval_g
    double t_eval_jac_g_; // time spent in eval_jac_g
    double t_eval_h_; // time spent in eval_h
    double t_callback_fun_;  // time spent in callback function
    double t_callback_prepare_; // time spent in callback preparation
    double t_mainloop_; // time spent in the main loop of the solver
  
    std::string formatStatus(int status) const;
  
    /// Pass the supplied options to Snopt
    void passOptions();
    
    /// Work arrays for SNOPT
    std::vector<char> snopt_cw_;
    std::vector<int> snopt_iw_;
    std::vector<double> snopt_rw_;
    
    void userfun(int & mode, int nnObj, int nnCon, int nnJac, int nnL, int neJac, double* x, double &fObj, double*gObj, double* fCon, double* gCon, int nState, char* cu, int lencu, int* iu, int leniu, double* ru, int lenru);
  
    int nnJac_;
    int nnObj_;
    int nnCon_;
    
    /// Classification arrays
    /// original variable index -> category w.r.t f
    std::vector<int> x_type_g_;
    /// original variable index -> category w.r.t g
    std::vector<int> x_type_f_;
    /// original constraint index -> category
    std::vector<int> g_type_;
    
    /// sorted variable index -> original variable index
    std::vector<int> x_order_;
    /// sorted constraint index -> original constraint index
    std::vector<int> g_order_;
    
    IMatrix A_structure_;
    std::vector<double> A_data_;
    
    
    std::vector<double> bl_;
    std::vector<double> bu_;
    std::vector<int> hs_;
    std::vector<double> x_;
    std::vector<double> pi_;
    std::vector<double> rc_;
    
    // Do detection of linear substructure
    bool detect_linear_;
    
    int m_;
    int iObj_;
    
    static void userfunPtr(int * mode, int* nnObj, int * nnCon, int *nJac, int *nnL, int * neJac, double *x, double *fObj, double *gObj, double * fCon, double* gCon, int* nState, char* cu, int* lencu, int* iu, int* leniu, double* ru, int *lenru);
    
    typedef std::map< std::string, std::pair< opt_type, std::string> > OptionsMap;
    
    OptionsMap optionsmap_;
    
    // Matrix A has a inear objective row
    bool jacF_row_;
    // Matrix A has a dummy row
    bool dummyrow_;
    
  private:
      void snInit(int iPrint, int iSumm);
      void snSeti(const std::string &snopt_name, int value);
      void snSetr(const std::string &snopt_name, double value);
      void snSet(const std::string &snopt_name, const std::string &value);
  };
  


} // namespace CasADi

#endif //SNOPT_INTERNAL_HPP
