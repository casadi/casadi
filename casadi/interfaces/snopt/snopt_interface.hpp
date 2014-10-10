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


#ifndef CASADI_SNOPT_INTERFACE_HPP
#define CASADI_SNOPT_INTERFACE_HPP

#include "casadi/core/function/nlp_solver_internal.hpp"
#include "casadi/interfaces/snopt/casadi_nlpsolver_snopt_export.h"
#include "casadi/interfaces/snopt/snopt.h"
#include "snoptProblem.hpp"

/** \defgroup plugin_NlpSolver_snopt
  SNOPT interface
*/

/** \pluginsection{NlpSolver,snopt} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{NlpSolver,snopt}
     @copydoc NlpSolver_doc
     @copydoc plugin_NlpSolver_snopt
  */
  class CASADI_NLPSOLVER_SNOPT_EXPORT SnoptInterface : public NlpSolverInternal {

  public:
    // Constructor
    explicit SnoptInterface(const Function& nlp);

    // Destructor
    virtual ~SnoptInterface();

    // Clone function
    virtual SnoptInterface* clone() const;

    /** \brief  Create a new NLP Solver */
    static NlpSolverInternal* creator(const Function& nlp)
    { return new SnoptInterface(nlp);}

    // (Re)initialize
    virtual void init();

    // Solve the NLP
    virtual void evaluate();


    virtual void setQPOptions();

    /// Read options from snopt parameter xml
    virtual void setOptionsFromFile(const std::string & file);

    /// Exact Hessian?
    bool exact_hessian_;

    std::map<int, std::string> status_;
    std::map<std::string, opt_type> ops_;

    std::string formatStatus(int status) const;

    void userfun(int* mode, int nnObj, int nnCon, int nnJac, int nnL, int neJac,
                 double* x, double* fObj, double*gObj, double* fCon, double* gCon,
                 int nState, char* cu, int lencu, int* iu, int leniu, double* ru, int lenru);

    void callback(int* iAbort, int* info, int HQNType, int* KTcond, int MjrPrt, int minimz,
        int m, int maxS, int n, int nb, int nnCon0, int nnCon, int nnObj0, int nnObj, int nS,
        int itn, int nMajor, int nMinor, int nSwap,
        double condHz, int iObj, double sclObj, double ObjAdd,
        double fMrt,  double PenNrm, double step,
        double prInf,  double duInf,  double vimax,  double virel, int* hs,
        int ne, int nlocJ, int* locJ, int* indJ, double* Jcol, int negCon,
        double* Ascale, double* bl, double* bu, double* fCon, double* gCon, double* gObj,
        double* yCon, double* pi, double* rc, double* rg, double* x,
        char* cu, int lencu, int* iu, int leniu, double* ru, int lenru,
        char* cw, int lencw, int* iw, int leniw, double* rw, int lenrw);

    int nnJac_;
    int nnObj_;
    int nnCon_;

    /// Classification arrays
    /// original variable index -> category w.r.t. f
    std::vector<int> x_type_g_;
    /// original variable index -> category w.r.t. g
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

    static void userfunPtr(int * mode, int* nnObj, int * nnCon, int *nJac, int *nnL, int * neJac,
                           double *x, double *fObj, double *gObj, double * fCon, double* gCon,
                           int* nState, char* cu, int* lencu, int* iu, int* leniu,
                           double* ru, int *lenru);
    static void snStopPtr(
        int* iAbort, int* info, int* HQNType, int* KTcond, int* MjrPrt, int* minimz,
        int* m, int* maxS, int* n, int* nb,
        int* nnCon0, int* nnCon, int* nnObj0, int* nnObj, int* nS,
        int* itn, int* nMajor, int* nMinor, int* nSwap,
        double * condHz, int* iObj, double * sclObj,  double *ObjAdd,
        double * fMrt,  double * PenNrm,  double * step,
        double *prInf,  double *duInf,  double *vimax,  double *virel, int* hs,
        int* ne, int* nlocJ, int* locJ, int* indJ, double* Jcol, int* negCon,
        double* Ascale, double* bl, double* bu, double* fCon, double* gCon, double* gObj,
        double* yCon, double* pi, double* rc, double* rg, double* x,
        char* cu, int* lencu, int* iu, int* leniu, double* ru, int *lenru,
        char* cw, int* lencw, int* iw, int *leniw, double* rw, int* lenrw);

    // Matrix A has a linear objective row
    bool jacF_row_;
    // Matrix A has a dummy row
    bool dummyrow_;

    /// A documentation string
    static const std::string meta_doc;

  private:
      // options
      std::map<std::string, std::string> intOpts_;
      std::map<std::string, std::string> realOpts_;
      std::map<std::string, std::pair<std::string, std::string> > strOpts_;

      /// Pass the supplied options to Snopt
      void passOptions(snoptProblemC &probC);

      // Accumulated time in last evaluate():
      double t_eval_grad_f_; // time spent in eval_grad_f
      double t_eval_jac_g_; // time spent in eval_jac_g
      double t_callback_fun_;  // time spent in callback function
      double t_mainloop_; // time spent in the main loop of the solver

      // Accumulated counts in last evaluate():
      int n_eval_grad_f_; // number of calls to eval_grad_f
      int n_eval_jac_g_; // number of calls to eval_jac_g
      int n_callback_fun_; // number of calls to callback function
      int n_iter_; // number of major iterations

  };

} // namespace casadi

/// \endcond
#endif // CASADI_SNOPT_INTERFACE_HPP
