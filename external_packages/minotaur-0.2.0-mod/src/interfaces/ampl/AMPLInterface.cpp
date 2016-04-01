//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file AMPLInterface.cpp
 * \brief Define the AMPLInterface class for reading problems from AMPL.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <sstream>
#include <stdint.h>
#include <iostream>
#include "opcode.hd"

#include "MinotaurConfig.h"
#include "AMPLInterface.h"
#include "AMPLNonlinearFunction.h"
#include "CGraph.h"
#include "CNode.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "Logger.h"
#include "LinearFunction.h"
#include "Option.h"
#include "PolynomialFunction.h"
#include "Problem.h"
#include "QuadraticFunction.h"
#include "Solution.h"
#include "SOS.h"
#include "Variable.h"


// 
// This is REALLY important. when passing error flags to AMPL functions,
// initialize their value (preferably to 0), even if AMPL is supposed to fill
// in the value and return it. 
//

using namespace MINOTAUR_AMPL;

const std::string AMPLInterface::me_ = "AMPLInterface: ";


// suf_sos_ASL implemented in ASL's suf_sos.c needs priority, ref, sos, sosno
// and sosref.
static SufDecl suftab[] =
{
  { const_cast<char*>("priority"), 0, ASL_Sufkind_var, 1 },
  { const_cast<char*>("ref"), 0, ASL_Sufkind_var | ASL_Sufkind_real, 1 },
  { const_cast<char*>("sos"), 0, ASL_Sufkind_var, 1 },
  { const_cast<char*>("sos"), 0, ASL_Sufkind_con, 1 },
  { const_cast<char*>("sosno"), 0, ASL_Sufkind_var | ASL_Sufkind_real, 1 },
  { const_cast<char*>("sosref"), 0, ASL_Sufkind_var | ASL_Sufkind_real, 1 }
};




// Constructor 
AMPLInterface::AMPLInterface(Minotaur::EnvPtr env, std::string solver) 
 : env_(env),
   intTol_(1e-8),
   myAsl_(NULL),
   nCons_(0),
   nDefVars_(0),
   nDefVarsBco_(0),
   nDefVarsCo1_(0),
   nVars_(0),
   zTol_(1e-8)
{
#ifndef USE_MINOTAUR_AMPL_INTERFACE
#error need to define USE_MINOTAUR_AMPL_INTERFACE
#endif
  std::string str="minotaurampl";
  logger_ = (Minotaur::LoggerPtr) new Minotaur::Logger((Minotaur::LogLevel)
      env_->getOptions()->findInt("ampl_log_level")->getValue());
  getOptionsFromEnv_(solver);
  addOptions_();
}


AMPLInterface::~AMPLInterface()
{
  freeASL();
  functionMap_.clear();
  vars_.clear();
}


void AMPLInterface::addDefinedVars_(Minotaur::ProblemPtr instance)
{
  std::string name;
  Minotaur::VariablePtr v;
  for (int i=0; i<nDefVars_; ++i) {
    std::stringstream name_stream;
    name_stream << "defvar" << i;
    name = name_stream.str();
    v = instance->newVariable(-INFINITY, INFINITY, Minotaur::Continuous, name);
    vars_.push_back(v);
  }
}


void AMPLInterface::addLinearConstraint_(int i, 
                                         Minotaur::ProblemPtr instance)
{
  std::string cName;
  Minotaur::LinearFunctionPtr lfPtr = Minotaur::LinearFunctionPtr(); // NULL
  Minotaur::FunctionPtr fPtr;
  
  assert (i <= myAsl_->i.n_con_);
  cName = std::string(con_name_ASL(myAsl_, i));
  addLinearTermsFromConstr_(lfPtr, i);
  fPtr = (Minotaur::FunctionPtr) new Minotaur::Function(lfPtr);

  instance->newConstraint(fPtr, myAsl_->i.LUrhs_[2*i], 
      myAsl_->i.LUrhs_[2*i+1], cName);
}


void AMPLInterface::addLinearObjective_(int i, Minotaur::ProblemPtr instance)
{
  Minotaur::FunctionPtr fPtr;
  Minotaur::LinearFunctionPtr lfPtr = Minotaur::LinearFunctionPtr(); // null
  Minotaur::ObjectiveType obj_sense;
  std::string oName;

  assert (i <= myAsl_->i.n_obj_);

  if (myAsl_->i.objtype_[i] == 1) {
    obj_sense = Minotaur::Maximize;
  } else {
    obj_sense = Minotaur::Minimize;
  }

  addLinearTermsFromObj_(lfPtr, i);

  oName = std::string(obj_name_ASL(myAsl_, i));
  fPtr = (Minotaur::FunctionPtr) new Minotaur::Function(lfPtr);
  instance->newObjective(fPtr, objconst_ASL(myAsl_, 0), 
      obj_sense, oName);
}


// if lfPtr is null, allocate memory and simply add terms from ASL. ASL can be
// trusted for no repititions. Otherwise, there are some existing terms in lf
// and we need to use lf->incTerm().
void AMPLInterface::addLinearTermsFromConstr_(Minotaur::LinearFunctionPtr & lf,
                                              int i)
{
  cgrad *cg; // for ampl
  if (lf) {
    for (cg = myAsl_->i.Cgrad_[i]; cg; cg = cg->next) {
      lf->incTerm(vars_[cg->varno], cg->coef);
    }
  } else {
    lf= (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    for (cg = myAsl_->i.Cgrad_[i]; cg; cg = cg->next) {
      lf->addTerm(vars_[cg->varno], cg->coef);
    }
  }
  if (lf->getNumTerms()==0) {
    lf.reset();
  }
}


// same as getLinearTermsFromConstr_() but for objective.
void AMPLInterface::addLinearTermsFromObj_(Minotaur::LinearFunctionPtr & lf, 
                                           int i)
{
  ograd *og; // for ampl
  if (lf) {
    for (og = myAsl_->i.Ograd_[i]; og; og = og->next) {
      lf->incTerm(vars_[og->varno], og->coef);
    }
  } else {
    lf= (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    for (og = myAsl_->i.Ograd_[i]; og; og = og->next) {
      lf->addTerm(vars_[og->varno], og->coef);
    }
  }
  if (lf->getNumTerms()==0) {
    lf.reset();
  }
}


void AMPLInterface::addOptions_()
{
  Minotaur::FlagOptionPtr f_option; 
  Minotaur::BoolOptionPtr b_option; 
  Minotaur::OptionDBPtr options = env_->getOptions();
  f_option = (Minotaur::FlagOptionPtr) new Minotaur::Option<bool>
    ("AMPL", "If given, then write .sol file for ampl.", true, false);
  options->insert(f_option, true);

  f_option = (Minotaur::FlagOptionPtr) new Minotaur::Option<bool>
    ("v", "If given, then write version information.", true, false);
  options->insert(f_option, true);

  f_option = (Minotaur::FlagOptionPtr) new Minotaur::Option<bool>
    ("=", "If given, then write all known options.", true, false);
  options->insert(f_option, true);

  f_option = (Minotaur::FlagOptionPtr) new Minotaur::Option<bool>
    ("?", "If given, then write help message.", true, false);
  options->insert(f_option, true);

  b_option = (Minotaur::BoolOptionPtr) new Minotaur::Option<bool>
    ("display_ampl_model", 
     "If true, write ampl model before creating the problem: <0/1>", 
     true, false);
  options->insert(b_option);

}


void AMPLInterface::addQuadraticConstraint_(int i, 
                                            Minotaur::ProblemPtr instance)
{
  {
    std::string cName;
    Minotaur::FunctionPtr fPtr = Minotaur::FunctionPtr();  //NULL
    Minotaur::LinearFunctionPtr lfPtr = Minotaur::LinearFunctionPtr();  //NULL
    Minotaur::QuadraticFunctionPtr qfPtr = Minotaur::QuadraticFunctionPtr();
    Minotaur::PolyFunPtr pfPtr = Minotaur::PolyFunPtr();  //NULL
    double c = 0;
    cde *constraint_cde;
    // cast myAsl_ from ASL into ASL_fg
    ASL_fg *asl_fg = (ASL_fg *)myAsl_;

    constraint_cde = asl_fg->I.con_de_+i;
    assert (constraint_cde);
    getPoly_(lfPtr, qfPtr, pfPtr, c, constraint_cde->e);
    addLinearTermsFromConstr_(lfPtr, i);
    if (logger_->getMaxLevel() >= Minotaur::LogDebug2) {
      if (lfPtr) {
        lfPtr->write(logger_->msgStream(Minotaur::LogDebug2));
      } 
      if (qfPtr) {
        qfPtr->write(logger_->msgStream(Minotaur::LogDebug2));
      }
      if (pfPtr) {
        pfPtr->write(logger_->msgStream(Minotaur::LogDebug2));
      }
      logger_->msgStream(Minotaur::LogDebug2) << c << std::endl;
    }

    // we are ready to add the constraint now.
    fPtr = (Minotaur::FunctionPtr) new Minotaur::Function(lfPtr, qfPtr, pfPtr);
    cName = std::string(con_name_ASL(myAsl_, i));
    instance->newConstraint(fPtr, myAsl_->i.LUrhs_[2*i] - c, 
                            myAsl_->i.LUrhs_[2*i+1] - c, cName); 
  }
}


void AMPLInterface::addQuadraticDefCons_(Minotaur::ProblemPtr instance)
{
  double c;
  Minotaur::QuadraticFunctionPtr qf;
  Minotaur::LinearFunctionPtr lf;
  Minotaur::FunctionPtr f;
  Minotaur::PolyFunPtr pf;
  cexp  *cstruct; // cexp is defined in nlp.h of asl
  cexp1 *cstruct1; // cexp1 is defined in nlp.h of asl
  int v_index;
  linpart *L;    // linpart is defined in asl.h

  // visit each 'defined variable' and add the constraint that is used to
  // define it. The defintion could be linear
  // or nonlinear.
  for (int i=0; i<nDefVarsBco_; ++i) {
    c = 0;
    qf = Minotaur::QuadraticFunctionPtr();  // NULL
    lf = Minotaur::LinearFunctionPtr();     // NULL
    cstruct = (((const ASL_fg *)myAsl_)->I.cexps_)+i;
    getPoly_(lf, qf, pf, c, cstruct->e);
    if (!lf) {
      lf = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    }
    L = cstruct->L;
    for (int j=0; j<cstruct->nlin; ++j) {
      v_index = ((uintptr_t)L->v.rp - (uintptr_t)((const ASL_fg *) myAsl_)->
          I.var_e_)/sizeof(expr_v);
      lf->incTerm(vars_[v_index], L->fac);
      ++L;
    }
    // the constraint is: lf + qf + c - defvar_[i] = 0
    lf->incTerm(vars_[nVars_+i],-1);
    f = (Minotaur::FunctionPtr) new Minotaur::Function(lf, qf, pf);
    instance->newConstraint(f, -c, -c);
  }
  for (int i=0; i<nDefVarsCo1_; ++i) {
    c = 0;
    qf = Minotaur::QuadraticFunctionPtr(); 
    lf = Minotaur::LinearFunctionPtr(); // NULL
    cstruct1 = (((const ASL_fg *)myAsl_)->I.cexps1_)+i;
    getPoly_(lf, qf, pf, c, cstruct1->e);
    if (!lf) {
      lf = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    }
    L = cstruct1->L;
    for (int j=0; j<cstruct1->nlin; ++j) {
      v_index = ((uintptr_t)L->v.rp - (uintptr_t)((const ASL_fg *) myAsl_)->
          I.var_e_)/sizeof(expr_v);
      lf->incTerm(vars_[v_index], L->fac);
      ++L;
    }
    // the constraint is: lf + qf + c - defvar_[i] = 0
    lf->incTerm(vars_[nVars_+nDefVarsBco_+i],-1);
    f = (Minotaur::FunctionPtr) new Minotaur::Function(lf, qf, pf);
    instance->newConstraint(f, -c, -c);
  }
}


void AMPLInterface::addQuadraticDefCons2_(Minotaur::ProblemPtr instance)
{
  double c;
  Minotaur::QuadraticFunctionPtr qf = Minotaur::QuadraticFunctionPtr();
  Minotaur::LinearFunctionPtr lf;
  Minotaur::FunctionPtr f;
  Minotaur::CGraphPtr cgraph;
  Minotaur::CNode *cnode;
  cexp  *cstruct; // cexp is defined in nlp.h of asl
  cexp1 *cstruct1; // cexp1 is defined in nlp.h of asl
  int v_index;
  linpart *L;    // linpart is defined in asl.h
  int err = 0;

  // visit each 'defined variable' and add the constraint that is used to
  // define it. The defintion could be linear
  // or nonlinear.
  for (int i=0; i<nDefVarsBco_; ++i) {
    c = 0;
    cstruct = (((const ASL_fg *)myAsl_)->I.cexps_)+i;
    lf = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    L = cstruct->L;
    for (int j=0; j<cstruct->nlin; ++j) {
      v_index = ((uintptr_t)L->v.rp - (uintptr_t)((const ASL_fg *) myAsl_)->
          I.var_e_)/sizeof(expr_v);
      lf->incTerm(vars_[v_index], L->fac);
      ++L;
    }
    // the constraint is: lf + nlf + c - defvar_[i] = 0
    lf->incTerm(vars_[nVars_+i],-1);

    cgraph = (Minotaur::CGraphPtr) new Minotaur::CGraph();
    cnode = getCGraph_(cstruct->e, cgraph, instance);
    cgraph->setOut(cnode);
    cgraph->finalize();
    if (Minotaur::Constant==cgraph->getType()) {
      c = cgraph->eval(0, &err); assert(0==err);
      cgraph.reset();
    } 
    f = (Minotaur::FunctionPtr) new Minotaur::Function(lf, qf, cgraph);
    instance->newConstraint(f, -c , -c);
  }


  for (int i=0; i<nDefVarsCo1_; ++i) {
    c = 0;
    cstruct1 = (((const ASL_fg *)myAsl_)->I.cexps1_)+i;
    L = cstruct1->L;
    lf = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    for (int j=0; j<cstruct1->nlin; ++j) {
      v_index = ((uintptr_t)L->v.rp - (uintptr_t)((const ASL_fg *) myAsl_)->
          I.var_e_)/sizeof(expr_v);
      lf->incTerm(vars_[v_index], L->fac);
      ++L;
    }
    // the constraint is: lf + nlf + c - defvar_[i] = 0
    lf->incTerm(vars_[nVars_+nDefVarsBco_+i],-1);

    cgraph = (Minotaur::CGraphPtr) new Minotaur::CGraph();
    cnode = getCGraph_(cstruct1->e, cgraph, instance);
    cgraph->setOut(cnode);
    cgraph->finalize();
    if (Minotaur::Constant==cgraph->getType()) {
      c = cgraph->eval(0, &err); assert(0==err);
      cgraph.reset();
    }
    f = (Minotaur::FunctionPtr) new Minotaur::Function(lf, qf, cgraph);
    instance->newConstraint(f, -c, -c);
  }
}


void AMPLInterface::addQuadraticObjective_(int i, 
                                           Minotaur::ProblemPtr instance)
{
  cde *obj_cde;
  // cast myAsl_ from ASL into ASL_fg
  ASL_fg *asl_fg = (ASL_fg *)myAsl_;

  if (myAsl_->i.n_obj_ > 0) {
    Minotaur::ObjectiveType obj_sense;
    if (myAsl_->i.objtype_[i] == 1) {
      obj_sense = Minotaur::Maximize;
    } else {
      obj_sense = Minotaur::Minimize;
    }

    if (myAsl_->i.nlo_ > 0) {
      std::string oName;
      Minotaur::LinearFunctionPtr lfPtr = Minotaur::LinearFunctionPtr();
      Minotaur::QuadraticFunctionPtr qfPtr = Minotaur::QuadraticFunctionPtr();
      Minotaur::PolyFunPtr pfPtr = Minotaur::PolyFunPtr();
      Minotaur::FunctionPtr fPtr;
      double c = 0;

      obj_cde = asl_fg->I.obj_de_+i;
      assert (obj_cde);
      getPoly_(lfPtr, qfPtr, pfPtr, c, obj_cde->e);
      addLinearTermsFromObj_(lfPtr, i);
      oName = std::string(obj_name_ASL(myAsl_, 0));
      fPtr = (Minotaur::FunctionPtr) new Minotaur::Function(lfPtr, qfPtr, 
          pfPtr);
      instance->newObjective(fPtr, objconst_ASL(myAsl_, 0) + c, 
          obj_sense, oName);
    } 
  }
}


void AMPLInterface::addVariablesFromASL_(Minotaur::ProblemPtr instance)
{
  Minotaur::UInt stop_index, start_index;
  Minotaur::VariablePtr vPtr;
  std::string vName;

  // add variables to the instance
  
  //
  // Ordering of variables by ampl:
  // http://www.ampl.com/REFS/HOOKING/index.html
  //
  // category           count
  // -------------------------------------------------------------------------
  // nonlinear          max(nlvc,nlvo); (see next table)
  // linear arcs        nwv;
  // other linear       n_var - (max{nlvc,nlvo} + niv + nbv + nwv);
  // binary             nbv;
  // other integers     niv;
  // -------------------------------------------------------------------------
  //
  // Ordering nonlinear variables
  // Nonlinear-category where appears   count
  // -------------------------------------------------------------------------
  // continuous         obj and constr  nlvb - nlvbi;
  // integer            obj and constr  nlvbi;
  // continuous         only constr     nlvc - (nlvb + nlvci);
  // integer            only constr     nlvci;
  // continuous         only obj        nlvo - (nlvc + nlvoi); **
  // integer            only obj        nlvoi;
  //
  // ** may seem to be a bug. you may think that nlvc should be replaced by
  // nlvo. it is NOT a bug.  
  // if nlvo > nlvc, the first nlvc variables in the objective may or may not
  // be nonlinear in objective, but the next nlvo-nlvc are definitely
  // nonlinear in objective. if nlvo <= nlvc all of the first nlvo variables
  // apper nonlinearly in the objective.
  // -------------------------------------------------------------------------
  //

  // first visit all nonlinear variables
  // continuous nonlinear variables in both obj and cons
  start_index = 0;
  stop_index = myAsl_->i.nlvb_ - myAsl_->i.nlvbi_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    vName = std::string(var_name_ASL(myAsl_, i));
    vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
        myAsl_->i.LUv_[2*i+1], Minotaur::Continuous, vName);
    vars_.push_back(vPtr);
  }
  // integer nonlinear variables in both obj and cons
  start_index = stop_index;
  stop_index += myAsl_->i.nlvbi_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    vName = std::string(var_name_ASL(myAsl_, i));
    vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
        myAsl_->i.LUv_[2*i+1], Minotaur::Integer, vName);
    vars_.push_back(vPtr);
  }
  // continuous nonlinear variables in cons only
  start_index = stop_index;
  stop_index += myAsl_->i.nlvc_ - (myAsl_->i.nlvb_ + myAsl_->i.nlvci_);
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    vName = std::string(var_name_ASL(myAsl_, i));
    vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
        myAsl_->i.LUv_[2*i+1], Minotaur::Continuous, vName);
    vars_.push_back(vPtr);
  }
  // integer nonlinear variables in cons only
  start_index = stop_index;
  stop_index += myAsl_->i.nlvci_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    vName = std::string(var_name_ASL(myAsl_, i));
    vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
        myAsl_->i.LUv_[2*i+1], Minotaur::Integer, vName);
    vars_.push_back(vPtr);
  }

  if (myAsl_->i.nlvo_ > myAsl_->i.nlvc_) {
    // there are some variables that are linear in constraints but nonlinear
    // in objective. the first 'nlvc' are already counted.
    // continuous nonlinear variables in obj only
    start_index = stop_index;
    stop_index += myAsl_->i.nlvo_ - (myAsl_->i.nlvc_ + myAsl_->i.nlvoi_);
    for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
      vName = std::string(var_name_ASL(myAsl_, i));
      vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
          myAsl_->i.LUv_[2*i+1], Minotaur::Continuous, vName);
      vars_.push_back(vPtr);
    }
    // integer nonlinear variables in obj only
    start_index = stop_index;
    stop_index += myAsl_->i.nlvoi_;
    for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
      vName = std::string(var_name_ASL(myAsl_, i));
      vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
          myAsl_->i.LUv_[2*i+1], Minotaur::Integer, vName);
      vars_.push_back(vPtr);
    }

    assert ((int) stop_index == myAsl_->i.nlvo_);
  } else {
    // all variables that are nonlinear in objective are also nonlinear in
    // constraints and hence we dont need to add anything more.
    assert ((int) stop_index == myAsl_->i.nlvc_);
  } // all nonlinear variables have been added.


  // we dont deal with linear arcs for now.
  assert(myAsl_->i.nwv_ == 0);

  // visit all linear continuous variables
  // continuous linear variables
  start_index = stop_index;
  stop_index  = myAsl_->i.n_var_ - (myAsl_->i.niv_ + myAsl_->i.nbv_);
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    vName = std::string(var_name_ASL(myAsl_, i));
    vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
        myAsl_->i.LUv_[2*i+1], Minotaur::Continuous, vName);
    vars_.push_back(vPtr);
  }

  // binary linear variables
  start_index = stop_index;
  stop_index += myAsl_->i.nbv_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    vName = std::string(var_name_ASL(myAsl_, i));
    vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
        myAsl_->i.LUv_[2*i+1], Minotaur::Binary, vName);
    vars_.push_back(vPtr);
  }

  // integer linear variables
  start_index = stop_index;
  stop_index += myAsl_->i.niv_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    vName = std::string(var_name_ASL(myAsl_, i));
    vPtr = instance->newVariable(myAsl_->i.LUv_[2*i], 
        myAsl_->i.LUv_[2*i+1], Minotaur::Integer, vName);
    vars_.push_back(vPtr);
  }

  assert ((int) stop_index == myAsl_->i.n_var_);
  // ALL variables have been added.
}



// Documentation of reading SOS information from the .nl files is available in
// the file README.suf in the asl directory.
void AMPLInterface::addSOS_(Minotaur::ProblemPtr instance)
{
  int flags = ASL_suf_sos_explict_free; // = caller will explicitly free
                                        // returned arrays.
  int   nsos    = 0;    // number of SOS constraints.
  int   nsosnz  = 0;    // total number of nonzeros in SOS constraints.
  char *sostype = 0;    // +1 Type-I, -1 Type-2
  int  *sospri  = 0;    
  int  *sosbeg  = 0;
  int  *sosind  = 0;
  real *sosref  = 0;
  int   copri[2];

  Minotaur::SOSType sostypem = Minotaur::SOS1;
  Minotaur::VarVector vars;

  logger_->msgStream(Minotaur::LogDebug2) << "Checking SOS information "
                                          << std::endl;
  copri[0] = 0;
  copri[1] = 0;
  nsos = suf_sos_ASL(myAsl_, flags, &nsosnz, &sostype, &sospri, copri, &sosbeg,
                     &sosind, &sosref);
  logger_->msgStream(Minotaur::LogDebug)  << "Number of SOS constraints = "
                                          << nsos << std::endl
                                          << "Number of SOS nonzeros = "
                                          << nsosnz << std::endl;
  
  for (int i=0; i<nsos; ++i) {
    if (sostype[i]=='1') {
      sostypem = Minotaur::SOS1;
    } else if (sostype[i]=='2') {
      sostypem = Minotaur::SOS2;
    } else {
      logger_->errStream() << "bad SOS type." << std::endl;
    }
    vars.clear();

    for (int j=sosbeg[i]; j<sosbeg[i+1]; ++j) {
      vars.push_back(instance->getVariable(sosind[j]));
    }

    instance->newSOS(sosbeg[i+1]-sosbeg[i], sostypem, sosref+sosbeg[i],
                     vars, sospri[i]);
  }

  if (0<nsos) {
    free(sosref);
    // All others are freed automatically.
    //free(sostype);
    //free(sospri);
    //free(sosbeg);
    //free(sosind);
  }
}


// copy the data desired data structures from ASL into Instance 
Minotaur::ProblemPtr AMPLInterface::copyInstanceFromASL_()
{
  Minotaur::UInt stop_index;
  Minotaur::UInt start_index;

  // new instance
  Minotaur::ProblemPtr instance = (Minotaur::ProblemPtr) 
                                      new Minotaur::Problem();
  addVariablesFromASL_(instance);
  addDefinedVars_(instance);

  // visit each constraint and copy the quadratic and linear parts.
  // add quadratic constraints 
  start_index = 0;
  stop_index = myAsl_->i.nlc_ - myAsl_->i.nlnc_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    addQuadraticConstraint_(i, instance);
  }

  // add constraints that are used to define 'defined variables'
  addQuadraticDefCons_(instance);

  // network constraints are not allowed.
  assert (myAsl_->i.nlnc_ == 0);
  assert (myAsl_->i.lnc_ == 0);

  // add linear constraints
  start_index = stop_index;
  stop_index  = myAsl_->i.n_con_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    addLinearConstraint_(i, instance);
  }

  // add objective. Do not accept more than one objective.
  assert (myAsl_->i.n_obj_ < 2);
  if (myAsl_->i.nlo_ > 0) {
    addQuadraticObjective_(0, instance);
  } else {
    addLinearObjective_(0, instance);
  }

  addSOS_(instance);
  //instance->calculateSize();
  //instance->writeSize();
  return instance;
}


Minotaur::ProblemPtr AMPLInterface::copyInstanceFromASL2_()
{
  Minotaur::UInt stop_index;
  Minotaur::UInt start_index;
  Minotaur::LinearFunctionPtr lf = Minotaur::LinearFunctionPtr();  //NULL
  Minotaur::QuadraticFunctionPtr qf = Minotaur::QuadraticFunctionPtr();  //NULL
  Minotaur::FunctionPtr f = Minotaur::FunctionPtr();  //NULL
  cde *constraint_cde = 0, *obj_cde = 0;
  ASL_fg *asl_fg = (ASL_fg *)myAsl_;
  Minotaur::CGraphPtr cgraph;
  Minotaur::CNode *cnode = 0;
  std::string name;
  Minotaur::ObjectiveType obj_sense = Minotaur::Minimize;
  double *x = 0, *grad = 0;
  int err = 0;

  // new instance
  Minotaur::ProblemPtr instance = (Minotaur::ProblemPtr) 
                                   new Minotaur::Problem();
  addVariablesFromASL_(instance);
  addDefinedVars_(instance);

  x = new double[instance->getNumVars()];
  grad = new double[instance->getNumVars()];
  memset(x, 0, instance->getNumVars()*sizeof(double));
  memset(grad, 0, instance->getNumVars()*sizeof(double));

  // visit each constraint and copy the linear parts and nonlinear parts.
  start_index = 0;
  stop_index = myAsl_->i.nlc_ - myAsl_->i.nlnc_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    cgraph = (Minotaur::CGraphPtr) new Minotaur::CGraph();
    lf = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    addLinearTermsFromConstr_(lf, i);
    constraint_cde = asl_fg->I.con_de_+i;
    cnode = getCGraph_(constraint_cde->e, cgraph, instance);
    cgraph->setOut(cnode);
    cgraph->finalize();
    assert(Minotaur::Constant!=cgraph->getType());
    if (Minotaur::Linear==cgraph->getType()) {
      if (!lf) {
        lf = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
      }
      cgraph->evalGradient(x, grad, &err);
      assert(0==err);
      for (Minotaur::UInt j=0; j<instance->getNumVars(); ++j) {
        if (fabs(grad[j])>1e-10) {
          lf->incTerm(instance->getVariable(j), grad[j]);
        }
      }
      memset(grad, 0, instance->getNumVars()*sizeof(double));
      cgraph.reset();
    }
    f = (Minotaur::FunctionPtr) new Minotaur::Function(lf, qf, cgraph);
    name = std::string(con_name_ASL(myAsl_, i));
    instance->newConstraint(f, myAsl_->i.LUrhs_[2*i], myAsl_->i.LUrhs_[2*i+1],
                            name); 
  }

  // add constraints that are used to define 'defined variables'
  addQuadraticDefCons2_(instance);

  // network constraints are not allowed.
  assert (myAsl_->i.nlnc_ == 0);
  assert (myAsl_->i.lnc_ == 0);

  // add linear constraints
  start_index = stop_index;
  stop_index  = myAsl_->i.n_con_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    addLinearConstraint_(i, instance);
  }

  assert (myAsl_->i.n_obj_ < 2);
  if (myAsl_->i.nlo_ > 0) {
    cgraph = (Minotaur::CGraphPtr) new Minotaur::CGraph();
    lf = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
    addLinearTermsFromObj_(lf, 0);
    obj_cde = asl_fg->I.obj_de_+0;
    assert (obj_cde);
    cnode = getCGraph_(obj_cde->e, cgraph, instance);
    cgraph->setOut(cnode);
    cgraph->finalize();
    name = std::string(obj_name_ASL(myAsl_, 0));
    f = (Minotaur::FunctionPtr) new Minotaur::Function(lf, qf, cgraph);
    if (myAsl_->i.objtype_[0] == 1) {
      obj_sense = Minotaur::Maximize;
    }
    instance->newObjective(f, objconst_ASL(myAsl_, 0), obj_sense, name);
  } else {
    addLinearObjective_(0, instance);
  }

  addSOS_(instance);

  delete [] grad;
  delete [] x;

  return instance;
}


void AMPLInterface::createFunctionMap_()
{
  // We want a map such that if an 'f' of type (efunc) is input, the
  // corresponding operation-code is output.

  // check if the map is already built.
  if (functionMap_.size() == N_OPS) {
    return;
  } else {
    // build the map.
    for (int i=0; i<N_OPS; ++i) {
      functionMap_[r_ops_ASL[i]] = i;
    }
  }
}


Minotaur::ProblemType AMPLInterface::findProblemType_()
{
  Minotaur::ProblemType problem_type = Minotaur::UnknownProblem;
  Minotaur::FunctionType obj_type  = Minotaur::UnknownFunction,
                         con_type  = Minotaur::UnknownFunction,
                         full_type = Minotaur::UnknownFunction;
  bool has_integers = false;

  // first check if it has integer constrained variables.
  if (myAsl_->i.nlvbi_ + myAsl_->i.nlvci_ + myAsl_->i.nlvoi_ + myAsl_->i.nbv_ 
      + myAsl_->i.niv_ > 0) {
    has_integers = true;
  }

  // check if any variable is nonlinear.
  if (myAsl_->i.nlvc_ + myAsl_->i.nlvo_ <= 0) {
    if (has_integers) {
      return Minotaur::MILP;
    } else {
      return Minotaur::LP;
    }
  }

  // Visit the objective function and check if it is a quadratic.
  if (myAsl_->i.n_obj_ == 0) {
    obj_type = Minotaur::Constant;
  } else {
    obj_type = getObjFunctionType_(0);
  }

  // Check if we can declare it a MINLP on the basis of objective function
  // alone and stop
  logger_->msgStream(Minotaur::LogDebug) << "Objective is " <<
    getFunctionTypeString(obj_type) << std::endl;
  if (obj_type == Minotaur::Nonlinear) {
    if (has_integers) {
      return Minotaur::MINLP;
    } else {
      return Minotaur::NLP;
    }
  }

  // Check constraints
  con_type = getConstraintsType_();
  full_type = funcTypesAdd(obj_type, con_type);
  switch (full_type) {
   case (Minotaur::Nonlinear):
    if (has_integers) {
      return Minotaur::MINLP;
    } else {
      return Minotaur::NLP;
    }
    break;
   case (Minotaur::Polynomial):
    if (has_integers) {
      return Minotaur::MIPOLYP;
    } else {
      return Minotaur::POLYP;
    }
    break;
   case (Minotaur::Quadratic):
    if (con_type==Minotaur::Linear) {
      if (has_integers) {
        return Minotaur::MIQP;
      } else {
        return Minotaur::QP;
      }
    } else if (has_integers) {
      return Minotaur::MIQCQP;
    } else {
      return Minotaur::QCQP;
    }
    break;
   case (Minotaur::Linear):
    if (has_integers) {
      return Minotaur::MILP;
    } else {
      return Minotaur::LP;
    }
    break;
   default:
    assert (!"Error finding constraints' type!");
  }

  return problem_type;
}


void AMPLInterface::findVars_(expr *e_ptr, std::set<int> & vars)
{
  int opcode = functionMap_.find(e_ptr->op)->second;
  switch (opcode) {
   case (OPPLUS):   // expr1 + expr2
   case (OPMINUS):  // expr1 - expr2
   case (OPMULT): // expr1*expr2
   case (OPDIV): // expr1/expr2
   case (OPREM): // remainder by dividing expr1/expr2
   case (OPPOW): // expr1 ^ expr2
     findVars_(e_ptr->L.e, vars);
     findVars_(e_ptr->R.e, vars);
     break;
   case (OPUMINUS):
     findVars_(e_ptr->L.e, vars);
     break;
   case (OPSUMLIST):
     expr **ep, **epe;
     ep = e_ptr->L.ep;
     epe = e_ptr->R.ep;
     while (ep < epe) {
       findVars_(*ep, vars);
       ++ep;
     }
     break;
   case (OP1POW): //  OPPOW for R = numeric constant
   case (OP2POW): //  expr^2
     findVars_(e_ptr->L.e, vars);
     break;
   case (OPCPOW): //  (constant)^expr
     findVars_(e_ptr->R.e, vars);
     break;
   case (OPNUM): //  numeric constant
       break;
   case (OPVARVAL): //  single variable
       assert((expr_v *)e_ptr - ((ASL_fg *)myAsl_)->I.var_e_ < nVars_);
       vars.insert((expr_v *)e_ptr - ((ASL_fg *)myAsl_)->I.var_e_);
       break;
   case (OP_tanh):  
   case (OP_tan):  
   case (OP_sqrt): 
   case (OP_sinh): 
   case (OP_sin):   
   case (OP_log10):
   case (OP_log):  
   case (OP_exp):   
   case (OP_cosh): 
   case (OP_cos):  
   case (OP_atanh):
   case (OP_atan2): 
   case (OP_atan): 
   case (OP_asinh):
   case (OP_asin): 
   case (OP_acosh): 
   case (OP_acos):
       findVars_(e_ptr->L.e, vars);
       break;
   case (OPLESS):      case (MINLIST):     case (MAXLIST):     case (FLOOR):
   case (CEIL):        case (ABS):         case (OPOR):        case (OPAND): 
   case (LT):          case (LE):          case (EQ):          case (GE):
   case (NE):          case (OPNOT):       case (OPIFnl):      case (OPintDIV):
   case (OPprecision): case (OPround):     case (OPtrunc):     case (OPCOUNT):
   case (OPNUMBEROF):  case (OPNUMBEROFs): case (OPATLEAST):   case (OPATMOST):
   case (OPPLTERM):    case (OPIFSYM):     case (OPEXACTLY):   
   case (OPNOTATLEAST):                    case (OPNOTATMOST): 
   case (OPNOTEXACTLY):                    case (ANDLIST):     case (ORLIST):
   case (OPIMPELSE):   case (OP_IFF):      case (OPALLDIFF):   case (OPFUNCALL):
       logger_->errStream() << "Operation code " << opcode
         << " Minotaur cannot handle this operation code."
         << std::endl;
     assert(!"cannot solve!");
     break;
  }
}


// free the myAsl_ data structure in ampl.
// it seems ampl does not free the initial point X0 itself. we have to free it
// ourselves.
void AMPLInterface::freeASL()
{
  if (myAsl_) {
    if (myAsl_->i.X0_) {
     free(myAsl_->i.X0_);
     myAsl_->i.X0_ = NULL;
    }
    ASL_free(&myAsl_);
  }
}


ASL * AMPLInterface::getAsl()
{
  return myAsl_;
}


Minotaur::CNode *AMPLInterface::getCGraph_(expr *e_ptr, 
                                           Minotaur::CGraphPtr cgraph,
                                           Minotaur::ProblemPtr instance)
{
  Minotaur::CNode *lchild = 0;
  Minotaur::CNode *rchild = 0;
  Minotaur::CNode **childr = 0;
  Minotaur::CNode *n = 0;
  int opcode = functionMap_[e_ptr->op];
  switch (opcode) {
  case (OPPLUS):   // expr1 + expr2
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    rchild = getCGraph_(e_ptr->R.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpPlus, lchild, rchild);
    break;
  case (OPMINUS):  
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    rchild = getCGraph_(e_ptr->R.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpMinus, lchild, rchild);
    break;
  case (OPMULT): // expr1*expr2
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    rchild = getCGraph_(e_ptr->R.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpMult, lchild, rchild);
    break;
  case (OPDIV): // expr1/expr2
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    rchild = getCGraph_(e_ptr->R.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpDiv, lchild, rchild);
    break;
  case (FLOOR):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpFloor, lchild, 0);
    break;
  case (CEIL):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpCeil, lchild, 0);
    break;
  case (ABS):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpAbs, lchild, 0);
    break;
  case (OPUMINUS):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpUMinus, lchild, 0);
    break;
  case (OP_tan):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpTan, lchild, 0);
    break;
  case (OP_tanh):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpTanh, lchild, 0);
    break;
  case (OP_sqrt):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpSqrt, lchild, 0);
    break;
  case (OP_sinh):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpSinh, lchild, 0);
    break;
  case (OP_sin):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpSin, lchild, 0);
    break;
  case (OP_log10):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpLog10, lchild, 0);
    break;
  case (OP_log):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpLog, lchild, 0);
    break;
  case (OP_exp):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpExp, lchild, 0);
    break;
  case (OP_cosh):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpCosh, lchild, 0);
    break;
  case (OP_cos):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpCos, lchild, 0);
    break;
  case (OP_atanh):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpTanh, lchild, 0);
    break;
  case (OP_atan2):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpAtan, lchild, 0);
    break;
  case (OP_atan):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpAtan, lchild, 0);
    break;
  case (OP_asinh):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpAsinh, lchild, 0);
    break;
  case (OP_asin):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpAsin, lchild, 0);
    break;
  case (OP_acosh):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpAcosh, lchild, 0);
    break;
  case (OP_acos):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpAcos, lchild, 0);
    break;
  case (OPSUMLIST):
    {
    //logger_->msgStream(Minotaur::LogNone) << " ++ " << std::endl;
    expr **ep, **epe;
    int i;
    ep = e_ptr->L.ep;
    epe = e_ptr->R.ep;
    i=0;
    while (ep < epe) {
      ++ep;
      ++i;
    }
    childr = new Minotaur::CNode *[i];
    ep = e_ptr->L.ep;
    epe = e_ptr->R.ep;
    i=0;
    while (ep < epe) {
      childr[i] = getCGraph_(*ep, cgraph, instance);
      ++ep;
      ++i;
    }
    n = cgraph->newNode(Minotaur::OpSumList, childr, i);
    delete []childr;
    }
    break;
  case (OPintDIV):
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    rchild = getCGraph_(e_ptr->R.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpIntDiv, lchild, rchild);
    break;
  case (OP1POW): 
    {
    //  OPPOW for R = numeric constant
    double power = ((expr_n *)e_ptr->R.e)->v;
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    rchild = cgraph->newNode(power);
    n = cgraph->newNode(Minotaur::OpPowK, lchild, rchild);
    }
    break;
  case (OP2POW): //  expr^2
    lchild = getCGraph_(e_ptr->L.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpSqr, lchild, 0);
    break;
  case (OPCPOW): //  (constant)^expr
    lchild = cgraph->newNode(((expr_n *)e_ptr->L.e)->v);
    rchild = getCGraph_(e_ptr->R.e, cgraph, instance);
    n = cgraph->newNode(Minotaur::OpCPow, lchild, rchild);
    break;
  case (OPNUM): //  numeric constant
    n = cgraph->newNode(((expr_n *)e_ptr)->v);
    break;
  case (OPVARVAL): //  single variable
    {
    int var_index = (expr_v *)e_ptr - ((ASL_fg *)myAsl_)->I.var_e_;
    //logger_->msgStream(Minotaur::LogNone) <<  "x" << var_index << std::endl;
    assert(var_index < nVars_+nDefVars_);
    n = cgraph->newNode(instance->getVariable(var_index));
    }
    break;
  default:
    unsupportedOp_(opcode);
  }
  return n;
}


Minotaur::FunctionType AMPLInterface::getConstraintsType_() 
{
  Minotaur::FunctionType function_type, overall_type;
  if (myAsl_->i.nlc_ == 0 && nDefVars_==0) {
    return Minotaur::Linear;
  }

  overall_type = Minotaur::Constant;
  // visit all nonlinear constraints that are not network-constraints and
  // check if all are Quadratic. 
  for (int i=0; i<myAsl_->i.nlc_ - myAsl_->i.nlnc_; ++i) {
    function_type = getConstraintType_(i);
    logger_->msgStream(Minotaur::LogDebug) << "Constraint (nonlin) is " <<
      getFunctionTypeString(function_type) << std::endl;
    overall_type = Minotaur::funcTypesAdd(overall_type, function_type);
  }


  // check the constraints that define the 'defined variables'
  for (int i=0; i<nDefVarsBco_; ++i) {
    function_type = getDefConstraintType_(i);
#if DEBUG
    logger_->msgStream(Minotaur::LogDebug) << "Constraint (def) " << i << 
      " is " << getFunctionTypeString(function_type) << std::endl;
#endif
    overall_type = Minotaur::funcTypesAdd(overall_type, function_type);
  }


  // check the constraints that define the 'defined variables' (only 1)
  for (int i=0; i<nDefVarsCo1_; ++i) {
    function_type = getDef1ConstraintType_(i);
#if DEBUG
    logger_->msgStream(Minotaur::LogDebug) << "Constraint (def1) " << i 
      << " is " << getFunctionTypeString(function_type) << std::endl;
#endif
    overall_type = Minotaur::funcTypesAdd(overall_type, function_type);
  }
  return overall_type;
}


Minotaur::FunctionType AMPLInterface::getConstraintType_(int i) 
{
  ASL_fg *asl_fg;
  cde *constraint_cde;

  if (i >= myAsl_->i.nlc_) {
    return Minotaur::Linear;
  }

  // first cast myAsl_ from ASL into ASL_fg
  asl_fg = (ASL_fg *)myAsl_;

  // get the constraint number i;
  constraint_cde = asl_fg->I.con_de_+i;

  assert (constraint_cde);
  // get the expression from this cde.
  return getExpressionType_(constraint_cde->e);
}


Minotaur::FunctionType AMPLInterface::getDef1ConstraintType_(int i) 
{
  Minotaur::FunctionType ftype;
  cexp1  *cstruct1; // cexp is defined in nlp.h of asl
  cstruct1 = (((const ASL_fg *)myAsl_)->I.cexps1_)+i;
  ftype = getExpressionType_(cstruct1->e);
  ftype = Minotaur::funcTypesAdd(ftype, Minotaur::Linear);
  return ftype;
}


Minotaur::FunctionType AMPLInterface::getDefConstraintType_(int i) 
{
  Minotaur::FunctionType ftype;
  cexp  *cstruct; // cexp is defined in nlp.h of asl
  cstruct = (((const ASL_fg *)myAsl_)->I.cexps_)+i;
  ftype = getExpressionType_(cstruct->e);
  ftype = Minotaur::funcTypesAdd(ftype, Minotaur::Linear);
  return ftype;
}


Minotaur::FunctionType AMPLInterface::getExpressionType_(expr *e_ptr)
{
  int opcode;
  Minotaur::FunctionType fun_type1, fun_type2;

  opcode = functionMap_[e_ptr->op];
  switch (opcode) {
   case (OPPLUS):   // expr1 + expr2
   case (OPMINUS):  // expr1 - expr2
     return getPMExpressionType_(e_ptr);
     break;
   case (OPMULT): // expr1*expr2
     return getMultExpressionType_(e_ptr);
     break;
   case (OPDIV): // expr1/expr2
     fun_type2 = getExpressionType_(e_ptr->R.e);
     if (fun_type2 == Minotaur::Constant) {
       fun_type1 = getExpressionType_(e_ptr->L.e);
       return fun_type1;
     } else {
       return Minotaur::Nonlinear;
     }
     break;
   case (OPREM): // remainder by dividing expr1/expr2
   case (OPPOW): // expr1 ^ expr2
     fun_type1 = getExpressionType_(e_ptr->L.e);
     fun_type2 = getExpressionType_(e_ptr->R.e);
     if (fun_type1 == Minotaur::Constant && fun_type2 == Minotaur::Constant) {
       return Minotaur::Constant;
     } else {
       return Minotaur::Nonlinear;
     }
     break;
   case (OPUMINUS):
     return getExpressionType_(e_ptr->L.e);
     break;
   case (OPSUMLIST):
     return getSumlistExpressionType_(e_ptr);
     break;
   case (OP1POW): //  OPPOW for R = numeric constant
     return getPow1ExpressionType_(e_ptr);
     break;
   case (OP2POW): //  expr^2
     fun_type1 = getExpressionType_(e_ptr->L.e);
     if (fun_type1 == Minotaur::Constant) {
       return Minotaur::Constant;
     } else if (fun_type1 == Minotaur::Linear) {
       return Minotaur::Quadratic;
     } else if (fun_type1 == Minotaur::Quadratic || 
                fun_type1 == Minotaur::Polynomial) {
       return Minotaur::Polynomial;
     } else {
       return Minotaur::Nonlinear;
     }
     break;
   case (OPCPOW): //  (constant)^expr
     fun_type2 = getExpressionType_(e_ptr->R.e);
     if (fun_type2 == Minotaur::Constant) {
       return Minotaur::Constant;
     } else {
       return Minotaur::Nonlinear;
     }
     break;
   case (OPNUM): //  numeric constant
       return Minotaur::Constant;
       break;
   case (OPVARVAL): //  single variable
       return Minotaur::Linear;
       break;
   case (OP_tanh):  case (OP_tan):   case (OP_sqrt):  case (OP_sinh): 
   case (OP_sin):   case (OP_log10): case (OP_log):   case (OP_exp):   
   case (OP_cosh):  case (OP_cos):   case (OP_atanh): case (OP_atan2): 
   case (OP_atan):  case (OP_asinh): case (OP_asin):  case (OP_acosh): 
   case (OP_acos):
     fun_type1 = getExpressionType_(e_ptr->L.e);
     if (fun_type1 == Minotaur::Constant) {
       return Minotaur::Constant;
     } else {
       return Minotaur::Nonlinear;
     }
     break;
   case (OPLESS):      case (MINLIST):     case (MAXLIST):     case (FLOOR):
   case (CEIL):        case (ABS):         case (OPOR):        case (OPAND): 
   case (LT):          case (LE):          case (EQ):          case (GE):
   case (NE):          case (OPNOT):       case (OPIFnl):      case (OPintDIV):
   case (OPprecision): case (OPround):     case (OPtrunc):     case (OPCOUNT):
   case (OPNUMBEROF):  case (OPNUMBEROFs): case (OPATLEAST):   case (OPATMOST):
   case (OPPLTERM):    case (OPIFSYM):     case (OPEXACTLY):   
   case (OPNOTATLEAST):                    case (OPNOTATMOST): 
   case (OPNOTEXACTLY):                    case (ANDLIST):     case (ORLIST):
   case (OPIMPELSE):   case (OP_IFF):      case (OPALLDIFF):   case (OPFUNCALL):
     logger_->errStream() << "Operation code " << opcode
                          << " Minotaur cannot handle this operation code."
                          << std::endl;
     return Minotaur::Nonlinear;
     break;
  }

  return Minotaur::Nonlinear;
}


const double * AMPLInterface::getInitialPoint() const
{
  if (myAsl_->i.X0_) {
    return myAsl_->i.X0_;
  } else {
    return NULL;
  }
}


// copy the data desired data structures from ASL into Instance 
// TODO: XXX: break this long function into smaller pieces.
Minotaur::ProblemPtr AMPLInterface::getInstanceFromASL_(
    std::vector<std::set<int> > &vids)
{
  Minotaur::LinearFunctionPtr lfPtr;
  AMPLNlfPtr nlfPtr;
  Minotaur::FunctionPtr fPtr;
  Minotaur::QuadraticFunctionPtr qfPtr = Minotaur::QuadraticFunctionPtr(); 
  std::vector<Minotaur::VariableSet> vars;
  Minotaur::VariableSet vset;

  Minotaur::ObjectiveType obj_sense = Minotaur::Minimize;
  std::string cName, oName;

  Minotaur::UInt stop_index;
  Minotaur::UInt start_index;

  Minotaur::UInt n = myAsl_->i.n_var_;
  Minotaur::VariablePtr v;


  // new instance
  Minotaur::ProblemPtr instance = (Minotaur::ProblemPtr) 
    new Minotaur::Problem();

  addVariablesFromASL_(instance);
  nDefVars_    = 0;
  nDefVarsBco_ = 0;
  nDefVarsCo1_ = 0; 

  // convert vids to vars
  for (std::vector<std::set<int> >::iterator it=vids.begin();
      it!=vids.end(); ++it) {
    for (std::set<int>::iterator it2=it->begin(); it2!=it->end(); ++it2) {
      v = instance->getVariable(*it2);
      vset.insert(v);
    }
    vars.push_back(vset);
    vset.clear();
  }




  // add constraints now
  
  //
  // Ordering of constraints by ampl:
  // http://www.ampl.com/REFS/HOOKING/index.html
  //
  // category           count
  // -------------------------------------------------------------------------
  // nonlinear-general  nlc - nlnc
  // nonlinear-network  nlnc
  // linear-network     lnc
  // linear-general     n_con - (nlc + lnc)
  //

  // add non-linear constraints 
  lfPtr = Minotaur::LinearFunctionPtr();  //NULL
  start_index = 0;
  stop_index = myAsl_->i.nlc_ - myAsl_->i.nlnc_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    cName = std::string(con_name_ASL(myAsl_, i));
    nlfPtr = (AMPLNlfPtr) new AMPLNonlinearFunction(i, n, myAsl_, false);
    if (vars.size()>0) {
      nlfPtr->setVars(vars[i+1].begin(), vars[i+1].end());
    }
    fPtr = (Minotaur::FunctionPtr) new Minotaur::Function(lfPtr, qfPtr, nlfPtr);
    instance->newConstraint(fPtr, myAsl_->i.LUrhs_[2*i], 
        myAsl_->i.LUrhs_[2*i+1], cName); 
  }

  assert (myAsl_->i.nlnc_ == 0);
  assert (myAsl_->i.lnc_ == 0);

  // add linear constraints
  start_index = stop_index;
  stop_index  = myAsl_->i.n_con_;
  for (Minotaur::UInt i=start_index; i<stop_index; ++i) {
    addLinearConstraint_(i, instance);
  }

  // add objective 
  assert (myAsl_->i.n_obj_ < 2);

  if (myAsl_->i.n_obj_ > 0) {
    if (myAsl_->i.objtype_[0] == 1) {
      obj_sense = Minotaur::Maximize;
    } else {
      obj_sense = Minotaur::Minimize;
    }

    if (myAsl_->i.nlo_ > 0) {
      oName = std::string(obj_name_ASL(myAsl_, 0));
      nlfPtr = (AMPLNlfPtr) new AMPLNonlinearFunction(0, n, myAsl_, true);
      if (vars.size()>0) {
        nlfPtr->setVars(vars[0].begin(), vars[0].end());
      }
      lfPtr = Minotaur::LinearFunctionPtr();  //NULL
      fPtr = (Minotaur::FunctionPtr) new 
        Minotaur::Function(lfPtr, qfPtr, nlfPtr);
      instance->newObjective(fPtr, objconst_ASL(myAsl_, 0), obj_sense, 
          oName);
    } else {
      // linear objective
      addLinearObjective_(0, instance);
    }
  }

  addSOS_(instance);

  // all done.
  return instance;
}


Minotaur::FunctionType AMPLInterface::getMultExpressionType_(expr *e_ptr) 
{
  Minotaur::FunctionType f1, f2;
  f1 = getExpressionType_(e_ptr->L.e);
  f2 = getExpressionType_(e_ptr->R.e);
  switch (f1) {
   case (Minotaur::Constant):
     return f2;
   case (Minotaur::Linear):
     if (f2 == Minotaur::Constant) {
       return Minotaur::Linear;
     } else if (f2 == Minotaur::Linear) {
       return Minotaur::Quadratic;
     } else if (f2 == Minotaur::Quadratic ||
                f2 == Minotaur::Polynomial) {
       return Minotaur::Polynomial;
     } else {
       return Minotaur::Nonlinear;
     }
   case (Minotaur::Quadratic):
     if (f2 == Minotaur::Constant) {
       return Minotaur::Quadratic;
     } else if (f2 == Minotaur::Linear || 
                f2 == Minotaur::Quadratic ||
                f2 == Minotaur::Polynomial) {
       return Minotaur::Polynomial;
     } else {
       return Minotaur::Nonlinear;
     }
   case (Minotaur::Polynomial):
     if (f2 == Minotaur::Constant ||
         f2 == Minotaur::Linear || 
         f2 == Minotaur::Quadratic ||
         f2 == Minotaur::Polynomial) {
       return Minotaur::Polynomial;
     }
   default:
     return Minotaur::Nonlinear;
  }
}


Minotaur::UInt AMPLInterface::getNumDefs() const
{
  return nDefVars_;
}


Minotaur::FunctionType 
AMPLInterface::getObjFunctionType_(Minotaur::UInt obj_index) 
{
  // first cast myAsl_ from ASL into ASL_fg
  ASL_fg *asl_fg = (ASL_fg *)myAsl_;
  // get the objective number obj_index
  cde *obj_cde = asl_fg->I.obj_de_+obj_index;

  assert (obj_cde);
  // get the expression from this cde.
  return getExpressionType_(obj_cde->e);
}


void AMPLInterface::getOptionsFromEnv_(std::string pre)
{
  std::string str;
  char *cval=0, *c=0, *cc=0;
  char **argv=0;
  std::vector<char *>argvec;
  size_t argc=1;

  str = pre+"_options";
  cval = getenv(str.c_str());
  if (cval) {
    // first word is just the name of the solver.
    cc = new char[pre.length()+1];
    strcpy (cc,pre.c_str());
    argvec.push_back(cc);

    // now add rest of the arguments.
    c = strtok (cval," ,\t\n");
    while (c != NULL) {
      ++argc;
      cc = new char[strlen(c)+1];
      memcpy (cc,c,strlen(c)+1);
      argvec.push_back(cc);
      c = strtok(NULL, "  ,\t\n");
    }

    // Set up argc, argv. Send them to env.
    assert(argc==argvec.size());
    argv = new char*[argc];
    for (size_t i=0; i<argc; ++i) {
      argv[i] = argvec[i];
    }
    env_->readOptions(argc, argv);

    for (size_t i=0; i<argc; ++i) {
      delete [] argv[i];
    }
    delete [] argv;
  }
}


Minotaur::FunctionType AMPLInterface::getPMExpressionType_(expr *e_ptr) 
{
  // In case of addition and subtraction, if any of the two operands 
  // is not quadratic then e_ptr is not quadratic
  Minotaur::FunctionType f1, f2;
  f1 = getExpressionType_(e_ptr->L.e);
  f2 = getExpressionType_(e_ptr->R.e);
  return Minotaur::funcTypesAdd(f1, f2);
}


void AMPLInterface::getPoly_(Minotaur::LinearFunctionPtr & lfPtr, 
                             Minotaur::QuadraticFunctionPtr & qfPtr, 
                             Minotaur::PolyFunPtr & pfPtr, 
                             double & c, expr *e_ptr)
{
  int opcode;
  Minotaur::LinearFunctionPtr lfPtr1 = Minotaur::LinearFunctionPtr();
  Minotaur::LinearFunctionPtr lfPtr2 = Minotaur::LinearFunctionPtr();
  Minotaur::QuadraticFunctionPtr qfPtr1 = Minotaur::QuadraticFunctionPtr();
  Minotaur::QuadraticFunctionPtr qfPtr2 = Minotaur::QuadraticFunctionPtr();
  Minotaur::PolyFunPtr pfPtr1 = Minotaur::PolyFunPtr();
  Minotaur::PolyFunPtr pfPtr2 = Minotaur::PolyFunPtr();
  lfPtr = Minotaur::LinearFunctionPtr();  //NULL
  qfPtr = Minotaur::QuadraticFunctionPtr();  //NULL
  double c1 = 0, c2 = 0;
  c = 0;
  int var_index;

  opcode = functionMap_[e_ptr->op];
  switch (opcode) {
   case (OPPLUS):   // expr1 + expr2
     //logger_->msgStream(Minotaur::LogNone) << " + " << std::endl;
     getPoly_(lfPtr1, qfPtr1, pfPtr1, c1, e_ptr->L.e);
     getPoly_(lfPtr2, qfPtr2, pfPtr2, c2, e_ptr->R.e);
     c = c1 + c2;
     lfPtr = lfPtr1 + lfPtr2;
     qfPtr = qfPtr1 + qfPtr2;
     if (pfPtr1) {
       pfPtr = pfPtr1;
       (*pfPtr) += pfPtr2;
     } else {
       pfPtr = pfPtr2;
     }
     break;
   case (OPMINUS):  // expr1 - expr2
     //logger_->msgStream(Minotaur::LogNone) << " - " << std::endl;
     getPoly_(lfPtr1, qfPtr1, pfPtr1, c1, e_ptr->L.e);
     getPoly_(lfPtr2, qfPtr2, pfPtr2, c2, e_ptr->R.e);
     c = c1 - c2;
     lfPtr = lfPtr1 - lfPtr2;
     qfPtr = qfPtr1 - qfPtr2;
     if (pfPtr2) {
       pfPtr = pfPtr2;
       (*pfPtr) *= -1;
       (*pfPtr) += pfPtr1;
     } else {
       pfPtr = pfPtr1;
     }
     break;
   case (OPMULT): // expr1*expr2
     // logger_->msgStream(Minotaur::LogNone) << " * " << std::endl;
     getPoly_(lfPtr1, qfPtr1, pfPtr1, c1, e_ptr->L.e);
     getPoly_(lfPtr2, qfPtr2, pfPtr2, c2, e_ptr->R.e);
     c = c1*c2;
     lfPtr = c1*lfPtr2 + c2*lfPtr1;
     qfPtr = lfPtr1*lfPtr2 + c1*qfPtr2 + c2*qfPtr1;
     pfPtr = c1*pfPtr2 + c2*pfPtr1 + qfPtr2*lfPtr1 + qfPtr1*lfPtr2 + 
        pfPtr2*lfPtr1 + pfPtr1*lfPtr2 + qfPtr1*qfPtr2 + pfPtr2*qfPtr1 + 
        pfPtr1*qfPtr2 + pfPtr1*pfPtr2;
     break;
   case (OPDIV): // expr1/expr2
     //logger_->msgStream(Minotaur::LogNone) << " / " << std::endl;
     getPoly_(lfPtr1, qfPtr1, pfPtr1, c1, e_ptr->L.e);
     c2 = ((expr_n *)e_ptr->R.e)->v;
     assert (fabs(c2)>zTol_);
     c     = c1/c2;
     c2    = 1/c2;
     lfPtr = c2*lfPtr1;
     qfPtr = c2*qfPtr1;
     pfPtr = c2*pfPtr1;
     break;
   case (OPREM): // remainder by dividing expr1/expr2
     assert(!"OPREM not implemented yet!");
     break;
   case (OPPOW): // expr1 ^ expr2
     assert(!"OPPOW not implemented yet!");
     break;
   case (OPUMINUS):
     getPoly_(lfPtr1, qfPtr1, pfPtr1, c1, e_ptr->L.e);
     if (lfPtr1) {
       lfPtr = -1.0 * lfPtr1;
     } 
     if (qfPtr1) {
       qfPtr = -1.0 * qfPtr1;
     }
     if (pfPtr1) {
       pfPtr = -1.0 * pfPtr1;
     }
     c = -1.0*c1;
     break;
   case (OPSUMLIST):
     //logger_->msgStream(Minotaur::LogNone) << " ++ " << std::endl;
     expr **ep, **epe;
     ep = e_ptr->L.ep;
     epe = e_ptr->R.ep;
     lfPtr = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
     qfPtr = (Minotaur::QuadraticFunctionPtr) new Minotaur::QuadraticFunction();
     pfPtr = (Minotaur::PolyFunPtr) new Minotaur::PolynomialFunction();
     while (ep < epe) {
       getPoly_(lfPtr2, qfPtr2, pfPtr2, c2, *ep);
       (*lfPtr) += lfPtr2;
       (*qfPtr) += qfPtr2;
       (*pfPtr) += pfPtr2;
              c += c2;
       ++ep;
     }
     if (lfPtr->getNumTerms() == 0) {
       lfPtr = Minotaur::LinearFunctionPtr(); //NULL
     }
     if (qfPtr->getNumTerms() == 0) {
       qfPtr = Minotaur::QuadraticFunctionPtr(); //NULL
     }
     if (pfPtr->isEmpty()) {
       pfPtr = Minotaur::PolyFunPtr(); //NULL
     }
     break;
   case (OP1POW): //  OPPOW for R = numeric constant
     {
       double power = ((expr_n *)e_ptr->R.e)->v;
       int ipower = -1;
       assert (power > 2 && fabs(power - floor(power+0.5)) < intTol_); 
       ipower = (int) power;
       getPoly_(lfPtr1, qfPtr1, pfPtr1, c1, e_ptr->L.e);
       pfPtr2 = (Minotaur::PolyFunPtr) new Minotaur::PolynomialFunction();
       (*pfPtr2) += c1;
       (*pfPtr2) += lfPtr1;
       (*pfPtr2) += qfPtr1;
       (*pfPtr2) += pfPtr1;
       pfPtr = pfPtr2->clone();
       for (int i=0; i<ipower-1; ++i) {
         (*pfPtr) *= pfPtr2;
       }
       break;
     }

     assert(!"OP1POW not implemented yet!");
   case (OP2POW): //  expr^2
     //logger_->msgStream(Minotaur::LogNone) <<  "^2" << std::endl;
     getPoly_(lfPtr1, qfPtr1, pfPtr1, c1, e_ptr->L.e);

     c = c1*c1;
     lfPtr = (2.0*c1) * lfPtr1;
     qfPtr = lfPtr1->square() + (2*c1)*qfPtr1;
     pfPtr = pfPtr1*pfPtr1 + qfPtr1*qfPtr1 + 2.*pfPtr1*qfPtr1 +
       2.*pfPtr1*lfPtr1 + 2.*qfPtr1*lfPtr1 + (2.*c1)*pfPtr1;
     break;
   case (OPCPOW): //  (constant)^expr
     assert(!"OPCPOW not implemented yet!");
   case (OPNUM): //  numeric constant
     c = ((expr_n *)e_ptr)->v;
     //logger_->msgStream(Minotaur::LogNone) << "constant = " <<  c  << std::endl;
     break;
   case (OPVARVAL): //  single variable
     // not sure if this var_index is correct. no documentation available!
     var_index = (expr_v *)e_ptr - ((ASL_fg *)myAsl_)->I.var_e_;
     //logger_->msgStream(Minotaur::LogNone) <<  "x" << var_index << std::endl;
     assert(var_index < nVars_+nDefVars_);
     lfPtr = (Minotaur::LinearFunctionPtr) new Minotaur::LinearFunction();
     lfPtr->addTerm(vars_[var_index], 1.0);
     break;
   default:
     logger_->msgStream(Minotaur::LogError) << "ASL opcode " << opcode 
       << std::endl;
     assert(!"can not recover function from ASL!");
     break;
  }
}


Minotaur::FunctionType AMPLInterface::getPow1ExpressionType_(expr *e_ptr) 
{
  Minotaur::FunctionType f1;
  double power = ((expr_n *)e_ptr->R.e)->v;
  f1 = getExpressionType_(e_ptr->L.e);
  if ( fabs(power) <= zTol_) {
    return Minotaur::Constant;
  } else if (fabs(power - 1.0) <= zTol_) {
    return f1;
  } else if (f1 == Minotaur::Constant) {
    return Minotaur::Constant;
  } else if (f1 == Minotaur::Linear) {
    if ( fabs(power - 2.0) <= zTol_ ) {
      return Minotaur::Quadratic;
    } else if ( fabs(power - floor(power+0.5)) < intTol_ ) {
      return Minotaur::Polynomial;
    } else {
      return Minotaur::Nonlinear;
    }
  } else if (f1 == Minotaur::Quadratic || f1 == Minotaur::Polynomial) {
    if (fabs(power - floor(power+0.5)) < intTol_) {
      return Minotaur::Polynomial;
    } else {
      return Minotaur::Nonlinear;
    }
  }
  return Minotaur::Nonlinear;
}


ReaderType AMPLInterface::getReaderType()
{
    return readerType_;
}


Minotaur::FunctionType AMPLInterface::getSumlistExpressionType_(expr *e_ptr) 
{
  Minotaur::FunctionType f1, f2;
  expr **ep, **epe;
  ep = e_ptr->L.ep;
  epe = e_ptr->R.ep;
  f2 = Minotaur::Constant;
  while (ep < epe) {
    f1 = getExpressionType_(*ep);
    f2 = Minotaur::funcTypesAdd(f1, f2);
    if (f2 == Minotaur::Nonlinear) {
      break;
    }
    ++ep;
  }
  return f2;
}


// Read the instance from a '.nl' file.
// XXX: TODO: Use AMPLMpsInterface derived class to read MPS/LP files 
void AMPLInterface::readFile_(std::string *fname, ReaderType readerType) 
{
  FILE *nl = NULL;
  char *fname_chars;

  assert(fname);

  fname_chars = (char *)malloc((fname->length()+1)*sizeof(char));
  strcpy(fname_chars, fname->c_str());

  // set flag for what reader type is now in use.
  // XXX: TODO: all asl routines called in the file should check what reader was
  // used.
  readerType_ = readerType;

  switch (readerType) {
   case (FReader):
     myAsl_ = ASL_alloc(ASL_read_f); 
     break;
   case (FGReader):
     myAsl_ = ASL_alloc(ASL_read_fg); 
     break;
   case (FGHReader):
     myAsl_ = ASL_alloc(ASL_read_fgh); 
     break;
   case (PFGReader):
     myAsl_ = ASL_alloc(ASL_read_pfg); 
     break;
   case (PFGHReader):
     myAsl_ = ASL_alloc(ASL_read_pfgh); 
     break;
  } 

  // initialization for linear program 
  nl = jac0dim_ASL(myAsl_, fname_chars, (fint) (fname->length()));

  free(fname_chars);

  // no networks allowed yet 
  assert(myAsl_->i.nwv_ == 0);

  // only one objective at most
  assert(myAsl_->i.n_obj_ < 2);

  // number of variables
  nVars_ = myAsl_->i.n_var_;

  // common expressions and defined variables
  nDefVarsBco_ = myAsl_->i.comb_ + myAsl_->i.comc_ + myAsl_->i.como_; 
  nDefVarsCo1_ = myAsl_->i.comc1_ + myAsl_->i.como1_; 
  nDefVars_    = nDefVarsBco_ + nDefVarsCo1_;

  // number of constraints
  nCons_ = myAsl_->i.n_con_;

  // setup space for initial guess
  myAsl_->i.X0_ = (real *)mymalloc_ASL(nVars_*sizeof(real));

  // Tell ASL about suffixes we want.
  suf_declare_ASL(myAsl_, suftab, sizeof(suftab) / sizeof(SufDecl));

  // read the full program 
  switch (readerType) {
   case (FReader):
     f_read_ASL(myAsl_,nl,0);
     break;
   case (FGReader):
     fg_read_ASL(myAsl_,nl,0);
     break;
   case (FGHReader):
     fgh_read_ASL(myAsl_,nl,0);
     break;
   case (PFGReader):
     pfg_read_ASL(myAsl_,nl,0);
     break;
   case (PFGHReader):
     pfgh_read_ASL(myAsl_,nl,0);
     break;
  }
}


Minotaur::ProblemPtr AMPLInterface::readInstance(std::string fname) 
{
  if (false==env_->getOptions()->findBool("use_native_cgraph")->getValue()) {
    return readInstanceASL_(fname);
  } 
  return readInstanceCG_(fname);
}


// read the file and build the instance 
Minotaur::ProblemPtr AMPLInterface::readInstanceASL_(std::string fname) 
{
  Minotaur::ProblemType problem_type = Minotaur::UnknownProblem;
  Minotaur::ProblemPtr instance;
  std::vector<std::set<int> > vars;
  bool do_poly = env_->getOptions()->findBool("expand_poly")->getValue();

  // We will read the stub using AMPL upto 2 times. In the first pass, read
  // the stub and find the problem type. In the next pass, we actually save
  // the data. If the problem is quadratic or polynomial, we save in native
  // format, otherwise the problem just keeps a pointer to ampl solver library
  // and uses information from there. 

  // ask AMPL to read the stub 
  readFile_(&fname, FGReader);

  // Load the array of opcodes
  createFunctionMap_();

  if (env_->getOptions()->findBool("display_ampl_model")->getValue()) {
    writeProblem(logger_->msgStream(Minotaur::LogNone));
  }

  problem_type = findProblemType_();

  // create a new instance 
  switch (problem_type) {
   case (Minotaur::LP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as LP" << std::endl;
     break;
   case (Minotaur::MILP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as MILP" << std::endl;
     break;
   case (Minotaur::QP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as QP" << std::endl;
     break;
   case (Minotaur::MIQP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as MIQP" << std::endl;
     break;
   case (Minotaur::QCQP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as QCQP" << std::endl;
     break;
   case (Minotaur::MIQCQP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as MIQCQP" << std::endl;
     break;
   case (Minotaur::POLYP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as POLYNOMIAL Program" << std::endl;
     break;
   case (Minotaur::MIPOLYP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as MI " << "POLYNOMIAL Program" << std::endl;
     break;
   case (Minotaur::NLP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as NLP" << std::endl;
     break;
   case (Minotaur::MINLP):
     logger_->msgStream(Minotaur::LogInfo) << me_ 
     << "Problem identified as MINLP" << std::endl;
     break;
   case (Minotaur::UnknownProblem):
    break;
  }

  if (problem_type==Minotaur::NLP || problem_type==Minotaur::MINLP ||
      ((problem_type==Minotaur::POLYP || problem_type==Minotaur::MIPOLYP) &&
       do_poly == false) ||
      ((problem_type==Minotaur::QP || problem_type==Minotaur::QCQP ||
       problem_type==Minotaur::MIQP || problem_type==Minotaur::MIQCQP) && 
      env_->getOptions()->findBool("expand_quad")->getValue()==false)){
    if (nDefVars_<1) {
      // identify variables in each nonlinear constraint and objective.
      saveNlVars_(vars);
    }
    // reset the reader and use PFGHReader
    freeASL();
    readFile_(&fname, PFGHReader);
    instance = getInstanceFromASL_(vars);
  } else {
    instance = copyInstanceFromASL_();
  }
  return instance;
}


Minotaur::ProblemPtr AMPLInterface::readInstanceCG_(std::string fname) 
{
  Minotaur::ProblemPtr instance;
  std::vector<std::set<int> > vars;

  // We will read the stub using AMPL upto 2 times. In the first pass, read
  // the stub and find the problem type. In the next pass, we actually save
  // the data. If the problem is quadratic or polynomial, we save in native
  // format, otherwise the problem just keeps a pointer to ampl solver library
  // and uses information from there. 

  // ask AMPL to read the stub 
  readFile_(&fname, FGReader);

  // Load the array of ASL opcodes
  createFunctionMap_();

  instance = copyInstanceFromASL2_();

  env_->getLogger()->msgStream(Minotaur::LogInfo) << me_ << "problem type is "
    << getProblemTypeString(instance->findType()) << std::endl;
  //instance->write(std::cout);
  return instance;
}


void AMPLInterface::saveNlVars_(std::vector<std::set<int> > &vars)
{
  std::set<int> vset;
  ASL_fg *asl_fg = (ASL_fg *)myAsl_;
  cgrad *cg; 
  ograd *og;

  findVars_((asl_fg->I.obj_de_+0)->e, vset);
  for (og = myAsl_->i.Ograd_[0]; og; og = og->next) {
    vset.insert(og->varno);
  }
  vars.push_back(vset);
  vset.clear();

  for (int i=0; i<myAsl_->i.nlc_ - myAsl_->i.nlnc_; ++i) {
    findVars_((asl_fg->I.con_de_+i)->e, vset);
    for (cg = myAsl_->i.Cgrad_[i]; cg; cg = cg->next) {
      vset.insert(cg->varno);
    }
    vars.push_back(vset);
    vset.clear();
  }
}


void AMPLInterface::unsupportedOp_(int opcode)
{
  switch(opcode) {
  case (OPREM):
    logger_->errStream() << "opcode OPREM is unsupported!" << std::endl;
    break;
  case (OPLESS):
    logger_->errStream() << "opcode OPLESS is unsupported!" << std::endl;
    break;
  case (MINLIST):
    logger_->errStream() << "opcode OPMINLIST is unsupported!" << std::endl;
    break;
  case (MAXLIST):
    logger_->errStream() << "opcode OPMAXLIST is unsupported!" << std::endl;
    break;
  case (OPOR):
    logger_->errStream() << "opcode OPOR is unsupported!" << std::endl;
    break;
  case (OPAND):
    logger_->errStream() << "opcode OPAND is unsupported!" << std::endl;
    break;
  case (LT):
    logger_->errStream() << "opcode OPLT is unsupported!" << std::endl;
    break;
  case (LE):
    logger_->errStream() << "opcode OPLE is unsupported!" << std::endl;
    break;
  case (EQ):
    logger_->errStream() << "opcode OPEQ is unsupported!" << std::endl;
    break;
  case (GE):
    logger_->errStream() << "opcode OPGE is unsupported!" << std::endl;
    break;
  case (GT):
    logger_->errStream() << "opcode OPGT is unsupported!" << std::endl;
    break;
  case (NE):
    logger_->errStream() << "opcode OPNE is unsupported!" << std::endl;
    break;
  case (OPNOT):
    logger_->errStream() << "opcode OPNOT is unsupported!" << std::endl;
    break;
  case (OPIFnl):
    logger_->errStream() << "opcode OPIFnl is unsupported!" << std::endl;
    break;
  case (OPprecision):
    logger_->errStream() << "opcode OPprecision is unsupported!" << std::endl;
    break;
  case (OPround):
    logger_->errStream() << "opcode OPround is unsupported!" << std::endl;
    break;
  case (OPtrunc):
    logger_->errStream() << "opcode OPtrunc is unsupported!" << std::endl;
    break;
  case (OPCOUNT):
    logger_->errStream() << "opcode OPCOUNT is unsupported!" << std::endl;
    break;
  case (OPNUMBEROFs):
    logger_->errStream() << "opcode OPNIMBEROFs is unsupported!" << std::endl;
    break;
  case (OPATLEAST):
    logger_->errStream() << "opcode OPATLEAST is unsupported!" << std::endl;
    break;
  case (OPATMOST):
    logger_->errStream() << "opcode OPATMOST is unsupported!" << std::endl;
    break;
  case (OPPLTERM):
    logger_->errStream() << "opcode OPLTERM is unsupported!" << std::endl;
    break;
  case (OPIFSYM):
    logger_->errStream() << "opcode OPIFSYM is unsupported!" << std::endl;
    break;
  case (OPEXACTLY):
    logger_->errStream() << "opcode OPEXACTLY is unsupported!" << std::endl;
    break;
  case (OPNOTATLEAST):
    logger_->errStream() << "opcode OPNOTATLEAST is unsupported!" << std::endl;
    break;
  case (OPNOTATMOST):
    logger_->errStream() << "opcode OPNOTATMOST is unsupported!" << std::endl;
    break;
  case (OPNOTEXACTLY):
    logger_->errStream() << "opcode OPNOTEXACTLY is unsupported!" << std::endl;
    break;
  case (ANDLIST):
    logger_->errStream() << "opcode OPANDLIST is unsupported!" << std::endl;
    break;
  case (ORLIST):
    logger_->errStream() << "opcode OPLIST is unsupported!" << std::endl;
    break;
  case (OPIMPELSE):
    logger_->errStream() << "opcode OPIMPELSE is unsupported!" << std::endl;
    break;
  case (OP_IFF):
    logger_->errStream() << "opcode OP_IFF is unsupported!" << std::endl;
    break;
  case (OPALLDIFF):
    logger_->errStream() << "opcode OPALLDIFF is unsupported!" << std::endl;
    break;
  case (OPHOL):
    logger_->errStream() << "opcode OPHOL is unsupported!" << std::endl;
    break;
  default: 
    logger_->errStream() << "opcode " << opcode << " is unknown!" << std::endl;
  }
  assert(!"unsupported opcode!");
}


void AMPLInterface::writeExpression_(int i, bool is_obj, std::ostream &out)
  const
{
  ASL_fg *asl_fg = (ASL_fg *)myAsl_;
  if (is_obj) {
    if (asl_fg->I.obj_de_->e) {
      writeExpression_(asl_fg->I.obj_de_->e, out);
    } else {
      out << "nonlinear function";
    }
  } else {
    if (i<=myAsl_->i.nlc_ && asl_fg->I.con_de_->e) {
      writeExpression_((asl_fg->I.con_de_+i)->e, out);
    } else {
      out << "nonlinear function";
    }
  }
}


void AMPLInterface::writeExpression_(expr *e_ptr, std::ostream &out) const
{
  int opcode = functionMap_.find(e_ptr->op)->second;
  switch (opcode) {
   case (OPPLUS):   // expr1 + expr2
     out << "(";
     writeExpression_(e_ptr->L.e, out);
     out << " + ";
     writeExpression_(e_ptr->R.e, out);
     out << ")";
     break;
   case (OPMINUS):  // expr1 - expr2
     out << "(";
     writeExpression_(e_ptr->L.e, out);
     out << " - ";
     writeExpression_(e_ptr->R.e, out);
     out << ")";
     break;
   case (OPMULT): // expr1*expr2
     out << "(";
     writeExpression_(e_ptr->L.e, out);
     out << " * ";
     writeExpression_(e_ptr->R.e, out);
     out << ")";
     break;
   case (OPDIV): // expr1/expr2
     out << "(";
     writeExpression_(e_ptr->L.e, out);
     out << " / ";
     writeExpression_(e_ptr->R.e, out);
     out << ")";
     break;
   case (OPREM): // remainder by dividing expr1/expr2
     writeExpression_(e_ptr->L.e, out);
     out << " % ";
     writeExpression_(e_ptr->R.e, out);
     break;
   case (OPPOW): // expr1 ^ expr2
     out << "(";
     writeExpression_(e_ptr->L.e, out);
     out << "^";
     writeExpression_(e_ptr->R.e, out);
     out << ")";
     break;
   case (OPUMINUS):
     out << "(";
     out << "-";
     writeExpression_(e_ptr->L.e, out);
     out << ")";
     break;
   case (OPSUMLIST):
     expr **ep, **epe;
     ep = e_ptr->L.ep;
     epe = e_ptr->R.ep;
     out << "(";
     while (ep < epe) {
       writeExpression_(*ep, out);
       ++ep;
       if (ep<epe) {
         out << " + ";
       }
     }
     out << ")";
     break;
   case (OP1POW): //  OPPOW for R = numeric constant
     writeExpression_(e_ptr->L.e, out);
     out << "^" << ((expr_n *)e_ptr->R.e)->v;
     break;
   case (OP2POW): //  expr^2
     writeExpression_(e_ptr->L.e, out);
     out << "^2";
     break;
   case (OPCPOW): //  (constant)^expr
     out << ((expr_n *)e_ptr->R.e)->v << "^";
     writeExpression_(e_ptr->R.e, out);
     break;
   case (OPNUM): //  numeric constant
       out << ((expr_n *)e_ptr)->v ;
       break;
   case (OPVARVAL): //  single variable
       out << std::string(var_name_ASL(myAsl_, 
             ((expr_v *)e_ptr - ((ASL_fg *)myAsl_)->I.var_e_)));
       break;
   case (OP_tanh):  
       out << "tanh(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_tan):  
       out << "tan(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_sqrt): 
       out << "sqrt(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_sinh): 
       out << "sinh(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_sin):   
       out << "sin(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_log10):
       out << "log10(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_log):  
       out << "log(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_exp):   
       out << "e^(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_cosh): 
       out << "acosh(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_cos):  
       out << "acos(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_atanh):
       out << "atanh(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_atan2): 
       out << "atan2(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_atan): 
       out << "atan(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_asinh):
       out << "asinh(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_asin): 
       out << "asin(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_acosh): 
       out << "acosh(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OP_acos):
       out << "acos(";
       writeExpression_(e_ptr->L.e, out);
       out << ") ";
       break;
   case (OPLESS):      case (MINLIST):     case (MAXLIST):     case (FLOOR):
   case (CEIL):        case (ABS):         case (OPOR):        case (OPAND): 
   case (LT):          case (LE):          case (EQ):          case (GE):
   case (NE):          case (OPNOT):       case (OPIFnl):      case (OPintDIV):
   case (OPprecision): case (OPround):     case (OPtrunc):     case (OPCOUNT):
   case (OPNUMBEROF):  case (OPNUMBEROFs): case (OPATLEAST):   case (OPATMOST):
   case (OPPLTERM):    case (OPIFSYM):     case (OPEXACTLY):   
   case (OPNOTATLEAST):                    case (OPNOTATMOST): 
   case (OPNOTEXACTLY):                    case (ANDLIST):     case (ORLIST):
   case (OPIMPELSE):   case (OP_IFF):      case (OPALLDIFF):   case (OPFUNCALL):
     out << "Operation code " << opcode
         << " Minotaur cannot handle this operation code."
         << std::endl;
     break;
  }
}


void AMPLInterface::writeLin_(Minotaur::UInt i, bool is_obj, 
    std::ostream &out) const
{
  if (is_obj) {
    ograd *og; 
    for (og = myAsl_->i.Ograd_[i]; og; og = og->next) {
      if (fabs(og->coef) > 1e-12) {
        if (og->coef>0) {
          out << " + ";
        } else {
          out << " - ";
        }
        out << fabs(og->coef) << "*" << var_name_ASL(myAsl_,og->varno);
      }
    }
    if (objconst_ASL(myAsl_, 0) > 0) {
      out << " + ";
    } else {
      out << " - ";
    }
    out << fabs(objconst_ASL(myAsl_, 0));
  } else {
    cgrad *cg; 
    for (cg = myAsl_->i.Cgrad_[i]; cg; cg = cg->next) {
      if (fabs(cg->coef) > 1e-12) {
        if (cg->coef>0) {
          out << " + ";
        } else {
          out << " - ";
        }
        out << fabs(cg->coef) << "*" << var_name_ASL(myAsl_,cg->varno);
      }
    }
  }
}


void AMPLInterface::writeProblem(std::ostream &out) const
{
  int mark=(myAsl_->i.nlvo_ > myAsl_->i.nlvc_)?myAsl_->i.nlvo_:myAsl_->i.nlvc_;
  out.precision(8);

  for (int i=0; i<myAsl_->i.n_var_; ++i) {
    out << "var " << std::string(var_name_ASL(myAsl_, i));
    if (myAsl_->i.LUv_[2*i]>-1e20) {
      out << ", >= " << myAsl_->i.LUv_[2*i];
    }
    if (myAsl_->i.LUv_[2*i+1]<1e20) {
      out << ", <= " << myAsl_->i.LUv_[2*i+1];
    }
    if (i >= myAsl_->i.nlvb_-myAsl_->i.nlvbi_ && i<myAsl_->i.nlvb_) {
      out << ", integer";
    } else if (i>=myAsl_->i.nlvc_-myAsl_->i.nlvci_ && i<myAsl_->i.nlvc_) {
      out << ", integer";
    } else if (myAsl_->i.nlvo_ > myAsl_->i.nlvc_ && 
        i>=myAsl_->i.nlvo_-myAsl_->i.nlvoi_ && i<myAsl_->i.nlvo_) {
      out << ", integer";
    } else if (i>=mark && i>=myAsl_->i.n_var_ - 
               (myAsl_->i.niv_ + myAsl_->i.nbv_)) {
      out << ", integer";
    }

    out << ";" << std::endl;
  }

  if (myAsl_->i.n_obj_ > 0) {
    if (myAsl_->i.objtype_[0] == 1) {
      out << "maximize ";
    } else {
      out << "minimize "; 
    }
    out << std::string(std::string(obj_name_ASL(myAsl_, 0))) << ": ";
    writeExpression_(0, true, out);
    writeLin_(0, true, out);
    out << ";" << std::endl;
  }

  for (int i=0; i<myAsl_->i.nlc_ - myAsl_->i.nlnc_; ++i) {
    out << "s.t. " << std::string(con_name_ASL(myAsl_, i)) << ": ";
    if (myAsl_->i.LUrhs_[2*i]>-INFINITY) {
      out << myAsl_->i.LUrhs_[2*i] << " <= ";
    }
    writeExpression_(i, false, out);
    writeLin_(i, false, out);
    if (myAsl_->i.LUrhs_[2*i+1]<INFINITY) {
      out << " <= " << myAsl_->i.LUrhs_[2*i+1];
    }
    out << ";" << std::endl;
  }
  for (int i=myAsl_->i.nlc_ - myAsl_->i.nlnc_; i<myAsl_->i.n_con_; ++i) {
    out << "s.t. " << std::string(con_name_ASL(myAsl_, i)) << ": ";
    if (myAsl_->i.LUrhs_[2*i]>-INFINITY) {
      out << myAsl_->i.LUrhs_[2*i] << " <= ";
    }
    writeLin_(i, false, out);
    if (myAsl_->i.LUrhs_[2*i+1]<INFINITY) {
      out << " <= " << myAsl_->i.LUrhs_[2*i+1];
    }
    out << ";" << std::endl;
  }
}


void AMPLInterface::writeSolution(Minotaur::ConstSolutionPtr sol, 
                                  Minotaur::SolveStatus status)
{
  Option_Info *option_info = NULL;
  double *y = NULL;
  double *x = NULL;
  char *cstr = new char[2];
  cstr[0] = ' ';
  cstr[1] = '\0';
  /*
  From README.suf in ASL code:
  Another detail not yet documented in "Hooking Your Solver..." is that
  you should assign a solve_result_num value to indicate success,
  failure, iteration limit, etc. before calling write_sol.  This is
  an integer value in one of the ranges indicated by AMPL's default
  $solve_result_table:
  
          0       solved
          100     solved?
          200     infeasible
          300     unbounded
          400     limit
          500     failure
  
  For successful solves, solve_result_num should thus be an integer
  in [0,99].  If the problem might be solved, but tolerances may have
  been too tight to satisfy all stopping tests, assign an integer in
  [100,199], etc. -- any value >= 500 indicates failure (such as
  insufficient memory).  Then AMPL's symbolic solve_result will be
  assigned one of the values in the second column of the above table,
  and scripts that need more detail can base tests on solve_result_num.
  
  Many of the sample solver interfaces appearing in subdirectories
  of netlib's ampl/solvers directory make use of suffixes and supply
  solve_result_num.
  */

  switch (status) {
  case (Minotaur::NotStarted):
    myAsl_->p.solve_code_ = 501;
    break;
  case (Minotaur::Started):
    myAsl_->p.solve_code_ = 502;
    break;
  case (Minotaur::Restarted):
    myAsl_->p.solve_code_ = 503;
    break;
  case (Minotaur::SolvedOptimal):
    myAsl_->p.solve_code_ = 0;
    break;
  case (Minotaur::SolvedInfeasible):
    myAsl_->p.solve_code_ = 200;
    break;
  case (Minotaur::SolvedUnbounded):
    myAsl_->p.solve_code_ = 300;
    break;
  case (Minotaur::SolvedGapLimit):
    myAsl_->p.solve_code_ = 401;
    break;
  case (Minotaur::SolvedSolsLimit):
    myAsl_->p.solve_code_ = 402;
    break;
  case (Minotaur::IterationLimitReached):
    myAsl_->p.solve_code_ = 403;
    break;
  case (Minotaur::Interrupted):
    myAsl_->p.solve_code_ = 501;
    break;
  case (Minotaur::TimeLimitReached):
    myAsl_->p.solve_code_ = 404;
    break;
  case (Minotaur::Finished):
    myAsl_->p.solve_code_ = 504;
    break;
  default:
    myAsl_->p.solve_code_ = 500;
    break;
  }

  if (sol) {
    const double *best_x = sol->getPrimal();
    x = new double[nVars_];

    // assume that sol has nVars_ variables.
    std::copy(best_x, best_x+nVars_, x);
    write_sol_ASL(myAsl_, cstr, x, y, option_info);
    delete [] x;
  } else {
    write_sol_ASL(myAsl_, cstr, x, y, option_info);
  }
  delete [] cstr;
}


// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
