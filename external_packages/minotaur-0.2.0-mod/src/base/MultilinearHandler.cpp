//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file MultilinearHandler.cpp
 * \brief Handles Multilinear Constraints
 * \author Mahdi Namizifar wrote the ugly parts.  Jeff Linderoth wrote the
 * good parts.
 */


#include <stack>

#include <cmath>
#include <stdio.h>
#include <iostream>
#include <cstdio>
#include "MinotaurConfig.h"
#include "Branch.h"
#include "BrVarCand.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "MultilinearHandler.h"
#include "Node.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "PolynomialFunction.h"
#include "QuadraticFunction.h"
#include "Relaxation.h"
#include "Solution.h"
#include "NLPEngine.h"
#include "Types.h"
#include "Variable.h"
#include "VarBoundMod.h"

using namespace Minotaur;
#undef DEBUG_MULTILINEAR_HANDLER
#undef DEBUG_MULTILINEAR_HANDLER2

MultilinearHandler::MultilinearHandler(EnvPtr env, ProblemPtr problem)
  : env_(env), problem_(problem)
{
  logger_  = (LoggerPtr) new Logger((LogLevel) 
                                    env->getOptions()->findInt("handler_log_level")->getValue());
  workingProblem_ = problem_->clone();

  linearizationCnt_ = 5;
  eTol_ = env->getOptions()->findDouble("ml_feastol")->getValue();

#if defined DEBUG_MULTILINEAR_HANDLER2
  std::cout << "Original problem variable pointers: " << std::endl;
  for(VariableConstIterator it = problem_->varsBegin(); it != problem_->varsEnd(); ++it) {
    ConstVariablePtr x = *it;
    std::cout << "ID: " << x->getId() << " Address: " << (  x.get()) << std::endl;
  }

  std::cout << "Working problem variable pointers: " << std::endl;
  for(VariableConstIterator it = workingProblem_->varsBegin(); it != workingProblem_->varsEnd(); ++it) {
    ConstVariablePtr x = *it;
    std::cout << "ID: " << x->getId() << " Address: " << (  x.get()) << std::endl;
  }

  std::cout << "Variable pointers in constraint functions: " << std::endl;

  for(ConstraintConstIterator it = workingProblem_->consBegin(); it != workingProblem_->consEnd(); ++it) {
    FunctionPtr constraintF = (*it)->getFunction();
    LinearFunctionPtr lcf = constraintF->getLinearFunction();
    QuadraticFunctionPtr qcf = constraintF->getQuadraticFunction();
    NonlinearFunctionPtr nlcf = constraintF->getNonlinearFunction();
    PolyFunPtr pcf = boost::dynamic_pointer_cast <PolynomialFunction> (nlcf);

    if(pcf) {
      pcf->removeLinear(lcf);
      pcf->removeQuadratic(qcf);
    }
    std::cout << "Constraint address: " << constraintF.get() << std::endl;
    if(lcf) {
      for(VariableGroupConstIterator vit = lcf->termsBegin(); vit != lcf->termsEnd(); ++vit) {
        std::cout << "  Address: " << vit->first.get() << std::endl;
      }
    }
    if (qcf) {
      for(VariablePairGroupConstIterator qit = qcf->begin(); qit != qcf->end(); ++qit) {
        ConstVariablePtr x1 = qit->first.first;
        ConstVariablePtr x2 = qit->first.second;
        std::cout << "  Address: " << x1.get() << std::endl;
        std::cout << "  Address: " << x2.get() << std::endl;
      }
    }
  }
      
#endif

}

bool MultilinearHandler::findOriginalVariable(ConstVariablePtr rv, 
                                              ConstVariablePtr & ov) const
{
  bool retval = false;
  std::map<ConstVariablePtr, ConstVariablePtr>::const_iterator it;
  it = rev_oVars_.find(rv);
  if (it != rev_oVars_.end()) {
    rv->write(std::cout);
    ov = it->second;
    ov->write(std::cout);
    retval = true;
  }
   
  return retval;
}


Branches MultilinearHandler::getBranches(BrCandPtr cand, DoubleVector & x,
                                         RelaxationPtr, SolutionPoolPtr)
{
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();

  
#if defined(DEBUG_MULTILINEAR_HANDLER)  
  std::cout << "Branching candidate (working problem) ID: " << v->getId() 
            << " address: " << (  v.get() );
  std::cout << " Relaxation ID: " << (oVars_[v])->getId() << std::endl;
  
#endif

  // x is a *relaxation* solution, while we have put the *original* (or
  // working) problem variables into the BrCandPtr, so we need to
  // update our value appropriately...
  double value = x[oVars_[v]->getId()];
  BranchPtr branch;
  Branches branches = (Branches) new BranchPtrVector();

  // can't branch on something that is at its bounds.
  if (!(value > v->getLb()+1e-8 && value < v->getUb()-1e-8)) {
    std::cerr << "Warning!  Branching on variable with bounds/value: [" << 
      v->getLb() << " , " << value << "  " << v->getUb() << " ]" << std::endl;
    //assert(value > v->getLb()+1e-8 && value < v->getUb()-1e-8);
  }

  // down branch
  VarBoundModPtr mod = (VarBoundModPtr) new VarBoundMod(v, Upper, value);
  branch = (BranchPtr) new Branch();
  branch->addMod(mod);
  branch->setActivity(0.5);// TODO: set this correctly
  branches->push_back(branch);

  // up branch
  mod    = (VarBoundModPtr) new VarBoundMod(v, Lower, value);
  branch = (BranchPtr) new Branch();
  branch->addMod(mod);
  branch->setActivity(0.5); // TODO: set this correctly
  branches->push_back(branch);

  logger_->msgStream(LogDebug2) << "branching on " << v->getName();
  logger_->msgStream(LogDebug2) << " <= " << value << " or " 
    << " >= " << value << std::endl;

#if defined(DEBUG_MULTILINEAR_HANDLER)  
  std::cout << "branching on " << v->getName();
  std::cout << " <= " << value << " or " << " >= " << value << std::endl;
#endif
  
  return branches;
}


void MultilinearHandler::getBranchingCandidates(RelaxationPtr, 
                                                const DoubleVector &x, ModVector &, 
                                                BrCandSet &cands, bool &isInf )
{

  // We use this just to check if we have already added a variable
  UIntSet cand_inds;

  std::pair<UIntSet::iterator, bool> ret;

  for(std::map<ConstVariablePtr, ConstVariablePair>::iterator iter = blterms_.begin();
      iter != blterms_.end(); ++iter) {
    ConstVariablePtr rx1 = iter->second.first;
    ConstVariablePtr rx2 = iter->second.second;
    ConstVariablePtr w = iter->first;
    double x1val = x[rx1->getId()];
    double x2val = x[rx2->getId()];
    double wval = x[w->getId()];
    
    // We put WORKING variable in as branching candidate.
    // (It work for what we want to do -- MaxViolbrancher only)
    // So we do a reverse lookup for the working variable pointer
    // here.  (I do it the 'safe' way, checking)

    std::map<ConstVariablePtr, ConstVariablePtr>::iterator wvar_it;
    wvar_it = rev_oVars_.find(rx1);
    if (wvar_it == oVars_.end()) {
      std::cerr << "WOAH.  Couldn't find working variable from relaxation variable" << std::endl;
      assert(0);
    }
    ConstVariablePtr x1 = wvar_it->second;
    wvar_it = rev_oVars_.find(rx2);
    if (wvar_it == oVars_.end()) {
      std::cerr << "WOAH.  Couldn't find working variable from relaxation variable" << std::endl;
      assert(0);
    }
    ConstVariablePtr x2 = wvar_it->second;

    double viol = fabs(wval-x1val*x2val);
    if (viol > eTol_) {
      // This may be a candidate
      ret = cand_inds.insert(x1->getId());
      if (ret.second == true) {
#if defined(DEBUG_MULTILINEAR_HANDLER)
        std::cout << "Candidate is: ";
        x1->write(std::cout);
        //std::cout << "This: " << (  x1.get()) << std::endl;
        std::cout << "Value: " << x1val << std::endl;
        std::cout << "ID: " << x1->getId() << std::endl;
#endif        
        BrVarCandPtr br_can = (BrVarCandPtr) new BrVarCand(x1, x1->getId(), viol, viol); 
        cands.insert(br_can);
      }
      ret = cand_inds.insert(x2->getId());
      if (ret.second == true) {
#if defined(DEBUG_MULTILINEAR_HANDLER)
        std::cout << "Candidate is: ";
        x2->write(std::cout);
        std::cout << "Value: " << x2val << std::endl;
        //std::cout << "This: " << (  x2.get()) << std::endl;
        std::cout << "ID: " << x2->getId() << std::endl;
#endif        
        BrVarCandPtr br_can = (BrVarCandPtr) new BrVarCand(x2, x2->getId(), viol, viol); 
        cands.insert(br_can);
      }
      
    }
  }

  //XXX: TODO Also need to do this for multilinear terms

  isInf = false;

}

ModificationPtr MultilinearHandler::getBrMod(BrCandPtr , DoubleVector &, 
                                             RelaxationPtr , BranchDirection )
{
  std::cout << "getBrMod() does nothing" << std::endl;
  return ModificationPtr();
}


bool
MultilinearHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr,
                               bool &should_prune)
{

  should_prune = false;

  const double *x = sol->getPrimal();

  bool feas = true;
  for(std::map<ConstVariablePtr, ConstVariablePair>::iterator iter = blterms_.begin();
      iter != blterms_.end(); ++iter) {
    ConstVariablePtr x1 = iter->second.first;
    ConstVariablePtr x2 = iter->second.second;
    ConstVariablePtr w = iter->first;
    double x1val = x[x1->getId()];
    double x2val = x[x2->getId()];
    double wval = x[w->getId()];

#if defined(DEBUG_MULTILINEAR_HANDLER2)
    std::cout << "x1: ";
    x1->write(std::cout);
    std::cout<< " val " << x1val << std::endl;
    std::cout << "x2: ";
    x2->write(std::cout);
    std::cout<< " val " << x2val << std::endl;
    std::cout << "w: ";
    w->write(std::cout);
    std::cout<< " val " << wval << std::endl;
#endif

    if (fabs(wval-x1val*x2val) > eTol_) {
      feas = false;
    }
  }

  //XXX Can break once feasible is false
  //XXX: TODO Also need to implement for multilinear terms


  return(feas);
}

void
MultilinearHandler::relax(RelaxationPtr relaxation, bool & should_prune)
{

  // Since this relaxation builds from scratch every time, we need to clear *everything* 
  //  out that may have held state.  This is a bad way to do it.
  clearAllContainers();

#if 0  
  // If the objective is nonlinear (Bilinear, Quadratic, or Multilinear) move it to the constraints
  // -take care of objective.
  ObjectivePtr oPtr;
  oPtr = workingProblem_->getObjective();
  FunctionPtr objF = oPtr->getFunction();
  FunctionType objFType = oPtr->getFunctionType();
#if defined(DEBUG_MULTILINEAR_HANDLER)
  std::cout << "Objective function is: " << getFunctionTypeString(objFType) << std::endl;
#endif
  objModified_ = 0;

  // --Quadratic objective
  if ((objFType != Linear) && (objFType == Quadratic || objFType == Bilinear)) {
    QuadraticFunctionPtr quadObjF = oPtr->getQuadraticFunction();

    VariablePtr z = workingProblem_->newVariable(-INFINITY, INFINITY, Continuous);
    LinearFunctionPtr lz = LinearFunctionPtr(new LinearFunction());
    lz->addTerm(z, -1.0);
    FunctionPtr objConF = FunctionPtr(new Function(lz, quadObjF));
    ConstraintPtr objCon = workingProblem_->newConstraint(objConF, 0.0, 0.0);

    LinearFunctionPtr mlz = LinearFunctionPtr(new LinearFunction());
    mlz->addTerm(z, 1.0);
    workingProblem_->removeQuadFromObj();
    workingProblem_->addToObj(mlz);
    objModified_ = 1;
  }
  
  // --Multilinear objective
  if ((oPtr->getFunctionType() != Linear && oPtr->getFunctionType() != Quadratic && 
       oPtr->getFunctionType() != Bilinear) && (oPtr->getFunctionType() == Polynomial)) {

    LinearFunctionPtr linObjF = oPtr->getLinearFunction();
    LinearFunctionPtr olf = LinearFunctionPtr(new LinearFunction());
    for(VariableGroupConstIterator vit = linObjF->termsBegin(); vit != linObjF->termsEnd(); ++vit) {
      olf->addTerm(vit->first, vit->second);
    }

    LinearFunctionPtr negLinObjF = -1*linObjF;
    objF->add(negLinObjF);

    VariablePtr z = workingProblem_->newVariable(-INFINITY, INFINITY, Continuous);
    olf->addTerm(z, +1.0);

    LinearFunctionPtr lz = LinearFunctionPtr(new LinearFunction());
    lz->addTerm(z, -1.0);
    objF->add(lz);
    ConstraintPtr objCon = workingProblem_->newConstraint(objF, 0.0, 0.0);

    workingProblem_->removeObjective();
    FunctionPtr newObjF = (FunctionPtr) new Function(olf);
    workingProblem_->newObjective(newObjF, oPtr->getConstant(), oPtr->getObjectiveType());
    objModified_ = 1;
  }

  // --Other nonlinear functions
  if (oPtr->getFunctionType() == Nonlinear)  {
    assert(!"can't handle nonlinear objective function yet!");
  }
#endif

  //STRG
  // 5: term by term
  // 4: 1 round of heuristic coverage grouping stratgy
  // 1 and 2: obsolete

  std::string s = env_->getOptions()->findString("ml_group_strategy")->getValue();
  if("ALL" == s) {
    makeGroupedConvexHull(relaxation, should_prune, 4, objModified_);
  }
  else if ("TT" == s) {
    makeGroupedConvexHull(relaxation, should_prune, 5, objModified_);
  }
  else if ("MC" == s) {
    makeMcCormick(relaxation, should_prune);
  }
  else {
    assert (!"ml_group_stretgy must be one of 'ALL', 'TT', 'MC'");
  }

}

void
MultilinearHandler::relaxNode(NodePtr node,
                              RelaxationPtr relaxation, bool & should_prune)
{


  ModificationConstIterator mod_it;
  ModificationPtr mod;

#if defined(DEBUG_MULTILINEAR_HANDLER)
  std::cout << "Relaxing node: " << std::endl;
  node->write(std::cout);
#endif

  // You need to reset the bounds on workingproblem from originalproblem
  for(VariableConstIterator it = problem_->varsBegin(); it != problem_->varsEnd(); ++it) {
    ConstVariablePtr v = *it;
    UInt id = v->getId();
    double lb = v->getLb();
    double ub = v->getUb();
    workingProblem_->changeBound(id, lb, ub);
  }
         

  // traceback to root and put in all modifications that need to go into the
  // relaxation and the engine.
  std::stack<NodePtr> predecessors;
  NodePtr t_node = node->getParent();
  
  while (t_node) {
    predecessors.push(t_node);
    t_node = t_node->getParent();
  }
  
  // Now apply them in order
   while (!predecessors.empty()) {
     t_node = predecessors.top();
     //t_node->write(std::cout);

     BranchPtr b = t_node->getBranch();
     if (b) {
       for(mod_it = b->modsBegin(); mod_it != b->modsEnd(); ++mod_it) {        
         mod = *mod_it;
         mod->applyToProblem(workingProblem_);
       }    
#if defined(DEBUG_MULTILINEAR_HANDLER)
       std::cout << "Branch: " << std::endl;
       b->write(std::cout);
#endif
     }

     // delete from stack
     predecessors.pop();
   }     

#if defined(DEBUG_MULTILINEAR_HANDLER)
   std::cout << "Node Mods: " << std::endl;
#endif

   node->applyMods(workingProblem_);
     
#if defined(DEBUG_MULTILINEAR_HANDLER)
  std::cout << "After applying mods problem is: " << std::endl;
  workingProblem_->write(std::cout);
#endif

  // Now you need to relax the problem into relaxation
  relax(relaxation, should_prune);

}

void 
MultilinearHandler::clearAllContainers()
{

#if defined(DEBUG_MULTILINEAR_HANDLER)
  std::cout << "Before, groups size: " << groups_.size() << std::endl;
#endif

  max_pow_.clear();

  oVars_.clear();
  rev_oVars_.clear();
  sqterms_.clear();
  rev_sqterms_.clear();

  blterms_cons_.clear();
  blterms_cons_coef_.clear();
  rev_blterms_cons_.clear();
  mlterms_cons_.clear();
  mlterms_cons_coef_.clear();
  rev_mlterms_cons_.clear();
  
  blterms_obj_.clear();
  blterms_obj_coef_.clear();
  rev_blterms_obj_.clear();
  mlterms_obj_.clear();
  mlterms_obj_coef_.clear();
  rev_mlterms_obj_.clear();
  
  blterms_.clear();

  blterms_coef_.clear();
  rev_blterms_.clear();
  mlterms_.clear();
  mlterms_coef_.clear();
  rev_mlterms_.clear();
  monomial_terms_.clear();
  groups_.clear();
  
  all_lambdas_.clear();
  newCopyVariables_.clear();

#if defined(DEBUG_MULTILINEAR_HANDLER)
  std::cout << "After, groups size: " << groups_.size() << std::endl;
#endif

}

void
MultilinearHandler::makeGroupedConvexHull(RelaxationPtr relaxation, bool &should_prune,
                                          int groupStrategy, bool objModified)
                                          
{

  // Loop 1: Find each (original) variable in the WORKING problem and add it to the relaxation
  for(ConstraintConstIterator it = workingProblem_->consBegin(); it != workingProblem_->consEnd(); ++it) {
    FunctionPtr constraintF = (*it)->getFunction();
    LinearFunctionPtr lcf = constraintF->getLinearFunction();
    QuadraticFunctionPtr qcf = constraintF->getQuadraticFunction();
    NonlinearFunctionPtr nlcf = constraintF->getNonlinearFunction();
    PolyFunPtr pcf = boost::dynamic_pointer_cast <PolynomialFunction> (nlcf);

    if(pcf) {
      pcf->removeLinear(lcf);
      pcf->removeQuadratic(qcf);
    }

    // In each term, if a binary variable takes a power bigger than 1 make that power 1
    makeBinaryVariablesPowers1(constraintF, pcf, qcf, lcf);


    // For each variable appearing in the constraints 
    // - have a map between the original variables and the new variables introduced for the "ralaxation"
    // - find the maximum exponent of the variable

    std::map<ConstVariablePtr, ConstVariablePtr>::iterator ovar_it;
    if(lcf) {
      for(VariableGroupConstIterator vit = lcf->termsBegin(); vit != lcf->termsEnd(); ++vit) {
        ovar_it = oVars_.find(vit->first);
        if(ovar_it == oVars_.end()) {
          VariablePtr newVar = relaxation->newVariable((vit->first)->getLb(), (vit->first)->getUb(), (vit->first)->getType());
          oVars_.insert(std::make_pair(vit->first, newVar));
          rev_oVars_.insert(std::make_pair(newVar, vit->first));
          max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
        }
      }
    }

    if(qcf) {
      for(VariablePairGroupConstIterator qit = qcf->begin(); qit != qcf->end(); ++qit) {
        ConstVariablePtr x1 = qit->first.first;
        ConstVariablePtr x2 = qit->first.second;
        ovar_it = oVars_.find(x1);
        if(ovar_it == oVars_.end()) {
          VariablePtr newVar;
          newVar = relaxation->newVariable(x1->getLb(), x1->getUb(), x1->getType());
          oVars_.insert(std::make_pair(x1, newVar));
 
          rev_oVars_.insert(std::make_pair(newVar, x1));
          max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
        }
        ovar_it = oVars_.find(x2);
        if(ovar_it == oVars_.end()) {
          VariablePtr newVar;
          newVar = relaxation->newVariable(x2->getLb(), x2->getUb(), x2->getType());
          oVars_.insert(std::make_pair(x2, newVar));

          rev_oVars_.insert(std::make_pair(newVar, x2));
          max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
        }
      }
    }

    if(pcf) {
      for(MonomialConstIter mit = pcf->termsBegin(); mit != pcf->termsEnd(); ++mit) {
        for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
          ovar_it = oVars_.find(tit->first);
          if(ovar_it == oVars_.end()) {
            VariablePtr newVar;
            newVar = relaxation->newVariable((tit->first)->getLb(), (tit->first)->getUb(), 
                                             (tit->first)->getType());
            oVars_.insert(std::make_pair(tit->first, newVar));
            rev_oVars_.insert(std::make_pair(newVar, tit->first));
            max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
          }
        }
      }
    } 
  }
#if defined(DEBUG_MULTILINEAR_HANDLER2)
  std::cout << "Original variables in problem: " << std::endl;
  for(std::map <ConstVariablePtr, ConstVariablePtr>::iterator it = oVars_.begin();
      it != oVars_.end(); ++it) {
    ConstVariablePtr v = it->first;
    std::cout << "Var ID: " << v->getId() << " ";
    v->write(std::cout);
    v = it->second;
    std::cout << "is mapped to relaxation var ID: " << v->getId() << " ";
    v->write(std::cout);
  }
#endif  
  std::vector<int> mlcid;
  int cix = 0;

  // Go through the multilinear constraints
  //  - for each variable find the highest power
  //  - for each monomial if it hasn't been seen before, add it to the list of monomials (monomial_terms_)
  for(ConstraintConstIterator it = workingProblem_->consBegin(); it != workingProblem_->consEnd(); ++it) {
    FunctionPtr constraintF = (*it)->getFunction();
    FunctionType ft = (*it)->getFunctionType();

    // Check the type of the function
    // - Linear: add the constraint to the relaxation
    // - Bilinear, Quadratic, and Polynomial: mark the constraint
    if (ft == Quadratic || ft == Bilinear || ft == Polynomial) {
      mlcid.push_back((*it)->getId());
      //c_list[cix] = true;
    }
    
    if (ft == Linear) {
      LinearFunctionPtr lf = LinearFunctionPtr(new LinearFunction());
      const LinearFunctionPtr olf = (*it)->getLinearFunction();
      for(VariableGroupConstIterator vit = olf->termsBegin(); vit != olf->termsEnd(); ++vit) {
        ConstVariablePtr tempVar = vit->first;
        lf->addTerm(oVars_.find(vit->first)->second, vit->second);
      }

      FunctionPtr originalLinearF = FunctionPtr(new Function(lf));
      ConstraintPtr originalLinearC = relaxation->newConstraint(originalLinearF, (*it)->getLb(), 
                                                                (*it)->getUb());
    }

    cix++;
  
    // If the function is Quadratic
    // - Find the maximum exponent of each variable
    if(ft == Quadratic) {
      LinearFunctionPtr lcf;
      QuadraticFunctionPtr qcf;
      lcf = constraintF->getLinearFunction();
      qcf = constraintF->getQuadraticFunction();
      if(qcf) {
        for(VariablePairGroupConstIterator qit = qcf->begin(); qit != qcf->end(); ++qit) {
          ConstVariablePtr x1 = qit->first.first;
          ConstVariablePtr x2 = qit->first.second;

          //Update max_pow_
          max_pow_it_ = max_pow_.find(oVars_.find(x1)->second); 
          if(x1->getId()==x2->getId()) {
            if(max_pow_it_->second < 2 ) {
              max_pow_it_->second = 2;
            }
          }
        }
      }
    }

    // If the function is Polynomial
    // - Find the maximum exponent of each variable
    // - Update the list of monomial terms
       
    if(ft == Polynomial) {
      LinearFunctionPtr lcf;
      QuadraticFunctionPtr qcf;
      NonlinearFunctionPtr nlcf = constraintF->getNonlinearFunction();
      PolyFunPtr pcf = boost::dynamic_pointer_cast <PolynomialFunction> (nlcf);
      lcf = constraintF->getLinearFunction();
      qcf = constraintF->getQuadraticFunction();
      
      // get the all quadratic part (no linear left) of the constraint function
      // update the max_pow_ map for the quadratic part
      if(qcf) {
        for(VariablePairGroupConstIterator qit = qcf->begin(); qit != qcf->end(); ++qit) {
          ConstVariablePtr x1 = qit->first.first;
          ConstVariablePtr x2 = qit->first.second;
          max_pow_it_ = max_pow_.find(oVars_.find(x1)->second); 
          if(x1->getId()==x2->getId()) {
            if(max_pow_it_->second < 2 ) {
              max_pow_it_->second = 2;
            }
          }
        }
      }
      
      // For the polynomial part update the max_pow_ and make a list of monomial terms
      if(pcf) {
        for(MonomialConstIter mit = pcf->termsBegin(); mit != pcf->termsEnd(); ++mit) {
          const VarIntMap *mit_terms = (*mit)->getTerms();
          if(monomial_terms_.find(*mit_terms) == monomial_terms_.end()) {
            // Find the bounds on the monomial
            double lb1 = 1;
            double ub1 = 1;
            double lb = 0;
            double ub = 0;
            for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
              for(UInt i = 1; i <= tit->second; i++) {
                BoundsOnProduct((tit->first)->getLb(), (tit->first)->getUb(), lb1, ub1, lb, ub);
                lb1 = lb;
                ub1 = ub;
              }
            }

            ConstVariablePtr y = relaxation->newVariable(lb, ub, Continuous);
            monomial_terms_.insert(std::pair<VarIntMap, ConstVariablePtr>(*mit_terms, y));
          }
          for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
            max_pow_it_ = max_pow_.find(oVars_.find(tit->first)->second);
            if(max_pow_it_->second < tit->second && tit->second >= 2) { 
              max_pow_it_->second = tit->second;
            }
          }
        }
      }
    }
  }

  // For each variable, based on its maximum power, build the right number of new variables
  // - For instance if max_pow_ for x1 is 3
  //   make a vector containing x1_1 and x1_2
  std::vector<ConstVariablePtr> relaxation_vars;
  for(VariableConstIterator vit = relaxation->varsBegin(); vit != relaxation->varsEnd(); ++vit) {
    relaxation_vars.push_back(*vit);
  }

  for(std::vector<ConstVariablePtr>::iterator vit = relaxation_vars.begin(); vit != relaxation_vars.end(); ++vit) {
    std::vector<ConstVariablePtr> copyV;
    max_pow_it_ = max_pow_.find(*vit);
    if(max_pow_it_ != max_pow_.end()) {
      ConstVariablePtr v = max_pow_it_->first;
      int vMaxPow = max_pow_it_->second;
      if(vMaxPow > 1) {
        for(int i=2; i<=vMaxPow; i++) {
          ConstVariablePtr newV = relaxation->newVariable(v->getLb(), v->getUb(), Continuous);
          copyV.push_back(newV);

          // Add a constraint that says v = newV
          LinearFunctionPtr slf = LinearFunctionPtr(new LinearFunction());
          slf->addTerm(v,-1.0);
          slf->addTerm(newV, 1.0);
          FunctionPtr sof = (FunctionPtr) new Function(slf);
          ConstraintPtr soc = relaxation->newConstraint(sof, 0.0, 0.0);
        }
      }
      newCopyVariables_.insert(std::pair<ConstVariablePtr, std::vector<ConstVariablePtr> > (v, copyV));
    }
  }
#if defined(DEBUG_MULTILINEAR_HANDLER)
  std::cout << "Original variables in problem: " << std::endl;
  for(std::map <ConstVariablePtr, ConstVariablePtr>::iterator it = oVars_.begin();
      it != oVars_.end(); ++it) {
    ConstVariablePtr v = it->first;
    std::cout << "Var ID: " << v->getId() << " ";
    v->write(std::cout);
    v = it->second;
    std::cout << "is mapped to relaxation var ID: " << v->getId() << " ";
    v->write(std::cout);
  }
#endif  
  // Get the objective function of the problem and add it to the relaxation
  LinearFunctionPtr rlf = LinearFunctionPtr(new LinearFunction());
  const ObjectivePtr obj = workingProblem_->getObjective();
  const LinearFunctionPtr objlf = obj->getLinearFunction();
  if (objlf){
    for(VariableGroupConstIterator it = objlf->termsBegin(); it != objlf->termsEnd();++it) {
      ConstVariablePtr tempVar;
      tempVar = it->first;
#if defined(DEBUG_MULTILINEAR_HANDLER)
      std::cout << "Orig Var ID:" << tempVar->getId() << " ";
      tempVar->write(std::cout);
      
      std::map <ConstVariablePtr, ConstVariablePtr>::iterator f;
      f = oVars_.find(tempVar);
      if (f == oVars_.end()) {
        std::cout << "NOT FOUND" << std::endl;
        assert(0);
      }
      else {
        ConstVariablePtr t2 = f->second;      
        std::cout << "Found in map: ";
        t2->write(std::cout);
      }
#endif
      rlf->addTerm(oVars_.find(tempVar)->second, it->second);
    }
  }
  FunctionPtr of = FunctionPtr(new Function(rlf));
  ObjectivePtr oo = relaxation->newObjective(of, 0, Minimize);

  // Go through multilinear rows.  Add constraints, and create maps for product vars
  //    We want to have the terms of the constraints and objective separately
  //  - If objective is modified, original constraints are all the mlcid except for the last one
  UInt mlConsCnt;
  if(objModified)
    mlConsCnt = mlcid.size()-1;
  else
    mlConsCnt = mlcid.size();

  // If the objective was modified, its function now is the last element of mlcid
  // For this function, form the maps and replace the multilinear terms
  if(objModified) {
    const ConstraintPtr objmlc = workingProblem_->getConstraint(mlcid[mlcid.size()-1]);
    const LinearFunctionPtr objlf = objmlc->getLinearFunction();
    const QuadraticFunctionPtr objqf = objmlc->getQuadraticFunction();
    const NonlinearFunctionPtr objnlf = objmlc->getNonlinearFunction();
    PolyFunPtr objpf = boost::dynamic_pointer_cast <PolynomialFunction> (objnlf);

    LinearFunctionPtr lf = LinearFunctionPtr(new LinearFunction());

    // Linear part of constraint remains the same
    if (objlf){
      for(VariableGroupConstIterator it = objlf->termsBegin(); it != objlf->termsEnd();++it) {
        ConstVariablePtr tempVar;
        tempVar = it->first;
        lf->addTerm(oVars_.find(tempVar)->second, it->second);
      }
    }

    // Quadratic part gets a new variable for every term
    if(objqf) {
      for(VariablePairGroupConstIterator it = objqf->begin(); it != objqf->end(); ++it) {

        ConstVariablePtr x1 = oVars_.find(it->first.first)->second;
        ConstVariablePtr x2 = oVars_.find(it->first.second)->second;

        ConstVariablePtr x_1;

        // Check to see if 'it' is a square term
        if(x1->getId() == x2->getId()) {
          newCopyVariables_it_ = newCopyVariables_.find(x1);
          x_1 = ((newCopyVariables_it_)->second)[0];

          // look to see if x_1^2 already exists...
          std::map <ConstVariablePtr, ConstVariablePair>::iterator pos ;
          ConstVariablePtr y;
          ConstVariablePair Z;
          if (sqterms_.find(x1) == sqterms_.end()) {
                        
            // Make a copy of the variable
            double lb = 0.0;
            double ub = 0.0;
            BoundsOnProduct(x1,x_1,lb,ub);
            y = relaxation->newVariable(lb, ub, Continuous);
            /*
            blterms_.insert(make_pair(y,ConstVariablePair(x1,x_1)));
            if(it->second >=0)
              blterms_sign_.insert(std::pair<ConstVariablePtr, bool> (y, 1));
            else
              blterms_sign_.insert(std::pair<ConstVariablePtr, bool> (y, 0));
            rev_blterms_.insert(make_pair(ConstVariablePair(x1,x_1), y));
            */
            sqterms_.insert(make_pair(x1,ConstVariablePair(y,x_1)));
            rev_sqterms_.insert(make_pair(ConstVariablePair(y,x_1), x1));
          }
          else {
            pos = sqterms_.find(x1);
            y = pos->second.first;
          }
          
          lf->addTerm(y, it->second);

          /*
          // Add this term to the list of square terms
          sqterms_.insert(make_pair(y,ConstVariablePair(x1,x_1)));
          rev_sqterms_.insert(make_pair(ConstVariablePair(x1,x_1), y));
          */
        }
        else {
          //Bounds on product depend on whether variable bounds are < 0, > 0
          double lb = 0.0;
          double ub = 0.0;
          BoundsOnProduct(x1,x2,lb,ub);
        
          // look to see if w var already exists...
          std::map <ConstVariablePair, ConstVariablePtr>::iterator pos ;
          ConstVariablePtr w;
          if (rev_blterms_obj_.find(ConstVariablePair(x1,x2)) == rev_blterms_obj_.end()) {
            w = relaxation->newVariable(lb, ub, Continuous);
            
            std::vector<double> temp;
            temp.push_back(it->second);

            blterms_obj_.insert(make_pair(w,ConstVariablePair(x1,x2)));
            blterms_obj_coef_.insert(make_pair(ConstVariablePair(x1,x2), temp));
            rev_blterms_obj_.insert(make_pair(ConstVariablePair(x1,x2), w));

            blterms_.insert(make_pair(w,ConstVariablePair(x1,x2)));
            blterms_coef_.insert(make_pair(ConstVariablePair(x1,x2), temp));
            rev_blterms_.insert(make_pair(ConstVariablePair(x1,x2), w));
          }
          else {
            pos = rev_blterms_obj_.find(ConstVariablePair(x1,x2));
            w = pos->second;
            std::map <ConstVariablePair, std::vector<double> >::iterator poscoef;
            poscoef = blterms_obj_coef_.find(ConstVariablePair(x1,x2));
            (poscoef->second).push_back(it->second);
            poscoef = blterms_coef_.find(ConstVariablePair(x1,x2));
            (poscoef->second).push_back(it->second);
          }
          
          lf->addTerm(w,it->second);
        }
      }
    }

    // Take care of the polynomial parts
    // - Go through the polynomial parts of the constraint
    //   form the linearizations
    
    if(objpf) {
      for(MonomialConstIter mit = objpf->termsBegin(); mit != objpf->termsEnd(); ++mit) {
        
        const VarIntMap* mit_terms = (*mit)->getTerms();
        double mit_coeff = (*mit)->getCoeff();
        monomial_terms_it_ = monomial_terms_.find(*mit_terms);
        if(monomial_terms_it_ != monomial_terms_.end()) {
          std::vector<ConstVariablePtr> mit_ml;
          for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
            ConstVariablePtr relax_tit_var = (oVars_.find(tit->first))->second;
            if(tit->second == 1) {
              mit_ml.push_back(relax_tit_var);
            }
            if(tit->second >= 2) {
              mit_ml.push_back(relax_tit_var);
              newCopyVariables_it_ = newCopyVariables_.find(relax_tit_var);
              for(UInt i=2; i <= tit->second; i++) {
                mit_ml.push_back((newCopyVariables_it_->second)[i-2]);
              }
            }
          }

          // See if the multilinear term 'mit_ml' is already in the list of multilinear terms
          std::map <std::vector<ConstVariablePtr>, ConstVariablePtr >::iterator rmcIt = rev_mlterms_obj_.find(mit_ml);
          if(rmcIt == rev_mlterms_obj_.end()) {
            std::vector<double> temp;
            temp.push_back(mit_coeff);
            mlterms_obj_.insert(std::pair<ConstVariablePtr, std::vector<ConstVariablePtr> > (monomial_terms_it_->second, mit_ml));
            rev_mlterms_obj_.insert(std::pair<std::vector<ConstVariablePtr>, ConstVariablePtr > (mit_ml, monomial_terms_it_->second));
            mlterms_obj_coef_.insert(std::pair<std::vector<ConstVariablePtr>, std::vector<double> > (mit_ml, temp));

            mlterms_.insert(std::pair<ConstVariablePtr, std::vector<ConstVariablePtr> > (monomial_terms_it_->second, mit_ml));
            rev_mlterms_.insert(std::pair<std::vector<ConstVariablePtr>, ConstVariablePtr > (mit_ml, monomial_terms_it_->second));
            mlterms_coef_.insert(std::pair<std::vector<ConstVariablePtr>, std::vector<double> > (mit_ml, temp));

          }
          else {
            std::map <std::vector<ConstVariablePtr>, std::vector<double> >::iterator mccIt = mlterms_obj_coef_.find(mit_ml);
            (mccIt->second).push_back(mit_coeff);
          }
        }

        bool termAdded = 0;
        for(VariableGroupConstIterator lf_it = lf->termsBegin(); lf_it != lf->termsEnd(); ++lf_it) {
          if(lf_it->first == monomial_terms_it_->second) {
            lf->incTerm(monomial_terms_it_->second, mit_coeff);
            termAdded = 1;
            break;
          }
        }
        if(!termAdded)
          lf->addTerm(monomial_terms_it_->second, mit_coeff);
      }
    }

    FunctionPtr of = (FunctionPtr) new Function(lf);
    ConstraintPtr oc = relaxation->newConstraint(of, objmlc->getLb(), objmlc->getUb());
  }

  // For the multilinear constraints, form the maps and replace the polynomial terms with the new variables
  for(UInt i = 0; i < mlConsCnt; i++) {
    const ConstraintPtr omlc = workingProblem_->getConstraint(mlcid[i]);
    const LinearFunctionPtr olf = omlc->getLinearFunction();
    const QuadraticFunctionPtr oqf = omlc->getQuadraticFunction();
    const NonlinearFunctionPtr onlf = omlc->getNonlinearFunction();
    PolyFunPtr opf = boost::dynamic_pointer_cast <PolynomialFunction> (onlf);

    LinearFunctionPtr lf = LinearFunctionPtr(new LinearFunction());
    // Linear part of constraint remains the same
    if (olf){
      for(VariableGroupConstIterator it = olf->termsBegin(); it != olf->termsEnd();++it) {
        ConstVariablePtr tempVar;
        tempVar = it->first;
        lf->addTerm(oVars_.find(tempVar)->second, it->second);
      }
    }
  
    // Quadratic part gets a new variable for every term
    if(oqf) {
      for(VariablePairGroupConstIterator it = oqf->begin(); it != oqf->end(); ++it) {
        ConstVariablePtr x1 = oVars_.find(it->first.first)->second;
        ConstVariablePtr x2 = oVars_.find(it->first.second)->second;
        //ConstVariablePtr x1 = oVars_.find(it->first.first)->second;
        //ConstVariablePtr x2 = oVars_.find(it->first.second)->second;
        ConstVariablePtr x_1;

        // Check to see if 'it' is a square term
        if(x1->getId() == x2->getId()) {
          newCopyVariables_it_ = newCopyVariables_.find(x1);
          x_1 = ((newCopyVariables_it_)->second)[0];

          // look to see if x_1^2 already exists...
          std::map <ConstVariablePtr, ConstVariablePair>::iterator pos ;
          ConstVariablePtr y;
          ConstVariablePair Z;
          if (sqterms_.find(x1) == sqterms_.end()) {
                        
            // Make a copy of the variable
            double lb = 0.0;
            double ub = 0.0;
            BoundsOnProduct(x1,x_1,lb,ub);
            y = relaxation->newVariable(lb, ub, Continuous);
            /*
            blterms_.insert(make_pair(y,ConstVariablePair(x1,x_1)));
            if(it->second >=0)
              blterms_sign_.insert(std::pair<ConstVariablePtr, bool> (y, 1));
            else
              blterms_sign_.insert(std::pair<ConstVariablePtr, bool> (y, 0));
            rev_blterms_.insert(make_pair(ConstVariablePair(x1,x_1), y));
            */
            sqterms_.insert(make_pair(x1,ConstVariablePair(y,x_1)));
            rev_sqterms_.insert(make_pair(ConstVariablePair(y,x_1), x1));
          }
          else {
            pos = sqterms_.find(x1);
            y = pos->second.first;
          }
          
          lf->addTerm(y, it->second);

          /*
          // Add this term to the list of square terms
          sqterms_.insert(make_pair(y,ConstVariablePair(x1,x_1)));
          rev_sqterms_.insert(make_pair(ConstVariablePair(x1,x_1), y));
          */
        }
        else {
          //Bounds on product depend on whether variable bounds are < 0, > 0
          double lb = 0.0;
          double ub = 0.0;
          BoundsOnProduct(x1,x2,lb,ub);
        
          // look to see if w var already exists...
          std::map <ConstVariablePair, ConstVariablePtr>::iterator pos ;
          ConstVariablePtr w;
          if (rev_blterms_cons_.find(ConstVariablePair(x1,x2)) == rev_blterms_cons_.end()) {
            w = relaxation->newVariable(lb, ub, Continuous);
            
            std::vector<double> temp;
            temp.push_back(it->second);

            blterms_cons_.insert(make_pair(w,ConstVariablePair(x1,x2)));
            blterms_cons_coef_.insert(make_pair(ConstVariablePair(x1,x2), temp));
            rev_blterms_cons_.insert(make_pair(ConstVariablePair(x1,x2), w));


          }
          else {
            pos = rev_blterms_cons_.find(ConstVariablePair(x1,x2));
            w = pos->second;
            std::map <ConstVariablePair, std::vector<double> >::iterator poscoef;
            poscoef = blterms_cons_coef_.find(ConstVariablePair(x1,x2));
            (poscoef->second).push_back(it->second);
          }

          // add the term to the map of all bilinear terms
          if(rev_blterms_.find(ConstVariablePair(x1,x2)) == rev_blterms_.end()) {
            std::vector<double> temp;
            temp.push_back(0.0); // the coefficient of this term in the objective is 0
            temp.push_back(it->second);
            
            blterms_.insert(make_pair(w,ConstVariablePair(x1,x2)));
            blterms_coef_.insert(make_pair(ConstVariablePair(x1,x2), temp));
            rev_blterms_.insert(make_pair(ConstVariablePair(x1,x2), w));
          }
          else {
            std::map <ConstVariablePair, std::vector<double> >::iterator bcIt;
            bcIt = blterms_coef_.find(ConstVariablePair(x1,x2));
            (bcIt->second).push_back(it->second);
          }
          
          lf->addTerm(w,it->second);
        }
      }
    }

    // Take care of the polynomial parts
    // - Go through the polynomial parts of the constraint
    //   form the linearizations
    
    if(opf) {
      for(MonomialConstIter mit = opf->termsBegin(); mit != opf->termsEnd(); ++mit) {
        
        const VarIntMap* mit_terms = (*mit)->getTerms();
        double mit_coeff = (*mit)->getCoeff();
        monomial_terms_it_ = monomial_terms_.find(*mit_terms);
        if(monomial_terms_it_ != monomial_terms_.end()) {
          std::vector<ConstVariablePtr> mit_ml;
          for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
            ConstVariablePtr relax_tit_var = (oVars_.find(tit->first))->second;
            if(tit->second == 1) {
              mit_ml.push_back(relax_tit_var);
            }
            if(tit->second >= 2) {
              mit_ml.push_back(relax_tit_var);
              newCopyVariables_it_ = newCopyVariables_.find(relax_tit_var);
              for(UInt i=2; i <= tit->second; i++) {
                mit_ml.push_back((newCopyVariables_it_->second)[i-2]);
              }
            }
          }

          // See if the multilinear term 'mit_ml' is already in the list of multilinear terms
          std::map <std::vector<ConstVariablePtr>, ConstVariablePtr >::iterator rmcIt = rev_mlterms_cons_.find(mit_ml);
          if(rmcIt == rev_mlterms_cons_.end()) {
            std::vector<double> temp;
            temp.push_back(mit_coeff);
            mlterms_cons_.insert(std::pair<ConstVariablePtr, std::vector<ConstVariablePtr> > (monomial_terms_it_->second, mit_ml));
            rev_mlterms_cons_.insert(std::pair<std::vector<ConstVariablePtr>, ConstVariablePtr > (mit_ml, monomial_terms_it_->second));
            mlterms_cons_coef_.insert(std::pair<std::vector<ConstVariablePtr>, std::vector<double> > (mit_ml, temp));
          }
          else {
            std::map <std::vector<ConstVariablePtr>, std::vector<double> >::iterator mccIt = mlterms_cons_coef_.find(mit_ml);
            (mccIt->second).push_back(mit_coeff);
          }

          // add the term to the map of all multilinear terms
          if(rev_mlterms_.find(mit_ml) == rev_mlterms_.end()) {
            std::vector<double> temp;
            temp.push_back(0.0); // the coefficient of this term in the objective is 0
            temp.push_back(mit_coeff);
            
            mlterms_.insert(std::pair<ConstVariablePtr, std::vector<ConstVariablePtr> > (monomial_terms_it_->second, mit_ml));
            rev_mlterms_.insert(std::pair<std::vector<ConstVariablePtr>, ConstVariablePtr > (mit_ml, monomial_terms_it_->second));
            mlterms_coef_.insert(std::pair<std::vector<ConstVariablePtr>, std::vector<double> > (mit_ml, temp));
          }
          else {
            std::map <std::vector<ConstVariablePtr>, std::vector<double> >::iterator mcIt;
            mcIt = mlterms_coef_.find(mit_ml);
            (mcIt->second).push_back(mit_coeff);
          }
        }

        bool termAdded = 0;
        for(VariableGroupConstIterator lf_it = lf->termsBegin(); lf_it != lf->termsEnd(); ++lf_it) {
          if(lf_it->first == monomial_terms_it_->second) {
            lf->incTerm(monomial_terms_it_->second, mit_coeff);
            termAdded = 1;
            break;
          }
        }
        if(!termAdded)
          lf->addTerm(monomial_terms_it_->second, mit_coeff);
      }
    }

    FunctionPtr of = (FunctionPtr) new Function(lf);
    ConstraintPtr oc = relaxation->newConstraint(of,omlc->getLb(), omlc->getUb());
  }

  // Add linearizations for square terms
  for(std::map <ConstVariablePtr, ConstVariablePair>::iterator sqiter = sqterms_.begin();
      sqiter != sqterms_.end(); ++sqiter) {
    ConstVariablePtr sqVar = sqiter->second.first;
    double sqVarUb = sqVar->getUb();
    double sqVarLb = sqVar->getLb();
    int numLin;
    if(sqVarUb - sqVarLb < 1+1e-6)
      numLin = linearizationCnt_;
    else
      numLin = floor((sqVarUb-sqVarLb)*linearizationCnt_) + 1;
    
    
    for(int i=0; i<linearizationCnt_; i++) {
      ConstVariablePtr y;
      ConstVariablePtr x1;
      x1 = sqiter->first;
      y = sqiter->second.first;
      LinearFunctionPtr sCutlf = LinearFunctionPtr(new LinearFunction());
      sCutlf = LinearFunctionPtr(new LinearFunction());
      double x1Val = x1->getLb()+i*(x1->getUb()-x1->getLb())/(linearizationCnt_-1);
      sCutlf->addTerm(x1, 2*x1Val);
      sCutlf->addTerm(y, -1);
      
      FunctionPtr sCutf = (FunctionPtr) new Function(sCutlf);
      ConstraintPtr sCutc = relaxation->newConstraint(sCutf, -INFINITY,
                                                      x1Val*x1Val );
    }

    // upper
    ConstVariablePtr y;
    ConstVariablePtr x1;
    x1 = sqiter->first;
    y = sqiter->second.first;
    LinearFunctionPtr sCutlf = LinearFunctionPtr(new LinearFunction());
    sCutlf->addTerm(x1, x1->getLb() + x1->getUb());
    sCutlf->addTerm(y, -1);
    FunctionPtr sCutf = (FunctionPtr) new Function(sCutlf);
    ConstraintPtr sCutc = relaxation->newConstraint(sCutf, x1->getLb() * x1->getUb(), INFINITY);
  }
  

  // Now call the module that groups the variables
  //  std::vector<std::vector<ConstVariablePtr> > groups;
  //  std::vector<std::vector<ConstVariablePtr> > all_lambdas;

#if defined(DEBUG_MULTILINEAR_HANDLER2)
  std::cout << "Before calling makeGroups.  Pointers of variables in blterms: " << std::endl;
  UInt jjj = 0;
  for(std::map<ConstVariablePtr, ConstVariablePair>::iterator iter = blterms_.begin();
      iter != blterms_.end(); ++iter) {
    std::cout << "Term: " << jjj++ << std::endl;
    ConstVariablePtr x1 = iter->second.first;
    ConstVariablePtr x2 = iter->second.second;
    x1->write(std::cout);
    std::cout << "Address: " << (  x1.get() ) << std::endl;
    x2->write(std::cout);
    std::cout << "Address: " << (  x2.get() ) << std::endl;
  }
#endif

  // Now we just make the groups...
  makeGroups(blterms_, 
             rev_blterms_,
             blterms_coef_,
             mlterms_, 
             rev_mlterms_,
             mlterms_coef_,
             blterms_obj_, 
             rev_blterms_obj_,
             blterms_obj_coef_,
             mlterms_obj_, 
             rev_mlterms_obj_,
             mlterms_obj_coef_,
             blterms_cons_, 
             rev_blterms_cons_,
             blterms_cons_coef_,
             mlterms_cons_, 
             rev_mlterms_cons_,
             mlterms_cons_coef_,
             groups_,
             groupStrategy);

  // Now we go through the groups and for each group we find the extreme points
  for(UInt it=0; it < groups_.size(); it++) {
    std::vector<double> lb;
    std::vector<double> ub;
    for (std::vector<ConstVariablePtr>::iterator it1 = groups_[it].begin(); it1 != groups_[it].end(); ++it1) {
      ConstVariablePtr v = *it1;
      lb.push_back(v->getLb());
      ub.push_back(v->getUb());
    }

    int numVars = groups_[it].size();
    std::vector<int> S(numVars);
    for(int i = 0; i < numVars; i++)
      S[i] = i;
    std::vector<std::vector<double> > V;
    allExtreme(S, lb, ub, V);

    // Now add the lambda variables
    std::vector<VariablePtr> lambdavars;
    for(UInt i = 0; i < V.size(); i++) {
      VariablePtr lambda = relaxation->newVariable(0.0, 1.0, Continuous);
      lambdavars.push_back(lambda);
    }
    all_lambdas_.push_back(lambdavars);

    // Add the linearization for bilinear terms
    // - Add the product constraints
    for(std::map<ConstVariablePtr, ConstVariablePair>::iterator iter = blterms_.begin();
        iter != blterms_.end(); ++iter) {
      ConstVariablePtr x1 = iter->second.first;
      ConstVariablePtr x2 = iter->second.second;
      bool x1InG = 0;
      bool x2InG = 0;
      int x1Pos = 0;
      int x2Pos = 0;
      int cntr = 0;
      for(std::vector<ConstVariablePtr>::iterator iter1 = groups_[it].begin();
          iter1 != groups_[it].end(); ++iter1) {
        if(x1->getId() == (*iter1)->getId()) {
          x1InG = 1;
          x1Pos = cntr;
        }
        if(x2->getId() == (*iter1)->getId()) {
          x2InG = 1;
          x2Pos = cntr;
        }
        ++cntr;
      }
      if(x1InG && x2InG) {
        LinearFunctionPtr lf = LinearFunctionPtr(new LinearFunction());
        lf->addTerm(iter->first, -1.0);
        for(UInt i=0; i<V.size(); i++) {
          lf->addTerm(lambdavars[i], V[i][x1Pos]*V[i][x2Pos]);
        }
        FunctionPtr f = (FunctionPtr) new Function(lf);
        ConstraintPtr c = relaxation -> newConstraint(f, 0.0, 0.0);
      }
    }

    // Add the linearization for multilinear terms
    // - Add the product constraints
    for(std::map <ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator iter = mlterms_.begin();
        iter != mlterms_.end(); ++iter) {

      bool term_in_group = 1;
      std::vector<int> pos;
      UInt term_member_cntr = 0;
      // Check to see if the multilinear term is in the group
      // - Iterate over variables of the term
      for(std::vector<ConstVariablePtr>::iterator term_it = iter->second.begin(); term_it != iter->second.end(); ++term_it) {
        term_member_cntr++;
        // - Iterate over members of the group
        int cntr = 0;
        for(std::vector<ConstVariablePtr>::iterator group_member_it = groups_[it].begin();
            group_member_it != groups_[it].end(); ++group_member_it) {
          if((*term_it)->getId() == (*group_member_it)->getId()) {
            pos.push_back(cntr);
          }
          cntr++;
        }
        if(pos.size() < term_member_cntr) {
          term_in_group = 0;
          break;
        }
      }

      if(term_in_group) {
        LinearFunctionPtr lf = LinearFunctionPtr(new LinearFunction());
        lf->addTerm(iter->first, -1.0);
        for(UInt i=0; i<V.size(); i++) {
          double lf_term_coeff = 1;
          for(UInt j=0; j<pos.size(); j++) {
            lf_term_coeff *= V[i][pos[j]];
          }
          lf->addTerm(lambdavars[i], lf_term_coeff);
        }
        FunctionPtr f = (FunctionPtr) new Function(lf);
        ConstraintPtr c = relaxation -> newConstraint(f, 0.0, 0.0);
      }
    }

    // - Add the convexity constraints
    LinearFunctionPtr convex_lf = LinearFunctionPtr(new LinearFunction());
    for(UInt i=0; i<V.size(); i++)
      convex_lf->addTerm(lambdavars[i], 1.0);
    FunctionPtr convex_f = (FunctionPtr) new Function(convex_lf);
    ConstraintPtr convex_c = relaxation->newConstraint(convex_f, 1.0, 1.0);
    
    // - Add the constraint x_i = \sum \lambda_k c_k
    int cntr = 0;
    for(std::vector<ConstVariablePtr>::iterator iter1 = groups_[it].begin();
        iter1 != groups_[it].end(); ++iter1) {
      LinearFunctionPtr sum_lf = LinearFunctionPtr(new LinearFunction());
      for(UInt i=0; i<V.size(); i++) {
        sum_lf->addTerm(lambdavars[i], V[i][cntr]);
      }
      sum_lf->addTerm(*iter1, -1.0);
      FunctionPtr sum_f = (FunctionPtr) new Function(sum_lf);
      ConstraintPtr sum_c = relaxation->newConstraint(sum_f, 0.0, 0.0);
      cntr++;
    }
  }
  should_prune = false;
  
}



void
MultilinearHandler::makeMcCormick(RelaxationPtr relaxation, bool &should_prune)
{
  // New try to initialize oVars_
  for(ConstraintConstIterator it = workingProblem_->consBegin(); it != workingProblem_->consEnd(); ++it) {
    FunctionPtr constraintF = (*it)->getFunction();
    LinearFunctionPtr lcf = constraintF->getLinearFunction();
    QuadraticFunctionPtr qcf = constraintF->getQuadraticFunction();
    NonlinearFunctionPtr nlcf = constraintF->getNonlinearFunction();
    PolyFunPtr pcf = boost::dynamic_pointer_cast <PolynomialFunction> (nlcf);

    if(pcf) {
      pcf->removeLinear(lcf);
      pcf->removeQuadratic(qcf);
    }

    // In each term, if a binary variable takes a power bigger than 1 make that power 1
    //    makeBinaryVariablesPowers1(constraintF, pcf, qcf, lcf);

    
    std::map<ConstVariablePtr, ConstVariablePtr>::iterator ovar_it;
    if(lcf) {
      for(VariableGroupConstIterator vit = lcf->termsBegin(); vit != lcf->termsEnd(); ++vit) {
        ovar_it = oVars_.find(vit->first);
        if(ovar_it == oVars_.end()) {
          VariablePtr newVar = relaxation->newVariable((vit->first)->getLb(), (vit->first)->getUb(), Continuous);
          oVars_.insert(std::make_pair(vit->first, newVar));
          rev_oVars_.insert(std::make_pair(newVar, vit->first));
          max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
        }
      }
    }

    if(qcf) {
      for(VariablePairGroupConstIterator qit = qcf->begin(); qit != qcf->end(); ++qit) {
        ConstVariablePtr x1 = qit->first.first;
        ConstVariablePtr x2 = qit->first.second;
        ovar_it = oVars_.find(x1);
        if(ovar_it == oVars_.end()) {
          VariablePtr newVar;
          newVar = relaxation->newVariable(x1->getLb(), x1->getUb(), Continuous);
          oVars_.insert(std::make_pair(x1, newVar));
          rev_oVars_.insert(std::make_pair(newVar, x1));
          max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
        }
        ovar_it = oVars_.find(x2);
        if(ovar_it == oVars_.end()) {
          VariablePtr newVar;
          newVar = relaxation->newVariable(x2->getLb(), x2->getUb(), Continuous);
          oVars_.insert(std::make_pair(x2, newVar));
          rev_oVars_.insert(std::make_pair(newVar, x2));
          max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
        }
      }
    }

    if(pcf) {
      for(MonomialConstIter mit = pcf->termsBegin(); mit != pcf->termsEnd(); ++mit) {
        for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
          ovar_it = oVars_.find(tit->first);
          if(ovar_it == oVars_.end()) {
            VariablePtr newVar;
            newVar = relaxation->newVariable((tit->first)->getLb(), (tit->first)->getUb(), Continuous);
            oVars_.insert(std::make_pair(tit->first, newVar));
            rev_oVars_.insert(std::make_pair(newVar, tit->first));
            max_pow_.insert(std::pair<ConstVariablePtr, int>(newVar, 0));
          }
        }
      }
    } 
  }

  std::vector<int> mlcid;
  int cix = 0;

  // Go through the multilinear constraints
  //  - for each variable find the highest power
  //  - for each monomial if it hasn't been seen before, add it to the list of monomials (monomial_terms_)
  for(ConstraintConstIterator it = workingProblem_->consBegin(); it != workingProblem_->consEnd(); ++it) {
    FunctionPtr constraintF = (*it)->getFunction();
    FunctionType ft = (*it)->getFunctionType();

    // Check the type of the function
    // - Linear: add the constraint to the relaxation
    // - Bilinear, Quadratic, and Polynomial: mark the constraint
    if (ft == Quadratic || ft == Bilinear || ft == Polynomial) {
      mlcid.push_back((*it)->getId());
      //c_list[cix] = true;
    }
    if (ft == Linear) {
      LinearFunctionPtr lf = LinearFunctionPtr(new LinearFunction());
      const LinearFunctionPtr olf = (*it)->getLinearFunction();
      for(VariableGroupConstIterator vit = olf->termsBegin(); vit != olf->termsEnd(); ++vit) {
        ConstVariablePtr tempVar = vit->first;
        lf->addTerm(oVars_.find(tempVar)->second, vit->second);
      }

      FunctionPtr originalLinearF = FunctionPtr(new Function(lf));
      ConstraintPtr originalLinearC = relaxation->newConstraint(originalLinearF, (*it)->getLb(), (*it)->getUb());
    }
    cix++;
  
    // If the function is Quadratic
    // - Find the maximum exponent of each variable
    if(ft == Quadratic) {
      LinearFunctionPtr lcf;
      QuadraticFunctionPtr qcf;
      lcf = constraintF->getLinearFunction();
      qcf = constraintF->getQuadraticFunction();
      if(qcf) {
        for(VariablePairGroupConstIterator qit = qcf->begin(); qit != qcf->end(); ++qit) {
          ConstVariablePtr x1 = qit->first.first;
          ConstVariablePtr x2 = qit->first.second;
          max_pow_it_ = max_pow_.find(oVars_.find(x1)->second); 
          if(x1->getId()==x2->getId()) {
            if(max_pow_it_->second < 2 ) {
              max_pow_it_->second = 2;
            }
          }
        }
      }
    }

    // If the function is Polynomial
    // - Find the maximum exponent of each variable
    // - Update the list of monomial terms
       
    if(ft == Polynomial) {
      LinearFunctionPtr lcf;
      QuadraticFunctionPtr qcf;
      NonlinearFunctionPtr nlcf = constraintF->getNonlinearFunction();
      PolyFunPtr pcf = boost::dynamic_pointer_cast <PolynomialFunction> (nlcf);
      lcf = constraintF->getLinearFunction();
      qcf = constraintF->getQuadraticFunction();
      
      // get the all quadratic part (no linear left) of the constraint function
      // update the max_pow_ map for the quadratic part
      if(qcf) {
        for(VariablePairGroupConstIterator qit = qcf->begin(); qit != qcf->end(); ++qit) {
          ConstVariablePtr x1 = qit->first.first;
          ConstVariablePtr x2 = qit->first.second;
          max_pow_it_ = max_pow_.find(oVars_.find(x1)->second); 
          if(x1->getId()==x2->getId()) {
            if(max_pow_it_->second < 2 ) {
              max_pow_it_->second = 2;
            }
          }
        }
      }
      
      // For the polynomial part update the max_pow_ and make a list of monomial terms
      if(pcf) {
        for(MonomialConstIter mit = pcf->termsBegin(); mit != pcf->termsEnd(); ++mit) {
          const VarIntMap *mit_terms = (*mit)->getTerms();
          if(monomial_terms_.find(*mit_terms) == monomial_terms_.end()) {
            // Find the bounds on the monomial
            double lb1 = 1;
            double ub1 = 1;
            double lb = 0;
            double ub = 0;
            for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
              for(UInt i = 1; i <= tit->second; i++) {
                BoundsOnProduct((tit->first)->getLb(), (tit->first)->getUb(), lb1, ub1, lb, ub);
                lb1 = lb;
                ub1 = ub;
              }
            }

            ConstVariablePtr y = relaxation->newVariable(lb, ub, Continuous);
            monomial_terms_.insert(std::pair<VarIntMap, ConstVariablePtr>(*mit_terms, y));
          }
          for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
            max_pow_it_ = max_pow_.find(oVars_.find(tit->first)->second);
            if(max_pow_it_->second < tit->second && tit->second >= 2) { 
              max_pow_it_->second = tit->second;
            }
          }
        }
      }
    }
  }

  // For each variable, based on its maximum power, build the right number of new variables
  // - For instance if max_pow_ for x1 is 3
  //   make a vector containing x1_1 and x1_2
  std::vector<ConstVariablePtr> relaxation_vars;
  for(VariableConstIterator vit = relaxation->varsBegin(); vit != relaxation->varsEnd(); ++vit) {
    relaxation_vars.push_back(*vit);
  }

  for(std::vector<ConstVariablePtr>::iterator vit = relaxation_vars.begin(); vit != relaxation_vars.end(); ++vit) {
    std::vector<ConstVariablePtr> copyV;
    max_pow_it_ = max_pow_.find(*vit);
    if(max_pow_it_ != max_pow_.end()) {
      ConstVariablePtr v = max_pow_it_->first;
      int vMaxPow = max_pow_it_->second;
      if(vMaxPow > 1) {
        for(int i=2; i<=vMaxPow; i++) {
          ConstVariablePtr newV = relaxation->newVariable(v->getLb(), v->getUb(), Continuous);
          copyV.push_back(newV);

          // Add a constraint that says v = newV
          LinearFunctionPtr slf = LinearFunctionPtr(new LinearFunction());
          slf->addTerm(v,-1.0);
          slf->addTerm(newV, 1.0);
          FunctionPtr sof = (FunctionPtr) new Function(slf);
          ConstraintPtr soc = relaxation->newConstraint(sof, 0.0, 0.0);
        }
      }
      newCopyVariables_.insert(std::pair<ConstVariablePtr, std::vector<ConstVariablePtr> > (v, copyV));
    }
  }

  // Get the objective function of the problem and add it to the relaxation
  LinearFunctionPtr rlf = LinearFunctionPtr(new LinearFunction());
  const ObjectivePtr obj = workingProblem_->getObjective();
  const LinearFunctionPtr objlf = obj->getLinearFunction();
  if (objlf){
    for(VariableGroupConstIterator it = objlf->termsBegin(); it != objlf->termsEnd();++it) {
      ConstVariablePtr tempVar;
      tempVar = it->first;
      rlf->addTerm(oVars_.find(tempVar)->second, it->second);
    }
  }
  FunctionPtr of = FunctionPtr(new Function(rlf));
  ObjectivePtr oo = relaxation->newObjective(of, 0, Minimize);

  // Go through multilinear rows.  Add constraints, and create maps for product vars
  for(UInt i = 0; i < mlcid.size(); i++) {
    const ConstraintPtr omlc = workingProblem_->getConstraint(mlcid[i]);
    const LinearFunctionPtr olf = omlc->getLinearFunction();
    const QuadraticFunctionPtr oqf = omlc->getQuadraticFunction();
    const NonlinearFunctionPtr onlf = omlc->getNonlinearFunction();
    PolyFunPtr opf = boost::dynamic_pointer_cast <PolynomialFunction> (onlf);

    LinearFunctionPtr lf = LinearFunctionPtr(new LinearFunction());
    // Linear part of constraint remains the same
    if (olf){
      for(VariableGroupConstIterator it = olf->termsBegin(); it != olf->termsEnd();++it) {
        ConstVariablePtr tempVar;
        tempVar = it->first;
        lf->addTerm(oVars_.find(tempVar)->second, it->second);
      }
    }
  
    // Quadratic part gets a new variable for every term
    if(oqf) {
      for(VariablePairGroupConstIterator it = oqf->begin(); it != oqf->end(); ++it) {
        
        ConstVariablePtr x1 = oVars_.find(it->first.first)->second;
        ConstVariablePtr x2 = oVars_.find(it->first.second)->second;
        ConstVariablePtr x_1;
        // Check to see if 'it' is a square term
        if(x1->getId() == x2->getId()) {
          newCopyVariables_it_ = newCopyVariables_.find(x1);
          x_1 = ((newCopyVariables_it_)->second)[0];

          // look to see if x1^2 already exists...
          std::map <ConstVariablePtr, ConstVariablePair>::iterator pos ;
          ConstVariablePtr y;
          ConstVariablePair Z;
          if (sqterms_.find(x1) == sqterms_.end()) {
            
            // Make a copy of the variable
            double lb = 0.0;
            double ub = 0.0;
            BoundsOnProduct(x1,x_1,lb,ub);
            y = relaxation->newVariable(lb, ub, Continuous);
            blterms_.insert(make_pair(y,ConstVariablePair(x1,x_1)));
            //            blterms_coef_.insert(make_pair(ConstVariablePair(x1,x_1), it->second));

            rev_blterms_.insert(make_pair(ConstVariablePair(x1,x_1), y));
            
            sqterms_.insert(make_pair(x1,ConstVariablePair(y,x_1)));
            rev_sqterms_.insert(make_pair(ConstVariablePair(y,x_1), x1));
          }
          else {
            pos = sqterms_.find(x1);
            y = pos->second.first;
          }
          
          lf->addTerm(y, it->second);
          
          /*
          // Add this term to the list of square terms
          sqterms_.insert(make_pair(y,ConstVariablePair(x1,x_1)));
          rev_sqterms_.insert(make_pair(ConstVariablePair(x1,x_1), y));
          */
        }
        else {
          //Bounds on product depend on whether variable bounds are < 0, > 0
          double lb = 0.0;
          double ub = 0.0;
          BoundsOnProduct(x1,x2,lb,ub);
          
          // look to see if w var already exists...
          std::map <ConstVariablePair, ConstVariablePtr>::iterator pos ;
          ConstVariablePtr w;
          if (rev_blterms_.find(ConstVariablePair(x1,x2)) == rev_blterms_.end()) {
            w = relaxation->newVariable(lb, ub, Continuous);
            
            blterms_.insert(make_pair(w,ConstVariablePair(x1,x2)));
            //            blterms_coef_.insert(make_pair(ConstVariablePair(x1,x2), it->second));

            rev_blterms_.insert(make_pair(ConstVariablePair(x1,x2), w));
          }
          else {
            pos = rev_blterms_.find(ConstVariablePair(x1,x2));
            w = pos->second;
          }
          
          lf->addTerm(w,it->second);
        }
      }
    }


    // Take care of the polynomial parts
    // - Go through the polynomial parts of the constraint
    //   form the linearizations
    
    if(opf) {
      for(MonomialConstIter mit = opf->termsBegin(); mit != opf->termsEnd(); ++mit) {
        
        const VarIntMap* mit_terms = (*mit)->getTerms();
        double mit_coeff = (*mit)->getCoeff();
        monomial_terms_it_ = monomial_terms_.find(*mit_terms);
        if(monomial_terms_it_ != monomial_terms_.end()) {
          std::vector<ConstVariablePtr> mit_ml;
          for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
            ConstVariablePtr relax_tit_var = (oVars_.find(tit->first))->second;
            if(tit->second == 1) {
              mit_ml.push_back(relax_tit_var);
            }
            if(tit->second >= 2) {
              mit_ml.push_back(relax_tit_var);
              newCopyVariables_it_ = newCopyVariables_.find(relax_tit_var);
              for(UInt i=2; i <= tit->second; i++) {
                mit_ml.push_back((newCopyVariables_it_->second)[i-2]);
              }
            }
          }

          mlterms_.insert(std::pair<ConstVariablePtr, std::vector<ConstVariablePtr> > (monomial_terms_it_->second, mit_ml));
          rev_mlterms_.insert(std::pair<std::vector<ConstVariablePtr>, ConstVariablePtr > (mit_ml, monomial_terms_it_->second));
        }
        
        bool termAdded = 0;
        for(VariableGroupConstIterator lf_it = lf->termsBegin(); lf_it != lf->termsEnd(); ++lf_it) {
          if(lf_it->first == monomial_terms_it_->second) {
            lf->incTerm(monomial_terms_it_->second, mit_coeff);
            termAdded = 1;
            break;
          }
        }
        if(!termAdded)
          lf->addTerm(monomial_terms_it_->second, mit_coeff);
      }
    }

    FunctionPtr of = (FunctionPtr) new Function(lf);
    ConstraintPtr oc = relaxation->newConstraint(of,omlc->getLb(), omlc->getUb());
  }

  // Add linearizations for square terms
  for(std::map <ConstVariablePtr, ConstVariablePair>::iterator sqiter = sqterms_.begin();
      sqiter != sqterms_.end(); ++sqiter) {
    for(int i=0; i<linearizationCnt_; i++) {
      ConstVariablePtr y;
      ConstVariablePtr x1;
      x1 = sqiter->first;
      y = sqiter->second.first;
      LinearFunctionPtr sCutlf = LinearFunctionPtr(new LinearFunction());
      sCutlf = LinearFunctionPtr(new LinearFunction());
      double x1Val = x1->getLb()+i*(x1->getUb()-x1->getLb())/(linearizationCnt_-1);
      sCutlf->addTerm(x1, 2*x1Val);
      sCutlf->addTerm(y, -1);
      
      FunctionPtr sCutf = (FunctionPtr) new Function(sCutlf);
      ConstraintPtr sCutc = relaxation->newConstraint(sCutf, -INFINITY,
                                                      x1Val*x1Val );
    }
  }

  // Now add all the constraints for each new bilinear term
  for(std::map<ConstVariablePtr, ConstVariablePair>::iterator it = blterms_.begin();
      it != blterms_.end(); ++it) {
    ConstVariablePtr w = it->first;
    ConstVariablePtr x1 = it->second.first;
    ConstVariablePtr x2 = it->second.second;
    LinearFunctionPtr lf1 = LinearFunctionPtr(new LinearFunction());
    LinearFunctionPtr lf2 = LinearFunctionPtr(new LinearFunction());
    LinearFunctionPtr lfw = LinearFunctionPtr(new LinearFunction());

    VariablePtr lamll = relaxation->newVariable(0.0, 1.0, Continuous);
    VariablePtr lamul = relaxation->newVariable(0.0, 1.0, Continuous);
    VariablePtr lamuu = relaxation->newVariable(0.0, 1.0, Continuous);
    VariablePtr lamlu = relaxation->newVariable(0.0, 1.0, Continuous);

    // Just enumerate extreme points yourself
    lf1->addTerm(x1,-1.0);
    lf1->addTerm(lamll, x1->getLb());
    lf1->addTerm(lamul, x1->getUb());
    lf1->addTerm(lamuu, x1->getUb());
    lf1->addTerm(lamlu, x1->getLb());

    // Just enumerate extreme points yourself
    lf2->addTerm(x2,-1.0);
    lf2->addTerm(lamll, x2->getLb());
    lf2->addTerm(lamul, x2->getLb());
    lf2->addTerm(lamuu, x2->getUb());
    lf2->addTerm(lamlu, x2->getUb());

    lfw->addTerm(w, -1.0);
    lfw->addTerm(lamll, x1->getLb()*x2->getLb());
    lfw->addTerm(lamul, x1->getUb()*x2->getLb());
    lfw->addTerm(lamuu, x1->getUb()*x2->getUb());
    lfw->addTerm(lamlu, x1->getLb()*x2->getUb());

    // Add the x1,x2,and w rows
    FunctionPtr f1 = (FunctionPtr) new Function(lf1);
    ConstraintPtr c1 = relaxation->newConstraint(f1, 0.0, 0.0);

    FunctionPtr f2 = (FunctionPtr) new Function(lf2);
    ConstraintPtr c2 = relaxation->newConstraint(f2, 0.0, 0.0);

    FunctionPtr fw = (FunctionPtr) new Function(lfw);
    ConstraintPtr cw = relaxation->newConstraint(fw, 0.0, 0.0);

    // Add the convexity constraint
    LinearFunctionPtr convex_lf =  LinearFunctionPtr(new LinearFunction());
    convex_lf->addTerm(lamll, 1.0);
    convex_lf->addTerm(lamul, 1.0);
    convex_lf->addTerm(lamuu, 1.0);
    convex_lf->addTerm(lamlu, 1.0);

    FunctionPtr convex_f = (FunctionPtr) new Function(convex_lf);
    ConstraintPtr convex_c = relaxation->newConstraint(convex_f, 1.0, 1.0);
  }

  // Now add all the constraints for multilinear terms
  // - Go through the multilinear terms
  for(mlterms_it_ = mlterms_.begin(); mlterms_it_ != mlterms_.end(); mlterms_it_++) {
    ConstVariablePtr x1 = (mlterms_it_->second)[0];
    ConstVariablePtr w;
    for(UInt i = 1; i < (mlterms_it_->second).size(); i++) {
      ConstVariablePtr x2 = (mlterms_it_->second)[i];
      double wl;
      double wu;
      BoundsOnProduct(x1, x2, wl, wu);
      w = relaxation->newVariable(wl, wu, Continuous); 

      LinearFunctionPtr lf1 = LinearFunctionPtr(new LinearFunction());
      LinearFunctionPtr lf2 = LinearFunctionPtr(new LinearFunction());
      LinearFunctionPtr lfw = LinearFunctionPtr(new LinearFunction());
      
      VariablePtr lamll = relaxation->newVariable(0.0, 1.0, Continuous);
      VariablePtr lamul = relaxation->newVariable(0.0, 1.0, Continuous);
      VariablePtr lamuu = relaxation->newVariable(0.0, 1.0, Continuous);
      VariablePtr lamlu = relaxation->newVariable(0.0, 1.0, Continuous);

      // Just enumerate extreme points yourself
      lf1->addTerm(x1,-1.0);
      lf1->addTerm(lamll, x1->getLb());
      lf1->addTerm(lamul, x1->getUb());
      lf1->addTerm(lamuu, x1->getUb());
      lf1->addTerm(lamlu, x1->getLb());
      
      // Just enumerate extreme points yourself
      lf2->addTerm(x2,-1.0);
      lf2->addTerm(lamll, x2->getLb());
      lf2->addTerm(lamul, x2->getLb());
      lf2->addTerm(lamuu, x2->getUb());
      lf2->addTerm(lamlu, x2->getUb());
      
      lfw->addTerm(w, -1.0);
      lfw->addTerm(lamll, x1->getLb()*x2->getLb());
      lfw->addTerm(lamul, x1->getUb()*x2->getLb());
      lfw->addTerm(lamuu, x1->getUb()*x2->getUb());
      lfw->addTerm(lamlu, x1->getLb()*x2->getUb());
      
      // Add the x1,x2,and w rows
      FunctionPtr f1 = (FunctionPtr) new Function(lf1);
      ConstraintPtr c1 = relaxation->newConstraint(f1, 0.0, 0.0);
      
      FunctionPtr f2 = (FunctionPtr) new Function(lf2);
      ConstraintPtr c2 = relaxation->newConstraint(f2, 0.0, 0.0);
      
      FunctionPtr fw = (FunctionPtr) new Function(lfw);
      ConstraintPtr cw = relaxation->newConstraint(fw, 0.0, 0.0);
      
      
      // Add the convexity constraint
      LinearFunctionPtr convex_lf =  LinearFunctionPtr(new LinearFunction());
      convex_lf->addTerm(lamll, 1.0);
      convex_lf->addTerm(lamul, 1.0);
      convex_lf->addTerm(lamuu, 1.0);
      convex_lf->addTerm(lamlu, 1.0);
      
      FunctionPtr convex_f = (FunctionPtr) new Function(convex_lf);
      ConstraintPtr convex_c = relaxation->newConstraint(convex_f, 1.0, 1.0);

      x1 = w;
    }
    // Add the constraint 'mlterm' = w
    LinearFunctionPtr mlterm_lf = LinearFunctionPtr(new LinearFunction());
    mlterm_lf->addTerm(mlterms_it_->first, -1.0);
    mlterm_lf->addTerm(w, 1.0);

    FunctionPtr mlterm_f = (FunctionPtr) new Function(mlterm_lf);
    ConstraintPtr mlterm_c = relaxation->newConstraint(mlterm_f, 0.0, 0.0);
  }

  should_prune = false;
}

void
MultilinearHandler::getMultilinearTerms(std::map <ConstVariablePtr, ConstVariablePair> blterms,
                             std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms,
                             UInt maxGroupSize,
                             std::vector<std::vector<ConstVariablePtr> > &terms)
{
  // Form a 2-D vector that has the bilinear and multilinear terms
  for(std::map<ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator ml_it = mlterms.begin(); ml_it != mlterms.end(); ++ml_it) {
    std::vector<ConstVariablePtr> temp;
    for(std::vector<ConstVariablePtr>::iterator t_it = ml_it->second.begin(); t_it != ml_it->second.end(); ++t_it) {
      temp.push_back(*t_it);
    }
    terms.push_back(temp);
    // if the term size is bigger than the max group size, quit
    assert(maxGroupSize > temp.size());
  }
  
  for(std::map<ConstVariablePtr, ConstVariablePair>::iterator bl_it = blterms.begin(); bl_it != blterms.end(); ++bl_it) {
    std::vector<ConstVariablePtr> temp;
    temp.push_back(bl_it->second.first);
    temp.push_back(bl_it->second.second);
    terms.push_back(temp);
  }
}












void
MultilinearHandler::termsAppearKtimes(std::vector<std::vector<ConstVariablePtr> > terms,
                                      std::vector <double> termsCoef,
                                      std::vector <std::vector<ConstVariablePtr> > groups,
                                      int k,
                                      std::vector<std::vector<ConstVariablePtr> > &termsk,                       
                                      std::vector <double> &termsKCoef,
                                      std::vector<int> &termRep)
{
  for(UInt m=0; m<terms.size(); m++) {
    int termCnt = 0;
    ConstVariablePtr x1 = terms[m][0];
    ConstVariablePtr x2 = terms[m][1];
    for(UInt i=0; i<groups.size(); i++) {
      bool x1InGroup = 0;
      bool x2InGroup = 0;
      for(std::vector<ConstVariablePtr>::iterator it1 = groups[i].begin(); it1 != groups[i].end(); ++it1) {
        if(x1->getId() == (*it1)->getId())
          x1InGroup = 1;
        if(x2->getId() == (*it1)->getId())
          x2InGroup = 1;        
      }
      if(x1InGroup && x2InGroup) {
        termCnt++;
      }
    }
    if(termCnt == k) {
      termsk.push_back(terms[m]);
      termsKCoef.push_back(termsCoef[m]);
    }
    termRep.push_back(termCnt);
  }  
}

void 
MultilinearHandler::countTermsAppearance(std::vector<std::vector<ConstVariablePtr> > terms,
                                         std::vector <std::vector<ConstVariablePtr> > groups,
                                         std::vector<int> &termRep)
{
  for(UInt i=0; i<terms.size(); i++) {
    int termCnt = 0;
    ConstVariablePtr x1 = terms[i][0];
    ConstVariablePtr x2 = terms[i][1];
    for(UInt i=0; i<groups.size(); i++) {
      bool x1InGroup = 0;
      bool x2InGroup = 0;
      for(std::vector<ConstVariablePtr>::iterator it1 = groups[i].begin(); it1 != groups[i].end(); ++it1) {
        if(x1->getId() == (*it1)->getId())
          x1InGroup = 1;
        if(x2->getId() == (*it1)->getId())
          x2InGroup = 1;        
      }
      if(x1InGroup && x2InGroup)
        termCnt++;
    }
    termRep.push_back(termCnt);
  }
}

void
MultilinearHandler::getMultilinearVariables(std::map <ConstVariablePtr, ConstVariablePair> blterms,
                                            std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms,
                                            std::vector<ConstVariablePtr> &mlVars)
{
  // Get the variables that appear in the bilinear terms
  for(std::map<ConstVariablePtr, ConstVariablePair>::iterator it = blterms.begin();
      it != blterms.end(); ++it) {
    
    ConstVariablePtr x1 = it->second.first;
    ConstVariablePtr x2 = it->second.second;
    bool x1InTemp = 0;
    bool x2InTemp = 0;
    for(std::vector<ConstVariablePtr>::iterator it1 = mlVars.begin(); it1 != mlVars.end(); ++it1) {
      if(x1->getId() == (*it1)->getId())
        x1InTemp = 1;
      if(x2->getId() == (*it1)->getId())
        x2InTemp = 1;
    }
    if(!x1InTemp) 
      mlVars.push_back(x1);

    if(!x2InTemp) 
      mlVars.push_back(x2);
  }
  
  // Get the variables that appear in multilinear terms
  for(std::map<ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator it = mlterms.begin();
      it != mlterms.end(); ++it) {
    for(std::vector<ConstVariablePtr>::iterator mlt_it = it->second.begin(); mlt_it != it->second.end(); ++mlt_it) {
      int varInTemp = 0;
      for(std::vector<ConstVariablePtr>::iterator gv_it = mlVars.begin(); gv_it != mlVars.end(); ++gv_it) {
        if((*gv_it)->getId() == (*mlt_it)->getId()) {
          varInTemp = 1;
          break;
        }
      }
      if(!varInTemp) 
        mlVars.push_back(*mlt_it);
    }
  }
}

void
MultilinearHandler::getMultilinearTermsCoef(std::map <ConstVariablePair, double> blterms_coef,
                                             std::vector<std::vector<ConstVariablePtr> > terms,
                                             std::vector<double> &termsCoef)
{
  for(UInt i=0; i<terms.size(); i++) {
    ConstVariablePtr x1 = terms[i][0];
    ConstVariablePtr x2 = terms[i][1];
    for(std::map <ConstVariablePair, double>::iterator coefIter=blterms_coef.begin();
        coefIter != blterms_coef.end(); ++coefIter) {
      bool x1In = 0;
      bool x2In = 0;
      if(x1->getId() == (coefIter->first.first)->getId())
        x1In = 1;
      if(x1->getId() == (coefIter->first.second)->getId())
        x1In = 1;
      if(x2->getId() == (coefIter->first.first)->getId())
        x2In = 1;
      if(x2->getId() == (coefIter->first.second)->getId())
        x2In = 1;

      if(x1In && x2In) {
        termsCoef.push_back(coefIter->second);
      }
    }
  }
}

void
MultilinearHandler::groupTermByTerm(std::map <ConstVariablePtr, ConstVariablePair> blterms,
                                    std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms,
                                    std::vector <std::vector<ConstVariablePtr> > &groups)
{
  for(std::map <ConstVariablePtr, ConstVariablePair>::iterator bit = blterms.begin(); 
      bit != blterms.end(); ++bit) {
    ConstVariablePtr x1 = (bit->second).first;
    ConstVariablePtr x2 = (bit->second).second;
    std::vector<ConstVariablePtr> temp;
    temp.push_back(x1);
    temp.push_back(x2);
    groups.push_back(temp);
  }
  
  for(std::map <ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator mit = mlterms.begin();
      mit != mlterms.end(); ++mit) {
    int vSize = (mit->second).size();
    std::vector<ConstVariablePtr> temp;
    for(int i=0; i<vSize; i++) {
      temp.push_back((mit->second)[i]);
    }
    groups.push_back(temp);
  }
}

void
MultilinearHandler::groupUsingDensest2ND(UInt gs, 
                                         std::map <ConstVariablePair, std::vector<double> > bl_coef,
                                         std::map <std::vector<ConstVariablePtr>, std::vector<double> > ml_coef,
                                         std::vector <ConstVariablePtr> vars,
                                         std::vector <std::vector<ConstVariablePtr> > &groups,
                                         UInt maxGroupsCnt)
{
  int varsCnt = vars.size();
  int *varsId = new int[varsCnt];
  for(int i=0; i<varsCnt; i++)
    varsId[i] = vars[i]->getId();

  std::vector <std::vector<ConstVariablePtr> > terms;
  std::vector <std::vector<double> > termsCoef;
  
  // put the bilinear and multilinear terms and their coef in vectors
  for(std::map <ConstVariablePair, std::vector<double> >::iterator bit = bl_coef.begin(); bit != bl_coef.end(); ++bit) {
    ConstVariablePtr x1 = (bit->first).first;
    ConstVariablePtr x2 = (bit->first).second;
    std::vector<ConstVariablePtr> tempVar;
    tempVar.push_back(x1);
    tempVar.push_back(x2);
    terms.push_back(tempVar);
    
    termsCoef.push_back(bit->second);
  }
  for(std::map <std::vector<ConstVariablePtr>, std::vector<double> >::iterator mit = ml_coef.begin(); mit != ml_coef.end(); ++mit) {
    terms.push_back(mit->first);
    termsCoef.push_back(mit->second);
  }

  int termsCnt = terms.size();  

  // terms variables coinsidence matrix
  double **term_var_matrix = new double* [termsCnt];
  for(int i=0; i<termsCnt; i++)
    term_var_matrix[i] = new double [varsCnt];
  
  // initialize the matrix elements to zero
  for(int i=0; i<termsCnt; i++)
    for(int j=0; j<varsCnt; j++)
      term_var_matrix[i][j] = 0;

  // for the elements of the bilinear terms map, put the coefficients in the matrix
  for(UInt i=0; i<terms.size(); i++) {
    for(UInt j=0; j<terms[i].size(); j++) {
      for(UInt k=0; k<vars.size(); k++) {
        if(vars[k]->getId() == terms[i][j]->getId()) {
          for(UInt l=0; l<termsCoef[i].size(); l++) {
            term_var_matrix[i][k] += fabs(termsCoef[i][l]);
          }
        }
      }
    }
  }

  while(groups.size() < maxGroupsCnt) {
    std::vector<std::vector<int> > comp;
    findGraphComponents(varsCnt, termsCnt, term_var_matrix, comp);

    if(comp.size() == 0)
      break;
    
    double maxDensity = 0;
    std::vector<ConstVariablePtr> temp;
    std::vector<int> tempInd;
    int emptyCompCnt = 0;
    for(UInt i=0; i<comp.size(); i++) {
      double density = 0;
      std::vector<ConstVariablePtr> densest;
      std::vector<int> densestInd;
      bool isCompEmpty = 0;
      findDensestSubgraph(gs, vars, terms, term_var_matrix, comp[i], densest, densestInd, isCompEmpty);
      emptyCompCnt =+ int(isCompEmpty);
      findSubgraphDensity(terms, term_var_matrix, densest, densestInd, density);
      if(density > maxDensity) {
        maxDensity = density;
        temp.clear();
        tempInd.clear();
        // temp = densest;
        for(UInt j=0; j<densest.size(); j++) {
          temp.push_back(densest[j]);
          tempInd.push_back(densestInd[j]);
        }
      }
    }

    groups.push_back(temp);
    reduceSubgraphEdges(termsCnt, term_var_matrix, terms, temp, tempInd);
  }

  // delete memory allocated
  delete [] varsId;

  for(int i = 0; i < termsCnt; i++) {
    delete [] term_var_matrix[i];
  }
  delete [] term_var_matrix;

}

void
MultilinearHandler::groupUsingDensest(UInt gs, 
                                      std::map <ConstVariablePair, std::vector<double> > bl_coef,
                                      std::map <std::vector<ConstVariablePtr>, std::vector<double> > ml_coef,
                                      std::vector <ConstVariablePtr> vars,
                                      std::vector <std::vector<ConstVariablePtr> > &groups)
{
  int varsCnt = vars.size();
  int *varsId = new int[varsCnt];
  for(int i=0; i<varsCnt; i++)
    varsId[i] = vars[i]->getId();

  std::vector <std::vector<ConstVariablePtr> > terms;
  std::vector <std::vector<double> > termsCoef;
  
  // put the bilinear and multilinear terms and their coef in vectors
  for(std::map <ConstVariablePair, std::vector<double> >::iterator bit = bl_coef.begin(); bit != bl_coef.end(); ++bit) {
    ConstVariablePtr x1 = (bit->first).first;
    ConstVariablePtr x2 = (bit->first).second;
    std::vector<ConstVariablePtr> tempVar;
    tempVar.push_back(x1);
    tempVar.push_back(x2);
    terms.push_back(tempVar);
    
    termsCoef.push_back(bit->second);
  }
  for(std::map <std::vector<ConstVariablePtr>, std::vector<double> >::iterator mit = ml_coef.begin(); mit != ml_coef.end(); ++mit) {
    terms.push_back(mit->first);
    termsCoef.push_back(mit->second);
  }

  int termsCnt = terms.size();  

  // terms variables coinsidence matrix
  double **term_var_matrix = new double* [termsCnt];
  for(int i=0; i<termsCnt; i++)
    term_var_matrix[i] = new double [varsCnt];
  
  // initialize the matrix elements to zero
  for(int i=0; i<termsCnt; i++)
    for(int j=0; j<varsCnt; j++)
      term_var_matrix[i][j] = 0;

  // for the elements of the bilinear terms map, put the coefficients in the matrix
  for(UInt i=0; i<terms.size(); i++) {
    for(UInt j=0; j<terms[i].size(); j++) {
      for(UInt k=0; k<vars.size(); k++) {
        if(vars[k]->getId() == terms[i][j]->getId()) {
          for(UInt l=0; l<termsCoef[i].size(); l++) {
            term_var_matrix[i][k] += fabs(termsCoef[i][l]);
          }
        }
      }
    }
  }

  bool isGraphEmpty = 0;
  while(!isGraphEmpty) {
    std::vector<std::vector<int> > comp;
    findGraphComponents(varsCnt, termsCnt, term_var_matrix, comp);
    
    if(comp.size() == 0)
      break;
    
    double maxDensity = 0;
    std::vector<ConstVariablePtr> temp;
    std::vector<int> tempInd;
    UInt emptyCompCnt = 0;
    for(UInt i=0; i<comp.size(); i++) {
      double density = 0;
      std::vector<ConstVariablePtr> densest;
      std::vector<int> densestInd;
      bool isCompEmpty = 0;
      findDensestSubgraph(gs, vars, terms, term_var_matrix, comp[i], densest, densestInd, isCompEmpty);
      emptyCompCnt =+ int(isCompEmpty);
      findSubgraphDensity(terms, term_var_matrix, densest, densestInd, density);
      if(density > maxDensity) {
        maxDensity = density;
        temp.clear();
        tempInd.clear();
        // temp = densest;
        for(UInt j=0; j<densest.size(); j++) {
          temp.push_back(densest[j]);
          tempInd.push_back(densestInd[j]);
        }
      }
    }
    //XXX
    if(emptyCompCnt == comp.size() && comp.size() == 1)
      isGraphEmpty = 1;
    groups.push_back(temp);
    removeSubgraphEdges(termsCnt, term_var_matrix, terms, temp, tempInd);
  }

  delete [] varsId;

  for(int i=0; i<termsCnt; i++)
    delete [] term_var_matrix[i];

  delete [] term_var_matrix;

}

void
MultilinearHandler::findSubgraphDensity(std::vector<std::vector<ConstVariablePtr> > terms,
                                        double** term_var_matrix,
                                        std::vector<ConstVariablePtr> nodes,
                                        std::vector<int> nodesInd,
                                        double &density)
{
  int termsCnt = terms.size();
  for(int i=0; i<termsCnt; i++) {
    UInt isTermIn = 0;
    for(UInt j=0; j<terms[i].size(); j++) {
      for(UInt k=0; k<nodes.size(); k++) {
        if(terms[i][j]->getId() == nodes[k]->getId()) {
          isTermIn++;
        }
      }
    }
    if(isTermIn >= terms[i].size()) {
      for(UInt j=0; j<terms[i].size(); j++) {
        for(UInt k=0; k<nodes.size(); k++) {
          if(nodes[k]->getId() == terms[i][j]->getId()) {
            density += term_var_matrix[i][nodesInd[k]];
          }
        }
      }
    }
  }
}

void
MultilinearHandler::removeSubgraphEdges(int termsCnt,
                                        double** &term_var_matrix,
                                        std::vector<std::vector<ConstVariablePtr> > terms,
                                        std::vector<ConstVariablePtr> nodes,
                                        std::vector<int> nodesInd)
{
  for(int i=0; i<termsCnt; i++) {
    UInt isTermIn = 0;
    for(UInt j=0; j<terms[i].size(); j++) {
      for(UInt k=0; k<nodes.size(); k++) {
        if(terms[i][j]->getId() == nodes[k]->getId()) {
          isTermIn++;
        }
      }
    }
    if(isTermIn >= terms[i].size()) {
      for(UInt j=0; j<terms[i].size(); j++) {
        for(UInt k=0; k<nodes.size(); k++) {
          if(nodes[k]->getId() == terms[i][j]->getId()) {
            term_var_matrix[i][nodesInd[k]] = 0;
          }
        }
      }
    }
  }
}

void
MultilinearHandler::reduceSubgraphEdges(int termsCnt,
                                        double** &term_var_matrix,
                                        std::vector<std::vector<ConstVariablePtr> > terms,
                                        std::vector<ConstVariablePtr> nodes,
                                        std::vector<int> nodesInd)
{
  for(int i=0; i<termsCnt; i++) {
    UInt isTermIn = 0;
    for(UInt j=0; j<terms[i].size(); j++) {
      for(UInt k=0; k<nodes.size(); k++) {
        if(terms[i][j]->getId() == nodes[k]->getId()) {
          isTermIn++;
        }
      }
    }
    if(isTermIn >= terms[i].size()) {
      for(UInt j=0; j<terms[i].size(); j++) {
        for(UInt k=0; k<nodes.size(); k++) {
          if(nodes[k]->getId() == terms[i][j]->getId()) {
            term_var_matrix[i][nodesInd[k]] *= 0.5;
          }
        }
      }
    }
  }
}

void
MultilinearHandler::findDensestSubgraph(UInt gs,
                                        std::vector<ConstVariablePtr> vars,
                                        std::vector<std::vector<ConstVariablePtr> > terms,
                                        double** term_var_matrix, 
                                        std::vector<int> component, 
                                        std::vector<ConstVariablePtr> &densest,
                                        std::vector<int> &densestInd,
                                        bool &isCompEmpty) 
{
  int varsCnt = vars.size();
  int termsCnt = terms.size();

  bool *termActive = new bool[termsCnt];
  for(int i=0; i<termsCnt; i++) {
    UInt termIsIn = 0;
    for(UInt j=0; j<terms[i].size(); j++) {
      for(UInt k=0; k<component.size(); k++) {
        if(vars[component[k]]->getId() == terms[i][j]->getId()) {
          termIsIn++;
          break;
        }
      }
    }
    if(termIsIn == terms[i].size())
      termActive[i] = 1;
    else
      termActive[i] = 0;
  }
  
  // Keep the index of the remaining variables
  bool* remainNodes = new bool [varsCnt];
  UInt remainNodesCnt = 0;
  for(int i=0; i<varsCnt; i++)
    remainNodes[i] = false;

  for(UInt i=0; i<component.size(); i++) {
    bool nodeEmpty = 1;
    for(int j=0; j<termsCnt; j++) {
      if(termActive[j]) {
        if(term_var_matrix[j][component[i]] > 1e-15) {
          nodeEmpty = 0;
          break;
        }
      }
    }
    if(!nodeEmpty) {
      remainNodes[component[i]] = true;
      remainNodesCnt++;
    }
  }

  bool groupAdded = 0;
  // if number of remaining nodes is less than or equal to the group size, just add them as the last group
  if(remainNodesCnt <= gs) {
    groupAdded = 1;
    isCompEmpty = 1;
    for(int i=0; i<varsCnt; i++) {
      if(remainNodes[i]) {
        densest.push_back(vars[i]);
        densestInd.push_back(i);
      }
    }
  }

  // local copy of term_var_matrix
  double **tempMatrix = new double* [termsCnt];
  for(int i=0; i<termsCnt; i++) 
    tempMatrix[i] = new double [varsCnt];
  for(int i=0; i<termsCnt; i++)
    for(int j=0; j<varsCnt; j++)
      tempMatrix[i][j] = term_var_matrix[i][j];

  while(!groupAdded) {
    // find the variable with the least weighted degree
    double minWD = 1e10;
    int minWDInd = -1;
    for(int i=0; i<varsCnt; i++) {
      double WD=0;
      if(remainNodes[i]) {
        for(int j=0; j<termsCnt; j++) {
          if(termActive[j]) {
            WD += tempMatrix[j][i];
          }
        }
        if(minWD > WD) {
          minWD = WD;
          minWDInd = i;
        }
      }
    }

    // zero out all the coefficients of this variable
    for(int i=0; i<termsCnt; i++) {
      if(termActive[i]) {
        if(tempMatrix[i][minWDInd] > 1e-15) {
          for(int j=0; j<varsCnt; j++) {
            tempMatrix[i][j] = 0;
          }
        }
      }
    }

    // update the list of remaining nodes
    remainNodesCnt = 0;
    for(int i=0; i<varsCnt; i++)
      remainNodes[i] = 0;
    
    for(UInt i=0; i<component.size(); i++) {
      bool nodeEmpty = 1;
      for(int j=0; j<termsCnt; j++) {
        if(termActive[j]) {
          if(tempMatrix[j][component[i]] > 1e-15) {
            nodeEmpty = 0;
            break;
          }
        }
      }
      if(!nodeEmpty) {
        remainNodes[component[i]] = 1;
        remainNodesCnt++;
      }
    }

    if(remainNodesCnt <= gs) {
      // add the remaining nodes as the densest subgraph
      for(int i=0; i<varsCnt; i++) {
        if(remainNodes[i]) {
          densest.push_back(vars[i]);
          densestInd.push_back(i);
        }
      }
      groupAdded = 1;
    }
  }

  delete [] termActive;
  delete [] remainNodes;

  for(int i=0; i<termsCnt; i++) {
    delete [] tempMatrix[i];
  }
  delete [] tempMatrix;

}

void
MultilinearHandler::findGraphComponents(int varsCnt,
                                        int termsCnt,
                                        double** term_var_matrix, 
                                        std::vector<std::vector<int> > &components)
{
  std::vector<int> visited;
  for(int i=0; i<varsCnt; i++) {
    bool nodeIsVisited = 0;
    for(UInt j=0; j<visited.size(); j++) {
      if(i == visited[j]) {
        nodeIsVisited = 1;
        break;
      }
    }
    if(!nodeIsVisited) {
      std::vector<int> temp;
      temp.push_back(i);
      visited.push_back(i);
      UInt iter = 0;
      while(iter < temp.size()) {
        for(int j=0; j<termsCnt; j++) {
          if(term_var_matrix[j][temp[iter]] > 1e-15) {
            for(int k=0; k<varsCnt; k++) {
              if(k!=temp[iter] && term_var_matrix[j][k] > 1e-15) {
                nodeIsVisited = 0;
                for(UInt l=0; l<visited.size(); l++) {
                  if(k == visited[l]) {
                    nodeIsVisited = 1;
                    break;
                  }
                }
                if(!nodeIsVisited) {
                  temp.push_back(k);
                  visited.push_back(k);
                }
              }
            }
          }
        }
        iter++;
      }
      if(temp.size() > 1)
        components.push_back(temp);
    }
  }
}


void
MultilinearHandler::findGroupDensity(std::vector<std::vector<ConstVariablePtr> > terms,
                                     std::vector <std::vector<ConstVariablePtr> > groups,
                                     std::vector <int> &density)
{
  for(UInt i=0; i<groups.size(); i++) {
    int numTermsInGroup = 0;
    for(UInt j=0; j<groups[i].size(); j++) {
      for(UInt k=j+1; k<groups[i].size(); k++) {
        ConstVariablePtr x1 = groups[i][j];
        ConstVariablePtr x2 = groups[i][k];
        for(UInt l=0; l<terms.size(); l++) {
          bool x1InTerm = 0;
          bool x2InTerm = 0;
          if(x1->getId() == terms[l][0]->getId() || x1->getId() == terms[l][1]->getId())
            x1InTerm = 1;
          if(x2->getId() == terms[l][0]->getId() || x2->getId() == terms[l][1]->getId())
            x2InTerm = 1;
          if(x1InTerm && x2InTerm) {
            numTermsInGroup ++;
            break;
          }          
        }
      }
    }
    density.push_back(numTermsInGroup);
  }
}

void
MultilinearHandler::findGroupCoef(std::vector<std::vector<ConstVariablePtr> > terms,
                                  std::vector <double> termsCoef,
                                  std::vector <std::vector<ConstVariablePtr> > groups,
                                  std::vector <double> &sumCoef)
{
  for(UInt i=0; i<groups.size(); i++) {
    double sumCoefTermsInGroup = 0;
    for(UInt j=0; j<groups[i].size(); j++) {
      for(UInt k=j+1; k<groups[i].size(); k++) {
        ConstVariablePtr x1 = groups[i][j];
        ConstVariablePtr x2 = groups[i][k];
        for(UInt l=0; l<terms.size(); l++) {
          bool x1InTerm = 0;
          bool x2InTerm = 0;
          if(x1->getId() == terms[l][0]->getId() || x1->getId() == terms[l][1]->getId())
            x1InTerm = 1;
          if(x2->getId() == terms[l][0]->getId() || x2->getId() == terms[l][1]->getId())
            x2InTerm = 1;
          if(x1InTerm && x2InTerm) {
            sumCoefTermsInGroup += termsCoef[l];
            break;
          }          
        }
      }
    }
    sumCoef.push_back(sumCoefTermsInGroup);
  }
}


void
MultilinearHandler::makeGroups(std::map <ConstVariablePtr, ConstVariablePair> blterms, 
                               std::map <ConstVariablePair, ConstVariablePtr> , 
                               std::map <ConstVariablePair, std::vector<double> > blterms_coef,
                               std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms, 
                               std::map <std::vector<ConstVariablePtr>, ConstVariablePtr> , 
                               std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_coef,
                               std::map <ConstVariablePtr, ConstVariablePair> , 
                               std::map <ConstVariablePair, ConstVariablePtr> , 
                               std::map <ConstVariablePair, std::vector<double> > blterms_obj_coef,
                               std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > , 
                               std::map <std::vector<ConstVariablePtr>, ConstVariablePtr> , 
                               std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_obj_coef,
                               std::map <ConstVariablePtr, ConstVariablePair> , 
                               std::map <ConstVariablePair, ConstVariablePtr> , 
                               std::map <ConstVariablePair, std::vector<double> > blterms_cons_coef,
                               std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > , 
                               std::map <std::vector<ConstVariablePtr>, ConstVariablePtr> , 
                               std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_cons_coef,
                               std::vector <std::vector<ConstVariablePtr> > &groups,                               
                               int groupStrategy)
{

  UInt gs = env_->getOptions()->findInt("ml_max_group_size")->getValue();
  

  // WTF is this: how is randomization used?  I am turning it off for now
  // srand ( time(NULL) );

  std::vector <ConstVariablePtr> vars;
  getMultilinearVariables(blterms, mlterms, vars);

  if(groupStrategy == 1)  {
    groupUsingDensest(gs, blterms_coef, mlterms_coef, vars, groups);
  }
  
  if(groupStrategy == 2)  {
    groupUsingDensest(gs, blterms_obj_coef, mlterms_obj_coef, vars, groups);
  }
  
  if(groupStrategy == 3) {
    groupUsingDensest(gs, blterms_obj_coef, mlterms_obj_coef, vars, groups);
    groupUsingDensest(gs, blterms_cons_coef, mlterms_cons_coef, vars, groups);
  }

  if(groupStrategy == 4) {
    groupUsingDensest(gs, blterms_coef, mlterms_coef, vars, groups);
    double aug_factor = env_->getOptions()->findDouble("ml_cover_augmentation_factor")->getValue();
    UInt maxGroupsCnt = (UInt) (groups.size()*aug_factor);

#if defined DEBUG_MULTILINEAR_HANDLER
    std::cout << "Second round, will augment from " << groups.size() << " to " << maxGroupsCnt << " groups." 
              << std::endl;
#endif

    if(groups.size() < maxGroupsCnt) {
      groupUsingDensest2ND(gs, blterms_coef, mlterms_coef, vars, groups, maxGroupsCnt);      
    }
  }

  if(groupStrategy ==5) {
    // term-by-term relaxation
    groupTermByTerm(blterms, mlterms, groups);
  }


#if defined DEBUG_MULTILINEAR_HANDLER
  printf("\n\ngroups:\n");
  for(UInt i = 0; i< groups.size(); i++) {
    for(UInt j = 0; j< groups[i].size(); j++) {
      printf("%d\t", groups[i][j]->getId());
    }
    printf("\n");
  }
  printf("\nNumber of groups: %d\n", int(groups.size()));
#endif
  
  //calculate the number of times each term appears in the grouping
  std::vector<std::vector<ConstVariablePtr> > terms;
  getMultilinearTerms(blterms, mlterms, gs, terms);
  std::vector<int>termRep;
  countTermsAppearance(terms, groups, termRep);
  
#if defined(DEBUG_MULTILINEAR_HANDLER)
  // Count the number of terms with k repetitions for k = 1,2,...
  // -- find the max rep
  int maxRep = 0;
  for(UInt i = 0; i<termRep.size(); i++) {
    if(termRep[i]>maxRep)
      maxRep = termRep[i];
  }

  int *kRepCnt = new int[maxRep+1];
  for(int i=0; i<maxRep+1; i++) 
    kRepCnt[i] = 0;
  for(UInt j=0; j<termRep.size(); j++) {
    kRepCnt[termRep[j]]++;
  }

  printf("maxRep: %d\n", maxRep);
  printf("\nrep count:\n\n");
  for(int i=0; i<maxRep+1; i++) {
    printf("%d\n", kRepCnt[i]);
  }
  delete [] kRepCnt;
#endif

}


void
MultilinearHandler::allExtreme(std::vector<int> &S, std::vector<double> &lb,
                               std::vector<double> &ub, std::vector<std::vector<double> > &E)
{
  std::vector<double> val(S.size());
  UInt ix = 0;
  visit(S, lb, ub, ix, val, E);

}

void
MultilinearHandler::visit(std::vector<int> &S, std::vector<double> &lb,
                       std::vector<double> &ub, UInt ix, std::vector<double> &val,
                       std::vector<std::vector<double> > &E)
{
  if (ix == S.size()) {
    E.push_back(val);
    return ;
  }
  val[ix] = lb[S[ix]];
  visit(S,lb,ub,ix+1,val,E);

  val[ix] = ub[S[ix]];
  visit(S,lb,ub,ix+1,val,E);

}

void 
MultilinearHandler::makeBinaryVariablesPowers1(FunctionPtr &f, PolyFunPtr &p, QuadraticFunctionPtr &q, LinearFunctionPtr &l) {
  // In each term, if a binary variable takes a power bigger than 1
  // make that power 1
  // - first look at the quadratic part
  bool no_f_l = 0;
  if(!l) {
    no_f_l = 1;
    l = LinearFunctionPtr (new LinearFunction);
  }
  
  if(q) {
    for(VariablePairGroupConstIterator qit = q->begin(); qit != q->end(); ++qit) {
      ConstVariablePtr x1 = qit->first.first;
      ConstVariablePtr x2 = qit->first.second;
      if(x1 == x2 && (x1->getType() == Binary || (x1->getType() == Integer && x1->getUb() == 1 && x1->getLb() == 0))) {
	int qitCoef = qit->second;
	q->incTerm(qit->first, -1*qit->second);
	
	l->incTerm(qit->first.first, qitCoef);
	if(no_f_l) 
	  f->add(l);
      }
    }
  }
  
  // - now look at the polynomial part
  if(p) {
    for(MonomialConstIter mit = p->termsBegin(); mit != p->termsEnd(); ++mit) {
      for(VarIntMapConstIterator tit = (*mit)->termsBegin(); tit != (*mit)->termsEnd(); ++tit) {
	if((tit->first)->getType() == Binary || ((tit->first)->getType() == Integer && (tit->first)->getUb() == 1 && (tit->first)->getLb() == 0)) {
	  if(tit->second > 1)
	    (*mit)->multiply(1.0, tit->first, -1*(tit->second - 1));
	}
      }
    }
    
    p->removeLinear(l);
    p->removeQuadratic(q);
    if(no_f_l) 
      f->add(l);
  }  
}


std::string MultilinearHandler::getName() const
{
   return "MultilinearHandler (Handling Multilinear terms).";
}




//XXX Jeff says -- under here, we may not need it?



void
MultilinearHandler::groupUsingIntersection(UInt gs, double rho,
                                           std::vector<std::vector<ConstVariablePtr> > terms,
                                           std::vector <std::vector<ConstVariablePtr> > &groups)
{
  for(std::vector<std::vector<ConstVariablePtr> >::iterator allTermsIt = terms.begin(); allTermsIt != terms.end(); ++allTermsIt) {
    UInt termSize = (*allTermsIt).size();
    bool termInCurrGroups = 0;
    UInt initialGroupsSize = groups.size();
    bool *isGroupActive = new bool [groups.size()];
    UInt *allIntersections = new UInt [groups.size()];
    UInt maxIntersection = 0;
    int maxIntersectionInd = -1;          

    // Calculate the intersection size of the term with each group
    for(UInt i = 0; i < groups.size(); i++) {
      UInt intersectSize = 0;
      // Iterate over the variables of the term
      for(std::vector<ConstVariablePtr>::iterator termIt = (*allTermsIt).begin(); termIt != (*allTermsIt).end(); ++termIt) {
        // Iterate over the variables of the group
        for(std::vector<ConstVariablePtr>::iterator groupIt = groups[i].begin(); groupIt != groups[i].end(); ++groupIt) {
          if((*groupIt)->getId() == (*termIt)->getId()) {
            intersectSize++;
          }
        }
      }
      if(termSize - intersectSize > 0 && termSize - intersectSize <= gs - groups[i].size() && groups[i].size() < gs) 
        isGroupActive[i] = 1;
      else
        isGroupActive[i] = 0;
      
      if(termSize - intersectSize == 0) 
        termInCurrGroups = 1;
      
      allIntersections[i] = intersectSize;
    }  

    // If the term is not in the current groups
    if(!termInCurrGroups) {
      // go through the active groups and find the max intersection
      for(UInt i = 0; i < groups.size(); i++) {
        if(isGroupActive[i] && allIntersections[i] >= maxIntersection) {
          maxIntersection = allIntersections[i];
          maxIntersectionInd = i;
        }
      }
      
      // If maximum intersection is 0 then add a group that has the variables of the term
      if(maxIntersection == 0) {
        groups.push_back((*allTermsIt));
      }
      else {
        // Add the variables of the term to the group with max intersection
        for(std::vector<ConstVariablePtr>::iterator termIt = (*allTermsIt).begin(); termIt != (*allTermsIt).end(); ++termIt) {
          bool varInGroup = 0;
          for(UInt j = 0; j < groups[maxIntersectionInd].size(); j++) {
            if((*termIt)->getId() == groups[maxIntersectionInd][j]->getId()) {
              varInGroup = 1;
              break;
            }
          }
          if(!varInGroup) {
            groups[maxIntersectionInd].push_back(*termIt);
          }
        }
      }
    }
    
    // Now go through the groups and by a probablity depending on the intersection
    // size of the group, add the term to the group
    int extraAdd = 0;
    for(UInt i = 0; i < initialGroupsSize; i++) {
      double r = (double) rand()/RAND_MAX;
      //if(r < (allIntersections[i])/((*allTermsIt).size()) && i != maxIntersectionInd && isGroupActive[i]) {
      if(r < 1-pow(1-rho, allIntersections[i]) && int(i) != maxIntersectionInd && isGroupActive[i]) {
        extraAdd++;
        for(std::vector<ConstVariablePtr>::iterator termIt = (*allTermsIt).begin(); termIt != (*allTermsIt).end(); ++termIt) {
          bool varInGroup = 0;
          for(UInt j = 0; j < groups[i].size(); j++) {
            if((*termIt)->getId() == groups[i][j]->getId()) {
              varInGroup = 1;
              break;
            }
          }
          if(!varInGroup) {
            groups[i].push_back(*termIt);
          }
        }
      }
    }
    //if(!extraAdd)
    //groups.push_back((*allTermsIt));
  }
}











void
MultilinearHandler::groupUsingIntersection2(UInt gs, 
                                           std::vector<std::vector<ConstVariablePtr> > terms,
                                           std::vector <std::vector<ConstVariablePtr> > &groups)
{
  for(UInt termsIt = 0; termsIt < terms.size(); termsIt++) {
    printf("%d---%d\n", terms[termsIt][0]->getId(), terms[termsIt][1]->getId());
  }
  int thrshd = 3;//int(gs/2);
  for(UInt termsIt = 0; termsIt < terms.size(); termsIt++) {
    UInt termSize = 2;
    bool termInCurrGroups = 0;
    bool *isGroupActive = new bool [groups.size()];
    UInt *allIntersections = new UInt [groups.size()];
    UInt maxIntersection = 0;
    int maxIntersectionInd = -1;        
    std::vector< std::vector<int> > allGroupsIntersection;
    bool termIsAdded = 0;
    
    //PIPI
    printf("\nterm in hand: %d --- %d\n\n", terms[termsIt][0]->getId(), terms[termsIt][1]->getId());
    printf("groups: \n");
    for(UInt ii=0; ii<groups.size(); ii++) {
      for(UInt jj = 0; jj<groups[ii].size(); jj++) {
        printf("%d\t", groups[ii][jj]->getId());
      }
      printf("\n");
    }



    // Calculate the intersection size of the term with each group
    for(UInt i = 0; i < groups.size(); i++) {
      UInt intersectSize = 0;
      // Iterate over the variables of the term
      for(UInt tIt = 0; tIt < terms[termsIt].size(); tIt++) {
        // Iterate over the variables of the group
        for(std::vector<ConstVariablePtr>::iterator groupIt = groups[i].begin(); groupIt != groups[i].end(); ++groupIt) {
          if((*groupIt)->getId() == terms[termsIt][tIt]->getId()) {
            intersectSize++;
          }
        }
      }
      if(termSize - intersectSize > 0 && termSize - intersectSize <= gs - groups[i].size() && 
         groups[i].size() < gs && intersectSize > 0) 
        isGroupActive[i] = 1;
      else
        isGroupActive[i] = 0;
      
      if(termSize - intersectSize == 0) 
        termInCurrGroups = 1;
      
      allIntersections[i] = intersectSize;
    }  


    // If the term is not in the current groups
    if(!termInCurrGroups) {
      // go through the active groups and find the max intersection
      for(UInt i = 0; i < groups.size(); i++) {
        if(isGroupActive[i] && allIntersections[i] >= maxIntersection) {
          maxIntersection = allIntersections[i];
          maxIntersectionInd = i;
        }
      }
      
      // If maximum intersection is 0 then add a group that has the variables of the term
      if(maxIntersection == 0) {
        groups.push_back(terms[termsIt]);
        termIsAdded = 1;
      }
      else {
        // go through the active groups and if adding the term to an active group
        // won't make the intersection of the active group with other groups more
        // than the threshhold, then add the term to the active group
        ConstVariablePtr x1 = terms[termsIt][0];
        ConstVariablePtr x2 = terms[termsIt][1];
        for(UInt i = 0; i < groups.size(); i++) {
          if(isGroupActive[i]) {
            int maxIntersectionWithGroups = 0;
            bool isX1InGroup = 0;
            bool isX2InGroup = 0;
            for(UInt j=0; j<groups[i].size(); j++) {
              if(groups[i][j]->getId() == x1->getId())
                isX1InGroup = 1;
              if(groups[i][j]->getId() == x2->getId())
                isX2InGroup = 1;
            }
            
            //PIPI
            printf("\n%d\tisX1: %d\t isX2: %d\n", i, isX1InGroup, isX2InGroup);

            if(isX1InGroup + isX2InGroup == 1) {
              if(isX1InGroup) {
                for(UInt k=0; k<groups.size(); k++) {
                  if(k!=i) {
                    int intersectionWithGroup = 0;
                    for(UInt l=0; l<groups[k].size(); l++) {
                      for(UInt m=0; m<groups[i].size(); m++) {
                        if(groups[k][l]->getId() == groups[i][m]->getId()) {
                          intersectionWithGroup++;
                        }
                      }
                      if(x2->getId() == groups[k][l]->getId())
                        intersectionWithGroup++;
                    }
                    if(maxIntersectionWithGroups < intersectionWithGroup)
                      maxIntersectionWithGroups = intersectionWithGroup;
                  }
                }
              }
              if(isX2InGroup) {
                for(UInt k=0; k<groups.size(); k++) {
                  if(k!=i) {
                    int intersectionWithGroup = 0;
                    for(UInt l=0; l<groups[k].size(); l++) {
                      for(UInt m=0; m<groups[i].size(); m++) {
                        if(groups[k][l]->getId() == groups[i][m]->getId()) {
                          intersectionWithGroup++;
                        }
                      }
                      if(x1->getId() == groups[k][l]->getId())
                        intersectionWithGroup++;
                    }
                    if(maxIntersectionWithGroups < intersectionWithGroup)
                      maxIntersectionWithGroups = intersectionWithGroup;
                  }
                }
              }
            }
            // if max intersections with groups is less than the max intersection threshhold
            // then add the term to the active group
            if(maxIntersectionWithGroups > 0 && maxIntersectionWithGroups <= thrshd) {
              termIsAdded = 1;
              if(isX1InGroup)
                groups[i].push_back(x2);
              if(isX2InGroup)
                groups[i].push_back(x1);
              break;
            }
          }
        }
        printf("term is Added: %d\n", termIsAdded);
        if(termIsAdded == 0) {
          groups.push_back(terms[termsIt]);
          termIsAdded = 1;
        }
      }
    }
  }
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
