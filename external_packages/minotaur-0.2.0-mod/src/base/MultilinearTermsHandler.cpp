//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file MultilinearTermsHandler.cpp
 * \brief Handles Multilinear Terms.  Each term must appear as equal to some
 * auxiliary variable.
 * \author Jeff Linderoth
 */

#include <algorithm>
#include <cmath>
#include <sstream>

#include "Branch.h"
#include "BrVarCand.h"
#include "Constraint.h"
#include "Environment.h"
#include "Function.h"
#include "LinearFunction.h"
#include "LinMods.h"
#include "Logger.h"
#include "MultilinearTermsHandler.h"
#include "Objective.h"
#include "Option.h"
#include "Node.h"
#include "Problem.h"
#include "Relaxation.h"
#include "Solution.h"
#include "Variable.h"

using namespace Minotaur;
using namespace std;

#undef DEBUG_MULTILINEARTERMS_HANDLER
#define PRINT_MULTILINEARTERMS_HANDLER_STATS

MultilinearTermsHandler::MultilinearTermsHandler(EnvPtr env, ProblemPtr problem)
  : env_(env), problem_(problem)
{
  logger_  = (LoggerPtr) new Logger((LogLevel) 
                                    env->getOptions()->findInt("handler_log_level")->getValue());
  eTol_ = env->getOptions()->findDouble("ml_feastol")->getValue();
  maxGroupSize_ = (UInt) env_->getOptions()->findInt("ml_max_group_size")->getValue();
  augmentCoverFactor_ = env_->getOptions()->findDouble("ml_cover_augmentation_factor")->getValue();

  initialTermCoverSize_ = 0;

}


void MultilinearTermsHandler::addConstraint(ConstraintPtr newcon, ConstVariablePtr ovar,
                                            set<ConstVariablePtr> ivars)
{
  Handler::addConstraint(newcon);
  termsO_.insert(make_pair(ovar, ivars));

}


Branches MultilinearTermsHandler::getBranches(BrCandPtr cand, DoubleVector & x,
                                         RelaxationPtr, SolutionPoolPtr)
{
  BrVarCandPtr vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();
  double value = x[v->getIndex()];

  // can't branch on something that is at its bounds.
  if (!(value > v->getLb()+1e-8 && value < v->getUb()-1e-8)) {
    cerr << "Warning!  Branching on variable with bounds/value: [" << 
      v->getLb() << " , " << value << "  " << v->getUb() << " ]" << endl;
    //assert(value > v->getLb()+1e-8 && value < v->getUb()-1e-8);
  }

  Branches branches = (Branches) new BranchPtrVector();
  BranchPtr branch;
  branch = doBranch_(DownBranch, v, value);
  branches->push_back(branch);
  branch = doBranch_(UpBranch, v, value);
  branches->push_back(branch);


  logger_->msgStream(LogDebug2) << "branching on " << v->getName();
  logger_->msgStream(LogDebug2) << " <= " << value << " or " 
    << " >= " << value << endl;

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)  
  cout << "branching on " << v->getName();
  cout << " <= " << value << " or " << " >= " << value << endl;
#endif

  return branches;
}


void MultilinearTermsHandler::getBranchingCandidates(RelaxationPtr, 
                                                const DoubleVector &x,
                                                ModVector &,
                                                BrVarCandSet &cands,
                                                BrCandVector &, bool &is_inf)
{
  // Implementation notes:
  //  (1) You must insert into the set cands
  //  (2) Super naive implmentation:  Just pick the variable from infeasible
  //      terms whose total 'bound product' is largest

  std::set<ConstVariablePtr> candidates;

  for (ConstTermIterator it = termsR_.begin(); it != termsR_.end(); ++it) {
    ConstVariablePtr zt = it->first;
    SetOfVars const &jt = it->second;

    if (allVarsBinary_(jt)) continue;

    double zval = x[zt->getIndex()];

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      std::cout << "Relaxation term variable: ";
      zt->write(std::cout);
      std::cout << " Has LP value: " << zval << std::endl;
#endif

    double termval = 1.0;
    double largest_score = eTol_;
    ConstVariablePtr largest_score_var;
    for(SetOfVars::const_iterator jt_it = jt.begin(); jt_it != jt.end(); ++jt_it) {
      ConstVariablePtr termvar = *jt_it;
      termval *= x[termvar->getIndex()];
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      std::cout << "Term variable: ";
      termvar->write(std::cout);
      std::cout <<  "  has value: " << x[termvar->getIndex()] << std::endl;
#endif

      double score = (termvar->getUb() - x[termvar->getIndex()])*
        (x[termvar->getIndex()] - termvar->getLb());

      if (score > largest_score) {
        largest_score = score;
        largest_score_var = termvar;
      }
    }
    if (fabs(zval - termval) > eTol_) {
      candidates.insert(largest_score_var);
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      largest_score_var->write(std::cout);
      std::cout << " will be a candidate" << std::endl;      
#endif
    }
    
  }

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "Branching candidates are: " << std::endl;
  for(SetOfVars::const_iterator it = candidates.begin(); it != candidates.end(); ++it) {
    (*it)->write(std::cout);
  }
#endif  

  for(SetOfVars::const_iterator it = candidates.begin(); it != candidates.end(); ++it) {
    ConstVariablePtr v = *it;
    BrVarCandPtr br_can = (BrVarCandPtr) new BrVarCand(v, v->getIndex(), 0.5, 0.5);
    cands.insert(br_can);
  }
  is_inf = false;
}

ModificationPtr MultilinearTermsHandler::getBrMod(BrCandPtr cand, DoubleVector &xval, 
                                                  RelaxationPtr , BranchDirection dir)
{
  LinModsPtr linmods;

  //XXX Put (bool init) back in handle{x,z}def...

  BrVarCandPtr  vcand = boost::dynamic_pointer_cast <BrVarCand> (cand);
  VariablePtr v = vcand->getVar();
  
  double branching_value = xval[v->getIndex()];
  BoundType lu;
  VariableType vtype = v->getType();

  // Change bounds on the x var (called v here)
  if (dir == DownBranch) { 
    lu = Upper;    
    if (vtype != Continuous)  branching_value = floor(branching_value);
  }
  else {
    lu = Lower;
    if (vtype != Continuous)  branching_value = ceil(branching_value);
  }

  linmods = (LinModsPtr) new LinMods();

  VarBoundModPtr vmod = (VarBoundModPtr) new VarBoundMod(v, lu, branching_value);
  linmods->insert(vmod);
  
  
  // This chunk of code changes the
  // x_{V_g} = \sum_{k=1}^{2 |V_g|} \lambda_k^g \chi^{k,g} \forall g \in G


  for (UInt gix = 0; gix < groups_.size(); ++gix) {
    for(SetOfVars::const_iterator it = groups_[gix].begin(); it != groups_[gix].end(); ++it) {
      ConstVariablePtr xvar = *it;
      if (v != xvar) continue;

      LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
      lf->addTerm(xvar, -1.0);

      UInt pix = 0;
      for (std::set<SetOfVars>::iterator it2 = points_[gix].begin(); it2 != points_[gix].end(); ++it2) {
        VariablePtr lam = lambdavars_[gix][pix];
        double val = -INFINITY;

        bool atLower = varIsAtLowerBoundAtPoint_(v, *it2);
        bool atUpper = !atLower;
        
        if (lu == Upper && atUpper) val = branching_value;
        else if (lu == Lower && atLower) val = branching_value;
        else val = (atLower ? v->getLb() : v->getUb());

        lf->addTerm(lam, val);
        ++pix;
      }       
      FunctionPtr f = (FunctionPtr) new Function(lf);

      IntVarPtrPairConstraintMap::iterator pos;
      pos = xConMap_.find(IntVarPtrPair(gix, xvar));
      if (pos == xConMap_.end()) {
        assert(0);
      }
      ConstraintPtr c = pos->second;

      LinConModPtr lcmod = (LinConModPtr) new LinConMod(c, lf, 0.0, 0.0);
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      std::cout << "getBrMod().  Will change 'x =' constraint to have linear function ";
      lf->write(std::cout);        
      std::cout << std::endl;
#endif
      linmods->insert(lcmod);
    
    }
  }

  // This will change the z_t = sum \sum_{k=1}^{2|V_g} \lambda_k^g \chi^{k,g}.
  //  Probably not very efficient way to do this...
  for(ConstTermIterator it = termsR_.begin(); it != termsR_.end(); ++it) {
    SetOfVars const &jt = it->second;

    for (UInt gix = 0; gix < groups_.size(); ++gix) {
      SetOfVars &vg = groups_[gix];

      std::set<ConstVariablePtr>::iterator pos1;
      pos1 = jt.find(v);
      if (pos1 == jt.end()) continue; // J_t does not contain v, go to next group
      // J_t is not in V_g, go to next group
      if (! std::includes(vg.begin(), vg.end(), jt.begin(), jt.end())) continue;

      ConstVariablePtr zvar = it->first;
      LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
      lf->addTerm(zvar, -1.0);

      // Get ConstraintToChange
      IntVarPtrPairConstraintMap::iterator pos2;
      pos2 = zConMap_.find(IntVarPtrPair(gix, zvar));
      if (pos2 == zConMap_.end()) {
        assert(0);
      }
      ConstraintPtr c = pos2->second;

      UInt pix = 0;
      for (std::set<SetOfVars>::iterator it2 = points_[gix].begin(); 
           it2 != points_[gix].end(); ++it2) {

        double prodval = 1.0;
        VariablePtr lam = lambdavars_[gix][pix];

        // Compute new extreme point value for this lambda
        for(SetOfVars::const_iterator jt_it = jt.begin(); jt_it != jt.end(); ++jt_it) {
          ConstVariablePtr xvar = *jt_it;
          double val = 0.0;
          bool atLower = varIsAtLowerBoundAtPoint_(xvar, *it2);
          bool atUpper = !atLower;
            
          if (xvar == v) {
            if (lu == Upper && atUpper) val = branching_value;
            else if (lu == Lower && atLower) val = branching_value;
            else val = (atLower ? xvar->getLb() : xvar->getUb());
          }
          else {
            val = atLower ? xvar->getLb() : xvar->getUb();
          }            
          prodval *= val;
        }
        lf->addTerm(lam, prodval);
        ++pix;
      }

      // Add new mod
      LinConModPtr lcmod = (LinConModPtr) new LinConMod(c, lf, 0.0, 0.0);

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      std::cout << "getBrMod(): Will change 'zt = ' constraint to have linear function: ";
      lf->write(std::cout);  
      std::cout << std::endl;
#endif
      linmods->insert(lcmod);
    }
          
  }

  return linmods;

  return ModificationPtr();
}


bool
MultilinearTermsHandler::isFeasible(ConstSolutionPtr sol, RelaxationPtr ,
                                    bool &, double &)
{
  const double *x = sol->getPrimal();
  bool is_feas = true;

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "Checking feasibility: " << std::endl;
#endif

  for (ConstTermIterator it = termsR_.begin(); it != termsR_.end(); ++it) {
    ConstVariablePtr zt = it->first;
    SetOfVars const &jt = it->second;

    if (allVarsBinary_(jt)) continue;

    double zval = x[zt->getIndex()];
    double xval = 1.0;
    for(SetOfVars::const_iterator jt_it = jt.begin(); jt_it != jt.end(); ++jt_it) {
      xval *= x[(*jt_it)->getIndex()];
    }
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
    std::cout << "Term for variable: " << std::endl;
    zt->write(std::cout);
    std::cout << " has error: " << fabs(zval - xval) << std::endl;
#endif
    if (fabs(zval - xval) > eTol_) {
      is_feas = false;
      break;
    }
  }
  
  return(is_feas);
}

void
MultilinearTermsHandler::relaxInitInc(RelaxationPtr relaxation, bool *)
{

  //  General notes...
  // Use Id to order variables
  // x[v->getIndex()] returns solution.
  
  /*
    0) Create "terms" map using the relaxation variables, not the original
       problem variables
    1) Create groups based on container of terms
    2) Add all lambda variables based on groups.
       (Be sure to keep a container of which constraints are associated with
       which groups)
   */

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  cout << "In MultilinearTermHandler::relaxInitInc()" << endl << " Current Relaxation: " << endl;
  relaxation->write(cout);
  cout << "And current terms: " << endl;

  for(ConstTermIterator it = termsO_.begin(); it != termsO_.end(); ++it) {
    std::cout << "zvar: ";
    it->first->write(std::cout);
    std::cout << "Contains vars: " << std::endl;
    for(SetOfVars::const_iterator it2 = it->second.begin(); 
        it2 != it->second.end(); ++it2) {
      (*it2)->write(cout);
    }
  }
#endif  

  // Copy the original variable pointers into a data structure consisting of 
  //  relaxation variables   (It just makes likfe easier to keep things in 
  //  terms of relaxation variables).
  for (ConstTermIterator it = termsO_.begin(); it != termsO_.end(); ++it) {
    ConstVariablePtr ov = it->first;
    SetOfVars const &ovset = it->second;
    ConstVariablePtr rv = relaxation->getRelaxationVar(ov);
    SetOfVars rvset;
    for(SetOfVars::const_iterator it2 = ovset.begin(); it2 != ovset.end(); ++it2) {
      rvset.insert( relaxation->getRelaxationVar(*it2));
    }
    termsR_.insert(make_pair(rv, rvset));
  }
      
    
  // First we make the groups.

  makeGroups_();

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  cout << "After making groups: " << endl;

  for (ConstGroupIterator it = groups_.begin(); it != groups_.end(); ++it) {  
    cout << "Group of: " << endl;
    for(set<ConstVariablePtr>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
      (*it2)->write(cout);
    }
  }

  // Just checking groups now.  Stop here.
  //exit(1);
#endif

  /* This chunk of code creates the (powerset) representation 
   * of the extreme points in the groups, adds the lambda variables,
   * and adds the convexity constraints
   */

  lambdavars_.resize(groups_.size());
  for(UInt gix = 0; gix < groups_.size(); ++gix) {    

    LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();

    std::set<SetOfVars> p = powerset_(groups_[gix]);
    int pix = 0;
    for (std::set<SetOfVars>::iterator it2 = p.begin(); it2 != p.end(); ++it2) {
      std::string name;
      std::stringstream name_stream;
      name_stream << "lam_" << gix << "_" << pix;
      name = name_stream.str();
      
      VariablePtr lam;
      lam = relaxation->newVariable(0.0, 1.0, Continuous, name);
      lf->addTerm(lam, 1.0);
#if defined(DEBUG_MULTILINEARTERMS_HANDLER2)
      std::cout << "Adding: " << name << std::endl;
#endif

      lambdavars_[gix].push_back(lam);
      pix++;
    }    
    points_.push_back(p);
    FunctionPtr f = (FunctionPtr) new Function(lf);
    relaxation->newConstraint(f, 1.0, 1.0);
  }  


  ModVector mods;  // Not used in initialization
  handleXDefConstraints_(relaxation, relaxInit_Call, mods);
  handleZDefConstraints_(relaxation, relaxInit_Call, mods);


#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "In MultilinearTermHandler::relaxInitInc()" << std::endl << " Final Relaxation: " << std::endl;
  relaxation->write(std::cout);
#endif
    
}


void
MultilinearTermsHandler::relaxNodeInc(NodePtr node,
                                      RelaxationPtr relaxation, bool *isInfeasible)
{
  *isInfeasible = false;

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "MultilinearTermsHandler::relaxNodeInc.  Node: ";
  node->write(std::cout);
  std::cout << "Initial Relaxation: ";
  relaxation->write(std::cout);
#endif

  ModVector mods;
  handleXDefConstraints_(relaxation, relaxNodeInc_Call, mods);
  handleZDefConstraints_(relaxation, relaxNodeInc_Call, mods);
  
  ModificationConstIterator it;
  for (it = mods.begin(); it != mods.end(); ++it) {
    node->addRMod(*it);
  }

#if defined(DEBUG_MULTILINEARTERMS_HANDLER) 
  std::cout << "MultilinearTermsHandler::relaxNodeInc.  Final relaxation:" << std::endl;
  relaxation->write(std::cout);
#endif

}


void
MultilinearTermsHandler::handleXDefConstraints_(RelaxationPtr relaxation, HandleCallingFunction wherefrom, ModVector &mods)
{
  for (UInt gix = 0; gix < groups_.size(); ++gix) {
    for(SetOfVars::const_iterator it = groups_[gix].begin(); it != groups_[gix].end(); ++it) {
      ConstVariablePtr xvar = *it;
      LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
      lf->addTerm(xvar, -1.0);
      
      int pix = 0;
      for (std::set<SetOfVars>::iterator it2 = points_[gix].begin(); it2 != points_[gix].end(); ++it2) {
        VariablePtr lam = lambdavars_[gix][pix];
        double val = varIsAtLowerBoundAtPoint_(xvar, *it2) ? xvar->getLb() : xvar->getUb();
        lf->addTerm(lam, val);
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
        std::cout << xvar->getName() << ", lam: " << gix << "," << pix << " value is: " 
                  << val << std::endl;
#endif        
        ++pix;
      }      
      FunctionPtr f = (FunctionPtr) new Function(lf);

      if (wherefrom == relaxInit_Call) {
        ConstraintPtr c = relaxation->newConstraint(f, 0.0, 0.0);
        xConMap_.insert(std::make_pair(IntVarPtrPair(gix, xvar), c));
      }
      else { 
        IntVarPtrPairConstraintMap::iterator pos;
        pos = xConMap_.find(IntVarPtrPair(gix, xvar));
        if (pos == xConMap_.end()) {
          assert(0);
        }
        ConstraintPtr c = pos->second;
        //XXX Here you should just check if the constraint really was going to
        //change and do nothing if it doesn't...  (will be faster).
        if (wherefrom == relaxNodeInc_Call) {
          relaxation->changeConstraint(c, lf, 0.0, 0.0);
        }
        else {
          assert(0);
        }
        LinConModPtr lcmod = (LinConModPtr) new LinConMod(c, lf, 0.0, 0.0); 
        mods.push_back(lcmod);
      }      
    }
  }  
  
}


// This chunk of code adds
// z_t = \sum_{k=1}^{2 |V_g|} \lambda_k^g, \Prod_{j \in J_t} \chi_j^{g,k}
//   \forall J_t \subseteq V_g
void
MultilinearTermsHandler::handleZDefConstraints_(RelaxationPtr relaxation, HandleCallingFunction wherefrom, ModVector &mods)
{
  for(ConstTermIterator it = termsR_.begin(); it != termsR_.end(); ++it) {
    ConstVariablePtr zt = it->first;
    SetOfVars const &jt = it->second;

    for (UInt gix = 0; gix < groups_.size(); ++gix) {
      SetOfVars &vg = groups_[gix];
       
      if (std::includes(vg.begin(), vg.end(), jt.begin(), jt.end())) { 
        // jt is a subset of vg, add constraint
        LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
        lf->addTerm(zt, -1.0);
        int pix = 0;

        for (std::set<SetOfVars>::iterator it2 = points_[gix].begin(); 
             it2 != points_[gix].end(); ++it2) {
          double prodval = 1.0;
          VariablePtr lam = lambdavars_[gix][pix];

          for(SetOfVars::const_iterator jt_it = jt.begin(); jt_it != jt.end(); ++jt_it) {
            ConstVariablePtr xvar = *jt_it;
            double tmp = varIsAtLowerBoundAtPoint_(xvar, *it2) ? xvar->getLb() : xvar->getUb();
            prodval *= tmp;
          }

          lf->addTerm(lam, prodval);
          ++pix;
        }        
        FunctionPtr f = (FunctionPtr) new Function(lf);

        if (wherefrom == relaxInit_Call) {
          ConstraintPtr c = relaxation->newConstraint(f, 0.0, 0.0);
          zConMap_.insert(std::make_pair(IntVarPtrPair(gix, zt), c));
        }
        else { 
          IntVarPtrPairConstraintMap::iterator pos;
          pos = zConMap_.find(IntVarPtrPair(gix, zt));
          if (pos == zConMap_.end()) {
            assert(0);
          }
          ConstraintPtr c = pos->second;
          //XXX Here you should just check if the constraint really was going to
          //change and do nothing if it doesn't...  (will be faster).
          if (wherefrom == relaxNodeInc_Call) {
            relaxation->changeConstraint(c, lf, 0.0, 0.0);
          }
          LinConModPtr lcmod = (LinConModPtr) new LinConMod(c, lf, 0.0, 0.0); 
          mods.push_back(lcmod);
        }
      }
      
    }
  } 


}

bool MultilinearTermsHandler::allVarsBinary_(SetOfVars const &s) const
{
  bool allbinary = true;
  for(SetOfVars::const_iterator it = s.begin(); it != s.end(); ++it) {
    VariableType vt = (*it)->getType();
    if ( vt == Continuous || vt == Integer || vt == ImplInt) {
      allbinary = false;
      break;
    }
  }
  return allbinary;
}

BranchPtr
MultilinearTermsHandler::doBranch_(BranchDirection UpOrDown, ConstVariablePtr v, 
                                   double bvalue)
{
  BranchPtr branch;
  BoundType lu;
  VariableType vtype = v->getType();

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "Branching: " << (UpOrDown == DownBranch ? "Down" : "Up")
            << " at value: " << bvalue << " on: " << std::endl;
  v->write(std::cout);
#endif

  branch = (BranchPtr) new Branch();

  double branching_value = bvalue;

  // Change bounds on the x var (called v here)
  if (UpOrDown == DownBranch) { 
    lu = Upper;    
    if (vtype != Continuous)  branching_value = floor(bvalue);
  }
  else {
    lu = Lower;
    if (vtype != Continuous)  branching_value = ceil(bvalue);
  }
 
  VarBoundModPtr vmod = (VarBoundModPtr) new VarBoundMod(v, lu, branching_value);
  assert(!"check whether this needs to be addRMod instead");
  branch->addPMod(vmod);

  branch->setActivity(0.5);// TODO: set this correctly
  return branch;

}

bool MultilinearTermsHandler::varIsAtLowerBoundAtPoint_(ConstVariablePtr &x,
                                                        SetOfVars const &p)
{
  SetOfVars::iterator it = p.find(x);
  return(it == p.end());
}

std::set<std::set<ConstVariablePtr> >
MultilinearTermsHandler::powerset_(SetOfVars const &s)
{
  typedef SetOfVars::const_iterator set_iter;
  typedef std::vector<set_iter> vec;

  struct local
  {
    static ConstVariablePtr dereference(set_iter v) { return *v; }
  };
 
  std::set<SetOfVars> result;
 
  vec elements;
  do
  {
    SetOfVars tmp;
    std::transform(elements.begin(), elements.end(),
                   std::inserter(tmp, tmp.end()),
                   local::dereference);
    result.insert(tmp);
    if (!elements.empty() && ++elements.back() == s.end())
    {
      elements.pop_back();
    }
    else
    {
      set_iter iter;
      if (elements.empty())
      {
        iter = s.begin();
      }
      else
      {
        iter = elements.back();
        ++iter;
      }
      for (; iter != s.end(); ++iter)
      {
        elements.push_back(iter);
      }
    }
  } while (!elements.empty());
 
  return result;
  
}


void 
MultilinearTermsHandler::makeGroups_()
{
  std::string s = env_->getOptions()->findString("ml_group_strategy")->getValue();

  if (s == "CONV") {
    std::cout << "Grouping: convex hull" << std::endl;
    SetOfVars vars;
    for(ConstTermIterator it = termsR_.begin(); it != termsR_.end(); ++it) {
      vars.insert(it->second.begin(),it->second.end());
    }
    groups_.push_back(vars);
  }  
  else if (s == "TT") {
    std::cout << "Grouping: term by term" << std::endl;
    for(ConstTermIterator it = termsR_.begin(); it != termsR_.end(); ++it) {
      groups_.push_back(it->second);
    }
  }
  else if (s == "TC") {
    std::cout << "Grouping: term cover" << std::endl;


    H_ = (HypergraphPtr) new Hypergraph(problem_);
    H_->create(termsR_);

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
    std::cout << "Hypergraph is: " << std::endl;
    H_->write(std::cout);
#endif    

    greedyDenseHeuristic_();

    // This is a stupid one
    //weightedDegreeHeuristic_();

  }
  else {
    assert(0);
  }

  // Always remove subsets.  (I don't think we need this, and it is slow. 
  //  May need it for random or something...)
  //removeSubsetsFromGroups_();

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "Final term cover has size: " << groups_.size() << std::endl;
#endif

}

void
MultilinearTermsHandler::addEdgeToGroups_(const SetOfVars &e, bool phaseOne)
{

  // In 'Phase 1', we do not create a new group if the edge is implied?
  bool room_for_edge = false;
  bool need_to_add_edge = true;
  int gixadd = -1;
  int ending_ix = phaseOne ? 0 : initialTermCoverSize_;
  int gix = 0;

  for(gix = ((int) groups_.size()) - 1; gix >= ending_ix; gix--) {

    const SetOfVars &g = groups_[gix];
    bool e_subsetof_g = edgeIsContainedInGroup_(e,g);
    if (phaseOne && e_subsetof_g) {
      need_to_add_edge = false;
      break;
    }
    if (!e_subsetof_g && edgeWillFitInGroup_(e,g)) {
      room_for_edge = true;
      gixadd = gix;
      break;
    }
  }

  if (need_to_add_edge) {
    if (room_for_edge) {
      groups_[gixadd].insert(e.begin(), e.end());
    }
    else {
      SetOfVars v;
      v.insert(e.begin(), e.end());
      groups_.push_back(v);
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      std::cout << "New group.  Size is now: " << groups_.size() << std::endl;
#endif
    }
  }
}

bool MultilinearTermsHandler::edgeIsContainedInGroup_(const SetOfVars &e,
                                                      const SetOfVars &g) const
{
  return (std::includes(g.begin(), g.end(), e.begin(), e.end()));
}

bool MultilinearTermsHandler::edgeWillFitInGroup_(const SetOfVars &e,
                                                  const SetOfVars &g) const
{
  UInt nUnique = 0;
  for(SetOfVars::const_iterator it = e.begin(); it != e.end(); ++it) {
    SetOfVars::iterator pos = g.find(*it);
    if (pos == g.end()) nUnique++;
  }

  bool fit = (g.size() + nUnique <= maxGroupSize_);
  return fit;

}

#if 0
MultilinearTermsHandler::WeightContainer::iterator 
MultilinearTermsHandler::findMaxWeight_()
{
  WeightContainer::iterator max_it = weights_.begin();
  double maxVal = -INFINITY;
  
  for(WeightContainer::iterator it = weights_.begin(); it != weights_.end(); ++it) {
    double val = it->second;
    if (val > maxVal) {
      maxVal = val;
      max_it = it;
    }
  }
  return max_it;
}
#endif


bool
MultilinearTermsHandler::varsAreGrouped_(SetOfVars const &termvars) const
{
  bool retval = false;
  for(UInt i = 0; i < groups_.size(); ++i) {
    if (std::includes(groups_[i].begin(), groups_[i].end(), termvars.begin(), termvars.end())) {
      retval = true;
      break;
    }
  }
  return retval;
}


void
MultilinearTermsHandler::removeSubsetsFromGroups_()
{
  GroupIterator it = groups_.begin();
  while(it != groups_.end()) {

    set<ConstVariablePtr> s1 = *it;
    bool removed = false;

    for(ConstGroupIterator it2 = groups_.begin(); it2 != groups_.end(); ++it2) {
      if (it == it2) continue;

      set<ConstVariablePtr> s2 = *it2;
      if (includes(s2.begin(), s2.end(), s1.begin(), s1.end())) {
        // s1 is subset of s2 erase it.
        // erase() returns iterator to next position
        it = groups_.erase(it);
        removed = true;
        break;
      }
    }
    if (!removed) ++it;
  }
}

#if 0
void 
MultilinearTermsHandler::randomCoverHeuristic_()
{
  bool positiveWeight = false;    
  SetOfVars e = H_->randomEdge(positiveWeight);
  while(positiveWeight) {
    double we = H_->getWeight(e);
    if (we <= 0.0) continue;  // If weight is not positive, we skip the edge
        
    addEdgeToGroups_(e, true);
    H_->setWeight(e, 0.0);
  }  

  initialTermCoverSize_ = groups_.size();
  
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "Initial term cover has size: " << groups_.size() << std::endl;
#endif
  
  // Do the second phase
  H_->resetWeights();
  UInt final_cover_size = (UInt) (augmentCoverFactor_ * initialTermCoverSize_);
  
  while(groups_.size() <= final_cover_size) {
    SetOfVars working_group;
    e = H_->randomEdge(positiveWeight);
    addEdgeToGroups_(e, true);
  }

  removeSubsetsFromGroups_();
  
}
#endif

void 
MultilinearTermsHandler::greedyDenseHeuristic_()
{

  // phase 1
  bool positiveWeight = false;    
  SetOfVars e = H_->heaviestEdge(positiveWeight);
  while(positiveWeight) {

    // Make a new group
    groups_.push_back(e);
    H_->setWeight(e, 0.0);

    // Add vertices to group (greedy).  If no incident vertex, we go to next group
    int gix = groups_.size()-1;
    bool incident_vertex = true;
    while(groups_[gix].size() < maxGroupSize_ && incident_vertex) {
      SetOfVars &g = groups_[gix];
      VariablePtr v = H_->heaviestIncidentVertex(g);
      if (v == 0) {
        incident_vertex = false;
      }
      else {
        groups_[gix].insert(v);      
        H_->adjustEdgeWeightsBetween(v,g,true);
      }
    }
    e = H_->heaviestEdge(positiveWeight);
  }  

  initialTermCoverSize_ = groups_.size();
  
#if defined(PRINT_MULTILINEARTERMS_HANDLER_STATS)
  std::cout << "Initial term cover has size: " << groups_.size() << std::endl;
#endif

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  cout << "Initial groups: " << endl;

  for (ConstGroupIterator it = groups_.begin(); it != groups_.end(); ++it) {  
    cout << "Group of: " << endl;
    for(set<ConstVariablePtr>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
      (*it2)->write(cout);
    }
  }
#endif

  // Do the second phase
  H_->resetWeights();
  UInt final_cover_size = (UInt) (augmentCoverFactor_ * initialTermCoverSize_);

  positiveWeight = true;

  while((groups_.size() <= final_cover_size) && positiveWeight) {
    e = H_->heaviestEdge(positiveWeight);
    double w = H_->getWeight(e);
    H_->setWeight(e, w/2.0);

    SetOfVars working_group;
    working_group.insert(e.begin(), e.end());

    bool incident_vertex = true;
    while( (working_group.size() < maxGroupSize_) && incident_vertex) {
      VariablePtr v = H_->heaviestIncidentVertex(working_group);
      if (v == 0) {
        incident_vertex = false;        
      }
      else {
        working_group.insert(v);      
        H_->adjustEdgeWeightsBetween(v,working_group,false);
      }
    }


    // Now see if it is a duplicate
    bool duplicate = false;
    for(UInt i = 0; i < groups_.size(); i++) {
      const SetOfVars &g = groups_[i];
      if (std::includes(g.begin(), g.end(), working_group.begin(), working_group.end())) {
#if  defined(DEBUG_MULTILINEARTERMS_HANDLER)  
        std::cout << "working_group is a duplicate" << std::endl;
        for(SetOfVars::const_iterator it = working_group.begin(); it != working_group.end(); ++it) {
           (*it)->write(std::cout);
        }
#endif
        duplicate = true;
        break;
      }
      
    }
    // Add it if not a duplicate
    if (!duplicate) {
      groups_.push_back(working_group);    
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)  
      std::cout << "Adding new element to groups.  Now have: " << groups_.size() << std::endl;
#endif
    }

    
    e = H_->heaviestEdge(positiveWeight);
  }  

#if defined(PRINT_MULTILINEARTERMS_HANDLER_STATS)
  std::cout << "Final term cover has size: " << groups_.size() << std::endl;
#endif
    
}
  

void 
MultilinearTermsHandler::weightedDegreeHeuristic_()
{
  // Cover all terms at least once
  bool positiveWeight = false;    
  do {
    VariablePtr heavyVar = H_->maxWeightedDegreeVertex(positiveWeight);      
    
    if (positiveWeight) {
      assert(heavyVar != 0);
      
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      std::cout << "Heaviest vertex is: ";
      heavyVar->write(std::cout);
#endif
      
      Hypergraph::ListOfSetOfVars edges = H_->incidentEdges(heavyVar);
      for(Hypergraph::ListOfSetOfVars::const_iterator e_it = edges.begin(); e_it != edges.end(); ++e_it) {
        const SetOfVars &e = *e_it;        
        double we = H_->getWeight(e);
        if (we <= 0.0) continue;  // If weight is not positive, we skip the edge
        
        addEdgeToGroups_(e, true);
        // Remove edge from H by setting weight to 0
        H_->setWeight(e, 0.0);
      }      
    }
  } while(positiveWeight);
  
  initialTermCoverSize_ = groups_.size();
  
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "Initial term cover has size: " << groups_.size() << std::endl;
#endif
  
  // Do the second phase
  H_->resetWeights();
  UInt final_cover_size = (UInt) (augmentCoverFactor_ * initialTermCoverSize_);
  
  while(groups_.size() <= final_cover_size) {
    VariablePtr heavyVar = H_->maxWeightedDegreeVertex(positiveWeight);      
    
    if (positiveWeight) {
      assert(heavyVar != 0);
      
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      std::cout << "Heaviest vertex is: ";
      heavyVar->write(std::cout);
#endif
      
      Hypergraph::ListOfSetOfVars edges = H_->incidentEdges(heavyVar);
      for(Hypergraph::ListOfSetOfVars::const_iterator e_it = edges.begin(); e_it != edges.end(); ++e_it) {
        const SetOfVars &e = *e_it;        
        double we = H_->getWeight(e);
        if (we <= 0.0) continue;  // If weight is not positive, we skip the edge
        
        addEdgeToGroups_(e, false);
        H_->setWeight(e, we/2.0);
      }      
    }
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
    std::cout << "Size of groups is now: " << groups_.size() << std::endl;
#endif
      
  }    
}

#if 0
void
MultilinearTermsHandler::setupWeights_()
{
  for(ConstTermIterator it = termsR_.begin(); it != termsR_.end(); ++it) {
    ConstVariablePtr zvar = it->first;
    double zweight = 0.0;
    for (ConstrSet::iterator it2 = zvar->consBegin(); it2 != zvar->consEnd(); ++it2) {
      const LinearFunctionPtr lf = (*it2)->getLinearFunction();
      double w = lf->getWeight(zvar);
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      zvar->write(std::cout);
      std::cout << "  has weight: " << w << " in constraint: ";
      (*it2)->write(std::cout);
#endif
      zweight += fabs(w);
    }
    weights_.insert(std::make_pair(zvar, zweight));
  }
  
}
#endif

void 
Hypergraph::adjustEdgeWeightsBetween(const VariablePtr v, const SetOfVars &g, 
                                     bool phaseOne)
{
  AdjListType::const_iterator pos = adjList_.find(v);
  assert(pos != adjList_.end());

  std::list<SetOfVars> e_list = pos->second;
  for(std::list<SetOfVars>::const_iterator e_it = pos->second.begin(); e_it != pos->second.end(); ++e_it) {
    SetOfVars new_e;
    new_e.insert(e_it->begin(), e_it->end());
    new_e.erase(v);
    if (std::includes(g.begin(), g.end(), new_e.begin(), new_e.end())) {
      if (phaseOne) {
        setWeight(*e_it, 0);
      }
      else {
        double w = getWeight(*e_it);
        setWeight(*e_it, w/2.0);        
      }
    }
  }
}


void 
Hypergraph::create(std::map<ConstVariablePtr, SetOfVars > const &terms)
{

  // Add the vertices
  std::map<ConstVariablePtr, SetOfVars >::const_iterator terms_it;
  for(terms_it = terms.begin(); terms_it != terms.end(); ++terms_it) {
    SetOfVars const &jt = terms_it->second;
    for(SetOfVars::iterator it = jt.begin(); it != jt.end(); ++it) {
      V_.insert(*it);
    }
  }

  // Add empty lists for the hyperedges
  for(SetOfVars::iterator v_it = V_.begin(); v_it != V_.end(); ++v_it) {
    std::list<SetOfVars> e_list;
    adjList_.insert(make_pair(*v_it, e_list));
  }

  // Now add the edges
  for(terms_it = terms.begin(); terms_it != terms.end(); ++terms_it) {
    SetOfVars const &jt = terms_it->second;
    E_.insert(jt);
    // Put it in the adjacency lists
    for (SetOfVars::const_iterator it = jt.begin(); it != jt.end(); ++it) {
      AdjListType::iterator pos = adjList_.find(*it);
      assert(pos != adjList_.end());
      std::list<SetOfVars> &e_list = pos->second;
      e_list.push_back(jt);
    }
    
    // Determine weight
    ConstVariablePtr zvar = terms_it->first;
    double zweight = 0.0;
    for (ConstrSet::iterator it2 = zvar->consBegin(); it2 != zvar->consEnd(); ++it2) {
      const LinearFunctionPtr lf = (*it2)->getLinearFunction();
      double w = lf->getWeight(zvar);
#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
      zvar->write(std::cout);
      std::cout << "  has weight: " << w << " in constraint: ";
      (*it2)->write(std::cout);
#endif
      zweight += fabs(w);
    }

    // Add objective weight
    const LinearFunctionPtr obj = problem_->getObjective()->getLinearFunction();
    if (obj != 0) {
      double w = obj->getWeight(zvar);
      zweight += fabs(w);
    }

    weights_.insert(make_pair(jt, zweight));
    originalWeights_.insert(make_pair(jt, zweight));
  }
}

SetOfVars 
Hypergraph::heaviestEdge(bool &positiveWeight) const
{
  positiveWeight = false;
  SetOfVarsDoubleMap::const_iterator max_pos;
  double maxWeight = 0.0;
  for (SetOfVarsDoubleMap::const_iterator it = weights_.begin(); it != weights_.end(); ++it) {
    double w = it->second;
    if (w > maxWeight) {
      positiveWeight = true;
      max_pos = it;
      maxWeight = w;
    }
  }
  SetOfVars e;
  if (positiveWeight) {
    e.insert(max_pos->first.begin(), max_pos->first.end());
  }
  return e;
}

VariablePtr
Hypergraph::heaviestIncidentVertex(const SetOfVars &g)
{

  VariablePtr bestv;
  double max_weight = 0.0;

  for(AdjListType::const_iterator it = adjList_.begin(); it != adjList_.end(); ++it) {
    VariablePtr v = it->first;    
    if (g.find(v) != g.end()) continue; // Skip it if in g

    double weight = 0.0;
    for (ListOfSetOfVars::const_iterator e_it = it->second.begin(); e_it != it->second.end(); ++e_it) {
      SetOfVars new_e;
      new_e.insert(e_it->begin(), e_it->end());
      new_e.erase(v);
      if (std::includes(g.begin(), g.end(), new_e.begin(), new_e.end())) {
        // This vertex is incident on the group with this edge
        weight += getWeight(*e_it);
      }
    }

    if (weight > max_weight) {
      max_weight = weight;
      bestv = v;
    }    
  }
  
  return bestv;
}

VariablePtr
Hypergraph::maxWeightedDegreeVertex(bool &positiveWeight) const
{
  VariablePtr heavyV;
  double max_deg = 0.0;
  positiveWeight = false;

  for(SetOfVars::iterator it = V_.begin(); it != V_.end(); ++it) {
    double w = weightedDegree(*it);    
    if (w > max_deg) {
      max_deg = w;
      positiveWeight = true;
      heavyV = *it;
    }
  }      
  return heavyV;
}

SetOfVars 
Hypergraph::randomEdge(bool &positiveWeight)
{
  positiveWeight = false;
  SetOfSetOfVars::iterator first_it = E_.begin();  
  int rix = (int) (drand48()*E_.size());
  std::advance(first_it, rix);

  SetOfSetOfVars::iterator e_it = first_it;
  do {
    const SetOfVars &e = *e_it;

    double w = getWeight(e);

    if (w > 0.0) {
      positiveWeight = true;
    }

    if (!positiveWeight) {
      e_it++;      
      if (e_it == E_.end()) e_it = E_.begin();
    }
    
  } while (!positiveWeight && e_it != first_it);
  
  return (*e_it);
}

void 
Hypergraph::resetWeights()
{

  for(SetOfVarsDoubleMap::iterator it = originalWeights_.begin(); it != originalWeights_.end(); ++it) {
    weights_[it->first] = it->second;    
  }

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  int ix = 0;
  for(SetOfVarsDoubleMap::iterator it = originalWeights_.begin(); it != originalWeights_.end(); ++it) {
    std::cout << "Resetting weight of edge " << (ix++) << " to: " << it->second << std::endl;
  }
#endif
}



double Hypergraph::getWeight(const SetOfVars &e) 
{
#if DEBUG
  // Just check that it exists (debugging)
  SetOfVarsDoubleMap::const_iterator pos = weights_.find(e);
  assert(pos != weights_.end());
#endif

  return weights_[e];

}

Hypergraph::ListOfSetOfVars
Hypergraph::incidentEdges(ConstVariablePtr v) const
{

  ListOfSetOfVars edges;
  AdjListType::const_iterator pos = adjList_.find(v);
  if (pos != adjList_.end()) {
    edges = pos->second;
  }          
  return edges;
}


void Hypergraph::setWeight(const SetOfVars &e, double w)
{
#if DEBUG
  // Just check that it exists (debugging)
  SetOfVarsDoubleMap::iterator pos = weights_.find(e);
  assert(pos != weights_.end());
#endif

  if (w < 0.0001) weights_[e] = 0.0;
  else weights_[e] = w;

#if defined(DEBUG_MULTILINEARTERMS_HANDLER)
  std::cout << "Setting weight of edge: " << std::endl;
  for(SetOfVars::const_iterator it = e.begin(); it != e.end(); ++it) {
    (*it)->write(std::cout);
  }
  std::cout << "to " << w << std::endl;
#endif  

}

double Hypergraph::weightedDegree(ConstVariablePtr v) const
{
  double val = 0.0;
  AdjListType::const_iterator pos = adjList_.find(v);

  if (pos != adjList_.end()) {
    std::list<SetOfVars> e_list = pos->second;
    for(std::list<SetOfVars>::const_iterator e_it = e_list.begin(); e_it != e_list.end(); ++e_it) {
      SetOfVarsDoubleMap::const_iterator it = weights_.find(*e_it);
      assert(it != weights_.end());
      val += it->second;
    }
  }
  return val;
}

void
Hypergraph::write(std::ostream &out) const
{
  for(AdjListType::const_iterator it = adjList_.begin(); it != adjList_.end(); ++it) {
    ConstVariablePtr v = it->first;
    std::cout << "Vertex: ";
    v->write(out);
    std::list<SetOfVars> e_list = it->second;    
    for(std::list<SetOfVars>::const_iterator e_it = e_list.begin(); e_it != e_list.end(); ++e_it) {
      std::cout << " has edge: " << std::endl;    
      for(SetOfVars::const_iterator v_it = e_it->begin(); v_it != e_it->end(); ++v_it) {
        (*v_it)->write(out);
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
