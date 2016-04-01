//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file Problem.cpp
 * \brief Define base class Problem.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Engine.h"
#include "Function.h"
#include "HessianOfLag.h"
#include "Jacobian.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "NonlinearFunction.h"
#include "Objective.h"
#include "Problem.h"
#include "ProblemSize.h"
#include "QuadraticFunction.h"
#include "SOS.h"
#include "Variable.h"

using namespace Minotaur;
const std::string Problem::me_ = "Problem: ";

Problem::Problem() 
: cons_(0), 
  consModed_(false),
  engine_(0),
  initialPt_(0), 
  nativeDer_(false),
  nextCId_(0),
  nextSId_(0),
  nextVId_(0),
  numDCons_(0),
  numDVars_(0),
  obj_(ObjectivePtr()), 
  size_(ProblemSizePtr()),
  vars_(0), 
  varsModed_(false)

{
  logger_ = (LoggerPtr) new Logger(LogInfo);
}


Problem::~Problem()
{
  clear();
  vars_.clear();
  cons_.clear();
  obj_.reset();
  sos1_.clear();
  sos2_.clear();

  if (initialPt_) {
    delete [] initialPt_;
    initialPt_ = 0;
  }
}


void Problem::addToObj(LinearFunctionPtr lf)
{
  assert(engine_ == 0 ||
      (!"Cannot change objective after loading problem to engine\n")); 
  if (obj_) {
    obj_->add_(lf);
  } else {
    assert(!"Cannot add lf to an empty objective!");
  }
  consModed_ = true;
}


void Problem::addToObj(double c)
{
  assert(engine_ == 0 ||
      (!"Cannot change objective after loading problem to engine\n")); 
  if (obj_) {
    obj_->add_(c);
  } else {
    assert(!"Cannot add c to an empty objective!");
  }
  consModed_ = true;
}


void Problem::addToCons(ConstraintPtr cons, double c) 
{
  cons->add_(c);
}


void Problem::calculateSize(bool shouldRedo)
{
  if (!size_) {
    shouldRedo=true;
    size_ = (ProblemSizePtr) new ProblemSize();
  }

  if (consModed_ || varsModed_) {
    shouldRedo = true;
  }

  if (shouldRedo) {
    size_->vars = vars_.size();
    size_->cons = cons_.size();
    size_->objs = (obj_)?1:0;

    countVarTypes_();
    countConsTypes_();
    countObjTypes_();
  }
  consModed_ = false;
  varsModed_ = false;
  return;
}


void Problem::changeBound(UInt index, BoundType lu, double new_val)
{

  assert(index < vars_.size() || 
      !"Problem::changeBound: index of variable exceeds no. of variables.");

  if (lu == Lower) {
    vars_[index]->setLb_(new_val);
  } else {
    vars_[index]->setUb_(new_val);
  }
  if (engine_) {
    engine_->changeBound(vars_[index], lu, new_val);
  }
}


void Problem::changeBound(UInt index, double new_lb, double new_ub)
{

  assert(index < vars_.size() || 
      !"Problem::changeBound: index of variable exceeds no. of variables.");

  vars_[index]->setLb_(new_lb);
  vars_[index]->setUb_(new_ub);
  if (engine_) {
    engine_->changeBound(vars_[index], new_lb, new_ub);
  }
}


void Problem::changeBound(VariablePtr var, BoundType lu, double new_val)
{

  assert(var == vars_[var->getIndex()] || 
      !"Problem: Bound of variable not in a problem can't be changed.");

  if (lu == Lower) {
    var->setLb_(new_val);
  } else {
    var->setUb_(new_val);
  }
  if (engine_) {
    engine_->changeBound(var, lu, new_val);
  }
}


void Problem::changeBound(VariablePtr var, double new_lb, double new_ub)
{

  assert(var == vars_[var->getIndex()] || 
      !"Problem: Bound of variable that is not in problem can't be changed.");

  var->setLb_(new_lb);
  var->setUb_(new_ub);
  if (engine_) {
    engine_->changeBound(var, new_lb, new_ub);
  }
}


void Problem::changeBound(ConstraintPtr con, double new_lb, double new_ub)
{

  assert(con == cons_[con->getIndex()] || 
      !"Problem: Bound of constraint that is not in problem can't be changed.");
  assert (engine_==0 || 
      (!"Cannot change constraint after loading problem to engine\n")); 

  con->setLb_(new_lb);
  con->setUb_(new_ub);
  consModed_ = true;
}


void Problem::changeBound(ConstraintPtr con, BoundType lu, double new_val)
{

  assert(con == cons_[con->getIndex()] || 
      !"Problem: Bound of constraint that is not in problem can't be changed.");

  if (engine_) {
    engine_->changeBound(con, lu, new_val);
  }
  if (lu == Lower) {
    con->setLb_(new_val);
  } else {
    con->setUb_(new_val);
  }
  consModed_ = true;
}


void Problem::changeConstraint(ConstraintPtr con, NonlinearFunctionPtr nlf)
{
  // simply replacing lf is sufficient to take care of jacobian and hessian as
  // well.

  FunctionPtr f = con->getFunction();

  assert(f);
  assert(con == getConstraint(con->getIndex()));

  // It is important to apply changes to engine first. Some engines use the
  // old constraint stored in problem to make changes.
  if (engine_) {
    engine_->changeConstraint(con, nlf);
  }


  for (VarSet::iterator vit=f->varsBegin(); vit!=f->varsEnd(); ++vit) {
    (*vit)->outOfConstraint_(con);
  }

  con->changeNlf_(nlf);

  f = con->getFunction();
  for (VarSet::iterator vit=f->varsBegin(); vit!=f->varsEnd(); ++vit) {
    (*vit)->inConstraint_(con);
  }
  consModed_ = true;
}


void Problem::changeConstraint(ConstraintPtr con, LinearFunctionPtr lf, 
                               double lb, double ub)
{
  // simply replacing lf is sufficient to take care of jacobian and hessian as
  // well.

  FunctionPtr f = con->getFunction();

  assert(f);
  assert(con == getConstraint(con->getIndex()));

  // It is important to apply changes to engine first. Some engines use the
  // old constraint stored in problem to make changes.
  if (engine_) {
    engine_->changeConstraint(con, lf, lb, ub);
  }


  for (VarSet::iterator vit=f->varsBegin(); vit!=f->varsEnd(); ++vit) {
    (*vit)->outOfConstraint_(con);
  }

  con->changeLf_(lf);
  con->setLb_(lb);
  con->setUb_(ub);

  f = con->getFunction();
  for (VarSet::iterator vit=f->varsBegin(); vit!=f->varsEnd(); ++vit) {
    (*vit)->inConstraint_(con);
  }
  consModed_ = true;
}


void Problem::changeObj(FunctionPtr f, double cb)
{
  int err = 0;
  FunctionPtr f2 = f ? f->cloneWithVars(vars_.begin(), &err) : 
    (FunctionPtr) new Function();
  std::string name = (obj_) ? obj_->getName() : "obj";
  assert (err==0);
  if (engine_) {
    engine_->changeObj(f2, cb);
  }
  obj_ = (ObjectivePtr) new Objective(f2, cb, Minimize, name);
  consModed_ = true;
}


int Problem::checkConVars() const
{
  ConstraintPtr c;
  VariablePtr v;
  int err = 0;

  for (ConstraintConstIterator citer=cons_.begin(); citer!=cons_.end();
       ++citer) {
    c = *citer;
    for (VarSet::iterator vit=c->getFunction()->varsBegin(); 
         vit!=c->getFunction()->varsEnd(); ++vit) {
      v = *vit;
      if (v != vars_[v->getIndex()]) {
        err = 1;
        logger_->msgStream(LogError) << me_ << "variable " << v->getName()
                                     << " in constraint " << c->getName()
                                     << " is not a variable of this problem."
                                     << std::endl;
        c->write(logger_->msgStream(LogError));
        logger_->msgStream(LogError) << std::endl;
      }
    }
  }
  return err;
}


void Problem::clear() 
{
  VariableIterator viter;
  ConstraintIterator citer;
  SOSIterator siter;
  VariablePtr v;
  ConstraintPtr c;

  for (viter=vars_.begin(); viter!=vars_.end(); viter++) {
    v = *viter;
    v->clearConstraints_();
  }

  for (citer=cons_.begin(); citer!=cons_.end(); ++citer) {
    c = *citer;
    c.reset();
  }
  for (siter=sos1_.begin(); siter!=sos1_.end(); siter++) {
    delete *siter;
    *siter = 0;
  }
  for (siter=sos2_.begin(); siter!=sos2_.end(); siter++) {
    delete *siter;
    *siter = 0;
  }
  if (engine_) {
    engine_->clear();
    engine_ = 0;
  }
  consModed_ = true;
  varsModed_ = true;
}


// Does not clone Jacobian and Hessian yet.
ProblemPtr Problem::clone() const
{
  ConstraintPtr c;
  ConstConstraintPtr cc;
  VariablePtr v;
  ConstVariablePtr cv;
  FunctionPtr f;
  ObjectivePtr oPtr;
  VariableConstIterator vit0;
  int err = 0;
  VarVector vvec;

  ProblemPtr clonePtr = (ProblemPtr) new Problem();

  // Copy the variables.
  for (VariableConstIterator it=vars_.begin(); it!=vars_.end(); ++it) {
    cv = *it;
    v = clonePtr->newVariable(cv->getLb(), cv->getUb(), cv->getType(),
                              cv->getName(), cv->getSrcType());
    v->setState_(cv->getState());
    v->setSrcType(cv->getSrcType());
    v->setFunType_(cv->getFunType());
    v->setId_(cv->getId());
  }
  
  vit0 = clonePtr->varsBegin();
  // add constraints
  for (ConstraintConstIterator it=cons_.begin(); it!=cons_.end(); ++it) {
    cc = *it;
    // clone the function.
    f = cc->getFunction()->cloneWithVars(vit0, &err);
    assert(err==0);
    c = clonePtr->newConstraint(f, cc->getLb(), cc->getUb(), cc->getName());
    c->setId_(cc->getId());
    c->setState_(cc->getState());
  }    

  // copy SOS1 constraints
  vvec.clear();
  for (SOSConstIterator it=sos1_.begin(); it!=sos1_.end(); ++it) {
    for (VariableConstIterator it2 = (*it)->varsBegin();
         it2!=(*it)->varsEnd(); ++it2) {
      vvec.push_back(*it2);
    }
    clonePtr->newSOS((*it)->getNz(), (*it)->getType(), (*it)->getWeights(),
                     vvec, (*it)->getPriority());
    vvec.clear();
  }

  // add objective
  oPtr = getObjective();
  if (oPtr) {
    f = oPtr->getFunction()->cloneWithVars(vit0, &err);
    assert(err==0);
    clonePtr->newObjective(f, oPtr->getConstant(),
                           oPtr->getObjectiveType(), oPtr->getName()); 
  } 

  // Now clone everything else...
  if (initialPt_) {
    clonePtr->initialPt_= new double[vars_.size()];
    std::copy(initialPt_, initialPt_+vars_.size(), clonePtr->initialPt_);
  }

  clonePtr->jacobian_  = JacobianPtr(); // NULL.
  clonePtr->nextCId_   = nextCId_;
  clonePtr->nextSId_   = nextSId_;
  clonePtr->nextVId_   = nextVId_;
  clonePtr->hessian_   = HessianOfLagPtr(); // NULL.
  clonePtr->logger_    = (LoggerPtr) new Logger(logger_->getMaxLevel());
  clonePtr->numDVars_  = numDVars_;	
  clonePtr->numDCons_  = numDCons_;	
  clonePtr->engine_    = 0;
  clonePtr->consModed_ = consModed_;
  clonePtr->varsModed_ = varsModed_;

  // clone size
  if (size_) {
    clonePtr->size_                   = (ProblemSizePtr) new ProblemSize();
    clonePtr->size_->vars             = size_->vars;
    clonePtr->size_->cons             = size_->cons;
    clonePtr->size_->objs             = size_->objs;
    clonePtr->size_->bins             = size_->bins;
    clonePtr->size_->fixed            = size_->fixed;
    clonePtr->size_->ints             = size_->ints;
    clonePtr->size_->conts            = size_->conts;
    clonePtr->size_->linCons          = size_->linCons;
    clonePtr->size_->SOS1Cons         = size_->SOS1Cons;
    clonePtr->size_->SOS2Cons         = size_->SOS2Cons;
    clonePtr->size_->bilinCons        = size_->bilinCons;
    clonePtr->size_->multilinCons     = size_->multilinCons;
    clonePtr->size_->quadCons         = size_->quadCons;
    clonePtr->size_->nonlinCons       = size_->nonlinCons;
    clonePtr->size_->consWithLin      = size_->consWithLin;
    clonePtr->size_->consWithBilin    = size_->consWithBilin;
    clonePtr->size_->consWithMultilin = size_->consWithMultilin;
    clonePtr->size_->consWithQuad     = size_->consWithQuad;
    clonePtr->size_->consWithNonlin   = size_->consWithNonlin;
    clonePtr->size_->linTerms         = size_->linTerms;
    clonePtr->size_->multiLinTerms    = size_->multiLinTerms;
    clonePtr->size_->quadTerms        = size_->quadTerms;
    clonePtr->size_->objLinTerms      = size_->objLinTerms;
    clonePtr->size_->objQuadTerms     = size_->objQuadTerms;
    clonePtr->size_->objType          = size_->objType;
  } else {
    clonePtr->size_ = ProblemSizePtr(); // NULL
  }
  clonePtr->nativeDer_ = nativeDer_; // NULL

  return clonePtr;
}


void Problem::countConsTypes_()
{
  ConstraintIterator citer;
  ConstraintPtr cPtr;
  LinearFunctionPtr lf;
  QuadraticFunctionPtr qf;
  UInt linCons=0, bilinCons=0, multilinCons=0, quadCons=0, nonlinCons=0;
  UInt consWithLin=0, consWithBilin=0, consWithMultilin=0, consWithQuad=0, 
       consWithNonlin=0;
  UInt linTerms=0, quadTerms=0;

  for (citer=cons_.begin(); citer!=cons_.end(); ++citer) {
    cPtr = *citer;
    switch (cPtr->getFunctionType()) {
     case Constant: // TODO: for now consider it linear
     case Linear:
       linCons++;
       break;
     case Bilinear:
       bilinCons++;
       break;
     case Multilinear:
       multilinCons++;
       break;
     case Quadratic:
       quadCons++;
       break;
     default:
       nonlinCons++;
       break;
    }
    lf = cPtr->getLinearFunction();
    if (lf) {
      consWithLin++;
      linTerms += lf->getNumTerms();
    }
    qf = cPtr->getQuadraticFunction();
    if (qf) {
      consWithQuad++;
      quadTerms += qf->getNumTerms();
    }
    if (cPtr->getNonlinearFunction()) {
      consWithNonlin++;
    }
  }

  size_->SOS1Cons     = sos1_.size();
  size_->SOS2Cons     = sos2_.size();
  size_->linCons      = linCons;
  size_->bilinCons    = bilinCons;
  size_->multilinCons = multilinCons;
  size_->quadCons     = quadCons;
  size_->nonlinCons   = nonlinCons;

  size_->consWithLin       = consWithLin;
  size_->consWithBilin     = consWithBilin;
  size_->consWithMultilin  = consWithMultilin;
  size_->consWithQuad      = consWithQuad;
  size_->consWithNonlin    = consWithNonlin;

  size_->linTerms = linTerms;
  size_->quadTerms = quadTerms;
  return;
}


void Problem::countObjTypes_()
{
  LinearFunctionPtr lf;
  QuadraticFunctionPtr qf;
  UInt linTerms=0, quadTerms=0;

  if (obj_) {
    size_->objType = obj_->getFunctionType();
    lf = obj_->getLinearFunction();
    if (lf) {
      linTerms += lf->getNumTerms();
    }
    qf = obj_->getQuadraticFunction();
    if (qf) {
      quadTerms += qf->getNumTerms();
    }
    if (obj_->getNonlinearFunction()) {
    }
  }

  size_->objLinTerms = linTerms;
  size_->objQuadTerms = quadTerms;
  return;
}


void Problem::countVarTypes_()
{
  VariableIterator viter;
  VariablePtr vPtr;
  UInt bins=0, ints=0, conts=0, fixed=0;

  for (viter=vars_.begin(); viter!=vars_.end(); viter++) {
    vPtr = *viter;
    switch (vPtr->getType()) {
     case (Binary):
       ++bins;
       break;
     case (Integer):
       ++ints;
       break;
     case (Continuous):
       ++conts;
       break;
     default:
       break;
    }
    if (fabs(vPtr->getUb() - vPtr->getLb()) < 1e-9) {
      ++fixed;
    }
  }
  size_->bins  = bins;
  size_->ints  = ints;
  size_->conts = conts;
  size_->fixed = fixed;
  findVarFunTypes_();
  return;
}


void Problem::delMarkedCons()
{
  if (numDCons_>0) {
    ConstraintPtr c;
    UInt i;
    std::vector<ConstraintPtr> copycons;
    std::vector<ConstraintPtr> delcons;

    for (ConstraintIterator it=cons_.begin(); it!=cons_.end(); ++it) {
      c = *it;
      if (c->getState() == DeletedCons) {
        delcons.push_back(c);
      } else {
        copycons.push_back(c);
      }
    }

    if (engine_) {
      engine_->removeCons(delcons);
    }

    for (ConstraintIterator it=delcons.begin(); it!=delcons.end(); ++it) {
      c = *it;
      for (VarSet::iterator vit=c->getFunction()->varsBegin(); 
           vit!=c->getFunction()->varsEnd(); ++vit) {
        (*vit)->outOfConstraint_(c);
      }
    }

    i = 0;
    for (ConstraintIterator it=copycons.begin(); it!=copycons.end(); ++it)
          {
      (*it)->setIndex_(i);
      ++i;
    }

    cons_ = copycons;
    consModed_ = true;
    numDCons_ = 0;
  }
}


void Problem::delMarkedVars()
{
  assert(engine_ == 0 ||
      (!"Cannot delete variables after loading problem to engine\n")); 
  if (numDVars_>0) {
    VariablePtr v;
    UInt i=0;
    std::vector<VariablePtr> copyvars;
    for (VariableIterator it=vars_.begin(); it!=vars_.end(); ++it) {
      v = *it;
      if (v->getState() == DeletedVar) {
        // std::cout << "deleting variable: ";
        // v->write(std::cout);
        for (ConstrSet::const_iterator cit=v->consBegin(); cit!=v->consEnd(); 
            ++cit) {
          (*cit)->delFixedVar_(v, v->getLb());
        }
        if (obj_) {
          obj_->delFixedVar_(v, v->getLb());
        }
      } else {
        v->setIndex_(i);
        ++i;
        copyvars.push_back(v);
      }
    }
    vars_ = copyvars;

    varsModed_ = true;
    numDVars_ = 0;
  }
}


ProblemType Problem::findType()
{
  calculateSize();

  if (size_->cons == size_->linCons && 
      (Constant == size_->objType || Linear == size_->objType)) { 
    return (size_->bins+size_->ints > 0) ? MILP : LP;
  } else if (size_->cons == size_->linCons && 
      (Quadratic == size_->objType || Bilinear == size_->objType)) { 
    return (size_->bins+size_->ints > 0) ? MIQP : QP;

  } else if (size_->cons == size_->linCons+size_->bilinCons+size_->quadCons &&
             (Quadratic == size_->objType || Bilinear == size_->objType)) {
    return (size_->bins+size_->ints > 0) ? MIQCQP : QCQP;
  } else if (isPolyp_()) {
    return (size_->bins+size_->ints > 0) ? MIPOLYP : POLYP;
  } else {
    return (size_->bins+size_->ints > 0) ? MINLP : NLP;
  }
}


void Problem::findVarFunTypes_()
{
  FunctionType ftype;
  FunctionPtr of = (obj_)?obj_->getFunction():FunctionPtr();
  FunctionPtr f;
  VariablePtr v;
  NonlinearFunctionPtr nlf;
  LinearFunctionPtr lf;
  QuadraticFunctionPtr qf;
  
#if 0
  for (VariableConstIterator it=vars_.begin(); it!=vars_.end(); ++it) {
    ftype = Constant;
    v = *it;
    for (ConstrSet::iterator cit=v->consBegin(); cit!=v->consEnd(); ++cit) {
      f = (*cit)->getFunction();
      if (f) {
       ftype = funcTypesAdd(ftype, f->getVarFunType(v));
      }
    }
    v->setFunType_(ftype);
  }

  if (of) {
    for (VarSet::iterator vit=of->varsBegin(); vit!=of->varsEnd(); ++vit) {
      v = *vit;
      ftype = funcTypesAdd(v->getFunType(), of->getVarFunType(v));
      v->setFunType_(ftype);
    }
  }
#else
  for (VariableConstIterator it=vars_.begin(); it!=vars_.end(); ++it) {
    (*it)->setFunType_(Constant);
  }
  for (ConstraintConstIterator it=cons_.begin(); it!=cons_.end(); ++it) {
    f = (*it)->getFunction();
    if (f) {
      lf = f->getLinearFunction();
      qf = f->getQuadraticFunction();
      nlf = f->getNonlinearFunction();
      if (lf) {
        for (VariableGroupConstIterator lit=lf->termsBegin();
             lit!=lf->termsEnd(); ++lit) {
          v = lit->first;
          ftype = funcTypesAdd(v->getFunType(), Linear);
          v->setFunType_(ftype);
        }
      }
      if (qf) {
        VarCountConstMap *vmap = qf->getVarMap();
        for (VarCountConstMap::const_iterator lit=vmap->begin();
             lit!=vmap->end(); ++lit) {
          v = lit->first;
          ftype = funcTypesAdd(v->getFunType(), Quadratic);
          v->setFunType_(ftype);
        }
      }
      if (nlf) {
        for(VariableSet::const_iterator lit=nlf->varsBegin();
            lit!=nlf->varsEnd(); ++lit) {
          v = *lit;
          v->setFunType_(Nonlinear);
        }
      }
    }
  }
  if (of) {
    for (VarSet::iterator vit=of->varsBegin(); vit!=of->varsEnd(); ++vit) {
      v = *vit;
      ftype = funcTypesAdd(v->getFunType(), of->getVarFunType(v));
      v->setFunType_(ftype);
    }
  }

#endif
}


ConstraintPtr Problem::getConstraint(UInt index) const
{
  return cons_[index];
}


HessianOfLagPtr Problem::getHessian() const
{ 
  return hessian_;
}


JacobianPtr Problem::getJacobian() const
{
  return jacobian_;
}


LoggerPtr Problem::getLogger() 
{
  return logger_;
}


UInt Problem::getNumHessNnzs() const
{
  if (hessian_) {
    return hessian_->getNumNz();
  } 
  return 0;
}


UInt Problem::getNumJacNnzs() const
{
  if (jacobian_) {
    return jacobian_->getNumNz();
  }
  return 0;
}


UInt Problem::getNumLinCons()
{
  return size_->linCons;
}


UInt Problem::getNumSOS1()
{
  return size_->SOS1Cons;
}


UInt Problem::getNumSOS2()
{
  return size_->SOS2Cons;
}


ObjectivePtr Problem::getObjective() const
{
  return obj_;
}


// evaluate the objective at a given x. x must be of dimension 'n', the
// number of variables in the problem
double Problem::getObjValue(const double *x, int *err) const
{
  if (obj_) {
    return obj_->eval(x, err);
  } else {
    return 0;
  }
}


ConstProblemSizePtr Problem::getSize() const
{
  return size_;
}


VariablePtr Problem::getVariable(UInt index) const
{
  return vars_[index];
}


bool Problem::hasNativeDer() const
{
  return nativeDer_;
}


bool Problem::isLinear()
{
  if (size_) {
    if (size_->cons == size_->linCons && 
        (Constant == size_->objType || Linear == size_->objType)) { 
      return true;
    } 
  }
  return false;
}


bool Problem::isMarkedDel(ConstConstraintPtr con)
{
  return (con->getState() == DeletedCons);
}


bool Problem::isMarkedDel(ConstVariablePtr var)
{
  return (var->getState() == DeletedVar);
}


bool Problem::isPolyp_()
{
  // assume that we have already check for linear, quadratic ...
  FunctionPtr f;
  if (obj_ && obj_->getFunction() && 
      (obj_->getFunction()->getType() == Nonlinear ||
       obj_->getFunction()->getType() == UnknownFunction)) {
    return false;
  }
  for (ConstraintIterator it=cons_.begin(); it!=cons_.end(); ++it) {
    f = (*it)->getFunction();
    if (f && (f->getType()==Nonlinear || f->getType()==UnknownFunction)) {
      return false;
    }
  }
  return false;
}


bool Problem::isQP()
{
  if (size_) {
    if (isLinear()) {
      return false;
    } else if ((size_->linCons == size_->cons) && 
        (Constant == size_->objType || Linear ==size_->objType || 
         Quadratic == size_->objType || Bilinear == size_->objType)) {
      return true;
    } 
  }
  return false;
}


bool Problem::isQuadratic()
{
  if (size_) {
    if (isLinear()) {
      return false;
    } else if ((size_->linCons + size_->quadCons + 
                size_->bilinCons == size_->cons) && 
               (Constant == size_->objType || Linear == size_->objType || 
                Quadratic == size_->objType || Bilinear == size_->objType)) {
      return true;
    } 
  }
  return false;
}


void Problem::markDelete(ConstraintPtr con)
{
  con->setState_(DeletedCons);
  ++numDCons_;
}


void Problem::markDelete(VariablePtr var)
{
  assert(engine_ == 0 ||
      (!"Cannot delete variables after loading problem to engine\n")); 
  var->setState_(DeletedVar);
  ++numDVars_;
}


void Problem::negateObj()
{
  if (engine_) {
    engine_->negateObj();
  }
  if (obj_) {
    obj_->negate_();
  }
  if (hessian_) {
    hessian_->negateObj();
  }
}


VariablePtr Problem::newBinaryVariable()
{
  assert(engine_ == 0 ||
      ("Cannot add variables after loading problem to engine\n")); 
  VariablePtr v;
  std::string name;
  std::stringstream name_stream;
  name_stream << "var" << vars_.size();
  name = name_stream.str();
  v = newVariable(0.0, 1.0, Binary, name);
  return v;
}


VariablePtr Problem::newBinaryVariable(std::string name)
{
  assert(engine_ == 0 ||
      (!"Cannot add variables after loading problem to engine\n")); 
  VariablePtr v = newVariable(0.0, 1.0, Binary, name);
  return v;
}


ConstraintPtr Problem::newConstraint(FunctionPtr f, double lb, double ub,
                                     std::string name)
{
  ConstraintPtr c = (ConstraintPtr) new Constraint(nextCId_, cons_.size(), f,
                                                   lb, ub, name);
  ++nextCId_;
  if (f) {
    for (VarSet::iterator vit=f->varsBegin(); vit!=f->varsEnd(); ++vit) {
      (*vit)->inConstraint_(c);
    }
  }
  cons_.push_back(c);
  if (engine_ != 0) {
    engine_->addConstraint(c);
  }
  consModed_ = true;
  return c;
}


ConstraintPtr Problem::newConstraint(FunctionPtr funPtr, double lb, double ub)
{
  // set a name and call newConstraint above.
  std::stringstream name_stream;
  std::string name;
  ConstraintPtr c;

  // build a name
  name_stream << "cons" << cons_.size();
  name = name_stream.str();

  // make a constraint,
  c = (ConstraintPtr) newConstraint(funPtr, lb, ub, name);
  if (engine_ != 0) {
    engine_->addConstraint(c);
  }

  return c;
}


ObjectivePtr Problem::newObjective(FunctionPtr f, double cb, 
                                   ObjectiveType otyp)
{
  assert(engine_ == 0 ||
      (!"Cannot add objective after loading problem to engine\n")); 
  // XXX: set name.
  std::string name;
  std::stringstream name_stream;
  ObjectivePtr o;
  name_stream << "obj";

  name = name_stream.str();
  o = newObjective(f, cb, otyp, name);
  return o;
}


ObjectivePtr Problem::newObjective(FunctionPtr f, double cb, 
                                   ObjectiveType otyp, std::string name)
{
  assert(engine_ == 0 ||
      (!"Cannot add objective after loading problem to engine\n")); 

  obj_ = (ObjectivePtr) new Objective(f, cb, otyp, name);
  consModed_ = true;
  return obj_;
}


SOSPtr Problem::newSOS(int n, SOSType type, const double *weights,
                       const VarVector &vars, int priority, std::string name)
{
  SOSPtr sos = new SOS(n, type, weights, vars, nextSId_, priority, name);
  ++nextSId_;
  if (SOS1 == type) {
    sos1_.push_back(sos);
  } else if (SOS2 == type) {
    sos2_.push_back(sos);
  }
  return sos;
}


SOSPtr Problem::newSOS(int n, SOSType type, const double *weights,
                       const VarVector &vars, int priority)
{
  std::string name;
  std::stringstream name_stream;
  if (SOS1==type) {
    name_stream <<  "sos1_" << (sos1_.size() + sos2_.size());
  } else if (SOS2==type) {
    name_stream <<  "sos2_" << (sos1_.size() + sos2_.size());
  } else {
    name_stream <<  "sos_" << (sos1_.size() + sos2_.size());
  }
  name = name_stream.str();
  SOSPtr sos = newSOS(n, type, weights, vars, priority, name);
  return sos;
}


VariablePtr Problem::newVariable(VarSrcType stype)
{
  assert(engine_ == 0 ||
      (!"Cannot add variables after loading problem to engine\n")); 
  VariablePtr v;
  std::string name;
  std::stringstream name_stream;
  name_stream <<  "var" << vars_.size();
  name = name_stream.str();
  v = newVariable(-INFINITY, INFINITY, Continuous, name, stype);
  return v;
}


VariablePtr Problem::newVariable(double lb, double ub, VariableType vtype,
                                 VarSrcType stype)
{
  assert(engine_ == 0 || 
      (!"Cannot add variables after loading problem to engine\n")); 
  VariablePtr v;
  std::string name;
  std::stringstream name_stream;
  name_stream <<  "var" << vars_.size();
  name = name_stream.str();
  v = newVariable(lb, ub, vtype, name, stype);
  return v;
} 



VariablePtr Problem::newVariable(double lb, double ub, VariableType vtype,
                                 std::string name, VarSrcType stype)
{
  assert(engine_ == 0 ||
      (!"Cannot add variables after loading problem to engine\n")); 
  VariablePtr v;
  v = (VariablePtr) new Variable(nextVId_, vars_.size(), lb, ub, vtype, name);
  v->setSrcType(stype);
  ++nextVId_;
  vars_.push_back(v);
  varsModed_ = true;
  return v;
}


void Problem::newVariables(VariableConstIterator v_begin, 
                           VariableConstIterator v_end, VarSrcType stype)
{
  assert(engine_ == 0 ||
      (!"Cannot add variables after loading problem to engine\n")); 
  VariableConstIterator v_iter;
  VariablePtr v;

  for (v_iter=v_begin; v_iter!=v_end; v_iter++) {
    v = newVariable((*v_iter)->getLb(), (*v_iter)->getUb(),
        (*v_iter)->getType(), (*v_iter)->getName(), stype);
  }
}


void Problem::prepareForSolve()
{
  bool reload = false;

  if (consModed_ || varsModed_) {
    reload = true;
  }
  calculateSize();
  if (nativeDer_ && (true == reload || !hessian_)) {
    setNativeDer();
  } 
}


void Problem::removeObjective()
{
  assert(engine_ == 0 ||
      (!"Cannot change objective after loading problem to engine\n")); 
  if (obj_) {
    obj_.reset();
  } 
  return;
}

QuadraticFunctionPtr Problem::removeQuadFromObj()
{
  assert(engine_ == 0 ||
      (!"Cannot change objective after loading problem to engine\n")); 
  if (obj_) {
    return obj_->removeQuadratic_();
  }
  consModed_ = true;
  return QuadraticFunctionPtr(); // NULL
}


void Problem::resetDer()
{
  jacobian_  = JacobianPtr(); // NULL.
  hessian_   = HessianOfLagPtr(); // NULL.
}


void Problem::reverseSense(ConstraintPtr cons) 
{
  cons->reverseSense_();
  consModed_ = true;
}


void Problem::setEngine(Engine* engine) 
{
  if (engine_) {
    engine_->clear();
  }
  engine_ = engine;
}


void Problem::setHessian(HessianOfLagPtr hessian) 
{ 
  hessian_ = hessian; 
}


void Problem::setIndex_(VariablePtr v, UInt i)
{
  v->setIndex_(i);
}


void Problem::setInitialPoint(const double *x) 
{
  // if x is null or if there are no variables, do nothing.
  if (!x || vars_.size() == 0) {
    return;
  }

  // if initial point hasnt been set before, allocate memory. otherwise just
  // use the old space.
  if (!initialPt_) {
    initialPt_ = new double[vars_.size()];
  }

  // copy
  std::copy(x, x+vars_.size(), initialPt_);
}


void Problem::setInitialPoint(const double *x, size_t k) 
{
  // if x is null or if there are no variables, do nothing.
  if (!x || vars_.size() == 0) {
    return;
  }

  // if initial point hasnt been set before, allocate memory. otherwise just
  // use the old space.
  if (!initialPt_) {
    initialPt_ = new double[vars_.size()];
  }

  // copy
  std::copy(x, x+k, initialPt_);
  std::fill(initialPt_+k, initialPt_+vars_.size(), 0.);
}


void Problem::setJacobian(JacobianPtr jacobian)
{
  jacobian_ = jacobian;
}


void  Problem::setLogger(LoggerPtr logger)
{ 
  logger_ = logger;
}


void Problem::setNativeDer()
{
  calculateSize();
  nativeDer_ = true;
  jacobian_ = (JacobianPtr) new Jacobian(cons_, vars_.size());
  hessian_ = (HessianOfLagPtr) new HessianOfLag(this);
}


void Problem::setVarType(VariablePtr var, VariableType type)
{
  assert(var == vars_[var->getIndex()] || 
      !"Problem: Type of variable that is not in problem can't be changed.");

  if (size_) {
    switch (var->getType()) {
    case (Binary):
    case (ImplBin):
      --(size_->bins);
      break;
    case (Integer):
    case (ImplInt):
      --(size_->ints);
      break;
    default:
      --(size_->conts);
    }

    switch (type) {
    case (Binary):
    case (ImplBin):
      ++(size_->bins);
      break;
    case (Integer):
    case (ImplInt):
      ++(size_->ints);
      break;
    default:
      ++(size_->conts);
    }
  }
  var->setType_(type);
  varsModed_ = true;
}


void Problem::subst(VariablePtr out, VariablePtr in, double rat)
{
  bool stayin;
  assert(engine_ == 0 ||
      (!"Cannot substitute variables after loading problem to engine\n")); 
  for (ConstrSet::const_iterator it=out->consBegin(); it!=out->consEnd(); 
      ++it) {
    (*it)->subst_(out, in, rat, &stayin);
    if (stayin) {
      in->inConstraint_(*it);
    } else {
      in->outOfConstraint_(*it);
    }
  }
  obj_->subst_(out, in, rat);
  consModed_ = varsModed_ = true;
}


void Problem::unsetEngine() 
{
  engine_ = 0;
}


void Problem::write(std::ostream &out, std::streamsize out_p) const 
{
  ConstraintConstIterator citer;
  ConstraintPtr cPtr;
  VariableConstIterator viter;
  std::streamsize old_p = out.precision();

  out.precision(out_p);
  if (size_) {
    writeSize(out);
  }

  for (viter = vars_.begin(); viter != vars_.end(); ++viter) {
    (*viter)->write(out);
  }

  if (obj_) {
    obj_->write(out);
    out << std::endl;
  }

  for (citer = cons_.begin(); citer != cons_.end(); ++citer) {
    (*citer)->write(out);
    out << std::endl;
  }
  out.precision(old_p);

  //for (citer = cons_.begin(); citer != cons_.end(); ++citer) {
  //  (*citer)->displayFunctionMap();
  //}

}


void Problem::writeSize(std::ostream &out) const
{
  out << "Problem size:" << std::endl
    << " Number of variables = " << size_->vars << std::endl
    << " Number of binary variables = " << size_->bins << std::endl
    << " Number of general integer variables = " << size_->ints << std::endl
    << " Number of continuous variables = " << size_->conts << std::endl
    << " Number of fixed variables = " << size_->fixed << std::endl
    << " Number of constraints = " << size_->cons << std::endl
    << " Number of linear constraints = " << size_->linCons << std::endl
    << " Number of SOS1 constraints = " << size_->SOS1Cons << std::endl
    << " Number of SOS2 constraints = " << size_->SOS2Cons << std::endl
    << " Number of bilinear constraints = " << size_->bilinCons << std::endl
    << " Number of multilinear constraints = " << size_->multilinCons 
    << std::endl
    << " Number of quadratic constraints = " << size_->quadCons << std::endl
    << " Number of nonlinear constraints = " << size_->nonlinCons << std::endl
    << " Number of constraints with linear terms = " << size_->consWithLin 
    << std::endl
    << " Number of constraints with bilinear terms = " << size_->consWithBilin 
    << std::endl
    << " Number of constraints with multilinear terms = " 
    << size_->consWithMultilin << std::endl
    << " Number of constraints with quadratic terms = " << size_->consWithQuad 
    << std::endl
    << " Number of linear terms in constraints = " << size_->linTerms 
    << std::endl
    << " Number of multilinear terms in constraints = " << size_->multiLinTerms 
    << std::endl
    << " Number of quadratic terms in constraints = " << size_->quadTerms 
    << std::endl
    << " Number of objectives = " << size_->objs << std::endl
    << " Number of linear terms in objective = " << size_->objLinTerms 
    << std::endl
    << " Number of quadratic terms in objective = " <<  size_->objQuadTerms
    << std::endl
    << " Type of objective = " <<  getFunctionTypeString(size_->objType)
    << std::endl;
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
