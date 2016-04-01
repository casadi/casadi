// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2010 - 2014 The MINOTAUR Team.
// 

/**
 * \file CutMan2.cpp
 * \brief Implement the methods of CutMan2 class. 
 * \author Mahdi Hamzeei, University of Wisconsin-Madison
 */

#include <cmath>
#include <iostream>

#include "MinotaurConfig.h"
#include "Cut.h"
#include "CutMan2.h"
#include "Constraint.h"
#include "Engine.h"
#include "Environment.h"
#include "Function.h"
#include "LinearFunction.h"
#include "Logger.h"
#include "Node.h"
#include "Problem.h"
#include "ProblemSize.h"
#include "Relaxation.h"
#include "Solution.h"
#include "SolutionPool.h"
#include "Timer.h"
#include "Types.h"
#include "Variable.h"

using namespace Minotaur;

typedef std::list<CutPtr>::const_iterator CCIter;
typedef std::vector<ConstraintPtr>::const_iterator ConstIter;

CutMan2::CutMan2()
  : env_(EnvPtr()),   // NULL
    hashVec_(0),
    p_(ProblemPtr()),  // NULL
    absTol_(5e-2),
    MaxInactiveInRel_(10000),
    PoolSize_(200),
    CtThrsh_(0),
    ctMngrtime_(0),
    PrntCntThrsh_(0),
    numCuts_(0)
{
  stats_ = new CutStat();
  logger_ = (LoggerPtr) new Logger(LogDebug2);
  stats_->numAddedCuts = 0;
  stats_->numDeletedCuts = 0;
  stats_->numPoolToRel = 0;
  stats_->numRelToPool = 0 ;
  stats_->numCuts= 0;

  ctmngrInfo_.t = ctMngrtime_;
  ctmngrInfo_.RelTr = MaxInactiveInRel_;
  ctmngrInfo_.PoolTr = PoolSize_;
  ctmngrInfo_.PrntActCnt = PrntCntThrsh_;
}

CutMan2::CutMan2(EnvPtr env, ProblemPtr p)
  : env_(env),
    p_(ProblemPtr()),
    absTol_(5e-2),
    MaxInactiveInRel_(100),
    PoolSize_(70),
    ctMngrtime_(0),
    PrntCntThrsh_(0)
{
  stats_ = new CutStat();
  timer_ = env_->getNewTimer();
  UInt n = p->getNumVars();
  hashVec_ = new double[n];
  CtThrsh_ = 6 * (n + p->getNumCons());

  for (UInt i = 0; i < n; i++) {
    hashVec_[i] = rand()/double (RAND_MAX);
  }

  updateTime_ = 0.0;
  checkTime_ = 0.0;
  processedTime_ = 0.0;
  branchedTime_ = 0.0;

  stats_->numAddedCuts = 0;
  stats_->numDeletedCuts = 0;
  stats_->numPoolToRel = 0;
  stats_->numRelToPool = 0 ;
  stats_->callNums = 0;
  stats_->PoolSize = 0;
  stats_->RelSize = 0;
  stats_->numCuts= 0;
  ctmngrInfo_.t = ctMngrtime_;
  ctmngrInfo_.PoolTr = PoolSize_;
  ctmngrInfo_.PrntActCnt = PrntCntThrsh_;
}

CutMan2::~CutMan2()
{
  if (stats_){
    writeStat();
    delete stats_;
  }
  if (hashVec_) {
    delete[] hashVec_;
  }
  env_.reset();
  pool_.clear();
  rel_.clear();
  p_.reset();
  delete timer_;
}

void CutMan2::updateRel(ConstSolutionPtr sol, ProblemPtr rel)
{
  if ( numCuts_ >= CtThrsh_){
    timer_->start();
    UInt m = rel->getNumCons();
    const double *y = new double[m];
    y = sol->getDualOfCons();
    int i;
    CutPtr cut;
    cutList temp;

    for (std::list<CutPtr>::iterator it = rel_.begin(); it != rel_.end();)
    {
      cut = *it;
      i =cut->getConstraint()->getId();
      if (y[i] > 1e-6 || y[i] < -1e-6)
      {
        ++(cut->getInfo()->cntSinceActive);
        cut->getInfo()->cntSinceActive = 0;
      } else {
        ++(cut->getInfo()->cntSinceActive);
      } 

      if ( (y[i] > -1e-6 && y[i] < 1e-6) && cut->getInfo()->parent_active_cnts <=PrntCntThrsh_ &&  
	   cut->getInfo()->cntSinceActive > MaxInactiveInRel_) 
      {
        temp.push_back(cut);
        it = rel_.erase(it);
      } else {
        ++it;
      }
    }  
    double al;

    if (stats_->numRelToPool != 0){
      al = 1.0 * stats_->numPoolToRel / stats_->numRelToPool;
      //am = 1.0 * stats_->RelSize / stats_->numCuts;
      if ( (al > 0.1 && MaxInactiveInRel_ < 250)){ // || (am < 1.0 &&  MaxInactiveInRel_ > 0 )){
        MaxInactiveInRel_++;
      } else if (MaxInactiveInRel_ > 0){
        MaxInactiveInRel_ = MaxInactiveInRel_-1;
      }
      ctmngrInfo_.RelTr = MaxInactiveInRel_;
    }
    if (temp.size() > 5){
      for (std::list<CutPtr>::iterator it = temp.begin(); 
  	   it != temp.end();++it){
	cut = *it;
        addToPool_(cut);
        cut->getInfo()->cntSinceViol = 0;
        rel->markDelete(cut->getConstraint());
        stats_->numRelToPool++;
      }
      rel->delMarkedCons();
    } else {
	rel_.merge(temp);        
    }
    double a1 = timer_->query();
    ctMngrtime_ += a1;
    updateTime_ += a1;
    timer_->stop();
  }
}

void CutMan2::updatePool(ProblemPtr rel, ConstSolutionPtr sol)
{
  if ( numCuts_ >= CtThrsh_){
    timer_->start();
    const double *x = sol->getPrimal();
    double viol;
    double score;
    int toRel = 0;
    int deleted = 0;
    CutPtr cut;
    for (std::list<CutPtr>::iterator it = pool_.begin(); it != pool_.end();)
    {
      cut = *it;
      cut->evalScore(x, &viol, &score);

      if (score < 1e-6){
        ++(cut->getInfo()->cntSinceViol);
        ++it;
      }else if (score < absTol_){ 
	++it;
      } else {    
	addToRel_(rel,cut,false);
	toRel++;
        stats_->numPoolToRel++;
  	it = pool_.erase(it);
      } 
    }
    stats_->numDeletedCuts += deleted;
    ctMngrtime_ += timer_->query();
    checkTime_ += timer_->query();
    timer_->stop(); 
  }
}

ConstraintPtr CutMan2::addCut(ProblemPtr rel,FunctionPtr fn, double lb, double ub, bool, bool neverDelete)
{
  CutPtr cut = (CutPtr) new Cut(rel,fn, lb, ub,neverDelete,false);
  addToRel_(rel,cut,true);
  stats_->numAddedCuts++;
  cut->getInfo()->inRel = true;
  stats_->numCuts++;
  return cut->getConstraint();
}

void CutMan2::addToRel_(ProblemPtr rel, CutPtr cut, bool newcut)
{
  rel_.push_back(cut);
  if (false == newcut)
    cut->applyToProblem(rel);

  ++(cut->getInfo()->numActive);
  cut->getInfo()->cntSinceActive = 0;
  cut->getInfo()->cntSinceViol = 0;
  cut->getInfo()->inRel = true;

}

void CutMan2::addToRel_(CutPtr cut)
{
  rel_.push_back(cut);
  ++(cut->getInfo()->numActive);
  cut->getInfo()->cntSinceActive = 0;
  cut->getInfo()->cntSinceViol = 0;
  cut->getInfo()->inRel = true;
}

void CutMan2::addToPool_(CutPtr cut)
{
/*
  std::list<CutPtr>::iterator it = pool_.begin();
  CutPtr c;
  c = *it;
  int n = c->getInfo()->parent_active_cnts;
  if (pool_.size() > PoolSize_ - 1){
    while ( n >  PrntCntThrsh_ && it != pool_.end() )
    {
      ++it;
      c = *it;
      n = c->getInfo()->parent_active_cnts;
    }
    if ( n <=  PrntCntThrsh_ ){
      it = pool_.erase(it);
    }
  }
*/
  if (pool_.size() > PoolSize_ - 1){
    std::list<CutPtr>::iterator it = pool_.begin();
    it = pool_.erase(it);
  }

  pool_.push_back(cut);
  cut->getInfo()->inRel = false;

}

void CutMan2::NodeIsBranched(NodePtr node, ConstSolutionPtr sol, int num)
{
  CutPtr cut;
  cutList cutlist;
  NodePtr child;
  const double *y = sol->getDualOfCons();
  int i;
  timer_->start();
    
  for (CCIter it=rel_.begin(); it != rel_.end(); ++it){
    cut = *it;
    i = cut->getConstraint()->getId();
    if (y[i] > 1e-6 || y[i] < -1e-6){
      cut->getInfo()->parent_active_cnts += num;
      cutlist.push_back(cut);
    }
  }
  NodeCutsMap_[node] = cutlist; 
  ChildNum_[node] = num;
  double a1 =  timer_->query();
  ctMngrtime_ += a1;
  branchedTime_ += a1;
  timer_->stop();
}

void CutMan2::NodeIsProcessed(NodePtr node)
{
  stats_->callNums++;
  CutPtr cut;

  timer_->start();
  if (NodeCutsMap_.count(node->getParent()) > 0) {
  	//cutList cutlist = NodeCutsMap_[node->getParent()];
  	//for (CCIter it=cutlist.begin(); it!=cutlist.end(); ++it){
	for (CCIter it = NodeCutsMap_[node->getParent()].begin(); 
		it != NodeCutsMap_[node->getParent()].end(); ++it){
	  cut = *it;
          --(cut->getInfo()->parent_active_cnts);
  	}
	ChildNum_[node->getParent()] += -1;
  }
  if (ChildNum_[node->getParent()] == 0){
    NodeCutsMap_.erase(node->getParent());
  }
   
  double a1 =  timer_->query();
  ctMngrtime_ += a1;
  processedTime_ += a1;
  timer_->stop();

  stats_->RelSize += rel_.size();
  stats_->PoolSize += pool_.size();

  ctmngrInfo_.RelSize = rel_.size();
  ctmngrInfo_.PoolSize = pool_.size();
  ctmngrInfo_.RelToPool = stats_->numRelToPool;
  ctmngrInfo_.PoolToRel = stats_->numPoolToRel;
  ctmngrInfo_.RelAve = (double)stats_->RelSize/stats_->callNums;
  ctmngrInfo_.PoolAve = (double)stats_->PoolSize/stats_->callNums;

}

void CutMan2::write(std::ostream &out) const {

  out << " CutManager " << std::endl;

}

void CutMan2::addCuts(CutVectorIter cbeg,CutVectorIter cend)
{
  for (CutVectorIter it=cbeg; it!=cend; ++it){
    addCut(*it);
  }
}


void CutMan2::addCut(CutPtr c)
{
  addToRel_(c);
}

void CutMan2::writeStats(std::ostream &out) const
{
  out << "nothing to do" << std::endl;
}
void CutMan2::writeStat()
{
  std::cout
    << "CutManager: number of cuts added........................ = " << stats_->numAddedCuts << std::endl
    << "CutManager: number of cuts deleted...................... = " << stats_->numDeletedCuts << std::endl
    << "CutManager: number of cuts moved from relaxation to pool = " << stats_->numRelToPool << std::endl
    << "CutManager: number of cuts moved from pool to relaxation = " << stats_->numPoolToRel << std::endl
    << "CutManager: number of calls............................. = " << stats_->callNums << std::endl
    << "CutManager: size of rel................................. = " << rel_.size() << std::endl  
    << "CutManager: size of pool................................ = " << pool_.size() << std::endl
    << "CutManager: average size of rel......................... = " << (double)stats_->RelSize/stats_->callNums << std::endl
    << "CutManager: average size of pool........................ = " << (double)stats_->PoolSize/stats_->callNums << std::endl
    << "CutManager: MaxInactiveInRel............................ = " << MaxInactiveInRel_ << std::endl
    << "CutManager: PrntActCnt.................................. = " << PrntCntThrsh_ << std::endl
    << "CutManager: time........................................ = " << ctMngrtime_ << std::endl
    << "CutManager: Map size.................................... = " << NodeCutsMap_.size() <<  "\n"
    << "CutManager: update cut.................................. = " << updateTime_ << "\n"
    << "CutManager: check cut................................... = " << checkTime_ << "\n"
    << "CutManager: processed................................... = " << processedTime_ << "\n"
    << "CutManager: branchedt................................... = " << branchedTime_ << "\n";

    
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
