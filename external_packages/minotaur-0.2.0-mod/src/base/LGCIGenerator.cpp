//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file LGCIGenerator.cpp 
 * \brief Declare base class LGCIGenerator. 
 * \author Serdar Yildiz, Argonne National Laboratory 
*/
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cmath>
using std::min;
#include <sstream>
using std::ostringstream;
#include <algorithm>
using std::max;
#include <cstring>

//#include "/sandbox/sey309/minotaur-0.1.1-dev-linux-x86_64/third-party/clp-1.14.6/include/coin/CoinPackedMatrix.hpp"
//#include "/sandbox/sey309/minotaur-0.1.1-dev-linux-x86_64/third-party/clp-1.14.6/include/coin/OsiClpSolverInterface.hpp"
//#include "/sandbox/sey309/minotaur-0.1.1-dev-linux-x86_64/third-party/clp-1.14.6/include/coin/CoinPackedVector.hpp"

#include "LGCIGenerator.h"
#include "LinearFunction.h"
#include "Function.h"
#include "Variable.h"
#include "Option.h"
#include "LPEngine.h"
#include "Variable.h"

using namespace Minotaur;

#define DEBUG_LEVEL 9

// Default constructor does nothing.
LGCIGenerator::LGCIGenerator() {}

LGCIGenerator::LGCIGenerator(ProblemPtr p, SolutionPtr s, 
                             EnvPtr env, LPEnginePtr lpengine)
  :  env_(env), p_(p), s_(s), lpengine_(lpengine)
{
  initialize();
  generateAllCuts();
}

LGCIGenerator::LGCIGenerator(RelaxationPtr rel, ConstSolutionPtr sol, 
                             EnvPtr env, LPEnginePtr lpengine)
{
  LGCIGenerator(rel, sol, env, lpengine);
}

LGCIGenerator::~LGCIGenerator()
{
  // deallocate heap.
  if (stats_) {
    delete stats_;
  }
}

void LGCIGenerator::initialize()
{
  // Integrality tolerance is assigned.
  intTol_ = env_->getOptions()->findDouble("int_tol")->getValue();
  stats_ = new LGCIGenStats();
  initStats();
  // We cannot take a copy of these objects, since there is no copy
  // constructor.
  numCons_ = 0;
  // Identify list of knapsack inequalities in the given problem.
  knapList_ = (KnapsackListPtr) new KnapsackList(p_);
  // Identify GUB inequalities in the given problem.
  // Serdar, We want to take the problem structure as a given data
  // For now, it is OK to create problem structure inside here just for test
  // purposes.
  // Later, we will add prob structure as a parameter to constructer.
  probstruct_ = (ConstProbStructPtr) new ProbStructure(p_, env_);
  if (DEBUG_LEVEL >= 0) {
    outfile_ = "GUBs.txt";
    output_.open(outfile_.c_str());
    // Print out the current solution.
    VariableConstIterator it;
    VariableConstIterator begin = p_->varsBegin();
    VariableConstIterator end   = p_->varsEnd();
    VarVector solution;
    for (it=begin; it!=end; ++it) {
      solution.push_back(*it);
    }
    output_.precision(1);
    output_ << "Given primal solution is: " << endl;
    s_->writePrimal(output_, &solution);
    output_ << "\nGiven dual solution is: " << endl;
    s_->writeDual(output_);
    output_.close();
    
  }
}

void LGCIGenerator::initStats()
{
  stats_->liftsubprobs = 0;
  stats_->totalcuts    = 0;
  stats_->cuts         = 0;
  stats_->knapcons     = 0;
  stats_->knapcov      = 0;
  stats_->knapgub      = 0;
  stats_->time         = 0;
}

void LGCIGenerator::generateAllCuts()
{
  // Implement knapsack identifier to use a map instead of a vector.
  // By this way it will be easy to identify if a Knapsack is not a GUB which
  // is necessary here.
  // We check if a knapsack is GUB at the time when we consider each knapsack.
  // If there exists any knapsack constraints in the problem, then proceed.
  if (knapList_->getNumKnaps() >= 1) {
    // Iterators for knapsack constraints.
    ConstraintConstIterator it;
    ConstraintConstIterator begin = knapList_->getListBegin();
    ConstraintConstIterator end   = knapList_->getListEnd();
    // Shows if the knapsack constraint has a cover.
    bool hascover    = false;
    // Shows if there is a set of GUBs that cover knapsack inequality.
    bool hasgubcover = false;
    // Consider each of the knapsack inequalities for cut generation.
    for (it=begin; it!=end; ++it) {
      // If debug mode is active, constraint is written to an output file.
      if (DEBUG_LEVEL >= 9) {
        output_.open(outfile_.c_str(), std::ios_base::app);
        output_ << "\n\n\nConstraint considered is: " << endl;
        (*it)->write(output_);
        output_.close();
      }
      // Increment the number of knapsack constraints considered until now.
      stats_->knapcons += 1;
      // If constraint has a cover then continue.
      hascover    = hasCover(it);
      if (hascover) {
        // Increment the number of knapsacks that has at least one cover.
        stats_->knapcov += 1;
        // GUB list that may cover knapsack ineqiuality is initialized as
        // empty.
        ConstConstraintVectorPtr gublist = 
          (ConstConstraintVectorPtr) new ConstConstraintVector;
        hasgubcover = hasGubCover(it, gublist);
        // If there is a GUB set that covers knapsack, then generate cuts.
        if (hasgubcover) {
          stats_->knapgub += 1;
          generateCuts(*it,gublist);
        } // end hasgubcover
      } // end hascover
    } // end for loop
    // If no knapsack constraint is used for cut generation print message.
  } else { 
    cerr << "LGCIGenerator::generateAllCuts. No constraint used for cut generation.\n";
  }
}

bool LGCIGenerator::hasCover(ConstraintConstIterator itcons)
{
  bool has = false;
  // The constraint considered.
  ConstraintPtr cons= *itcons; 
  double b = cons->getUb();
  // Iterators for variables in function.
  LinearFunctionPtr lf = cons->getLinearFunction();
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  // Sum of coefficients of the variables.
  double sum = 0.0;
  // double coefficient of current variable.
  double coeff = 0.0;
  // Iterate through all the variables in the linear function.
  for (it=begin; it!=end; ++it) {
    coeff = it->second;
    // Check if the sum is greater than rhs, i.e. \sum a_j > b.
    sum += coeff;
    if (sum > b) {
      has = true;
      break;
    }
  }
  // No cover exists.
  return has;
}


// This function has to be rewwritten by using max and min maps.
// This is really inefficient.
bool LGCIGenerator::hasGubCover(ConstraintConstIterator itcons, 
                                ConstConstraintVectorPtr gublist)
{
  // Initialize how many times the variable has not been covered.
  CoverSetPtr numuncovered = (CoverSetPtr) new CoverSet();
  // The fractional values of variable.
  CoverSetPtr fractions  = (CoverSetPtr) new CoverSet();
  // Variable is covered or not.
  CoverSetPtr covered    = (CoverSetPtr) new CoverSet();
  // Obtain the linear function.
  const LinearFunctionPtr lf = (*itcons)->getLinearFunction();
  // Iterators for variable.
  VariableGroupConstIterator it0;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  // Primal solution.
  const double * x = s_->getPrimal();
  // Current variable.
  ConstVariablePtr var;
  // Index of current variable.
  UInt index = 0;
  // fractional solution of current variable.
  double value = 0.0;
  // how many times current variable is uncovered.
  double curnumuncov = 0.0;
  // Fractional value of current variable.
  double fraction = 0.0;
  // Initialize the vectors that shows the degree of uncoveredness of the variable,
  // is it covered or not, what is the fractional part of solution value of variable.
  // Iterate through all variables.
  for(it0=begin; it0!=end; ++it0) {
    var = it0->first;
    index = var->getIndex();
    value = x[index];
    // Assumes 0 <= value <= 1.
    fraction = min(value, 1.0-value);
    curnumuncov = fraction;
    VariableValuePair uncov(var,curnumuncov); 
    numuncovered->push_back(uncov);
    VariableValuePair frac(var,fraction);
    fractions->push_back(frac);
    VariableValuePair varcover(var,0);
    covered->push_back(varcover);
  }

  // While there is some variables uncovered add more GUB's.
  UInt mincovered = 0;
  // Variable to be covered.
  ConstVariablePtr vartocover;
  // How to compare, use the decreasing order. Thus, thw minimum element is
  // the maximum element.
  CompareValueVariablePair compare;
  // The variable pair identified to cover.
  CoverSetIterator varpair;
  
  // This outer loop runs until all elements in the knapsack are covered.
  do {

    // For debug, print the current fractionals, covered and  number of
    // uncoveredness.
    if (DEBUG_LEVEL >= 9) {
      printIneq(fractions, covered->size()-1, Cons, "Fractions.");
      printIneq(covered, covered->size()-1, Cons, "Covered vector.");
      printIneq(numuncovered, numuncovered->size()-1, Cons, "NumUncovered vector.");
    }
    
    // Check if there exist an uncovered element and get the maximum uncovered one.
    bool uncovered = false;
    while (true) {      
      varpair = std::min_element(numuncovered->begin(), numuncovered->end(), compare);
      if (varpair->second == -1) {
        break;
      }
      vartocover = varpair->first;
      for (CoverSetIterator it=covered->begin(); it!=covered->end(); ++it) {
        if ((it->first->getId() == vartocover->getId()) &&
            (it->second == 0) ) {
          //covered->erase(it);
          it->second = 1;
          uncovered = true;
          break;
        }
      }
      if (uncovered) {
        break;
      }
    }

    // Print out variable to cover.
    if (DEBUG_LEVEL >= 9) {
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "Variable to cover is " << vartocover->getName() << endl;
      output_.close();
    }

    
    // The GUB that covers the variable and that has the largest summation of
    // total uncoveredness of all elements it include is selected.
    ConstConstraintVectorPtr vargubs = probstruct_->getVarGUBs(vartocover);
    ConstConstraintPtr maxgub;
    double maxtotaluncovered = -10000000000; // Serdar Change this to infinity.
    for (ConstConstraintVector::iterator it = vargubs->begin(); it!=vargubs->end(); ++it) {
      double totaluncovered = 0.0;
      for (VariableGroupConstIterator it2=(*it)->getLinearFunction()->termsBegin();
           it2!=(*it)->getLinearFunction()->termsEnd(); ++it2) {
        for (CoverSetIterator it3=numuncovered->begin(); it3!=numuncovered->end(); ++it3) {
          if (it2->first->getId() == it3->first->getId()) {
            totaluncovered += it3->second;
            break;
          }
        }
      }
      if (totaluncovered >= maxtotaluncovered) {
        maxgub = *it;
        maxtotaluncovered = totaluncovered;
      }
    }

    // Print out the GUB selected to cover the variable.
    if (DEBUG_LEVEL >= 9) {
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "GUB selected to cover the variable is Constraint " 
              << maxgub->getName() << " with id " << maxgub->getId() << endl;
      output_.close();
    }

    // Now, we cover the variables in the added GUB, so we can take it out
    // from consideration in the proceeding loops.
    VariableGroupConstIterator it3;
    VariableGroupConstIterator begin3 = maxgub->getLinearFunction()->termsBegin();
    VariableGroupConstIterator end3   = maxgub->getLinearFunction()->termsEnd();
    VariablePtr curvarlf;
    VariablePtr curvargub;
    for (it3=begin3; it3!=end3; ++it3) {
      curvargub = it3->first;
      for (CoverSetIterator it4=covered->begin(); it4!=covered->end(); ++it4) {
        if (it4->first == curvargub) {
          it4->second = 1;
          break;
        }
      }

      for (CoverSetIterator it5=numuncovered->begin(); it5!=numuncovered->end(); ++it5) {
        if (it5->first == curvargub) {
          it5->second = -1;
          break;
        }
      }
    }

    if (DEBUG_LEVEL >= 9) {
      printIneq(covered, covered->size()-1, Cons, "Covered vector after covering variable.");
      printIneq(numuncovered, numuncovered->size()-1, Cons, "NumUncovered vector af ter covering variable.");
    }

    // Add the GUB to the list.
    gublist->push_back(maxgub);
    
    // Check if there is any variable uncovered. 
    mincovered = (UInt) (std::max_element(covered->begin(), covered->end(), compare) )->second;
  } while (mincovered == 0);

  return true;
}

void LGCIGenerator::generateCuts(ConstConstraintPtr cons, 
                                 ConstConstraintVectorPtr gublist)
{
  // Generate a GNS LGCI cut.
  GNS(cons,gublist);
}

bool LGCIGenerator::GNS(ConstConstraintPtr cons, ConstConstraintVectorPtr gublist)
{
  // Initial cover is empty and filled later.
  CoverSetPtr cover = (CoverSetPtr) new CoverSet();
  bool covgenerated = coverSetGeneratorGNS(cons,cover);
  // If no cover generated terminate.
  if (covgenerated == false) {
    return false;
  }
  
  // Obtain the minimal cover.
  minimalCover(cover, cons);

  // Initialize sets C, C1, C2, F and R as described in GNS paper.
  CoverSetPtr cone  = (CoverSetPtr) new CoverSet();
  CoverSetPtr ctwo  = (CoverSetPtr) new CoverSet();
  CoverSetPtr fset  = (CoverSetPtr) new CoverSet();
  CoverSetPtr rset  = (CoverSetPtr) new CoverSet();
  // Construct the partitions. 
  coverPartitionGNS(cons, cover, cone, ctwo, fset, rset);

  // Cover inequality is initialized as empty.
  CoverSetPtr ineq = (CoverSetPtr) new CoverSet();
  double b = 0.0;
  bool generated = liftingGNS(cone,ctwo,fset,rset,ineq,cons,gublist,b);
  if (generated == true) {
    CutFail failtype;
    bool cutgen = addCut(ineq, b, Gns, failtype);
    return cutgen;
  } else {
    stats_->noinitcov += 1;
    return false;
  }

  return generated;
}

/** Assumption: An empty cover is given.
    Initial cover will be generated from scratch.
 */
bool LGCIGenerator::coverSetGeneratorGNS(ConstConstraintPtr cons, 
                                         CoverSetPtr cover)
{
  // This is the knapsack constraint to be covered.
  const LinearFunctionPtr lf = cons->getLinearFunction();
  // Obtain nonzero variables in the given solution.
  CoverSetPtr nonzerovars = (CoverSetPtr) new CoverSet();
  CoverSetPtr zerovars    = (CoverSetPtr) new CoverSet();
  nonzeroVars(lf, nonzerovars, zerovars);
  // Sort variables in nonincreasing order of thir values in the given solution.
  sortNonIncreasing(nonzerovars);

  // Cover is being obtained from knapsack inequality.
  // Iterators for nonzeros vector.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = nonzerovars->begin();
  CoverSetConstIterator end   = nonzerovars->end(); 
  // The partial summation of variable coefficients in knapsack constraint. 
  double sum = 0;
  // right hnad side "b"
  double b = cons->getUb();
  // Current pair.
  VariableValuePair pair;
  // Current variable.
  VariablePtr variable; 
  // Coefficient of current variable.
  double coeff = 0.0;
  // Assumption is that "b" is positive. Let's add an assert later.
  for (it=begin; it!=end; ++it) {
    // Add variables until sum exceeds right hand side "b".
    if (sum <= b) {
      pair = *it;
      variable = pair.first;
      coeff = lf->getWeight(variable);
      VariableValuePair newpair(variable,coeff);
      cover->push_back(newpair);
      sum += coeff;
    } else {
      break;
    }
  }
  
  // Check if we obtain a cover.
  if (sum <= b) {
    // No cover generated.
    if (DEBUG_LEVEL >= 9) {
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "No initial cover obtained from GNS." << endl;
      output_.close();
    }
    return false;
  } else { 
    // Cover generated.
    if (DEBUG_LEVEL >= 9) {
      printIneq(cover,cover->size()-1, Cover, "GNS initial cover generated.");
    }
    return true;
  }
}

void LGCIGenerator::nonzeroVars(const ConstLinearFunctionPtr lf,
                                CoverSetPtr nonzerovars,
                                CoverSetPtr zerovars)
{
  // Iterators to iterate on variables in knapsack constraint.
  VariableGroupConstIterator it;
  const VariableGroupConstIterator begin = lf->termsBegin();
  const VariableGroupConstIterator end   = lf->termsEnd();
  
  // Current variable.
  ConstVariablePtr variable;
  // Index of variable in solution.
  UInt varindex = 0;
  // Variable value in solution.
  double varvalue = 0.0;
  // Solution vector.
  const double * const x = s_->getPrimal();
  // Fill in nonzeros and zeros vectors.
  // Iterate through all the variables in knapsack constraint.
  for (it=begin; it!=end; ++it) {
    variable = it->first;
    varindex = variable->getIndex();
    varvalue = x[varindex];
    // Variable to be added to list nonzeros or zeros.
    VariableValuePair currentvar(variable, varvalue);
    if (varvalue != 0) { // add a tolerance here.
      nonzerovars->push_back(currentvar);
    } else {
      zerovars->push_back(currentvar);
    }
  }
}

void LGCIGenerator::sortNonIncreasing(CoverSetPtr nonzeros)
{
  /**
   * We sort in nonincreasing order the variables according to their values in
   * fractional solution.
   */
  CompareValueVariablePair compare;
  CoverSetIterator begin = nonzeros->begin();
  CoverSetIterator end   = nonzeros->end();
  sort(begin, end, compare);
}

/** Assumption: 
 *  cover is an minimal cover.
 *  cone, ctwo, fset, and fbar are empty sets.*/
void LGCIGenerator::coverPartitionGNS(const ConstConstraintPtr cons,
                                      const ConstCoverSetPtr cover,
                                      CoverSetPtr cone,
                                      CoverSetPtr ctwo,
                                      CoverSetPtr fset,
                                      CoverSetPtr rset)
{
  // Partition cover ser C into C1 and C2.
  // Iterators for cover set C.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = cover->begin();
  CoverSetConstIterator end   = cover->end();
  // Get the primal solution.
  const double * x = s_->getPrimal();
  // Current variable.
  ConstVariablePtr variable;
  // Index of variable. 
  UInt varindex = 0;
  // Variable value in solution.
  double value = 0.0;
  // Since |C1| >= 1, we change it as |C1| >= 2. 
  // Serdar explain this, if only one then trivial liftings!
  // This is the number of remaining variables in cover for C1.
  // This can be relaxed to |C1| >= 1 in the future.
  UInt numremaining = cover->size();
  // Iterate through elements of C.
  for (it=begin; it!=end; ++it) {
    variable = it->first;
    varindex = variable->getIndex();
    value    = x[varindex];
    // If x_j^* == 1 then add it to C2.
    if (value == 1 && (numremaining >= 3)) {
      ctwo->push_back(*it);
      numremaining -= 1;
    } else {
      // If x_j^* != 1 then add it to C1.
      cone->push_back(*it);
    }
  }

  // Partition set cbar (N-C) into sets F and R as described by GNS.
  // Construct N-C.
  CoverSetPtr cbar = (CoverSetPtr) new CoverSet();
  cBar(cover, cbar, cons); 

  // Divide cbar into sets F and R.
  // Iterators for cbar.
  CoverSetConstIterator itcbar;
  CoverSetConstIterator begincbar = cbar->begin();
  CoverSetConstIterator endcbar   = cbar->end();
  // Current variable.
  ConstVariablePtr varcbar;
  // Index of variable.
  UInt varcbarind = 0;
  // Variable value in solution.
  double valuecbar = 0.0;
  // Iterate through elements of C-N.
  for (itcbar=begincbar; itcbar!=endcbar; ++itcbar) {
    varcbar    = itcbar->first;
    varcbarind = varcbar->getIndex();
    valuecbar  = x[varcbarind];
    // If x_j^* > 0 then add it to F.
    if (valuecbar > 0) {
      fset->push_back(*itcbar);
    } else {
      // If x_j^* <= 0 then add it to R.
      rset->push_back(*itcbar);
    }
  }

  // Print out the cover debugging.
  if (DEBUG_LEVEL >= 9) {
    printIneq(cover, cover->size()-1, Set, "Minimal cover for GNS cover partition.\n C");
    printIneq(cbar,  cbar->size()-1,  Set, "N-C");
    printIneq(cone,  cone->size()-1,  Set, "C1");
    printIneq(ctwo,  ctwo->size()-1,  Set, "C2");
    printIneq(fset,  fset->size()-1,  Set, "F");
    printIneq(rset,  rset->size()-1,  Set, "R");
  }

}

/** Lifting strategy of GNS.
 * Assumptions:
 * Sets C1, C2, F and R are given as described in paper of GNS (not empty).
 * For lifting sequence, we do not consider GNS for GUBwise lifting.
 * The reason is we do not solve lifting problem exactly.
 * The GUBwise lifting decreases computational effort to solve consequent
 * lifting problems.
 * We simply use GNS lifting sequence as in given order order in their set.
 * This order can be changed later.
 * Inequality is given as empty.
 * Vectors are not sorted or anything else.
 */
bool LGCIGenerator::liftingGNS(const ConstCoverSetPtr cone,
                               const ConstCoverSetPtr ctwo,
                               const ConstCoverSetPtr fset,
                               const ConstCoverSetPtr rset,
                               CoverSetPtr ineq,
                               const ConstConstraintPtr cons,
                               ConstConstraintVectorPtr gublist,
                               double & rhs)
{
  // Generate cover set version of GUBs from constraint type.
  // This is exactly the same as GUB constraints but in cover set data ype.
  boost::shared_ptr< std::vector<CoverSetPtr> > guborigcoeffs = 
    (boost::shared_ptr<std::vector<CoverSetPtr> >)   new std::vector<CoverSetPtr>;
  // This one will be modified to be non-overlapping.
  boost::shared_ptr< std::vector<CoverSetPtr> > gubcoverlist =
    (boost::shared_ptr< std::vector<CoverSetPtr> >)  new std::vector<CoverSetPtr>;
  // GUB lists are generated.
  generateGubMaps(gublist, gubcoverlist, guborigcoeffs);
  // Get the non-overlapping GUB set that covers knapsack inequality.
  nonOverlap(gubcoverlist);
  nonOverlap(guborigcoeffs);

  // Initialize needed data elements.
  // Initialize inequality as C1 set and rhs of inequality as |C1|-1.
  rhs = cone->size() - 1;
  // Coefficients of variables in the initial cover inequality is assigned to be ones.
  initCoverIneq(cone, ineq);
  
  // Construct objective function.
  CoverSetPtr obj = (CoverSetPtr) new CoverSet(*ineq);
  // Construct initial constraints of lifting problem.
  // First constraint is similar to C1.
  CoverSetPtr consknap = (CoverSetPtr) new CoverSet(*cone);
  // remaining constraints are related to GUB constraints.
  boost::shared_ptr< std::vector<CoverSetPtr> > gubcons =  
    (boost::shared_ptr< std::vector<CoverSetPtr> >) new std::vector<CoverSetPtr>;
  initGubCons(cone,gubcoverlist, gubcons);
  
  // Construct rhs of lifting problem, i.e. for knapsack it is b, 
  // and initialize rhs of GUB inequalities as ones.
  double initialbknap = cons->getUb();
  UInt numgubcons = gubcons->size();
  double * initialbgub = new double[numgubcons];
  std::fill_n(initialbgub, numgubcons, 1.0);
  // Change the rhs according to downlifted variables.
  if (ctwo->size() >= 1) {
    // Iterators for elements of C2.
    CoverSetConstIterator it;
    CoverSetConstIterator begin = ctwo->begin();
    CoverSetConstIterator end   = ctwo->end();
    // This can be simplified too much, not efficient yet.
    for (it=begin; it!=end; ++it) {
      // Update rhs of knapsack.
      initialbknap -= it->second;
      // Update GUB rhs.
      for (UInt index=0; index<numgubcons; ++index) {
        // Iterators for elements of GUB elements
        CoverSetConstIterator itgub;
        CoverSetConstIterator begingub = (*gubcoverlist)[index]->begin();
        CoverSetConstIterator endgub   = (*gubcoverlist)[index]->end();
        // Iterate thorugh elements GUB.
        for (itgub=begingub; itgub!=endgub; ++itgub) {
          // Check if element is the same.
          if (it->first->getId() == itgub->first->getId()) {
            initialbgub[index] -= 1;
          }
        }
      }
    }
  }
  
  // Print the rhs of knapsack and GUBs.
  if (DEBUG_LEVEL >= 9) {
    output_.open(outfile_.c_str(), std::ios_base::app);
    output_ << "Rhs of knapsack inequality for initial problem is " << initialbknap << endl;
    output_ << "Rhs of GUBs for initial problem." << endl;
    for (UInt i = 0; i<gubcons->size(); ++i) {
      output_ << initialbgub[i] << " ";
    }
    output_ << endl;
    output_.close();
    // Print out knapsack constraint.
    printIneq(cone, initialbknap, Cons, "Initial knapsack constraint");
  }

  bool liftup = true;
  // We are going to uplift variables in set F.
  if (fset->size() >= 1) {
    // Uplift if set F is nonempty.
    liftSet(obj, guborigcoeffs, consknap, gubcons, fset, ineq, rhs, initialbknap, initialbgub, liftup);
  }

  if (ctwo->size() >= 1) {
    // Variables in C2 are going to be down lifted.
    liftup = false;
    // Lift down the variables.
    liftSet(obj, guborigcoeffs, consknap, gubcons, ctwo, ineq, rhs, initialbknap, initialbgub, liftup);
  }
  
  if (rset->size() >= 1) {
    // Variables in set R are going to be up-lifted.
    liftup = true;
    // Lift up variables.
    liftSet(obj, guborigcoeffs, consknap, gubcons, rset, ineq, rhs, initialbknap, initialbgub, liftup);
  }
  
 

  delete [] initialbgub;

  return true;
}

/** Assumptions: Following data elements should be suitable for lifitng
 * problem construction.
 * Objective, lifting problem constraints, variable set to be lifted,
 * inequality to be lifted, rhs of inequality, rhs of lifting problem
 * constraints.
 * Under these assumptions, this function determine the coefficient of lifted
 * variable, 
 * adds the new variable to inequality being lifted,
 * lifting problem objective function,
 * and lifting problem constraints.
 * obj, consknap, gubcons, rhs, initialbknap, initialbgub are changed at the
 * end.
 * varset and liftup does not change.
 */
void LGCIGenerator::liftSet(CoverSetPtr obj,
                            boost::shared_ptr<std::vector<CoverSetPtr> > origgubs,
                            CoverSetPtr consknap,
                            boost::shared_ptr<std::vector<CoverSetPtr> > gubcons,
                            const ConstCoverSetPtr varset,
                            CoverSetPtr coverineq,
                            double & rhs,
                            double & initialbknap,
                            double * initialbgub,
                            bool liftup)
{
  // Check if the set being lifted is empty.
    if (varset->empty() == false) {
      // Iterators for the set of variables to be lifted.
      CoverSetConstIterator it;
      CoverSetConstIterator begin = varset->begin();
      CoverSetConstIterator end   = varset->end();
      // Initialize coefficient of variable to be lifted.
      double alpha = 0.0;
      // Lift the variables one by one in the given order.
      for (it=begin; it!=end; ++it) {
        // Lift the variable.
        alpha = lift(obj, origgubs, consknap, gubcons, it, rhs, initialbknap, initialbgub, liftup);
        
        if (DEBUG_LEVEL >= 9) {
          output_.open(outfile_.c_str(), std::ios_base::app);
          output_ << "Alpha obtained is " << alpha << endl;
        }
        
        // Construct new pair for objective.
        VariableValuePairPtr newobjvar = 
          (VariableValuePairPtr) new VariableValuePair(*it);
        // Update coefficient of variable with alpha.
        newobjvar->second = alpha;
        // Update objective function.
        obj->push_back(*newobjvar);
        // Update cover inequality.
        coverineq->push_back(*newobjvar);
        // Update knapsack problem constraint.
        consknap->push_back(*it);

        // Update GUB constraints.
        // Iterator for GUB constraints.
        std::vector<CoverSetPtr>::iterator itgub;
        std::vector<CoverSetPtr>::iterator begingub = origgubs->begin();
        std::vector<CoverSetPtr>::iterator endgub   = origgubs->end();
        // Iterate thorugh all GUB constraints.
        UInt index = 0;
        for (itgub=begingub; itgub!=endgub; ++itgub) {
          // Iterators for elements of GUB.
          CoverSetIterator itel;
          CoverSetIterator beginel = (*itgub)->begin();
          CoverSetIterator endel   = (*itgub)->end();
          // Iterate through elements of GUB.
          for (itel=beginel; itel!=endel; ++itel) {
            // If initially element is included in GUB, then add it again.
            if (itel->first == it->first) {
              VariableValuePair newvar = VariableValuePair(*it);
              newvar.second = 1.0; 
              (*gubcons)[index]->push_back(newvar);
            }
          }
          index += 1;
        }

        if (DEBUG_LEVEL >= 9) {
          output_.close();
          output_.open(outfile_.c_str(), std::ios_base::app);
          printIneq(coverineq, rhs, Cons, "Updated cover inequality.");
          output_.close();
        }
      }
    }
}

/**
 * This is used to call OsiClpInterface.
 * This function prepares the problem data for lifting problem solver.
 * Updates rhs of lifting problem and rhs of cover inequality.
 */
// double LGCIGenerator::lift(CoverSetPtr obj,
//                          CoverSetPtr consknap,
//                          boost::shared_ptr<std::vector<CoverSetPtr> > gubcons,
//                          const CoverSetConstIterator variable,
//                          double & rhs,
//                          double & initialbknap,
//                          double * initialbgub,
//                          bool liftup)
// {
  // // Increment number of relaxations solved.
  // stats_->liftsubprobs += 1;

  // // Write the problem as Clp expects.
  // OsiClpSolverInterface * si = new OsiClpSolverInterface();
  // // Number of variables.
  // UInt n_cols = p_->getNumVars(); // number of variables.
  // double * objective = new double[n_cols]; // objective coefficients
  // double * col_lb    = new double[n_cols]; // variable lower bounds
  // double * col_ub    = new double[n_cols]; // variable upper bounds
  // // Define objective coefficients.
  // // Iterators for objective function.
  // CoverSetConstIterator it;
  // CoverSetConstIterator begin = obj->begin();
  // CoverSetConstIterator end   = obj->end();
  // UInt index = 0;
  // for (it=begin; it!=end; ++it) {
  //   objective[index] = it->second;
  //   index += 1;
  // }

  // // Define variable lower upper bounds.
  // std::fill_n(col_lb,n_cols,0.0);
  // std::fill_n(col_ub,n_cols,1.0);

  // // Number of constraints , i.e. 1(knapsack)+nuumber of GUB constraints.
  // UInt n_rows = 1 + gubcons->size(); // number of constraints.
  // double * row_lb = new double[n_rows]; // constraint lower bounds.
  // double * row_ub = new double[n_rows]; // constraint upper bounds.
  
  // // Define constraint matrix.
  // CoinPackedMatrix matrix(false,0,0);
  // matrix.setDimensions(0,n_cols);
  // // The list of constraints that includes variables.
  // // Add knapsack constraint to the matrix.
  // CoinPackedVector knapsack;
  // // Iterators for variables of knapsack.
  // CoverSetConstIterator itknap;
  // CoverSetConstIterator beginknap = consknap->begin();
  // CoverSetConstIterator endknap   = consknap->end();
  // // Iterate through all elements of knapsack.
  // UInt varindex = 0;
  // for  (itknap=beginknap; itknap!=endknap; ++itknap) {
  //   if (itknap->second != 0) {
  //     varindex = itknap->first->getIndex();
  //     knapsack.insert(varindex, itknap->second);
  //   }
  // }
  // // Serdar change this to -inf later. This may result in an infeasible
  // // lifting problem.
  // //row_lb[0] = 0.0; 
  // row_lb[0] = si->getInfinity();
  // row_ub[0] = initialbknap;
  // matrix.appendRow(knapsack);
  // // Add the GUB constraints to the matrix.
  // // Iterators for GUB constraints.
  // std::vector<CoverSetPtr>::const_iterator itgublist;
  // std::vector<CoverSetPtr>::const_iterator begingublist = gubcons->begin();
  // std::vector<CoverSetPtr>::const_iterator endgublist   = gubcons->end();
  // // Iterate through all GUB constraints.
  // UInt consindex = 1;
  // for (itgublist=begingublist; itgublist!=endgublist; ++itgublist) {
  //   CoinPackedVector gubnew;
  //   // Iterators for elements of GUB.
  //   CoverSetConstIterator itgubel;
  //   CoverSetConstIterator begingubel = (*itgublist)->begin();
  //   CoverSetConstIterator endgubel   = (*itgublist)->end();  
  //   // Iterate through all elements of the GUB constraint.
  //   for (itgubel=begingubel; itgubel!=endgubel; ++itgubel) {
  //     if (itgubel->second != 0) {
  //       varindex = itgubel->first->getIndex();
  //       gubnew.insert(varindex,itgubel->second);
  //     }
  //   }
  //   row_lb[consindex] = 0;
  //   row_ub[consindex] = 1;
  //   matrix.appendRow(gubnew);
    
  //   consindex += 1;
  // }
  
  // // Load problem to OSI.
  // si->loadProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);
  // // we want to maximize the objective function.
  // si->setObjSense(-1);
  // // solve the linear program.
  // si->initialSolve();

  // double alpha = roundHeur(si);

  // // To be continued, change rhs etc, update for uplift and downlift.


  // delete si;
  // return alpha;

  // Delete this later.
//   return 0;
// }

void LGCIGenerator::addCons(CoverSetPtr obj,
                            CoverSetPtr consknap,
                            double bknap,
                            boost::shared_ptr< std::vector<CoverSetPtr> > gubcons,
                            double * bgubs,
                            OrigLiftVarsPtr varmap,
                            ProblemPtr liftprob)
{
  LinearFunctionPtr lfobj = (LinearFunctionPtr) new LinearFunction();
  ObjectiveType otype     = Maximize;
  double cb = 0;
  // Construct objective and its variables.
  CoverSetIterator it;
  CoverSetIterator begin = obj->begin();
  CoverSetIterator end   = obj->end();
  // Iterate through each element of objective.
  for (it=begin; it!=end; ++it) {
    VariablePtr origvar = it->first;
    // check if variable exits in the map
    VariablePtr liftvar;
    liftvar = addVar(origvar, varmap, liftprob);
    double coeff = it->second;
    if (coeff != 0) {
      lfobj->addTerm(liftvar, coeff);
    }
  } // end of for. 

  // Generate objecive function.
  FunctionPtr funobj = (FunctionPtr) new Function(lfobj);
  // Add the objective function.
  liftprob->newObjective(funobj, cb, otype);


  // Add knapsack constraint.
  double ub = bknap;
  LinearFunctionPtr lfknap = (LinearFunctionPtr) new LinearFunction();
  // Iterators for elements of knapsack.
  CoverSetIterator itknap;
  CoverSetIterator beginknap = consknap->begin();
  CoverSetIterator endknap   = consknap->end();
  // Iterate through elements of knapsack constraint.
  for (itknap=beginknap; itknap!=endknap; ++itknap) {
    VariablePtr origknapvar = itknap->first;
    double coeff = itknap->second;
    VariablePtr liftknapvar = addVar(origknapvar, varmap, liftprob);
    if (coeff != 0) {
      lfknap->addTerm(liftknapvar, coeff);
    }
  } // end of for.
  // Generate knapsack function.
  FunctionPtr funknap = (FunctionPtr) new Function(lfknap);
  double lb = 0.0;
  liftprob->newConstraint(funknap, lb, ub);

  // Add the gub constraints.
  // Iterators for GUBs.
  std::vector<CoverSetPtr>::iterator itgubs;
  std::vector<CoverSetPtr>::iterator begingubs = gubcons->begin();
  std::vector<CoverSetPtr>::iterator endgubs   = gubcons->end();
  UInt index = 0;
  for (itgubs=begingubs; itgubs!=endgubs; ++itgubs) {
    double ubgub = bgubs[index];
    LinearFunctionPtr lfgub = (LinearFunctionPtr) new LinearFunction();
    // Iterators for elements of GUB.
    CoverSetIterator itel;
    CoverSetIterator beginel = (*itgubs)->begin();
    CoverSetIterator endel   = (*itgubs)->end();
    for (itel=beginel; itel!=endel; ++itel) {
      VariablePtr origgubvar = itel->first;
      double coeff = itel->second;
      VariablePtr liftgubvar = addVar(origgubvar, varmap, liftprob);
      if (coeff != 0) {
        lfgub->addTerm(liftgubvar, coeff);
      }
    } // end of elements of gub for.
    // Generate gub function.
    FunctionPtr fungub = (FunctionPtr) new Function(lfgub);
    // We assume every coefficient is non-negative.
    double lbgub = 0.0;
    liftprob->newConstraint(fungub, lbgub, ubgub); 

    index += 1;
  } // end of gub for.
    
}

/** Checks if the variable is already added. 
 * If added it returns a pointer to the variable.
 * Otherwise, it creates a new variable and sends its pointer.
 */
VariablePtr LGCIGenerator::addVar(VariablePtr var, OrigLiftVarsPtr varmap, 
                                  ProblemPtr liftprob)
{
  // Iterators for map that maps variables in the original problem to the lifting problem.
  OrigLiftVars::iterator it;
  OrigLiftVars::iterator end;
  // Check if variable is added already.
  it  = varmap->find(var);
  end = varmap->end();
  // If variable exists send its pointer.
  if (it!=end) {
    return (it->second);
  } else {
    // Variable does not exist.
    // Create a new variable.
    double lb = 0.0;
    double ub = 1.0;
    VariableType vtype = Continuous;
    // Add new variable to lifting problem.
    VariablePtr newvar = liftprob->newVariable(lb, ub, vtype);
    // Add new variable to variable map.
    varmap->insert(std::pair<ConstVariablePtr, ConstVariablePtr>(var,newvar));
    return newvar;
  }
  
}

double LGCIGenerator::lift(CoverSetPtr obj,
                           boost::shared_ptr<std::vector<CoverSetPtr> > origgubs,
                           CoverSetPtr consknap,
                           boost::shared_ptr<std::vector<CoverSetPtr> > gubcons,
                           const CoverSetConstIterator variable,
                           double & rhs,
                           double & initialbknap,
                           double * initialbgub,
                           bool liftup)
{
  if (DEBUG_LEVEL >= 9) {
    output_.close();
    output_.open(outfile_.c_str(), std::ios_base::app);
    if (output_.fail()) {
      cerr << "Error!";
    }
    output_ << "Variable being lifted is " << variable->first->getName() << endl;
    output_.close();
  }

  // Create lifting problem.
  ProblemPtr liftprob = (ProblemPtr) new Problem();
  // Clean the lp engine before use.
  lpengine_->clear();
  OrigLiftVarsPtr varmap = (OrigLiftVarsPtr) new OrigLiftVars();

  double knapb = 0.0;
  // Uplift specific adjustment.
  if (liftup == true) {
    knapb = initialbknap - variable->second;
  } else {
    // Downlifting.
    knapb = initialbknap + variable->second;
  }

  // Adjust GUB rhs.
  // Iterator for GUB constraints.
  std::vector<CoverSetPtr>::iterator it;
  std::vector<CoverSetPtr>::iterator begin = origgubs->begin();
  std::vector<CoverSetPtr>::iterator end   = origgubs->end();
  UInt index = 0;
  double * bgub;
  double * copybgub = 0;
  if (liftup == false) {
    for (it=begin ; it!=end; ++it) {
      // Iterators for variables of GUB.
      CoverSetIterator itel;
      CoverSetIterator beginel = (*it)->begin();
      CoverSetIterator endel   = (*it)->end();
      for (itel=beginel; itel!=endel; ++itel) {
        if (variable->first == itel->first) {
          if (liftup == true) {
            initialbgub[index] -= 1;
          } else {
            initialbgub[index] += 1;
          }
          break;
        }
      } // end of for.
      index += 1;
    } // end of for.
    bgub = initialbgub;
  } else {
    copybgub = new double[gubcons->size()];
    std::memcpy(copybgub, initialbgub, sizeof(double)*(gubcons->size()));
    // Uplifting case.
    for (it=begin ; it!=end; ++it) {
      // Iterators for variables of GUB.
      CoverSetIterator itel;
      CoverSetIterator beginel = (*it)->begin();
      CoverSetIterator endel   = (*it)->end();
      for (itel=beginel; itel!=endel; ++itel) {
        if (variable->first == itel->first) {
          if (liftup == true) {
            copybgub[index] -= 1;
          } else {
            copybgub[index] += 1;
          }
          break;
        }
      } // end of for.
      index += 1;
    } // end of for.
    bgub = copybgub;
  }

  // Define problem.
  addCons(obj, consknap, knapb, gubcons, bgub, varmap, liftprob); 

  // Prepare for solve    
  liftprob->prepareForSolve();

  // write problem.
  if (DEBUG_LEVEL >= 9) {
    output_.open(outfile_.c_str(), std::ios_base::app);
    liftprob->write(output_);
    output_.close();
  }

  // Load the problem to solver.
  lpengine_->load(liftprob);
  
  // Solve problem.
  lpengine_->solve();

  // Solution of lifting problem is obtained.
  double gamma = roundHeur(liftprob);
  
  double alpha = 0.0;
  double ksi   = 0.0;
  if (liftup == true) {
    // If uplifting.
    alpha = rhs - gamma;
    delete [] copybgub;
    return alpha;
  } else {
    // If downlifting.
    ksi   = gamma - rhs;
    rhs += ksi;
    initialbknap = knapb; 
    
    if (DEBUG_LEVEL >= 9) {
      output_.close();
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "gamma " << gamma << endl;
      output_ << "rhs " << rhs << endl;
      output_ << "ksi " << ksi << endl; 
      output_.close();
    }

    return ksi;
  }
}

// Assumes that lifting problem is solved already.
double LGCIGenerator::roundHeur(ProblemPtr prob)
{
  std::vector<UInt> indices;
  std::vector<double> values;
  EngineStatus status = lpengine_->getStatus();
  // Check if problem is solved optimally or else.
  if ((status == ProvenOptimal) || (status == FailedFeas)) { // Here check if we have a feasible solution.  
    double gamma = 0.0;
    // Serdar check that, since it considers minimize it gives negative result.
    double gammainitial = -lpengine_->getSolutionValue();
    ConstSolutionPtr sol = lpengine_->getSolution();
    const double * solution = sol->getPrimal();
    UInt n_cols = prob->getNumVars();
    
    if (DEBUG_LEVEL >= 9) {
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "Gamma initial obtained is " << gammainitial << endl;
      output_ << "Solution of lifting subproblem is: " << endl;
      for (UInt i=0; i<n_cols; ++i) {
        output_ << solution[i] << " " << std::flush;
      }
      output_ << endl;
      output_.close();
    }

    //UInt n_rows = prob->getNumCons();
    // Get list of nonzero values in solution.
    for (UInt i=0; i<n_cols; ++i){
      if (solution[i] != 0) {
        indices.push_back(i);
        values.push_back(solution[i]);
      }
    }
    
    std::vector<UInt> fracindices;
    std::vector<double> fracvalues;
    // Check how many fractionals do we have.
    // Iterators for nonzero values.
    std::vector<double>::iterator it;
    std::vector<double>::iterator begin = values.begin();
    std::vector<double>::iterator end   = values.end();
    UInt index = 0;
    for (it=begin; it!=end; ++it) {
      if ((0.00 != *it) && (*it != 1.0)) {
        fracindices.push_back(indices[index]);
        fracvalues.push_back(*it);
      }
      index += 1;
    }
    
    // Check how many fractionals do we have.
    switch(fracindices.size()) {
    case 0: 
      gamma = gammainitial;
      return gamma;
      break;
    case 1: 
      gamma = gammainitial - fracvalues[0]*solution[fracindices[0]];
      return gamma;
      break;
    case 2: 
    {
      // Get row index of first fractional variable.
      UInt indexfirst = fracindices[0];
      VariablePtr var1 = *(prob->varsBegin() + indexfirst);
      UInt numrows = var1->getNumCons();
      UInt rowindex1 = 0;
      ConstraintPtr cons1 = *(var1->consBegin());
      if (numrows == 1) {
        rowindex1 = cons1->getIndex();
      }
      LinearFunctionPtr lf = cons1->getLinearFunction(); 
      //double coeff1 = lf->getWeight(var1);
     
 
      // Row index of second fractional.
      UInt indexsecond = fracindices[1];
      VariablePtr var2 = *(prob->varsBegin() + indexsecond);
      UInt numrows2 = var2->getNumCons();
      UInt rowindex2 = 0;
      ConstraintPtr cons2 = *(var2->consBegin());
      if (numrows2 == 1) {
        rowindex2 = cons2->getIndex();
      }
      LinearFunctionPtr lf2 = cons2->getLinearFunction();
      // double coeff2 = lf->getWeight(var2);
      
      if (rowindex1 != rowindex2) {
        cerr << "LGCIGenerator.cpp:roundHeur. Two fractional values are not in the same constraint. " << endl;
      } 

      // Serdar this is wrong check again.
      //double kfirstcoeff = (*prob->consBegin())->getLinearFunction()->getWeight(var1);
      //double ksecondcoeff =(*prob->consBegin())->getLinearFunction()->getWeight(var2);
      double kfirstcoeff = 1.0;
      double ksecondcoeff = 1.0;

      if (kfirstcoeff <= ksecondcoeff) {
        return (gammainitial - (ksecondcoeff*fracvalues[0]) + (kfirstcoeff*(1-fracvalues[0])) ); 
      } else {
        return (gammainitial + (ksecondcoeff*(1-fracvalues[0])) + (kfirstcoeff*fracvalues[0]) ); 
      }
      break;
    }
    default:
      cerr << "LGCIGenerator.cpp: More than two fractional values in the solution." << endl;
    }

  }else if (status == ProvenInfeasible){
    cerr << "LGCIGenerator.cpp::roundHeur. Lifting problem is infeasible." << endl;
    return 0;
  } else {
    cerr << "LGCIGenerator.cpp::roundHeur. Lifting problem is not solved to optimality." << endl;
  }
  
  // Serdar This should never return anything, make the upper part as if-else etc. a
  // complete structure.
  return 0.0;
}


void LGCIGenerator::initGubCons(const ConstCoverSetPtr cone,
                                boost::shared_ptr<std::vector<CoverSetPtr> > gubcoverlist, 
                                boost::shared_ptr<std::vector<CoverSetPtr> > gubcons)
{
  // Generate empty GUB constraints.
  UInt numgubs = gubcoverlist->size();
  for (UInt index=0; index<numgubs; ++index) {
    CoverSetPtr gubconstraint = (CoverSetPtr) new CoverSet();
    gubcons->push_back(gubconstraint);
  }

  // Iterators for each element of C1.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = cone->begin();
  CoverSetConstIterator end   = cone->end();
  // Iterators for unchanged GUB constraints.
  std::vector<CoverSetPtr>::iterator itgublist;
  std::vector<CoverSetPtr>::iterator begingublist = gubcoverlist->begin();
  std::vector<CoverSetPtr>::iterator endgublist   = gubcoverlist->end();
  // Iterate through each element of C1.
  for (it=begin; it!=end; ++it) {
    // Iterate through GUB constraint.
    UInt whichcons = 0;
    for (itgublist=begingublist; itgublist!=endgublist; ++itgublist) {
      // Iterators for each element of a GUB constraint.
      CoverSetIterator itgubel;
      CoverSetIterator begingubel = (*itgublist)->begin();
      CoverSetIterator endgubel   = (*itgublist)->end();
      // Iterate through each element of constraint.
      // Check if element is included in the current GUBconstraint
      for (itgubel=begingubel; itgubel!=endgubel; ++itgubel) {
        if(it->first->getId() == itgubel->first->getId()) {
          (*gubcons)[whichcons]->push_back(*itgubel);
        }
      }
      // end for elements of a GUB constraint. 
      whichcons += 1;
    } // end for GUB constraints.
  } // end for C1.

  // Print the initialized GUBs for debugging.
  if (DEBUG_LEVEL >= 9) {
    // Iterators for the GUBs.
    std::vector<CoverSetPtr>::const_iterator itdebug;
    std::vector<CoverSetPtr>::const_iterator begindebug = gubcons->begin();
    std::vector<CoverSetPtr>::const_iterator enddebug   = gubcons->end();
    // Iterate through all GUB constraints.
    bool first = true;
    for (itdebug=begindebug; itdebug!=enddebug; ++itdebug) {
      if (first) {
        printIneq((*itdebug), 1, Cons, "GUB constrants for the initial problem.");
        first = false;
      } else {
        printIneq((*itdebug), 1, Cons, "");
      }
    } // end of for.
  } // end of debug.
}

void LGCIGenerator::generateGubMaps(ConstConstraintVectorPtr gublist, 
                                    boost::shared_ptr< std::vector<CoverSetPtr> > gubcoverlist, 
                                    boost::shared_ptr< std::vector<CoverSetPtr> > guborigcoeffs)
{
  // Number of GUB constraints to be converted to coversets.
  ConstConstraintIterator it;
  ConstConstraintIterator begin = gublist->begin();
  ConstConstraintIterator end   = gublist->end();
  for (it=begin; it!=end; ++it) {
    // Construct cover set for coreesponding gub.
    CoverSetPtr guborigcoeff = (CoverSetPtr) new CoverSet();
    // Here both of these cover sets should be the same.
    // Copy one of them to the other one.
    varCoeff(*it,guborigcoeff);
    guborigcoeffs->push_back(guborigcoeff);
    CoverSetPtr gubcoeffones = (CoverSetPtr) new CoverSet();
    initCoverIneq(guborigcoeff, gubcoeffones);
    gubcoverlist->push_back(gubcoeffones);
  }
  
  // In debug mode, print out the GUB list generated.
  if (DEBUG_LEVEL >= 9) {
    std::vector<CoverSetPtr>::const_iterator it2;
    std::vector<CoverSetPtr>::const_iterator begin2 = guborigcoeffs->begin();
    std::vector<CoverSetPtr>::const_iterator end2   = guborigcoeffs->end();
    bool  first = true;
    // Print out original gub lists generated.
    for (it2=begin2; it2!=end2; ++it2) {
      if (first) {
        printIneq((*it2), 1, Cons, "GUB list original");
        first = false;
      } else {
        printIneq((*it2), 1, Cons, "");
      }
    } // end of for.
    // Print out gub listts that will be modified to be non-overlapping later.
    first = true;
    for (it2=gubcoverlist->begin(); it2!=gubcoverlist->end(); ++it2) {
      if (first) {
        printIneq((*it2), 1, Cons, "GUB list to be made non-overlapping.");
        first = false;
      } else {
        printIneq((*it2), 1, Cons, "");
      }
    } // end of for
  } // end of debug.

}

void LGCIGenerator::nonOverlap(boost::shared_ptr< std::vector<CoverSetPtr> > gubcoverlist)
{
  // Iterators for the gubs in the list.
  std::vector<CoverSetPtr>::iterator it;
  std::vector<CoverSetPtr>::iterator begin = gubcoverlist->begin();
  std::vector<CoverSetPtr>::iterator end   = gubcoverlist->end();
  // We iterate through each set.
  for (it=begin; it!=end; ++it) {
    CoverSetIterator itrefset;
    CoverSetIterator beginrefset = (*it)->begin();
    CoverSetIterator endrefset   = (*it)->end();
    // We iterate through each element of each set.
    for (itrefset=beginrefset; itrefset!=endrefset; ++itrefset) {
      // We iterate through each remaining set.
      std::vector<CoverSetPtr>::iterator itremset;
      for (itremset=(it+1); itremset!=end; itremset++) {
        // Iterate through each element of remaining set.
        CoverSetIterator itremsetel;
        CoverSetIterator beginremsetel = (*itremset)->begin();
        CoverSetIterator endremsetel   = (*itremset)->end();
        for (itremsetel=beginremsetel; itremsetel!=endremsetel; ++itremsetel) {
          // Check if element is included in this set.
          if (itremsetel->first->getId() == itrefset->first->getId()) {
            (*itremset)->erase(itremsetel);
          }
        } // end each element of a remaining set.
      } // end remaining set loop
    } // end each element loop.
  } // end each set loop
  
  // In debug mode, print out the non-overlapping GUB list.
  if (DEBUG_LEVEL >= 9) {
    std::vector<CoverSetPtr>::const_iterator it2;
    std::vector<CoverSetPtr>::const_iterator begin2 = gubcoverlist->begin();
    std::vector<CoverSetPtr>::const_iterator end2   = gubcoverlist->end();
    bool first = true;
    // Print out original gub lists generated.
    for (it2=begin2; it2!=end2; ++it2) {
      if (first) {
        printIneq((*it2), 1, Cons, "GUB list non-overlapping.");
        first = false;
      } else {
        printIneq((*it2), 1, Cons, "");
      }
    } // end of for
  } // end of debug.


}


void LGCIGenerator::initCoverIneq(const ConstCoverSetPtr cover,
                                  CoverSetPtr ineq)
{
  // Iterators for cover inequality.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = cover->begin();
  CoverSetConstIterator end   = cover->end();
  // Iterate through all variables of cone or any other initial set of variables.
  for (it=begin; it!=end; ++it) {
    // New variable for inequality.
    VariableValuePair newvar(*it);
    // Change the coefficient as one.
    newvar.second = 1.0;
    // Add variable to inequality.
    ineq->push_back(newvar);
  }
}

/**
   This function generates set N-C, the variables outdide of cover set.
 */
void LGCIGenerator::cBar(const ConstCoverSetPtr cover,
                         CoverSetPtr cbar,
                         const ConstConstraintPtr cons)
{
  // Now, define cbar the set of N-C.
  // Get the linear function.
  LinearFunctionPtr lf = cons->getLinearFunction();
  // Take a copy of cover set.
  CoverSetPtr covertodelete = (CoverSetPtr) new CoverSet(*cover);
  // Iterators for whole constraint.
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  // Add the variables which are not in the cover set.
  // Current variable of linear function.
  ConstVariablePtr curlfvar;
  // Current variable of cover set.
  ConstVariablePtr curcovvar;
  for (it=begin; it!=end; ++it) {
    bool in = false;
    curlfvar = it->first;
    // Iterators for cover elements, this changes dynamically thus it is
    // inside the loop.
    CoverSetIterator itcov;
    CoverSetIterator begincov = covertodelete->begin();
    CoverSetIterator endcov   = covertodelete->end();
    for (itcov=begincov; itcov!=endcov; ++itcov) {
      curcovvar = itcov->first;
      // Check if the element is in cover set.       
      // For now, this comparison works but may be we should check variable
      // id.
      // Pointers may be changed later. 
      // Assumption is Minotaur will not change these pointer values. 
      if (curlfvar == curcovvar) {
        in = true;
        // take out the element of copy of cover.
        covertodelete->erase(itcov);
        // Go to the new variable in linear function.
        break;
      }
    }
    // If variable is not in cover C then it is cbar N\C.
    if (in == false) {
      cbar->push_back(*it);
    }

  }
}

void LGCIGenerator::minimalCover(CoverSetPtr cover, ConstConstraintPtr cons)
{
  // Sort in nonincreasing order of coefficients.
  sortNonIncreasing(cover);
  // Initialize loop to make cover minimal (Take out elements).
  // The partial summation of coefficients.
  double sum = 0.0;
  // Right hand side of knapsack constraint considered.
  double b = cons->getUb();
  // Difference between the sum and b.
  double difference = 0.0;
  // Reverse iterartor used to get the minimum coefficient variable from the
  // end of coverset vector.
  CoverSet::reverse_iterator minpair;
  // minimum coefficient value.
  double min = 0.0;
  // Current sum is the summation of coefficients of all variables in cover.
  sum = sumCoeffs(cover);
  // Loop that removes the variables until minimal cover is obtained.
  do {
    difference = sum-b;
    minpair = cover->rbegin(); // Reverse iterator!
    min = minpair->second;
    if (min < difference) {
      cover->pop_back();
      sum -= min;
    } else {
      break;
    }
  } while(sum > b); // make sure that we do not erase variables such that
                    // sum=<b
  // Print out the minimal cover.
  if (DEBUG_LEVEL >= 9) {
    printIneq(cover, cover->size()-1, Cover, "Minimal cover cut for GNS.");
  } 
}

double LGCIGenerator::sumCoeffs(CoverSetPtr cover)
{
  // Iterators for the vector considered.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = cover->begin();
  CoverSetConstIterator end   = cover->end();
  // Coefficient of current variable.
  double coeff = 0.0;
  // Summation of coefficients in the vector.
  double sum = 0.0;
  for (it=begin; it!=end; ++it) {
    coeff = it->second;
    sum += coeff;
  }
  
  return sum;
}

void LGCIGenerator::varCoeff(ConstConstraintPtr cons, CoverSetPtr cover)
{
  // Iterators for the constraint considered.
  LinearFunctionPtr lf = cons->getLinearFunction();
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = cons->getLinearFunction()->termsEnd();
  // Iterate through elements of constraints.
  for (it=begin; it!=end; ++it) {
    VariableValuePairPtr currentvar = (VariableValuePairPtr) new VariableValuePair();
    currentvar->first  = it->first;
    currentvar->second = it->second;
    cover->push_back(*currentvar);
  }
}


// Add the cut to the lits and insert the violation of cut to violation list.
// Adds the cut to the corresponding list without checking anything else
// beside violation.
void LGCIGenerator::addCut(CutPtr cut)
{
  // Add cut to the cut list.
  cutVec_.push_back(cut);
  // Calculate violation.
  double viol = violation(cut);
  // If cut is violated by current relaxation solution, increment the number
  // of violated cuts.
  if (viol > 0) {
    if (DEBUG_LEVEL >= 9) {
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "Cut is a violated cut." << endl;
      output_ << "Violation of the cut is " << viol << endl;
      output_.close();
    }

    // Add the cut to the list of violated cuts.
    violatedCuts_.push_back(cut);
    // Increment the number of cuts violating current relaxation solution.
    stats_->violated += 1;
  } else {
    // Cut is not violated.
    if (DEBUG_LEVEL >= 9) {
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "Cut is not violated." << endl;
      output_.close();
    }
  }
  // Add the corresponding violation to violation list.
  violList_.push_back(viol);
}

bool LGCIGenerator::addCut(CoverSetPtr cov, double rhs, UInt cuttype, CutFail&)
{
  // In debug mode we write all cuts generated.
  if (DEBUG_LEVEL >= 9) {
    // These lines coverts the number to a string.
    int numcuts = stats_->totalcuts + 1;
    string numc;
    ostringstream convert;
    convert << numcuts;
    numc = convert.str();
    // until here.
    
    switch (cuttype) {
    case Gns:
      printIneq(cov, rhs, Cons, "Cut " + numc + "gns"); break;
    default:
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "LGCIGenerator.cpp::addCut. Undefined cut type." << endl;
      output_.close();
    } // end of switch.
  }

  // Total number of cuts generated is increased by one.
  stats_->totalcuts += 1;
  // Check violation of the coefficients.
  bool cutexists = checkExists(cov, rhs);
  if (cutexists == false){
    CutPtr cut;
    generateCut(cov, rhs, cut);
    addCut(cut);
    // Increment the number of cuts for the type of cuts generated.
    stats_->cuts += 1;
    switch (cuttype) {
    case Gns: stats_->gns += 1; break;
    default:
      output_.open(outfile_.c_str(), std::ios_base::app);
      output_ << "LGCIGenerator.cpp::addCut. Undefined cut type." << endl;
      output_.close();
    }
    // Cut is not duplicate.
    return true;
  } else { // Cut exists, it is a duplicate.
    return false;
  } 

}

// Assumption a cut pointer is given and the cut will be returned by this pointer.
void LGCIGenerator::generateCut(const ConstCoverSetPtr inequality, double rhs, CutPtr cut)
{
  // Generate linear function for cut.
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  // create iterators for cover set elements.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = inequality->begin();
  CoverSetConstIterator end   = inequality->end();
  // One by one add the variables and their coefficients.
  for (it=begin; it!=end; ++it) {
    // If coefficients if zero, do not add to the cut.
    if (it->second != 0) {
      lf->addTerm(it->first, it->second);
    }
  } // end of for.
  
  // generate function.
  FunctionPtr f = (FunctionPtr) new Function(lf);
  // Assumption: coefficients are positive lower bound is zero.
  double lb = 0.0;
  // Generate cover cut.
  cut = (CutPtr) new Cut(p_->getNumVars(), f, lb, rhs, false, false);
}

/**
 * TODO: Create a new index for variables in the knapsack only.
 * Do not consider all variables for coefficients.
 */
bool LGCIGenerator::checkExists(CoverSetPtr inequality, double rhs)
{
  // Iterators for variables in cut.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = inequality->begin();
  CoverSetConstIterator end   = inequality->end();
  // Number of variables in the problem.
  // Later change this to number of variables in the linear function.
  UInt numvars = p_->getNumVars();
  // Vector to include coefficients.
  std::vector<double> coeffs(numvars,0);
  // Index of the variable in the problem.
  UInt index = 0;
  // Current variable considered.
  VariablePtr var;
  // Coefficient divided by rhs of cover inequality.
  double dividedcoeff;
  // Add coefficients to coefficient vector one by one.
  for (it=begin; it!=end; ++it) {
    var           = it->first;
    index         = var->getIndex();
    dividedcoeff  = it->second / rhs;
    coeffs[index] = dividedcoeff;
  }
  
  // Check if the cut already exists.
  std::map< std::vector<double>, UInt>::iterator found  = cutmap_.find(coeffs);
  std::map< std::vector<double>, UInt>::iterator endmap = cutmap_.end();
  if (found == endmap) {
    cutmap_.insert(std::pair< std::vector<double>, UInt > (coeffs,stats_->totalcuts));
    return false; 
  } else {    
    return true;
  }
}

double LGCIGenerator::violation(CutPtr cut)
{
  FunctionPtr f = cut->getFunction();
  const double * x = s_->getPrimal();
  int error = 0;
  double evaluation = f->eval(x, &error);
  double violub     = max(0.0, evaluation - cut->getUb());
  double viollb     = max(0.0, cut->getLb() - evaluation);
  double viol       = max(violub, viollb);

  return viol;
}


void LGCIGenerator::printIneq(const ConstCoverSetPtr cover, double rhs, 
                                  PrintType type, string message)
{
  // Iterators for variables.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = cover->begin();
  CoverSetConstIterator end   = cover->end();
  
  // Current variable.
  ConstVariablePtr var;
  // Coefficient of variable.
  double coeff;
  // Name of variable.
  string name;
  // Close the file in any case it is open.
  output_.close();
  // Reopen the file.
  output_.open(outfile_.c_str(), std::ios_base::app);
  output_.precision(1);
  output_ << endl;

  // Print out the initial message.
  if (message.empty() == false) {
    output_ << message << endl;
  }
  
  // Check  if the given cover, constraint, objective or set is empty or not.
  if (cover->empty()) {
    switch(type) {
    case Cover: 
      output_ << "Empty cover."      << endl; break;
    case Cons:
      output_ << "Empty constraint." << endl; break;
    case Obj:
      output_ << "Empty objective."  << endl; break;
    case Set:
      output_ << "Empty set."        << endl; break;
    default:
      output_ << "CoverCutGenerator.cpp::printIneq. Invalid print type." << endl;
    } // end of switch.
  } else { // end of if.
    bool first = true;
    if (type == Obj) {
      output_ << "max " << std::flush;
    }

    for (it=begin; it!=end; ++it) {
      var   = it->first;
      coeff = it->second;
      name  = var->getName();
      switch (type) {
      case Obj: // Same for objective and constraint.        
      case Cons:
        if (true/*coeff != 0*/) {
          if (first == true) { // The first variable is different to write.
            output_ << coeff << "*" << name << std::flush;           
          } else {
            if (coeff >= 0) {
              output_ << " + ";
            }
            output_ << coeff << "*" << name  << std::flush;
          }
        }
        break;
      case Cover:
        if (first == true) {
          output_ << name << std::flush;
        } else {
          output_ << " + " << name << std::flush;
        }
        break;
      case Set:
        if (first == true) {
          output_ << name << std::flush;
        } else {
          output_ << ", " << name << std::flush;
        }
        break;
      default:
        output_ << "CoverCutGenerator.cpp::printIneq. Invalid print type." << endl;
      } // end of switch.
      first = false;
    } // end of for.
    // Print the sense and rhs of constraint.
    if (type == Cover || type == Cons) {
      // Assumption: All constraints etc. are in sense "<=".
      output_ << " <= " << rhs << endl;
    }
  } // end of else. 

  // Close output file.
  output_.close();
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
