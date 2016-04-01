//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2009 - 2014 The MINOTAUR Team.
//

/**
 * \file CoverCutGenerator.cpp
 * \brief Define base class CoverCutGenerator.
 * \author Serdar Yildiz, Argonne National Laboratory
 */
#include <cmath>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <ostream>
#include <cstring>
#include <limits>
using std::numeric_limits;
#include <algorithm>
using std::max;
#include <string>
using std::string;
#include <sstream>
using std::ostringstream;

#include "CoverCutGenerator.h"
#include "LinearFunction.h"
#include "Function.h"
#include "Types.h"
#include "Variable.h"
#include "Relaxation.h"
#include "Option.h"


using namespace Minotaur;

# define DEBUG_LEVEL -1

CoverCutGenerator::CoverCutGenerator() {}

// Currently unused, probably removed later.
CoverCutGenerator::CoverCutGenerator(ProblemPtr , SolutionPtr , EnvPtr )
// : stats_(0)
{
  // Check if initialization is successful.
  //bool successinit = false;
  //successinit = initialize(p,s, env);
  //if (successinit) {
  //  generateAllCuts();
  // }
}

CoverCutGenerator::CoverCutGenerator(RelaxationPtr rel, ConstSolutionPtr sol, EnvPtr env)
// : stats_(0)
{
  //ProblemPtr p;
  //p = rel;
  //ProblemPtr p = boost::static_pointer_cast<Problem>(rel);
  bool successinit = false;
  successinit = initialize(rel,sol, env);
  if (successinit) {
    generateAllCuts();
  }
}

CoverCutGenerator::~CoverCutGenerator()
{
  // deallocate heap.
  if (stats_) {
    delete stats_;
  }
}

bool CoverCutGenerator::initialize(RelaxationPtr p, ConstSolutionPtr s, EnvPtr env)
{
  // Environment is obtained for parameters.
  env_ = env;
  // Integrality tolerance is assigned.
  intTol_ = env_->getOptions()->findDouble("int_tol")->getValue();
  // Objective value change tolerance.
  objtol_ = 1e-6;
  // If no cut generation occured this will return 0s for each statistics.
  stats_ = new CovCutGenStats();
  initStats();
  numCons_ = 0;
  // We cannot take a copy of these objects, since there is no copy constructor.
  p_ = p;
  s_ = s;
  generateKnapList();
   
  // Check if the solution is non-integer. If the given solution is an integer
  // feasible solution then there is no need to generate cuts.   
  bool integral = false;
  integral = checkIntegral(p, s);
  // If the solution is not integral, then generate cuts.
  if (integral == false) {
    // If debugging is active create an output file.
    if (DEBUG_LEVEL >= 0) {
      outfile_ = "Liftings.txt";
      output_= (OfstreamPtr) new ofstream();
      output_->open(outfile_.c_str());
      // Print out the current solution.
      VariableConstIterator it;
      VariableConstIterator begin = p->varsBegin();
      VariableConstIterator end   = p->varsEnd();
      VarVector solution; 
      for (it=begin; it!=end; ++it) {
        solution.push_back(*it);
      }
      //(*output_).fixed;
      (*output_).precision(1);
      (*output_) << "Primal given solution" << endl;
      s->writePrimal((*output_),&solution);
      (*output_) << "\nDual solution is written" << endl;;
      s->writeDual((*output_));
      output_->close();        
    }
    
    // Initialization is succesfull.
    return true;
  } else {
    // Initialization unsuccessfull. No need to generate cut since solution is
    // integral already.
    return false;
  }
}

void CoverCutGenerator::initStats()
{
  stats_->knaps = 0;
  stats_->cuts = 0;
  stats_->violated = 0;
  stats_->extended = 0;
  stats_->simple = 0;
  stats_->gns = 0;
  stats_->singlectwo = 0;
  stats_->basic = 0;
}

void CoverCutGenerator::generateKnapList()
{
  knapsackListPtr_ = (KnapsackListPtr) new KnapsackList(p_);
}

void CoverCutGenerator::generateAllCuts()
{
  if (knapsackListPtr_->getNumKnaps() >= 1) {
    ConstraintIterator it;
    ConstraintIterator begin = knapsackListPtr_->getListBegin();
    ConstraintIterator end   = knapsackListPtr_->getListEnd(); 
    for (it=begin; it != end; ++it) {
      // If debug option is given write the constraint considered to output
      // file.
      if (DEBUG_LEVEL >= 9) {
        output_->open(outfile_.c_str(), std::ios_base::app);
        // output_->open("/sandbox/sey309/Liftings.txt");
        (*output_) << "\n\n\nConstraint considered is:" << endl;
        (*it)->write(*output_);
        output_->close();
      }
      
      // Check if constraint is GUB.
      bool isGUB = GUB(it); 
      if (isGUB == true) {
        continue;
      }
      
      // Check if constraint has a cover.
      bool has = hasCover(it);
      if (has) {
        numCons_ += 1;
        generateCuts(*it);
      }
    }
    // Check if no constraint is considered for cut generation.
    // may be use aseert.
  }else {
    cerr << "CoverCutGenerator::generateAllCuts. No constraint used for cut generation.\n";
  }
}

// Check if it is a GUB. If it is a GUB, we do not generate any cover cuts
// by using it. May be we should eliminate such constraints as well.
// x_1 + x_2 + x_3 <= 5
bool CoverCutGenerator::GUB(ConstraintIterator itcons)
{
  // Constraint considered.
  ConstraintPtr cons = *itcons;
  double b = cons->getUb();
  // Iterators for variables in function.
  LinearFunctionPtr lf = cons->getLinearFunction();
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  // Check if all the coefficients divided by rhs is one.
  // Iterate through all variables in linear function.
  double division = 0.0;
  double gendivision = 0.0;
  bool first = true;
  for (it=begin; it!=end; ++it) {
    division = (it->second) / b; 
    if (first == true){
      gendivision = division;
      first = false;
    }

    if (division != gendivision) {
      // It is not similar to GUB.
      return false;
    }
  }

  // It is similar to GUB.
  return true;
}

bool CoverCutGenerator::hasCover(ConstraintIterator itcons)
{
  // Variable that shows if the constraint has a cover.
  bool has = false;
  // Constraint considered.
  ConstraintPtr cons = *itcons;
  double b = cons->getUb();
  // Iterators for variables in function.
  LinearFunctionPtr lf = cons->getLinearFunction();
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd(); 
  // Sum of coefficients of the variables. 
  double sum = 0.0;
  // Iterate through all the variables in the linear function.
  for (it=begin; it!=end; ++it) {
    // Check if the sum is greater than rhs, i.e. \sum a_j > b.
    sum += it->second;
    if (sum > b) {
      has = true;
      break;
    }
  }
  // No cover exists.
  return has;
}

void CoverCutGenerator::nonzeroVars(LinearFunctionPtr lf,
                                    CoverSetPtr nonzerovars,
                                    CoverSetPtr zerovars)
{
  // Iterators to iterate on variables in knapsack constraint.
  VariableGroupConstIterator it;
  const VariableGroupConstIterator begin = lf->termsBegin();
  const VariableGroupConstIterator end   = lf->termsEnd();
  
  // Fill in nonzeros and zeros vectors.
  // Iterate through all the variables in knapsack constraint.
  for (it=begin; it != end; ++it) {
    // Index of variable in solution.
    UInt varindex = (it->first)->getIndex();
    // Variable value in solution.
    double varvalue = s_->getPrimal()[varindex];
    // Serdar check this, this should take value of variable as second element.!!!!
    VariableValuePair currentvar(it->first,varvalue);
    if (varvalue != 0) { // add a tolerance value here, if needed.
      nonzerovars->push_back(currentvar);
    } else {
      zerovars->push_back(currentvar);
    }
  }
}

void CoverCutGenerator::sortNonIncreasing(CoverSetPtr nonzeros)
{
  /**
   * We sort in nonincreasing order the variables according to their values in
   * fractional solution.
   */
  CompareValueVariablePair compare; 
  sort(nonzeros->begin(), nonzeros->end(), compare);
}

void CoverCutGenerator::sortReducedCosts(CoverSetPtr & vars)
{
  if (vars->empty() == false) {
    // Iterators for variables.
    CoverSetConstIterator it;
    CoverSetConstIterator begin = vars->begin();
    CoverSetConstIterator end   = vars->end();
    // First construct an array of index-reduced cost pair.
    std::vector<id> ordered;
    // An array of reduced costs of variables. 
    const double * reducedcosts = s_->getDualOfVars();
    // The index of variable.
    UInt index = 0;
    // Reduced cost of variable.
    double reducedcost = 0.0;
    // Index of variable in vars vector.
    UInt varsindex = 0;
    for (it=begin; it!=end; ++it) {
      index = it->first->getIndex();
      reducedcost = reducedcosts[index];
      // vars index-reduced cost to be added.
      // Very important!!! negative of reduced cost is used to not to construct
      // another compare for sort.
      // Better to construct a separate compares.
      id indreduced(varsindex, -reducedcost);
      ordered.push_back(indreduced);
      varsindex += 1;
    }
  
    // Order the variables in nondecreasing order of reduced costs.
    CompareIntDouble compare;
    sort(ordered.begin(), ordered.end(), compare);
  
    // Construct the oredered vars vector.
    CoverSetPtr varssorted = (CoverSetPtr) new CoverSet();
    UInt numvars = vars->size();
    UInt sortedindex = 0;
    for(UInt i=0; i<numvars; ++i) {
      sortedindex = ordered.at(i).first;
      varssorted->push_back(vars->at(sortedindex));
   }
    //Change the vars as the sorted vars vector.
    vars = varssorted;
  }
}

CoverSetPtr CoverCutGenerator::varCoeff(LinearFunctionPtr lf)
{
  // The cover set to be returned.
  CoverSetPtr varcoeff = (CoverSetPtr) new CoverSet();
  // Iterators for the vector considered.
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  // Summation of coefficients in the vector.
  for (it=begin; it!=end; ++it) {
    // For now this cannot be taken out.
    VariableValuePairPtr currentvar = 
      (VariableValuePairPtr) new VariableValuePair();
    currentvar->first  = it->first;
    currentvar->second = it->second;
    varcoeff->push_back(*currentvar);
  }

  return varcoeff; 
}


CoverSetPtr CoverCutGenerator::coverSetGeneratorDefault(ConstConstraintPtr cons)
{
  const LinearFunctionPtr lf = cons->getLinearFunction();
  // Copy all variable-coeffs from linear function to a vector.
  CoverSetPtr vars = varCoeff(lf);
  // Sort in nonincreasing order of coeffs.
  sortNonIncreasing(vars);
  CoverSetPtr def = (CoverSetPtr) new CoverSet();
  
  // Iterators are prepared.
  CoverSetIterator it;
  CoverSetIterator begin = vars->begin();
  CoverSetIterator end   = vars->end();
  double b = cons->getUb();
  // Add variables until sum > b.
  double sum = 0;
  for (it=begin; it!=end; ++it) {
    // Check if cover obtained.
    if (sum > b) {
      break;
    }
    def->push_back(*it);
    sum += it->second;
  }

  if(DEBUG_LEVEL >= 9) {
    printIneq(def, def->size()-1, Set, "Default cover generated.");
  }

  return (def); 
}

/** Assumption the cover set is empty at the beginning.
 * Initial cpver will be generated from scratch.
 */
bool CoverCutGenerator::coverSetGeneratorGNS(ConstConstraintPtr cons,
                                                    CoverSetPtr cover)
{
  // This is the knapsack constraint to be covered.
  const LinearFunctionPtr lf = cons->getLinearFunction();
  // Obtain nonzero variables in the given solution.
  CoverSetPtr nonzerovars = (CoverSetPtr) new CoverSet();
  CoverSetPtr zerovars    = (CoverSetPtr) new CoverSet();
  nonzeroVars(lf, nonzerovars, zerovars);
  // Sort variables in nonincreasing order of their values in the given solution.  
  sortNonIncreasing(nonzerovars);

  // Cover is being obtained from knapsack inequality.
  // Iterators for nonzeros vector.
  CoverSetIterator it;
  CoverSetIterator begin = nonzerovars->begin();
  CoverSetIterator end   = nonzerovars->end();
  // The partial summation of variable coefficients in knapsack constraint.
  double sum = 0; 
  // right hand side "b"
  double b = cons->getUb();
  // Assumption is that ub is positive. Let's add an assert later.
  for (it=begin; it!=end; ++it) {
    // Add variables until sum exceeds right hand side "b".
    if (sum <= b) {
      double weight = lf->getWeight(it->first);
      VariableValuePair newpair(it->first, weight);
      cover->push_back(newpair);
      sum += lf->getWeight(it->first);
    } else { 
      break;
    }
  }
  
  // Check if we obtained an initial cover by using GNS method.
  if (sum <= b) {
    // No cover generated.
    // New
    if (DEBUG_LEVEL >= 9) {
      output_->open(outfile_.c_str(), std::ios_base::app);
      (*output_) << "No initial cover obtained from GNS." << endl;
      output_->close();
    }
    return false;
  } else {
    // Print the initial cover if debugging needed.
    if (DEBUG_LEVEL >= 9) {
      printIneq(cover,cover->size()-1, Cover, "GNS initial cover generated.");
    }

    // Cover generated.
    return true;
  }

}

// Not well tested.
CoverSetPtr CoverCutGenerator::coverSetGenGNSModified(ConstConstraintPtr cons)
{
  const LinearFunctionPtr lf = cons->getLinearFunction();
  // Obtain nonzero variables in the given solution.
  VariableValuePairVectorPtr nonzerovars = 
    (VariableValuePairVectorPtr) new VariableValuePairVector();
  VariableValuePairVectorPtr zerovars = 
    (VariableValuePairVectorPtr) new VariableValuePairVector();
  nonzeroVars(lf, nonzerovars, zerovars);
  // Sort variables in nonincreasing order of their values in the given solution.  
  sortNonIncreasing(nonzerovars);

  // Cover is being obtained from knapsack inequality.
  CoverSetPtr coverset = (CoverSetPtr) new CoverSet();
  // Iterators for nonzeros vector.
  VariableValuePairVectorConstIterator it;
  VariableValuePairVectorConstIterator begin = nonzerovars->begin();
  VariableValuePairVectorConstIterator end   = nonzerovars->end();
  // The summation of variable coefficients in knapsack constraint.
  int sum = 0; 
  // right hand side "b
  int b = (int) cons->getUb();
  // Assumption is that ub is positive. Let's add an assert later.
  for (it=begin; it!=end; ++it) {
    // Check if the right hand side of knapsack is exceeded.
    if (sum <= b) {
      coverset->push_back(*it);
      sum += (int) lf->getWeight(it->first);
    }
    else
      break;
  }

  // Iterators for zeros vector.
  VariableValuePairVectorConstIterator itzero;
  VariableValuePairVectorConstIterator beginzero = zerovars->begin();
  VariableValuePairVectorConstIterator endzero   = zerovars->end();
  // If sum <= b, then we continue to add more variables. 
  if (sum <= b) {
    for (itzero=beginzero; itzero!=endzero; ++itzero) {
      // Check if cover set obtained.
      if (sum <= b) {
        coverset->push_back(*it);
        sum += (int) lf->getWeight(it->first);
      }
    }
  }
 
  // Add an assesrt here to check if it is a cover. 
  return coverset;
}


double CoverCutGenerator::sumCoeffs(CoverSetPtr cover)
{
  // Iterators for the vector considered.
  CoverSetIterator it;
  CoverSetIterator begin = cover->begin();
  CoverSetIterator end   = cover->end();
  // Summation of coefficients in the vector.
  double sum = 0.0;
  for (it=begin; it!=end; ++it) {
    sum+=it->second;
  }
  return sum;
}

void CoverCutGenerator::minimalCover(CoverSetPtr cover, 
                                     ConstConstraintPtr cons)
{
  // Sort in nonincreasing order of coefficients.
  sortNonIncreasing(cover);  
  // Initialization of minimizer loop. 
  // variable that holds the summation of coefficients.
  double sum = 0.0; 
  // rhs of the knapsack constraint. 
  double b = cons->getUb();
  // difference between sum and b.
  double difference = 0.0;
  // Reverse iterator used to get the minimum coefficient of variable from
  // the end of nonincreasing ordered vector.
  CoverSet::reverse_iterator minpair;
  double min = 0.0; // minimum coefficient value.
  // Loop that removes the variables until the minimal cover is obtained.
  sum = sumCoeffs(cover);
  // Loop that reamoves the variables until minimal cover is obtained.
  do {
    difference = sum-b;
    minpair = cover->rbegin(); // Reverse iterator!
    min = minpair->second;
    if (min < difference) {
      cover->pop_back();
      sum = sum - min;
    } else {
      break;
    }
  }while(sum > b); // make sure that we do not erase variables such that sum
                    // =< b
}

// add implementation details to here.
// Assumes that the given cover set is minimal.
// Assumes that cone, ctwo, fset, and fbar are empty sets.
void CoverCutGenerator::coverPartitionGNS(const ConstConstraintPtr cons,
                                          const ConstCoverSetPtr cover,
                                          CoverSetPtr cone,
                                          CoverSetPtr ctwo,
                                          CoverSetPtr fset,
                                          CoverSetPtr rset)
{
  // Partition cover set C into C1 and C2.
  // Iterators for cover set.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = cover->begin();
  CoverSetConstIterator end   = cover->end();
  // Get the primal solution.
  const double * x = s_->getPrimal();
  // Since |C1| >= 1, we change it as |C1| >= 2. 
  // Serdar explain this, if only one then trivial liftings!
  // This is the number of remaining variables in cover for C1.
  UInt numremaining = cover->size();
  // Check if x_j^* == 1 or not.
  for (it=begin; it!=end; ++it) {
    UInt index = it->first->getIndex(); 
    // Check if |C1| >= 2.
    if (x[index] == 1 && (numremaining >= 3)) {
      ctwo->push_back(*it);
      numremaining -= 1;
    } else {
      cone->push_back(*it);
    }
  }
  
  // Partition set cbar (N-C) into sets F and R as described in the paper of
  // Gu Nemhasuser and Savelsbergh.
  // Construct N-C.
  CoverSetPtr cbar = (CoverSetPtr) new CoverSet();
  cBar(cover, cbar, cons);

  // Divide cbar into sets F and R.
  // Iterators for cbar.
  CoverSetConstIterator itcbar;
  CoverSetConstIterator begincbar = cbar->begin();
  CoverSetConstIterator endcbar   = cbar->end();
  // Check if x_j^* > 0 or not.
  for (itcbar=begincbar; itcbar!=endcbar; ++itcbar) {
    UInt index = itcbar->first->getIndex();
    if (x[index] > 0) {
      fset->push_back(*itcbar);
    } else {
      rset->push_back(*itcbar);
    }
  }

  // Print out the cover for debugging.
  if (DEBUG_LEVEL >= 9) {
    printIneq(cover, cover->size()-1, Set,"Minimal cover for GNS cover partition.\n C");
    printIneq(cbar, cbar->size()-1, Set,"N-C");
    printIneq(cone, cone->size()-1, Set,"C1");
    printIneq(ctwo, ctwo->size()-1, Set,"C2");
    printIneq(fset, fset->size()-1, Set,"F");
    printIneq(rset, rset->size()-1, Set,"R");
  }
}

/* This function generates set N\C, the variables outside of cover set. 
 */
void CoverCutGenerator::cBar(const ConstCoverSetPtr cover, 
                             CoverSetPtr cbar,
                             const ConstConstraintPtr cons)
{
  // Now, define cBar the set of N\C.
  LinearFunctionPtr lf = cons->getLinearFunction();
  // Iterators for whole constraint.
  VariableGroupConstIterator it;
  VariableGroupConstIterator begin = lf->termsBegin();
  VariableGroupConstIterator end   = lf->termsEnd();
  // Iterators for cover set.
  CoverSetConstIterator itcov;
  CoverSetConstIterator begincov = cover->begin();
  CoverSetConstIterator endcov   = cover->end();
  // Construct C bar when initialCover is obtained and update when minimal
  // cover is obtained !!! This is just for efficiency.

  // Delete the variables in the cover set.
  // This part can be made much more efficient. 
  for (it=begin; it!=end; ++it) {
    bool in = false;
    for (itcov=begincov; itcov != endcov; ++itcov) {
      if (it->first == itcov->first) {
        in = true;
        break;
      }
    }
    // Since variable is not in C then it is in C bar.
    if (in == false) {
      cbar->push_back(*it);
    }
  }  
}

/*
  We have to check if these coefficients are integer, however they can be
  noninteger as well !!!!!
 */


/* This function prepares the problem data for lifting problem solver.
   It updates the rhs of lifting problem, rhs of lifting inequality, 
 */
double CoverCutGenerator::lift(const ConstCoverSetPtr obj,
                               const ConstCoverSetPtr constraint,
                               const CoverSetConstIterator variable,
                               double & rhs,
                               double & initialb,
                               bool uplift)
{
  // Increment number of knapsacks solved.
  stats_->knaps += 1;

  // Number of variables.
  UInt n = obj->size();
  // Objective functions coefficients.
  double * c = new double[n];
  CoverSetConstIterator it;
  CoverSetConstIterator begin = obj->begin();
  CoverSetConstIterator end   = obj->end();
  UInt i = 0;
  for (it=begin; it!=end; ++it) {
    c[i] = it->second;
    i += 1;
  }
  
  // Constraint Coefficients.
  double * a = new double[n];
  CoverSetConstIterator itCoef;
  CoverSetConstIterator beginCoef = constraint->begin();
  CoverSetConstIterator endCoef   = constraint->end();
  i = 0;
  for (itCoef=beginCoef; itCoef!=endCoef; ++itCoef) {
    a[i] = itCoef->second;
    i += 1;
  }
  
  // First I have to order variables suitable for knapsack solver!!!.
  double * temp = new double[n];
  for(i = 0; i < n; ++i) {
    if (a[i] != 0) {
      temp[i] = c[i]/a[i];
    } else {
      temp[i] = 0;
    }
  }
 
  std::vector<id> ordered;
  for (i = 0; i < n; ++i) {
    ordered.push_back(id(i,temp[i]));
  }
  // order the vector
  // This can be made more efficient.
  CompareIntDouble compare;
  sort(ordered.begin(), ordered.end(), compare);
  double * tempa = new double[n];
  double * tempc = new double[n];

  for (i = 0; i<n; ++i) {
    // The index of variable in initial a and c vectors.
    int index = ordered[i].first;
    tempa[i] = a[index];
    tempc[i] = c[index];
  }
  // Now, we do not need a and c since we sopied them into tempa and tempc.
  delete [] a;
  delete [] c;
  a = tempa;
  c = tempc;

  // Solution.
  double gamma = 0.0;

  // Solution vector.
  int * x = new int[n];

  // Rhs value of constraint (b-a_i) where a_i is the variable to lift up.
  if (uplift == true) {
    double b = initialb - variable->second;
    // Solve binary knapsack solver.
    binaryKnapsackSolver(n,b,c,a,gamma,x);    
    // Alpha is obtained.
    double alpha = rhs-gamma;   
    // No need to update rhs of inequality for up-ifting.
    
    // If debug is activated write the lifting problem to output file.
    if (DEBUG_LEVEL >= 9) {
      printLiftProb(obj,constraint,variable,rhs,initialb,uplift,b,gamma,alpha);
    }

    // Clean-up, deallocate memory.
    delete [] tempa;
    delete [] tempc;
    delete [] temp;
    delete [] x;
 
    return alpha;

  } else {
    // Order of these statements are important.

    // Rhs value of constraint b.
    double  b = initialb + variable->second;
    // Solve binary knapsack solver.
    binaryKnapsackSolver(n,b,c,a,gamma,x);
    // ksi is obtained.
    double ksi = gamma-rhs;
    // Update rhs of inequality.
    rhs = rhs + ksi;
    // If debug is activated write the lifting problem to output file.
    if (DEBUG_LEVEL >= 9) {
      printLiftProb(obj,constraint,variable,rhs,initialb,uplift,b,gamma,ksi);
    }
    // Update initial bound of constraint.
    initialb += variable->second;
    
    // Clean-up, deallocate memory.
    delete [] tempa;
    delete [] tempc;
    delete [] temp;
    delete [] x;
  
    return ksi;
  
  } 
}

// Simple lifted cover cut
void CoverCutGenerator::simple(const ConstCoverSetPtr cover,
                               const ConstCoverSetPtr cbar,
                               const ConstConstraintPtr cons)
{
    // that is the upper bound of constraint.
    if ((cbar->empty() == false)) {
      double initialb = cons->getUb();
      CoverSetPtr consknap = (CoverSetPtr) new CoverSet(*cover);
      CoverSetPtr coverineq = (CoverSetPtr) new CoverSet();
      initCoverIneq(cover, coverineq);
      double rhs = coverineq->size() -1;
      CoverSetPtr obj = (CoverSetPtr) new CoverSet(*coverineq);
      liftSet(obj,consknap, cbar, coverineq, rhs, initialb, true);      
      addCut(coverineq, rhs, Simple);
    }
    
}

void CoverCutGenerator::allCTwo(const ConstCoverSetPtr cover,
                                const ConstCoverSetPtr cone,
                                const ConstCoverSetPtr cbar,
                                const ConstConstraintPtr cons)
{
  // Initialize C2 as an empty set.
  CoverSetPtr ctwo = (CoverSetPtr) new CoverSet();
  // initialize constraint.
  CoverSetPtr constraint = (CoverSetPtr) new CoverSet();
  initCoverIneq(cover,constraint);
  // Generate all possible |C2|=1.
  // Up lift Cbar, downlift C2.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = cone->begin();
  CoverSetConstIterator end   = cone->end();
  for (it=begin; it!=end; ++it) {
    ctwo->push_back(*it);
    // // Take out the element from cone
    CoverSetPtr conein = (CoverSetPtr) new CoverSet(*cover);
    CoverSetPtr lifted = (CoverSetPtr) new CoverSet(*constraint);

    CoverSetIterator itconein;
    CoverSetIterator beginconein = conein->begin();
    CoverSetIterator endconein   = conein->end();
    CoverSetIterator itlifted;
    CoverSetIterator beginlifted = lifted->begin();
    CoverSetIterator endlifted   = lifted->end();
    for (itconein = beginconein; itconein!=endconein; ++itconein) {
      if (itconein->first->getId() == it->first->getId()) {
        conein->erase(itconein);
        break;
      }
    }

    for (itlifted=beginlifted; itlifted!=endlifted; ++itlifted) {
      if(itlifted->first->getId() == it->first->getId()) {
        lifted->erase(itlifted);
        break;
      }
    }
    
    double initialb = cons->getUb() - it->second;
    double rhs = lifted->size() -1;
    CoverSetPtr obj = (CoverSetPtr) new CoverSet(*lifted);
    if (cbar->empty() == false) {
      // conein have to be updated.
      liftSet(obj, conein, cbar, lifted, rhs, initialb, true);
    }
    if (ctwo->empty() == false) {
      liftSet(obj,conein, ctwo, lifted, rhs, initialb, false);
    }
    addCut(lifted, rhs, Singlectwo);
    ctwo->pop_back();
  }
}

CutPtr CoverCutGenerator::generateCut(const ConstCoverSetPtr constraint, const double ub)
{
  // generate linear function from the cut.
  LinearFunctionPtr lf = (LinearFunctionPtr) new LinearFunction();
  // create iterators.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = constraint->begin();
  CoverSetConstIterator end   = constraint->end();
  // One by one add the variables and their coefficients.
  for (it=begin; it!=end; ++it) {
    // If coefficient is 0 do not add to the cut.
    if (it->second != 0) {
      lf->addTerm(it->first, it->second);
    }
  }
  // generate function.
  FunctionPtr f = (FunctionPtr) new Function(lf);
  // Since coefficients are positive lower bound is zero.
  double lb = 0.0;

  // Generate cover cut.
  CutPtr cut = (CutPtr) new Cut(p_->getNumVars(), f, lb, ub, false, false);
  
  return cut;
}

// Add the cut to the list and insert the violation of cut to violation list.
void CoverCutGenerator::addCut(CutPtr cut)
{
  // Add the cut to cut list.
  cutVec_.push_back(cut);
  // Calculate the violation
  double viol = violation(cut);
  // If cut is violated by solution of current relaxation, then increment
  // number of violated cuts.
  if (viol > objtol_) {
    // Add the cut to the list of violated cuts.
    violatedCuts_.push_back(cut);
    // Increment the number of cuts violating solution of current relaxation.
    stats_->violated += 1;
    // Add the corresponding violation to violation list.
    violList_.push_back(viol);
  }
  // // Add the corresponding violation to violation list.
  // violList_.push_back(viol);
}

bool CoverCutGenerator::addCut(CoverSetPtr cov, double rhs, UInt cuttype)
{
  // In debug mode we write all cuts generated.
  if (DEBUG_LEVEL >= 9) {
    // These lines converts the number to a string.
    int numcuts = stats_->totalcuts + 1;
    string numc;
    ostringstream convert;
    convert << numcuts;
    numc = convert.str();
    // until here
    
    switch (cuttype) {
    case Extended:   printIneq(cov,rhs,Cons,"Cut " +numc +" extended");   break;
    case Simple:     printIneq(cov,rhs,Cons,"Cut " +numc +" simple");     break;
    case Gns:        printIneq(cov,rhs,Cons,"Cut " +numc +" gns");        break;
    case Singlectwo: printIneq(cov,rhs,Cons,"Cut " +numc +" singlwctwo"); break;
    case Basic:      printIneq(cov,rhs,Cons,"Cut " +numc +" basic");      break;
    default:
      // New
      output_->open(outfile_.c_str(), std::ios_base::app);
      (*output_) << "CoverCutGenerator.cpp::addCut. Undefined cut type." << endl;
      output_->close();
    }
  }

  // Total number of cuts generated increased by one.
  stats_->totalcuts += 1;
  // I have to check violation and check integrality of coefficients as well.
  bool cutexists = checkExists(cov, rhs);
  if (DEBUG_LEVEL >= 10){
    cerr << "EXISTS = " << cutexists << endl; 
  }
  if (cutexists == false) {
    CutPtr cut = generateCut(cov, rhs);
    addCut(cut);
    // Increment the number of cuts for the type of cut generated.
    stats_->cuts += 1;
    switch(cuttype) {
    case Extended:    stats_->extended += 1;    break;
    case Simple:      stats_->simple += 1;      break;
    case Gns:         stats_->gns += 1;         break;
    case Singlectwo:  stats_->singlectwo += 1;  break;
    case Basic:       stats_->basic += 1;       break;
    default:
      if (DEBUG_LEVEL >= 10) {
        // New
        output_->open(outfile_.c_str(), std::ios_base::app);
        (*output_) << "CoverCutGenerator.cpp::addCut. Undefined cut type." <<  endl;
        output_->close();
      }
    }
    // New
    // Cut is not a duplicate.
    return true;  
  } else {
    // Cuty is duplicate.
    return false;
  }

}


// Careful! It checks values according to indices of variable not to ID's.
// This may cause a problem at the higher level when we check duplicacy in Cut manager.
bool CoverCutGenerator::checkExists(CoverSetPtr cov, double rhs)
{
  // Iterators for variables in cut.
  CoverSetIterator it;
  CoverSetIterator begin = cov->begin();
  CoverSetIterator end = cov->end();
  // Number of variables in the linear function.
  UInt numvars = p_->getNumVars();
  // Vector to include the coefficients.
  std::vector<double> coeffs(numvars,0);
  // Index of the variable in the problem.
  UInt index;
  // Current variable considered.
  VariablePtr var;
  // Coefficient divided by rhs of cover inequality.
  double dividedcoeff;
  // Add the coefficients to coefficient vector one by one.
  if (DEBUG_LEVEL >= 10) {
    cerr << "Coeffs: ";
  }
  for (it=begin; it!=end; ++it) {
    var = it->first;
    index  = var->getIndex();
    // Serdar check this int, what it does, probably it should not be here.
    dividedcoeff = double(it->second) / rhs;
    coeffs[index] = dividedcoeff;
    if (DEBUG_LEVEL >=  10) {
      cerr << dividedcoeff << " ";
    } 
  }
  if (DEBUG_LEVEL >= 10) {
    cerr << endl;
  }

  // Check if the cut already exists.
  std::map< std::vector<double>, UInt>::iterator found = cutmap.find(coeffs);
  std::map< std::vector<double>, UInt>::iterator endmap = cutmap.end();
  if (found == endmap) {
    cutmap.insert(std::pair< std::vector<double>, UInt > (coeffs,stats_->totalcuts));
    return false;
  } else {
    return true;
  }

}


/* This function generates various cover inequalities.
 */
void CoverCutGenerator::generateCuts(ConstConstraintPtr cons)
{
  // Add all these cuts to the list !!!.
  // Initial cover set is generated (possibly not mimimal cover).
  CoverSetPtr cover = coverSetGeneratorDefault(cons);
  // If it is needed minimal cover cut is added.
  CoverSetPtr basiccoverineq = (CoverSetPtr) new CoverSet();
  initCoverIneq(cover, basiccoverineq);
  double ub = cover->size() - 1;
  addCut(basiccoverineq, ub, Basic);
  // Drop some of the variables so that the cover set becomes minimal.
  minimalCover(cover, cons);
  if (DEBUG_LEVEL >= 9) {
    printIneq(cover,cover->size()-1,Cover,"Minimal cover for simple lifted cover inequality.");
  }
  // Requires a minimum cover.
  CoverSetPtr cone = (CoverSetPtr) new CoverSet(*cover);
  CoverSetPtr cbar = (CoverSetPtr) new CoverSet();
  cBar(cover, cbar, cons);
  
  // Those functions below assumes minimal cover as an input.
  // Generate a simple lifted cut.
  if (DEBUG_LEVEL >= 9) {
    printIneq(cbar, 0, Set, "N-C for simple lifted cover inequality.");
  }
  simple(cover, cbar, cons);

  // Different lifted cuts.
  // I guess this does not generate different cuts from initial cover inequalities.
  //allCTwo(cover, cone, cbar, cons);

  // Generates Gu, Nemhauser, Savelsbergh cover. 
  GNS(cons);
}


void CoverCutGenerator::initCoverIneq(const ConstCoverSetPtr coverset, 
                                      CoverSetPtr coverineq)
{
  // Iterators for cover inequality.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = coverset->begin();
  CoverSetConstIterator end   = coverset->end();
  for (it=begin; it!=end; ++it) {
    VariableValuePair newvar(*it);
    newvar.second = 1.0;
    coverineq->push_back(newvar);
  }
}

// Lifting strategy of Gu, Nemhauser and Savelsbergh.
// Add implementation details here.
// Assumption sets F, R, C1, and C2 are given. 
// They are already constructed as described in paper.
// However, vectors are not sorted or anything else.
bool CoverCutGenerator::liftingGNS(const ConstCoverSetPtr cone,
                                   const ConstCoverSetPtr ctwo,
                                   const ConstCoverSetPtr fset,
                                   const ConstCoverSetPtr rset,
                                   CoverSetPtr constraint,
                                   const ConstConstraintPtr cons,
                                   double & rhs)
  
{
  // Initialize the needed data elements.
  // rhs of cover inequality.
  rhs = cone->size() - 1;
  // Initial constraint which is basically cover inequality obtained from C1
  // by changing coefficients to 1s.
  initCoverIneq(cone,constraint);
  // Initial constraint of lifting problem which is exactly C1. 
  CoverSetPtr consknap = (CoverSetPtr) new CoverSet(*cone);
  // Initial objective of lifting problem which is exactly same as C1 with
  // coefficients 1s.
  CoverSetPtr obj = (CoverSetPtr) new CoverSet(*constraint);

  // Calculate the initial b.
  double initialb = cons->getUb();
  // Iterators for C2.
  CoverSetConstIterator it;
  CoverSetConstIterator begin = ctwo->begin();
  CoverSetConstIterator end   = ctwo->end();
  // Decrease initial b by the coefficients of C2.
  for (it=begin; it!=end; ++it) {
    initialb -= it->second;
  }

  bool liftup = true;
  // We are going to up-lift variables in set F. 
  if (fset->size() >= 1) { 
    liftup = true;
    // Uplift if set F is non-empty.
    liftSetF(obj, consknap, fset, constraint, rhs, initialb, liftup);
  }

  // If there is no violation, what to do?
  // is the current half lifted cut is valid for the problem?
  // Check if there is violation or not.
  CutPtr tempcut = generateCut(constraint,rhs);
  double viol = violation(tempcut);
  // If there is no violation stop the algorithm.
  if (viol == 0) {
    stats_->noviol += 1;
    return false;
  }

  if (ctwo->size() >= 1) {
    // Variables in set C2 are going to be down-lifted.
    liftup = false;
    // First order the variables in nondecreasing magnitude of reduced costs.
    CoverSetPtr sortedctwo = (CoverSetPtr) new CoverSet(*ctwo);
    sortReducedCosts(sortedctwo);
    // Lift down the variables.
    liftSet(obj, consknap, sortedctwo, constraint, rhs, initialb, liftup);
  }

  if (rset->size() >= 1) {
    // Variables in set R are going to be up-lifted.
    liftup = true;
    // First order the variables in nondecreasing magnitude of reduced costs.
    CoverSetPtr sortedrset = (CoverSetPtr) new CoverSet(*rset);
    sortReducedCosts(sortedrset);
    // Lift up the variables.
    liftSet(obj, consknap, sortedrset, constraint, rhs, initialb, liftup);
  }
  
  return true;
}

/* Assumptions: Following data elements should be suitable for lifting problem
 * construction.
 * Objective, lifting problem constraint, variable set to be lifted,
 * inequality to be lifted, rhs of inequality, the rhs of lifitng problem.
 * Under these assumptions, this function determine the coefficient of
 * lifted variable, 
 * adds the new variable to inequality being lifted, lifting problem objective
 * function, and lifting problem constraint.
 * obj, consknap, constraint, rhs, initialb are changed at the end.
 * varset and liftup does not change.
 */
void CoverCutGenerator::liftSet(CoverSetPtr obj,
                                     CoverSetPtr consknap,
                                     const ConstCoverSetPtr varset,
                                     CoverSetPtr constraint,
                                     double & rhs,
                                     double & initialb,
                                     bool liftup)
{
  // Check if the set being liftes is empty or not.
  if (varset->empty() == false) {
    // Iterators for the set of variables to be lifted.
    CoverSetConstIterator it;
    CoverSetConstIterator begin = varset->begin();
    CoverSetConstIterator end   = varset->end();

    // Initialize the coefficient of variable to be lifted.
    double alpha = 0.0; 
    // Lift the variables one by one in the given order.
    for (it = begin; it!=end; ++it) {
      // Lift the variable.
      alpha = lift(obj, consknap, it, rhs, initialb, liftup);
      // Construct new pair for objective.
      VariableValuePairPtr newobjvar =
        (VariableValuePairPtr) new VariableValuePair(*it);
      // Update coefficient of variable with alpha.
      newobjvar->second = alpha;
      // Update objective function.
      obj->push_back(*newobjvar);
      // Update cover inequality.
      constraint->push_back(*newobjvar);
      // Update knapsack problem constraint.
      consknap->push_back(*it);
    }
  }
}

/* This function uplifts the variables in the set F as described in Gu,
 * Nemhauser, Savelsbergh paper.
 * Assumptions: The following data elemenst are assumed to be suitable to
 * construct lifting problem:
 * ojective and constraint of lifting problem, variable set to be uplifted
 * (set F), lifting inequality, rhs of lifting inequality, rhs of lifting
 * problem constraint.
 * obj, consknap, constraint, rhs and initial b will change at the end.
 * varset and liftup remains the same.
 */
void CoverCutGenerator::liftSetF(CoverSetPtr obj,
                                 CoverSetPtr consknap, 
                                 const ConstCoverSetPtr setf,
                                 CoverSetPtr coverineq,
                                 double & rhs,
                                 double & initialb,
                                 const bool liftup)
{
  // If set F is empty, no need for lifting.
  if (setf->empty() ==  false) {
    // Get a copy of set F.
    CoverSetPtr fin = (CoverSetPtr) new CoverSet(*setf);
    // Get primal solution x^*.
    const double * x = s_->getPrimal();

    // Lift up the variables in the given set in a greedy way.
    while (fin->empty() == false) {
      // Iterators for given set.
      CoverSetIterator it;
      CoverSetIterator begin = fin->begin();
      CoverSetIterator end   = fin->end();
    
      // Initialize the best variable to be uplifted.
      CoverSetIterator bestvar = begin;
      double bestalpha = 0.0;
      double bestcontribution = 0.0;
      // At each iteration calculate the alpha for each element in F.
      for (it=begin; it!=end; ++it) {
        // Prepare parameters. 
        CoverSetPtr objlocal = (CoverSetPtr) new CoverSet(*obj);
        CoverSetPtr constknaplocal = (CoverSetPtr) new CoverSet(*consknap);
        CoverSetIterator varlocal = it;
        double initblocal = initialb;
        double rhslocal = rhs;
        // Check if lift changes anything in the problem that should not be
        // changed!!!.
        double alpha = lift(objlocal, constknaplocal, varlocal, rhslocal, initblocal, liftup);
        // Get x_j^* from current solution.
        UInt index = it->first->getIndex();
        double xjvalue = x[index]; 
        double contribution = alpha * xjvalue;
        // Check if the best variable changed.      
        // Select the variable in F with largest alpha_j*x_j.  
        if (contribution > bestcontribution) {
          bestvar = it;
          // Update best variable and best alpha information.
          bestalpha = alpha;
          bestcontribution = contribution;
        }
      }

      // Up-lift the best variable.
      // Construct the new variable-alpha pair for objective.
      VariableValuePairPtr newobjvar = 
        (VariableValuePairPtr) new VariableValuePair(*bestvar);
      // Change the coefficient to alpha.
      newobjvar->second = bestalpha;
      // Update objective function.
      obj->push_back(*newobjvar);
      // Update the cover inequality.
      coverineq->push_back(*newobjvar);
      // Update the constraint of the lifting problem.
      consknap->push_back(*bestvar);
      // Take out the best variable that is lifted from out of F.
      fin->erase(bestvar);
    }
  }
}


bool CoverCutGenerator::GNS(const ConstConstraintPtr cons)  
{
  // Initial cover is empty and filled later.
  CoverSetPtr cover = (CoverSetPtr) new CoverSet();
  bool covgenerated = coverSetGeneratorGNS(cons, cover);
  // Check if cover generated.
  if(covgenerated == false) {
    return false;
  }

  // Obtain the minimal cover.
  // Serdar print this out for debug as well!!!.
  minimalCover(cover, cons);

  CoverSetPtr cone = (CoverSetPtr) new CoverSet();
  CoverSetPtr ctwo = (CoverSetPtr) new CoverSet();
  CoverSetPtr fset = (CoverSetPtr) new CoverSet();
  CoverSetPtr rset = (CoverSetPtr) new CoverSet();
  coverPartitionGNS(cons, cover, cone, ctwo, fset, rset);

  // Cover inequality is initialized as empty.
  CoverSetPtr ineq = (CoverSetPtr) new CoverSet();
  double ub = 0.0;
  bool generated = liftingGNS(cone, ctwo, fset, rset, ineq, cons, ub);
  if (generated == true) {
    bool cutgen = addCut(ineq, ub, Gns);
    return cutgen;
  } else {
    stats_->noinitcov += 1;
    return false;
  }
}

// Not very well tested.
void CoverCutGenerator::extendedCover(CoverSetPtr cover, ConstConstraintPtr cons)
{
  const LinearFunctionPtr lf = cons->getLinearFunction();
  // The whole constraint is transferred to a vector.
  CoverSetPtr vars = varCoeff(lf);
  // Variables ordered for extending.
  sortNonIncreasing(vars);
  // This is the cover set to be extended.
  CoverSetPtr extended = (CoverSetPtr) new CoverSet(*cover);
  // Rhs of cover cut does not change.
  double ub = cover->size()-1;
  // Determine the variable that has the largest coefficient in cover set.
  CoverSetIterator begin = cover->begin();
  CoverSetIterator end   = cover->end();
  CompareValueVariablePair compare;
  // Since ordering rule is reverse, I do min, indeed this hsould be changed
  // to  max by creating a new compare type.
  CoverSetIterator maxincover = min_element(begin, end, compare);
  // Add each element that has higher coefficient than cover elements. 
  CoverSetIterator it;
  CoverSetIterator beginvars = vars->begin();
  CoverSetIterator endvars   = vars->end();
  // Here, we add the variables according to nonincreasing order of coefficients.
  // We can use CBAR, by this way it may work faster.
  double maxcoeff = maxincover->second;
  for (it=beginvars; it!=endvars; ++it) {
    double currcoeff = it->second;
    if (currcoeff > maxcoeff) {
      extended->push_back(*it);
    } else if (currcoeff == maxcoeff) {
      // If the coefficients are the same then check if the variables are the
      // same or not
      UInt idmaxcoeff  = maxincover->first->getId();
      UInt idcurrcoeff = it->first->getId();
      if (idcurrcoeff != idmaxcoeff) {
        extended->push_back(*it);
      }
    } else {
      break;
    }
  }

  addCut(extended, ub, Extended);
}


// This code is based on CglKnapsackCover Inequalities.
// Assumption: Items are sorted in the order of c_1/a_1 >= c_2/a_2 >= ... c_n/a_n  
UInt CoverCutGenerator::binaryKnapsackSolver(UInt n, double b, 
                                             double const * c, double const * a, 
                                             double & z, int * x)
{
  // Set the solution as a vector of zeros. // change this to memset.
  memset(x, 0, (n)*sizeof(int));
  // Current solution vector.
  UInt * xhat = new UInt[n+1];
  // Set the solution as a vector of zeros.
  memset(xhat,0,(n+1)*sizeof(int));
  UInt j = 0;
 
  // set up: adding extra elements
  double * cIn = new double[n+2];
  double * aIn = new double[n+2];
  UInt ii = 0;
  for (ii=1; ii<n+1; ii++) {
    cIn[ii]=c[ii-1];
    aIn[ii]=a[ii-1];
  }

  // 1. initialize
  double zhat = 0.0;
  z = 0.0;
  double bhat = b+0.000001;
  cIn[n+1] = 0.0;
  aIn[n+1] = numeric_limits<double>::infinity(); // put infinity here
  j = 1;

  while(1) {
    // 2. compute upper bound u
    // "find r = min {i:sum k=j, i a_k>bhat};"
    ii=j;
    double aSemiSum = aIn[j];
    double cSemiSum = cIn[j];
    while (aSemiSum <= bhat && ii<n+2) {
      ii++;
      aSemiSum += aIn[ii];
      cSemiSum += cIn[ii];
    }
    if (ii==n+2) {
      cout << "Exceeded iterator limit. Aborting...\n";
      abort();
    }
  
    // r = ii at this point
    aSemiSum -= aIn[ii];
    cSemiSum -= cIn[ii];
    double u = cSemiSum + floor((bhat - aSemiSum)*cIn[ii]/aIn[ii]);
    // "if (z >= zhat + u) goto 5: backtrack;"
    if (!(z >= zhat + u)) {
      do {
        // 3. perform a forward step.
        while(aIn[j] <= bhat) {
          bhat -= aIn[j];
          zhat += cIn[j];
          xhat[j]=1;
          j+=1;
        }
        if (j<=n) {
          xhat[j]=0;
          j+=1;
        }
      } while (j==n);
      // "if (j<n) goto 2: compute_ub;"
      if (j<n) {
        continue;
      }
      // 4. update the best solution so far.
      if (zhat > z) {
        z=zhat;
        UInt k;
        for (k=0; k<n; k++) {
          x[k]=xhat[k+1];
        }
      }
      j=n;
      if (xhat[n] == 1) {
        bhat += aIn[n];
        zhat -= cIn[n];
        xhat[n]=0;
      }
    }
    // 5. backtrack
    // "find i=max{k<j:xhat[k]=1};" 
    UInt i=j-1;
    while (!(xhat[i]==1)&& i>0) {
      i--;
    }
    // "if (no such i exists) return;"
    if (i==0) {
      // clean up, deallocate memory.
      delete [] cIn;
      delete [] aIn;
      delete [] xhat;

      return 1;
    }
    bhat += aIn[i];
    zhat -= cIn[i];
    xhat[i] = 0;
    j=i+1;
    // "goto 2: compute_ub;"
  }

}

/* This function evaluates the function of the cut given for the solution
 * given. Then, it calculates the violation to the lower bound and upper bound
 * of the constraint. 
 * Takes the maximum of these violations.
 * This violation measure is different form GNS, since they only consider
 * upper bound violation of cut.
 */
double CoverCutGenerator::violation(CutPtr cut)
{
  FunctionPtr f = cut->getFunction();
  const double * x = s_->getPrimal();
  int error = 0;
  double evaluation = f->eval(x, &error);
  double violub = max(0.0, evaluation - cut->getUb());
  double viollb = max (0.0, cut->getLb() - evaluation); 
  double viol = max(violub, viollb);
  return viol;
}

bool CoverCutGenerator::checkIntegral(RelaxationPtr p, ConstSolutionPtr s)
{
  // Iterators for variables in problem.
  VariableConstIterator it;
  VariableConstIterator begin = p->varsBegin();
  VariableConstIterator end   = p->varsEnd();

  // Primal solution.
  const double * x = s->getPrimal(); 
  // Current variable.
  ConstVariablePtr var;
  // Index of variable.
  UInt index = 0;
  // Value of variable.
  double value = 0.0;
  // Absolute of fractional part of variable value.
  double fraction = 0.0;
  // Type of variable.
  VariableType type;
  // Iterate through all variables in the problem.
  for (it=begin; it!=end; ++it) {
    var = (*it);
    type = var->getType();
    //Check if variable is type integer.
    if (type == Binary || type == Integer) {
      index = var->getIndex();
      value = x[index];
      fraction = fabs(value - floor(value+0.5));
      if (fraction > intTol_) {
        // Check if the variable has fractional value.
        return false;
      }
    }
  }
  
  return true;
}


void CoverCutGenerator::printIneq(const ConstCoverSetPtr cover, double rhs, PrintType type, string message)
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
  output_->close();
  output_->open(outfile_.c_str(), std::ios_base::app);
  //(*output_).fixed;
  (*output_).precision(1); 
  (*output_) << endl;

  // Print out initial message.
  if (message.empty() == false) {
    (*output_) << message << endl;
  }  
  
  // Check if the given cover, constraint, objective or set is empty or not.
  if (cover->empty()) { 
    switch (type) {
    case Cover:
      (*output_) << "Empty cover."      << endl; break;
    case Cons:
      (*output_) << "Empty constraint." << endl; break;
    case Obj:
      (*output_) << "Empty objective."  << endl; break;
    case Set:
      (*output_) << "Empty set."        << endl; break;
    default:
      (*output_) << "CoverCutGenerator.cpp::printIneq. Invalid print type." << endl; 
    }  
  } else {
    // Iterate through all variables in the inequality.
    bool first = true;
    if (type == Obj) {
      (*output_) << "max " << std::flush;
    }
    for (it=begin; it!=end; ++it) {
      var   = it->first;
      coeff = it->second;
      name  = var->getName();
      switch (type) {
      case Obj: // Same for objective and constraint.
      case Cons:
        if (true/*coeff != 0*/) {
          if (first == true) {
            (*output_) << coeff << "*" << name << std::flush;
          } else {
            (*output_) << " + " << coeff << "*" << name << std::flush;
          }
        }
        break;
      case Cover:
        if (first == true) {
          (*output_) << name  << std::flush;
        } else {
          (*output_) << " + " << name << std::flush;
        }
        break;
      case Set:
        if (first == true) {
          (*output_) << name << std::flush;
        } else {
          (*output_) << ", " << name << std::flush;
        }
        break;
      default:
        cerr << "CoverCutGenerator.cpp::printIneq. Invalid print type." << endl; 
      } // end of switch.
      first = false;
    } // end of for.
    // Print the sense and rhs of constraint.
    if (type == Cover || type == Cons) {
      (*output_) << " <= " << rhs << endl;
    }
  } // end of else.
   
  // Close output file.
  output_->close();
}

void CoverCutGenerator::printLiftProb(const ConstCoverSetPtr obj,
                                      const ConstCoverSetPtr consknap,
                                      const CoverSetConstIterator variable,
                                      double rhs,
                                      double ,
                                      bool uplift,
                                      double b,
                                      double gamma,
                                      double alpha)
{
  output_->close();
  output_->open(outfile_.c_str(), std::ios_base::app);
  // Add two lines empty.
  (*output_) << endl << endl;;
  // Print out for uplifting.
  if (uplift) {
    (*output_) << "Uplifting variable " << variable->first->getName() 
               << " with coefficient " << variable->second << endl; 
  } 
  
  if (uplift == false) {
    (*output_) << "Downlifting variable " << variable->first->getName() 
               << " with coefficient " << variable->second << endl;;
  }
  
  // Objective is printed.
  (*output_) << "Objective max:" << endl;
  output_->open(outfile_.c_str(), std::ios_base::app);  
  printIneq(obj,0,Obj,"");
  // Knapsack constraint.
  output_->open(outfile_.c_str(), std::ios_base::app);  
  printIneq(consknap, b, Cons,"");

  //Solution of lifting problem.
  output_->open(outfile_.c_str(), std::ios_base::app);  
  (*output_) << "\nObj. value (gamma) = " << gamma << endl;
  (*output_) << "Right hand side of inequality = " << rhs << endl;
  (*output_) << "Alpha = " << alpha << endl;
  
  // print the inequality
  CoverSetPtr ineq = (CoverSetPtr) new CoverSet(*obj);
  VariableValuePair newvar(*variable);
  newvar.second = alpha;
  ineq->push_back(newvar);
  output_->open(outfile_.c_str(), std::ios_base::app);  
  printIneq(ineq, rhs, Cons, "One variable more uplifted cover inequality.");
  
  output_->close();
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
