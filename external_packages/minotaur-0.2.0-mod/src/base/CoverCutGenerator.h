//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file CoverCutGenerator.h 
 * \brief Declare base class CoverCutGenerator. 
 * \author Serdar Yildiz, Argonne National Laboratory 
 */

#ifndef MINOTAURCOVERCUTGENERATOR_H
#define MINOTAURCOVERCUTGENERATOR_H

#include <map>
#include <fstream>
using std::ofstream;
#include <string>
using std::string;

#include "KnapsackList.h"
#include "Problem.h"
#include "Solution.h"
#include "Types.h"
#include "Cut.h"
#include "Relaxation.h"
#include "Environment.h"

namespace Minotaur {

  typedef enum {
    Cover = 0,
    Cons,
    Obj,
    Set
  } PrintType;

  typedef enum {
    Totalcuts = 0,
    Cuts,
    Violated,
    Extended,
    Simple,
    Gns,
    Singlectwo,
    Basic,
    Noviol,
    Noinitcover
  } KnapCovType;

  struct CovCutGenStats
  {
    UInt knaps; /// Number of total knapsacks solved.
    UInt totalcuts; /// Number of all cuts i.e. included duplicates as well.
    UInt cuts;  /// Number of total cover cuts generated.
    UInt violated; /// Number of violated cuts.
    UInt extended; /// Number of extended cuts generated.
    UInt simple; /// Number of simple lifted cover cuts generated.
    UInt gns; /// Number of Gu, Nemhauser, Savelsbergh cuts generated.
    UInt singlectwo; /// Number of general lifted only a single element for C2
    /// is downlifted.
    UInt basic; /// Number of cuts generated from initial cover set.
    // Time is not calculated yet.
    UInt noviol; // GNS procedure terminated since there is no violation after
    // upliftin F
    UInt noinitcov; // GNS procedure terminated since there is no initial cover generated.
    double time; // Total time used by generator.
  };



  typedef boost::shared_ptr<const Cut> ConstCutPtr;
  typedef CovCutGenStats* CovCutGenStatsPtr; 
  typedef CovCutGenStats const * ConstCovCutGenStatsPtr;
  typedef boost::shared_ptr<ofstream> OfstreamPtr;

  /**
   * The CoverCutGenerator class generates a set of minimal covers 
   * for each knapsack constraint.
   * First, the knapsack list is created.
   * Second, Cover cuts for a chosen minimal cover for a corresponding
   * knapsack inequality is created
   * the inequalites are added to a cut list.
   */
  class CoverCutGenerator {
  public:
    // Default constructor
    CoverCutGenerator();

    // Constructor that uses a problem and a solution given.
    CoverCutGenerator(ProblemPtr p, SolutionPtr s, EnvPtr env);

    // Constructor that uses a relaxation and a solution given.
    CoverCutGenerator(RelaxationPtr rel, ConstSolutionPtr sol, EnvPtr env);

    // Destructor
    ~CoverCutGenerator();

    // Initialize data elements.
    bool initialize(RelaxationPtr p, ConstSolutionPtr s, EnvPtr env);

    // Obtain the knapsack constraint list from problem.
    void generateKnapList();

    // Checks if the constraint has a cover set.
    bool hasCover(ConstraintIterator it);

    // Check if it is a GUB.
    bool GUB(ConstraintIterator itcons);

    // Constructs a vector of nonzero variables in the given solution.
    void nonzeroVars(LinearFunctionPtr lf,
                     CoverSetPtr nonzerovars,
                     CoverSetPtr zerovars);

    // Sort nonzero variables array in nonincreasing order.
    void sortNonIncreasing(CoverSetPtr nonzeros);

    // Sort the variables in nondecreasing order of their reduced costs.
    void sortReducedCosts(CoverSetPtr & vars);

    /** Generates a cover set from a vector of variables and their coefficients
     * in the knapsack constraint.
     * It uses Gu, Nemhauser, Savelsbergh approach.
     */
    bool coverSetGeneratorGNS(ConstConstraintPtr cons, CoverSetPtr cover);



    /* Modified GNS cover set generator such that it always generates a cover.
       This function generates a cover set by using different strategies.
     * For now, it generates a modified version of Gu, Nemhauser, Savelsbergh 
     * such that it generates a cover set even though the number of nonzero
     * elements in solution vector given is not enough for cover set generation 
     * s.t sum(a_i) <= b for i s.t. x^*_i != 0.
     */
    CoverSetPtr coverSetGenGNSModified(ConstConstraintPtr cons);


    /** This is the default cover set generator.
     * It simply adds the variables to the set until the sum > b.
     * It orders the variables in order of their coefficients in nonincreasing
     * order. 
     * By this way, the cover will probably be minimal cover.
     */
    CoverSetPtr coverSetGeneratorDefault(ConstConstraintPtr cons);


    // Obtain initial cover.
    //CoverSetPtr initialCover(ConstConstraintPtr constraintPtr);

    /**
     * This generates the variable-coefficient pair vector 
     * from a given linear function.
     */
    CoverSetPtr varCoeff(LinearFunctionPtr lf);

    /** Constructs the variable-value pair vector from a given vector and a
     * linear function that inclued coefficients.
     */
    void variableCoeffPair(CoverSetPtr cover,LinearFunctionPtr lf);

    /** Calculates the sum of coefficients of a given vector of
     * variable-coefficient pairs
     */
    double sumCoeffs(CoverSetPtr cover);

    /** Drops some of the variables and obtains a minimal cover form a given
     * cover.
     */
    void minimalCover(CoverSetPtr cover,
                      ConstConstraintPtr cons);

    // Generates all the cover cuts from all knapsack constraints.
    void generateAllCuts();

    // Generate cover partitions C1, C2 and Cbar according to Gu, Nemhauser,
    // Savelsbergh.
    void coverPartitionGNS(const ConstConstraintPtr cons,
                           const ConstCoverSetPtr cover,
                           CoverSetPtr cone,
                           CoverSetPtr ctwo,
                           CoverSetPtr fset,
                           CoverSetPtr rset);

    // Generates set N\C, the variables outside of cover set.
    void cBar(const ConstCoverSetPtr coverSetPtr, 
              CoverSetPtr cBar,
              const ConstConstraintPtr constraint);

    // Lifts a variables up or down as it is specified by uplift.
    double lift(const ConstCoverSetPtr obj,
                const ConstCoverSetPtr constraint,
                const CoverSetConstIterator variable,
                double & rhs,
                double & inititalb,
                bool uplift);

    // Simple lifted cover
    void simple(const ConstCoverSetPtr cover,
                const ConstCoverSetPtr cbar,
                const ConstConstraintPtr cons);

    // Simple lifted cover
    void allCTwo(const ConstCoverSetPtr cover,
                 const ConstCoverSetPtr cone,
                 const ConstCoverSetPtr cbar,
                 const ConstConstraintPtr cons);  

    // Initialize the cover inequality by changing the coefficients of cover set
    // by 1s.
    void initCoverIneq(const ConstCoverSetPtr coverset, CoverSetPtr coverineq);

    // Lifting strategy of Gu, Nemhauser, and Savelsbergh.i
    bool liftingGNS(const ConstCoverSetPtr cone,
                    const ConstCoverSetPtr ctwo,
                    const ConstCoverSetPtr fset,
                    const ConstCoverSetPtr rset,
                    CoverSetPtr constraint,
                    const ConstConstraintPtr cons,
                    double & ub);


    // This lifts the variables in a given set according to Gu, Nemhauser and
    // Savelsbergh algorithm by considering the amount of contribution.
    void liftSetF(CoverSetPtr obj,
                  CoverSetPtr consknap,
                  const ConstCoverSetPtr setf,
                  CoverSetPtr coverineq,
                  double & ub,
                  double & initialb,
                  const bool liftup);

    // This function lifts up and lifts down the variables similar to liftSet
    // but the assumprions are changed as the same as liftSetGNS.
    void liftSet(CoverSetPtr obj,
                 CoverSetPtr consknap,
                 const ConstCoverSetPtr varset,
                 CoverSetPtr constraint,
                 double & ub,
                 double & initialb,
                 bool liftup);

    // Generates a Gu, Nemhauser, Savelsbergh lifted cover inequality.
    bool GNS(const ConstConstraintPtr cons);

    // Generates a cover cut from a cover set.
    CutPtr generateCut(const ConstCoverSetPtr constraint, const double ub);

    // Add the cut to the list and insert the corresponding violation to
    // violation list.
    void addCut(CutPtr cut);

    // Check if the same cut is already included.
    bool checkExists(CoverSetPtr cov, double rhs);

    // Generates all the cover cuts from a given knapsack constraint.
    void generateCuts(ConstConstraintPtr constraint);

    // Generates an extended cover from a given cover set.
    void extendedCover(CoverSetPtr cover, ConstConstraintPtr cons); 

    // Implementation based on Horowitz-Shahni algorithm implementation in
    // CglKnapsackCover
    UInt binaryKnapsackSolver(UInt n, double b, double const * c,
                              double const *a, double & z, int * x);

    // Calculates the violation for the given cut.
    double violation(CutPtr cut);

    // Return const pointer for problem, check if this works!!!.
    ConstProblemPtr getProblem() const {return ConstProblemPtr(p_);}    

    // Return const knapsack list.
    ConstKnapsackListPtr getKnapsackList() const 
    {return ConstKnapsackListPtr(knapsackListPtr_);}

    // Return const solution.
    ConstSolutionPtr getSolution() const
    {return ConstSolutionPtr(s_);}

    // Return const cut list.
    CutVector getCutList() const
    {return cutVec_;}

    // Return violation list.
    DoubleVector getViolList() const
    {return violList_;}

    // Return number of constraints considered.
    UInt getNumCons() const {return numCons_;}

    // Return statistics of cover cut generator.
    ConstCovCutGenStatsPtr getStats() const {return ConstCovCutGenStatsPtr(stats_);}

    // Initialize statistics
    void initStats();

    // Adds a cut from a given cover set and rhs by checking integrality and if
    // the cut already exists.
    bool addCut(CoverSetPtr cov, double rhs, UInt cuttype);

    // Return only the violated cuts.
    CutVector getViolatedCutList() const
    {return violatedCuts_;};

    // Check if the given solution satisfies integrality for given problem.
    bool checkIntegral(RelaxationPtr p, ConstSolutionPtr s);

    // Print inequality
    void printIneq(const ConstCoverSetPtr cov, double rhs, 
                   PrintType type, string message);

    // Print lifting problem.
    void printLiftProb(const ConstCoverSetPtr obj,
                       const ConstCoverSetPtr consknap,
                       const CoverSetConstIterator variable,
                       double rhs,
                       double initialb,
                       bool uplift,
                       double b,
                       double gamma,
                       double alpha);

  private:
    // Environment.
    EnvPtr env_;

    // Problem that cover cuts will be generated for.
    ProblemPtr p_;

    // List of cuts generated.
    CutVector cutVec_;

    // List of violated cuts.
    CutVector violatedCuts_;

    // List of violations obtained from cuts generated in the same order.
    DoubleVector violList_;

    // List of knapsack inequalities in the problem.
    KnapsackListPtr knapsackListPtr_;

    /**
     * Given (possibly fractional) solution. 
     * Cut will be designed to violate this solution.
     */
    ConstSolutionPtr s_;

    // Number of knapsack inequalities considered.
    UInt numCons_;

    // Statistics for cover cut generator.
    CovCutGenStatsPtr stats_;

    // Hash map that is used to check if a cut is already created or not.
    std::map< std::vector<double>, UInt> cutmap;  

    // Integer tolerance.
    double intTol_;

    // Objective tolerance.
    double objtol_;

    // Output file.
    OfstreamPtr  output_;

    // Output file name
    string  outfile_;
  };
}

#endif // MINOTAURCOVERCUTGENERATOR_H


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
