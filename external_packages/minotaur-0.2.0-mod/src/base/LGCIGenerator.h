//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file LGCIGenerator.h 
 * \brief Declare base class KnapsackList. 
 * \author Serdar Yildiz, Argonne National Laboratory 
 */

#ifndef MINOTAURLGCIGENERATOR_H
#define MINOTAURLGCIGENERATOR_H

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
#include "ProbStructure.h"
#include "LPEngine.h"

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
  Gns,
  Noviol,
  Noinitcover
} LGCIType;
  
// Shows why cut is not added to the list.
typedef enum {
  Duplicate = 0,
  NotViolated
}CutFail;

struct LGCIGenStats
{
  UInt liftsubprobs; /// Number of total subproblems solved for lifting.
  UInt totalcuts; /// Number of all cuts i.e. including duplicates as well.
  UInt cuts; /// Number of total cover cuts generated excluding duplicates.
  UInt violated; /// Number of violated cuts.
  UInt knapcons; // Number of knapsack constraints.
  UInt knapcov; /// Number of knapsack inequalities considered that has cover
                /// inequalities.
  UInt knapgub; /// Number of knapsack inequalities that has a set of GUBs
                /// that cover its variables.
  UInt gns;     /// Number of Gu, Nemhauser, Savelsbergh cuts generated.
  UInt noviol;  /// Number cuts that does not violate current relaxation
                /// solution.
  UInt noinitcov; /// GNS procedure terminated since no inital GNS cover constructed.
  double time; /// Total time used by generator.
};

// Typedefs for ease of coding.
class LGCIGenerator;
typedef boost::shared_ptr<LGCIGenerator> LGCIGeneratorPtr;
typedef boost::shared_ptr<const LGCIGenerator> ConstLGCIGeneratorPtr;
typedef LGCIGenStats* LGCIGenStatsPtr;
typedef LGCIGenStats const * ConstLGCIGenStatsPtr;
typedef boost::shared_ptr<const LinearFunction> ConstLinearFunctionPtr;
typedef std::map<ConstVariablePtr, ConstVariablePtr> OrigLiftVars;
typedef boost::shared_ptr<OrigLiftVars> OrigLiftVarsPtr;
  
/**
 * LGCIGenerator class generates lifted GUB cover inequalities from a 
 * given problem and a solution.
 * It generates the list of knapsack inequalities suitable for cut generation.
 * It checks if the considered knapsack constraint has at least one cover.
 * If not then it checks another knapsack inequality.
 * It generates an overlapping GUB set that covers knapsack variables.
 * It by using a simple elimination rule, i.e., take out duplicates, obtains a non-overlapping GUB set.
 */
class LGCIGenerator {
public:
  // Default constructor.
  LGCIGenerator();
    
  // Constructor that uses a problem and a solution.
  LGCIGenerator(ProblemPtr p, SolutionPtr s, 
                EnvPtr env, LPEnginePtr lpengine);

  // Constructor that uses a relaxation and a solution given.
  LGCIGenerator(RelaxationPtr rel, ConstSolutionPtr sol, 
                EnvPtr env, LPEnginePtr lpengine);
    
  // Destructor.
  ~LGCIGenerator();

  // Initialize data elements.
  void initialize();

  // Initialize statistics.
  void initStats();

  // Checks if the constraint has a cover set.
  bool hasCover(ConstraintConstIterator it);

  // Checks if there is a set of GUBs that covers the knapsack inequality.
  bool hasGubCover(ConstraintConstIterator it,
                   ConstConstraintVectorPtr gublist);

  // Generates initial cover by using GNS algorithm.
  bool coverSetGeneratorGNS(ConstConstraintPtr cons, 
                            CoverSetPtr cover);

  /* Modified GNS cover set generator such that it always generates a cover.
   * This function generates a cover set by using different strategies.
   * For now, it generates a modified version of Gu, Nemhauser, Savelsbergh 
   * such that it generates a cover set even though the number of nonzero
   * elements in solution vector given is not enough for cover set generation 
   * s.t sum(a_i) <= b for i s.t. x^*_i != 0.
   */
  bool coverSetGenGNSModified(ConstConstraintPtr cons);

  /** Constructs the variable-value pair vector from a given vector and a
   * linear function that inclued coefficients.
   */
  void varCoeff(ConstConstraintPtr cons, CoverSetPtr cover);

  /** Calculates the sum of coefficients of a given vector of
   * variable-coefficient pairs
   */
  double sumCoeffs(CoverSetPtr cover);
    
  /** Drops some of the variables and obtains a minimal cover form a given
   * cover.
   */
  void minimalCover(CoverSetPtr cover, ConstConstraintPtr cons);

  // Generates all LGCI cover cuts from all constraints.
  void generateAllCuts();
    
  // Generate cover partitions C1, C2, F and R
  void coverPartitionGNS(const ConstConstraintPtr cons,
                         const ConstCoverSetPtr cover,
                         CoverSetPtr cone,
                         CoverSetPtr ctwo,
                         CoverSetPtr fset,
                         CoverSetPtr rset);

  // Lifts a variable up or down as it is specified by uplift.
  double lift(const ConstCoverSetPtr obj,
              const ConstCoverSetPtr inequality,
              const CoverSetConstIterator variable,
              double & rhs,
              double & initialb,
              bool uplift);

  // Initialize cover inequality by changing the coefficients of cover set
  // by ones.
  void initCoverIneq(const ConstCoverSetPtr cover, CoverSetPtr coverineq);

  /** Lifting strategy of GNS for LGCI lifting. Different from Knapsack
   * cover genereator.
   * Generates Gu, Nemhauser, Savelsbergh LGCI generation algorithm.
   */
  bool liftingGNS(const ConstCoverSetPtr cone,
                  const ConstCoverSetPtr ctwo,
                  const ConstCoverSetPtr fset,
                  const ConstCoverSetPtr rset,
                  CoverSetPtr constraint,
                  const ConstConstraintPtr cons,
                  ConstConstraintVectorPtr gublist,
                  double & ub);


  /** This lifts the variables in a given set in the given order
   */
  bool liftSet(CoverSetPtr obj,
               CoverSetPtr constraintset,
               const ConstCoverSetPtr varset,
               CoverSetPtr inequality,
               double & ub,
               double & initialb,
               bool liftup);

  // Generates all the cover cuts from a given knapsack constraint.
  void generateCuts(ConstConstraintPtr cons, ConstConstraintVectorPtr gublist);

  // Generate a cut from a given variable set and rhs.
  void generateCut(const ConstCoverSetPtr inequality, double ub, CutPtr cut);

  // Generates a GNS LGCI cut.
  bool GNS(const ConstConstraintPtr cons, ConstConstraintVectorPtr gublist);

  // Add the cut to the list and insert corresponding violation to violation
  // list.
  void addCut(CutPtr cut);
  
  bool addCut(CoverSetPtr cov, double rhs, UInt cuttype, CutFail& failtype);

  // Check if the same cut is already included in the list.
  bool checkExists(CoverSetPtr inequality, double rhs);

  // Solves the lifting problem.
  double liftingProbsolver();

  // Calculates the violation for the given cut.
  double violation(CutPtr cut);

  // Constructs a vector of nonzero variables in the given solution.
  void nonzeroVars(const ConstLinearFunctionPtr lf,
                   CoverSetPtr nonzerovars,
                   CoverSetPtr zerovars);

  // Sort nonzero variab;e array in nonincreasing order.
  void sortNonIncreasing(CoverSetPtr nonzeros);

  // This function generates set N\C, the variables outside of cover set.
  void cBar(const ConstCoverSetPtr cover,
            CoverSetPtr cbar,
            const ConstConstraintPtr cons);

  // Generate GUBs from constraints in map data type.
  void generateGubMaps(ConstConstraintVectorPtr gublist, 
                       boost::shared_ptr< std::vector<CoverSetPtr> > gubcoverlist,
                       boost::shared_ptr< std::vector<CoverSetPtr> > guborigcoeffs);
  
  // Generate non overlapping GUB constraint set.
  void nonOverlap(boost::shared_ptr< std::vector<CoverSetPtr> > gubcoverlist);

  // Initialize the GUB constraints in lifting problem.
  void initGubCons(const ConstCoverSetPtr cone,
                   boost::shared_ptr<std::vector<CoverSetPtr> > gubcoverlist,
                   boost::shared_ptr<std::vector<CoverSetPtr> > gubcons);

  // Lift elements in the set.
  void liftSet(CoverSetPtr obj, 
               boost::shared_ptr<std::vector<CoverSetPtr> >origgubs,
               CoverSetPtr consknap,
               boost::shared_ptr<std::vector<CoverSetPtr> > gubcons,
               const ConstCoverSetPtr varset,
               CoverSetPtr coverineq,
               double & rhs,
               double & initialknap,
               double * initialgub,
               bool liftup);

  // Lift the current variable.
  double lift(CoverSetPtr obj,
              boost::shared_ptr<std::vector<CoverSetPtr> > origgubs,
            CoverSetPtr consknap,
            boost::shared_ptr<std::vector<CoverSetPtr> > gubcons,
            const CoverSetConstIterator variable,
            double & rhs,
            double & initialbknap,
            double * initialgub,
            bool liftup);

  // Return const GUB identifier.
  ConstProbStructPtr getProbStruct() const
    {return probstruct_;}

  // Return const knapsack list.
  ConstKnapsackListPtr getKnapsackList() const
    {return ConstKnapsackListPtr(knapList_);}

  // Return statistics of LGCIGenerator.
  ConstLGCIGenStatsPtr getStats() const
    {return ConstLGCIGenStatsPtr(stats_);}

  // Print inequality.
  void printIneq(const ConstCoverSetPtr cov, double rhs, 
                 PrintType type, string message);


  void addCons(CoverSetPtr obj,
               CoverSetPtr consknap,
               double bknap,
               boost::shared_ptr< std::vector<CoverSetPtr> > gubcons,
               double * bgubs,
               OrigLiftVarsPtr varmap,
               ProblemPtr liftprob);

  VariablePtr addVar(VariablePtr var, OrigLiftVarsPtr varmap, ProblemPtr liftprob);
  

  double roundHeur(ProblemPtr prob);

  // Return const pointer for problem.

private:
  // Environment.
  EnvPtr env_;
  // Problem that cover cuts will be generated for.
  ProblemPtr p_;   
  // List of cuts generated.
  CutVector cutVec_;
  // List of violated cut.
  CutVector violatedCuts_;
  // List of violations for cuts corresponding in cut list.
  DoubleVector violList_; 
  // List of knapsack inequalities.
  KnapsackListPtr knapList_; 
  /**
   * Given (possibly fractional) solution. 
   * Cut will be designed to violate this solution.
   */
  SolutionPtr s_;
  // Problem structure to obtaib GUB list for variables in knapsack inequality.
  ConstProbStructPtr probstruct_;  
  // Number of knapsack ineqs used for cut generation.
  UInt numCons_; 
  // Statistics for LGCI generator.
  LGCIGenStatsPtr stats_;
  // Hash map that is used to check if a cut is already created or not.
  std::map< std::vector<double>, UInt> cutmap_;
  // Integer tolerance.
  double intTol_;

  // Output file
  ofstream output_;
  // Output file name
  string outfile_;
  // LP solver.
  LPEnginePtr lpengine_;

};
  
}




#endif // MINOTAURLGCIGENERATOR_H





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
 
