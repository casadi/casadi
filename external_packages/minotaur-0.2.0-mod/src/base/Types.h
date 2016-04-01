// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Types.h
 * \brief Declare important 'types' used in Minotaur.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURTYPES_H
#define MINOTAURTYPES_H

#include <boost/shared_ptr.hpp>
#include <list>
#include <map>
#include <deque>
#include <set>
#include <stack>
#include <vector>

//#ifndef MINOTAUR_RUSAGE
//# error MINOTAUR_RUSAGE is NOT defined
//#endif
namespace Minotaur {

  /// Standard types.
  typedef unsigned int UInt;

  /// Containers for standard types
  typedef std::deque  < UInt >   UIntQ;
  typedef std::vector < bool >   BoolVector;
  typedef std::vector < double > DoubleVector;
  typedef std::vector < int >    IntVector;
  typedef std::vector < UInt >   UIntVector;
  typedef std::set    < UInt >   UIntSet;

  /// The different classes of problems
  typedef enum {
    LP,
    MILP,
    QP,
    MIQP,
    QCQP,
    MIQCQP,
    POLYP,
    MIPOLYP,
    NLP,
    MINLP,
    UnknownProblem
  } ProblemType;

  /**
   * Objective sense. Minotaur minimizes. Always. If the sense is to Maximize,
   * then Minotaur converts the objective.
   */
  typedef enum {
    Minimize,
    Maximize
  } ObjectiveType;

  /// Different types of functions in Minotaur.
  typedef enum {
    Constant,
    Linear,
    Bilinear,
    Multilinear,
    Quadratic,
    Nonlinear,
    Polynomial,
    UnknownFunction
  } FunctionType;

  /// Different types of variables.
  typedef enum {
    Binary,     ///< Variable is constrained to be binary.
    Integer,    ///< Variable is constrained to be integer.
    ImplBin,    ///< Variable is continuous, but will take binary values only.
    ImplInt,    ///< Variable is continuous, but will take integer values only.
    Continuous  ///< Variable is continuous.
  } VariableType;

  /// Where did the variable originate from?
  typedef enum {
    VarOrig,    ///> Variable was part of the original problem
    VarTran,    ///> Variable created by Transformer
    VarHand,    ///> Variable created by Handler
    VarOtherSrc ///> Variable created by some other source
  } VarSrcType;
  
  /// Different types of variable-bounds.
  typedef enum {
    Lower,
    Upper
  } BoundType;

  /// SOS types
  typedef enum {
    SOS1,
    SOS2
  } SOSType;

  /// Different states a variable can be in.
  typedef enum {
    DeletedVar, ///< Marked for deletion.
    FixedVar,   ///< Fixed to a value.
    FreeVar,    ///< Doesn't have bounds. Same as superbasic in LP
    NormalVar   ///< Not in any other category
  } VarState;

  /// Different states a constraint can be in.
  typedef enum {
    DeletedCons,  ///< Marked for deletion.
    FreeCons,     ///< Doesn't have bounds. Implies redundant.
    NormalCons    ///< Not in any other category
  } ConsState;

  /// Different states an objective can be in.
  typedef enum {
    DeletedObj,  ///< Marked for deletion.
    NormalObj    ///< Not in any other category
  } ObjState;

  /// Different states an algorithm like branch-and-bound can be in.
  typedef enum {
    NotStarted,
    Started,
    Restarted,
    SolvedOptimal,
    SolvedInfeasible,
    SolvedUnbounded,
    SolvedGapLimit,
    SolvedSolsLimit,
    IterationLimitReached,
    Interrupted,
    TimeLimitReached,
    SolLimitReached,
    SolveError,
    Finished
  } SolveStatus;

  /// Different status that an external engine may report.
  typedef enum {
    ProvenOptimal,
    ProvenLocalOptimal,
    ProvenInfeasible,
    ProvenLocalInfeasible, 
    ProvenUnbounded,
    ProvenObjectiveCutOff,
    EngineIterationLimit,
    ProvenFailedCQFeas,
    ProvenFailedCQInfeas,
    FailedFeas,
    FailedInfeas,
    EngineError,
    EngineUnknownStatus
  } EngineStatus;     

  /// What can a brancher do to a node in branch-and-bound.
  typedef enum {
    ModifiedByBrancher,
    PrunedByBrancher,
    NotModifiedByBrancher
  } BrancherStatus;     

  /// Two directions for branching
  typedef enum {
    DownBranch,
    UpBranch
  } BranchDirection;

  /// Status of a node in branch-and-bound.
  typedef enum {
    NodeNotProcessed, //< has been created but not yet been processed.
    NodeInfeasible,   //< can be pruned
    NodeHitUb,        //< can be pruned
    NodeDominated,    //< can be pruned
    NodeOptimal,      //< The solution obtained is integral, feasible and optimal 
                      //< for the subtree. can prune.
    NodeContinue,     //< solution is neither optimal, nor the relaxation 
                      //< provably infeasible.
    //NodeFrac,       //< fractional solution but optimal for the relaxation
    //NodeInt,        //< integer solution found at this node,
                      //< but may not be feasible/optimal. branch more.
    NodeStopped 
  } NodeStatus;     


  /// Status from separation routine:
  typedef enum {
    SepaContinue, ///< Separation routine wants nothing.
    SepaResolve,  ///< Separation routine wants the relaxation resolved.
    SepaPrune,    ///< No need to further solve the subproblem.
    SepaNone,     ///< No separation found.
    SepaError     ///< Problem separating a point.
  } SeparationStatus;


  /// Levels of verbosity.
  typedef enum {
    LogNone,          //< absolutely no messages
    LogError,         //< only error messages and warnings if any
    LogInfo,          //< errors, timed status, final solution
    LogExtraInfo,     //< some extra info, stats etc
    LogDebug,         //< more verbose
    LogDebug1,        //< more verbose
    LogDebug2         //< more verbose
  } LogLevel;

  /// Order of tree search.
  typedef enum {
    DepthFirst,   
    BestFirst,
    BestThenDive     /// First find the best bound, then dive until pruned.
  } TreeSearchOrder;


  /// Type of algorithms that can be used
  typedef enum {
    DefaultAlgo,
    QPNLPBnb,
    QG,
    NLPBnb
  } AlgoType;

  /// Logging
  class Logger;
  typedef boost::shared_ptr<Logger> LoggerPtr;

  /// What is the function type on adding two functions f1 and f2?
  FunctionType funcTypesAdd(FunctionType f1, FunctionType f2);

  /// What is the function type on multiplying two functions f1 and f2?
  FunctionType funcTypesMult(FunctionType f1, FunctionType f2);

  /// Get the string equivalent of ProblemType.
  std::string getProblemTypeString(ProblemType p);

  /// Get the string equivalent of FunctionType.
  std::string getFunctionTypeString(FunctionType f);

  /// Get the string equivalent of SolveStatus.
  std::string getSolveStatusString(SolveStatus s);

  class Constraint;
  typedef boost::shared_ptr<Constraint> ConstraintPtr;
  typedef boost::shared_ptr<const Constraint> ConstConstraintPtr;
  typedef std::vector<ConstraintPtr>::iterator ConstraintIterator;
  typedef std::vector<ConstraintPtr>::const_iterator ConstraintConstIterator;
  typedef std::vector<ConstraintPtr> ConstraintVector;

  class SOS;
  typedef std::vector<SOS *> SOSVector;
  typedef SOSVector::iterator SOSIterator;
  typedef SOSVector::const_iterator SOSConstIterator;

  // Serdar defined
  typedef boost::shared_ptr<ConstraintVector> ConstraintVectorPtr;
  typedef std::vector<ConstConstraintPtr> ConstConstraintVector;
  typedef boost::shared_ptr<ConstConstraintVector> ConstConstraintVectorPtr;
  typedef ConstConstraintVector::iterator ConstConstraintIterator;
  typedef std::set<ConstraintPtr> ConstrSet;
  typedef std::deque<ConstConstraintPtr> ConstrQ;


  class Variable;
  typedef boost::shared_ptr<Variable> VariablePtr;
  typedef boost::shared_ptr<Variable> ConstVariablePtr;
  typedef std::vector<VariablePtr> VarVector;
  typedef std::deque<VariablePtr> VarQueue;
  typedef VarVector::const_iterator VariableConstIterator;
  typedef VarVector::iterator VariableIterator;
  typedef std::vector<ConstVariablePtr>::iterator ConstVarIter;
  typedef VarQueue::const_iterator VarQueueConstIter;
  typedef VarQueue::iterator VarQueueIter;
  struct CompareVariablePtr {
    bool operator()(ConstVariablePtr v1, ConstVariablePtr v2) const;
  };
  typedef std::set<ConstVariablePtr, CompareVariablePtr> VariableSet;
  typedef VariableSet::const_iterator VarSetConstIterator;
  typedef std::set<VariablePtr> VarSet;

  // Serdar added.
  typedef boost::shared_ptr<VarSet> VarSetPtr;
  typedef boost::shared_ptr<const VarSet> ConstVarSetPtr;
  typedef std::map<ConstVariablePtr, UInt, CompareVariablePtr> VarIntMap;
  typedef VarIntMap::const_iterator VarIntMapConstIterator;
  typedef VarIntMap::iterator VarIntMapIterator;

  typedef VarSet::iterator           VarSetIter;
  typedef VarSet::const_iterator     VarSetConstIter;
  typedef std::set<ConstVariablePtr> ConstVarSet;
  typedef VarSet::const_iterator     ConstVarSetIter;

  class Node;
  typedef boost::shared_ptr<Node> NodePtr;
  typedef std::vector<NodePtr> NodePtrVector;
  typedef std::vector<NodePtr>::iterator NodePtrIterator;

  class Handler;
  typedef boost::shared_ptr<Handler> HandlerPtr;
  typedef boost::shared_ptr<const Handler> ConstHandlerPtr;
  typedef std::vector<HandlerPtr> HandlerVector;
  typedef HandlerVector::iterator HandlerIterator;
  typedef HandlerVector::const_iterator HandlerConstIterator;

  class   Environment;
  typedef boost::shared_ptr<Environment> EnvPtr;
  class   Problem;
  typedef boost::shared_ptr<Problem> ProblemPtr;
  typedef boost::shared_ptr<const Problem> ConstProblemPtr;

  class   Modification;
  typedef boost::shared_ptr<Modification> ModificationPtr;
  typedef boost::shared_ptr<const Modification> ConstModificationPtr;  
  typedef std::vector<ModificationPtr> ModVector;  
  typedef std::vector<ModificationPtr>::const_iterator 
          ModificationConstIterator;
  typedef std::vector<ModificationPtr>::const_reverse_iterator 
          ModificationRConstIterator;
  typedef std::deque<ModificationPtr> ModQ;
  typedef std::stack<ModificationPtr> ModStack;

  class   Branch;
  typedef boost::shared_ptr<Branch> BranchPtr;
  typedef boost::shared_ptr<const Branch> ConstBranchPtr;  
  typedef std::vector<BranchPtr> BranchPtrVector;
  typedef boost::shared_ptr<BranchPtrVector> Branches;
  typedef std::vector<BranchPtr>::const_iterator BranchConstIterator;

  class   BrCand;
  typedef boost::shared_ptr<BrCand> BrCandPtr;
  typedef boost::shared_ptr<const BrCand> ConstBrCandPtr;  
  typedef std::vector< BrCandPtr > BrCandVector;
  typedef BrCandVector::iterator BrCandVIter;

  class   BrVarCand;
  typedef boost::shared_ptr<BrVarCand> BrVarCandPtr;
  typedef boost::shared_ptr<const BrVarCand> ConstBrVarCandPtr; 
  struct CompareVarBrCand {
    bool operator()(ConstBrVarCandPtr c1, ConstBrVarCandPtr c2) const;
  };
  typedef std::set<BrVarCandPtr, CompareVarBrCand> BrVarCandSet;
  typedef BrVarCandSet::iterator BrVarCandIter;
  typedef std::vector< BrVarCandPtr >::iterator BrVarCandVIter;

  // Serdar added this block 
  // Compare some pairs.
  typedef std::pair<int,double> id;
  struct CompareIntDouble {
    bool operator() (id id1, id id2) const;
  };

  // Compare rule to compare the value of the pair.
  typedef std::pair<ConstVariablePtr, double> VariableValuePair;
  typedef boost::shared_ptr<VariableValuePair> VariableValuePairPtr;
  struct CompareValueVariablePair {
    bool operator() (VariableValuePair v1,VariableValuePair v2) const;
  };

  // Map to determine GUB cover.
  struct CompareValueVarInc{
    bool operator() (VariableValuePair v1, VariableValuePair v2) const;
  };
  typedef std::map<ConstVariablePtr, double, CompareValueVarInc>
          VariableValueMap;


  // Vector of pair<ConstVariablePtr,double value> which is VariableValuePair.
  typedef std::vector<VariableValuePair> VariableValuePairVector;
  typedef boost::shared_ptr<VariableValuePairVector> VariableValuePairVectorPtr;
  typedef std::vector<VariableValuePair>::iterator
          VariableValuePairVectorIterator;
  typedef std::vector<VariableValuePair>::const_iterator 
          VariableValuePairVectorConstIterator;

  // Cover is a VariableSet
  typedef VariableValuePairVector CoverSet;
  typedef boost::shared_ptr<CoverSet> CoverSetPtr;
  typedef boost::shared_ptr<const CoverSet> ConstCoverSetPtr;
  typedef CoverSet::iterator CoverSetIterator;
  typedef CoverSet::const_iterator CoverSetConstIterator; 

  // Pointers to classes added.
  class KnapsackList;
  typedef boost::shared_ptr<KnapsackList> KnapsackListPtr;
  typedef boost::shared_ptr<const KnapsackList> ConstKnapsackListPtr;

  // Cut class vectors etc.
  class Cut;
  typedef boost::shared_ptr<Cut> CutPtr;
  typedef std::vector<CutPtr> CutVector;
  typedef CutVector::iterator CutVectorIter;
  typedef CutVector::const_iterator CutVectorConstIter;
  typedef std::list<CutPtr> CutList;
  typedef CutList::iterator CutListIter;

  //LiftingProblem which is a knapsack problem type of Problem.
  typedef Problem LiftingProblem;
  typedef ProblemPtr LiftingProblemPtr;
  // CoverCutGenerator pointer
  class CoverCutGenerator;
  typedef boost::shared_ptr<CoverCutGenerator> CoverCutGeneratorPtr;
  // LGCIGenerator pointer
  class LGCIGenerator;
  typedef boost::shared_ptr<LGCIGenerator> LGCIGeneratorPtr;

  class Function;
  typedef boost::shared_ptr<const Function> ConstFunctionPtr; 
  // Serdar ended.

  // declare options database
  class OptionDB;
  typedef boost::shared_ptr <OptionDB> OptionDBPtr;

  // declare options
  template <class T> class Option;
  typedef boost::shared_ptr< Option<bool>        > BoolOptionPtr;
  typedef boost::shared_ptr< Option<int>         > IntOptionPtr;
  typedef boost::shared_ptr< Option<double>      > DoubleOptionPtr;
  typedef boost::shared_ptr< Option<std::string> > StringOptionPtr;
  typedef BoolOptionPtr FlagOptionPtr;

  // define sets of pointers
  typedef std::set<BoolOptionPtr> BoolOptionSet;
  typedef std::set<IntOptionPtr> IntOptionSet;
  typedef std::set<DoubleOptionPtr> DoubleOptionSet;
  typedef std::set<StringOptionPtr> StringOptionSet;
  typedef BoolOptionSet FlagOptionSet;

  // define iterators of pointers
  typedef BoolOptionSet::iterator BoolOptionSetIter;
  typedef IntOptionSet::iterator IntOptionSetIter;
  typedef DoubleOptionSet::iterator DoubleOptionSetIter;
  typedef StringOptionSet::iterator StringOptionSetIter;
  typedef BoolOptionSetIter FlagOptionSetIter;


  /// Variables should always be constant within a group
  typedef std::map<ConstVariablePtr, double, CompareVariablePtr> VariableGroup;
  typedef std::map<ConstVariablePtr, double>::iterator VariableGroupIterator;
  typedef std::map<ConstVariablePtr, double>::const_iterator 
    VariableGroupConstIterator;

  /// Pairs of variables are used in quadratic functions.
  typedef std::pair<ConstVariablePtr, ConstVariablePtr> VariablePair;
  typedef std::pair<ConstVariablePtr, ConstVariablePtr> ConstVariablePair;
  
  struct CompareVariablePair {
    bool operator()(ConstVariablePair tv1, ConstVariablePair tv2) const;
  };

  /// Again, Variables should always be constant within a group
  typedef std::map<ConstVariablePair, double, CompareVariablePair> 
    VariablePairGroup;
  typedef std::map<ConstVariablePair, UInt, CompareVariablePair> VarPairIntMap;
  typedef std::map<VariablePair, double>::iterator VariablePairGroupIterator;
  typedef std::map<ConstVariablePair, double>::const_iterator 
    VariablePairGroupConstIterator;
  typedef std::map<ConstVariablePtr, ConstVariablePtr, 
          CompareVariablePtr> VarVarMap;
  typedef VarVarMap::const_iterator VarVarMapConstIterator;

  typedef const std::map<ConstVariablePtr, UInt, CompareVariablePtr> 
    VarCountConstMap;

  class   Heuristic;
  typedef boost::shared_ptr<Heuristic> HeurPtr;
  typedef std::vector<HeurPtr> HeurVector;
}

#endif

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
