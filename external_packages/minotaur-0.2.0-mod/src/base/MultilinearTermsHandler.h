//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file MultilinearTermsHandler.h
 * \brief A Handler for a general collection of MultilinearTerms
 */

#ifndef MINOTAURMULTILINEARTERMSHANDLER_H
#define MINOTAURMULTILINEARTERMSHANDLER_H

#include <algorithm>

#include "Handler.h"
#include "LPEngine.h"


namespace Minotaur {

//typedef boost::shared_ptr<const Problem> ConstProblemPtr;
//typedef boost::shared_ptr<const Variable> ConstVariablePtr;

// These should maybe go in Terms.h and be public -- they seem pretty generic
typedef std::set<ConstVariablePtr> SetOfVars;
typedef std::set<SetOfVars> SetOfSetOfVars;


/*
 * This class is needed to store information about term cover
 *  It is not a horribly efficient implementation, as I didn't even bother
 *  to make an adjacency list (yet)
 */

class Hypergraph {

public:

class CompareSetsOfVars
{
public:
  bool operator()(const SetOfVars &lhs, const SetOfVars &rhs ) const
    {
      return lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
    }
};
  
  typedef std::map<SetOfVars, double, CompareSetsOfVars> SetOfVarsDoubleMap;
  typedef std::list<SetOfVars> ListOfSetOfVars;
  typedef std::map<ConstVariablePtr, ListOfSetOfVars> AdjListType;
  
public:
  Hypergraph(ConstProblemPtr p) : problem_(p) {} ;
  virtual ~Hypergraph () {};

  void adjustEdgeWeightsBetween(const VariablePtr v, const SetOfVars &g, 
                                bool phaseOne);
  void create(std::map<ConstVariablePtr, SetOfVars > const &terms);
  double getWeight(const SetOfVars &e);
  SetOfVars heaviestEdge(bool &maxWeightPositive) const;
  VariablePtr heaviestIncidentVertex(const SetOfVars &g);
  ListOfSetOfVars incidentEdges(ConstVariablePtr v) const;

  VariablePtr maxWeightedDegreeVertex(bool &maxWeightPositive) const;
  int numEdges() const { return E_.size(); }
  int numVertices() const { return V_.size(); }

  SetOfVars randomEdge(bool &maxWeightPositive);
  void resetWeights();
  void setWeight(const SetOfVars &e, double w);
  double weightedDegree(ConstVariablePtr v) const;
  void write(std::ostream &out) const;

private:

  ConstProblemPtr problem_;

  SetOfVars V_;
  AdjListType adjList_;

  SetOfSetOfVars E_;
  SetOfVarsDoubleMap weights_;
  SetOfVarsDoubleMap originalWeights_;
  
};


typedef boost::shared_ptr<Hypergraph> HypergraphPtr;
typedef boost::shared_ptr<const Hypergraph> HypergraphConstPtr;

/** A MultilinearTermsHandler handles all (single) multilinear term
 * constraints
 */
class MultilinearTermsHandler : public Handler {

public:
  typedef enum { 
    relaxInit_Call, 
    relaxNodeInc_Call, 
    getBrMod_Call 
  } HandleCallingFunction;


public:
    
  /// Default constructor.
  MultilinearTermsHandler(EnvPtr env, ProblemPtr problem); 
        
  /// Destroy.
  virtual ~MultilinearTermsHandler() {};

  /**
   * Add constraint
   */
  void addConstraint(ConstraintPtr newcon, ConstVariablePtr ovar,
                     std::set<ConstVariablePtr> ivars);
  
private:
  void addConstraint(ConstraintPtr) {};

public:

  // base class method
  virtual void getBranchingCandidates(RelaxationPtr rel, 
                                      const DoubleVector &x, ModVector &mods, 
                                      BrVarCandSet &cands, BrCandVector &,
                                      bool &is_inf);
    
  /**
   * Check if each multilinear term is satisfied. Stops on the first
   * violated constraint.
   */
  bool isFeasible(ConstSolutionPtr, RelaxationPtr, bool &should_prune,
                  double &inf_meas );
    
  // Does nothing.
  void relaxInitFull(RelaxationPtr, bool *) { assert(0); } ;

  // Build initial relaxation
  void relaxInitInc(RelaxationPtr, bool *);

  // Does nothing.
  void relaxNodeFull(NodePtr, RelaxationPtr, bool *) {assert(0); };

  // All changes are in the branches to the relaxation, so this need not do anything
  void relaxNodeInc(NodePtr n, RelaxationPtr r , bool *is_inf);

    
  /// Can not return any cuts for this case.
  void separate(ConstSolutionPtr, NodePtr , RelaxationPtr , CutManager *, 
                SolutionPoolPtr, bool *, SeparationStatus *) {};
    
  virtual ModificationPtr getBrMod(BrCandPtr , DoubleVector &, 
                                   RelaxationPtr , BranchDirection );
    
  virtual Branches getBranches(BrCandPtr cand, DoubleVector &x , 
                               RelaxationPtr rel, SolutionPoolPtr s_pool);

  // presolve.
  virtual SolveStatus presolve(PreModQ *, bool *) {return Finished;};
    
  // Implement Handler::presolveNode()
  virtual bool presolveNode(RelaxationPtr, NodePtr, SolutionPoolPtr, ModVector &,
                    ModVector &)
  {return false;};

    
  // Write name
  std::string getName() const { return std::string("Multilinear Term Handler"); }
    
protected:

  /// Environment
  EnvPtr env_;
  
  /// The handler 'handles' constraints of this problem
  ConstProblemPtr problem_;
    
private:
    
  LPEnginePtr lp_engine_;
    
  // (Absolute?) error tolerance for branching.
  double eTol_;
  LoggerPtr logger_;

  // Parameters for term cover
  UInt maxGroupSize_;
  double augmentCoverFactor_;
  int initialTermCoverSize_;

private:
  typedef std::set<ConstVariablePtr> SetOfVars;
  typedef std::set<SetOfVars> SetOfSetOfVars;

  // This contains map in *original* variables (not relaxation variables)
  typedef std::map<ConstVariablePtr, SetOfVars > TermContainer;
  TermContainer termsO_;

  TermContainer termsR_;

  // This contains relaxation variable pointers
  typedef std::vector<SetOfVars> GroupContainer;
  GroupContainer groups_;

  typedef std::vector<SetOfSetOfVars> PointContainer;
  PointContainer points_;  

  // Stores the lambda variables for group <groupIx,pointIx>
  //  These (of course) will be variables in the relaxation
  std::vector<std::vector<ConstVariablePtr> > lambdavars_;

#if 0
  //XXX Don't need these anymore?
  // Multimap between <x_j, Constraint expressing x_j as convex comb>
  // The map is between <relaxation var, relaxation constraint>
  typedef std::multimap<ConstVariablePtr, ConstraintPtr> VariableConstraintMMap;
  VariableConstraintMMap xConMMap_;

  typedef std::map<ConstConstraintPtr, int> ConstraintIntMap;
  ConstraintIntMap conGroupMap_;

  // Map between <z_t, Constraint expressing z_t as convex comb of endpoints>
  // The map is between <relaxation var, relaxation constraint>
  std::map<ConstVariablePtr, ConstraintPtr> zConMap_;
#endif

  // We want (int, VarPtr) as key.  ConstraintPtr as value
  typedef std::pair<int, ConstVariablePtr> IntVarPtrPair;
  typedef std::map<IntVarPtrPair, ConstraintPtr> IntVarPtrPairConstraintMap;
  
  IntVarPtrPairConstraintMap xConMap_;
  IntVarPtrPairConstraintMap zConMap_;


  // Hypergraph for termcover
  HypergraphPtr H_;
         
public:
  typedef TermContainer::const_iterator ConstTermIterator;
  typedef TermContainer::iterator TermIterator;

  typedef GroupContainer::const_iterator ConstGroupIterator;
  typedef GroupContainer::iterator GroupIterator;


private:

  // A bunch of helper functions
  void addEdgeToGroups_(const SetOfVars &e, bool phaseOne);
  bool allVarsBinary_(const SetOfVars &v) const;
  bool edgeIsContainedInGroup_(const SetOfVars &e, const SetOfVars &g) const;
  bool edgeWillFitInGroup_(const SetOfVars &e, const SetOfVars &g) const;

  // Helper function to do branch
  BranchPtr doBranch_(BranchDirection UpOrDown, ConstVariablePtr v, 
                      double bvalue);

  // A greedy dense term covering heuristic
  void greedyDenseHeuristic_();

  /* This chunk of code adds (or modifies)
     x_{V_g} = \sum_{k=1}^{2 |V_g|} \lambda_k^g \chi^{k,g} \forall g \in G
  */
  void handleXDefConstraints_(RelaxationPtr relaxation, HandleCallingFunction wherefrom, ModVector &mods);

  /* This chunk of code adds
   z_t = \sum_{k=1}^{2 |V_g|} \lambda_k^g, \Prod_{j \in J_t} \chi_j^{g,k}
                            \forall J_t \subseteq V_g
  */         
  void handleZDefConstraints_(RelaxationPtr relaxation, HandleCallingFunction wherefrom, ModVector &mods); 
  
  // Helper function to make the groups
  void makeGroups_();

  // Returns the powerset of s
  std::set<SetOfVars> powerset_(SetOfVars const &s);

  // Removes subsets that appear in the groups
  void removeSubsetsFromGroups_();

  // Setup the (initial) weights for weighted term cover approach
  //void setupWeights_();

  bool varsAreGrouped_(SetOfVars const &termvars) const;

  //WeightContainer::iterator findMaxWeight_();

  bool varIsAtLowerBoundAtPoint_(ConstVariablePtr &x, SetOfVars const &p);
  bool varIsAtUpperBoundAtPoint_(ConstVariablePtr &x, SetOfVars const &p) {
    return !(varIsAtLowerBoundAtPoint_(x, p));
  }

  // A stoopid weighted degree heuristic
  void weightedDegreeHeuristic_();

};
typedef boost::shared_ptr<MultilinearTermsHandler> MultilinearTermsHandlerPtr;
typedef boost::shared_ptr<const MultilinearTermsHandler> MultilinearConstHandlerPtr;

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
