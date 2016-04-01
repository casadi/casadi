//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file MultilinearHandler.h
 * \brief Declare the MultilinearHandler class for handling Multilinear 
 * constraints.
 */

#ifndef MINOTAURMULTILINEARHANDLER_H
#define MINOTAURMULTILINEARHANDLER_H

#include "Handler.h"
#include "LPEngine.h"


namespace Minotaur {

class   Function;
class   PolynomialFunction;
class   QuadraticFunction;
class   Problem;
typedef boost::shared_ptr<Function> FunctionPtr;
typedef boost::shared_ptr<PolynomialFunction> PolyFunPtr;
typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
typedef boost::shared_ptr<const Problem> ConstProblemPtr;

/// A MultilinearHandler handles terms like \f$x_0x_1^2\f$.
class MultilinearHandler : public Handler {

public:
    
  /// Default constructor.
  MultilinearHandler(EnvPtr env, ProblemPtr problem); 
        
  /// Destroy.
  virtual ~MultilinearHandler() {};

    
  bool findOriginalVariable(ConstVariablePtr rv, ConstVariablePtr &ov) const;

  /**
   * Both \f$x_0, x_1\f$ are branching candidates. The McCormick
   * inequalities must be updated after branching.
   */
  virtual void getBranchingCandidates(RelaxationPtr rel, 
                                      const DoubleVector &x, ModVector &mods, 
                                      BrVarCandSet &cands,
                                      BrCandVector &gencands, bool &is_inf );
    
  /**
   * Check if each multilinear constraint is satisfied. Stops on the first
   * violated constraint.
   */
  bool isFeasible(ConstSolutionPtr, RelaxationPtr, bool &should_prune,
                  double &inf_meas);
    
  void relax(RelaxationPtr relaxation, bool & isInfeasible);

  // Does nothing.
  void relaxInitFull(RelaxationPtr, bool *) {} ;

  // Does nothing.
  void relaxInitInc(RelaxationPtr, bool *) {};

  // Does nothing.
  void relaxNodeFull(NodePtr, RelaxationPtr, bool *)  {};

  // Does nothing.
  void relaxNodeInc(NodePtr, RelaxationPtr, bool *) {};

  void relaxNode(NodePtr node, RelaxationPtr relaxation, bool & isInfeasible);

    
  /**
   * Return rev_mlterms_
   */    
  std::map <std::vector<ConstVariablePtr>, ConstVariablePtr > getRevMlterms() 
  {return rev_mlterms_; }

  /**
   * Return mlterms_
   */    
  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > getMlterms() 
  {return mlterms_; }

  /**
   * Retrun the max_pow_
   */
  std::map <ConstVariablePtr, UInt> getMaxPow() {return max_pow_; }
    
  /**
   * Retrun the blterms_
   */
  std::map <ConstVariablePtr, ConstVariablePair> getBlterms() {return blterms_; }
    
  /**
   * Retrun the rev_blterms_
   */
  std::map <ConstVariablePair, ConstVariablePtr> getRevBlterms() {return rev_blterms_; }
    
  /**
   * Retrun the monomial_terms__
   */
  std::map <VarIntMap, ConstVariablePtr> getMonomialterms() {return monomial_terms_; }
    
  /**
   * Retrun the sqterms_
   */
  std::map <ConstVariablePtr, ConstVariablePair> getSqterms() {return sqterms_; }    
    
  /**
   * Retrun the rev_sqterms_
   */
  std::map <ConstVariablePair, ConstVariablePtr> getRevSqterms() {return rev_sqterms_; }    
    
  /**
   * Return the groups of variables that are made for the Grouped Convex Hull relaxation  
   */
  std::vector<std::vector<ConstVariablePtr> > getGroups() {return groups_; }
    
  /**
   * Return the lambda variables introduced for the Grouped Convex Hull relaxation
   */
  std::vector<std::vector<ConstVariablePtr> > getAllLambdas() {return all_lambdas_; }
    
  /**
   * Return the map between the original variables and the variables in the relaxation
   */
  std::map <ConstVariablePtr, ConstVariablePtr> getOriginalVariablesMap() {return oVars_; }
    
  /**
   * Return the REVERSE map between the original variables and the variables in the relaxation
   */
  std::map <ConstVariablePtr, ConstVariablePtr> getRevOriginalVariablesMap() {return rev_oVars_; }
    
  /**
   * Return the map between the original variables and the copy variables for
   * its exponents
   */
  std::map<ConstVariablePtr, std::vector<ConstVariablePtr> > getNewCopyVariables() 
    {return newCopyVariables_; }
    
  /**
   * Retrun the map of monomial terms
   */
  std::map <VarIntMap, ConstVariablePtr> getMonomialTerms() {return monomial_terms_; }
    
  /**
   * Return the map of bilinear terms (map of a variable pair and a
   * substiture variable)
   */
  std::map <ConstVariablePair, ConstVariablePtr> getRevBilinearTerms() 
    {return rev_blterms_; }


  /// Can not return any cuts for this case.
  void separate(ConstSolutionPtr, NodePtr , RelaxationPtr , CutManager *,  
                SolutionPoolPtr, bool *, SeparationStatus *) {};
    
  virtual ModificationPtr getBrMod(BrCandPtr , DoubleVector &, 
                                   RelaxationPtr , BranchDirection );
    
  virtual Branches getBranches(BrCandPtr cand, DoubleVector &x , 
                               RelaxationPtr rel, SolutionPoolPtr s_pool);

  // presolve.
  virtual SolveStatus presolve(PreModQ *, bool *changed) {return Finished;};
    
  // Implement Handler::presolveNode()
  virtual bool presolveNode(RelaxationPtr, NodePtr, SolutionPoolPtr, ModVector &,
                    ModVector &)
  {return false;};

    
  // Write name
  std::string getName() const;
    
protected:

  /// Environment
  EnvPtr env_;

  // /**
  // The problem for which the handler was created.
  //XXX This should be const!
  // */
  ConstProblemPtr problem_;

  ProblemPtr workingProblem_;
    
    
    
private:
    
  LPEnginePtr lp_engine_;
  int linearizationCnt_;
    
  // Whether the objective was nonlinear and was moved to the constraints or not
  bool objModified_;

  // (Absolute?) error tolerance for branching.
  double eTol_;

  // max_pow maps each variable with its highest power that appears 
  // in the polynomial or quadratic part
  std::map <ConstVariablePtr, UInt> max_pow_;
  std::map <ConstVariablePtr, UInt>::iterator max_pow_it_;

  LoggerPtr logger_;

  // The map for the original variables and their copy in the relaxation (x, newVar)
  std::map <ConstVariablePtr, ConstVariablePtr> oVars_;
    
  // The map for the original variables and their copy in the relaxation -- reverse map (newVar, x)
  std::map <ConstVariablePtr, ConstVariablePtr> rev_oVars_;

  // The square terms <x1,(y,x_1)>
  std::map <ConstVariablePtr, ConstVariablePair> sqterms_;
    
  // The square terms -- reverse map
  std::map <ConstVariablePair, ConstVariablePtr> rev_sqterms_;
    
  // ******
  // The bilinear terms of the original constraints
  std::map <ConstVariablePtr, ConstVariablePair> blterms_cons_;
  std::map <ConstVariablePair, std::vector<double> > blterms_cons_coef_;

  // The bilinear terms of the original constraints -- reverse map
  std::map <ConstVariablePair, ConstVariablePtr> rev_blterms_cons_;

  // The multilinear terms of the original constraints
  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms_cons_;
  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator mlterms_cons_it_;
  std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_cons_coef_;
    
  // The multilinear terms of the original constraints -- reverse map
  std::map <std::vector<ConstVariablePtr>, ConstVariablePtr > rev_mlterms_cons_;
  // ******


  // ******
  // The bilinear terms of the original objective
  std::map <ConstVariablePtr, ConstVariablePair> blterms_obj_;
  std::map <ConstVariablePair, std::vector<double> > blterms_obj_coef_;

  // The bilinear terms of the original constraints -- reverse map
  std::map <ConstVariablePair, ConstVariablePtr> rev_blterms_obj_;

  // The multilinear terms of the original constraints
  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms_obj_;
  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator mlterms_obj_it_;
  std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_obj_coef_;
    
  // The multilinear terms of the original constraints -- reverse map
  std::map <std::vector<ConstVariablePtr>, ConstVariablePtr > rev_mlterms_obj_;
  // ******


  // All the bilinear terms 
  std::map <ConstVariablePtr, ConstVariablePair> blterms_;
  std::map <ConstVariablePair, std::vector<double> > blterms_coef_;

  // All the bilinear terms -- reverse map
  std::map <ConstVariablePair, ConstVariablePtr> rev_blterms_;
    
  // All the multilinear terms
  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms_;
  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator mlterms_it_;
  std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_coef_;

  // All the multilinear terms -- reverse map
  std::map <std::vector<ConstVariablePtr>, ConstVariablePtr > rev_mlterms_;

  // The monomial terms
  std::map <VarIntMap, ConstVariablePtr> monomial_terms_;
  std::map <VarIntMap, ConstVariablePtr>::iterator monomial_terms_it_;

  // The groups
  std::vector<std::vector<ConstVariablePtr> > groups_;

  // All the lambda variables
  std::vector<std::vector<ConstVariablePtr> > all_lambdas_;

  // The copy variables for the variables with exponent
  std::map<ConstVariablePtr, std::vector<ConstVariablePtr> > newCopyVariables_;
  std::map<ConstVariablePtr, std::vector<ConstVariablePtr> >::iterator newCopyVariables_it_;


  /// A (horrible) helper function to clean all the containers from the last time relax() was called
  void clearAllContainers();

  ///  A helper function to make McCormick.  This works for multilinear -- (recursive) as well

  void makeMcCormick(RelaxationPtr relaxation, bool &isInfeasible);
    
  ///  A helper function to make Convex Hull
  void makeGroupedConvexHull(RelaxationPtr relaxation, bool &isInfeasible, 
                             int groupStrategy, bool objModified);
    
  ///  A function that groups the variables to form the convex hull
  void makeGroups(std::map <ConstVariablePtr, ConstVariablePair> blterms, 
                  std::map <ConstVariablePair, ConstVariablePtr> rev_blterms, 
                  std::map <ConstVariablePair, std::vector<double> > blterms_coef,
                  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms, 
                  std::map <std::vector<ConstVariablePtr>, ConstVariablePtr> rev_mlterms, 
                  std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_coef,
                  std::map <ConstVariablePtr, ConstVariablePair> blterms_obj, 
                  std::map <ConstVariablePair, ConstVariablePtr> rev_blterms_obj, 
                  std::map <ConstVariablePair, std::vector<double> > blterms_obj_coef,
                  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms_obj, 
                  std::map <std::vector<ConstVariablePtr>, ConstVariablePtr> rev_mlterms_obj, 
                  std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_obj_coef,
                  std::map <ConstVariablePtr, ConstVariablePair> blterms_cons, 
                  std::map <ConstVariablePair, ConstVariablePtr> rev_blterms_cons, 
                  std::map <ConstVariablePair, std::vector<double> > blterms_cons_coef,
                  std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms_cons, 
                  std::map <std::vector<ConstVariablePtr>, ConstVariablePtr> rev_mlterms_cons, 
                  std::map <std::vector<ConstVariablePtr>, std::vector<double> > mlterms_cons_coef,
                  std::vector <std::vector<ConstVariablePtr> > &groups,                               
                  int groupStrategy);


  /// A function that gets the map of multilinear and bilinear terms and returns all the terms
  void getMultilinearTerms(std::map <ConstVariablePtr, ConstVariablePair> blterms,
                           std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms,
                           UInt maxGroupSize,
                           std::vector<std::vector<ConstVariablePtr> > &terms);

  /**
   * get the coefficients of the terms
   */
  void getMultilinearTermsCoef(std::map <ConstVariablePair, double> blterms_coef,
                               std::vector<std::vector<ConstVariablePtr> > terms,
                               std::vector<double> &termsCoeff);


  /**
   * Get the variables that appear in the multilinear terms
   */
  void getMultilinearVariables(std::map <ConstVariablePtr, ConstVariablePair> blterms,
                               std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms,
                               std::vector<ConstVariablePtr> &mlVars);

  /**
   * This function takes the terms that appear in the problem
   * and groups them using the following strategy:
   * at each iteration it takes a term
   * if the term does not already appear in the current groups
   *   it's added it to the group which has the highest intersection with it
   *   if all intersections are empty, then it's added as a new group.
   * the term is also added to other groups with nonempty intersection with the
   * term based on a probablity depending of rho
   */
  void groupUsingIntersection(UInt gs, double rho,
                              std::vector<std::vector<ConstVariablePtr> > terms,
                              std::vector <std::vector<ConstVariablePtr> > &groups);

  /**
   *
   */
  void groupUsingIntersection2(UInt gs,
                               std::vector<std::vector<ConstVariablePtr> > terms,
                               std::vector <std::vector<ConstVariablePtr> > &groups);

  void groupTermByTerm(std::map <ConstVariablePtr, ConstVariablePair> blterms,
                       std::map <ConstVariablePtr, std::vector<ConstVariablePtr> > mlterms,
                       std::vector <std::vector<ConstVariablePtr> > &groups);

  void groupUsingDensest(UInt gs, 
                         std::map <ConstVariablePair, std::vector<double> > bl_coef,
                         std::map <std::vector<ConstVariablePtr>, std::vector<double> > ml_coef,
                         std::vector <ConstVariablePtr> vars,
                         std::vector <std::vector<ConstVariablePtr> > &groups);
      
  void groupUsingDensest2ND(UInt gs, 
                            std::map <ConstVariablePair, std::vector<double> > bl_coef,
                            std::map <std::vector<ConstVariablePtr>, std::vector<double> > ml_coef,
                            std::vector <ConstVariablePtr> vars,
                            std::vector <std::vector<ConstVariablePtr> > &groups,
                            UInt maxGroupsCnt);

  /**
   * Groups by adding variables which make the sum of the abs value of coefficients of
   * the terms appearing in the group maximum
   */
  void groupUsingCoef(UInt gs,
                      std::map <ConstVariablePair, std::vector<double> > bl_coef,
                      std::map <std::vector<ConstVariablePtr>, std::vector<double> > ml_coef,
                      std::vector <ConstVariablePtr> graphVars,
                      std::vector <std::vector<ConstVariablePtr> > &groups);

  void findSubgraphDensity(std::vector<std::vector<ConstVariablePtr> > terms,
                           double** term_var_matrix,
                           std::vector<ConstVariablePtr> nodes,
                           std::vector<int> nodesInd,
                           double &density);

  void findDensestSubgraph(UInt gs,
                           std::vector<ConstVariablePtr> vars,
                           std::vector<std::vector<ConstVariablePtr> > terms,
                           double** term_var_matrix,
                           std::vector<int> component,
                           std::vector<ConstVariablePtr> &densest,
                           std::vector<int> &densestInd,
                           bool &isCompEmpty);

  void findGraphComponents(int varsCnt,
                           int termsCnt,
                           double** term_var_matrix, 
                           std::vector<std::vector<int> > &components);

  void findGroupDensity(std::vector<std::vector<ConstVariablePtr> > terms,
                        std::vector <std::vector<ConstVariablePtr> > groups,
                        std::vector <int> &density);

  void findGroupCoef(std::vector<std::vector<ConstVariablePtr> > terms,
                     std::vector <double> termsCoef,
                     std::vector <std::vector<ConstVariablePtr> > groups,
                     std::vector <double> &sumCoef);

  void removeSubgraphEdges(int termsCnt,
                           double** &term_var_matrix,
                           std::vector<std::vector<ConstVariablePtr> > terms,
                           std::vector<ConstVariablePtr> nodes,
                           std::vector<int> nodesInd);

  void reduceSubgraphEdges(int termsCnt,
                           double** &term_var_matrix,
                           std::vector<std::vector<ConstVariablePtr> > terms,
                           std::vector<ConstVariablePtr> nodes,
                           std::vector<int> nodesInd);

  /**
   * This function takes a vector of terms and a grouping and returns
   * the terms that have been covered in the grouping k times
   */
  void termsAppearKtimes(std::vector<std::vector<ConstVariablePtr> > terms,
                         std::vector <double> termsCoef,
                         std::vector <std::vector<ConstVariablePtr> > groups,
                         int k,
                         std::vector<std::vector<ConstVariablePtr> > &termsk,
                         std::vector <double> &termsKCoef,
                         std::vector<int> &termRep);

  /**
   * Counts the number of appearance of each term in the grouping
   */
  void countTermsAppearance(std::vector<std::vector<ConstVariablePtr> > terms,
                            std::vector <std::vector<ConstVariablePtr> > groups,
                            std::vector<int> &termRep);

  /**
   * In each term of a polynomial function, if a binary variable takes a 
   * power bigger than 1 make that power 1
   */
  void makeBinaryVariablesPowers1(FunctionPtr &f, PolyFunPtr &p, 
                                  QuadraticFunctionPtr &q, LinearFunctionPtr &l);    
    
  /*
 ///  A helper function to make Convex Hull
 void makeGroupedConvexHull(RelaxationPtr relaxation, bool &isInfeasible, 
 std::vector<bool> & c_list, int groupStrategy);
  */
    
  /// A function that finds all the extreme points for a set of variables
  void allExtreme(std::vector<int> &S, std::vector<double> &lb, 
                  std::vector<double> &ub, 
                  std::vector<std::vector<double> > &E);
    
  /// 
  void visit(std::vector<int> &S, std::vector<double> &lb, 
             std::vector<double> &ub, UInt ix, std::vector<double> &val, 
             std::vector<std::vector<double> > &E);
    
    
};
typedef boost::shared_ptr<MultilinearHandler> MultilinearHandlerPtr;
typedef boost::shared_ptr<const MultilinearHandler> MultilinearConstHandlerPtr;
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
