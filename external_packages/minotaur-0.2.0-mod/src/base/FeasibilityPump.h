//
//     MINOTAUR -- It's only 1/2 bull  
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//
/** 
 * \file FeasibilityPump.h
 * \brief Declare the class FeasibilityPump derived from base class Heuristic.
 * \author Jayash Koshal, Argonne National Laboratory
 */

#ifndef FEASIBILITYPUMP_H
#define FEASIBILITYPUMP_H

#include "Heuristic.h"

namespace Minotaur {
  class Engine;
  class Problem;
  typedef boost::shared_ptr<Engine> EnginePtr;
  typedef boost::shared_ptr<Problem> ProblemPtr;

  /// statistics for Feasibility Pump heuristic
  struct FeasPumpStats {
    UInt numNLPs;         /// Number of NLPs solved in the heuristic
    UInt errors;          /// Number of errors in NLP solved
    UInt numCycles;       /// Number of time the same solution was obtained
    double time;          /// Total time taken by the heuristic
    double bestObjValue;  /// Best objective value of the feasible solution
  };

  /**
   * \brief Feasibility Pump for MINLPs
   *
   * A Feasibility Pump heuristic used to find solutions for Mixed
   * Integer NLPs by solving a relaxed NLP using an NLP engine. After an
   * initial solve a sequence of points are found which are closest to
   * the nearest integer solution obtained by rounding the previous
   * point.
   */

  class FeasibilityPump : public Heuristic {

  public:

    /// default constructor
    FeasibilityPump(EnvPtr env, ProblemPtr p, EnginePtr e);

    /// constructor for derived class
    FeasibilityPump(EnvPtr env, ProblemPtr p, EnginePtr nlpe, EnginePtr
                    e);

    /// default destructor
    virtual ~FeasibilityPump();

    /// call to the heuristic
    void solve(NodePtr node, RelaxationPtr rel, SolutionPoolPtr s_pool);

    /// write statistic to the logger
    void writeStats(std::ostream &out) const;

  protected:

    /// Message name for the heuristic
    const static std::string me_;

    /// Binary/integer variables present in the problem
    VarVector bins_;

    /// Pointer to the engine to be used to solve the problem
    EnginePtr e_;

    /// Pointer to the environment
    EnvPtr env_;

    /// Vector for hash value of solutions
    DoubleVector hashVal_;

    /// Tolerance for a number to be considered as an integer
    double intTol_;

    /// Pointer to the logger
    LoggerPtr logger_;

    /// Number of variables to be flipped if cycling is detected
    UInt nToFlip_;

    /// Pointer to the problem being solved
    ProblemPtr p_;

    /// A random vector for inner product with the solution
    DoubleVector random_;

    /// Vector of rounded solution
    DoubleVector roundedSol_;

    /// Statistics for the Feasibility Pump heuristic
    FeasPumpStats* stats_;

    /// Timer of the heuristic
    Timer* timer_;

    /** 
     * \brief Function to construct/update the objective function
     *
     * \param[in] prob Pointer to the cloned problem
     * \param[in] sol Pointer to the solution of relaxation that was
     *            previously solved.
     *
     * This function selects the variable which are fractional and
     * construct the objective function out of such variables. The
     * new objective function replaces the objective of the problem 
     * passed
     */
    virtual void constructObj_(ProblemPtr prob, ConstSolutionPtr sol);

    /**
     * \brief Function to convert a solution of the cloned problem to 
     * that of an original problem
     *
     * \param[in] s_pool Pointer to solution pool of original problem
     * \param[in] sol Solution pointer to the modified (cloned) 
     * problem
     *
     * The binary variables are fixed by changing their bounds in the
     * original problem and then the original problem is resolved to
     * obtain a feasible solution. Bounds are relaxed at the end. 
     */
    void convertSol_(SolutionPoolPtr s_pool, ConstSolutionPtr sol); 

    /**
     * \brief A search function for detection of cycling
     *
     * \param[in] find_value Value to be located in a vector containing
     * hash values of already visited rounded solutions
     *
     * \return true If the point already exist implying cycling 
     * and false otherwise
     */
    bool cycle_(double find_value);

    /** 
     * \brief A function to hash the solutions
     *
     * \return The hash value of the solution
     *
     * A hashing function which calculates the inner product of vector
     * x with a FIXED random vector in (0,1]^n. This function uses the
     * roundedSol_ vector and random_ vector to calculate the hash value
     */
    double hash_();

    /** 
     * \brief Function to implement the Feasibility Pump method
     *
     * \param[in] Constant pointer to primal solution
     * \param[in] Pointer to solution pool
     *
     */
    virtual void implementFP_(const double* x, SolutionPoolPtr s_pool);

    /** 
     * \brief Function to check the integrality of a variable
     *
     * \param[in] Constant pointer to the primal solution
     *
     * return true if there is any fractional variable else false
     */
    bool isFrac_(const double* x);


    /** 
     * \brief A function to perturb the rounded solution in case of 
     * cycling
     *
     * \param[in] hash_val Hash value of the current rounded solution
     * \param[in] n_to_flip Number of variables to be flipped
     *
     * This function uses selectToFlip_ to get candidates for flipping
     * till a new rounded solution is obtained which is not visited 
     * earlier. Cycling_ is used to detect earlier existence. The
     * parameter n_to_flip is passed by the calling method is usually
     * an indicator of integer infeasibilities of the solution.
     */
    void perturb_(double hash_val, UInt n_to_flip);

    /**
     * \brief A function to restore the upper and lower bounds of the
     * problem.
     *
     * \param[in] LB_copy. Pointer to an array of lower bound. 
     * \param[in] UB_copy. Pointer to an array of upper bound. 
     * \param[in] numvars. Number of variables in the problem.
     *
     */
    void restoreBounds_(double* LB_copy, double* UB_copy, UInt vars);

    /**
     * \brief A function to save the upper and lower bounds of the
     * problem.
     *
     * \param[in] LB_copy. Pointer to an array of lower bound. Space has
     * to be allocated.
     * \param[in] UB_copy. Pointer to an array of upper bound. Space has
     * to be allocated.
     * \param[in] numvars. Number of variables in the problem.
     *
     */
    void saveBounds_(double* LB_copy, double* UB_copy, UInt vars);

    /**
     * \brief A funtion to randomly select "n" binary/integer 
     * variable to be flipped.
     *
     * \param[in] n_to_flip Number of variables to be flipped
     *
     * \return The vector of pointer to the variables selected as
     * candidate
     * 
     * This function uses random sampling for selection of variables  as
     * a candidate for be flipped from 0 to 1 or 1 to 0.
     */
    VarVector selectToFlip_(UInt n_to_flip);

    /**
     * \brief Function to decide whether to use Feasibility Pump
     *
     * return true or false. 
     *
     * We decide not to use Feasibility Pump is problem has integer
     * variables and/or if the problem is nonlinear (constraints or
     * objective)
     */
    virtual bool shouldFP_();

  };

  typedef boost::shared_ptr<FeasibilityPump> FeasPumpPtr;

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
