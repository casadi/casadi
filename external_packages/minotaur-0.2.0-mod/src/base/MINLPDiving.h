//
//     MINOTAUR -- It's only 1/2 bull  
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//
 
/**
 * \file MINLPDiving.h
 * \brief Define the class MINLPDiving derived from base class Heuristics.
 * \author Jayash Koshal, Argonne National Laboratory
 *
 * Declares the class MINLPDiving.
 */
 
 #ifndef MINOTAURMINLPDIVING_H
 #define MINOTAURMINLPDIVING_H
 
 #include "Heuristic.h"
 #include <vector>
 #include <stack>


 namespace Minotaur {
   class Engine;
   class LinearHandler;
   class Problem;
   class Solution;
   class Timer;
   class VarBoundMod;
   typedef boost::shared_ptr<const Solution> ConstSolutionPtr;
   typedef boost::shared_ptr<Engine> EnginePtr;
   typedef boost::shared_ptr<Problem> ProblemPtr;
   typedef boost::shared_ptr<VarBoundMod> VarBoundModPtr;

   /// Direction of rounding
   typedef enum {
     Floor,        /// Round to the smallest integer
     Ceil,         /// Round to the largest integer
     Nearest,      /// Round to the nearest integer
     Farthest      /// Round to the farthest integer
   } Direction;

   /// Order of rounding: least fractional or most fractional
   typedef enum {
     Least,     /// Select the variable with Least fractional part
     Most       /// Select the variable with Most fractional part
   } Order;

   /// Type of score evaluation for fractional variable
   typedef enum {
     Fractional,     /// Score of variable is the fractional part`
     VectorLength,   /// Score of variable is the Vector Length
     LexBound,       /// Score of variable is the index of the variable
     ReducedCost     /// Score of variable is the reduced cost
   } Scoretype;

   /**
    * \brief A statistic struct for MINLP Diving heuristic
    */
   struct DivingheurStats {
     UInt numNLPs[4];       /// NLPs solved in each selection method.
     UInt numInfeas[4];     /// Infeasible NLPs solved. 
     UInt errors[4];        /// NLPs that encountered errors.
     UInt numSol[4];        /// Solutions obtained.
     double time[4];        /// Time for each selection method.
     UInt iterations[4];    /// NLP iterations.
     UInt totalNLPs;        /// NLPs solved.
     UInt totalSol;         /// Solutions found.
     UInt numLocal;         /// Local optimal solutions obtained.
     double best_obj_value; /// Best objective value of feasible solution.
     double totalTime;      /// Total time taken in Diving heuristics.
   };

   /**
    * \brief Diving heuristif for MINLPs.
    *
    * A Diving Heuristic used to find solutions for Mixed Integer NLPs 
    * by solving the Relaxed NLP using an NLP engine. The engine is 
    * called once initially to generate a solution  which is rounded 
    * and used for diving. 
    */

   class MINLPDiving : public Heuristic {
   public:

     /// default constructor
     MINLPDiving(EnvPtr env, ProblemPtr p, EnginePtr e);

     /// default destructor
     ~MINLPDiving();

     /// call to heuristic
     void solve(NodePtr node, RelaxationPtr rel, SolutionPoolPtr s_pool); 

     /// writing the statistics to the logger
     void writeStats(std::ostream &out) const;

   protected:

     const static std::string me_;

     /// Average of dual variables from the previous iterations
     /// which is to be used to reduced cost diving
     DoubleVector avgDual_;

     /// Engine being used to solve the problem
     EnginePtr e_;

     /// Environment
     EnvPtr env_;

     /// Gradient of objective function for vector length diving
     double* gradientObj_;

     /// If a value is this far from an integer, it is considered
     /// integer. 
     double intTol_;

     /// Mods implied by bound changes in the previous node's presolve. 
     /// We only store mods of one node only.
     ModVector lastNodeMods_;

     /// Linear Handler for presolving.
     LinearHandler *lh_;

     /// Logger
     LoggerPtr logger_;

     /// Maximum number of NLPs allowed for diving heuristic
     UInt maxNLP_;

     /// Maximum number of feasible solution required from heuristic
     UInt maxSol_;

     /** A stack of modification pointer. All modification are stored
      * in stack and for backtracking the stack is unwinded to restore
      * feasibility
      */
     std::stack<VarBoundModPtr> mods_;

     /// Number of method for selection of variables
     UInt nSelector_;

     /// Problem to be solved
     ProblemPtr p_;

     /// Violated variable fraction part list
     DoubleVector score_;

     /// Statistics for the heuristic
     DivingheurStats *stats_;

     /// Timer for this heuristic
     Timer* timer_;

     /// violated variable index list
     UIntVector violated_;

     typedef UInt (MINLPDiving::*FuncPtr) (UInt numfrac, 
                                           const double* x, Direction d, Order o);

     /**
      * \brief Backtracking method
      *
      * \param[in] n_flipped Number of bound changes made in previous
      * dive
      *
      * All the changes are made by unwinding the modification stack
      * n_flipped number of times
      */
     void backtrack_(UInt n_flipped);

     /** 
      * \brief Fractional selection method for fractional variable
      * 
      * \param[in] numfrac Number of fractional variables in current
      * solution
      * \param[in] x Constant pointer to primal solution
      * \param[in] d Direction of rounding
      * \param[in] o Order for selection of fractional variables
      */
     UInt FracBounds_(UInt numfrac, const double* x, 
                      Direction d, Order o);

     /**
      * \brief Get the score of solution
      *
      * \param[in] x Const pointer to the primal solution
      * \param[in] s Scoretype
      *
      * The score of the solution is determined based on one of the
      * predefined scoretypes.
      */
     void getScore_(const double* x, Scoretype s);

     /**
      * \brief Function to implement diving
      *
      * \param[in] i Index of the method number
      * \param[in] x Const pointer to the root node primal solution
      * \param[in] s_pool Pointer to the solution pool
      *
      * Method for selection of fractional variable as candidate for
      * rouding are chosen for diving. Changes made to the problem are
      * stored in the stack of modification.
      */
     void implementDive_(int i, const double*x, SolutionPoolPtr s_pool);

     /** 
      * \brief The number of fractional variables in current solution
      *
      * \param[in] x Const pointer to primal solution
      *
      * \return Number of fractional variables
      */
     UInt isFrac_(const double* x);

     /** 
      * \brief Lexicographic selection method for fractional variable
      * 
      * \param[in] numfrac Number of fractional variables in current
      * solution
      * \param[in] x Constant pointer to primal solution
      * \param[in] d Direction of rounding
      * \param[in] o Order for selection of fractional variables
      */
     UInt LexBounds_(UInt numfrac, const double* x, 
                     Direction d, Order o);

     /** 
      * \brief Reduced cost diving selection method for fractional variable
      * 
      * \param[in] numfrac Number of fractional variables in current
      * solution
      * \param[in] x Constant pointer to primal solution
      * \param[in] d Direction of rounding
      * \param[in] o Order for selection of fractional variables
      */
     UInt ReducedCost_(UInt numfrac, const double* x, 
                       Direction d, Order o);

     /**
      * \brief Restore the bounds for the problem
      *
      * \param[in] LB_copy Pointer to saved bounds
      * \param[in] UB_copt Pointer to saved bounds
      * \param[in] vars Number of variables in the problem
      */
     void restoreBounds_(double* LB_copy, double* UB_copy, UInt vars);

     /**
      * \brief Rounding a value in a given direction
      *
      * \param[in] value Fractional value to be rounded
      * \param[in] d Direction to be used for rounding
      *
      * \return rounded value of the fractional variable
      */
     double rounding_(double value, Direction d);

     /** 
      * \brief Save bounds of the problem
      *
      * \param[in] LB_copy Pointer to bounds, space has to be allocated
      * \param[in] UB_copt Pointer to bounds, space has to be allocated
      * \param[in] vars Number of variables in the problem
      */
     void saveBounds_(double* LB_copy, double* UB_copy, UInt vars);

     /**
      * \brief Select the method, ordering and direction.
      *
      * \param[in] i Index of the method number
      * \param     d Reference to direction of rouding
      * \param     o Reference to order for selecting variable
      *
      * \return    FuncPtr Address of the selected Diving method
      */
     FuncPtr selectHeur_(int i, Direction &d, Order &o);

     /** 
      * \brief Function to decide on diving or not
      *
      * return true or false 
      *
      * We decide not to dive is the problem size is large and hessian
      * is dense
      */
     bool shouldDive_();

     /** 
      * \brief Sort the score
      * 
      * \param[in] left Beginning of the container
      * \param[in] right End of the container
      *
      * If the number of non-zero in hessian or number of bin+int <20
      * don't dive
      */
     void sort_(UInt left, UInt right);

     /**
      * \brief Update the average of dual multiplier
      * 
      * \param[in] sol Constant Pointer to the solution
      */ 
     void updateAvgDual_(ConstSolutionPtr sol);

     /** 
      * \brief Vector Length selection method for fractional variable
      * 
      * \param[in] numfrac Number of fractional variables in current
      * solution
      * \param[in] x Constant pointer to primal solution
      * \param[in] d Direction of rounding
      * \param[in] o Order for selection of fractional variables
      */
     UInt VectorLength_(UInt numfrac, const double* x, 
                        Direction d, Order o);

     /// Function to decide on vector length diving
     bool vectorFlag_(UInt min_vlength);

   };

   typedef boost::shared_ptr<MINLPDiving> MINLPDivingPtr;
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
