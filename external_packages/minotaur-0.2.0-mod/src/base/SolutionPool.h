// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file SolutionPool.h
 * \brief Declare a container for storing solutions and their qualities.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */



#ifndef MINOTAURSOLUTIONPOOL_H
#define MINOTAURSOLUTIONPOOL_H

#include "Problem.h"
#include "Solution.h"

namespace Minotaur {

  class Environment;
  class Timer;
  typedef boost::shared_ptr<Environment> EnvPtr;

  class SolutionPool {
  public:
    /// Default constructor.
    SolutionPool();

    /// Construct a solution pool of a given size for a given problem
    SolutionPool(EnvPtr env, ProblemPtr problem, UInt limit=100);

    /// Add Solution to the pool
    void addSolution(ConstSolutionPtr);

    /// Get number of solutions in the pool
    UInt getNumSols() const;

    /// Get number of solutions in the pool
    UInt getNumSolsFound() const;

    /// Get the limit on the number of solutions in the pool
    UInt getSizeLimit() const;

    /// Put a limit on the number of solutions in the pool
    void setSizeLimit(UInt limit);

    /// Get iterator for the first solution ...
    SolutionIterator solsBegin() { return sols_.begin(); }

    /// ... and the end.
    SolutionIterator solsEnd() { return sols_.end(); }

    /// Create a solution from a double array and add Solution to the pool.
    void addSolution(const double *x, double obj_value);

    /**
     * Get a solution with the best objective function value. Return NULL if
     * the pool is empty.
     */
    SolutionPtr getBestSolution();

    /// Get the best objective function value
    double getBestSolutionValue() const;

    /// Write statistics to the outstream.
    void writeStats(std::ostream &out) const; 

  private:
    /// The solutions are stored in a vector. 
    std::vector<SolutionPtr> sols_;

    /**
     * The best solution in terms of objective function value. In case of tie,
     * the most recently found one.
     */
    SolutionPtr bestSolution_;

    /// For logging.
    const static std::string me_;

    /// The number of solutions added to the pool.
    UInt numSolsFound_;

    /// Problem for which we are saving solutions
    ProblemPtr problem_;

    /// The limit on number of solutions in the pool.
    UInt sizeLimit_;

    /// Time when the best solution is found.
    double timeBest_;

    /// Time when the first solution is found.
    double timeFirst_;

    /// Global timer.
    const Timer* timer_;

  };

  typedef boost::shared_ptr<SolutionPool> SolutionPoolPtr;
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
