// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file Relaxation.h
 * \brief Declare the class Relaxation for storing and manipulating
 * relaxations.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURRELAXATION_H
#define MINOTAURRELAXATION_H

#include "Problem.h"

namespace Minotaur {


/**
 * Relaxation is a derived class of Problem. A relaxation is what is
 * actually solved at each iteration of an algorithm. 
 * 
 * A Relaxation need not be a relaxation of the original problem. It could
 * be a relaxation of the problem being solved at a given node. Since a
 * relaxation could be as big as the 
 * the original problem (and even much bigger), we should take care to have
 * as few copies of a relaxation around as possible. 
 * 
 * A relaxation is created from a Problem by using Handlers. It is assumed
 * that a Problem is in its standard form. Each constraint is relaxed by one
 * or more handlers. A handler may break the constraint into two parts,
 * adding a new variable to the relaxed problem. For instance: a constraint
 * \f$x^2 + y^2 - z^2 \leq 5\f$ can be broken as \f$x^2 + y^2 - v \leq 0\f$,
 * \f$-z^2 + w \leq 0\f$, \f$v - w \leq 5\f$.
 * 
 * Besides having constraints, variables, functions, jacobians and hessians
 * much like the Problem, a relaxation can have additional data structures
 * that may be used to solve the problem, e.g. reduced costs, conflict
 * graphs, implications, valid inequalities, symmetry groups etc.
 * Typically one would input a Problem and ask Minotaur to solve it using,
 * say, branch-and-bound. 
 * 
 * Since a Relaxation is obtained by reformulating the original Problem, it
 * should thus have methods to translate and return the solution values of
 * the variables of the original problem. It should also keep a pointer to
 * the original Problem from which it is derived. Some information stored in
 * a relaxation may be global, i.e. valid for the original problem. Other
 * information may be local, i.e. valid for the current iteration of
 * branch-and-bound only.
 */
class Relaxation : public Problem {
    
public:
  /// Default constructor.
  Relaxation();

  /**
   * Construct a relaxation from an Problem. In the default
   * implementation, as many new variables and
   * constraints as in the Problem are created. If nonlinear function can't be
   * cloned, their pointers are saved.  Everything else in the constraint
   * (bounds, sense, map etc.) are copied.
   */
  Relaxation(ProblemPtr problem);

  /// Destructor. No need yet. Use ~Problem().
  ~Relaxation() {};

  VariablePtr getOriginalVar(VariablePtr r_var);
  
  VariablePtr getRelaxationVar(VariablePtr p_var);

  void setProblem(ProblemPtr p);

protected:
  /// Pointer to the original problem.
  ConstProblemPtr p_;
};

typedef boost::shared_ptr<Relaxation> RelaxationPtr;
typedef boost::shared_ptr<const Relaxation> ConstRelaxationPtr;  
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
