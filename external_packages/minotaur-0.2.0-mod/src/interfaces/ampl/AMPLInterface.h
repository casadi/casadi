// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file AMPLInterface.h
 * \brief Declare the AMPLInterface class fo reading problems from AMPL.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURAMPLINTERFACE_H
#define MINOTAURAMPLINTERFACE_H

#include "Types.h"

#include "asl.h"
#include "nlp.h"
//#undef filename

namespace Minotaur {
class   CGraph;
class   CNode;
class   Environment;
class   LinearFunction;
class   PolynomialFunction;
class   Problem;
class   QuadraticFunction;
class   Solution;
typedef boost::shared_ptr<CGraph> CGraphPtr;
typedef boost::shared_ptr<Environment> EnvironmentPtr;
typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;
typedef boost::shared_ptr<PolynomialFunction> PolyFunPtr;
typedef boost::shared_ptr<Problem> ProblemPtr;
typedef boost::shared_ptr<QuadraticFunction> QuadraticFunctionPtr;
typedef boost::shared_ptr<const Solution> ConstSolutionPtr;
}

namespace MINOTAUR_AMPL {

class AMPLNonlinearFunction;
typedef boost::shared_ptr<AMPLNonlinearFunction> AMPLNlfPtr;

/// What kind of ASL reader is used to read the .nl file.
typedef enum {
  FReader,    /// No derivatives, linear objective and constraints only.
  FGReader,   /// first derivatives. This reader will be used to read the 
              /// expression tree.
  FGHReader,  /// First derivatives and Hessian-vector products. Need a 
              /// loop for computing full hessian.
  PFGReader,  /// First derivatives and partially separable structure.
  PFGHReader  /// First and second derivatives and partially separable
              /// structure. Can compute full hessian. This reader is used
              /// when NLPs are solved by calling ASL's evaluation routines.
} ReaderType;


class AMPLInterface;
typedef AMPLInterface * AMPLInterfacePtr;


// ------------------------------------------------------------------------ //
// ------------------------------------------------------------------------ //

/**
 * \brief Interface to read ampl models using AMPL Solver Library.
 * AMPLInterface class provides methods to read and evaluate an instance
 * from AMPL generated .nl files. It also provides methods to evaluate
 * functions, gradients, jacobians and hessians.
 */
class AMPLInterface {
public:
  /// Constructor.
  AMPLInterface(Minotaur::EnvPtr env, std::string solver="minotaurampl");

  /// Destroy.
  ~AMPLInterface();

  /// Free the ASL data structure myAsl_;
  void freeASL();

  /// Get pointer to the ASL structure myAsl_
  ASL * getAsl();

  /// Get the initial point provided in the .nl file
  const double * getInitialPoint() const;

  /// Get the number of defined variables.
  Minotaur::UInt getNumDefs() const;

  /// What kind of reader from ASL was used to read the .nl file.
  ReaderType getReaderType();

  /// Read an instance from a .nl file 'fname'.
  Minotaur::ProblemPtr readInstance(std::string fname);

  /// Write the solution to the AMPL acceptable .sol file.
  void writeSolution(Minotaur::ConstSolutionPtr sol,
                     Minotaur::SolveStatus status);

  void writeProblem(std::ostream &out) const;
private:
  /// Log manager
  Minotaur::EnvPtr env_;

  /**
   * A map describing what kind of an AMPL OPCODE a given function has.
   * See comments for createFunctionMap_().
   */
  std::map <efunc*, int> functionMap_;

  /// Tolerance to check integrality (of powers for example).
  const double intTol_;

  /// Log manager
  Minotaur::LoggerPtr logger_;

  /// For logging.
  static const std::string me_;

  /// Pointer to ASL.
  ASL *myAsl_;

  /// Number of constraints. Needs to be stored for jacobian evaluations.
  int nCons_;

  /**
   * The number of defied variables. AMPL stores the number of 5 kinds of
   * defined variables: 
   * ncomo = number of defined variables in more than one objectives only
   * ncomc = number of defined variables in more than one constraints only
   * ncomb = number of defined variables in both objectives and constraints
   * ncomo1 = number of defined variables in only one objective
   * ncomc1 = number of defined variables in only one constraint.
   */
  int nDefVars_;

  /// ncomo+ncomc+ncomb
  int nDefVarsBco_;

  /// ncomo1+ncomc1
  int nDefVarsCo1_;

  /**
   * Total number of variables. It does not include the number of "defined
   * variables" or common expressions. 
   */
  int nVars_;

  /// Reader that was used to read the .nl file.
  ReaderType readerType_;

  /// A vector of variables that are in this instance.
  std::vector<Minotaur::VariablePtr> vars_;

  /// Tolerance to check if zero.
  const double zTol_;

  /// Add 'defined variables', i.e., variables for common-expressions.
  void addDefinedVars_(Minotaur::ProblemPtr instance);

  /// Add the i-th linear constraint from ASL.
  void addLinearConstraint_(int i, Minotaur::ProblemPtr problemPtr);

  /// Add the i-th linear objective from ASL.
  void addLinearObjective_(int i, Minotaur::ProblemPtr problemPtr);

  /**
   * Get linear terms in the i-th constraint and put them into lf. If lf
   * is NULL, new space is allocated.
   */
  void addLinearTermsFromConstr_(Minotaur::LinearFunctionPtr &lf, int i);

  /**
   * Get linear terms in the i-th objective and put them into lf. If lf
   * is NULL, new space is allocated.
   */
  void addLinearTermsFromObj_(Minotaur::LinearFunctionPtr &lf, int i);

  /// Add ampl specific options to Minotaur options.
  void addOptions_();

  /// Add the i-th non-linear constraint from ASL as a quadratic constraint
  void addQuadraticConstraint_(int i, Minotaur::ProblemPtr problemPtr);

  void addQuadraticDefCons_(Minotaur::ProblemPtr instance);
  void addQuadraticDefCons2_(Minotaur::ProblemPtr instance);

  /// Add the i-th non-linear objective from ASL as a quadratic objective
  void addQuadraticObjective_(int i, Minotaur::ProblemPtr problemPtr);

  /// Add variables in ASL to the instance 
  void addVariablesFromASL_(Minotaur::ProblemPtr problemPtr);

  /// Add 'SOS' constrants of type 1 and 2.
  void addSOS_(Minotaur::ProblemPtr instance);

  /**
   * Call ASL routines and create an instance. The expression tree is
   * traversed and the instance is stored in Minotaur data structures.
   * ASL is not needed once the instance is constructed.
   */
  Minotaur::ProblemPtr copyInstanceFromASL_();
  Minotaur::ProblemPtr copyInstanceFromASL2_();

  /**
   * AMPL has "N_OPS" different operations defined in opcode.hd in the ASL
   * folder. Each operation has an associated function (efunc), say f_i, i
   * = 0...(N_OPS-1). The data structure myAsl_->I->r_ops_ is an array of
   * these functions. Thus myAsl_->I->r_ops_[0] is a pointer to a function
   * corresponding to OPPLUS in opcode.hd. Similarly myAsl_->I->r_ops_[1]
   * for OPMINUS and so on.
   * 
   * Given an operation, say OPMINUS, it is easy to find the function that
   * performs that operation. However, AMPL expression-trees do not store
   * these opcodes. They just store function pointers. It is difficult to
   * know what operation is associated with a given function. So we need
   * to create a map from each such function to operation-code.
   */
  void createFunctionMap_();

  /// Find the type of input problem by querying ASL expression trees.
  Minotaur::ProblemType findProblemType_();

  /// Find variables in an expression and add them to the set 'vars'.
  void findVars_(expr *e_ptr, std::set<int> & vars);


  /// Get computational graph from an expression.
  Minotaur::CNode* getCGraph_(expr *e_ptr, Minotaur::CGraphPtr cgraph, 
                              Minotaur::ProblemPtr instance);

  /// Get the most general function type that describes all constraints.
  Minotaur::FunctionType getConstraintsType_();

  /// Get the function type of i-th constraint.
  Minotaur::FunctionType getConstraintType_(int i);

  /**
   * Get the function type that describes the defined variable of type
   * '1'
   */
  Minotaur::FunctionType getDef1ConstraintType_(int i);

  /// Get the function type that describes the defined variable 
  Minotaur::FunctionType getDefConstraintType_(int i);

  /**
   * Return function type of the ASL expression e_ptr. This function is
   * recursive as an expression can have another expression as its
   * operand.
   */
  Minotaur::FunctionType getExpressionType_(expr *e_ptr);

  /**
   * Call ASL routines and create an instance. Any nonlinear functions
   * will be evaluated by calling ASL.
   */
  Minotaur::ProblemPtr getInstanceFromASL_(
                                           std::vector<std::set<int> > &vids);
  /**
   * Return function type of the ASL sub-expresstion e_ptr when the root 
   * node has a plus or a minus. This function is
   * called when an opcode OPMULT is found.
   */
  Minotaur::FunctionType getMultExpressionType_(expr *e_ptr);

  /// Return the function type of objective.
  Minotaur::FunctionType getObjFunctionType_(Minotaur::UInt obj_index=0);

  /** 
   * When a solver is called from ampl, options are passed as
   * environment variables. For e.g. if the solver name is qg, and
   * option "presolve" is set to 1, then environment variable qg_presolve
   * is set to 1.
   */
  void getOptionsFromEnv_(std::string pre);

  /**
   * Return function type of the ASL sub-expresstion e_ptr when the root 
   * node has a plus or a minus. This function is
   * called when an opcode OPPLUS or OPMINUS is encountered.
   */
  Minotaur::FunctionType getPMExpressionType_(expr *e_ptr);

  void getPoly_(Minotaur::LinearFunctionPtr & lfPtr, 
                Minotaur::QuadraticFunctionPtr & qfPtr, 
                Minotaur::PolyFunPtr & pfPtr, 
                double & c, expr *e_ptr);

  /**
   * Return function type of the ASL sub-expresstion e_ptr when the root 
   * node has a (expr)^(constant) operation. This function is
   * called when an opcode OP1POW is found.
   */
  Minotaur::FunctionType getPow1ExpressionType_(expr *e_ptr);

  /**
   * Return function type of the ASL sub-expresstion e_ptr when the root 
   * node has a SUMLIST operation. This function is
   * called when an opcode SUMLIST is found.
   */
  Minotaur::FunctionType getSumlistExpressionType_(expr *e_ptr);

  /// Read a stub file using one of several readers provided by ASL
  void readFile_(std::string *fname, ReaderType readerType);

  Minotaur::ProblemPtr readInstanceASL_(std::string fname);
  Minotaur::ProblemPtr readInstanceCG_(std::string fname);

  void saveNlVars_(std::vector<std::set<int> > &vars);

  /**
   * \brief Complain about an unsupported operation.
   * \param [in] opcode Minotaur opcode that we are complaining about.
   */
  void unsupportedOp_(int opcode);

  /**
   * \brief Write an expression.
   * \param [in] e_ptr ASL expression that should be written.
   * \param [in] out output stream that we write to.
   */
  void writeExpression_(expr *e_ptr, std::ostream &out) const;

  /**
   * \brief Write the expression (function) of a constraint or objective.
   * \param [in] i The index of constraint that needs to be written. It is
   * ignored if is_obj is true.
   * \param [in] is_obj True if we want to write objective.
   * \param [in] out output stream that we write to.
   */
  void writeExpression_(int i, bool is_obj, std::ostream &out) const;

  /**
   * \brief Write linear part of a function.
   * \param [in] i The index of constraint whose linear part we want to write.
   * It is ignored if is_obj is true.
   * \param [is_obj] True if objective, false otherwise.
   * \param [in] stream to which to write.
   */
  void writeLin_(Minotaur::UInt i, bool is_obj, std::ostream &out) const;
};

} // namespace MINOTAUR_AMPL

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
