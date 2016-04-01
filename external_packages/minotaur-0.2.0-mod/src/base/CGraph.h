//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file CGraph.h
 * \brief Declare class CGraph for storing computational graph of a nonlinear
 * function.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURCGRAPH_H
#define MINOTAURCGRAPH_H

#include <stack>

#include "Types.h"
#include "NonlinearFunction.h"
#include "OpCode.h"

namespace Minotaur {

class CGraph;
class CNode;
typedef boost::shared_ptr<CGraph> CGraphPtr;
typedef std::deque<CNode *> CNodeQ;
typedef std::vector<CNode *> CNodeVector;
typedef std::map<ConstVariablePtr, CNode*, CompareVariablePtr> VarNodeMap;

class CGraph : public NonlinearFunction {
public:
  /// Default constructor.
  CGraph();

  /// Default constructor.
  ~CGraph();

  void addConst(const double eps, int &err); 

  NonlinearFunctionPtr clone(int *err) const;

  NonlinearFunctionPtr cloneWithVars(VariableConstIterator, int *) const;

  // base class method.
  NonlinearFunctionPtr getPersp(VariablePtr z, double eps, int *err) const;

  // base class method.
  void computeBounds(double *lb, double *ub, int *error);

  // Evaluate at a given array.
  double eval(const double *x, int *err);

  // Evaluate gradient at a given array.
  void evalGradient(const double *x, double *grad_f, int *error);

  // Evaluate hessian of at a given vector.
  void evalHessian(double mult, const double *x, 
                   const LTHessStor *stor, double *values, 
                   int *error);

  // Fill hessian sparsity.
  void fillHessStor(LTHessStor *stor);

  // Finalize hessian offsets, if needed.
  void finalHessStor(const LTHessStor *stor);

  // Add gradient values to sparse Jacobian
  void fillJac(const double *x, double *values, int *error);

  /**
   * After adding all the nodes of the graph, finalize is called to create the
   * forward and backward traversal queues, and related book-keeping.
   */
  void finalize(); 

  // base class method.
  double getFixVarOffset(VariablePtr v, double val);

  // get type of function
  FunctionType getType() const;

  UInt getNumNodes();

  UInt getHessNz();

  /// Get the output node.
  const CNode* getOut() const;

  VariablePtr getVar(const CNode *cnode) const;

  // Get variables that exist in this function.
  void getVars(VariableSet *vars);

  bool isIdenticalTo(CGraphPtr cg);

  // base class method
  bool isSumOfSquares() const;

  // multiply by a constant.
  void multiply(double c);

  /**
   * \brief Create a new node with one or two children, and add it to the
   * graph. The children should already be a nodes of the graph.
   *
   * \param [in] op OpCode of the new node
   * \param [in] lchild The left child of the new node. It should
   * never be NULL.
   * \param [in] rchild If the opcode needs two children, e.g. OpPlus, OpMult
   * etc, then rchild should be the right child. For univariate functions like
   * OpSqrt, it should be NULL.
   * \return The new node created in this function.
   */
  CNode* newNode(OpCode op, CNode *lchild, CNode *rchild);

  /**
   * \brief Create a new node with more than two children, and add it to the
   * graph. The children should already be a nodes of the graph.
   *
   * \param [in] op OpCode of the new node
   * \param [in] child An array of pointers to children. It should be size 'n'.
   * \param [in] n The size of array 'child'.
   * \return The new node created in this function.
   */
  CNode* newNode(OpCode op, CNode **child, UInt n);

  /**
   * \brief Create a new node with constant real value. It does not have any
   * children.
   *
   * \param [in] d The value.
   * \return The new node created in this function. Its opcode is OpNum
   */
  CNode* newNode(double d);

  /**
   * \brief Create a new node with constant integer value. It does not have any
   * children.
   *
   * \param [in] i The value.
   * \return The new node created in this function. Its opcode is OpInt
   */
  CNode* newNode(int i);

  /**
   * \brief Create a new node denoting an input variable. It does not have any
   * children. 
   *
   * The function checks if the variable already exists in the graph. If so,
   * it returns the pointer to that node. Otherwise it creates a new node.
   *
   * \param [in] v The variable.
   * \return The new node found or created in this function. Its opcode is
   * OpVar
   */
  CNode* newNode(VariablePtr v);

  // base class method.
  void prepJac(VarSetConstIter vbeg, VarSetConstIter vend);

  // base class method.
  void removeVar(VariablePtr v, double val);

  /**
   * \brief Set the node that should be the output of this graph. This node
   * should already be a part of this graph (created by newNode() function).
   *
   * \param [in] node The node that is the output.
   */
  void setOut(CNode *node);

  // base class method.
  void sqrRoot(int &err);

  // base class method.
  void subst(VariablePtr out, VariablePtr in, double rat);

  // base class method.
  void varBoundMods(double lb, double ub, VarBoundModVector &mods,
                    SolveStatus *status);

  // display.
  void write(std::ostream &out) const;

private:
  /// All nodes of the graph.
  CNodeVector aNodes_; 

  bool changed_;

  /// All dependent nodes, i.e. nodes with OpCode different from OpVar, OpInt
  /// and OpNum.
  CNodeQ dq_;

  UIntVector hInds_;
  UInt hNnz_;
  UIntVector hOffs_;
  UIntVector hStarts_;
  UIntVector gOffs_;


  /// Topmost node or output node. We assume only one is present.
  CNode *oNode_;

  /// A map that tells which node corresponds to a given variable.
  VarNodeMap varNode_;

  /// All nodes with OpCode OpVar.
  CNodeQ vq_;

  CGraphPtr clone_(int *err) const;

  void fwdGrad_(CNode *node);
  void fwdGrad2_(std::stack<CNode *> *st2, CNode *node);

  void fillHessInds_(CNode *node, UIntQ *inds);
  void fillHessInds2_(CNode *node, UIntQ *inds);

  /// Recursive function to check whether CGraph represents a sum of squares.
  bool isSOSRec_(CNode *node) const;

  void revHess_(int *error);
  void revHess2_(std::stack<CNode *> *st2, double mult, UInt vind,
                 double *values, UInt *nz, int *error);

  /**
   *  Routine to propagate gradient by a reverse mode traversal.
   *
   *  \param [out] error Set to a nonzero if an error occurs. Otherwise, leave
   *  it unchanged.
   */
  void grad_(int *error);

  void setupHess_(VariablePtr v, CNode *node, std::set<ConstVariablePair, 
                  CompareVariablePair> & vps);

  void simplifyDq_();

};
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
