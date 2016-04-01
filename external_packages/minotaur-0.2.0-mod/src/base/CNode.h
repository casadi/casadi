//
//    MINOTAUR -- It's only 1/2 bull
//
//    (C)opyright 2008 - 2014 The MINOTAUR Team.
//


/**
 * \file CNode.h
 * \brief Declare class CNode to represent node of a computational graph of a
 * nonlinear function.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURCNODE_H
#define MINOTAURCNODE_H

#include "OpCode.h"
#include "Types.h"

namespace Minotaur {

class Variable;
class CNode;

struct CompareCNodes {
  bool operator()(const CNode* n1, const CNode *n2) const;
};
struct CompareCNodesR {
  bool operator()(const CNode* n1, const CNode *n2) const;
};
typedef boost::shared_ptr<Variable> VariablePtr;
typedef std::set<CNode *, CompareCNodes> CNodeSet;
typedef std::set<CNode *, CompareCNodesR> CNodeRSet;

class CQIter2 {
public:
  CNode *node;
  CQIter2 *next;
  CQIter2 *prev;
};

class CNode {
public:

  /// Default constructor
  CNode();

  /**
   * Constructor with a specific opcode and children. Children can be
   * NULL, and num_child can be zero.
   */
  CNode(OpCode op, CNode **children, UInt num_child);
  CNode(OpCode op, CNode *lchild, CNode *rchild);
  ~CNode();


  CNode *clone() const;
  void addPar(CNode *node);

  /// Replace one of the children of a given node by another node.
  void changeChild(CNode *out, CNode *in);

  void copyParChild(CNode *out, std::map<const CNode*, CNode*> *nmap) const;
  void eval(const double *x, int *error);
  double eval(double x, int *error) const;
  FunctionType findFType();
  void fwdGrad();
  bool getB() const { return b_; };
  double getG() const { return g_; };
  double getH() const {return h_;};
  UInt getId() const {return id_;};
  CNode* getL() const { return l_; };
  double getLb() const { return lb_; };
  CNode** getListL() const { return child_; };
  CNode** getListR() const { return child_+numChild_; };
  OpCode getOp() const { return op_; };
  CQIter2* getParB() const { return parB_; };
  CNode* getR() const { return r_; };
  int getTempI() const { return ti_; };
  FunctionType getType() const { return fType_; } ;
  double getUb() const { return ub_; };
  CNode* getUPar() const { return uPar_; };
  const Variable* getV() const { return v_; };
  double getVal() const { return val_; };
  void grad(int *error);
  UInt numChild() const;
  UInt numPar() const { return numPar_; };
  void propHessSpa();
  void propHessSpa2(CNodeRSet *nset);
  void propBounds(bool *is_inf, int *error);
  void hess(int *error);
  void hess2(CNodeRSet *nset, int *error);
  void setB(bool b) {b_ = b;};
  void setBounds(double lb, double ub) {lb_ = lb; ub_ = ub; };
  void setDouble(double d) {d_ = d;};
  void setG(double g) {g_ = g;};
  void setGi(double gi) {gi_ = gi;};
  void setH(double h) {h_ = h;};
  void setId(UInt i) {id_ = i;};
  void setInt(int i) {i_ = i;};
  void setL(CNode *n) {l_ = n;};
  void setOp(OpCode op) {op_ = op;};
  void setR(CNode *n) {r_ = n;};
  void setTempI(int i) { ti_ = i; };
  void setType(FunctionType t);
  void setVal(double v);
  void setV(VariablePtr v) {v_ = v.get();};
  void updateBnd(int *error);

  void write(std::ostream &out) const;
  void writeSubExp(std::ostream &out) const;

protected:
  bool b_;
  CNode **child_; // array of size numChild_ + 1
  double d_;
  FunctionType fType_;
  double g_;
  double gi_;
  double h_;
  int i_;
  UInt id_;
  CNode *l_;
  double lb_;
  UInt numChild_;
  UInt numPar_;
  OpCode op_;
  CQIter2 *parB_;
  CQIter2 *parE_;
  CNode *r_;
  int ti_;
  double ub_;

  /// Unique parent of this node. NULL if the node has multiple parents
  CNode *uPar_;
  const Variable* v_;
  double val_;
  void propBounds_(double lb, double ub, bool *is_inf);

};
typedef std::vector<CNode*> CNodeVector;
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
