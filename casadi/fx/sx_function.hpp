/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef SX_FUNCTION_HPP
#define SX_FUNCTION_HPP

#include "fx.hpp"
#include "../mx/mx_node.hpp"
#include "../sx/sx_matrix.hpp"
#include <map>
#include <vector>
#include <set>
#include <stack>

namespace CasADi{

/** \brief  Forward declaration of internal class */
class SXFunctionNode;

/** \brief  Function mapping from and to SX or SXMatrix 
  \author Joel Andersson 
  \date 2010
*/
class SXFunction : public FX{

public:

/** \brief  CONSTRUCTORS: */
/** \brief  default constructor */
  SXFunction();

  /** \brief  Single (scalar/matrix/vector valued) input, single (scalar/matrix/vector valued) output */
  SXFunction(const SXMatrix& arg, const SXMatrix& res);

  /** \brief  Multiple (vector valued) input, single (scalar/vector/matrix valued) output */
  SXFunction(const vector< vector<SX> >& arg, const SXMatrix& res);

  /** \brief  Multiple (matrix valued) input, single (scalar/vector/matrix valued) output */
  SXFunction(const vector< SXMatrix>& arg, const SXMatrix& res);

  /** \brief  Multiple (vector valued) input, multiple (vector valued) output */
  SXFunction(const vector< vector<SX> >& arg, const vector< vector<SX> >& res);

  /** \brief  Multiple (matrix valued) input, multiple (matrix valued) output */
  SXFunction(const vector< SXMatrix>& arg, const vector<SXMatrix>& res);
  
/** \brief  Access functions of the node */
  SXFunctionNode* operator->();
  const SXFunctionNode* operator->() const;
  
  /** \brief  evaluate symbolically */
  std::vector<SXMatrix> eval(const std::vector<SXMatrix>& arg);

  /** \brief  evaluate symbolically, single input, single output */
  SXMatrix eval(const SXMatrix& arg);

  /** \brief  evaluate symbolically (pass and get non-zero entries) */
  std::vector< std::vector<SX> > eval(const std::vector< std::vector<SX> >& arg);

  /** \brief  evaluate symbolically, single input, single output (pass and get non-zero entries) */
  std::vector<SX> eval(const std::vector<SX>& arg);
  
  /** \brief Jacobian of output oind with respect to input iind */
  SXFunction jacobian(int iind=0, int oind=0);
  
  /** \brief Hessian of output oind with respect to input iind */
  SXFunction hessian(int iind=0, int oind=0);

  /// Jacobian via source code transformation
  SXMatrix jac(int iind=0, int oind=0);

  /// Gradient via source code transformation
  SXMatrix grad(int iind=0, int oind=0);
  
  /// Hessian (forward over adjoint) via source code transformation
  SXMatrix hess(int iind=0, int oind=0);

  /** \brief get the an input argument symbolically */
  SXMatrix getArgumentIn(int iind=0) const;

  /** \brief get the an output argument symbolically */
  SXMatrix getArgumentOut(int iind=0) const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

};
  
template<int n>
struct AlgElData{
    // Partial derivatives
    double d[n+1];
};

/** \brief  Internal node class for SXFunction
  \author Joel Andersson 
  \date 2010
 A regular user should never work with any Node class. Use SXFunction directly.
*/
class SXFunctionNode : public FXNode{
  friend class SXFunction;
  
  protected:
    /** \brief  Constructor (only to be called from SXFunction, therefore protected) */
    SXFunctionNode(const std::vector<SXMatrix>& inputv, const std::vector<SXMatrix>& outputv);

  public:

  /** \brief  Make a deep copy */
  virtual SXFunctionNode* clone() const;
    
/** \brief  Destructor */
  virtual ~SXFunctionNode();

/** \brief  Clear the memory */
  virtual void clear(int ord=0);

/** \brief  Evaluate the function with partial derivatives up to order ord */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  evaluate symbolically */
  void eval(const std::vector<SXMatrix>& input_s, std::vector<SXMatrix>& output_s);
  void eval(const SXMatrix &x, SXMatrix &res) const; 

/** \brief  evaluate symbolically, replacing nodes */
  void eval(const SXMatrix &x, SXMatrix &res, const std::map<int,SX>& replace, SXMatrix &repres) const;

/** \brief  Check if smooth */
  bool isSmooth() const;

/** \brief  Print the algorithm */
  void printAlgorithm(std::ostream &stream=std::cout) const;
//   void printValues(std::ostream &stream=std::cout);

  /** \brief Jacobian of output oind with respect to input iind */
  virtual FX jacobian(int iind=0, int oind=0);
  
  /** \brief Hessian of output oind with respect to input iind */
  virtual FX hessian(int iind=0, int oind=0);
  
  /** \brief  Print */
  virtual void print(std::ostream &stream) const;

  /// Jacobian via source code transformation
  SXMatrix jac(int iind=0, int oind=0);

  /// Gradient via source code transformation
  SXMatrix grad(int iind=0, int oind=0);
  
  /// Hessian (forward over adjoint) via source code transformation
  SXMatrix hess(int iind=0, int oind=0);

/** \brief  DATA MEMBERS */
  
/** \brief  Indices of the nodes corresponding to the inputs */
  std::vector<std::vector<int> > input_ind;
  
/** \brief  Indices of the nodes corresponding the non-zeros of the outputs */
  std::vector<std::vector<int> > output_ind;

/** \brief  An elemenent of the algorithm, namely a binary operation */
  struct AlgEl{
    unsigned short op; // operator
    int ind; // index of the binary operaton to be evaluated
    int ch[2]; // indices of the arguments
  };
    
/** \brief  all binary nodes of the tree in the order of execution */
  std::vector<AlgEl> algorithm;
  std::vector<AlgElData<1> > pder1;
  std::vector<AlgElData<2> > pder2;
  
/** \brief  All nodes */
  std::vector<SXNode*> tree;

  /** \brief  Working vector for numeric calculation */
  std::vector< std::vector<double> > work;        // work array during the evaluation
  int worksize;

  /// work vector for symbolic calculations (allocated first time)
  std::vector<SX> work_sym;
  
/** \brief  Initialize */
  virtual void init();

  /** Maximal order of the automatic differentiation*/
  int maxorder;

/** \brief  Print to a c file */
  void generateCode(const std::string& filename) const;
  
/** \brief  Topological sorting of the nodes based on Depth-First Search (DFS) */
  static void sort_depth_first(std::stack<SXNode*>& s, std::vector<BinarySXNode*>& algnodes);

/** \brief  Topological (re)sorting of the nodes based on Bredth-First Search (BFS) (Kahn 1962) */
  static void resort_bredth_first(std::vector<BinarySXNode*>& algnodes);

/** \brief  Topological (re)sorting of the nodes with the purpose of postponing every calculation as much as possible, as long as it does not influence a dependent node */
  static void resort_postpone(std::vector<BinarySXNode*>& algnodes, std::vector<int>& lind);
  
/** \brief  Inputs of the function (needed for symbolic calculations) */
  std::vector<SXMatrix> inputv;

/** \brief  Outputs of the function (needed for symbolic calculations) */
  std::vector<SXMatrix> outputv;


};


} // namespace CasADi

#endif // SX_FUNCTION_HPP
