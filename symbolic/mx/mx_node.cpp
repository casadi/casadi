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

#include "mx_node.hpp"
#include "mx_tools.hpp"
#include "../stl_vector_tools.hpp"
#include <cassert>
#include <typeinfo> 
#include "../matrix/matrix_tools.hpp"
#include "transpose.hpp"
#include "reshape.hpp"
#include "multiplication.hpp"
#include "subref.hpp"
#include "subassign.hpp"
#include "getnonzeros.hpp"
#include "setnonzeros.hpp"
#include "densification.hpp"
#include "solve.hpp"
#include "unary_mx.hpp"
#include "binary_mx.hpp"


// Template implementations
#include "setnonzeros_impl.hpp"
#include "multiplication_impl.hpp"
#include "solve_impl.hpp"

using namespace std;

namespace CasADi{

  MXNode::MXNode(){
    temp = 0;
  }

  MXNode::~MXNode(){

    // Start destruction method if any of the dependencies has dependencies
    for(vector<MX>::iterator cc=dep_.begin(); cc!=dep_.end(); ++cc){
      // Skip if null
      if(cc->isNull()) continue;
    
      // Check if there are other "owners" of the node
      if(cc->getCount()!= 1){
      
	// Replace with a null pointer
	*cc = MX();
      
      } else {
	// Stack of experssions to be deleted
	std::stack<MX> deletion_stack;
        
	// Move the child to the deletion stack
	deletion_stack.push(*cc);
	*cc = MX();
      
	// Process stack
	while(!deletion_stack.empty()){
        
	  // Top element
	  MX t = deletion_stack.top();
          
	  // Check if the top element has dependencies with dependencies
	  bool found_dep = false;
        
	  // Start destruction method if any of the dependencies has dependencies
	  for(vector<MX>::iterator ii=t->dep_.begin(); ii!=t->dep_.end(); ++ii){
          
	    // Skip if null
	    if(ii->isNull()) continue;
          
	    // Check if this is the only reference to the element
	    if(ii->getCount()==1){
            
	      // Remove and add to stack
	      deletion_stack.push(*ii);
	      *ii = MX();
	      found_dep = true;
	      break;
	    } else {
	      // Replace with an element without dependencies
	      *ii = MX();
	    }
	  }
        
	  // Pop from stack if no dependencies found
	  if(!found_dep){
	    deletion_stack.pop();
	  }
	}
      }
    }
  }

  const string& MXNode::getName() const{
    throw CasadiException(string("MXNode::getName() not defined for class ") + typeid(*this).name());
  }

  bool MXNode::__nonzero__() const {
    casadi_error("Can only determine truth value of a numeric MX.");

  }

  const MX& MXNode::dep(int ind) const{
    return dep_.at(ind);
  }
  
  MX& MXNode::dep(int ind){
    return dep_.at(ind);
  }
  
  int MXNode::ndep() const{
    return dep_.size();
  }

  void MXNode::setSparsity(const CRSSparsity& sparsity){
    sparsity_ = sparsity;
  }

  void MXNode::setDependencies(const MX& dep){
    dep_.resize(1);
    dep_[0] = dep;
  }
    
  void MXNode::setDependencies(const MX& dep1, const MX& dep2){
    dep_.resize(2);
    dep_[0] = dep1;
    dep_[1] = dep2;
  }
    
  void MXNode::setDependencies(const MX& dep1, const MX& dep2, const MX& dep3){
    dep_.resize(3);
    dep_[0] = dep1;
    dep_[1] = dep2;
    dep_[2] = dep3;
  }

  int MXNode::addDependency(const MX& dep){
    dep_.push_back(dep);
    return dep_.size()-1;
  }

  void MXNode::assign(const MX& d, const std::vector<int>& inz, bool add){
    casadi_assert(0);
  }

  void MXNode::assign(const MX& d, const std::vector<int>& inz, const std::vector<int>& onz, bool add){
    casadi_assert(0);
  }

  void MXNode::setDependencies(const std::vector<MX>& dep){
    dep_ = dep;
  }

  int MXNode::numel() const{
    return sparsity_.numel();
  }

  int MXNode::size() const{
    return sparsity_.size();
  }

  int MXNode::size1() const{
    return sparsity_.size1();
  }

  int MXNode::size2() const{
    return sparsity_.size2();
  }

  const CRSSparsity& MXNode::sparsity() const{
    return sparsity_;
  }

  const CRSSparsity& MXNode::sparsity(int oind) const{
    casadi_assert_message(oind==0, "Index out of bounds");
    return sparsity_;
  }

  void MXNode::repr(std::ostream &stream) const{
    stream << "MX(";
    print(stream);
    stream << ")";
  }

  void MXNode::print(std::ostream &stream) const{
    long remaining_calls = MX::getMaxNumCallsInPrint();
    print(stream,remaining_calls);
  }

  void MXNode::print(std::ostream &stream, long& remaining_calls) const{
    if(remaining_calls>0){
      remaining_calls--;
      printPart(stream,0);
      for(int i=0; i<ndep(); ++i){
	if (dep(i).isNull()) {
	  stream << "MX()";
	} else {
	  dep(i)->print(stream,remaining_calls);
	}
	printPart(stream,i+1);
      }
    } else {
      stream << "...";
    }
  }

  void MXNode::printPart(std::ostream &stream, int part) const{
    casadi_assert(ndep()>1);
    casadi_assert(part>0);
    casadi_assert(part<ndep());
    stream << ",";
  }

  FX& MXNode::getFunction(){
    throw CasadiException(string("MXNode::getFunction() not defined for class ") + typeid(*this).name());
  }

  int MXNode::getFunctionOutput() const{
    throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") + typeid(*this).name());
  }

  int MXNode::getFunctionInput() const{
    throw CasadiException(string("MXNode::getFunctionOutput() not defined for class ") + typeid(*this).name());
  }

  void MXNode::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp){
    DMatrixPtrVV fwdSeed, fwdSens, adjSeed, adjSens;
    evaluateD(input,output,fwdSeed, fwdSens, adjSeed, adjSens, itmp, rtmp);
  }
  
  void MXNode::evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, 
			 const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, 
			 const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens){
    throw CasadiException(string("MXNode::evaluateD not defined for class ") + typeid(*this).name());
  }

  void MXNode::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, std::vector<int>& itmp, std::vector<SX>& rtmp){
    SXMatrixPtrVV fwdSeed, fwdSens, adjSeed, adjSens;
    evaluateSX(input,output,fwdSeed, fwdSens, adjSeed, adjSens, itmp, rtmp);
  }
  
  void MXNode::evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, 
			  const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, 
			  const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens){
    throw CasadiException(string("MXNode::evaluateSX not defined for class ") + typeid(*this).name());
  }

  void MXNode::evaluateMX(const MXPtrV& input, MXPtrV& output){
    MXPtrVV fwdSeed, fwdSens, adjSeed, adjSens;
    evaluateMX(input,output,fwdSeed, fwdSens, adjSeed, adjSens,false);
  }

  void MXNode::propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd){
    throw CasadiException(string("MXNode::propagateSparsity not defined for class ") + typeid(*this).name());
  }

  void MXNode::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    SharedObjectNode::deepCopyMembers(already_copied);
    dep_ = deepcopy(dep_,already_copied);
  }

  MX MXNode::getOutput(int oind) const{
    casadi_assert_message(oind==0,"Output index out of bounds");
    return shared_from_this<MX>();
  }

  void MXNode::generateOperation(std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const{
    stream << "#error " <<  typeid(*this).name() << ": " << arg << " => " << res << endl;
  }

  double MXNode::getValue() const{
    throw CasadiException(string("MXNode::getValue not defined for class ") + typeid(*this).name());    
  }

  Matrix<double> MXNode::getMatrixValue() const{
    throw CasadiException(string("MXNode::getMatrixValue not defined for class ") + typeid(*this).name());    
  }

  MX MXNode::getTranspose() const{
    if(sparsity().dense()){
      return MX::create(new DenseTranspose(shared_from_this<MX>()));
    } else {
      return MX::create(new Transpose(shared_from_this<MX>()));
    }
  }

  MX MXNode::getReshape(const CRSSparsity& sp) const{
    return MX::create(new Reshape(shared_from_this<MX>(),sp));
  }
  
  MX MXNode::getMultiplication(const MX& y) const{
    // Transpose the second argument
    MX trans_y = trans(y);
    CRSSparsity sp_z = sparsity().patternProduct(trans_y.sparsity());
    MX z = MX::zeros(sp_z);
    if(sparsity().dense() && y.dense()){
      return MX::create(new DenseMultiplication<false,true>(z,shared_from_this<MX>(),trans_y));
    } else {
      return MX::create(new Multiplication<false,true>(z,shared_from_this<MX>(),trans_y));    
    }
  }
    
  MX MXNode::getSolve(const MX& r, bool tr) const{
    if(tr){
      return MX::create(new Solve<true>(r,shared_from_this<MX>()));
    } else {
      return MX::create(new Solve<false>(r,shared_from_this<MX>()));
    }
  }

  MX MXNode::getGetNonzeros(const CRSSparsity& sp, const std::vector<int>& nz) const{
    if(nz.size()==0){
      return MX::zeros(sp);
    } else {
      MX ret;
      if(Slice::isSlice(nz)){
	ret = MX::create(new GetNonzerosSlice(sp,shared_from_this<MX>(),Slice(nz)));
      } else if(Slice::isSlice2(nz)){
	Slice outer;
	Slice inner(nz,outer);
	ret = MX::create(new GetNonzerosSlice2(sp,shared_from_this<MX>(),inner,outer));
      } else {
	ret = MX::create(new GetNonzerosVector(sp,shared_from_this<MX>(),nz));
      }
      simplify(ret);
      return ret;
    }
  }

  MX MXNode::getSetNonzeros(const MX& y, const std::vector<int>& nz) const{
    if(nz.size()==0){
      return y;
    } else {
      MX ret;
      if(Slice::isSlice(nz)){
	ret = MX::create(new SetNonzerosSlice<false>(y,shared_from_this<MX>(),Slice(nz)));
      } else if(Slice::isSlice2(nz)){
	Slice outer;
	Slice inner(nz,outer);
	ret = MX::create(new SetNonzerosSlice2<false>(y,shared_from_this<MX>(),inner,outer));
      } else {
	ret = MX::create(new SetNonzerosVector<false>(y,shared_from_this<MX>(),nz));
      }
      simplify(ret);
      return ret;
    }
  }


  MX MXNode::getAddNonzeros(const MX& y, const std::vector<int>& nz) const{
    if(nz.size()==0){
      return y;
    } else {
      MX ret;
      if(Slice::isSlice(nz)){
	ret = MX::create(new SetNonzerosSlice<true>(y,shared_from_this<MX>(),Slice(nz)));
      } else if(Slice::isSlice2(nz)){
	Slice outer;
	Slice inner(nz,outer);
	ret = MX::create(new SetNonzerosSlice2<true>(y,shared_from_this<MX>(),inner,outer));
      } else {
	ret = MX::create(new SetNonzerosVector<true>(y,shared_from_this<MX>(),nz));
      }
      simplify(ret);
      return ret;
    }
  }

  MX MXNode::getDensification(const CRSSparsity& sp) const{
    return MX::create(new Densification(shared_from_this<MX>(),sp));
  }

  
  MX MXNode::getSubRef(const Slice& i, const Slice& j) const{
    return MX::create(new SubRef(shared_from_this<MX>(),i,j));
  }

  MX MXNode::getSubAssign(const MX& y, const Slice& i, const Slice& j) const{
    return MX::create(new SubAssign(shared_from_this<MX>(),y,i,j));
  }

  MX MXNode::getUnary(int op) const{
    if(operation_checker<F0XChecker>(op) && isZero()){
      // If identically zero
      return MX::sparse(size1(),size2());
    } else {
      // Create a new node
      return MX::create(new UnaryMX(Operation(op),shared_from_this<MX>()));
    }
  }

  MX MXNode::getBinary(int op, const MX& y) const{
    // Make sure that dimensions match
    casadi_assert_message((sparsity().scalar() || y.scalar() || (sparsity().size1()==y.size1() && size2()==y.size2())),"Dimension mismatch." << "lhs is " << sparsity().dimString() << ", while rhs is " << y.dimString());
  
    // Quick return if zero
    if((operation_checker<F0XChecker>(op) && isZero()) || 
       (operation_checker<FX0Checker>(op) && y->isZero())){
      return MX::sparse(std::max(size1(),y.size1()),std::max(size2(),y.size2()));
    }
  
    // Create binary node
    if(sparsity().scalar())
      return getScalarMatrix(op,y);
    else if(y.scalar())  
      return getMatrixScalar(op,y);
    else
      return getMatrixMatrix(op,y);
  }

  MX MXNode::getMatrixMatrix(int op, const MX& y) const{
    // Check if we can carry out the operation only on the nonzeros
    if((sparsity().dense() && y.dense()) ||
       (operation_checker<F00Checker>(op) && sparsity()==y.sparsity())){
      // Loop over nonzeros only
      return MX::create(new NonzerosNonzerosOp(Operation(op),shared_from_this<MX>(),y)); 
    } else {
      // Sparse matrix-matrix operation necessary
      return MX::create(new SparseSparseOp(Operation(op),shared_from_this<MX>(),y)); 
    }
  }

  MX MXNode::getScalarMatrix(int op, const MX& y) const{
    // If the scalar is sparse, then replace with dense
    if(size()==0){
      return MX(0)->getScalarMatrix(op,y);
    } else {
      // Check if it is ok to loop over nonzeros only
      if(y.dense() || operation_checker<FX0Checker>(op)){
	// Loop over nonzeros
	return MX::create(new ScalarNonzerosOp(Operation(op),shared_from_this<MX>(),y));
      } else {
	// Put a densification node in between
	return getScalarMatrix(op,densify(y));
      }
    }
  }


  MX MXNode::getMatrixScalar(int op, const MX& y) const{
    // Check if the scalar is sparse, then replace with dense
    if(y.size()==0){
      return getMatrixScalar(op,0);
    } else {
      // Check if it is ok to loop over nonzeros only
      if(sparsity().dense() || operation_checker<F0XChecker>(op)){
	// Loop over nonzeros
	return MX::create(new NonzerosScalarOp(Operation(op),shared_from_this<MX>(),y));
      } else {
	// Put a densification node in between
	return densify(shared_from_this<MX>())->getMatrixScalar(op,y);
      }
    }
  }

  Matrix<int> MXNode::mapping() const{
    throw CasadiException(string("MXNode::mapping not defined for class ") + typeid(*this).name());
  }

} // namespace CasADi

