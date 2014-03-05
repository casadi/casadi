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

#include "fx_internal.hpp"
#include "../mx/call_fx.hpp"
#include "../fx/mx_function.hpp"
#include <typeinfo> 
#include "../stl_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "parallelizer.hpp"

using namespace std;

namespace CasADi{
  
  FX::FX(){
  }

  FX::~FX(){
  }

  FX FX::create(FXInternal* node){
    FX ret;
    ret.assignNode(node);
    return ret;
  }

  const FXInternal* FX::operator->() const{
    return static_cast<const FXInternal*>(OptionsFunctionality::operator->());
  }

  FXInternal* FX::operator->(){
    return static_cast<FXInternal*>(OptionsFunctionality::operator->());
  }

  vector<MX> FX::call(const MX &arg){
    vector<MX> xflatten(1,arg);
    return call(xflatten);
  }

  vector<MX> FX::call(const vector<MX> &arg){
    MXVectorVector dummy;
    MXVector res;
    call(arg,res,dummy,dummy,dummy,dummy);
    return res;
  }

  void FX::call(const MXVector& arg, MXVector& res,  const MXVectorVector& fseed, MXVectorVector& fsens, 
                const MXVectorVector& aseed, MXVectorVector& asens){
    casadi_assert_message(arg.size()==getNumInputs(),"FX::call: dimension mismatch. You supplied " << arg.size() << " arguments instead of suspected " << getNumInputs() << ".");
    (*this)->call(arg,res,fseed,fsens,aseed,asens,false,true);
  }

  vector<vector<MX> > FX::call(const vector<vector<MX> > &x, const Dictionary& paropt){
    assertInit();
  
    // Make sure not empty
    casadi_assert_message(x.size()>1,"FX: call(vector<vector<MX> >): argument must be of length > 1. You supplied length " << x.size() << ".");
  
    // Return object
    vector<vector<MX> > ret(x.size());
    
    // Check if we are bypassing the parallelizer
    Dictionary::const_iterator ii=paropt.find("parallelization");
    if(ii!=paropt.end() && ii->second=="expand"){
      for(int i=0; i<x.size(); ++i){
        ret[i] = call(x[i]);
      }
      return ret;
    }
      
    // Create parallelizer object and initialize it
    Parallelizer p(vector<FX>(x.size(),*this));
    p.setOption(paropt);
    p.init();
  
    // Concatenate the arguments
    vector<MX> p_in;
    p_in.reserve(x.size() * getNumInputs());
    for(int i=0; i<x.size(); ++i){
      p_in.insert(p_in.end(),x[i].begin(),x[i].end());
      p_in.resize(p_in.size()+getNumInputs()-x[i].size());
    }
  
    // Call the parallelizer
    vector<MX> p_out = p.call(p_in);
    casadi_assert(p_out.size() == x.size() * getNumOutputs());

    // Collect the outputs
    vector<MX>::const_iterator it=p_out.begin();
    for(int i=0; i<x.size(); ++i){
      ret[i].insert(ret[i].end(),it,it+getNumOutputs());
      it += getNumOutputs();
    }
    return ret;
  }

  void FX::evaluate(){
    assertInit();
    (*this)->evaluate();
  }

  void FX::solve(){
    evaluate();
  }

  int FX::getNumInputNonzeros() const{
    return (*this)->getNumInputNonzeros();
  }

  int FX::getNumOutputNonzeros() const{
    return (*this)->getNumOutputNonzeros();
  }

  int FX::getNumInputElements() const{
    return (*this)->getNumInputElements();
  }

  int FX::getNumOutputElements() const{
    return (*this)->getNumOutputElements();
  }

  FX FX::jacobian(int iind, int oind, bool compact, bool symmetric){
    assertInit();
    return (*this)->jacobian(iind,oind,compact,symmetric);
  }

  void FX::setJacobian(const FX& jac, int iind, int oind, bool compact){
    (*this)->setJacobian(jac,iind,oind,compact);
  }

  FX FX::gradient(int iind, int oind){
    assertInit();
    return (*this)->gradient(iind,oind);
  }

  FX FX::tangent(int iind, int oind){
    assertInit();
    return (*this)->tangent(iind,oind);
  }

  FX FX::hessian(int iind, int oind){
    assertInit();
    return (*this)->hessian(iind,oind);  
  }

  FX FX::fullJacobian(){
    assertInit();
    return (*this)->fullJacobian();
  }

  void FX::setFullJacobian(const FX& jac){
    (*this)->full_jacobian_ = jac;
  }

  bool FX::checkNode() const{
    return dynamic_cast<const FXInternal*>(get())!=0;
  }

  void FX::addMonitor(const string& mon){
    (*this)->monitors_.insert(mon);
  }

  void FX::removeMonitor(const string& mon){
    (*this)->monitors_.erase(mon);
  }

  const Dictionary & FX::getStats() const{
    return (*this)->getStats();
  }

  GenericType FX::getStat(const string& name) const{
    return (*this)->getStat(name);
  }

  Sparsity& FX::jacSparsity(int iind, int oind, bool compact, bool symmetric){
    return (*this)->jacSparsity(iind,oind,compact, symmetric);
  }

  void FX::setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact){
    (*this)->setJacSparsity(sp,iind,oind,compact);
  }

  std::vector<MX> FX::symbolicInput() const{
    return (*this)->symbolicInput();
  }

  std::vector<SX> FX::symbolicInputSX() const{
    return (*this)->symbolicInputSX();
  }

  void FX::setInputScheme(const IOScheme &scheme) {
    return (*this)->setInputScheme(scheme);
  }

  void FX::setOutputScheme(const IOScheme &scheme) {
    return (*this)->setOutputScheme(scheme);
  }

  IOScheme FX::getInputScheme() const {
    return (*this)->getInputScheme();
  }

  IOScheme FX::getOutputScheme() const {
    return (*this)->getOutputScheme();
  }

  FX FX::operator[](int k) const {

    // Argument checking
    if (k<0) k+=getNumOutputs();
    casadi_assert_message(k<getNumOutputs(),"FX[int k]:: Attempt to select the k'th output with k=" << k << ", but should be smaller or equal to number of outputs (" << getNumOutputs() << ").");
  
    // Get the inputs in MX form
    std::vector< MX > in = symbolicInput();
  
    // Clone such that we can use const.
    FX clone = *this;
  
    // Get the outputs
    std::vector< MX > result = clone.call(in);
  
    // Construct an MXFunction with only the k'th output
    MXFunction ret(in,result[k]);
  
    ret.setInputScheme(getInputScheme());
  
    // Initialize it
    ret.init();
  
    // And return it, will automatically cast to FX
    return ret;
  }

  vector<SX> FX::evalSX(const vector<SX>& arg){
    casadi_assert_message(arg.size()==getNumInputs(),"FX::evalSX: dimension mismatch. You supplied " << arg.size() << " arguments instead of suspected " << getNumInputs() << ".");
    vector<SX> res;
    vector<vector<SX> > dummy;
    (*this)->evalSX(arg,res,dummy,dummy,dummy,dummy);
    return res;
  }

  vector<MX> FX::evalMX(const vector<MX>& arg){
    casadi_assert_message(arg.size()==getNumInputs(),"FX::evalMX: dimension mismatch. You supplied " << arg.size() << " arguments instead of suspected " << getNumInputs() << ".");
    vector<MX> res;
    vector<vector<MX> > dummy;
    (*this)->evalMX(arg,res,dummy,dummy,dummy,dummy);
    return res;
  }

  void FX::evalSX(const std::vector<SX>& arg, std::vector<SX>& res, 
                  const std::vector<std::vector<SX> >& fseed, std::vector<std::vector<SX> >& fsens, 
                  const std::vector<std::vector<SX> >& aseed, std::vector<std::vector<SX> >& asens){
    casadi_assert_message(arg.size()==getNumInputs(),"FX::evalSX: dimension mismatch. You supplied " << arg.size() << " arguments instead of suspected " << getNumInputs() << ".");
    (*this)->evalSX(arg,res,fseed,fsens,aseed,asens);
  }

  void FX::evalMX(const std::vector<MX>& arg, std::vector<MX>& res, 
                  const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens, 
                  const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens){
    casadi_assert_message(arg.size()==getNumInputs(),"FX::evalMX: dimension mismatch. You supplied " << arg.size() << " arguments instead of suspected " << getNumInputs() << ".");
    (*this)->evalMX(arg,res,fseed,fsens,aseed,asens);
  }
                        
  void FX::eval(const std::vector<SX>& arg, std::vector<SX>& res, 
                const std::vector<std::vector<SX> >& fseed, std::vector<std::vector<SX> >& fsens, 
                const std::vector<std::vector<SX> >& aseed, std::vector<std::vector<SX> >& asens){
    casadi_assert_message(arg.size()==getNumInputs(),"FX::eval: dimension mismatch. You supplied " << arg.size() << " arguments instead of suspected " << getNumInputs() << ".");
    (*this)->evalSX(arg,res,fseed,fsens,aseed,asens);
  }

  void FX::eval(const std::vector<MX>& arg, std::vector<MX>& res, 
                const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens, 
                const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens){
    casadi_assert_message(arg.size()==getNumInputs(),"FX::eval: dimension mismatch. You supplied " << arg.size() << " arguments instead of suspected " << getNumInputs() << ".");
    (*this)->evalMX(arg,res,fseed,fsens,aseed,asens);
  }

  void FX::spEvaluate(bool fwd){
    (*this)->spEvaluate(fwd);
  }

  bool FX::spCanEvaluate(bool fwd){
    return (*this)->spCanEvaluate(fwd);
  }

  void FX::spInit(bool fwd){
    (*this)->spInit(fwd);
  }

  FX FX::derivative(int nfwd, int nadj){
    return (*this)->derivative(nfwd,nadj);
  }

  void FX::setDerivative(const FX& fcn, int nfwd, int nadj){
    (*this)->setDerivative(fcn,nfwd,nadj);
  }

  void FX::generateCode(const string& filename){
    (*this)->generateCode(filename);
  }

  const IOScheme& FX::inputScheme() const{
    return (*this)->inputScheme();
  }
  
  const IOScheme& FX::outputScheme() const{
    return (*this)->outputScheme();
  }
  
  IOScheme& FX::inputScheme(){
    return (*this)->inputScheme();
  }
  
  IOScheme& FX::outputScheme(){
    return (*this)->outputScheme();
  }
  
  const IOSchemeVector<DMatrix>& FX::input_struct() const{
    return (*this)->input_struct();
  }
  
  const IOSchemeVector<DMatrix>& FX::output_struct() const{
    return (*this)->output_struct();
  }
  
  IOSchemeVector<DMatrix>& FX::input_struct(){
    return (*this)->input_struct();
  }
  
  IOSchemeVector<DMatrix>& FX::output_struct(){
    return (*this)->output_struct();
  }  

  void FX::checkInputs() const {
    return (*this)->checkInputs();
  }

} // namespace CasADi

