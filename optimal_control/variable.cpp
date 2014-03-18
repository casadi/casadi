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

#include "variable_internal.hpp"
#include "../symbolic/casadi_exception.hpp"

using namespace std;
namespace CasADi{

  Variable::Variable(){
  }

  Variable::Variable(const string& name){
    assignNode(new VariableInternal(name));
  }

  Variable::~Variable(){
  }

  VariableInternal* Variable::operator->(){
    return static_cast<VariableInternal*>(SharedObject::operator->());
  }

  const VariableInternal* Variable::operator->() const{
    return static_cast<const VariableInternal*>(SharedObject::operator->());
  }
  
  SX Variable::der() const{
    return (*this)->der_;
  }

  SX Variable::var() const{
    return (*this)->var_;
  }

  SX Variable::binding(bool derivative) const{
    if(derivative){
      return (*this)->der_binding_;
    } else {
      return (*this)->binding_;
    }
  }

  string Variable::getName() const{
    return (*this)->var_.getName();
  }

  Variability Variable::getVariability() const{
    return (*this)->variability_;
  }

  void Variable::setVariability(Variability variability){
    (*this)->variability_ = variability;
  }

  Causality Variable::getCausality() const{
    return (*this)->causality_;
  }

  void Variable::setCausality(Causality causality){
    (*this)->causality_ = causality;
  }

  Category Variable::getCategory() const{
    return (*this)->category_;
  }

  void Variable::setCategory(Category category){
    (*this)->category_ = category;
  }

  Alias Variable::getAlias() const{
    return (*this)->alias_;
  }

  void Variable::setAlias(Alias alias){
    (*this)->alias_ = alias;
  }
    
  const string& Variable::getDescription() const{
    return (*this)->description_;
  }

  void Variable::setDescription(const string& description){
    (*this)->description_ = description;
  }
    
  int Variable::getValueReference() const{
    return (*this)->valueReference_;
  }

  void Variable::setValueReference(int valueReference){
    (*this)->valueReference_ = valueReference;
  }
    
  double Variable::getMin() const{
    return (*this)->min_;
  }

  double& Variable::min(){
    return (*this)->min_;
  }

  void Variable::setMin(double min){
    (*this)->min_ = min;
  }
    
  double Variable::getMax() const{
    return (*this)->max_;
  }

  double& Variable::max(){
    return (*this)->max_;
  }

  void Variable::setMax(double max){
    (*this)->max_ = max;
  }
    
  double Variable::getNominal() const{
    return (*this)->nominal_;
  }

  void Variable::setNominal(double nominal){
    (*this)->nominal_ = nominal;
  }

  double& Variable::nominal(){
    return (*this)->nominal_;
  }
    
  double Variable::getStart() const{
    return (*this)->start_;
  }

  void Variable::setStart(double start){
    (*this)->start_ = start;
  }

  double& Variable::start(){
    return (*this)->start_;
  }

  double Variable::getDerivativeStart() const{
    return (*this)->derivative_start_;
  }

  double& Variable::derivativeStart(){
    return (*this)->derivative_start_;
  }

  double Variable::getInitialGuess() const{
    return (*this)->initial_guess_;
  }

  double& Variable::initialGuess(){
    return (*this)->initial_guess_;
  }

  void Variable::setInitialGuess(double initial_guess){
    (*this)->initial_guess_ = initial_guess;
  }

  void Variable::setDerivativeStart(double start){
    (*this)->derivative_start_ = start;
  }
    
  const string& Variable::getUnit() const{
    return (*this)->unit_;
  }

  string& Variable::unit(){
    return (*this)->unit_;
  }

  void Variable::setUnit(const string& unit){
    (*this)->unit_ = unit;
  }
    
  const string& Variable::getDisplayUnit() const{
    return (*this)->displayUnit_;
  }

  string& Variable::displayUnit(){
    return (*this)->displayUnit_;
  }

  void Variable::setDisplayUnit(const string& displayUnit){
    (*this)->displayUnit_ = displayUnit;
  }

  void Variable::setExpression(const SX& v){
    (*this)->var_ = v;
  }

  void Variable::setDerivative(const SX& d){
    (*this)->der_ = d;
  }

  void Variable::setBinding(const SX& binding, bool derivative){
    if(derivative){
      (*this)->der_binding_ = binding;
    } else {
      (*this)->binding_ = binding;
    }
  }

  bool Variable::checkNode() const{
    return dynamic_cast<const VariableInternal*>(get())!=0;
  }

  SX Variable::atTime(double t, bool allocate) const{
    return (*this)->atTime(t,allocate);
  }

  SX Variable::atTime(double t, bool allocate){
    return (*this)->atTime(t,allocate);
  }

  int Variable::index() const{
    return (*this)->index_;
  }

  void Variable::setIndex(int ind){
    (*this)->index_ = ind;
  }
    
  bool Variable::isDifferential() const{
    return (*this)->is_differential_;
  }

  void Variable::setDifferential(bool is_differential){
    (*this)->is_differential_ = is_differential;
  }

  bool Variable::getFree() const{
    return (*this)->free_;
  }

  bool& Variable::free(){
    return (*this)->free_;
  }

  void Variable::setFree(bool free){
    (*this)->free_ = free;
  }


} // namespace CasADi
