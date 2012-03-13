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

#include "flat_ocp_internal.hpp"
#include "variable_tools.hpp"

using namespace std;
namespace CasADi{

FlatOCP::FlatOCP(){
}
    
FlatOCP::FlatOCP(const std::string& filename){
  assignNode(new FlatOCPInternal(filename));
}

void FlatOCP::parse(){
  (*this)->parse();
}

FlatOCPInternal* FlatOCP::operator->(){
  return (FlatOCPInternal*)(OptionsFunctionality::operator->());
}

const FlatOCPInternal* FlatOCP::operator->() const{
  return (const FlatOCPInternal*)(OptionsFunctionality::operator->());
}

bool FlatOCP::checkNode() const{
  return dynamic_cast<const FlatOCPInternal*>(get())!=0;
}

Variable& FlatOCP::variable(const std::string& name){
  return (*this)->variable(name);
}

void FlatOCP::addVariable(const std::string& name, const Variable& var){
  (*this)->addVariable(name,var);
}

void FlatOCP::makeAlgebraic(const std::string& name){
  (*this)->makeAlgebraic(variable(name));
}

std::vector<Variable>& FlatOCP::x(){
  return (*this)->x_;
}

std::vector<Variable>& FlatOCP::xd(){
  return (*this)->xd_;
}

std::vector<Variable>& FlatOCP::xa(){
  return (*this)->xa_;
}

std::vector<Variable>& FlatOCP::xq(){
  return (*this)->xq_;
}

std::vector<Variable>& FlatOCP::y(){
  return (*this)->y_;
}

std::vector<Variable>& FlatOCP::p(){
  return (*this)->p_;
}

std::vector<Variable>& FlatOCP::u(){
  return (*this)->u_;
}

#ifdef NEW_FLAT_OCP

SXMatrix FlatOCP::dae() const{
  return (*this)->dae_;
}
    
SXMatrix FlatOCP::ode() const{
  return (*this)->ode_;
}
    
SXMatrix FlatOCP::alg() const{
  return (*this)->alg_;
}
    
SXMatrix FlatOCP::quad() const{
  return (*this)->quad_;
}
    
SXMatrix FlatOCP::dep() const{
  return (*this)->dep_;
}
    
SXMatrix FlatOCP::initial() const{
  return (*this)->initial_;
}

SXMatrix FlatOCP::mterm() const{
  return (*this)->mterm_;
}

SXMatrix FlatOCP::lterm() const{
  return (*this)->lterm_;
}

SXMatrix FlatOCP::path() const{
  return (*this)->path_;
}

DMatrix FlatOCP::path_min() const{
  return (*this)->path_min_;
}

DMatrix FlatOCP::path_max() const{
  return (*this)->path_max_;
}

SXMatrix FlatOCP::t() const{
  return (*this)->t_;
}

#else // NEW_FLAT_OCP

std::vector<SX>& FlatOCP::dae(){
  return (*this)->dae_;
}
    
std::vector<SX>& FlatOCP::ode(){
  return (*this)->ode_;
}
    
std::vector<SX>& FlatOCP::alg(){
  return (*this)->alg_;
}
    
std::vector<SX>& FlatOCP::quad(){
  return (*this)->quad_;
}
    
std::vector<SX>& FlatOCP::dep(){
  return (*this)->dep_;
}
    
std::vector<SX>& FlatOCP::initial(){
  return (*this)->initial_;
}

std::vector<SX>& FlatOCP::mterm(){
  return (*this)->mterm_;
}

std::vector<SX>& FlatOCP::lterm(){
  return (*this)->lterm_;
}

std::vector<SX>& FlatOCP::path(){
  return (*this)->path_;
}
  
std::vector<double>& FlatOCP::path_min(){
  return (*this)->path_min_;
}

std::vector<double>& FlatOCP::path_max(){
  return (*this)->path_max_;
}

SX FlatOCP::t() const{
  return (*this)->t_;
}

#endif // NEW_FLAT_OCP

double FlatOCP::t0() const{
  return (*this)->t0_;
}

bool FlatOCP::t0_free() const{
  return (*this)->t0_free_;
}

double FlatOCP::tf() const{
  return (*this)->tf_;
}

bool FlatOCP::tf_free() const{
  return (*this)->tf_free_;
}
  
void FlatOCP::set_t0(double t){
  (*this)->t0_ = t;
}

void FlatOCP::set_tf(double t){
  (*this)->tf_ = t;
}

void FlatOCP::set_t0_free(bool free){
  (*this)->t0_free_ = free;
}

void FlatOCP::set_tf_free(bool free){
  (*this)->tf_free_ = free;
}

void FlatOCP::eliminateDependent(){
  (*this)->eliminateDependent();
}

void FlatOCP::eliminateInterdependencies(){
  (*this)->eliminateInterdependencies();
}

void FlatOCP::sortDAE(){
  (*this)->sortDAE();
}

void FlatOCP::makeExplicit(){
  (*this)->makeExplicit();
}

std::vector<Variable> FlatOCP::x_all() const{
  return (*this)->x_all();
}

std::vector<SXMatrix> FlatOCP::daeArg() const{
  return (*this)->daeArg();
}
    
FX FlatOCP::daeFcn() const{
  return (*this)->daeFcn();
}

vector<SXMatrix> FlatOCP::substituteDependents(const vector<SXMatrix>& x) const{
  return (*this)->substituteDependents(x);
}

void FlatOCP::generateMuscodDatFile(const std::string& filename, const Dictionary& mc2_ops) const{
  (*this)->generateMuscodDatFile(filename, mc2_ops);
}

void FlatOCP::eliminateLagrangeTerms(){
  (*this)->eliminateLagrangeTerms();
}

void FlatOCP::eliminateQuadratureStates(){
  (*this)->eliminateQuadratureStates();
}

} // namespace CasADi
