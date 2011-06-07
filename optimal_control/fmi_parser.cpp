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

#include "fmi_parser_internal.hpp"

using namespace std;
namespace CasADi{
namespace OptimalControl{

FMIParser::FMIParser(){
}
    
FMIParser::FMIParser(const std::string& filename){
  assignNode(new FMIParserInternal(filename));
}

OCP& FMIParser::parse(){
  return (*this)->parse();
}

OCP& FMIParser::ocp(){
  return (*this)->ocp();
}

const OCP& FMIParser::ocp() const{
  return (*this)->ocp();
}

FMIParserInternal* FMIParser::operator->(){
  return (FMIParserInternal*)(SharedObject::operator->());
}

const FMIParserInternal* FMIParser::operator->() const{
  return (const FMIParserInternal*)(SharedObject::operator->());
}

bool FMIParser::checkNode() const{
  return dynamic_cast<const FMIParserInternal*>(get())!=0;
}

} // namespace OptimalControl
} // namespace CasADi
