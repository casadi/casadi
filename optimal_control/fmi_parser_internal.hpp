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

#ifndef FMI_PARSER_INTERNAL_HPP
#define FMI_PARSER_INTERNAL_HPP

#include "fmi_parser.hpp"
#include "optimica_ocp.hpp"
#include "xml_node.hpp"

/** \brief  Forward declarations */
class TiXmlElement;
class TiXmlNode;

namespace CasADi{
namespace OptimalControl{

class FMIParserInternal : public SharedObjectNode{
  public:
    /// Constructor
    explicit FMIParserInternal(const std::string& filename);
  
    /// destructor
    virtual ~FMIParserInternal(); 

    /// Parse from XML to C++ format
    OCP& parse();

    /// Get the OCP
    OCP& ocp();

    /// Get the OCP (const ref)
    const OCP& ocp() const;

    ///  Print representation
    virtual void repr(std::ostream &stream=std::cout) const;

    /// Print description 
    virtual void print(std::ostream &stream=std::cout) const;

  protected:
    
    /// Add model variables */
    void addModelVariables();

    /// Add binding equations */
    void addBindingEquations();

    /// Add dynamic equations
    void addDynamicEquations();

    /// Read an equation
    SX readExpr(const XMLNode& odenode);

    /// Read a variable
    Variable& readVariable(const XMLNode& node);

    /// Add initial equations
    void addInitialEquations();

    /// Add optimization */
    void addOptimization();
    void addObjectiveFunction(const XMLNode& onode);
    void addIntegrandObjectiveFunction(const XMLNode& onode);
    void addConstraints(const XMLNode& onode);
    void addIntervalStartTime(const XMLNode& onode);
    void addIntervalFinalTime(const XMLNode& onode);

    // NOTE 1: Joel: The FMIParserInternal class will later have to be changed to work with the MX class instead of SX, 
    //               therefore I had to change the implementation so that it is more generic

    // NOTE 2: Joel: Will there really ever be so many functions that it will motivate a binary search of the functions rather than a simple linear search?

    /// Look-up table mapping XML names to SX unary functions
    std::map<std::string,SX (*)(const SX&)> unary_;

    /// Look-up table mapping XML names to SX binary functions
    std::map<std::string,SX (*)(const SX&,const SX&)> binary_;

    /// The optimal control problem representation -- keep synchronized with the XML representation!
    OCP ocp_;

    /// Parsed XML document
    XMLNode document_;

    /// Filename
    std::string filename_;
};

} // namespace OptimalControl
} // namespace CasADi

#endif //FMI_PARSER_INTERNAL_HPP
