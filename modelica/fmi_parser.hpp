#ifndef FMI_PARSER_HPP
#define FMI_PARSER_HPP

#include "xml_parser.hpp"
#include "optimica_ocp.hpp"

namespace CasADi{
namespace Modelica{

class FMIParser : public XMLParser{

public:
FMIParser(const std::string& filename);
virtual ~FMIParser(); // destructor

/** \brief Parse from XML to C++ format */
OCP& parse();

/** \brief Get the OCP */
OCP& ocp();

/** \brief Get the OCP (const ref)*/
const OCP& ocp() const;

protected:

/** \brief  Add model variables */
void addModelVariables();

/** \brief  Add binding equations */
void addBindingEquations();

/** \brief  Add dynamic equations */
void addDynamicEquations();

/** \brief  Read an equation */
SX readExpr_new(const XMLNode& odenode);

/** \brief  Read a variable */
Variable readVariable(const XMLNode& node) const;

/** \brief  Add initial equations */
void addInitialEquations();

/** \brief  Add optimization */
void addOptimization();
void addObjectiveFunction(const XMLNode& onode);
void addConstraints(const XMLNode& onode);
void addIntervalStartTime(const XMLNode& onode);
void addIntervalFinalTime(const XMLNode& onode);

// NOTE 1: Joel: The FMIParser class will later have to be changed to work with the MX class instead of SX, 
//               therefore I had to change the implementation so that it is more generic

// NOTE 2: Joel: Will there really ever be so many functions that it will motivate a binary search of the functions rather than a simple linear search?

/// Look-up table mapping XML names to SX unary functions
std::map<std::string,SX (*)(const SX&)> unary_;

/// Look-up table mapping XML names to SX binary functions
std::map<std::string,SX (*)(const SX&,const SX&)> binary_;

/** \brief  The optimal control problem representation -- keep synchronized with the XML representation! */
OCP ocp_;
};

} // namespace Modelica
} // namespace CasADi

#endif //FMI_PARSER_HPP
