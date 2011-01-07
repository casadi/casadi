%{
#include "modelica/variable.hpp"
#include "modelica/ocp_variables.hpp"
#include "modelica/optimica_ocp.hpp"
#include "modelica/fmi_parser.hpp"
%}

%include "modelica/variable.hpp"
%include "modelica/ocp_variables.hpp"
%include "modelica/optimica_ocp.hpp"

namespace CasADi{
  namespace Modelica{

class FMIParser{
public:
  FMIParser(const std::string& filename);    // constructor
  virtual ~FMIParser(); // destructor
  OCP& parse();
};

%extend FMIParser {
std::string __str__()  { return $self->getDescription(); }


}


} // namespace Modelica
} // namespace CasADi

