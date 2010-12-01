#include "simulator.hpp"
#include "simulator_internal.hpp"
#include "integrator_internal.hpp"
#include <cassert>

namespace CasADi{

Simulator::Simulator(){
}

Simulator::Simulator(const Integrator& integrator, const FX& output_fcn, const vector<double>& grid){
  assignNode(new SimulatorInternal(integrator,output_fcn,grid));
}

Simulator::Simulator(const Integrator& integrator, const vector<double>& grid){
  // Create a dummy function (returns the whole state)
  SXMatrix t("t");
  SXMatrix x("x",integrator->nx_);
  SXMatrix p("p",integrator->np_);
  vector<SXMatrix> arg(OUTPUT_NUM_IN);  arg[OUTPUT_T] = t; arg[OUTPUT_X] = x; arg[OUTPUT_P] = p;
  
  // Create the output function
  SXFunction output_fcn(arg,vector<SXMatrix>(1,x));
  assignNode(new SimulatorInternal(integrator,output_fcn,grid));
}

Simulator::Simulator(const Simulator& ref) : FX((const FX&)ref){
}

SimulatorInternal* Simulator::operator->(){
  return (SimulatorInternal*)(FX::operator->());
}

const SimulatorInternal* Simulator::operator->() const{
   return (const SimulatorInternal*)(FX::operator->()); 
}

void Simulator::assertNode() const{
  if(!dynamic_cast<const SimulatorInternal*>(get()))
    throw CasadiException("Simulator::assertNode");
}


} // namespace CasADi

