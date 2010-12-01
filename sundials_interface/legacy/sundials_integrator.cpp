#include "sundials_integrator.hpp"

namespace OPTICON{

SundialsIntegrator::SundialsIntegrator(OCP_old &ocp_) : IntegratorNode_old(ocp_){
  addOption("mupper",                    OT_INTEGER,  Option()); // upper band-width of banded jacobians
  addOption("mlower",                    OT_INTEGER,  Option()); // lower band-width of banded jacobians

  addOption("linear_solver",             OT_STRING, "dense"); // "dense", "band" or "sparse"
  addOption("sparse_solver",             OT_STRING, "spgmr"); // "spgmr", "spbcg" or  "sptfqmr"
  addOption("pretype",                   OT_STRING, "none"); // "none", "left", "right", "both"

}

SundialsIntegrator::~SundialsIntegrator(){

}

void SundialsIntegrator::init(){
  IntegratorNode_old::init();
}

void SundialsIntegrator::sundialsAssert(int flag) const{
  if(flag<0) throw "Error in Sundials function";
}

} // namespace OPTICON

