#include "muscod_internal.hpp"

#include <cassert>
#include <limits>

using namespace std;

namespace CasADi{

MuscodInterface::MuscodInterface(){ 
}

MuscodInterface::MuscodInterface(muscodSetupFcn setupFcn){
  assignNode(new MuscodInternal(setupFcn));
}

MuscodInternal* MuscodInterface::operator->(){
  return (MuscodInternal*)(OptionsFunctionality::operator->());
}

const MuscodInternal* MuscodInterface::operator->() const{
   return (const MuscodInternal*)(OptionsFunctionality::operator->()); 
}

void MuscodInterface::assertNode() const{
  if(!dynamic_cast<const MuscodInternal*>(get()))
    throw CasadiException("MuscodInterface::assertNode");
}

void MuscodInterface::solve(){
  (*this)->solve();
}


} // namespace CasADi

