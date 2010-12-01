#include "printable_object.hpp"
#include <typeinfo>
#include <sstream>

using namespace std;
namespace CasADi{

ostream& operator<<(ostream &stream, const PrintableObject& obj){
  obj.repr(stream);
  return stream;
}

void PrintableObject::repr(std::ostream &stream) const{
  // Print description by default
  print(stream);
}

void PrintableObject::print(std::ostream &stream) const{
  // Print name by default
  stream << typeid(this).name();
}

string PrintableObject::getRepresentation() const{
  stringstream ss;
  repr(ss);
  return ss.str();
}

string PrintableObject::getDescription() const{
  stringstream ss;
  print(ss);
  return ss.str();  
}

    
} // namespace CasADi
    


