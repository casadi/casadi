#ifndef MUSCOD_AUX_HPP
#define MUSCOD_AUX_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

/**
  Auxilliary functions for reading and writing DAT/RES files.
  Joel Andersson, 2010
*/

namespace CasADi{

/** Write a vector valued entry to file */
template<typename A>
void print_dat(std::ostream &stream,const std::string& name,const std::vector<A>& v){
  stream << name << std::endl;
  for(unsigned int i=0; i<v.size(); ++i)
    stream << i << ": " << v[i] << std::endl;
  stream << std::endl;
}

/** Write a scalar valued entry to file */
template<typename A>
void print_dat(std::ostream &stream,const std::string& name,const A& v){
  stream << name << std::endl << v << std::endl << std::endl;
}

/** Write a scalar/vector valued entry with an index argument to file */
template<typename A,typename I>
void print_dat(std::ostream &stream,const std::string& name,const A& v,I i){
  std::stringstream ss;
  ss << name << "(" << i << ")";
  print_dat(stream,ss.str(),v);
}

/** Write a scalar/vector valued entry with two index arguments to file */
template<typename A,typename I1,typename I2>
void print_dat(std::ostream &stream,const std::string& name,const A& v,I1 i1,I2 i2){
  std::stringstream ss;
  ss << name << "(" << i1 << "," << i2 << ")";
  print_dat(stream,ss.str(),v);
}

/** Structure to hold the information about an entry */
struct ResEntry{
  /** Type of argument */
  enum {DOUBLE, INT, STRING} type;
  
  /** Is the argument vector valued */
  bool isVector;
  
  /** Index of the argument */
  int index;
  
};


} // namespace CasADi

#endif //MUSCOD_AUX_HPP
