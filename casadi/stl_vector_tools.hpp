#ifndef STL_VECTOR_TOOLS_HPP
#define STL_VECTOR_TOOLS_HPP

#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <limits>
#include "casadi_exception.hpp"

namespace std{

/** \brief Enables flushing an STL vector to a stream
  \author Joel Andersson 
  \date 2010
*/

template<typename T>
ostream& operator<<(ostream &stream, const vector<T> &v){

  // check if empty
  if(v.empty()) return stream << "<empty vector>";

  // print vector style
  stream << "[" << v.size() << "]"; // Print dimension 

  // Print elements
  stream << "(";
  for(unsigned int i=0; i<v.size()-1; ++i)
    stream << v[i] << ",";
  if(!v.empty()) stream << v.back();
  stream << ")";
  
  return stream;
}

// template<typename T>
// vector<T> operator+(const vector<T> &v1, const vector<T> &v2){
//   assert(v1.size() == v2.size());
//   vector<T> res(v1.size());
//   for(unsigned int i=0; i<v1.size(); ++i)
//     res[i] = v1[i] + v2[i];
//   return res;
// }
// 
// template<typename T>
// vector<T>& operator+=(vector<T> &v1, const vector<T> &v2){
//   assert(v1.size() == v2.size());
//   for(unsigned int i=0; i<v1.size(); ++i)
//     v1[i] += v2[i];
//   return v1;
// }

// template<typename T>
// vector<T> operator-(const vector<T> &v1, const vector<T> &v2){
//   assert(v1.size() == v2.size());
//   vector<T> res(v1.size());
//   for(unsigned int i=0; i<v1.size(); ++i)
//     res[i] = v1[i] - v2[i];
//   return res;
// }
// 
// template<typename T>
// vector<T>& operator-=(vector<T> &v1, const vector<T> &v2){
//   assert(v1.size() == v2.size());
//   for(unsigned int i=0; i<v1.size(); ++i)
//     v1[i] -= v2[i];
//   return v1;
// }

template<typename T>
void write_matlab(ostream &stream, const vector<T> &v){
  copy(v.begin(), v.end(), ostream_iterator<T>(stream, " "));
}

template<typename T>
void write_matlab(ostream &stream, const vector<vector<T> > &v){
  for(unsigned int i=0; i<v.size(); ++i){    
    copy(v[i].begin(), v[i].end(), ostream_iterator<T>(stream, " "));
    stream << endl;
  }
}

template<typename T>
void read_matlab(istream &stream, vector<T> &v){
  v.clear();
  
  while(!stream.eof()) {
    T val;
    stream >> val;
    if(stream.fail()){
      stream.clear();
      string s;
      stream >> s;
      if(s.compare("inf") == 0)
	val = numeric_limits<T>::infinity();
      else
	break;
    }
    v.push_back(val);
  }
}

template<typename T>
void read_matlab(ifstream &file, vector<vector<T> > &v){
  v.clear();
  string line;
    while(!getline(file, line, '\n').eof()) {
        istringstream reader(line);
        vector<T> lineData;

        string::const_iterator i = line.begin();

        while(!reader.eof()) {
            T val;
            reader >> val;
            if(reader.fail()){
		reader.clear();
		string s;
		reader >> s;
		if(s.compare("inf") == 0)
		  val = numeric_limits<T>::infinity();
		else
		  break;
	    }

            lineData.push_back(val);
        }
        v.push_back(lineData);
    }
}

// Push back using the stream operator (convenient when adding many elements)
// template<typename T, typename E>
// vector<T> & operator<<(vector<T> &v, const E &el){
//   v.push_back(el);
//   return v;
// }

} // namespace std

namespace CasADi{
  
// Linspace
template<typename T, typename F, typename L>
void linspace(std::vector<T> &v, const F& first, const L& last){
  if(v.size()<2) throw CasadiException("std::linspace: vector must contain at least two elements");
  
  // Increment
  T increment = (last-first)/T(v.size()-1);
  
  v[0] = first;
  for(unsigned i=1; i<v.size()-1; ++i)
    v[i] = v[i-1] + increment;
  v[v.size()-1] = last;
}

} // namespace CasADi



#endif // STL_VECTOR_TOOLS_HPP
