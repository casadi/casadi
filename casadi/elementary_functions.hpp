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

#ifndef ELEMENTARY_FUNCTIONS_HPP
#define ELEMENTARY_FUNCTIONS_HPP

#include <iostream>
#include <string>

namespace CasADi{
/// Enum for quick access to any node
enum OPERATION2{
  ADD,  SUB,  MUL,  DIV,
  NEG,  EXP,  LOG,  POW,  
  SQRT,  SIN,  COS,  TAN,  
  ASIN,  ACOS,  ATAN,  
  STEP,  
  FLOOR,  CEIL,  
  EQUALITY,  ERF,  FMIN,  FMAX
};

/// Unary function
template<int I>
class UnFun{
  public:
    static void print(std::ostream &stream, const std::string& x);
};


/// Printing unary function
template<int I>
void uprint(std::ostream &stream, const std::string& x);

/// Printing binary function
template<int I>
void bprint(std::ostream &stream, const std::string& x, const std::string& y){ uprint<I>(stream,x); }

/// Unary function
template<typename T, int I>
T ufcn(const T& x);

/// Binary function
template<typename T, int I>
T bfcn(const T& x, const T& y){ return ufcn<T,I>(x);}

/// Partial derivatives - unary function
template<typename T, int I>
void uder(const T& x, T& f, T* d);

/// Partial derivatives - binary function
template<typename T, int I>
void bder(const T& x, const T& y, T& f, T* d){ uder<T,I>(x,f,d); d[1]=0; } 

/// Printing - template specialization
void bprint<ADD>(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "+" << y << ")";}
void bprint<SUB>(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "-" << y << ")";}
void bprint<MUL>(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "*" << y << ")";}
void bprint<DIV>(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "/" << y << ")";}
void uprint<EXP>(std::ostream &stream, const std::string& x){ stream << "exp(" << x << ")";}
void uprint<LOG>(std::ostream &stream, const std::string& x){ stream << "log(" << x << ")";}
void bprint<POW>(std::ostream &stream, const std::string& x, const std::string& y){stream << "pow(" << x << "," << y << ")";}
void uprint<SQRT>(std::ostream &stream, const std::string& x){ stream << "sqrt(" << x << ")";}
void uprint<SIN>(std::ostream &stream, const std::string& x){ stream << "sin(" << x << ")";}
void uprint<COS>(std::ostream &stream, const std::string& x){ stream << "cos(" << x << ")";}
void uprint<TAN>(std::ostream &stream, const std::string& x){ stream << "tan(" << x << ")";}

/// Elementary functions - partial specialization
template<typename T> T bfcn<T,ADD>(const T& x, const T& y){ return x+y;}
template<typename T> T bfcn<T,SUB>(const T& x, const T& y){ return x-y;}
template<typename T> T bfcn<T,MUL>(const T& x, const T& y){ return x*y;}
template<typename T> T bfcn<T,DIV>(const T& x, const T& y){ return x/y;}
template<typename T> T ufcn<T,EXP>(const T& x){ return exp(x);}
template<typename T> T ufcn<T,LOG>(const T& x){ return log(x);}
template<typename T> T bfcn<T,POW>(const T& x, const T& y){ return pow(x,y);}
template<typename T> T ufcn<T,SQRT>(const T& x){ return sqrt(x);}
template<typename T> T ufcn<T,SIN>(const T& x){ return sin(x);}
template<typename T> T ufcn<T,COS>(const T& x){ return cos(x);}
template<typename T> T ufcn<T,TAN>(const T& x){ return tan(x);}

/// Derivatives - partial specialization
template<typename T> void bder<T,ADD>(const T& x, const T& y, T& f, T* d){f = x+y; d[0]=1; d[1]=1; }
template<typename T> void bder<T,SUB>(const T& x, const T& y, T& f, T* d){f = x-y; d[0]=1; d[1]=-1; }
template<typename T> void bder<T,ADD>(const T& x, const T& y, T& f, T* d){f = x*y; d[0]=1; d[1]=y; }
template<typename T> void bder<T,SUB>(const T& x, const T& y, T& f, T* d){f = x/y; d[0]=1; d[1]=x; }
template<typename T> void uder<T,EXP>(const T& x, T& f, T* d){f = exp(x); d[0]=f;}
template<typename T> void uder<T,LOG>(const T& x, T& f, T* d){f = log(x); d[0]=1/x;}
template<typename T> void bder<T,POW>(const T& x, const T& y, T& f, T* d){f = pow(x,y); d[0]=y*f/x; d[1]=log(x)*f; }
template<typename T> void uder<T,SQRT>(const T& x, T& f, T* d){f = sqrt(x); d[0]=1/(2*f);}
template<typename T> void uder<T,SIN>(const T& x, T& f, T* d){f = sin(x); d[0]=cos(x);}
template<typename T> void uder<T,COS>(const T& x, T& f, T* d){f = cos(x); d[0]=-sin(x);}
template<typename T> void uder<T,TAN>(const T& x, T& f, T* d){f = cos(x); T cosx = cos(x); d[0] = 1/(cosx*cosx);}

} // namespace CasADi

#endif //ELEMENTARY_FUNCTIONS_HPP
