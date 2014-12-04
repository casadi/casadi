/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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

%{
#include <streambuf>
#include <ostream>
  namespace casadi {
    /** Stream buffer to allow printing to mex conveniently
        Adapted from answer to stackoverflow question 243696.
    */
    class mex_buf : public std::streambuf {
    public:
      mex_buf() {}
    protected:
      virtual int_type overflow(int_type ch) {
        if(ch != traits_type::eof()) {
          mexPrintf("%c", static_cast<char>(ch));
        }
        return ch;
      }
      /* virtual int sync() { // needed?
         mexEvalString("drawnow;");
         return 0;
      } */
      virtual std::streamsize xsputn(const char* s, std::streamsize num) {
        mexPrintf("%.*s", static_cast<int>(num), s);
        return num;
      }
    };

    // Corresponding output stream
    class mexostream : public std::ostream {
    protected:
      mex_buf buf;
    public:
      mexostream() : std::ostream(&buf) {}
    };

    // Instantiation (cf. cout)
    static mexostream mexout;
  } // namespace casadi
%}

%inline %{
/// int
template <>
int meta< int >::toCpp(GUESTOBJECT * p, int *m, swig_type_info *type) {
  if (is_a(p, type)) {
    int *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, type, 0) == -1) return false;
    *m=*mp;
    return true;
  }
  if (is_a(p, *meta< casadi::Matrix<int> >::name)) {
    casadi::Matrix<int> *temp;
    SWIG_ConvertPtr(p, (void **) &temp, *meta< casadi::Matrix<int> >::name, 0 );
    if (temp->numel()==1 && temp->size()==1) {
      *m = temp->data()[0];
      return true;
    }
    return false;
  } else {
    return false;
  }
}

/// double
template <>
int meta< double >::toCpp(GUESTOBJECT * p, double *m, swig_type_info *type) {
  if (is_a(p, *meta< double >::name)) {
    double *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< double >::name, 0) == -1)
      return false;
    *m=*mp;
    return true;
  }
  if (is_a(p, *meta< casadi::Matrix<double> >::name)) {
    casadi::Matrix<double> *temp;
    SWIG_ConvertPtr(p, (void **) &temp, *meta< casadi::Matrix<double> >::name, 0 );
    if (temp->numel()==1 && temp->size()==1) {
      *m = temp->data()[0];
      return true;
    }
    return false;
  } else if (is_a(p, *meta< casadi::Matrix<int> >::name)) {
    casadi::Matrix<int> *temp;
    SWIG_ConvertPtr(p, (void **) &temp, *meta< casadi::Matrix<int> >::name, 0 );
    if (temp->numel()==1 && temp->size()==1) {
      *m = temp->data()[0];
      return true;
    }
    return false;
  } else {
    return false;
  }
}

// Explicit intialization of these two member functions, so we can use them in meta< casadi::SXElement >
template<> int meta< casadi::SX >::toCpp(GUESTOBJECT *p, casadi::SX *, swig_type_info *type);

/// casadi::SX
template <>
int meta< casadi::SXElement >::toCpp(GUESTOBJECT *p, casadi::SXElement *s, swig_type_info *type) {
  if (is_a(p, *meta< casadi::SXElement >::name)) {
    casadi::SXElement *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< casadi::SXElement >::name, 0) == -1)
      return false;
    *s=*mp;
    return true;
  }
  if (meta< double >::couldbe(p)) {
    double res;
    if (!meta< double >::toCpp(p, &res, *meta< double >::name)) return false;
    *s=casadi::SXElement(res);
  } else {
    return false;
  }
  return true;
}

/// casadi::Matrix<int>
template <>
int meta< casadi::Matrix<int> >::toCpp(GUESTOBJECT * p,casadi::Matrix<int> *m, swig_type_info *type) {
  if (is_a(p, *meta< casadi::Matrix<int> >::name)) {
    casadi::Matrix<int> *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< casadi::Matrix<int> >::name, 0) == -1)
      return false;
    *m=*mp;
    return true;
  }
  if (meta< int >::couldbe(p)) {
    int t;
    int res = meta< int >::toCpp(p, &t, *meta< int >::name);
    *m = t;
    return res;
  } else {
    //SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

/// casadi::Matrix<double>
template <>
int meta< casadi::Matrix<double> >::toCpp(GUESTOBJECT * p, casadi::Matrix<double> *m, swig_type_info *type) {
  if (is_a(p, *meta< casadi::Matrix<double> >::name)) {
    casadi::Matrix<double> *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< casadi::Matrix<double> >::name, 0) == -1)
      return false;
    *m=*mp;
    return true;
  }
  if (is_a(p, *meta< casadi::Matrix<int> >::name)) {
    casadi::Matrix<int> *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< casadi::Matrix<int> >::name, 0) == -1)
      return false;
    *m=*mp;
    return true;
  }
  if (meta< double >::couldbe(p)) {
    double t;
    int res = meta< double >::toCpp(p, &t, *meta< double >::name);
    *m = t;
    return res;
  } else {
    //SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

/// casadi::SX
template <>
int meta< casadi::SX >::toCpp(GUESTOBJECT * p,casadi::SX *m, swig_type_info *type) {
  if (is_a(p, *meta< casadi::SX >::name)) {
    casadi::SX *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< casadi::SX >::name, 0) == -1)
      return false;
    *m=*mp;
    return true;
  }
  if (is_a(p, *meta< casadi::SXElement >::name)) {
    casadi::SXElement *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< casadi::SXElement >::name, 0) == -1)
      return false;
    *m=*mp;
    return true;
  }
  casadi::DMatrix mt;
  if(meta< casadi::Matrix<double> >::toCpp(p, &mt, *meta< casadi::Matrix<double> >::name)) {
    *m = casadi::SX(mt);
  } else {
    //SWIG_Error(SWIG_TypeError, "asSX: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
 }

meta_vector(std::vector<casadi::SXElement>);
meta_vector(casadi::SXElement);
meta_vector(casadi::Matrix< casadi::SXElement >);
meta_vector(std::vector< casadi::Matrix< casadi::SXElement > >);

/// casadi::MX
template <>
int meta< casadi::MX >::toCpp(GUESTOBJECT * p,casadi::MX *m, swig_type_info *type) {
  if (is_a(p, *meta< casadi::MX >::name)) {
    casadi::MX *mp;
    if (SWIG_ConvertPtr(p, (void **) &mp, *meta< casadi::MX >::name, 0) == -1)
      return false;
    *m=*mp;
    return true;
  }
  casadi::DMatrix mt;
  if(meta< casadi::Matrix<double> >::toCpp(p, &mt, *meta< casadi::Matrix<double> >::name)) {
    *m = casadi::MX(mt);
    return true;
  }
  return false;
}

meta_vector(casadi::MX);
meta_vector(std::vector< casadi::MX >);
%}
