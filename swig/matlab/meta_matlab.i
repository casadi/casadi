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

