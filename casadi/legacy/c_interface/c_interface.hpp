/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef C_INTERFACE_HPP
#define C_INTERFACE_HPP

#include <vector>
#include <iostream>
#include <sstream>

// C++ to C error handling
#define c_fcn_wrapper(name,x) \
  try{ \
    x \
    return 0; \
  } catch(exception &e){ \
    cerr << "Exception thrown in module \"" << name << "\": " << e.what() << endl; \
  } catch(const char* str){ \
    cerr << "Error (old style!) in module \"" << name << "\": " << str << endl; \
  } catch(...){ \
    cerr << "Unknown error in module \"" << name << "\"" << endl; \
  } \
  return 1;

// For constructors
#define c_constructor_wrapper(name,x) \
try{ \
  x \
  } catch(exception &e){ \
    cerr << "Exception thrown in module \"" << name << "\": " << e.what() << endl; \
  } catch(const char* str){ \
    cerr << "Error (old style!) in module \"" << name << "\": " << str << endl; \
  } catch(...){ \
    cerr << "Unknown error in module \"" << name << "\"" << endl; \
  } \
  return 0; \

  
#endif // C_INTERFACE_HPP

