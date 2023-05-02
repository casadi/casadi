#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
import re
import os

import casadi


with open('compiler.sh','w') as compilerscript:
    compilerscript.write('#!/bin/bash\n')
    compilerscript.write(f"""g++ "$1" -D_GLIBCXX_USE_CXX11_ABI=0 {casadi.CasadiMeta.compiler_flags()} -I{casadi.GlobalOptions.getCasadiIncludePath()} -L{casadi.GlobalOptions.getCasadiPath()} -Wl,-rpath,{casadi.GlobalOptions.getCasadiPath()} -lcasadi -o "$2" """)

