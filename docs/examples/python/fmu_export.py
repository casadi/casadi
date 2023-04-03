#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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
# -*- coding: utf-8 -*-
from casadi import *

# Example how to export FMUs from CasADi
# Joel Andersson, joel@jaeandersson.com

# Start with an empty DaeBuilder instance
dae = DaeBuilder('vdp')

# States
x1 = dae.add_x('x1')
x2 = dae.add_x('x2')
q = dae.add_q('q')

# Input
u = dae.add_u('u')

# Set ODE right-hand-sides
dae.set_ode('x1', (1 - x2 * x2)*x1 - x2 + u)
dae.set_ode('x2', x1)
dae.set_ode('q', x1**2 + x2**2 + u**2)

# Specify initial conditions
dae.set_start('x1', 0)
dae.set_start('x2', 0)

# Add bounds
dae.set_min('u', -0.75)
dae.set_max('u', 1)
dae.set_start('u', 0)

# Print DAE
dae.disp(True)

# Export FMU
funcs = dae.export_fmu()
print('Generated files: {}'.format(funcs))

# Compile DLL
import os
os.system('gcc --shared -fPIC -I../../../external_packages/FMI-Standard-3.0/headers/ vdp.c vdp_wrap.c -o vdp.so')
print('Compiled vdp.so')

# Package into an FMU
import zipfile
with zipfile.ZipFile('vdp_generated.fmu', 'w') as fmufile:
    # Add generated files to the archive
    for f in funcs:
      arcname = f if f == 'modelDescription.xml' else 'sources/' + f
      fmufile.write(f, arcname = arcname)
      os.remove(f)
    # Add compile DLL to the archive (assume Linux 64 but)
    fmufile.write('vdp.so', arcname = 'binaries/x86-linux/vdp.so')
    os.remove('vdp.so')
print('Created FMU: vdp_generated.fmu')