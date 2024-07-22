#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
dae.set_ode('x1', (1 - x2 * x2) * x1 - x2 + u)
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
