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
#
# -*- coding: utf-8 -*-
import casadi as ca
from matplotlib import pyplot as plt
import numpy as np
from zipfile import ZipFile
from pathlib import Path
import os

# Simulating a bouncing ball with DaeBuilder and event handling
# Joel Andersson, 2025

# Start with an empty DaeBuilder instance
dae = ca.DaeBuilder('bouncing_ball')

# Model variables
t = dae.add('t', 'independent')
h = dae.add('h', 'output', dict(start = 5, initial = 'exact'))
v = dae.add('v', 'output', dict(start = 0, initial = 'exact'))

# Dynamic equations
dae.eq(dae.der(h), v)
dae.eq(dae.der(v), -9.81)
dae.disp(True)

# Event dynamics: When h < 0, reinitialize v to -0.8*v
dae.when(h < 0, [dae.reinit('v', -0.8*dae.pre(v))])
dae.disp(True)

# Default experiment
dae.set_start_time(0)
dae.set_stop_time(7)

# Simulate
tgrid = np.linspace(dae.start_time(), dae.stop_time(), 100)
sim = ca.integrator('sim', 'cvodes', dae.create(), 0, tgrid,
                    dict(transition = dae.transition()))
simres = sim(x0 = dae.start(dae.x()))

# Visualize the solution
plt.figure(1)
plt.clf()
plt.plot(tgrid, simres['xf'][0, :].T)
plt.grid()
plt.show()

# Export FMU
fmu_files = dae.export_fmu()
print('Generated files: {}'.format(fmu_files))

# Compile DLL
fmi_headers = Path(__file__).parent.parent.parent.parent \
  / 'external_packages' / 'FMI-Standard-3.0' / 'headers'
cfiles = " ".join([f for f in fmu_files if f.endswith('.c')])
sofile = dae.name() + '.so'
os.system(f'gcc --shared -fPIC -I{fmi_headers} {cfiles} -o {sofile}')
print(f'Compiled {sofile}')
fmu_files[sofile] = 'binaries/x86_64-linux'

# Package into an FMU
fmuname = dae.name() + '.fmu'
with ZipFile(fmuname, 'w') as fmufile:
    for f, arcpath in fmu_files.items():
      fmufile.write(f, arcname = arcpath + '/' + f)
      os.remove(f)
print(f'Created FMU: {fmuname}')

# Load the FMU in FMPy
try:
    import fmpy
    fmpy.dump(fmuname)
    # Simulate the generated FMU
    res = fmpy.simulate_fmu(fmuname)
    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.clf()
    plt.plot(res['time'], res['h'],'-', label = 'h')
    plt.xlabel('time')
    plt.legend()
    plt.grid()
    plt.show()

except ImportError as e:
   print('FMPy not installed. Skipping FMU simulation.')
