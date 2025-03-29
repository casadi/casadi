#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
from casadi import *

import subprocess
import os
import unittest
from helpers import *

class Daebuildertests(casadiTestCase):

  def test_reference_fmus(self):
    if "ghc-filesystem" not in CasadiMeta.feature_list(): return
    for name in ["VanDerPol2","VanDerPol3"]:
        fmu_file = "../data/" + name + ".fmu"
        if not os.path.exists(fmu_file):
            print("Skipping test_reference_fmus, resource not available")
            return
        dae = DaeBuilder("car",fmu_file)
        dae.disp(True)
        x0 = SX.sym('x0')
        x1 = SX.sym('x1')
        mu = SX.sym("mu")
        
        x = vertcat(x0,x1)
        c = mu
        f = dae.create('f',['x'],['ode'])
        f_ref = Function('f',[x],[vertcat(x1,1 * ((1 - x0 * x0) * x1) - x0)])
        test_point = [vertcat(1.1,1.3)]
        self.checkfunction(f,f_ref,inputs=[vertcat(1.1,1.3)],digits=7)
        if not name.endswith("3"):
            self.check_serialize(f,inputs=test_point)
  
  def test_cstr(self):
    fmu_file = "../data/cstr.fmu"
    if not os.path.exists(fmu_file):
        print("Skipping test_fmu_zip, resource not available")
        return
    unzipped_name = "cstr"
    unzipped_path = os.path.join(os.getcwd(), unzipped_name)
    import shutil
    if os.path.isdir(unzipped_path): shutil.rmtree(unzipped_path)
    import zipfile
    with zipfile.ZipFile(fmu_file, 'r') as zip_ref:
        zip_ref.extractall(unzipped_name)
    dae = DaeBuilder("cstr",unzipped_name)
    f = dae.create('f',['x','u'],['ode'])
    with self.assertInException("No stats available"):
        f.stats()
    dae = None
    
    C_A = SX.sym('C_A')
    C_B = SX.sym('C_B')
    V = 1
    k = 0.5
    
    q_in = SX.sym('q_in')
    C_A_in = SX.sym('C_A_in')
    C_B_in = SX.sym('C_B_in')

    x = vertcat(C_A,C_B)
    u = vertcat(C_A_in,C_B_in,q_in)    
    ode = vertcat((q_in*(C_A_in - C_A) - V*k*C_A*C_B)/V,(q_in*(C_B_in - C_B) + V*k*C_A*C_B)/V)


    f_ref = Function('f',[x,u],[ode],['x','u'],['ode'])
    
    test_point = [vertcat(1.1,1.3),vertcat(1.7,1.11,1.13)]
    self.checkfunction(f,f_ref,inputs=test_point,digits=4,hessian=False,evals=1)
    print(f.stats())
          
  def test_rumoca(self):
    if "rumoca" not in CasadiMeta.feature_list(): return
    rumoca = os.path.join(GlobalOptions.getCasadiPath(),'rumoca')
    rumoca_exe = os.path.join(GlobalOptions.getCasadiPath(),'rumoca.exe')
    if os.path.exists(rumoca) or os.path.exists(rumoca_exe):
        p = subprocess.run([rumoca,"-t","../assets/casadi_daebuilder.jinja","-m", "../assets/hello_world.mo"]) 

  def test_fmu_zip(self):
    fmu_file = "../data/VanDerPol2.fmu"
    if not os.path.exists(fmu_file):
        print("Skipping test_fmu_zip, resource not available")
        return
    unzipped_name = "VanDerPol2"
    unzipped_path = os.path.join(os.getcwd(), unzipped_name)
    for serialize_mode in ["link","embed"]:
        for use_zip in [True]:#,False]:
            if use_zip:
                if "ghc-filesystem" in CasadiMeta.feature_list():
                    dae = DaeBuilder("cstr",fmu_file,{"resource_serialize_mode": serialize_mode})
                else:
                    with self.assertInException("passing fmu files to DaeBuilder is unsupported"):
                        dae = DaeBuilder("cstr",fmu_file)
                    continue
            else:
                import shutil
                if os.path.isdir(unzipped_path): shutil.rmtree(unzipped_path)
                import zipfile
                with zipfile.ZipFile(fmu_file, 'r') as zip_ref:
                    zip_ref.extractall(unzipped_name)
                print('Unzipped %s into %s' % (fmu_file, unzipped_path))
                dae = DaeBuilder("cstr",unzipped_path)
            dae.disp(True)
            f = dae.create('f',['x'],['ode'])
            dae = None
            
            x0 = SX.sym('x0')
            x1 = SX.sym('x1')
            mu = SX.sym("mu")
            
            x = vertcat(x0,x1)
            c = mu
            f_ref = Function('f',[x],[vertcat(x1,1 * ((1 - x0 * x0) * x1) - x0)])
            test_point = [vertcat(1.1,1.3)]

            self.checkfunction_light(f,f_ref,inputs=test_point,digits=7)
            
            f.save('f.casadi')
            
            f = None
            
            f = Function.load("f.casadi")
            self.checkfunction_light(f,f_ref,inputs=test_point,digits=7)
            
            # Type decay, so test twice
            f.save('f.casadi')
            self.checkfunction_light(f,f_ref,inputs=test_point,digits=7)
            f = None
            f = Function.load("f.casadi")
            self.checkfunction(f,f_ref,inputs=test_point,hessian=False,digits=7,evals=1)
            
if __name__ == '__main__':
    unittest.main()
