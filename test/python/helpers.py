#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
import casadi
from numpy import *
import unittest
import sys
from math import isnan, isinf
import itertools

import argparse
import struct

platform_arch = 8 * struct.calcsize("P")

parser = argparse.ArgumentParser()
parser.add_argument('--known_bugs', help='Run with known bugs', action='store_true')
parser.add_argument('--ignore_memory_heavy', help='Skip those tests that have a high memory footprint', action='store_true')
parser.add_argument('--ignore_memory_light', help='Skip those tests that have a lightweight memory footprint', action='store_true')
parser.add_argument('unittest_args', nargs='*')

args = parser.parse_args()

import sys
sys.argv[1:] = ['-v'] + args.unittest_args

from StringIO import StringIO

class TeeString(StringIO):
  def __init__(self,stream):
    StringIO.__init__(self)
    self.stream = stream
        
  def write(self, data):
    StringIO.write(self,data)
    self.stream.write(data)

class Stdout():
    def __init__(self):
        self.stream = TeeString(sys.stdout)

    def __enter__(self):
        sys.stdout = self.stream
        return self.stream

    def __exit__(self, type, value, traceback):
        sys.stdout = self.stream.stream
        
class Stderr():
    def __init__(self):
        self.stream = TeeString(sys.stderr)

    def __enter__(self):
        sys.stderr = self.stream
        return self.stream

    def __exit__(self, type, value, traceback):
        sys.stderr = self.stream.stream
        
        
class FunctionPool:
  def __init__(self):
    self.numpyoperators=[]
    self.casadioperators=[]
    self.names=[]
    self.flags=[]
  def append(self,cas,num,name="",flags=set()):
    self.casadioperators.append(cas)
    self.numpyoperators.append(num)
    self.names.append(name)
    self.flags.append(flags)

class casadiTestCase(unittest.TestCase):

  def __init__(self,*margs,**kwargs):
    fun = getattr(getattr(self,margs[0]),'im_func')
    if not hasattr(fun,'tag_memory_heavy'):
      fun.tag_memory_heavy = False
    
    
    if args.ignore_memory_heavy and fun.tag_memory_heavy:
      fun.__unittest_skip__ = True
      fun.__unittest_skip_why__ = "Ignoring memory_heavy tests (--ignore_memory_heavy)"
    if args.ignore_memory_light and not(fun.tag_memory_heavy):
      fun.__unittest_skip__ = True
      fun.__unittest_skip_why__ = "Ignoring memory_light tests (--ignore_memory_light)"
      
    unittest.TestCase.__init__(self,*margs,**kwargs)

  def randDMatrix(self,n,m=1,sparsity=1,valuegenerator=lambda : random.normal(0,1),symm=False ):
    if sparsity < 1:
      spp = self.randDMatrix(n,m,sparsity=1,valuegenerator=lambda : random.uniform(0,1) ,symm=symm)
      spm = (spp < sparsity)
      makeSparse(spm)
      ret = DMatrix(spm.sparsity(),[valuegenerator() for i in range(spm.size())])
      if symm:
        return (ret + ret.T)/2
      else:
        return ret
    else:
      ret = DMatrix([valuegenerator() for i in range(n*m)],n,m)
      if symm:
        return (ret + ret.T)/2
      else:
        return ret
  
  def message(self,s):
      print s
      sys.stdout.flush()

  def assertAlmostEqual(self,first, second, places=7, msg=""):
      msg+= " %.16e <-> %.16e"  % (first, second)
      n =  max(abs(first),abs(second))
      if n>1e3:
        n = 10**floor(log10(n))
      else:
        n = 1.0
      unittest.TestCase.assertAlmostEqual(self,first/n,second/n,places=places,msg=msg + "  scaled by %d" %n)

  def checkarray(self,zr,zt,name="",failmessage="",digits=10):
      """
      Checks for equality of two numpy matrices.
      The check uses dense form.
      
      zr - reference
      zt - to test
      
      name - a descriptor that will be included in error messages
      
      """
      if isinstance(zr,tuple):
        zr=array([list(zr)])
      if isinstance(zt,tuple):
        zt=array([list(zt)])
      if isinstance(zr,list):
        zr=array([zr])
      if isinstance(zt,list):
        zt=array([zt])
      if not(hasattr(zt,'shape')) or len(zt.shape)==0:
        zt=array([[zt]])
      if not(hasattr(zr,'shape')) or len(zr.shape)==0:
        zr=array([[zr]])
      if len(zr.shape)==1:
        zr=array(zr,ndmin=2).T
      if len(zt.shape)==1:
        zt=array(zt,ndmin=2).T
        

      self.assertEqual(zt.shape[0],zr.shape[0],"In %s: %s dimension error. Got %s, expected %s. %s <-> %s" % (name,failmessage,str(zt.shape),str(zr.shape),str(zt),str(zr)))
      self.assertEqual(len(zt.shape),len(zr.shape),"In %s: %s dimension error. Got %s, expected %s. %s <-> %s" % (name,failmessage,str(zt.shape),str(zr.shape),str(zt),str(zr)))
      self.assertEqual(zt.shape[1],zr.shape[1],"In %s: %s dimension error. Got %s, expected %s. %s <-> %s" % (name,failmessage,str(zt.shape),str(zr.shape),str(zt),str(zr)))
      for i in range(zr.shape[0]):
        for j in range(zr.shape[1]):
          try:
            float(zt[i,j])
            float(zr[i,j])
          except:
            self.assertTrue(isEqual(zt[i,j],zr[i,j]),"Expressions (%s,%s) are not equal \n %s <-> \n %s at elem(%d,%d): %s <-> %s" % (type(zt),type(zr),str(zt),str(zr),i,j,str(zt[i,j]),str(zr[i,j])))
            continue
          if zt[i,j]==zr[i,j]:
            continue
          #if (isnan(zt[i,j]) or isinf(zt[i,j])) and  (isinf(zt[i,j]) or isnan(zt[i,j])):
          #  continue
          self.assertAlmostEqual(zt[i,j],zr[i,j],digits,"In %s: %s evaluation error.\n %s <->\n %s\n [digits=%d] at elem(%d,%d): " % (name,failmessage,str(zt),str(zr), digits, i,j, ))

  def evaluationCheck(self,yt,yr,x,x0,name="",failmessage="",fmod=None,setx0=None):
    """ General unit test for checking casadi evaluation against a reference solution.
    
        Checks if yr is the same as S/MXFunction(x,yt) , evaluated for x0
    
        x: the symbolic seed for the test. It should be of a form accepted as first argument of SXFunction/MXFunction
        yt: the test expression.
        yr: the reference solution: a numpy matrix.
        
        name - a descriptor that will be included in error messages
    """
    self.message(":"+ name)
    if (type(x)==list):
      sample = x[0]
    else :
      sample = x
      

    if isinstance(sample,SX) or isinstance(sample,SXMatrix):
      f = SXFunction(x,yt)
    else:
      f = MXFunction(x,yt)
    
    f.init()
    if (not(fmod is None)):
      f=fmod(f,x)
    if not(type(x0)==list):
      x0=[x0]
    
    if setx0 is None:
      setx0=x0
      
    if not(type(setx0)==list):
      setx0=[setx0]
      
    for i in range(len(x0)):
      try:
        f.input(i).set(setx0[i])
      except Exception as e:
         print f.input(i).shape
         raise e
         raise Exception("ERROR! Tried to set input with %s which is of type  %s \n%s" %(str(x0[i]), str(type(x0[i])),name))
    f.evaluate()
    zt = f.output(0).toArray()
    self.checkarray(yr,zt,name,failmessage)
    
  def numpyEvaluationCheck(self,ft,fr,x,x0,name="",failmessage="",fmod=None,setx0=None):
    """ General unit test for checking casadi versus numpy evaluation.
    
        Checks if 'fr(x0)' yields the same as S/MXFunction(x,[ft(x)]) , evaluated for x0
    
        x: the symbolic seed for the test. It should be of a form accepted as first argument of SXFunction/MXFunction
        ft: the test function. This function should operate on the casadi matrix x and return MX or SXMatrix.
        fr: the reference function. This function works on the numpy array x0.
        
        name - a descriptor that will be included in error messages
    """
    self.message(":"+ name)
    fx=None
    frx=None
    try:
      fx=ft(x)
      frx=fr(x0)
    except Exception as e:
      print "Error calling functions in %s" % name
      raise e
    self.evaluationCheck([fx],frx,x,x0,name,failmessage,fmod=fmod,setx0=setx0)

  
  def numpyEvaluationCheckPool(self,pool,x,x0,name="",fmod=None,setx0=None,excludeflags=set()):
    """ Performs a numpyEvaluationCheck for all members of a function pool"""
    for i in range(len(pool.numpyoperators)):
      if len(excludeflags.intersection(pool.flags[i]))>0:
        continue
      self.numpyEvaluationCheck(pool.casadioperators[i],pool.numpyoperators[i],x,x0,"%s:%s" % (name,pool.names[i]),"\n I tried to apply %s (%s) from test case '%s' to numerical value %s. But the result returned: " % (str(pool.casadioperators[i]),pool.names[i],name, str(x0)),fmod=fmod,setx0=setx0)

  def checkfx(self,trial,solution,fwd=True,adj=True,jacobian=True,gradient=True,hessian=True,sens_der=True,evals=True,digits=9,digits_sens=None,failmessage="",allow_empty=True,verbose=True,indirect=False):

    if indirect:
      ins = trial.symbolicInput()
      extra_trial = MXFunction(ins,trial.call(ins))
      extra_trial.init()
      for i in range(trial.getNumInputs()):
        extra_trial.setInput(trial.input(i),i)
      self.checkfx(extra_trial,solution,fwd,adj,jacobian,gradient,hessian,sens_der,evals,digits=digits,digits_sens=digits_sens,failmessage=failmessage,allow_empty=allow_empty,verbose=verbose,indirect=False)
      for i in range(trial.getNumInputs()):
        trial.setInput(extra_trial.input(i),i)

    if digits_sens is None:
      digits_sens = digits
     
    for i in range(trial.getNumInputs()):
      if (allow_empty and (trial.input(i).empty() or solution.input(i).empty() )): continue
      message = "input(%d)" % i
      if verbose: print message + ": " + str(trial.getInput(i))
      self.checkarray(trial.getInput(i),solution.getInput(i),"",digits=digits,failmessage=failmessage+": "+ message)

    trial_inputs    = [ DMatrix(trial.getInput(k)) for k in range(trial.getNumInputs())]
    solution_inputs = [ DMatrix(solution.getInput(k)) for k in range(solution.getNumInputs())] 

    try:
      trial.evaluate(0,0)
      solution.evaluate(0,0)
    except Exception as e:
      raise Exception(str(e) + "\nThis occured for simple evaluate(%d,%d) for: %s" % (0,0,failmessage) )

    for i in range(trial.getNumInputs()):
      message = "input(%d) modified by evaluate" % i
      self.checkarray(trial.getInput(i),trial_inputs[i],"",digits=digits,failmessage=failmessage+": "+ message)
      self.checkarray(solution.getInput(i),solution_inputs[i],"",digits=digits,failmessage=failmessage+": "+ message)

    self.assertEqual(trial.getNumOutputs(),solution.getNumOutputs(),failmessage+": trial has %d outputs while solution has %d." % (trial.getNumOutputs(),solution.getNumOutputs()) )
    self.assertEqual(trial.getNumInputs(),solution.getNumInputs(),failmessage+": trial has %d inputs while solution has %d." % (trial.getNumInputs(),solution.getNumInputs()) )

    for i in range(trial.getNumOutputs()):
      message = "output(%d)" % i
      if verbose: print message + ": " + str(trial.getOutput(i))
      if (allow_empty and (trial.output(i).empty() or solution.output(i).empty() )): continue
      self.checkarray(trial.getOutput(i),solution.getOutput(i),"",digits=digits,failmessage=failmessage+": "+message)
      
    try:
      trial.evaluate(0,0)
      solution.evaluate(0,0)
    except Exception as e:
      raise Exception(str(e) + "\nThis occured for repeated evaluate(%d,%d) for: %s" % (0,0,failmessage) )

    for i in range(trial.getNumInputs()):
      message = "input(%d) modified by repeated evaluate" % i
      self.checkarray(trial.getInput(i),trial_inputs[i],"",digits=digits,failmessage=failmessage+": "+ message)
      self.checkarray(solution.getInput(i),solution_inputs[i],"",digits=digits,failmessage=failmessage+": "+ message)

    self.assertEqual(trial.getNumOutputs(),solution.getNumOutputs(),failmessage+": trial has %d outputs while solution has %d." % (trial.getNumOutputs(),solution.getNumOutputs()) )
    self.assertEqual(trial.getNumInputs(),solution.getNumInputs(),failmessage+": trial has %d inputs while solution has %d." % (trial.getNumInputs(),solution.getNumInputs()) )

    for i in range(trial.getNumOutputs()):
      message = "output(%d)" % i
      if verbose: print message + ": " + str(trial.getOutput(i))
      if (allow_empty and (trial.output(i).empty() or solution.output(i).empty() )): continue
      self.checkarray(trial.getOutput(i),solution.getOutput(i),"",digits=digits,failmessage=failmessage+": "+message)
    
    if fwd:
      fsm = 1.7
      for i in range(trial.getNumInputs()):
        for j in range(min(trial.input(i).size(),solution.input(i).size())):
          trial.setFwdSeed(0,i)
          solution.setFwdSeed(0,i)
          trial.fwdSeed(i)[j]=fsm
          solution.fwdSeed(i)[j]=fsm
          
          try:
            trial.evaluate(fwd,0)
            solution.evaluate(fwd,0)
          except Exception as e:
            raise Exception(str(e) + "\nThis occured for simple evaluate(%d,%d) for fwdSeed(%d)[%d]=1 for: %s" % (fwd,0,i,j,failmessage) )
          
          for k in range(trial.getNumOutputs()):
            if (allow_empty and (trial.output(k).empty() or solution.output(k).empty() )): continue
            message="fwdSeed(%d)[%d]=1 => fwdSens(%d)" % (i,j,k)
            if verbose: print message + ": " + str(trial.getFwdSens(k))
            self.checkarray(trial.getFwdSens(k),solution.getFwdSens(k),"",digits=digits_sens,failmessage=failmessage+": "+message)
            
          trial.setFwdSeed(0,i)
          solution.setFwdSeed(0,i)
        
    if adj:
      asm = 1.7
      for i in range(trial.getNumOutputs()):
        for j in range(min(trial.output(i).size(),solution.output(i).size())):
          trial.setAdjSeed(0,i)
          solution.setAdjSeed(0,i)
          trial.adjSeed(i)[j]=asm
          solution.adjSeed(i)[j]=asm
          
          try:
            trial.evaluate(0,adj)
            solution.evaluate(0,adj)
          except Exception as e:
            raise Exception(str(e) + "\nThis occured for simple evaluate(%d,%d) for adjSeed(%d)[%d]=1 for: %s" % (0,adj,i,j,failmessage) )

          for k in range(trial.getNumInputs()):
            if (allow_empty and (trial.input(k).empty() or solution.input(k).empty() )): continue
            message="adjSeed(%d)[%d]=1 => adjSens(%d)" % (i,j,k)
            if verbose: print message + ": " + str(trial.getAdjSens(k))
            self.checkarray(trial.getAdjSens(k),solution.getAdjSens(k),"",digits=digits_sens,failmessage=failmessage+": "+message)
            
          trial.setAdjSeed(0,i)
          solution.setAdjSeed(0,i)

    for i in range(trial.getNumInputs()):
      message = "input(%d) modified by evaluate" % i
      self.checkarray(trial.getInput(i),trial_inputs[i],"",digits=digits,failmessage=failmessage+": "+ message)
      self.checkarray(solution.getInput(i),solution_inputs[i],"",digits=digits,failmessage=failmessage+": "+ message)
    
    if jacobian:
      for i in range(trial.getNumInputs()):
        if (allow_empty and (trial.input(i).empty() or solution.input(i).empty() )): continue
        for j in range(trial.getNumOutputs()):
          trialjac = trial.jacobian(i,j)
          trialjac.init()
          self.assertEqual(trialjac.getNumInputs(),trial.getNumInputs())
          self.assertEqual(trialjac.getNumOutputs(),trial.getNumOutputs()+1)
          for k in range(trial.getNumInputs()): trialjac.setInput(trial_inputs[k],k)
          solutionjac = solution.jacobian(i,j)
          solutionjac.init()
          self.assertEqual(solutionjac.getNumInputs(),solution.getNumInputs())
          self.assertEqual(solutionjac.getNumOutputs(),solution.getNumOutputs()+1)
          for k in range(solution.getNumInputs()): solutionjac.setInput(solution_inputs[k],k)
          
          self.checkfx(trialjac,solutionjac,fwd=fwd if sens_der else False,adj=adj if sens_der else False,jacobian=False,gradient=False,hessian=False,evals=False,digits=digits_sens,failmessage="(%s).jacobian(%d,%d)" % (failmessage,i,j),allow_empty=allow_empty,verbose=verbose)

    if gradient:
      for i in range(trial.getNumInputs()):
        if (allow_empty and (trial.input(i).empty() or solution.input(i).empty() )): continue
        for j in range(trial.getNumOutputs()):
          if trial.output(j).scalar() and solution.output(j).scalar():
            trialgrad = trial.gradient(i,j)
            trialgrad.init()
            self.assertEqual(trialgrad.getNumInputs(),trial.getNumInputs())
            self.assertEqual(trialgrad.getNumOutputs(),trial.getNumOutputs()+1)
            for k in range(trial.getNumInputs()): trialgrad.setInput(trial_inputs[k],k)
            solutiongrad = solution.gradient(i,j)
            solutiongrad.init()
            self.assertEqual(solutiongrad.getNumInputs(),solution.getNumInputs())
            self.assertEqual(solutiongrad.getNumOutputs(),solution.getNumOutputs()+1)
            for k in range(solution.getNumInputs()): solutiongrad.setInput(solution_inputs[k],k)
            self.checkfx(trialgrad,solutiongrad,fwd=fwd  if sens_der else False,adj=adj if sens_der else False,jacobian=False,gradient=False,hessian=False,evals=False,digits=digits_sens,failmessage="(%s).gradient(%d,%d)" % (failmessage,i,j),allow_empty=allow_empty,verbose=verbose)

    if hessian:
      for i in range(trial.getNumInputs()):
        if (allow_empty and (trial.input(i).empty() or solution.input(i).empty() )): continue
        for j in range(trial.getNumOutputs()):
          if trial.output(j).scalar() and solution.output(j).scalar():
            trialhess = trial.hessian(i,j)
            trialhess.init()
            self.assertEqual(trialhess.getNumInputs(),trial.getNumInputs())
            self.assertEqual(trialhess.getNumOutputs(),trial.getNumOutputs()+2)
            for k in range(trial.getNumInputs()): trialhess.setInput(trial_inputs[k],k)
            solutionhess = solution.hessian(i,j)
            solutionhess.init()
            self.assertEqual(solutionhess.getNumInputs(),solution.getNumInputs())
            self.assertEqual(solutionhess.getNumOutputs(),solution.getNumOutputs()+2)
            for k in range(solution.getNumInputs()): solutionhess.setInput(solution_inputs[k],k)
            self.checkfx(trialhess,solutionhess,fwd=fwd  if sens_der else False,adj=adj  if sens_der else False,jacobian=False,gradient=False,hessian=False,evals=False,digits=digits_sens,failmessage="(%s).hessian(%d,%d)" % (failmessage,i,j),allow_empty=allow_empty,verbose=verbose)     

    for k in range(trial.getNumInputs()):
      trial.setInput(trial_inputs[k],k)
      solution.setInput(solution_inputs[k],k)
      
    if evals is True:
      evals = 2

    if evals:

      def remove00(x):
        ret = DMatrix(x)
        ret[0,0] = DMatrix.sparse(1,1)
        return ret
        
      spmods = [lambda x: x , remove00]
      spmods = [lambda x: x]
      
      sym = msym
      Function = MXFunction
      
      storage2 = {}
      storage = {}
      
      ndir = 2
      trial.requestNumSens(ndir,ndir)
      solution.requestNumSens(ndir,ndir)
      
      def flatten(l):
        ret = []
        for i in l:
          ret.extend(i)
        return ret
        
      vf_reference = None
        
      for f,values in [(solution,solution_inputs),(trial,trial_inputs)]:
      
        # dense
        for spmod,spmod2 in itertools.product(spmods,repeat=2):
          fseeds = [[sym("f",spmod(f.getInput(i)).sparsity()) for i in range(f.getNumInputs())]  for d in range(ndir)]
          aseeds = [[sym("a",spmod2(f.getOutput(i)).sparsity())  for i in range(f.getNumOutputs())] for d in range(ndir)]
          inputss = [sym("i",f.input(i).sparsity()) for i in range(f.getNumInputs())]
      
          res,fwdsens,adjsens = f.eval(inputss,fseeds,aseeds)
          
          vf = Function(inputss+flatten([fseeds[i]+aseeds[i] for i in range(ndir)]),list(res) + flatten([list(fwdsens[i])+list(adjsens[i]) for i in range(ndir)]))
          
          vf.init()

          for i,v in enumerate(values):
            vf.setInput(v,i)
        
          # Complete random seeding
          random.seed(1)
          for i in range(vf.getNumInputs()):
            vf.setInput(DMatrix(vf.input(i).sparsity(),random.random(vf.input(i).size())),i)
          
          vf.evaluate()
          storagekey = (spmod,spmod2)
          if not(storagekey in storage):
            storage[storagekey] = []
          storage[storagekey].append([vf.getOutput(i) for i in range(vf.getNumOutputs())])
          
          if vf_reference is None:
            vf_reference = vf
            
          if evals>1:

            # Second order sensitivities
            for spmod_2,spmod2_2 in itertools.product(spmods,repeat=2):
              fseeds2 = [[sym("f",vf_reference.input(i).sparsity()) for i in range(vf.getNumInputs())] for d in range(ndir)]
              aseeds2 = [[sym("a",vf_reference.output(i).sparsity())  for i in range(vf.getNumOutputs()) ] for d in range(ndir)]
              inputss2 = [sym("i",vf_reference.input(i).sparsity()) for i in range(vf.getNumInputs())]
           
              res2,fwdsens2,adjsens2 = vf.eval(inputss2,fseeds2,aseeds2)

              vf2 = Function(inputss2+flatten([fseeds2[i]+aseeds2[i] for i in range(ndir)]),list(res2) + flatten([list(fwdsens2[i])+list(adjsens2[i]) for i in range(ndir)]))
              vf2.init()
                
              random.seed(1)
              for i in range(vf2.getNumInputs()):
                vf2.setInput(DMatrix(vf2.input(i).sparsity(),random.random(vf2.input(i).size())),i)
              
              vf2.evaluate()
              storagekey = (spmod,spmod2)
              if not(storagekey in storage2):
                storage2[storagekey] = []
              storage2[storagekey].append([vf2.getOutput(i) for i in range(vf2.getNumOutputs())])

      # Remainder of eval testing
      for store,order in [(storage,"first-order"),(storage2,"second-order")][:evals]:
        for stk,st in store.items():
          for i in range(len(st)-1):
            for k,(a,b) in enumerate(zip(st[0],st[i+1])):
              if b.numel()==0 and sparse(a).size()==0: continue
              if a.numel()==0 and sparse(b).size()==0: continue
              self.checkarray(IMatrix(a.sparsity(),1),IMatrix(b.sparsity(),1),("%s, output(%d)" % (order,k))+str(vf.getInput(0)),digits=digits_sens)
              self.checkarray(a,b,("%s, output(%d)" % (order,k))+str(vf.getInput(0)),digits=digits_sens)
              
    for k in range(trial.getNumInputs()):
      trial.setInput(trial_inputs[k],k)
      solution.setInput(solution_inputs[k],k)

      

class run_only(object):
  def __init__(self, args):
    self.args = []
    for a in args:
      self.args.append(a)

  def __call__(self, c):
    print "run_only:"
    for i in dir(c):
      if i.startswith('test_'):
        n = i[5:]
        if not n in self.args:
          delattr(c,i)
        else:
          print i
    return c

class requires(object):
  def __init__(self,att):
    self.att = att
  
  def __call__(self,c):
    if hasattr(casadi,self.att):
      return c
    else:
      print "Not available %s, skipping unittests" % self.att
      return None

class skip(object):
  def __init__(self, skip=True):
    self.skip = skip

  def __call__(self, c):
    if not self.skip: return c
    if isinstance(c,unittest.TestCase):
      if i.startswith('test_'):
        delattr(c,i)
      return c
    else:
      print self.skiptext(c.__name__)
      return None
   
  def skiptext(self,name):
    return "Skipping test '%s'" % name

def known_bug():
  return unittest.skipIf(not(args.known_bugs),"known bug (run with --known_bugs to include it)")
    
class memory_heavy(object):
  def __init__(self):
    pass
    
  def __call__(self, c):
    print c
    c.tag_memory_heavy = True
    return c
