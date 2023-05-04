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
import casadi
import numpy
from numpy import random, array
import unittest
import sys
from math import isnan, isinf
import itertools
import time
from contextlib import contextmanager
from casadi.tools import capture_stdout

import argparse
import struct

if sys.version_info >= (3, 0):
  import builtins
else:
  import __builtin__
  builtins = __builtin__

x = MX.sym("x",2,3)
y = MX.sym("x",3,4)

try:
  x+y
except TypeError:
  swig4 = True
except RuntimeError:
  swig4 = False

try:
  nlpsol(123)
except Exception as e:
  if "SXDict" in str(e):
    systemswig = True
  elif "dict:SX" in str(e):
    systemswig = False
  else:
    systemswig = True

# Pending deprecation in numpy
check_matrix = False

platform_arch = 8 * struct.calcsize("P")

parser = argparse.ArgumentParser()
parser.add_argument('--known_bugs', help='Run with known bugs', action='store_true')
parser.add_argument('--ignore_memory_heavy', help='Skip those tests that have a high memory footprint', action='store_true')
parser.add_argument('--ignore_memory_light', help='Skip those tests that have a lightweight memory footprint', action='store_true')
parser.add_argument('--run_slow', help='Skip those tests that take a long time to run', action='store_true')
parser.add_argument('unittest_args', nargs='*')

args = parser.parse_args()

import sys
sys.argv[1:] = ['-v'] + args.unittest_args

from io import StringIO


class LazyString(object):
  def __init__(self,f):
     self.f = f
     self.context = dict(sys._getframe(1).f_locals)

  def __add__(self, other):
     return LazyString('str(self) + other')
  def __radd__(self, other):
     return LazyString('other + str(self)')

  def __str__(self):
    d = self.context
    exec("ret = " + self.f.replace("\n","\\n"), d)
    return str(d["ret"])

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
  def zip(self):
    return list(zip(self.casadioperators,self.numpyoperators,self.names,self.flags))


def toSX_fun(fun):
  ins = fun.sx_in()
  return Function("f",ins,fun.call(ins))

def toMX_fun(fun):
  ins = fun.mx_in()
  return Function("f",ins,fun(ins))

def jacobian_old(f, i, j):
    return f.factory(f.name() + '_jac', f.name_in(),
        ['jac:' + f.name_out(j) + ':' + f.name_in(i)] + f.name_out())

def hessian_old(f, i, j):
    return f.factory(f.name() + '_hess', f.name_in(),
        ['hess:' + f.name_out(j) + ':' + f.name_in(i) + ":" + f.name_in(i),
        'grad:' + f.name_out(j) + ':' + f.name_in(i)] + f.name_out())

class casadiTestCase(unittest.TestCase):

  @classmethod
  def tearDownClass(cls):
    print("STATUS_RAN_ALL_TESTS")

  @contextmanager
  def assertInAnyOutput(self,s):
    e = ""
    with capture_stdout() as result:
        try:
          yield
        except Exception as err:
          e = str(err)
    print(result[0])
    self.assertTrue(s in e or s in result[0] or s in result[1])

  @contextmanager
  def assertInException(self,s):
    e = None
    try:
      yield
    except Exception as err:
      e = str(err)
    self.assertFalse(e is None)
    self.assertTrue(s in e,msg=e + "<->" + s)

  @contextmanager
  def assertOutput(self,included,excluded):
    with capture_stdout() as result:
      yield
    if not(isinstance(included,list)):
      included = [included]
    if not(isinstance(excluded,list)):
      excluded = [excluded]
    for e in included:
      self.assertTrue(e in result[0],msg=result[0] + "<->" + e)
    for e in excluded:
      self.assertFalse(e in result[0],msg=result[0] + "<->" + e)
  def tearDown(self):
    t = time.time() - self.startTime
    print("deltaT %s: %.3f" % ( self.id(), t))

  def __init__(self,*margs,**kwargs):
    self.startTime = time.time()
    if sys.version_info >= (3, 0):
      fun = getattr(self,margs[0])
      if not hasattr(fun,'tag_memory_heavy'):
        fun.__dict__["tag_memory_heavy"] = False
      if not hasattr(fun,'tag_slow'):
        fun.__dict__["tag_slow"] = False

      if args.ignore_memory_heavy and fun.tag_memory_heavy:
        fun.__dict__["__unittest_skip__"] = True
        fun.__dict__["__unittest_skip_why__"] = "Ignoring memory_heavy tests (--ignore_memory_heavy)"
      if args.ignore_memory_light and not(fun.tag_memory_heavy):
        fun.__dict__["__unittest_skip__"] = True
        fun.__dict__["__unittest_skip_why__"] = "Ignoring memory_light tests (--ignore_memory_light)"

      if not(args.run_slow) and fun.tag_slow:
        fun.__dict__["__unittest_skip__"] = True
        fun.__dict__["__unittest_skip_why__"] = "Ignoring slow tests (--run_slow)"

    else:
      fun = getattr(getattr(self,margs[0]),'im_func')
      if not hasattr(fun,'tag_memory_heavy'):
        fun.tag_memory_heavy = False
      if not hasattr(fun,'tag_slow'):
        fun.tag_slow = False

      if args.ignore_memory_heavy and fun.tag_memory_heavy:
        fun.__unittest_skip__ = True
        fun.__unittest_skip_why__ = "Ignoring memory_heavy tests (--ignore_memory_heavy)"
      if args.ignore_memory_light and not(fun.tag_memory_heavy):
        fun.__unittest_skip__ = True
        fun.__unittest_skip_why__ = "Ignoring memory_light tests (--ignore_memory_light)"

      if not(args.run_slow) and fun.tag_slow:
        fun.__unittest_skip__ = True
        fun.__unittest_skip_why__ = "Ignoring slow tests (--run_slow)"

    unittest.TestCase.__init__(self,*margs,**kwargs)

  def randDM(self,n,m=1,sparsity=1,valuegenerator=lambda : random.normal(0,1),symm=False ):
    if sparsity < 1:
      spp = self.randDM(n,m,sparsity=1,valuegenerator=lambda : random.uniform(0,1) ,symm=symm)
      spm = (spp < sparsity)
      spm = sparsify(spm)
      ret = DM(spm.sparsity(),array([valuegenerator() for i in range(spm.nnz())]))
      if symm:
        return (ret + ret.T)/2
      else:
        return ret
    else:
      ret = DM([valuegenerator() for i in range(n*m)]).reshape((n,m))
      if symm:
        return (ret + ret.T)/2
      else:
        return ret

  def message(self,s):
      print(s)
      sys.stdout.flush()

  def assertAlmostEqual(self,first, second, places=7, msg=""):
      if isnan(first ) and isnan(second): return
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
            self.assertTrue(is_equal(zt[i,j],zr[i,j]),LazyString('"Expressions (%s,%s) are not equal \n %s <-> \n %s at elem(%d,%d): %s <-> %s" % (type(zt),type(zr),str(zt),str(zr),i,j,str(zt[i,j]),str(zr[i,j]))'))
            #self.assertTrue(is_equal(zt[i,j],zr[i,j]),"Expressions (%s,%s) are not equal \n %s <-> \n %s at elem(%d,%d): %s <-> %s" % (type(zt),type(zr),str(zt),str(zr),i,j,str(zt[i,j]),str(zr[i,j])))
            continue
          if zt[i,j]==zr[i,j]:
            continue
          #if (isnan(zt[i,j]) or isinf(zt[i,j])) and  (isinf(zt[i,j]) or isnan(zt[i,j])):
          #  continue
          self.assertAlmostEqual(zt[i,j],zr[i,j],digits,LazyString('"In %s: %s evaluation error.\n %s <->\n %s\n [digits=%d] at elem(%d,%d): " % (name,failmessage,str(zt),str(zr), digits, i,j, )'))
          #self.assertAlmostEqual(zt[i,j],zr[i,j],digits,"In %s: %s evaluation error.\n %s <->\n %s\n [digits=%d] at elem(%d,%d): " % (name,failmessage,str(zt),str(zr), digits, i,j, ))

  def evaluationCheck(self,yt,yr,x,x0,name="",failmessage="",fmod=None,setx0=None):
    """ General unit test for checking casadi evaluation against a reference solution.

        Checks if yr is the same as Function(name,x,yt) , evaluated for x0

        x: the symbolic seed for the test. It should be of a form accepted as first argument of Function
        yt: the test expression.
        yr: the reference solution: a numpy matrix.

        name - a descriptor that will be included in error messages
    """
    self.message(":"+ name)
    if (type(x)==list):
      sample = x[0]
    else :
      sample = x


    if isinstance(sample,SX):
      f = Function("f", x, yt)
    else:
      f = Function("f", x, yt)

    if (not(fmod is None)):
      f=fmod(f,x)
    if not(type(x0)==list):
      x0=[x0]

    if setx0 is None:
      setx0=x0

    if not(type(setx0)==list):
      setx0=[setx0]

    f_in = [0]*f.n_in()
    for i in range(len(x0)):
      try:
        f_in[i]=setx0[i]
      except Exception as e:
         print(f.size_in(i))
         raise e
         raise Exception("ERROR! Tried to set input with %s which is of type  %s \n%s" %(str(x0[i]), str(type(x0[i])),name))
    f_out = f.call(f_in)
    zt = f_out[0].full()
    self.checkarray(yr,zt,name,failmessage)

  def numpyEvaluationCheck(self,ft,fr,x,x0,name="",failmessage="",fmod=None,setx0=None):
    """ General unit test for checking casadi versus numpy evaluation.

        Checks if 'fr(x0)' yields the same as Function(name,x,[ft(x)]) , evaluated for x0

        x: the symbolic seed for the test. It should be of a form accepted as first argument of Function
        ft: the test function. This function should operate on the casadi matrix x and return MX or SX.
        fr: the reference function. This function works on the array x0.

        name - a descriptor that will be included in error messages
    """
    self.message(":"+ name)
    function=None
    frx=None
    try:
      function=ft(x)
      frx=fr(x0)
    except Exception as e:
      print("Error calling functions in %s" % name)
      raise e
    self.evaluationCheck([function],frx,x,x0,name,failmessage,fmod=fmod,setx0=setx0)


  def numpyEvaluationCheckPool(self,pool,x,x0,name="",fmod=None,setx0=None,excludeflags=set()):
    """ Performs a numpyEvaluationCheck for all members of a function pool"""
    for i in range(len(pool.numpyoperators)):
      if len(excludeflags.intersection(pool.flags[i]))>0:
        continue
      self.numpyEvaluationCheck(pool.casadioperators[i],pool.numpyoperators[i],x,x0,"%s:%s" % (name,pool.names[i]),"\n I tried to apply %s (%s) from test case '%s' to numerical value %s. But the result returned: " % (str(pool.casadioperators[i]),pool.names[i],name, str(x0)),fmod=fmod,setx0=setx0)
      
  def check_eval_mx(self, mx):
    if isinstance(mx,list):
        return self.check_eval_mx(vvcat(mx))
    assert isinstance(mx,MX)
    args = symvar(mx)
    f = Function("f",args,[mx])
    d = MX.sym("d") # dummy
    #mx_eval = f.call(args,True,False)[0]
    new_args = [MX.sym(e.name(),e.sparsity()) for e in args]
    mx_eval = f.call(new_args,True,False)[0]
    if isinstance(mx_eval,DM):
        mx = evalf(mx)
        self.checkarray(mx,mx_eval)
        return
    mx_eval = substitute([mx_eval],new_args,args)[0]
    if isinstance(mx_eval,MX): # In all sparse case, might falsely look like SX
        self.check_identical_mx(mx_eval,mx)

  def check_identical_fun(self, a, b):
    s = a.serialize()
    s2 = b.serialize()
    a.save('a.casadi',{"debug":True})
    b.save('b.casadi',{"debug":True})
    self.assertTrue(s==s2)
        
  def check_identical_mx(self, a, b):
    assert isinstance(a,MX)
    assert isinstance(b,MX)
    arg_a = symvar(a)
    arg_b = symvar(b)
    fa = Function('f',arg_a,[a])
    fb = Function('f',arg_b,[b])
    self.check_identical_fun(fa,fb)    

  def checkfunction_light(self,trial,solution,inputs=None,**kwargs):
    self.checkfunction(trial,solution,inputs,fwd=False,adj=False,jacobian=False,gradient=False,hessian=False,sens_der=False,evals=False,**kwargs)
  def checkfunction(self,trial,solution,inputs=None,fwd=True,adj=True,jacobian=True,gradient=True,hessian=True,sens_der=True,evals=True,digits=9,digits_sens=None,failmessage="",allow_empty=True,verbose=True,indirect=False,sparsity_mod=True,allow_nondiff=False):

    if isinstance(inputs,dict):
      d = inputs
      inputs = [0]*trial.n_in()
      ns = trial.name_in()
      for k,v in list(d.items()):
        inputs[ns.index(k)] = v

    if indirect:
      ins = trial.mx_in()
      extra_trial = Function("extra_trial", ins,trial(ins))
      self.checkfunction(extra_trial,solution,fwd,adj,jacobian,gradient,hessian,sens_der,evals,inputs=inputs,digits=digits,digits_sens=digits_sens,failmessage=failmessage,allow_empty=allow_empty,verbose=verbose,indirect=False)

    if digits_sens is None:
      digits_sens = digits

    for i in range(trial.n_in()):
      if (allow_empty and (trial.sparsity_in(i).is_empty() or solution.sparsity_in(i).is_empty() )): continue
      message = "input(%d: '%s')" % (i, trial.name_in(i))

    for i in range(2): # repeated evaluation
      try:
        trial_outputs = trial.call(inputs)
        solution_outputs = solution.call(inputs)
      except Exception as e:
        raise Exception(str(e) + "\nThis occured for simple evaluate(%d,%d) for: %s" % (0,0,failmessage) )

      self.assertEqual(trial.n_out(),solution.n_out(),failmessage+": trial has %d outputs while solution has %d." % (trial.n_out(),solution.n_out()) )
      self.assertEqual(trial.n_in(),solution.n_in(),failmessage+": trial has %d inputs while solution has %d." % (trial.n_in(),solution.n_in()) )

      for i in range(trial.n_out()):
        message = "output(%d: '%s')" % (i, trial.name_out(i))
        if (allow_empty and (trial.sparsity_out(i).is_empty() or solution.sparsity_out(i).is_empty() )): continue
        if (allow_nondiff and (trial.sparsity_out(i).nnz()==0 or solution.sparsity_out(i).nnz()==0 )): continue
        self.checkarray(trial_outputs[i],solution_outputs[i],"",digits=digits,failmessage=failmessage+": "+message)

    if jacobian:
      for i in range(trial.n_in()):
        if (allow_empty and (trial.sparsity_in(i).is_empty() or solution.sparsity_in(i).is_empty() )): continue
        for j in range(trial.n_out()):
          trialjac = jacobian_old(trial, i, j)
          self.assertEqual(trialjac.n_in(),trial.n_in())
          self.assertEqual(trialjac.n_out(),trial.n_out()+1)
          solutionjac = jacobian_old(solution, i, j)
          self.assertEqual(solutionjac.n_in(),solution.n_in())
          self.assertEqual(solutionjac.n_out(),solution.n_out()+1)

          self.checkfunction(trialjac,solutionjac,inputs=inputs,fwd=fwd if sens_der else False,adj=adj if sens_der else False,jacobian=False,gradient=False,hessian=False,evals=False,digits=digits_sens,failmessage="(%s).jacobian_old(%d,%d)" % (failmessage,i,j),allow_empty=allow_empty,verbose=verbose,allow_nondiff=allow_nondiff)

    if hessian:
      for i in range(trial.n_in()):
        if (allow_empty and (trial.sparsity_in(i).is_empty() or solution.sparsity_in(i).is_empty() )): continue
        for j in range(trial.n_out()):
          if trial.sparsity_out(j).is_scalar() and solution.sparsity_out(j).is_scalar():
            trialhess = hessian_old(trial, i, j)
            self.assertEqual(trialhess.n_in(),trial.n_in())
            self.assertEqual(trialhess.n_out(),trial.n_out()+2)
            solutionhess = hessian_old(solution, i, j)
            self.assertEqual(solutionhess.n_in(),solution.n_in())
            self.assertEqual(solutionhess.n_out(),solution.n_out()+2)
            self.checkfunction(trialhess,solutionhess,inputs=inputs,fwd=fwd  if sens_der else False,adj=adj  if sens_der else False,jacobian=False,gradient=False,hessian=False,evals=False,digits=digits_sens,failmessage="(%s).hessian_old(%d,%d)" % (failmessage,i,j),allow_empty=allow_empty,verbose=verbose,allow_nondiff=allow_nondiff)

    if evals is True:
      evals = 2

    if evals:

      def remove_first(x):
        ret = DM(x)
        if ret.numel()>0:
          ret[0,0] = DM(1,1)
          return ret.sparsity()
        else:
          return ret.sparsity()

      def remove_last(x):
        ret = DM(x)
        if ret.nnz()>0:
          ret[ret.sparsity().row()[-1],ret.sparsity().get_col()[-1]] = DM(1,1)
          return ret.sparsity()
        else:
          return ret.sparsity()

      spmods = [lambda x: x , remove_first, remove_last]
      #spmods = [lambda x: x]
      #spmods = [lambda x: x , remove_first]

      sym = MX.sym

      storage2 = {}
      storage = {}

      def vec(l):
        ret = []
        for i in l:
          ret.extend(i)
        return ret

      vf_reference = None

      for f in [solution,trial]:

        # dense
        for spmod,spmod2 in itertools.product(spmods,repeat=2):
          fseeds = [sym("f",spmod(f.sparsity_in(i))) for i in range(f.n_in())]
          aseeds = [sym("a",spmod2(f.sparsity_out(i)))  for i in range(f.n_out())]
          inputss = [sym("i",f.sparsity_in(i)) for i in range(f.n_in())]
          res = f.call(inputss,True)
          #print res, "sp", [i.sparsity().dim(True) for i in fseeds]
          opts = {"helper_options": {"is_diff_in": f.is_diff_in(), "is_diff_out": f.is_diff_out()}}
          [fwdsens] = forward(res, inputss, [fseeds],opts)
          [adjsens] = reverse(res, inputss, [aseeds],opts)

          vf = Function("vf", inputss+vec([fseeds+aseeds]),list(res) + vec([list(fwdsens)+list(adjsens)]),{"is_diff_in": f.is_diff_in()+f.is_diff_in()+f.is_diff_out(), "is_diff_out": f.is_diff_out()+f.is_diff_out()+f.is_diff_in()})

          vf_in = list(inputs)
          # Complete random seeding
          for i in range(f.n_in(),vf.n_in()):
            random.seed(i)
            vf_in.append(DM(vf.sparsity_in(i),random.random(vf.nnz_in(i))))

          vf_out = vf.call(vf_in)
          storagekey = (spmod,spmod2)
          if not(storagekey in storage):
            storage[storagekey] = []
          storage[storagekey].append(vf_out)

          if vf_reference is None:
            vf_reference = vf

          if evals>1:

            # Second order sensitivities
            for spmod_2,spmod2_2 in itertools.product(spmods,repeat=2):
              fseeds2 = [sym("f",vf_reference.sparsity_in(i)) for i in range(vf.n_in())]
              aseeds2 = [sym("a",vf_reference.sparsity_out(i))  for i in range(vf.n_out()) ]
              inputss2 = [sym("i",vf_reference.sparsity_in(i)) for i in range(vf.n_in())]

              res2 = vf.call(inputss2)
              opts = {"helper_options": {"is_diff_in": vf.is_diff_in(), "is_diff_out": vf.is_diff_out()}}
              [fwdsens2] = forward(res2, inputss2, [fseeds2],opts)
              [adjsens2] = reverse(res2, inputss2, [aseeds2],opts)

              vf2 = Function("vf2", inputss2+vec([fseeds2+aseeds2]),list(res2) + vec([list(fwdsens2)+list(adjsens2)]),{"is_diff_in": vf.is_diff_in()+vf.is_diff_in()+vf.is_diff_out(), "is_diff_out": vf.is_diff_out()+vf.is_diff_out()+vf.is_diff_in()})

              vf2_in = list(inputs)

              for i in range(f.n_in(),vf2.n_in()):
                random.seed(i)
                vf2_in.append(DM(vf2.sparsity_in(i),random.random(vf2.nnz_in(i))))

              vf2_out = vf2.call(vf2_in)
              storagekey = (spmod,spmod2)
              if not(storagekey in storage2):
                storage2[storagekey] = []
              storage2[storagekey].append(vf2_out)

      # Remainder of eval testing
      for store,order in [(storage,"first-order"),(storage2,"second-order")][:evals]:
        for stk,st in list(store.items()):
          for i in range(len(st)-1):
            for k,(a,b) in enumerate(zip(st[0],st[i+1])):
              if b.numel()==0 and sparsify(a).nnz()==0: continue
              if a.numel()==0 and sparsify(b).nnz()==0: continue

              if (allow_nondiff and (a.nnz()==0 or b.nnz()==0 )): continue
              #self.checkarray(IM(a.sparsity(),1),IM(b.sparsity(),1),("%s, output(%d)" % (order,k))+str(vf.getInput(0))+failmessage,digits=digits_sens)
              self.checkarray(a,b,("%s, output(%d)" % (order,k))+failmessage,digits=digits_sens)

  def check_sparsity(self, a,b):
    self.assertTrue(a==b, msg=str(a) + " <-> " + str(b))

  def check_codegen(self,F,inputs=None, opts=None,std="c89",extralibs="",check_serialize=False,extra_options=None,main=False,definitions=None,with_jac_sparsity=False,external_opts=None):

    if args.run_slow:
      import hashlib
      name = "codegen_%s" % (hashlib.md5(("%f" % np.random.random()+str(F)+str(time.time())).encode()).hexdigest())
      if opts is None: opts = {}
      if main: opts["main"] = True
      cg = CodeGenerator(name,opts)
      cg.add(F,with_jac_sparsity)
      cg.generate()
      import subprocess

      libdir = GlobalOptions.getCasadiPath()
      includedir = GlobalOptions.getCasadiIncludePath()
      includedirs = [includedir,os.path.join(includedir,"highs")]

      if isinstance(extralibs,list):
        extralibs_clean = []
        for lib in extralibs:
            if os.name=='nt':
                if "." in lib:
                    if lib.endswith(".dll"):
                        extralibs_clean.append(lib[:-4]+".lib")
                    else:
                        extralibs_clean.append(lib)
                else:
                    extralibs_clean.append(lib+".lib")
            else:
                if "." in lib:
                    extralibs_clean.append(lib)
                else:
                    extralibs_clean.append("-l"+lib)
        extralibs = " " + " ".join(extralibs_clean)

      if isinstance(extra_options,bool) or extra_options is None:
        extra_options = ""
      if isinstance(extra_options,list):
        extra_options = " " + " ".join(extra_options)
      if definitions is None:
        definitions = []

      def get_commands(shared=True):
        if os.name=='nt':
          defs = " ".join(["/D"+d for d in definitions])
          commands = "cl.exe {shared} {definitions} {includedir} {name}.c {extra} /link  /libpath:{libdir}".format(shared="/LD" if shared else "",std=std,name=name,libdir=libdir,includedir=" ".join(["/I" + e for e in includedirs]),extra=extralibs + extra_options + extralibs + extra_options,definitions=defs)
          if shared:
            output = "./" + name + ".dll"
          else:
            output = name + ".exe"
          return [commands, output]
        else:
          defs = " ".join(["-D"+d for d in definitions])
          output = "./" + name + (".so" if shared else "")
          commands = "gcc -pedantic -std={std} -fPIC {shared} -Wall -Werror -Wextra {includedir} -Wno-unknown-pragmas -Wno-long-long -Wno-unused-parameter -O3 {definitions} {name}.c -o {name_out} -L{libdir} -Wl,-rpath,{libdir} -Wl,-rpath,.".format(shared="-shared" if shared else "",std=std,name=name,name_out=name+(".so" if shared else ""),libdir=libdir,includedir=" ".join(["-I" + e for e in includedirs]),definitions=defs) + (" -lm" if not shared else "") + extralibs + extra_options
          if sys.platform=="darwin":
            commands+= " -Xlinker -rpath -Xlinker {libdir}".format(libdir=libdir)
            commands+= " -Xlinker -rpath -Xlinker .".format(libdir=libdir)
          return [commands, output]

      [commands, libname] = get_commands(shared=True)

      print("compile library",commands)
      p = subprocess.Popen(commands,shell=True).wait()
      if sys.platform=="darwin":
        subprocess.run(["otool","-l",libname])
      if external_opts is None: external_opts = {}
      F2 = external(F.name(), libname,external_opts)

      if main:
        [commands, exename] = get_commands(shared=False)
        print("here",commands)
        env = os.environ
        if os.name=='nt':
            env["PATH"] = env["PATH"]+";"+libdir
        else:
            env["LD_LIBRARY_PATH"] = libdir
        p = subprocess.Popen(commands,shell=True,env=env).wait()
        inputs_main = inputs
        if isinstance(inputs,dict):
          inputs_main = F.convert_in(inputs)
        F.generate_in(F.name()+"_in.txt", inputs_main)

      Fout = F.call(inputs)
      Fout2 = F2.call(inputs)

      if main:
        with open(F.name()+"_out.txt","w") as stdout:
          with open(F.name()+"_in.txt","r") as stdin:
            commands = exename+" "+F.name()
            print(commands+" < " + F.name()+"_in.txt")
            p = subprocess.Popen(commands,shell=True,stdin=stdin,stdout=stdout)
            out = p.communicate()
        assert p.returncode==0
        outputs = F.generate_out(F.name()+"_out.txt")
        print(outputs)
        if isinstance(inputs,dict):
          outputs = F.convert_out(outputs)
          for k in F.name_out():
            self.checkarray(Fout[k],outputs[k],digits=15)
        else:
          for i in range(F.n_out()):
            self.checkarray(Fout[i],Fout2[i],digits=15)

      if isinstance(inputs, dict):
        self.assertEqual(F.name_out(), F2.name_out())
        for k in F.name_out():
          self.checkarray(Fout[k],Fout2[k],digits=15)
      else:
        for i in range(F.n_out()):
          self.checkarray(Fout[i],Fout2[i],digits=15)

      if self.check_serialize:
        self.check_serialize(F2,inputs=inputs)
        
      return F2, libname

  def check_thread_safety(self,F,inputs=None,N=20):

    FP = F.map(N, 'thread',2)
    FS = F.map(N, 'thread')
    self.checkfunction_light(FP, FS, inputs)


  def check_serialize(self,F,inputs=None):
      F2 = Function.deserialize(F.serialize({"debug":True}))

      Fout = F.call(inputs)
      Fout2 = F2.call(inputs)

      if isinstance(inputs, dict):
        self.assertEqual(F.name_out(), F2.name_out())
        for k in F.name_out():
          self.checkarray(Fout[k],Fout2[k],digits=16)
      else:
        for i in range(F.n_out()):
          self.checkarray(Fout[i],Fout2[i],digits=16)


      if False:
          import hashlib
          h = hashlib.md5(F.serialize({"debug":True}).encode("ascii")).hexdigest()
          path = "serialize_"+casadi.__version__
          os.makedirs(path, exist_ok=True)
          F.save(os.path.join(path,h+".casadi"),{"debug":True})
          DM.set_precision(16)
          h_in = hashlib.md5(str(inputs).encode("ascii")).hexdigest()
          DM.set_precision(7)
          if isinstance(inputs,dict):
            inputs = F.convert_in(inputs)
            Fout = F.convert_out(Fout)
          F.generate_in(os.path.join(path,h+"_"+h_in+"_in.txt"),inputs)
          F.generate_out(os.path.join(path,h+"_"+ h_in +"_out.txt"),Fout)


  def check_pure(self,F,inputs=None):
      Fout = F.call(inputs)
      Fout2 = F.call(inputs)

      if isinstance(inputs, dict):
        for k in F.name_out():
          self.checkarray(Fout[k],Fout2[k],digits=16)
      else:
        for i in range(F.n_out()):
          self.checkarray(Fout[i],Fout2[i],digits=16)

class run_only(object):
  def __init__(self, args):
    self.args = []
    for a in args:
      self.args.append(a)

  def __call__(self, c):
    print("run_only:")
    for i in dir(c):
      if i.startswith('test_'):
        n = i[5:]
        if not n in self.args:
          delattr(c,i)
        else:
          print(i)
    return c

class requires(object):
  def __init__(self,att):
    self.att = att

  def __call__(self,c):
    if hasattr(casadi,self.att):
      return c
    else:
      print("Not available %s, skipping unittests" % self.att)
      return None

class requires_conic(object):
  def __init__(self,n):
    self.n = n

  def __call__(self,c):
    try:
      load_conic(self.n)
      return c
    except:
      print("Not available QP plugin %s, skipping unittests" % self.n)
      return None

class requires_linsol(object):
  def __init__(self,n):
    self.n = n

  def __call__(self,c):
    try:
      load_linsol(self.n)
      return c
    except:
      print("Not available linsol plugin %s, skipping unittests" % self.n)
      return None

class requires_nlpsol(object):
  def __init__(self,n):
    self.n = n

  def __call__(self,c):
    try:
      load_nlpsol(self.n)
      return c
    except:
      print("Not available NLP plugin %s, skipping unittests" % self.n)
      return None

class requires_expm(object):
  def __init__(self,n):
    self.n = n

  def __call__(self,c):
    try:
      load_expm(self.n)
      return c
    except:
      print("Not available Expm plugin %s, skipping unittests" % self.n)
      return None

class requires_integrator(object):
  def __init__(self,n):
    self.n = n

  def __call__(self,c):
    try:
      load_integrator(self.n)
      return c
    except:
      print("Not available integrator plugin %s, skipping unittests" % self.n)
      return None

class requires_rootfinder(object):
  def __init__(self,n):
    self.n = n

  def __call__(self,c):
    try:
      load_rootfinder(self.n)
      return c
    except:
      print("Not available RFP plugin %s, skipping unittests" % self.n)
      return None

class requiresPlugin(object):
  def __init__(self,att,n):
    self.att = att
    self.n = n

  def __call__(self,c):
    try:
      self.att.load_plugin(self.n)
      return c
    except:
      print("Not available %s plugin %s, skipping unittests" % (str(self.att),self.n))
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
      print(self.skiptext(c.__name__))
      return None

  def skiptext(self,name):
    return "Skipping test '%s'" % name

def known_bug():
  return unittest.skipIf(not(args.known_bugs),"known bug (run with --known_bugs to include it)")

class memory_heavy(object):
  def __init__(self):
    pass

  def __call__(self, c):
    print(c)
    c.__dict__["tag_memory_heavy"] = True
    return c

class slow(object):
  def __init__(self):
    pass

  def __call__(self, c):
    print("slow", c)
    c.__dict__["tag_slow"] = True
    return c
