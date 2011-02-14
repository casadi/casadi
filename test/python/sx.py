from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *

class SXtests(casadiTestCase):

  def setUp(self):
    self.pool=FunctionPool()
    self.pool.append(lambda x: sqrt(x[0]),sqrt,"sqrt")
    self.pool.append(lambda x: sin(x[0]),sin,"sin")
    self.pool.append(lambda x: cos(x[0]),cos,"cos")
    self.pool.append(lambda x: tan(x[0]),tan,"tan")
    self.pool.append(lambda x: arctan(x[0]),arctan,"arctan")
    self.pool.append(lambda x: arcsin(x[0]),arcsin,"arcsin")
    self.pool.append(lambda x: arccos(x[0]),arccos,"arccos")
    self.pool.append(lambda x: exp(x[0]),exp,"exp")
    self.pool.append(lambda x: log(x[0]),log,"log")
    self.pool.append(lambda x: x[0]**0,lambda x : x**0,"x^0")
    self.pool.append(lambda x: x[0]**1,lambda x : x**1,"^1")
    self.pool.append(lambda x: x[0]**(-2),lambda x : x**(-2),"^-2")
    self.pool.append(lambda x: x[0]**(0.3),lambda x : x**(0.3),"^0.3")
    self.pool.append(lambda x: floor(x[0]),floor,"floor")
    self.pool.append(lambda x: ceil(x[0]),ceil,"ceil")
    self.Jpool=FunctionPool()
    self.Jpool.append(lambda x: sqrt(x[0]),lambda x:diag(1/(2.0*sqrt(x))),"sqrt")
    self.Jpool.append(lambda x: sin(x[0]),lambda x:diag(cos(x)),"sin")
    self.Jpool.append(lambda x: cos(x[0]),lambda x:diag(-sin(x)),"cos")
    self.Jpool.append(lambda x: tan(x[0]),lambda x:diag(1.0/cos(x)**2),"tan")
    self.Jpool.append(lambda x: arctan(x[0]),lambda x:diag( 1.0/(x**2+1)),"arctan")
    self.Jpool.append(lambda x: arcsin(x[0]),lambda x:diag( 1.0/sqrt(1-x**2)),"arcsin")
    self.Jpool.append(lambda x: arccos(x[0]),lambda x: diag(-1.0/sqrt(1-x**2)),"arccos")
    self.Jpool.append(lambda x: exp(x[0]),lambda x: diag(exp(x)),"exp")
    self.Jpool.append(lambda x: log(x[0]),lambda x: diag(1.0/x),"log")
    self.Jpool.append(lambda x: x[0]**0,lambda x :diag(zeros(x.shape)),"x^0")
    self.Jpool.append(lambda x: x[0]**1,lambda x : diag(ones(x.shape)),"^1")
    self.Jpool.append(lambda x: x[0]**(-2),lambda x : diag(-2.0/x**3),"^-2")
    self.Jpool.append(lambda x: x[0]**(0.3),lambda x :diag( 0.3/x**0.7),"^0.3")
    self.matrixpool=FunctionPool()
    self.matrixpool.append(lambda x: norm_2(x[0]),linalg.norm,"norm_2")
    self.matrixbinarypool=FunctionPool()
    self.matrixbinarypool.append(lambda a: a[0]+a[1],lambda a: a[0]+a[1],"Matrix+Matrix")
    self.matrixbinarypool.append(lambda a: a[0]-a[1],lambda a: a[0]-a[1],"Matrix-Matrix")
    self.matrixbinarypool.append(lambda a: a[0]*a[1],lambda a: a[0]*a[1],"Matrix*Matrix")
    #self.matrixbinarypool.append(lambda a: inner_prod(a[0],trans(a[1])),lambda a: dot(a[0].T,a[1]),name="inner_prod(Matrix,Matrix)") 
    self.matrixbinarypool.append(lambda a: c.dot(a[0],trans(a[1])),lambda a: dot(a[0],a[1].T),"dot(Matrix,Matrix.T)")

    #self.pool.append(lambda x: erf(x[0]),erf,"erf") # numpy has no erf
    
    
  
  def test_scalarSX(self):
      x=symbolic("x")
      x0=0.738
      
      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarSX")
      
  def test_gradient(self):
      self.message("jacobian of SX**number")
      x=symbolic("x");
      x0=1;
      p=3 # increase to 20 to showcase ticket #56
      y=x**p;
      dx=jacobian(y,x);
      dxr=p;
      self.evaluationCheck([dx],dxr,[x],x0,name="jacobian");
      dxr=1
      for i in list(range(p)):
        y=jacobian(y,x)
        dxr=dxr*(p-i)
      
      
      self.evaluationCheck([y],dxr,[x],x0,name="recursive jacobian");

  def test_gradient2(self):
      self.message("jacobian of SX**SX")
      x=symbolic("x");
      p=symbolic("p");
      x0=1;
      p0=3 # increase to 20 to showcase ticket #56
      y=x**p;
      dx=jacobian(y,x);
      #print dx
      dxr=p0;
      self.evaluationCheck([dx],dxr,[x,p],[x0,p0],name="jacobian");

      dxr=1
      for i in list(range(p0)):
        y=jacobian(y,x)
        dxr=dxr*(p0-i)
      
      self.evaluationCheck([y],dxr,[x,p],[x0,p0],name="jacobian");
      
  def test_SXMatrixJacobian(self):
      self.message("SXMatrix(1,1) unary operation, jacobian")
      x=symbolic("x")
      x0=array([[0.738]])

      def fmod(f,x):
        J=f.jacobian()
        J.init()
        return J
      
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SXMatrix unary operations, jacobian",fmod=fmod)
      
  def test_SXMatrixJac(self):
      self.message("SXMatrix(1,1) unary operation, jac")
      x=symbolic("x")
      x0=array([[0.738]])

      def fmod(f,x):
        j=f.jac()
        J=SXFunction(x,[j])
        J.init()
        return J
      
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SXMatrix unary operations, jac",fmod=fmod)
      
  def test_SXMatrixJacobians(self):
      self.message("SXMatrix(3,1) unary operation, jacobian")
      x=symbolic("x",3)
      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        J=f.jacobian()
        J.init()
        return J
      
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SXMatrix unary operations, jacobian",fmod=fmod)
      
  def test_SXMatrix(self):
      self.message("SXMatrix unary operations")
      x=symbolic("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      
      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="SXMatrix")
      
      x=symbolic("x",3,3)
      x0=array([[0.738,0.2,0.3],[ 0.1,0.39,-6 ],[0.99,0.999999,-12]])
      #self.numpyEvaluationCheck(lambda x: c.det(x[0]), lambda   x: linalg.det(x),[x],x0,name="det(SXMatrix)")
      self.numpyEvaluationCheck(lambda x: SXMatrix([c.det(x[0])]), lambda   x: linalg.det(x),[x],x0,name="det(SXMatrix)")
      self.numpyEvaluationCheck(lambda x: c.inv(x[0]), lambda   x: linalg.inv(x),[x],x0,name="inv(SXMatrix)")
        
  def test_SXMatrixSparse(self):
      self.message("SXMatrix unary operations, sparse")
      from scipy.sparse import csr_matrix
      x=SX("x")
      y=SX("y")
      z=SX("z")
      x=SXMatrix(3,4,[1,2,1],[0,2,2,3],[x,y,z])
      x0=DMatrix(3,4,[1,2,1],[0,2,2,3],[0.738,0.1,0.99]).toCsr_matrix()
      
      self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="SXMatrix",setx0=x0)
      
  def test_SXMatrixbinary(self):
      self.message("SXMatrix binary operations")
      x=symbolic("x",3,2)
      y=symbolic("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="SXMatrix")
      self.assertRaises(RuntimeError, lambda : c.dot(x,y))
      
  def test_SXMatrixslicing(self):
      self.message("SXMatrix slicing")
      x=symbolic("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])

      self.numpyEvaluationCheck(lambda x: SXMatrix(x[0][0,0]), lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]")
      self.numpyEvaluationCheck(lambda x: SXMatrix(x[0][1,0]), lambda x: matrix(x)[1,0],[x],x0,name="x[1,0]")
      self.numpyEvaluationCheck(lambda x: SXMatrix(x[0][0,1]), lambda x: matrix(x)[0,1],[x],x0,name="x[1,0]")
      self.numpyEvaluationCheck(lambda x: x[0][:,0], lambda x: matrix(x)[:,0],[x],x0,name="x[:,0]")
      self.numpyEvaluationCheck(lambda x: x[0][:,1], lambda x: matrix(x)[:,1],[x],x0,name="x[:,1]")
      self.numpyEvaluationCheck(lambda x: x[0][1,:], lambda x: matrix(x)[1,:],[x],x0,name="x[1,:]")
      self.numpyEvaluationCheck(lambda x: x[0][0,:], lambda x: matrix(x)[0,:],[x],x0,name="x[0,:]")
      self.numpyEvaluationCheck(lambda x: x[0][-1,:], lambda x: matrix(x)[-1,:],[x],x0,name="x[-1,:]")
      self.numpyEvaluationCheck(lambda x: x[0][:,-2], lambda x: matrix(x)[:,-2],[x],x0,name="x[:,-2]")
      self.numpyEvaluationCheck(lambda x: x[0][0:-2,0:-1], lambda x: matrix(x)[0:-2,0:-1],[x],x0,name="x[0:-2,0:-1]")
      self.numpyEvaluationCheck(lambda x: x[0][0:2,0:2], lambda x: matrix(x)[0:2,0:2],[x],x0,name="x[0:2,0:2]")
  
  def test_SX1(self):
    self.message("SXFunction evaluation")
    fun=lambda x,y: [x+y,x*y,x**2+y**3]
    x=SX("x")
    y=SX("y")
    f=SXFunction([[x,y]],[fun(x,y)])
    f.init()
    L=[2,3]
    f.setInput(L)
    f.evaluate()
    z=f.output(0).toArray()
    zr=fun(*L)
    for i in range(3):
      self.assertAlmostEqual(z[i], zr[i],10,'SXfunction output in correct')
    
    J=f.jacobian()
    J.init()
    J.setInput(L)
    J.evaluate()
    Jr=matrix([[1,1],[3,2],[4,27]])
    self.checkarray(J.output(0),Jr,"SXfunction jacobian evaluates incorrectly")
          
  def test_SX2(self):
    self.message("SXFunction evalution + forward and adjoint seeds")
    fun = lambda x,y: [3-sin(x*x)-y, sqrt(y)*x]
    # variables
    x = SX("x")
    y = SX("y")

    # Create function
    f = fun(x,y)
    self.assertEqual(str(f),'[((3-sin((x*x)))-y), (sqrt(y)*x)]','SX representation is wrong')
    fcn = SXFunction([[x,y]],[f])

    # Set some options
    fcn.setOption("name","f")

    self.assertEqual(repr(fcn),'f','SX representation is wrong')
    # Initialize the function for numerical calculations
    fcn.init()

    # Pass inputs
    L=[2,3]
    fcn.setInput(L)

    # Pass forward seed
    sF=[1,0]
    fcn.setFwdSeed(sF);

    sA=[1,0]
    # Pass adjoint seed
    fcn.setAdjSeed(sA)

    # Evaluate numerically
    fcn.evaluate(1,1)

    # Get the results
    res = tuple(fcn.output())
    self.assertAlmostEqual(res[0], fun(*L)[0],10,'SXfunction evaluation wrong')
    self.assertAlmostEqual(res[1], fun(*L)[1],10,'SXfunction evaluation wrong')

    fsens = tuple(fcn.fwdSens())
    e=1e-8
    p=fun((L[0]+e*sF[0]),(L[1]+e*sF[1]))
    fsensn=[(p[0]-res[0])/e,(p[1]-res[1])/e]
    self.assertAlmostEqual(fsensn[0], fsens[0],6,'SXfunction forward mode evaluation wrong')
    self.assertAlmostEqual(fsensn[1], fsens[1],6,'SXfunction forward mode evaluation wrong')
    
    J=fcn.jacobian()
    J.init()
    J.setInput(L)
    J.evaluate()
    J=J.output(0).toArray()
    
    fsensJ=dot(J,array(sF))
    self.assertAlmostEqual(fsensJ[0], fsens[0],10,'SXfunction forward mode evaluation wrong')
    self.assertAlmostEqual(fsensJ[1], fsens[1],10,'SXfunction forward mode evaluation wrong')
    
    asens = tuple(fcn.adjSens())
    asensJ=dot(J.T,array(sA))
    self.assertAlmostEqual(asensJ[0], asens[0],10,'SXfunction adjoint mode evaluation wrong')
    self.assertAlmostEqual(asensJ[1], asens[1],10,'SXfunction adjoint mode evaluation wrong')

    # evaluate symbolically
    fsubst = fcn.eval([[1-y,1-x]])
    #print "fsubst = ", fsubst
    
  def test_SXFunctionc(self):
    self.message("SXFunction constructors")
    x=SX("x")
    y=symbolic("y",2,3)
    f=SXFunction(([x,x,x],[x,x,x,x]),[[x]])
    self.checkarray(f.input(0).shape,(3,1),"SXFunction constructors")
    self.checkarray(f.input(1).shape,(4,1),"SXFunction constructors")
    self.checkarray(f.output(0).shape,(1,1),"SXFunction constructors")
    
    f=SXFunction(((x,x),(x,x)),[(x,x)])
    self.checkarray(f.input(0).shape,(2,1),"SXFunction constructors")
    self.checkarray(f.input(1).shape,(2,1),"SXFunction constructors")
    self.checkarray(f.output(0).shape,(2,1),"SXFunction constructors")
    
    f=SXFunction([y],[y])
    self.checkarray(f.input(0).shape,(2,3),"SXFunction constructors")
    self.checkarray(f.output(0).shape,(2,3),"SXFunction constructors")
    
    f=SXFunction([y,[x,x]],[y,y])
    self.checkarray(f.input(0).shape,(2,3),"SXFunction constructors")
    self.checkarray(f.input(1).shape,(2,1),"SXFunction constructors")
    self.checkarray(f.output(0).shape,(2,3),"SXFunction constructors")
    self.checkarray(f.output(1).shape,(2,3),"SXFunction constructors")

    self.assertRaises(NotImplementedError,lambda: SXFunction(y,[y,y]))
    self.assertRaises(NotImplementedError,lambda: SXFunction(x,[x,x]))

  def test_evalfail(self):
    self.message("eval fail test")
    x = symbolic("x",2,2)
    f = SXFunction([x], [x])
    self.assertRaises(TypeError,lambda: f.eval(x))

  def test_SXconversion(self):
    self.message("Conversions from and to SXMatrix")
    y=SX("y")
    x=symbolic("x",3,3)
    SXMatrix(y)
    SXMatrix(x)
    c.det(x)
    y=array(x)
    c.det(y)

  def test_SXineq(self):
    self.message("Test (in)equality operators")
    
    list= {"x>y": lambda x,y: x>y,
              "x<y": lambda x,y: x<y,
              "x<=y": lambda x,y: x<=y,
              "x>=y": lambda x,y: x>=y,
              "x==y": lambda x,y: x==y,
              "x!=y": lambda x,y: x!=y
    }
    
    
    group = { "x SXMatrix, y SXMatrix": (symbolic("x"),symbolic("y")),
                    "x SX, y SXMatrix": (SX("x"),symbolic("y")),
                    "x SXMatrix, y SX": (SX("x"),SX("y")),
                    "x SX, y SX": (SX("x"),SX("y"))
                  }
    for gname, e in group.items():
      self.message(":"+ gname)
      for name, op in list.items():
        self.message("::" + name)
        self.assertTrue(isinstance(op(*e),SXMatrix),gname + "/" + name)
        
  def test_SXFunctionc2(self):
    self.message("SXmatrix typemaps constructors")
    #simplify(SX("x"))                 
    isEmpty(array([[SX("x")]]))
    list = [ ("SX" ,SX("x"),(1,1)),
                ("number",2.3, (1,1)),
                ("list(SX)", [SX("x"),SX("y")], (2,1)),
                ("list(SX,number)", [SX("x"),2.3], (2,1) ),
                ("tuple(SX)", (SX("x"),SX("y")), (2,1)),
                ("tuple(SX,number)", (SX("x"),2.3), (2,1)),
                ("SXMatrix", symbolic("x"), (1,1)),
                ("numpy.ndarray1D(SX)", array([SX("x"),SX("y")]), (2,1)),
                ("numpy.ndarray(SX)", array([[SX("x"),SX("y")],[SX("w"),SX("z")]]), (2,2)),
                ("numpy.ndarray(SX,number)", array([[SX("x"),2.3]]), (1,2))
    ];
    for name, arg,shape in list:
      self.message(":" + name)
      i=trans(trans(arg))
      self.assertEqual(i.shape[0],shape[0],"shape mismatch")
      self.assertEqual(i.shape[1],shape[1],"shape mismatch")
      isEmpty(arg)
    
  def test_SXFunctionc3(self):
    self.message("vector(SXmatrix) typemaps constructors")
    y=SX("y")
    x=symbolic("x",3,1)
    vertcat([x,x])
    vertcat([y,y])
    #vertcat([x,[y]])
    
  def test_eval(self):
    self.message("SXFunction eval")
    x=symbolic("x",2,2)
    y=symbolic("y",2,2)
    f  = SXFunction([x,y], [x*y])
    f.eval([x,y])


    
if __name__ == '__main__':
    unittest.main()

