from numpy import *
from casadi import *

DMatrix.setPrecision(16)

solver = None
  
def sdqp_sol(h=None,c=None,a=None,uba=None,f=None,g=None):
  global solver
  if solver is None:
    solver = SDPSDQPSolver(sdqpStruct(h=h.sparsity(),f=f.sparsity(),g=g.sparsity(),a=a.sparsity()))
    solver.setOption("sdp_solver",DSDPSolver)
    solver.setOption("sdp_solver_options",{"_printlevel": 0})
    solver.init()
    
  solver.setInput(h,"h")
  solver.setInput(f,"f")
  solver.setInput(g,"g")
  solver.setInput(a,"a")
  solver.setInput(c,"c")
  solver.setInput(-Inf,"lba")
  solver.setInput(uba,"uba")
  solver.evaluate()
  
  return solver.output("x"),solver.output("lam_a"), solver.output("dual")

x = ssym("x",2)

f = (1-x[0])**2+100*(x[1]-x[0]**2)**2
nsd = blockcat([[-x[0],2],[2,-x[1]**2]])  # <=0

g = eig_symbolic(nsd)

nlp = SXFunction(nlpIn(x=x),nlpOut(f=f,g=g))
nlp.init()

# Find a refence solution with another
ipopt = IpoptSolver(nlp)
ipopt.init()
ipopt.setInput(-Inf,"lbg")
ipopt.setInput(0,"ubg")
ipopt.solve()


print "reference sol= ", ipopt.output("x")

g = DMatrix(0,1)

lambd = ssym("lambda",g.shape)
Lambd = ssym("lambda",nsd.sparsity())

lag = f+mul(lambd.T,g)+trace(mul(Lambd,nsd))

oracle = SXFunction(customIO(x=x,lambd=lambd,Lambd=Lambd),customIO(f=f,g=g,nsd=nsd,hess=hessian(lag,x), gradF=gradient(f,x), jacG= jacobian(g,x),jac_nsd=jacobian(vec(nsd),x)))
oracle.init()

lambda_k = DMatrix([0])
Lambda_k = DMatrix([0])
x_k = DMatrix([2,3])

for i in range(25):
  print i, x_k
  oracle.setInput(x_k,"x")
  oracle.setInput(lambda_k,"lambd")
  oracle.setInput(Lambda_k,"Lambd")
  oracle.evaluate()
  
  step, lambda_k, Lambda_k = sdqp_sol(h=oracle.output("hess"),c=oracle.output("gradF"),a=oracle.output("jacG"),uba=-oracle.output("g"),f=vertcat([ oracle.output("jac_nsd")[:,i].reshape(oracle.output("nsd").shape) for i in range(x_k.size())]),g=-oracle.output("nsd"))
   
  x_k+= step
  
print linalg.eig(oracle.output("nsd"))[0]
