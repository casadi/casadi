from casadi import *
import casadi as c
from numpy import *
import unittest
from types import *
from helpers import *
from time import time
import sys
from scipy import linalg, std, mean
from scipy.stats import t

student = t

class ComplexityTests(casadiTestCase):
  maxt = 4   # [s] Aim for execution times up to maxt
  mint = 0.5 # [s] Only trust execution times larger than mint
  rejectat = 0.01  # [-] Reject the null hypothesis if p-score is smaller than rejectat
  testorders = [0,1,2,3] # Test these orders
  
  def checkOrders(self,Ns,ts):
    """
    Test the hypothesis that the slope of Ns is order.
    
    
    An order -N means that the order is between (N-1) and N
    
    """
    # Test statistic http://stattrek.com/AP-Statistics-4/Test-Slope.aspx?Tutorial=Stat
    # http://mathworld.wolfram.com/LeastSquaresFitting.html
    # http://en.wikipedia.org/wiki/Student%27s_t-test#Slope_of_a_regression_line

    orders   = self.testorders
    rejectat = self.rejectat
    
    Ns = log(Ns)
    ts = log(ts)
    
    m = Ns.shape[0]

    sigmaN = std(Ns)
    sigmaT = std(ts)
        
    # covariance
    cov = mean(ts*Ns)-mean(ts)*mean(Ns)
    
    # correlation coefficient
    rho = (mean(ts*Ns)-mean(ts)*mean(Ns))/sigmaN/sigmaT
    
    # sqrt(Variance) estimate 
    s = m/sqrt(m-2)*sqrt(sigmaT**2 - cov**2/sigmaN**2)
    
    # Estimate for b
    b = sigmaT/sigmaN * rho
    
    a = mean(ts) - b * mean(Ns)
    
    # Standard error on estimate for b
    sigmab = s/sqrt(m*sigmaN**2)
    
    # Standard deviation
    sigmaa = s * sqrt(1/m + mean(ts)**2/(m*sigmaT**2))
    
    conf = 0.05
    
    print "O(f) = %.3e N^(%.3f) [s]" % (exp(a),b)
     
    results = []
    for order in orders:
      # Null hypothesis: the slope of the regression is equal to order
      # test statistic
      t = abs(b-order)/sigmab

      p = (1-student.cdf(t,m))
      print "If the order were really %d, then the measurements have a p-score of %.8f" % (order,p)
      if p >= rejectat:
        results.append(order)

    
    if len(results)==1:
      order = results[0]
      a = mean(ts - order*Ns)
      sigmaa = std(ts - order*Ns)/sqrt(m)
      print "O(f) = %.3e N^%d [s]    | 95%% confidence:   [%.3e , %.3e] N^%d [s]" % (exp(a),order,exp(student.ppf(conf, m, loc=a, scale=sigmaa)),exp(student.ppf(1-conf, m, loc=a, scale=sigmaa)),order)
    
    
    for i, order in zip(range(len(orders))[1:],orders[1:]):
      if b < order and b > order-1 and not(order in results) and not(order -1 in results):
        results.append(-order)
    
    return results
      
  
  def complexity(self,setupfun,fun, order, depth = 0):
    
    N = 1
    
    Ns = []
    ts = []
    
    dt = 0
    while dt<self.maxt:
      sys.stdout.write("\r%d..." % N)
      sys.stdout.flush()
      p = setupfun(self,N) # Setup
      t = time()
      fun(self,N,p) # Run the function
      dt = time()-t
      Ns.append(N)
      ts.append(dt)
      N=int(N*1.1)
      N+=1
    print ""

    Ns = array(Ns)
    ts = array(ts)
    valid = ts > self.mint
    orders = self.checkOrders(Ns[valid],ts[valid])
    if len(orders)!=1 or orders[0]!=order:
      if (depth<3):
        return self.complexity(setupfun,fun, order, depth+1 )
      else:
        self.assertTrue(False,"We expected order %d, but found %s" % (order,str(orders)))
    
    
  def test_DMatrixadd(self):
    self.message("DMatrix add column vectors")
    def setupfun(self,N):
      return {'A': DMatrix(N,1,0), 'B': DMatrix(N,1,0)}
    def fun(self,N,setup):
      setup['A'] + setup['B']
    
    self.complexity(setupfun,fun, 1)
    
    self.message("DMatrix add rows vectors")
    def setupfun(self,N):
      return {'A': DMatrix(1,N,0), 'B': DMatrix(1,N,0)}
    
    self.complexity(setupfun,fun, 1)
    
    
if __name__ == '__main__':
    unittest.main()

