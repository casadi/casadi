from casadi import *

# work-around
def fix_getitem(C):
  if not hasattr(C,"getitem_old"):
    C.getitem_old = C.__getitem__
    
    def getitem_new(self,s):
      if isinstance(s,tuple) and len(s)==2 and isinstance(s[0],slice) and isinstance(s[1],slice):
        I1 = s[0].indices(self.size1())
        I2 = s[1].indices(self.size2())
        
        J1 = []
        for i in range(I1[0],I1[1],I1[2]):
          J1.append(i)
        
        J2 = []
        for i in range(I2[0],I2[1],I2[2]):
          J2.append(i)
          
        return self.getitem_old((J1,J2))
      else:
        return self.getitem_old(s)
      
    C.__getitem__ = getitem_new

fix_getitem(SXMatrix)
fix_getitem(DMatrix)

# Test
A = symbolic("A",2,3)
A.resize(5,4)

print "A = ", A
print "A[0:2] = ", list(A[0:2])
print "A[0,2] = ", A[0,2]
print "A[:2,1:] = ", A[:2,1:]

B = DMatrix(4,5,1)
B[:] = range(1,21)
print "B[:,0:1] = ", B[:,0:1]
print "B[:,2:] = ", B[:,2:]
