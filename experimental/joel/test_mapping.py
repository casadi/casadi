from casadi import *
#row = [0,0,1,1,2]
#col = [0,0,1,1,2]
#row = [2,1,1,0,0]
#col = [2,0,0,0,0]
row = [0,1,2]
col = [0,0,0]


mapping = IVector()
print "row = ", row
print "col = ", col

sp = sp_triplet(60,1,row,col,mapping,False,True)
print "mapping = ", mapping
print "sp = ", sp
print "sp.col() = ", sp.col()
print "sp.getRow() = ", sp.getRow()

a = DMatrix(10,1)
a[1,0] = 1
nz = IVector([0,1,2])
print "nz = ", nz
a.sparsity().getNZInplace(nz)
print "nz = ", nz
