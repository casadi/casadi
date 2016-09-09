import casadi.*

x = SX.sym('x')
y = SX.sym('y',2)
z = SX.sym('z',2,2)
a = SX.sym('a',Sparsity.upper(2))


f = Function('f',{x,y,z,a},{x,y,z,a})
f.generate('fmex',struct('mex',true))

mex fmex.c -largeArrayDims
Anum = 1;
Bnum = [2;3];
Cnum = [4 5; 6 7];
Dnum = sparse([8 9; 0 10]);
[A,B,C,D] = fmex('f',1,[2;3],[4 5; 6 7],sparse([8 9; 0 10]))

assert(all(A==Anum))
assert(all(B==Bnum))
assert(all(all(C==Cnum)))
assert(all(all(D==Dnum)))
