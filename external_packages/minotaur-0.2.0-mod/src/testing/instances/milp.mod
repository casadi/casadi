var x0  binary;
var x1  binary;
var x2  binary;
var x3  binary;
var x4  binary;

minimize obj: x4;

subject to cons0: 2*x0 + 2*x1 + 2*x2 + 2*x3 + x4 = 1;
