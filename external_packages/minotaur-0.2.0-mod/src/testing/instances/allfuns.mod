var x{0..73};

minimize obj: x[0]; 

s.t. opacos:  x[0] + x[1]/(1+acos(x[2]+x[3])) = 0;
s.t. opacosh: x[0] + x[4]/(1+acosh(x[5]+x[6])) = 0;
s.t. opasin:  x[0] + x[7]/(1+asin(x[8]+x[9])) = 0;
s.t. opasinh: x[0] + x[10]/(1+asinh(x[11]+x[12])) = 0;
s.t. opatan:  x[0] + x[13]/(1+atan(x[14]+x[15])) = 0;
s.t. opatanh: x[0] + x[16]/(1+atanh(x[17]+x[18])) = 0;
s.t. opcos:   x[0] + x[19]/(1+cos(x[20]+x[21])) = 0;
s.t. opcosh:  x[0] + x[22]/(1+cosh(x[23]+x[24])) = 0;
s.t. opcpow:  x[0] + x[25]/(1+3.4^(x[26]+x[27])) = 0;
s.t. opdiv:   x[0] + x[28]/(1+x[29]/x[30]) = 0;
s.t. opexp:   x[0] + x[31]/(1+exp(x[32]+x[33])) = 0;
s.t. oplog:   x[0] + x[34]/(1+log(x[35]+x[36])) = 0;
s.t. oplog10: x[0] + x[37]/(1+log10(x[38]+x[39])) = 0;
s.t. opminus: x[0] + x[40]/(1-(x[41]+x[42])) = 0;
s.t. opmult:  x[0] + x[43]/(1+(x[44]*x[45])) = 0;
s.t. oppowK:  x[0] + x[46]/(1+(x[47]+x[48])^3.4) = 0;
s.t. opsin:   x[0] + x[49]/(1+sin(x[50]+x[51])) = 0;
s.t. opsinh:  x[0] + x[52]/(1+sinh(x[53]+x[54])) = 0;
s.t. opsqr:   x[0] + x[55]/(1+(x[56]+x[57])^2) = 0; 
s.t. opsqrt:  x[0] + x[58]/(1+sqrt(x[59]+x[60])) = 0;
s.t. opsum:   x[0] + x[61]*x[62] + x[63]*x[64] + x[65]*x[66]*x[67] = 0;
s.t. optan:   x[0] + x[68]/(1+tan(x[69]+x[70])) = 0;
s.t. optanh:  x[0] + x[71]/(1+tanh(x[72]+x[73])) = 0;

