model vdp
  Real x1(start=0, fixed = true);
  Real x2(start=1, fixed = true);
  input Real u(min = -1, max = 1);
equation
  der(x1) = (1-x2^2)*x1 - x2 + u;
  der(x2) = x1;
end vdp;
