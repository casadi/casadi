model rocket
  Real s(start=0, unit="m");
  Real v(start=0, unit="m/s");
  Real m(start=1, unit="kg");
  parameter Real a = 0.05;
  parameter Real b = 0.1;
  input Real u(start = 0.4, min = -0.5, max = 0.5);
equation
  der(s) = v;
  der(v) = (u - a*v^2)/m;
  der(m) = -b*u^2;
end rocket;
