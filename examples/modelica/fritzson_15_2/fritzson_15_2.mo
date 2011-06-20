model BasicVolume "Conservation of Mass"
  import Modelica.SIunits.*;
  parameter SpecificHeatCapacity R = 287;
  Pressure P;
  Volume V;
  Mass m(start=0.00119);
  Temperature T;
  MassFlowRate mdot_in;
  MassFlowRate mdot_out;
equation
  // Boundary equations
  V=1e-3;
  T=293;
  mdot_in=0.1e-3;
  mdot_out=0.01e-3;
  // Conservation of mass
  der(m)=mdot_in-mdot_out;
  // Equation of state (ideal gas)
  P*V=m*R*T;
end BasicVolume;
