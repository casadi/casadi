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

model BasicVolume2
  import Modelica.SIunits.*;
  parameter SpecificInternalEnergy u_0 = 209058;
  parameter SpecificHeatCapacity c_v = 717;
  parameter Temperature T_0 = 293;
  parameter Mass m_0 = 0.00119;
  parameter SpecificHeatCapacity R = 287;
  Pressure P;
  Volume V;
  Mass m(start=m_0);
  Temperature T;
  MassFlowRate mdot_in;
  MassFlowRate mdot_out;
  SpecificEnthalpy h_in, h_out;
  SpecificEnthalpy h;
  Enthalpy H;
  SpecificInternalEnergy u;
  InternalEnergy U(start=u_0*m_0);
equation

  // Boundary equations
  V=1e-3;
  T=293;
  mdot_in=0.1e-3;
  mdot_out=0.01e-3;
  h_in = 300190;
  h_out = h;

  // Conservation of mass
  der(m) = mdot_in-mdot_out;

  // Conservation of energy
  der(U) = h_in*mdot_in - h_out*mdot_out;

  // Specific internal energy (ideal gas)
  u = U/m;
  u = u_0+c_v*(T-T_0);

  // Specific enthalpy
  H = U+P*V;
  h = H/m;

  // Equation of state (ideal gas)
  P*V=m*R*T;
end BasicVolume2;

