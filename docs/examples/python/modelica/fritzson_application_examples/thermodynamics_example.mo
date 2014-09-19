model BasicVolumeMassConservation "Conservation of Mass"
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

end BasicVolumeMassConservation;

model BasicVolumeEnergyConservation
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
end BasicVolumeEnergyConservation;

model BasicVolumeHeatTransfer
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
  HeatFlowRate Qdot;
  Power Wdot_e, Wdot_s;

equation

  // Boundary equations
  V=1e-3 + 0.1e-3*(if time>0.5 then time else 0);
  mdot_in=0.1e-3;
  mdot_out=0.01e-3;
  h_in = 300190;
  h_out = h;

  // Conservation of mass
  der(m) = mdot_in-mdot_out;

  // Conservation of energy
  der(U) = h_in*mdot_in - h_out*mdot_out + Qdot - Wdot_e;
  Wdot_e - Wdot_s = P*der(V);

  // Specific internal energy (ideal gas)
  u = U/m;
  u = u_0+c_v*(T-T_0);

  // Specific enthalpy
  H = U+P*V;
  h = H/m;

  // Equation of state (ideal gas)
  P*V=m*R*T;
end BasicVolumeHeatTransfer;

connector FlowConnector
  import Modelica.SIunits.*;
  Pressure P;
  SpecificEnthalpy h;
  flow MassFlowRate mdot;
end FlowConnector;

model IdealGas "Property model for a simple ideal gas"
  import Modelica.SIunits.*;
  parameter SpecificInternalEnergy u_0=209058; // Air at T=293
  parameter SpecificHeatCapacity c_v = 717;
  parameter SpecificHeatCapacity R = 287;
  parameter Temperature T_0 = 293;
  Pressure P;
  SpecificVolume v;
  Temperature T;
  SpecificInternalEnergy u;
equation
  u=u_0+c_v*(T-T_0); // Linear model for inner energy
  P*v = R*T; // Equation of state: ideal gas
end IdealGas;

model BasicVolume "Control Volume with Mass Flow Interfaces"
  import Modelica.SIunits.*;
  extends IdealGas;
  FlowConnector inlet, outlet;
  parameter Mass m_0 = 0.00119;
  Mass m(start=m_0);
  Volume V;
  MassFlowRate mdot_in, mdot_out;
  SpecificEnthalpy h_in, h_out;
  SpecificEnthalpy h;
  Enthalpy H;
  InternalEnergy U(start=u_0*m_0);
  HeatFlowRate Qdot;
  Power Wdot_e, Wdot_s;
equation
  // Interface relations
  mdot_in = inlet.mdot;
  h_in = inlet.h;
  P = inlet.P;
  mdot_out = -outlet.mdot;
  h_out = outlet.h;
  P = outlet.P;
  h = h_out;
  // Conservation of mass
  der(m) = mdot_in-mdot_out;
  // Conservation of energy
  der(U) = h_in*mdot_in-h_out*mdot_out+Qdot-Wdot_e;
  /* Wdot_e=P*der(V)+Wdot_s; */ // Original formulation
  Wdot_e=P*(0.0001*((1/(1+exp((-200*(time-0.5)))))+(((-((time-0.5)/(1+exp((-200*(time-0.5))))))/(1+exp((-200*(time-0.5)))))*(exp((-200*(time-0.5)))*-200))))+Wdot_s; // Workaround
  H = U+P*V;
  u = U/m; // Specific internal energy (ideal gas)
  v = V/m; // Specific volume
  h = H/m; // Specific enthalpy
end BasicVolume;

model BasicVolumeTest "testing our general volume"
  extends BasicVolume;
equation
  // Boundary conditions: they are normally defined by the neighboring models
  inlet.h = 300190;
  // Specifying the remaining 3 degrees of freedom
  /*  V= 1e-3 + 0.1*(if time > 0.5 then time - 0.5 else 0);*/ // Original formulation
  V= 1e-3 + 0.1e-3*((time-0.5)/(1+exp(-2*100*(time-0.5)))); // Smooth approximation
  Qdot = 0.0;
  Wdot_s = 0.0;
end BasicVolumeTest;

model PressureEnthalpySource "Constant pressure and enthalpy source"
  import Modelica.SIunits.*;
  FlowConnector outlet;
  parameter Pressure P_0=1.0e5;
  parameter SpecificEnthalpy h_0 = 300190;
equation
  outlet.P=P_0;
  outlet.h=h_0;
end PressureEnthalpySource;

model PressureSink "Pressure sink, constant pressure"
  FlowConnector inlet;
  parameter Modelica.SIunits.Pressure P_0=1.0e5;
equation
  inlet.P=P_0;
end PressureSink;

model SimpleValveFlow
  import Modelica.SIunits;
  parameter SIunits.Area A=1e-4;
  parameter Real beta = 5.0e-5;
  SIunits.Pressure P_in, P_out;
  SIunits.MassFlowRate mdot;
equation // Boundary conditions
  P_in = 1.2e5;
  P_out = 1.0e5;
  // Constitutive relation
  beta*A^2*(P_in-P_out) = mdot^2;
end SimpleValveFlow;

model ValveFlow
  import Modelica.SIunits;
  FlowConnector inlet, outlet;
  parameter SIunits.Area A = 1e-4;
  parameter Real beta = 5.0e-5;
  SIunits.Pressure P_in, P_out;
  SIunits.MassFlowRate mdot;
equation // Interface relations
  P_in = inlet.P;
  P_out = outlet.P;
  mdot = inlet.mdot;
  0 = inlet.mdot + outlet.mdot; // Conservation of mass (no storage)
  // Conservation of energy (no storage)
  0 = inlet.mdot*inlet.h + outlet.mdot*outlet.h;
  beta*A^2*(P_in-P_out) = mdot^2; // Constitutive relation
  assert(P_in-P_out>0, "Error: model not valid for backflow");
end ValveFlow;


model ControlledValveFlow
  import Modelica.SIunits;
  import Modelica.Blocks.Interfaces;
  FlowConnector inlet, outlet;
  Interfaces.RealInput u;
  parameter SIunits.Area A = 1e-4;
  parameter Real beta = 5.0e-5;
  SIunits.Pressure P_in, P_out;
  SIunits.MassFlowRate mdot;
equation // Interface relations
  P_in = inlet.P;
  P_out = outlet.P;
  mdot = inlet.mdot;
  0 = inlet.mdot + outlet.mdot; // Conservation of mass (no storage)
  // Conservation of energy (no storage)
  0 = inlet.mdot*inlet.h + outlet.mdot*outlet.h;
  beta*(A*u)^2*(P_in-P_out) = mdot^2; // Constitutive relation
  assert(P_in-P_out>0, "Error: model not valid for backflow");
end ControlledValveFlow;

model SmoothStep "A smooth version of Modelica.Blocks.Sources.Step"
  parameter Real offset=0;
  parameter Real startTime=0;
  parameter Real height=1;
  parameter Real K = 100;
  Real y;
equation
  y = offset + height/(1+exp(-2*K*(time-startTime)));
end SmoothStep;

model CtrlFlowSystem
  import Modelica.Blocks.Sources;
  PressureEnthalpySource pressureEnthalpySource(P_0=1.2e5);
  SmoothStep step(offset=1,startTime=0.5,height=-0.5);
  ControlledValveFlow ctrlFlowModel; // negative step signal
  PressureSink pressureSink;
equation 
  connect(pressureEnthalpySource.outlet,ctrlFlowModel.inlet);
  connect(ctrlFlowModel.outlet,pressureSink.inlet);
  connect(step.y, ctrlFlowModel.u);
end CtrlFlowSystem;


  
  


