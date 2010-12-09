within ;
package CSTRLib

 package Components
 model Two_CSTRs_Series
   // Import Modelica SI unit library
   import SI = Modelica.SIunits;

   // The model has to be rescaled in order to enable
   // use of SI units.

   // Model parameters
   parameter Real c1=10 "Constant 1";
   parameter Real c2=48.1909 "Constant 2";
   parameter Real qf=100 "Feed Volumetric Flow Rate (L/min)";
   parameter Real CAf=1.0 "Feed Concentreation of A (mol/L)";
   parameter SI.Temperature Tf=350 "Feed Temperature (K)";
   parameter Real ko=7.2e10
        "Pre-exponential Factor, Arrhenius Equation (min^-1)";
   parameter Real EoverR=1e4 "Activation Energy/Gas Constant (K)";
   parameter Real rho=1000 "Density of Fluid (g/L)";
   parameter Real Cp=0.239 "Heat Capacity of Fluid (J/(g*K))";
   parameter Real dH=4.78e4 "Heat of Reaction (J/mol)";

   // Algebraic variables
   Real q1(start=100,nominal=100) "Volumetric Flow Rate Leaving Reactor 1 (L/min)";
   Real q2(start=100,nominal=100) "Volumetric Flow Rate Leaving Reactor 2 (L/min)";
   Real Qc(start=-50,nominal=50) "Energy Flow from Reactor #1 - Cooling (J/min)";

  // Initial values for the states
   parameter Real V1_0=200;
   parameter Real CA1_0=0.03;
   parameter SI.Temperature T1_0=450;
   parameter Real V2_0=100;
   parameter Real CA2_0=0.002;
   parameter SI.Temperature T2_0=450;

   //Guess Values for steady-state solution to be
   //used as initial conditions
   Real V1(start=V1_0,nominal=V1_0) "Volume in Reactor 1 (L)";
   Real CA1(start=CA1_0,nominal=CA1_0) "Concentration of Reactor 1 (mol/L)";
   SI.Temperature T1(start=T1_0,nominal=T1_0) "Temperature in Reactor 1 (K)";
   Real V2(start=V2_0,nominal=V2_0) "Volume in Reactor 2 (L)";
   Real CA2(start=CA2_0,nominal=CA2_0) "Concentration of Reactor 2 (mol/L)";
   SI.Temperature T2(start=T2_0,nominal=T2_0) "Temperature in Reactor 2 (K)";

   // Model inputs
   Modelica.Blocks.Interfaces.RealInput u1(start=1) annotation (Placement(
            transformation(extent={{-100,20},{-60,60}}), iconTransformation(
              extent={{-100,20},{-60,60}})));
   Modelica.Blocks.Interfaces.RealInput u2(start=1)
        "Valve position at the outlet of the second reactor" annotation (
          Placement(transformation(extent={{-100,-60},{-60,-20}}),
            iconTransformation(extent={{-100,-60},{-60,-20}})));

 equation
   Qc = -c2*u2;
   q1 = c1*((V1 - V2)^.5);
   q2 = c1*u1*(V2^.5);

   der(V1) = qf - q1;
   der(CA1) = ((qf*CAf)/V1) - (ko*CA1*exp(-EoverR/T1)) - ((q1*CA1)/V1) - ((CA1*
     der(V1))/V1);
 //  der(T1) = ((qf*Tf)/V1) + ((dH*ko*CA1*exp(-EoverR/T1))/(rho*Cp)) - ((q1*T1)/V1)
 //     + Qc - ((T1*der(V1))/V1);
   der(T1) = ((qf*Tf)/V1) + ((dH*ko*CA1/exp(EoverR/400)*exp(-EoverR/T1+EoverR/400))/(rho*Cp)) - ((q1*T1)/V1)
      + Qc - ((T1*der(V1))/V1);

   der(V2) = q1 - q2;
   der(CA2) = ((q1*CA1)/V2) - (ko*CA2*exp(-EoverR/T2)) - ((q2*CA2)/V2) - ((CA2*
     der(V2))/V2);
   der(T2) = ((q1*T1)/V2) + ((dH*ko*CA2*exp(-EoverR/T2))/(rho*Cp)) - ((q2*T2)/V2)
      - ((T2*der(V2))/V2);

      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{
                -100,-100},{100,100}}), graphics), Icon(coordinateSystem(
              preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
            graphics={
            Rectangle(extent={{-40,28},{18,-22}}, lineColor={0,0,255}),
            Rectangle(extent={{10,-46},{72,-94}}, lineColor={0,0,255}),
            Text(
              extent={{-66,72},{104,40}},
              lineColor={0,0,255},
              textString="2 CSTRs")}));
 end Two_CSTRs_Series;

   model Two_CSTRs_stat_init
     extends Two_CSTRs_Series;
   initial equation
     der(V1) = 0;
     der(CA1) = 0;
     der(T1) = 0;
     der(V2) = 0;
     der(CA2) = 0;
     der(T2) = 0;

   end Two_CSTRs_stat_init;

 end Components;
  annotation (uses(Modelica(version="3.0.1")));

  package Examples
    model SimulationExperiment
      Components.Two_CSTRs_Series two_CSTRs_Series
        annotation (Placement(transformation(extent={{6,6},{26,26}})));
      Modelica.Blocks.Sources.Step step(
        startTime=1,
        height=1,
        offset=1)
        annotation (Placement(transformation(extent={{-60,28},{-40,48}})));
      Modelica.Blocks.Sources.Step step1(
        startTime=1,
        height=-0.5,
        offset=1)
        annotation (Placement(transformation(extent={{-50,-24},{-30,-4}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{
                -100,-100},{100,100}}), graphics));
      Components.Two_CSTRs_stat_init two_CSTRs_stat_init(
        CA2_0=1,
        T1_0(displayUnit="K"),
        T2_0(displayUnit="K"))
        annotation (Placement(transformation(extent={{8,-20},{28,0}})));
    equation
      connect(step.y, two_CSTRs_Series.u1) annotation (Line(
          points={{-39,38},{-16,38},{-16,20},{8,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(step1.y, two_CSTRs_Series.u2) annotation (Line(
          points={{-29,-14},{-12,-14},{-12,12},{8,12}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(two_CSTRs_stat_init.u1, step.y) annotation (Line(
          points={{10,-6},{-16,-6},{-16,38},{-39,38}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(two_CSTRs_stat_init.u2, step1.y) annotation (Line(
          points={{10,-14},{-29,-14}},
          color={0,0,127},
          smooth=Smooth.None));
    end SimulationExperiment;

    model CSTRs_Opt
      Components.Two_CSTRs_Series two_CSTRs_Series
        annotation (Placement(transformation(extent={{0,2},{20,22}})));
      Modelica.Blocks.Interfaces.RealInput u1
        annotation (Placement(transformation(extent={{-120,20},{-80,60}})));
      Modelica.Blocks.Interfaces.RealInput u2
        annotation (Placement(transformation(extent={{-120,-40},{-80,0}})));
      annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{
                -100,-100},{100,100}}), graphics));
    equation
      connect(u1, two_CSTRs_Series.u1) annotation (Line(
          points={{-100,40},{-50,40},{-50,16},{2,16}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(u2, two_CSTRs_Series.u2) annotation (Line(
          points={{-100,-20},{-50,-20},{-50,8},{2,8}},
          color={0,0,127},
          smooth=Smooth.None));
    end CSTRs_Opt;
  end Examples;
end CSTRLib;
