package QuadTank_pack
	        
  model QuadTank
    // Process parameters
	parameter Modelica.SIunits.Area A1=2.8e-3, A2=3.2e-3, A3=2.8e-3, A4=3.2e-3;
	parameter Modelica.SIunits.Area a1=7.1e-6, a2=5.7e-6, a3=7.1e-6, a4=5.7e-6;
	parameter Modelica.SIunits.Acceleration g=9.81;
	parameter Real k1_nmp(unit="m/s/V") = 3.14e-6, k2_nmp(unit="m/s/V") = 3.29e-6;
	parameter Real g1_nmp=0.70, g2_nmp=0.70;

    // Initial tank levels
	parameter Modelica.SIunits.Length x1_0 = 0.04102638;
	parameter Modelica.SIunits.Length x2_0 = 0.06607553;
	parameter Modelica.SIunits.Length x3_0 = 0.00393984;
	parameter Modelica.SIunits.Length x4_0 = 0.00556818;
	
    // Tank levels
	Modelica.SIunits.Length x1(start=x1_0,min=0.0001/*,max=0.20*/);
	Modelica.SIunits.Length x2(start=x2_0,min=0.0001/*,max=0.20*/);
	Modelica.SIunits.Length x3(start=x3_0,min=0.0001/*,max=0.20*/);
	Modelica.SIunits.Length x4(start=x4_0,min=0.0001/*,max=0.20*/);

	// Inputs
	input Modelica.SIunits.Voltage u1;
	input Modelica.SIunits.Voltage u2;

  equation
    	der(x1) = -a1/A1*sqrt(2*g*x1) + a3/A1*sqrt(2*g*x3) +
					g1_nmp*k1_nmp/A1*u1;
	der(x2) = -a2/A2*sqrt(2*g*x2) + a4/A2*sqrt(2*g*x4) +
					g2_nmp*k2_nmp/A2*u2;
	der(x3) = -a3/A3*sqrt(2*g*x3) + (1-g2_nmp)*k2_nmp/A3*u2;
	der(x4) = -a4/A4*sqrt(2*g*x4) + (1-g1_nmp)*k1_nmp/A4*u1;

end QuadTank;

optimization QuadTank_Opt (objective = cost(finalTime),
                     startTime = 0,
                     finalTime = 50)

    extends QuadTank(u1(initialGuess=u1_r),u2(initialGuess=u2_r),
                     x1(initialGuess=x1_0,fixed=true),
    	         x2(initialGuess=x2_0,fixed=true),
                     x3(initialGuess=x3_0,fixed=true),
                     x4(initialGuess=x4_0,fixed=true));

// Reference values
    parameter Modelica.SIunits.Length x1_r = 0.06410371;
    parameter Modelica.SIunits.Length x2_r = 0.10324302;
    parameter Modelica.SIunits.Length x3_r = 0.006156;
    parameter Modelica.SIunits.Length x4_r = 0.00870028;
    parameter Modelica.SIunits.Voltage u1_r = 2.5;
    parameter Modelica.SIunits.Voltage u2_r = 2.5;

    Real cost(start=0,fixed=true);

    equation
      der(cost) = 40000*((x1_r - x1))^2 + 
                        40000*((x2_r - x2))^2 + 
                        40000*((x3_r - x3))^2 + 
                        40000*((x4_r - x4))^2 + 
                        ((u1_r - u1))^2 + 
                        ((u2_r - u2))^2;

end QuadTank_Opt;

optimization QuadTank_Static(objective=(x1_meas-x1)^2 + (x2_meas-x2)^2 + 
                                       (x3_meas-x3)^2 + (x4_meas-x4)^2, static=true)
  
  extends QuadTank(a1(free=true),a2(free=true));

 parameter Real x1_meas = x1_0+0.01;
 parameter Real x2_meas = x2_0-0.01;
 parameter Real x3_meas = x3_0-0.01;
 parameter Real x4_meas = x4_0+0.01;

initial equation
 der(x1) = 0;
 der(x2) = 0;
 der(x3) = 0;
 der(x4) = 0;

end QuadTank_Static;

end QuadTank_pack;
