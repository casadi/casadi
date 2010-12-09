optimization CSTR2_Opt(startTime=0,finalTime=10,objective=cost(finalTime))

  CSTRLib.Examples.CSTRs_Opt cstr(u1(min=0.1,max=3,start=1),u2(min=0.1,max=3,start=1),
                                     two_CSTRs_Series.V1(fixed=true),two_CSTRs_Series.T1(fixed=true),two_CSTRs_Series.CA1(fixed=true),
                                     two_CSTRs_Series.V2(fixed=true),two_CSTRs_Series.T2(fixed=true),two_CSTRs_Series.CA2(fixed=true));

  input Real der_u2;

  Real cost(start=0,fixed=true);

  parameter Real CA1_ref = 0.03;
  parameter Real CA2_ref = 0.001;
  parameter Real u1_ref = 1;
  parameter Real u2_ref = 1;
  parameter Real alpha = 1e5;
  parameter Real beta = 5e1;

equation
  // Fix u1, since this input has little effect on the outputs
  cstr.u1 = 1.1;
  // Derivative of input u2
  der(cstr.u2) = der_u2;
  der(cost) = alpha*(CA1_ref - cstr.two_CSTRs_Series.CA1)^2 + 
alpha*(CA2_ref - cstr.two_CSTRs_Series.CA2)^2 + beta*(u2_ref-cstr.u2)^2;

end CSTR2_Opt;