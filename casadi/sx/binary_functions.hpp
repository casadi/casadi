#ifndef BINARY_FUNCTIONS_HPP
#define BINARY_FUNCTIONS_HPP

#include <iostream>
#include <string>

#define RES(i,j) res[i*(i+1)/2+j]
#define F RES(0,0)
#define FX RES(1,0)
#define FY RES(1,1)
#define FXX RES(2,0)
#define FXY RES(2,1)
#define FYY RES(2,2)

namespace CasADi{
using namespace std;

//@{
/** \brief  Prototypes for a binary function */
typedef SX (*sfun_ptr)(const SX&, const SX&); // symbolic function
typedef SX (*sder_ptr)(const SX&, const SX&, const SX&); // symbolic partial derivative
typedef void (*print_ptr)(std::ostream &stream, const std::string& x, const std::string& y); // print as a c command
typedef void (*nfun_ptr)(double, double, double*);
//@}
	
/** \brief  Zero derivative */
static SX zero_sder(const SX &f, const SX& x, const SX& y){ return 0; }

//@{
/** \brief  Addition */
static SX add_sfcn(const SX& x, const SX& y){ return x+y; }
static SX add_sder1(const SX &f, const SX& x, const SX& y){ return 1; }
static SX add_sder2(const SX &f, const SX& x, const SX& y){ return 1; }
static void add_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "+" << y << ")";}
static void add_nfcn0(double x, double y, double *res){ F = x+y;}
static void add_nfcn1(double x, double y, double *res){ add_nfcn0(x,y,res); FX = 1;  FY = 1;}
static void add_nfcn2(double x, double y, double *res){ add_nfcn1(x,y,res); FXX = 0; FXY = 0; FYY = 0;}
//@}

//@{
/** \brief  Subtraction */
static SX sub_sfcn(const SX& x, const SX& y){ return x-y; }
static SX sub_sder1(const SX &f, const SX& x, const SX& y){ return 1; }
static SX sub_sder2(const SX &f, const SX& x, const SX& y){ return -1; }
static void sub_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "-" << y << ")";}
static void sub_nfcn0(double x, double y, double *res){ F = x-y;}
static void sub_nfcn1(double x, double y, double *res){ sub_nfcn0(x,y,res); FX = 1; FY = -1;}
static void sub_nfcn2(double x, double y, double *res){ sub_nfcn1(x,y,res); FXX = 0; FXY = 0; FYY = 0;}
//@}

//@{
/** \brief  Multiplication */
static SX mul_sfcn(const SX& x, const SX& y){ return x*y; }
static SX mul_sder1(const SX &f, const SX& x, const SX& y){ return y; }
static SX mul_sder2(const SX &f, const SX& x, const SX& y){ return x; }
static void mul_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "*" << y << ")";}
static void mul_nfcn0(double x, double y, double *res){ F = x*y;}
static void mul_nfcn1(double x, double y, double *res){ mul_nfcn0(x,y,res); FX = y; FY = x;}
static void mul_nfcn2(double x, double y, double *res){ mul_nfcn1(x,y,res); FXX = 0; FXY = 1; FYY = 0;}
//@}

//@{
/** \brief  Negation */
static SX neg_sfcn(const SX& x, const SX& y){ return -x; }
static SX neg_sder1(const SX &f, const SX& x, const SX& y){ return -1; }
static void neg_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "(-" << x << ")";}
static void neg_nfcn0(double x, double y, double *res){ F = -x;}
static void neg_nfcn1(double x, double y, double *res){ neg_nfcn0(x,y,res); FX = -1; FY = 0;}
static void neg_nfcn2(double x, double y, double *res){ neg_nfcn1(x,y,res); FXX = 0; FXY = 0; FYY = 0;}
//@}

//@{
/** \brief  Division */
static SX div_sfcn(const SX& x, const SX& y){ return x/y; }
static SX div_sder1(const SX &f, const SX& x, const SX& y){ return 1/y; }
static SX div_sder2(const SX &f, const SX& x, const SX& y){ return -f/y; }
static void div_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "/" << y << ")";}
static void div_nfcn0(double x, double y, double *res){ F = x/y;}
static void div_nfcn1(double x, double y, double *res){ div_nfcn0(x,y,res); FX = 1/y; FY = -F/y;}
static void div_nfcn2(double x, double y, double *res){ div_nfcn1(x,y,res); FXX = 0;  FXY = -FX*FX; FYY = -2*FY/y;}
//@}

//@{
/** \brief  Natural exponential */
static SX exp_sfcn(const SX& x, const SX& y){ return exp(x); }
static SX exp_sder1(const SX &f, const SX& x, const SX& y){ return f; }
static void exp_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "exp(" << x << ")";}
static void exp_nfcn0(double x, double y, double *res){ F = exp(x);}
static void exp_nfcn1(double x, double y, double *res){ exp_nfcn0(x,y,res); FX = F; FY = 0;}
static void exp_nfcn2(double x, double y, double *res){ exp_nfcn1(x,y,res); FXX = F; FXY = 0;  FYY = 0;}
//@}

//@{
/** \brief  Natural logarithm */
static SX log_sfcn(const SX& x, const SX& y){ return log(x); }
static SX log_sder1(const SX &f, const SX& x, const SX& y){ return 1/x; }
static void log_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "log(" << x << ")";}
static void log_nfcn0(double x, double y, double *res){ F = log(x);}
static void log_nfcn1(double x, double y, double *res){ log_nfcn0(x,y,res); FX = 1/x; FY = 0;}
static void log_nfcn2(double x, double y, double *res){ log_nfcn1(x,y,res); FXX = -1/(x*x);  FXY = 0;  FYY = 0;}
//@}

//@{
/** \brief  Power */
static SX pow_sfcn(const SX& x, const SX& y){ return pow(x,y); }
static SX pow_sder1(const SX &f, const SX& x, const SX& y){ return y*f/x; }
static SX pow_sder2(const SX &f, const SX& x, const SX& y){ return log(x)*f; }
static void pow_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "pow(" << x << "," << y << ")";}
static void pow_nfcn0(double x, double y, double *res){ F = pow(x,y);}
static void pow_nfcn1(double x, double y, double *res){ pow_nfcn0(x,y,res); FX = y*F/x; FY = log(x)*F;}
static void pow_nfcn2(double x, double y, double *res){ pow_nfcn1(x,y,res); FXX = FX/x*(y-1); FXY = F/x + FY/x*y;  FYY = log(x)*FY;}
//@}

//@{
/** \brief  Square root */
static SX sqrt_sfcn(const SX& x, const SX& y){ return sqrt(x); }
static SX sqrt_sder1(const SX &f, const SX& x, const SX& y){ return 1/(2*f);}
static void sqrt_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "sqrt(" << x << ")";}
static void sqrt_nfcn0(double x, double y, double *res){ F = sqrt(x);}
static void sqrt_nfcn1(double x, double y, double *res){ sqrt_nfcn0(x,y,res); FX = 0.5/F; FY = 0;}
static void sqrt_nfcn2(double x, double y, double *res){ sqrt_nfcn1(x,y,res); FXX = -0.5*FX/F;  FXY = 0;  FYY = 0;}
//@}

//@{
/** \brief  Sine */
static SX sin_sfcn(const SX& x, const SX& y){ return sin(x); }
static SX sin_sder1(const SX &f, const SX& x, const SX& y){ return cos(x); }
static void sin_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "sin(" << x << ")";}
static void sin_nfcn0(double x, double y, double *res){ F = sin(x);}
static void sin_nfcn1(double x, double y, double *res){ sin_nfcn0(x,y,res); FX = cos(x); FY = 0;}
static void sin_nfcn2(double x, double y, double *res){ sin_nfcn1(x,y,res); FXX = -F;  FXY = 0;  FYY = 0;}
//@}

//@{
/** \brief  Cosine */
static SX cos_sfcn(const SX& x, const SX& y){ return cos(x); }
static SX cos_sder1(const SX &f, const SX& x, const SX& y){ return -sin(x); }
static void cos_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "cos(" << x << ")";}
static void cos_nfcn0(double x, double y, double *res){ F = cos(x);}
static void cos_nfcn1(double x, double y, double *res){ cos_nfcn0(x,y,res); FX = -sin(x); FY = 0;}
static void cos_nfcn2(double x, double y, double *res){ cos_nfcn1(x,y,res); FXX = -F;  FXY = 0;  FYY = 0;}
//@}

//@{
/** \brief  Tangens */
static SX tan_sfcn(const SX& x, const SX& y){ return tan(x); }
static SX tan_sder1(const SX &f, const SX& x, const SX& y){ SX cosx = cos(x); return 1/(cosx*cosx);}
static void tan_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "tan(" << x << ")";}
static void tan_nfcn0(double x, double y, double *res){ F = tan(x);}
static void tan_nfcn1(double x, double y, double *res){ tan_nfcn0(x,y,res); double cosx = cos(x); FX = 1/(cosx*cosx); FY = 0;}
static void tan_nfcn2(double x, double y, double *res){ tan_nfcn1(x,y,res); FXX = 2*F*FX;  FXY = 0;  FYY = 0;}
//@}

//@{
/** \brief  arcsin */
static SX asin_sfcn(const SX& x, const SX& y){ return asin(x); }
static SX asin_sder1(const SX &f, const SX& x, const SX& y){ return 1/sqrt(1-x*x); }
static void asin_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "asin(" << x << ")";}
static void asin_nfcn0(double x, double y, double *res){ F = asin(x);}
static void asin_nfcn1(double x, double y, double *res){ asin_nfcn0(x,y,res); FX = 1/sqrt(1-x*x); FY = 0;}
static void asin_nfcn2(double x, double y, double *res){ asin_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

//@{
/** \brief  arccos */
static SX acos_sfcn(const SX& x, const SX& y){ return acos(x); }
static SX acos_sder1(const SX &f, const SX& x, const SX& y){ return -1/sqrt(1-x*x);}
static void acos_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "acos(" << x << ")";}
static void acos_nfcn0(double x, double y, double *res){ F = acos(x);}
static void acos_nfcn1(double x, double y, double *res){ acos_nfcn0(x,y,res); FX = -1/sqrt(1-x*x); FY = 0;}
static void acos_nfcn2(double x, double y, double *res){ acos_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

//@{
/** \brief  arctan */
static SX atan_sfcn(const SX& x, const SX& y){ return atan(x); }
static SX atan_sder1(const SX &f, const SX& x, const SX& y){ return 1/(1+x*x); }
static void atan_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "atan(" << x << ")";}
static void atan_nfcn0(double x, double y, double *res){ F = atan(x);}
static void atan_nfcn1(double x, double y, double *res){ atan_nfcn0(x,y,res); FX = 1/(1+x*x); FY = 0;}
static void atan_nfcn2(double x, double y, double *res){ atan_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

//@{
/** \brief  Step node */
static SX step_sfcn(const SX& x, const SX& y){ return x >= 0; }
static void step_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << ">=0)";}
static void step_nfcn0(double x, double y, double *res){ F = x >= 0;}
static void step_nfcn1(double x, double y, double *res){ step_nfcn0(x,y,res); FX = FY = 0;}
static void step_nfcn2(double x, double y, double *res){ step_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

//@{
/** \brief  Floor */
static SX floor_sfcn(const SX& x, const SX& y){ return floor(x); }
static void floor_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "floor(" << x << ")";}
static void floor_nfcn0(double x, double y, double *res){ F = floor(x);}
static void floor_nfcn1(double x, double y, double *res){ floor_nfcn0(x,y,res); FX = FY = 0;}
static void floor_nfcn2(double x, double y, double *res){ floor_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

//@{
/** \brief  Ceil */
static SX ceil_sfcn(const SX& x, const SX& y){ return ceil(x); }
static void ceil_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "ceil(" << x << ")";}
static void ceil_nfcn0(double x, double y, double *res){ F = ceil(x);}
static void ceil_nfcn1(double x, double y, double *res){ ceil_nfcn0(x,y,res); FX = FY = 0;}
static void ceil_nfcn2(double x, double y, double *res){ ceil_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

//@{
/** \brief  Min */
static SX fmin_sfcn(const SX& x, const SX& y){ return fmin(x,y); }
static SX fmin_sder1(const SX &f, const SX& x, const SX& y){ return x<=y; }
static SX fmin_sder2(const SX &f, const SX& x, const SX& y){ return y<=x; }
static void fmin_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "fmin(" << x << "," << y << ")";}
static void fmin_nfcn0(double x, double y, double *res){ F = fmin(x,y);}
static void fmin_nfcn1(double x, double y, double *res){ fmin_nfcn0(x,y,res); FX = x<=y; FY = !FX;}
static void fmin_nfcn2(double x, double y, double *res){ fmin_nfcn1(x,y,res); FXX = 0;   FXY = 0; FYY = 0;}
//@}

//@{
/** \brief  Max */
static SX fmax_sfcn(const SX& x, const SX& y){ return fmax(x,y); }
static SX fmax_sder1(const SX &f, const SX& x, const SX& y){ return x>=y; }
static SX fmax_sder2(const SX &f, const SX& x, const SX& y){ return y>=x; }
static void fmax_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "fmax(" << x << "," << y << ")";}
static void fmax_nfcn0(double x, double y, double *res){ F = fmax(x,y);}
static void fmax_nfcn1(double x, double y, double *res){ fmax_nfcn0(x,y,res); FX = x>=y; FY = !FX;}
static void fmax_nfcn2(double x, double y, double *res){ fmax_nfcn1(x,y,res); FXX = 0;   FXY = 0;  FYY = 0;}
//@}

//@{
/** \brief  Equality */
static SX equality_sfcn(const SX& x, const SX& y){ return x == y; }
static void equality_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "(" << x << "==" << y << ")";}
static void equality_nfcn0(double x, double y, double *res){ F = x==y;}
static void equality_nfcn1(double x, double y, double *res){ equality_nfcn0(x,y,res); FX = FY = 0;}
static void equality_nfcn2(double x, double y, double *res){ equality_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

//@{
/** \brief  Error function */
static SX erf_sfcn(const SX& x, const SX& y){ return erf(x); }
static SX erf_sder1(const SX &f, const SX& x, const SX& y){ return (2/sqrt(M_PI))*exp(-x*x);}
static void erf_print(std::ostream &stream, const std::string& x, const std::string& y){stream << "erf(" << x << ")";}
static void erf_nfcn0(double x, double y, double *res){ F = erf(x);}
static void erf_nfcn1(double x, double y, double *res){ erf_nfcn0(x,y,res); FX = (2/sqrt(M_PI))*exp(-x*x); FY = 0;}
static void erf_nfcn2(double x, double y, double *res){ erf_nfcn1(x,y,res); FXX = FXY = FYY = 0;}
//@}

/// Enum for quick access to any node
enum OPERATION{ // 
  ADD_NODE,  SUB_NODE,  MUL_NODE,  DIV_NODE,  NEG_NODE,  EXP_NODE,  LOG_NODE,  POW_NODE,  SQRT_NODE,  SIN_NODE,  COS_NODE,  TAN_NODE,  ASIN_NODE,  ACOS_NODE,  ATAN_NODE,  STEP_NODE,  FLOOR_NODE,  CEIL_NODE,  EQUALITY_NODE,  ERF_NODE,  FMIN_NODE,  FMAX_NODE
};
/// Vector of symbolic functions
static sfun_ptr sfcn[] = { // 
  add_sfcn,  sub_sfcn,  mul_sfcn,  div_sfcn,  neg_sfcn,  exp_sfcn,  log_sfcn,  pow_sfcn,  sqrt_sfcn,  sin_sfcn,  cos_sfcn,  tan_sfcn,  asin_sfcn,  acos_sfcn,  atan_sfcn,  step_sfcn,  floor_sfcn,  ceil_sfcn,  equality_sfcn,  erf_sfcn,  fmin_sfcn,  fmax_sfcn
};
///Array of symbolic partial derivatives
static sder_ptr sder1[] = { // 
  add_sder1, sub_sder1, mul_sder1, div_sder1, neg_sder1, exp_sder1, log_sder1, pow_sder1, sqrt_sder1, sin_sder1, cos_sder1, tan_sder1, asin_sder1, acos_sder1, atan_sder1, zero_sder,  zero_sder,   zero_sder,  zero_sder,      erf_sder1, fmin_sder1, fmax_sder1
};
/// Array of symbolic partial derivatives
static sder_ptr sder2[] = { 
  add_sder2, sub_sder2, mul_sder2, div_sder2, zero_sder, zero_sder, zero_sder, pow_sder2, zero_sder,  zero_sder, zero_sder, zero_sder, zero_sder,  zero_sder,  zero_sder,  zero_sder,  zero_sder,   zero_sder,  zero_sder,      zero_sder, fmin_sder2, fmax_sder2
};
/// Array of printing functions
static print_ptr print_c[] = { 
  add_print, sub_print, mul_print, div_print, neg_print, exp_print, log_print, pow_print, sqrt_print, sin_print, cos_print, tan_print, asin_print, acos_print, atan_print, step_print, floor_print, ceil_print, equality_print, erf_print, fmin_print, fmax_print
};
/// Partial derivatives up to order 0
static nfun_ptr nfun0[] = {   
  add_nfcn0, sub_nfcn0, mul_nfcn0, div_nfcn0, neg_nfcn0, exp_nfcn0, log_nfcn0, pow_nfcn0, sqrt_nfcn0, sin_nfcn0, cos_nfcn0, tan_nfcn0, asin_nfcn0, acos_nfcn0, atan_nfcn0, step_nfcn0, floor_nfcn0, ceil_nfcn0, equality_nfcn0, erf_nfcn0, fmin_nfcn0, fmax_nfcn0
};
/// Partial derivatives up to order 1
static nfun_ptr nfun1[] = {   
  add_nfcn1, sub_nfcn1, mul_nfcn1, div_nfcn1, neg_nfcn1, exp_nfcn1, log_nfcn1, pow_nfcn1, sqrt_nfcn1, sin_nfcn1, cos_nfcn1, tan_nfcn1, asin_nfcn1, acos_nfcn1, atan_nfcn1, step_nfcn1, floor_nfcn1, ceil_nfcn1, equality_nfcn1, erf_nfcn1, fmin_nfcn1, fmax_nfcn1
};
/// Partial derivatives up to order 2
static nfun_ptr nfun2[] = {   
  add_nfcn2, sub_nfcn2, mul_nfcn2, div_nfcn2, neg_nfcn2, exp_nfcn2, log_nfcn2, pow_nfcn2, sqrt_nfcn2, sin_nfcn2, cos_nfcn2, tan_nfcn2, asin_nfcn2, acos_nfcn2, atan_nfcn2, step_nfcn2, floor_nfcn2, ceil_nfcn2, equality_nfcn2, erf_nfcn2, fmin_nfcn2, fmax_nfcn2
};


} // namespace CasADi

#undef RES
#undef F
#undef FX
#undef FY
#undef FXX
#undef FXY
#undef FYY

#endif //BINARY_FUNCTIONS_HPP
