#include <casadi/casadi.hpp>

// Templated blackbox function
template<typename T>
T testfun(T x, T y) {
  // Conditionals
  T z = x > y ? x + sin(y) : y + sin(x);
  // While loop
  while (z < 20) {
    z = z + 4;
  }
  // Done
  return z;
}

// Class with forward derivative calculation
class FwdAD {
 public:
  // Nondifferentiated value
  double v;
  // Sensitivities
  double s;
  // Constructor with implicit type conversion
  FwdAD(double v, double s = 0) : v(v), s(s) {}
  // Addition
  inline friend FwdAD operator+(FwdAD x, FwdAD y) {
    return FwdAD(x.v + y.v, x.s + y.s);
  }
  // Sine
  inline friend FwdAD sin(FwdAD x) {
    return FwdAD(std::sin(x.v), std::cos(x.v) * x.s);
  }
  // Comparisons return booleans
  inline friend bool operator==(FwdAD x, FwdAD y) {
    return x.v == y.v;
  }
  inline friend bool operator!=(FwdAD x, FwdAD y) {
    return x.v != y.v;
  }
  inline friend bool operator<(FwdAD x, FwdAD y) {
    return x.v < y.v;
  }
  inline friend bool operator>(FwdAD x, FwdAD y) {
    return x.v > y.v;
  }
  inline friend bool operator<=(FwdAD x, FwdAD y) {
    return x.v <= y.v;
  }
  inline friend bool operator>=(FwdAD x, FwdAD y) {
    return x.v >= y.v;
  }
};


int main(){
  FwdAD x = 2;
  x.s = 1;
  FwdAD z = testfun(x, x);



  std::cout << "here: z.v = " << z.v << ", z.s = " << z.s << std::endl;

  return 0;
}
