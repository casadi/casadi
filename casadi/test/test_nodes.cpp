#include "../sx/sx_matrix.hpp"

#include <iostream>

#ifdef HAVE_UBLAS
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;
#endif

using namespace std;
using namespace CasADi;

class Base{
  public:
  Base(){   cout << "base def constr\n"; }
  Base(const Base& obj){   cout << "base copy constr\n"; }
  virtual ~Base(){   cout << "base destr\n"; }
};

class Derived : public Base{
  public:
  Derived(){ cout << "derived constr\n"; }
  Derived(const Derived& obj) : Base(obj) {   cout << "derived copy constr\n"; }
  Derived(int i){   cout << "int constr\n"; }
  virtual ~Derived(){   cout << "derived destr\n"; }
};


main(){

  
 // Derived der;
  Derived intder(5);
  Derived copy(intder);


  SX x("x");
  SX y("y");
  SX z("z");
  
  SX a = 2;
  SX b = 3.2;

  SX n1 = a + x;
  SX n2 = n1 * y;
  SX n3 = z - n1;
  SX n4 = n2 / n3 - b;

  cout << "n1 = " << n1 << endl;
  cout << "n2 = " << n2 << endl;
  cout << "n3 = " << n3 << endl;
  cout << "n4 = " << n4 << endl;

  #ifdef HAVE_UBLAS
  // Work with ublas
  ublas::matrix<SX> A(2,2);
  A(0,1) = x;
  A(1,1) = 6;
  A(1,0) = n3;
  cout << A << endl;

  ublas::matrix<SX> A2 = prod(A,A);
  cout << A2 << endl;
  
   Matrix A3 = A2;
   cout << A3 << endl;

  ublas::matrix<double> B(2,2);
  B(0,1) = 4;
  B(1,1) = 6;
  B(1,0) = 2;


  cout << B << endl;
  ublas::matrix<SX> Bexp = B;
  cout << Bexp << endl;

  Matrix Bexp2 = B;
  cout << Bexp2 << endl;
#endif // HAVE_UBLAS


  // Limits 
  SX two = 2;
  SX three = 3.0;

  SX ssnan = SX::nan;
  SX ssinf = SX::inf;
  SX ssminf = SX::minf;


  cout << two    << " is nan? " << two->isNan()    << " is inf? " << two->isInf()    << " is minusinf? " << two->isMinusInf()    << endl;
  cout << three  << " is nan? " << three->isNan()  << " is inf? " << three->isInf()  << " is minusinf? " << three->isMinusInf()  << endl;
  cout << ssnan  << " is nan? " << ssnan->isNan()  << " is inf? " << ssnan->isInf()  << " is minusinf? " << ssnan->isMinusInf()  << endl;
  cout << ssinf  << " is nan? " << ssinf->isNan()  << " is inf? " << ssinf->isInf()  << " is minusinf? " << ssinf->isMinusInf()  << endl;
  cout << ssminf << " is nan? " << ssminf->isNan() << " is inf? " << ssminf->isInf() << " is minusinf? " << ssminf->isMinusInf() << endl;

  // Comparisons
  cout << (x>2.0) << endl;
  cout << (x>2) << endl;


}
