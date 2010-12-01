#ifndef MATRIX_SIZE_HPP
#define MATRIX_SIZE_HPP

#include <iostream>

namespace CasADi{

/** \brief  Container class for size of a matrix
  \author Joel Andersson 
  \date 2010	
*/
class MatrixSize{
  public:
  MatrixSize() : nrow(0), ncol(0){ }  // default constructor 0-by-0 matrix
  MatrixSize(int nrow_, int ncol_) : nrow(nrow_), ncol(ncol_){ }  // size n-by-m

/*  bool operator==(const MatrixSize &x) const{
     return nrow == x.nrow && ncol == x.ncol;
  }

  bool operator!=(const MatrixSize &x) const{
      return nrow != x.nrow || ncol != x.ncol;
  }*/
 
 friend std::ostream& operator<<(std::ostream &stream, const MatrixSize &x) {
      return stream << "[" << x.nrow << "," << x.ncol << "]"; 
  }

  int nrow, ncol;
};

} // namespace CasADi

#endif // MATRIX_SIZE_HPP
