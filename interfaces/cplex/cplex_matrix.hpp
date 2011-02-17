#ifndef CPLEX_MATRIX_HPP
#define CPLEX_MATRIX_HPP

#include "casadi/matrix/matrix.hpp"

namespace CasADi{

/// This class is used to convert csr matrix format used in casadi to the cplex sparse format (almost identical to csc)
class CplexMatrix{
  public:
    /// default contructor
    CplexMatrix();
    /// default destructor
    ~CplexMatrix();
    /// frees allocated memory
    void free();
    /// reads matrix in casadi format
    void read(const Matrix<double>& mat);
    /// returns non-zero values
    double* matval() const;
    /// returns indices of the beginning of columns
    int* matbeg() const;
    /// returns number of entries per column
    int* matcnt() const;
    /// returns row numbers
    int* matind() const;
  private:
    /// contains non-zero values column by column
    double* matval_; 
    /// matbeg[j] contains the index of the beginning of column j
    int* matbeg_;
    /// matcnt[j] contains the number of entries in column j
    int* matcnt_;
    /// matind[k] specifies the row number of the corresponding coefficient, matval[k]
    int* matind_;
};

} // namespace CasADi

#endif //CPLEX_MATRIX_HPP