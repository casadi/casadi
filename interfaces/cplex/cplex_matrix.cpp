#include "cplex_matrix.hpp"

using namespace CasADi;

CplexMatrix::CplexMatrix(){
  matval_ = NULL;
  matbeg_ = NULL;
  matcnt_ = NULL;
  matind_ = NULL;
}

CplexMatrix::~CplexMatrix(){
  this->free();
}

void CplexMatrix::free(){
  if (matval_){
    delete[] matval_;
    matval_ = NULL;
  }
  if (matbeg_){
    delete[] matbeg_;
    matbeg_ = NULL;
  }
  if (matcnt_){
    delete[] matcnt_;
    matcnt_ = NULL;
  }
  if (matind_){
    delete[] matind_;
    matind_ = NULL;
  }
}

void CplexMatrix::read(const Matrix<double>& mat){
  int nnz = mat.size();
  int nc = mat.size2();
  int nr = mat.size1();
  
  double* val; // non zero elements in csr format
  int* col_ind; // column indices in csr format
  int* row_ptr; // row pointer in csr format
  
  val = &mat[0];
  col_ind = &mat.col()[0];
  row_ptr = &mat.rowind()[0];
  
  matval_ = new double[nnz];
  matbeg_ = new int[nc+1];
  matcnt_ = new int[nc];
  matind_ = new int[nnz];
  
  e_cnt = 0; // element counter
  for(int c_ind=0 ; c_ind<nc; ++c_ind){
    int c_cnt = 0; // column element counter
    for(int r_ind=0; r_ind<nr; ++r_ind){
      int r_start = row_ptr[r_ind];
      int r_end = row_ptr[r_ind+1];
      for(int el=r_start; el<r_end; ++el){
        if(col==c_ind){
          matval_[e_cnt] = val[el];
          matind_[e_cnt] = r_ind;
          ++c_cnt;
          ++e_cnt;
          break;
        } else if (col>c_ind) {
          break;
        }
      }
    }
    matcnt_[c_ind]=c_cnt;
  }
  matbeg_[nc] = nnz+1;
}

double* CplexMatrix::matval() const{
  return matval_;
}

int* CplexMatrix::matbeg() const{
  return matebeg_;
}

int* CplexMatrix::matcnt() const{
  return matcnt_;
}

int* CplexMatrix::matind() const{
  return matind_;
}

