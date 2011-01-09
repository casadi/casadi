%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_matrix.hpp"
#include "casadi/sx/sx_tools.hpp"
%}

%include "typemaps.i"
%include "casadi/matrix/crs_sparsity.hpp"
%include "casadi/matrix/matrix.hpp"
// %include "casadi/sx/sx.hpp"

#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

%include "casadi/sx/sx.hpp"


// Template instantiations
%template(vector_PyObject) std::vector<PyObject*>;
%template(vector_sx) std::vector<CasADi::SX>;
%template(vector_vector_sx) std::vector< std::vector<CasADi::SX> >;
%template(vector_vector_vector_sx) std::vector< std::vector< std::vector<CasADi::SX> > >;

%template(matrix_sx) CasADi::Matrix<CasADi::SX>;
%template(matrix_double) CasADi::Matrix<double>;
%template(vector_matrix_sx) std::vector<CasADi::Matrix<CasADi::SX> >;
%template(vector_matrix_double) std::vector<CasADi::Matrix<double> >;

namespace CasADi{

#ifdef WITH_IMPLICITCONV
%implicitconv SXMatrix;
#endif WITH_IMPLICITCONV

// NOTE NOTE NOTE: SXMatrix should be deleted!!!
/*class SXMatrix : public std::vector<SX>{*/ // PROBLEM WITH INHERITANCE!
class SXMatrix{
public:

/** \brief  constructors */
SXMatrix();                               // empty 0-by-0 matrix
SXMatrix(int n, int m);                   // empty n-by-m matrix
SXMatrix(int n, int m, const SX& val);    // dense n-by-m matrix filled with val

/** \brief  These constructors enable implicit type conversion */
SXMatrix(const SX &scalar);      // create a variable from a scalar
SXMatrix(double val);            // constant
SXMatrix(const std::vector<SX>& x);
SXMatrix(const std::vector<SX>& x,  int n, int m);

/** \brief  Create a matrix of symbolic variables  */
explicit SXMatrix(const std::string& name, int n=1, int m=1);   // variable

/** \brief  destructor */
~SXMatrix();
  
void clear();
void resize(int n, int m);
void reserve(int nnz);

int numel() const;       // get the number of elements
int size1() const;       // get the first dimension
int size2() const;       // get the second dimension
int size() const;        // number of non-zero elements
bool empty() const; // is the matrix empty
bool scalar() const; // is the matrix scalar
bool vector() const; // is the matrix a vector

// Get sparsity in compressed row storage (CRS) format
const std::vector<int>& col() const; // vector of length nnz containing the columns for all the indices of the non-zero elements
const std::vector<int>& rowind() const; // vector of length n+1 containing the index of the last non-zero element up till each row 
};

%extend SXMatrix {
std::string __str__()  { return $self->getRepresentation(); }
std::string __repr__() { return $self->getRepresentation(); }

// Matrix product (quick fix)
SXMatrix dot(const SXMatrix &y){ return prod(*$self,y); }

// Get and set elements
SX __getitem__(int i){ return $self->getElement(i);}
SX __getitem__(const std::vector<int> &I ){ if(I.size()!=2) throw CasADi::CasadiException("__getitem__: not 2D"); return $self->getElement(I[0],I[1]);}



SXMatrix __getitem__(const std::vector<PyObject*> &I ) {
      if(I.size()!=2) throw CasADi::CasadiException("__getitem__(vector<PyObject*>): not 2D");
      int  	i[2]; // numpy - style:  i:n:k
		  int  	n[2];
		  int  	k[2];
		  
		  int s[2]={(*$self).size1(),(*$self).size2()};
		  
		  for (int j=0;j<2;j++) {
		    PyObject *q = I[j];
        if( PySlice_Check(q) ) {
          PySliceObject* slice=(PySliceObject*) q;
          if (slice->start==Py_None) {i[j]=0; } else {i[j]=PyInt_AsLong(slice->start);}
          if (slice->stop==Py_None)  {n[j]=s[j];} else {n[j]=PyInt_AsLong(slice->stop);}
          if (slice->step==Py_None)  {k[j]=1; } else {k[j]=PyInt_AsLong(slice->step);}
          if (i[j]<0)
            i[j]+=s[j];
          if (n[j]<0)
            n[j]+=s[j];
        } else if (PyInt_Check(q)) {
          i[j]=PyInt_AsLong(q);
          if (i[j]<0)
            i[j]+=s[j];
          n[j]=i[j]+1;
          k[j]=1;
        } else {
          SWIG_Error(SWIG_TypeError, "__getitem__: Slice or int expected");
        }

		  }
		  
		  CasADi::SXMatrix res;
      getSub(res,*$self,i[0],i[1],n[0]-i[0],n[1]-i[1],k[0],k[1]);
      return res;
}

#define SETTERS(T) \
void __setitem__(int i, T el){ $self->getElementRef(i) = el;}\
void __setitem__(const std::vector<int> &I, T el){ if(I.size()!=2) throw CasADi::CasadiException("__setitem__: not 2D"); $self->getElementRef(I[0],I[1]) = el;}
SETTERS(const SX&)
SETTERS(double)
#undef SETTERS

binops(SXMatrix, const SX&)
binops(SXMatrix, const SXMatrix&)
binops(SXMatrix, double)
binops(SXMatrix, const std::vector<SX>&)
unops(SXMatrix)

}

} // namespace CasADi

// Template instantiations
%template(vector_sx_matrix) std::vector<CasADi::SXMatrix>;

#undef unops
#undef binops




//namespace CasADi{

// This simple typemap will convert a returned SXMAtrix to a numpy array of SX's

//%typemap(out) SXMatrix
//{
 //    #PyObject * 
 //    PyErr_SetString(PyExc_TypeError,"Hoorah");
 //    #$result = PySXMatrix2Numpy($1);
//}
//}



