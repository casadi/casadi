%module casadi

// Turn off the warnings that certain methods are effectively ignored, this seams to be a false warning, 
// for example vertcat(SXMatrixVector), vertcat(DMatrixVector) and vertcat(MXVector) appears to work fine
#pragma SWIG nowarn=509

#ifdef WITH_DOXDOC
%include "../doc/doc.i"
#endif 

#ifndef WITH_DOXDOC
//#warning "Not using doxygen. Run /doc/doc2swig.sh to include it. A make rebuild-cache may be necessary."
#endif

%feature("autodoc", "1");

// STL
%include "std_string.i"
%include "std_vector.i"

// Move to C++
#ifndef SWIG
namespace CasADi{
  typedef std::vector<bool> BVector;
  typedef std::vector<std::vector<bool> > BVectorVector;
  typedef std::vector< std::vector<std::vector<bool> > > BVectorVectorVector;
  
  typedef std::vector<int> IVector;
  typedef std::vector<std::vector<int> > IVectorVector;
  typedef std::vector< std::vector<std::vector<int> > > IVectorVectorVector;
  
  typedef std::vector<double> DVector;
  typedef std::vector<std::vector<double> > DVectorVector;
  typedef std::vector< std::vector<std::vector<double> > > DVectorVectorVector;
} // namespace CasADi
#else // SWIG
%template(BVector)             std::vector<bool>;
%template(BVectorVector)       std::vector<std::vector<bool> > ;
%template(BVectorVectorVector) std::vector< std::vector<std::vector<bool> > > ;

%template(IVector)             std::vector<int>;
%template(IVectorVector)       std::vector<std::vector<int> > ;
%template(IVectorVectorVector) std::vector< std::vector<std::vector<int> > > ;

%template(DVector)             std::vector<double>;
%template(DVectorVector)       std::vector<std::vector<double> > ;
%template(DVectorVectorVector) std::vector< std::vector<std::vector<double> > > ;

#endif // SWIG

// The following is a work-around since it appears not possible to use the standard print functions from stl_vector tools,
// nor the std::stringstream class, since these are included _after_ std::vector in the C++ generated wrapper code
%extend std::vector<double>{  
  std::string __repr__(){
    char buffer[32];
    std::string ret;
    ret += "[";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%g", $self->at(i));
      ret += buffer;
    }
    ret += "]";
    return ret;
  }
  std::string __str__(){
    char buffer[32];
    std::string ret;
    
    // Add dimension
    snprintf(buffer, 32, "[%d]", $self->size());
    ret += buffer; 
    ret += "(";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%g", $self->at(i));
      ret += buffer;
    }
    ret += ")";
    return ret;
  }
};

// The same workaround for vector<double>
%extend std::vector<int>{  
  std::string __repr__(){
    char buffer[32];
    std::string ret;
    ret += "[";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%d", $self->at(i));
      ret += buffer;
    }
    ret += "]";
    return ret;
  }
  std::string __str__(){
    char buffer[32];
    std::string ret;
    
    // Add dimension
    snprintf(buffer, 32, "[%d]", $self->size());
    ret += buffer; 
    ret += "(";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%d", $self->at(i));
      ret += buffer;
    }
    ret += ")";
    return ret;
  }
};

#ifndef WITH_NUMPY
//#warning "Not using numpy. option(WITH_NUMPY = OFF)"
#endif

#ifdef WITH_NUMPY
//#warning "Using numpy. option(WITH_NUMPY = ON)"
%{
#define SWIG_FILE_WITH_INIT
%}
// Get the NumPy typemaps
%include "numpy.i"
%init %{
  import_array();
%}
#endif

// Auxilliary casadi functions
%include "casadi_aux.i"

// Matrix typemaps class
%include "matrixtypemaps.i"

// SX class
%include "sx.i"

// MX class
%include "mx.i"

// FX
%include "fx.i"

// Matrix tools
%include "casadi/matrix/matrix_tools.hpp"

// Instansiate the functions
MATRIX_TOOLS_TEMPLATES(double)
MATRIX_TOOLS_TEMPLATES(CasADi::SX)

// Sparsity tools
%{
#include "casadi/matrix/sparsity_tools.hpp"
%}
%include "casadi/matrix/sparsity_tools.hpp"

// SX tools
%include "casadi/sx/sx_tools.hpp"

// Optimal control
%include "optimal_control.i"

// IPOPT
#ifdef WITH_IPOPT
%include "ipopt_interface.i"
#endif

// ACADO
#ifdef WITH_ACADO
%include "acado_interface.i"
#endif

// Sundials
#ifdef WITH_SUNDIALS
%include "sundials_interface.i"
#endif

// SuperLU
#ifdef WITH_SUPERLU
%include "superlu.i"
#endif

#ifdef WITH_LAPACK
%include "lapack_interface.i"
#endif

#ifdef WITH_KNITRO
%{ 
  #include "interfaces/knitro/knitro_solver.hpp"
%}
%include "interfaces/knitro/knitro_solver.hpp"
#endif

#ifdef WITH_CPLEX
%{ 
  #include "interfaces/cplex/cplex_solver.hpp"
%}
#include "interfaces/cplex/cplex_solver.hpp"
#endif

#ifdef WITH_LIFTOPT
%{
  #include "interfaces/liftopt/liftopt_solver.hpp"
%}
%include "interfaces/liftopt/liftopt_solver.hpp"
#endif

#ifdef WITH_CSPARSE
%{
#include "interfaces/csparse/csparse.hpp"
%}
%include "interfaces/csparse/csparse.hpp"
#endif

