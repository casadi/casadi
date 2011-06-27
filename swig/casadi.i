#ifdef SWIGOCTAVE
%module casadi_interface
#else
%module casadi
#endif

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
%include "std_pair.i"

%template(StringVector) std::vector<std::string>;

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

%template(Pair_Int_Int) std::pair<int,int>;
%template(VectorPair_Int_Int) std::vector< std::pair<int,int> >;
#endif // SWIG

// Lower value means wil be checked first
#define PRECEDENCE_DVector 98
#define PRECEDENCE_IVector 99

#define PRECEDENCE_DMatrix 100
#define PRECEDENCE_DMatrixVector 101
#define PRECEDENCE_SXMatrix 102
#define PRECEDENCE_SX 103
#define PRECEDENCE_SXMatrixVector 103
#define PRECEDENCE_MX 104
#define PRECEDENCE_MXVector 105
#define PRECEDENCE_PAIR_SLICE_SLICE 204
#define PRECEDENCE_SLICE 205
#define PRECEDENCE_IndexVector 210
#define PRECEDENCE_PAIR_IVector_IVector 206
#define PRECEDENCE_GENERICTYPE 300
#define PRECEDENCE_DICTIONARY 301


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
    snprintf(buffer, 32, "[%lu]", (unsigned long) $self->size());
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
    snprintf(buffer, 32, "[%lu]", (unsigned long) $self->size());
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

#ifdef SWIGOCTAVE
%inline %{
#include <octave/dim-vector.h>
%}
#endif

#ifdef SWIGPYTHON
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

#ifdef WITH_PYTHON_INTERRUPTS
#include <pythonrun.h>

%{
void SigIntHandler(int) {
  std::cerr << "Keyboard Interrupt" << std::endl;
	signal(SIGINT, SIG_DFL);
	kill(getpid(), SIGINT);
}
%}

%init %{
PyOS_setsig(SIGINT, SigIntHandler);
%}
#endif // WITH_PYTHON_INTERRUPTS
#endif // SWIGPYTHON

// typemaphelpers
%include "typemaphelpers.i"

// Auxilliary casadi functions
%include "casadi_aux.i"

// These methods must be added since the implicit type cast does not work
# define binopsFull(argtype,argCast,selfCast,returntype) \
returntype __pow__ (argtype) const{ return selfCast(*$self).__pow__(argCast(b));} \
returntype __rpow__(argtype) const{ return argCast(b).__pow__(selfCast(*$self));} \
returntype __add__ (argtype) const{ return selfCast(*$self) + argCast(b);} \
returntype __radd__(argtype) const{ return argCast(b) + selfCast(*$self);} \
returntype __sub__ (argtype) const{ return selfCast(*$self) - argCast(b);} \
returntype __rsub__(argtype) const{ return argCast(b) - selfCast(*$self);} \
returntype __mul__ (argtype) const{ return selfCast(*$self) * argCast(b);} \
returntype __rmul__(argtype) const{ return argCast(b) * selfCast(*$self);} \
returntype __div__ (argtype) const{ return selfCast(*$self) / argCast(b);} \
returntype __rdiv__(argtype) const{ return argCast(b) / selfCast(*$self);} \
returntype prod (argtype) const{ return prod(selfCast(*$self) , argCast(b));} \
returntype rprod (argtype) const{ return prod(argCast(b) , selfCast(*$self));} \
returntype __mrdivide__  (argtype) const { return selfCast(*$self).__mrdivide__(argCast(b));} \
returntype __rmrdivide__ (argtype) const { return argCast(b).__mrdivide__(selfCast(*$self));} \
returntype __ldivide__   (argtype) const { return selfCast(*$self).__mrdivide__(argCast(b));} \
returntype __rmldivide__ (argtype) const { return argCast(b).__mrdivide__(selfCast(*$self));} \
returntype __mpower__    (argtype) const{ return selfCast(*$self).__mpower__(argCast(b));} \
returntype __rmpower__   (argtype) const{ return argCast(b).__mpower__(selfCast(*$self));}

// Matrix typemaps class
%include "matrix.i"

// SX class
%include "sx.i"

// MX class
%include "mx.i"

// Matrix tools
%include "casadi/matrix/matrix_tools.hpp"

// Instansiate the functions
MATRIX_TOOLS_TEMPLATES(double)
MATRIX_TOOLS_TEMPLATES(CasADi::SX)

// FX
%include "fx.i"

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
#include "interfaces/csparse/csparse_tools.hpp"
%}
%include "interfaces/csparse/csparse.hpp"
%include "interfaces/csparse/csparse_tools.hpp"
#endif

#ifdef WITH_GSL
%{
#include "interfaces/gsl/gsl_integrator.hpp"
%}
%include "interfaces/gsl/gsl_integrator.hpp"
#endif
