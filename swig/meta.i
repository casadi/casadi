%{
#include <sstream>
#include "casadi/stl_vector_tools.hpp"
#include "casadi/printable_object.hpp"
#include "casadi/shared_object.hpp"
#include "casadi/generic_type.hpp"
#include "casadi/options_functionality.hpp"
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/slice.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"
#include "casadi/fx/fx.hpp"
%}


%inline %{
template<> swig_type_info** meta< double >::name = &SWIGTYPE_p_double;
template<> swig_type_info** meta< int >::name = &SWIGTYPE_p_int;
template<> swig_type_info** meta< std::vector<double> >::name = &SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t;
template<> swig_type_info** meta< std::vector<int> >::name = &SWIGTYPE_p_std__vectorT_int_std__allocatorT_int_t_t;
template<> swig_type_info** meta< CasADi::GenericType >::name = &SWIGTYPE_p_CasADi__GenericType;
template<> swig_type_info** meta< CasADi::GenericType::Dictionary >::name = &SWIGTYPE_p_std__mapT_std__string_CasADi__GenericType_t;
template<> swig_type_info** meta< CasADi::SX >::name = &SWIGTYPE_p_CasADi__SX;
template<> swig_type_info** meta< CasADi::Matrix<CasADi::SX> >::name = &SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<CasADi::SX> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_CasADi__SX_t_std__allocatorT_CasADi__MatrixT_CasADi__SX_t_t_t;
template<> swig_type_info** meta< std::vector< CasADi::SX > >::name = &SWIGTYPE_p_std__vectorT_CasADi__SX_std__allocatorT_CasADi__SX_t_t;
template<> swig_type_info** meta< std::vector< std::vector< CasADi::SX > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__SX_std__allocatorT_CasADi__SX_t_t_std__allocatorT_std__vectorT_CasADi__SX_std__allocatorT_CasADi__SX_t_t_t_t;

template<> swig_type_info** meta< std::vector< CasADi::Matrix<double> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<int> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_int_t_std__allocatorT_CasADi__MatrixT_int_t_t_t;
template<> swig_type_info** meta< CasADi::CRSSparsity >::name = &SWIGTYPE_p_CasADi__CRSSparsity;
template<> swig_type_info** meta< CasADi::Matrix<double> >::name = &SWIGTYPE_p_CasADi__MatrixT_double_t;
template<> swig_type_info** meta< CasADi::Matrix<int> >::name = &SWIGTYPE_p_CasADi__MatrixT_int_t;
template<> swig_type_info** meta< CasADi::MX >::name = &SWIGTYPE_p_CasADi__MX;
template<> swig_type_info** meta< std::vector< CasADi::MX> >::name = &SWIGTYPE_p_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t;
	template<> swig_type_info** meta< std::vector< std::vector< CasADi::MX> > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__MX_p_std__allocatorT_CasADi__MX_p_t_t_std__allocatorT_std__vectorT_CasADi__MX_p_std__allocatorT_CasADi__MX_p_t_t_t_t;
template<> swig_type_info** meta< std::vector<std::string> >::name = &SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t;
%}



#ifdef SWIGOCTAVE
%include "meta_octave.i"
#endif

#ifdef SWIGPYTHON
%include "meta_python.i"
#endif
