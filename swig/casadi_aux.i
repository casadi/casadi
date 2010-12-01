%{
#include "casadi/options_functionality.hpp"
%}

// Exceptions handling
%include "exception.i"
%exception {
try {
$action
} catch (const std::exception& e) {
SWIG_exception(SWIG_RuntimeError, e.what());

} catch (const char* e) { // depreciated!!
SWIG_exception(SWIG_RuntimeError, e);
}
}

namespace CasADi {

class PrintableObject{
};

%extend PrintableObject{
std::string __str__()  { return $self->getDescription(); }
std::string __repr__()  { return $self->getRepresentation(); }
}

class SharedObject : public PrintableObject{
};

class OptionsFunctionality : public SharedObject{
	public:
		/** \brief  Print options to a stream */
		void printOptions(std::ostream &stream=std::cout) const;
};

%extend OptionsFunctionality {
void setOption(const std::string &name, const std::string& val){$self->setOption(name,val);} 
void setOption(const std::string &name, const std::vector<double>& val){$self->setOption(name,val);} 
void setOption(const std::string &name, double val){$self->setOption(name,val);} 

}

} // namespace CasADi
