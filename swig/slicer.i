%{
#include "casadi/mx/slicer.hpp"
%}

namespace CasADi{
#ifdef WITH_IMPLICITCONV
%implicitconv Slicer;
#endif WITH_IMPLICITCONV

class SlicerPrimitive {};

class SlicerPrimitiveFromTo : public SlicerPrimitive {
public:
	SlicerPrimitiveFromTo(int from, int to, int step=1);
};

class SlicerPrimitiveAll : public SlicerPrimitiveFromTo {
public:
	SlicerPrimitiveAll();
};

class Slicer {
public:

	Slicer(const char *s_);
	Slicer(const std::string &s);
	Slicer(std::vector<int> &i);
	Slicer(int i);
	Slicer(const SlicerPrimitiveAll &p_);
	Slicer(const Slicer &s);
	Slicer(const SlicerPrimitiveFromTo &p_);
};

//%typemap(in) SlicerPrimitiveFromTo {
//        $1 = PySliceObjectToSlicerPrimitiveFromTo($input);
//}


} // namespace CasADi

%inline %{
CasADi::SlicerPrimitiveFromTo PySliceObjectToSlicerPrimitiveFromTo(PySliceObject* slice) {
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
	return CasADi::SlicerPrimitiveFromTo(0,-1);
      }
      PyObject *startptr = slice->start;
      PyObject *stopptr = slice->stop;
      PyObject *stepptr = slice->step;
      int start=0;
      int stop=-1;
      int step=1;
      if (PyInt_Check(startptr)) {start=PyInt_AsLong(startptr);}
      if (PyInt_Check(stopptr))  {stop=PyInt_AsLong(stopptr);}
      if (PyInt_Check(stepptr))  {step=PyInt_AsLong(stepptr);}

      return CasADi::SlicerPrimitiveFromTo(start,stop, step);
}
%}
