// Data structure in the target language holding data
#define GUESTOBJECT SwigRV
#define CASADI_NOT_NONE(p) !p || p->kind != SWIG_RV_NIL

// C# types of the casadi types without a proxy class of their own
%typemap(cstype) casadi::Matrix<casadi_int> "DM"
%typemap(cstype) casadi::SXElem "SX"
%typemap(cstype) casadi::GenericType "object"
// Host values (DM.full(), DM.sparse()) as their natural C# value
%typemap(cstype) SwigRV "object"
