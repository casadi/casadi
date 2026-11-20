// Data structure in the target language holding data
#define GUESTOBJECT SWIG_WASMJS_Object

// Own std::string, bool, casadi_int and double typemaps: double passes through the wasm ABI
#define CASADI_TARGET_SCALAR_TYPEMAPS
