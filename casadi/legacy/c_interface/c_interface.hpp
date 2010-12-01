#ifndef C_INTERFACE_HPP
#define C_INTERFACE_HPP

#include <vector>
#include <iostream>
#include <sstream>

// C++ to C error handling
#define c_fcn_wrapper(name,x) \
  try{ \
    x \
    return 0; \
  } catch(exception &e){ \
    cerr << "Exception thrown in module \"" << name << "\": " << e.what() << endl; \
  } catch(const char* str){ \
    cerr << "Error (old style!) in module \"" << name << "\": " << str << endl; \
  } catch(...){ \
    cerr << "Unknown error in module \"" << name << "\"" << endl; \
  } \
  return 1;

// For constructors
#define c_constructor_wrapper(name,x) \
try{ \
  x \
  } catch(exception &e){ \
    cerr << "Exception thrown in module \"" << name << "\": " << e.what() << endl; \
  } catch(const char* str){ \
    cerr << "Error (old style!) in module \"" << name << "\": " << str << endl; \
  } catch(...){ \
    cerr << "Unknown error in module \"" << name << "\"" << endl; \
  } \
  return 0; \

  
#endif // C_INTERFACE_HPP

