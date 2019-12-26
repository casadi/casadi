/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "function.hpp"
#include "../casadi_c.h"
#include "serializer.hpp"
#include <deque>

using namespace casadi;

static std::vector<Function> casadi_c_loaded_functions;
static std::deque<int> casadi_c_load_stack;
static int casadi_c_active = -1;

int casadi_c_int_width() {
  return sizeof(casadi_int);
}
int casadi_c_real_width() {
  return sizeof(double);
}

int casadi_c_id(const char* funname) {
  int ret = -1;
  std::string fname = funname;
  for (int i=0;i<casadi_c_loaded_functions.size();++i) {
    if (fname==casadi_c_loaded_functions.at(i).name()) {
      if (ret!=-1) {
        std::cerr << "Ambiguous function name '" << fname << "'" << std::endl;
        return -2;
      } else {
        ret = i;
      }
    }
  }
  if (ret==-1) {
    std::cerr << "Could not find function named '" << fname << "'." << std::endl;
    std::cerr << "Available functions: ";
    for (const auto& f : casadi_c_loaded_functions) {
      std::cerr << f.name() << " ";
    }
    std::cerr << std::endl;
    return -1;
  }
  return ret;
}

int casadi_c_n_loaded() { return casadi_c_loaded_functions.size(); }

inline int casadi_c_push_file_internal(const char *filename) {
  try {
    FileDeserializer fs(filename);
    auto type = fs.pop_type();
    if (type==SerializerBase::SerializationType::SERIALIZED_FUNCTION) {
      casadi_c_loaded_functions.push_back(fs.blind_unpack_function());
      return 0;
    } else if (type==SerializerBase::SerializationType::SERIALIZED_FUNCTION_VECTOR) {
      for (const Function& f : fs.blind_unpack_function_vector()) {
        casadi_c_loaded_functions.push_back(f);
      }
      return 0;
    } else {
      std::cerr << "Serializer file should contain a 'function' or 'function_vector'. "
                  "Got '" + SerializerBase::type_to_string(type) + "' instead." << std::endl;
      return -1;
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -2;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return -3;
  }
}

int casadi_c_push_file(const char *filename) {
  int before = casadi_c_loaded_functions.size();
  int ret = casadi_c_push_file_internal(filename);
  int after = casadi_c_loaded_functions.size();
  casadi_c_load_stack.push_back(after-before);
  return ret;
}

void casadi_c_clear(void) {
  casadi_c_load_stack.clear();
  casadi_c_loaded_functions.clear();
  casadi_c_active = -1;
}

void casadi_c_pop(void) {
  int count = casadi_c_load_stack.back();
  casadi_c_load_stack.pop_back();
  casadi_c_loaded_functions.erase(
    casadi_c_loaded_functions.begin()+(casadi_c_loaded_functions.size()-count),
    casadi_c_loaded_functions.end());
}

inline int sanitize_id(int id) {
  if (id<0 || id>=casadi_c_loaded_functions.size()) {
    std::cerr << "id " << id << " is out of range: must be in [0, ";
    std::cerr << casadi_c_loaded_functions.size() << "[" << std::endl;
    return 1;
  }
  return 0;
}

int casadi_c_activate(int id) {
  if (sanitize_id(id)) return -1;
  casadi_c_active = id;
  return 0;
}

void casadi_c_incref(void) {}
void casadi_c_decref(void) {}
void casadi_c_incref_id(int id) {}
void casadi_c_decref_id(int id) {}

int casadi_c_checkout(void) {
  return casadi_c_checkout_id(casadi_c_active);
}
int casadi_c_checkout_id(int id) {
  if (sanitize_id(id)) return -1;
  try {
    return casadi_c_loaded_functions.at(id).checkout();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -2;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return -3;
  }
}

void casadi_c_release(int mem) {
  casadi_c_release_id(casadi_c_active, mem);
}
void casadi_c_release_id(int id, int mem) {
  sanitize_id(id);
  try {
    casadi_c_loaded_functions.at(id).release(mem);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
  }
}

double casadi_c_default_in(casadi_int i) {
  return casadi_c_default_in_id(casadi_c_active, i);
}
double casadi_c_default_in_id(int id, casadi_int i) {
  if (sanitize_id(id)) return -1;
  try {
    return casadi_c_loaded_functions.at(id).default_in(i);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -2;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return -3;
  }
}

casadi_int casadi_c_n_in(void) {
  return casadi_c_n_in_id(casadi_c_active);
}
casadi_int casadi_c_n_in_id(int id) {
  if (sanitize_id(id)) return -1;
  try {
    return casadi_c_loaded_functions.at(id).n_in();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -2;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return -3;
  }
}

casadi_int casadi_c_n_out(void) {
  return casadi_c_n_out_id(casadi_c_active);
}
casadi_int casadi_c_n_out_id(int id) {
  if (sanitize_id(id)) return -1;
  try {
    return casadi_c_loaded_functions.at(id).n_out();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -2;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return -3;
  }
}

const char* casadi_c_name() {
  return casadi_c_name_id(casadi_c_active);
}
const char* casadi_c_name_id(int id) {
  if (sanitize_id(id)) return "";
  return casadi_c_loaded_functions.at(id).name().c_str();
}

const char* casadi_c_name_in(casadi_int i) {
  return casadi_c_name_in_id(casadi_c_active, i);
}
const char* casadi_c_name_in_id(int id, casadi_int i) {
  if (sanitize_id(id)) return "";
  try {
    return casadi_c_loaded_functions.at(id).name_in(i).c_str();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return "";
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return "";
  }
}

const char* casadi_c_name_out(casadi_int i) {
  return casadi_c_name_out_id(casadi_c_active, i);
}
const char* casadi_c_name_out_id(int id, casadi_int i) {
  if (sanitize_id(id)) return "";
  try {
    return casadi_c_loaded_functions.at(id).name_out(i).c_str();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return "";
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return "";
  }
}

const casadi_int* casadi_c_sparsity_in(casadi_int i) {
  return casadi_c_sparsity_in_id(casadi_c_active, i);
}
const casadi_int* casadi_c_sparsity_in_id(int id, casadi_int i) {
  if (sanitize_id(id)) return nullptr;
  try {
    return casadi_c_loaded_functions.at(id).sparsity_in(i);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return nullptr;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return nullptr;
  }
}

const casadi_int* casadi_c_sparsity_out(casadi_int i) {
  return casadi_c_sparsity_out_id(casadi_c_active, i);
}
const casadi_int* casadi_c_sparsity_out_id(int id, casadi_int i) {
  if (sanitize_id(id)) return nullptr;
  try {
    return casadi_c_loaded_functions.at(id).sparsity_out(i);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return nullptr;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return nullptr;
  }
}

int casadi_c_work(casadi_int *sz_arg, casadi_int* sz_res,
    casadi_int *sz_iw, casadi_int *sz_w) {
  return casadi_c_work_id(casadi_c_active, sz_arg, sz_res, sz_iw, sz_w);
}

int casadi_c_work_id(int id, casadi_int *sz_arg, casadi_int* sz_res,
    casadi_int *sz_iw, casadi_int *sz_w) {
  if (sanitize_id(id)) return -1;
  try {
    *sz_arg = casadi_c_loaded_functions.at(id).sz_arg();
    *sz_res = casadi_c_loaded_functions.at(id).sz_res();
    *sz_iw = casadi_c_loaded_functions.at(id).sz_iw();
    *sz_w = casadi_c_loaded_functions.at(id).sz_w();
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -2;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return -3;
  }
  return 0;
}

int casadi_c_eval(const double** arg, double** res, casadi_int* iw, double* w, int mem) {
  return casadi_c_eval_id(casadi_c_active, arg, res, iw, w, mem);
}

int casadi_c_eval_id(int id, const double** arg, double** res, casadi_int* iw, double* w, int mem) {
  if (sanitize_id(id)) return -1;
  try {
    return casadi_c_loaded_functions.at(id)(arg, res, iw, w, mem);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -2;
  } catch (...) {
    std::cerr << "Uncaught exception" << std::endl;
    return -3;
  }
  return 0;
}
