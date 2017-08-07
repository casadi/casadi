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



#include "code_generator.hpp"
#include "function_internal.hpp"
#include <iomanip>
#include <casadi_runtime_str.h>

using namespace std;
namespace casadi {

  CodeGenerator::CodeGenerator(const std::string& name, const Dict& opts) {
    // Default options
    this->verbose = true;
    this->mex = false;
    this->cpp = false;
    this->main = false;
    this->real_t = "double";
    this->codegen_scalars = false;
    this->with_header = false;
    this->with_mem = false;
    this->with_export = true;
    indent_ = 2;

    // Read options
    for (auto&& e : opts) {
      if (e.first=="verbose") {
        this->verbose = e.second;
      } else if (e.first=="mex") {
        this->mex = e.second;
      } else if (e.first=="cpp") {
        this->cpp = e.second;
      } else if (e.first=="main") {
        this->main = e.second;
      } else if (e.first=="real_t") {
        this->real_t = e.second.to_string();
      } else if (e.first=="codegen_scalars") {
        this->codegen_scalars = e.second;
      } else if (e.first=="with_header") {
        this->with_header = e.second;
      } else if (e.first=="with_mem") {
        this->with_mem = e.second;
      } else if (e.first=="with_export") {
        this->with_export = e.second;
      } else if (e.first=="indent") {
        indent_ = e.second;
        casadi_assert(indent_>=0);
      } else {
        casadi_error("Unrecongnized option: " << e.first);
      }
    }

    // Start at new line with no indentation
    newline_ = true;
    current_indent_ = 0;

    // Divide name into base and suffix (if any)
    string::size_type dotpos = name.rfind('.');
    if (dotpos==string::npos) {
      this->name = name;
      this->suffix = this->cpp ? ".cpp" : ".c";
    } else {
      this->name = name.substr(0, dotpos);
      this->suffix = name.substr(dotpos);
    }

    // Make sure that the base name is sane
    casadi_assert(Function::check_name(this->name));

    // Includes needed
    if (this->main) addInclude("stdio.h");

    // Mex and main need string.h
    if (this->mex || this->main) {
      addInclude("string.h");
    }

    // Memory struct entry point
    if (this->with_mem) {
      addInclude("casadi/mem.h", false);
      this->header << "#include <casadi/mem.h>" << endl;
    }

    // Mex file?
    if (this->mex) {
      addInclude("mex.h", false, "MATLAB_MEX_FILE");
      // Define printf (note file should be compilable both with and without mex)
      this->auxiliaries
        << "#ifdef MATLAB_MEX_FILE" << endl
        << "#define PRINTF mexPrintf" << endl
        << "#else" << endl
        << "#define PRINTF printf" << endl
        << "#endif" << endl;
    } else {
      // Define printf as standard printf from stdio.h
      this->auxiliaries << "#define PRINTF printf" << endl;
    }

    if (this->with_export) {
      this->auxiliaries
        << "#ifndef CASADI_SYMBOL_EXPORT" << endl
        << "#if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)" << endl
        << "#if defined(STATIC_LINKED)" << endl
        << "#define CASADI_SYMBOL_EXPORT" << endl
        << "#else /* defined(STATIC_LINKED) */" << endl
        << "#define CASADI_SYMBOL_EXPORT __declspec(dllexport)" << endl
        << "#endif /* defined(STATIC_LINKED) */" << endl
        << "#elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)" << endl
        << "#define CASADI_SYMBOL_EXPORT __attribute__ ((visibility (\"default\")))" << endl
        << "#else /* defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__) */"  << endl
        << "#define CASADI_SYMBOL_EXPORT" << endl
        << "#endif /* defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__) */" << endl
        << "#endif /* CASADI_SYMBOL_EXPORT */" << endl;
    }
  }

  void CodeGenerator::add(const Function& f) {
    f->generateFunction(*this, f.name(), false);
    if (this->with_header) {
      if (this->cpp) this->header << "extern \"C\" " ; // C linkage
      this->header << f->signature(f.name()) << ";" << endl;
    }
    f->generateMeta(*this, f.name());
    this->exposed_fname.push_back(f.name());
  }

  std::string CodeGenerator::dump() const {
    stringstream s;
    dump(s);
    return s.str();
  }

  void CodeGenerator::file_open(std::ofstream& f, const std::string& name) const {
    // Open a file for writing
    f.open(name);

    // Print header
    f << "/* This file was automatically generated by CasADi.\n"
      << "   The CasADi copyright holders make no ownership claim of its contents. */\n";

    // C linkage
    if (!this->cpp) {
      f << "#ifdef __cplusplus\n"
        << "extern \"C\" {\n"
        << "#endif\n\n";
    }
  }

  void CodeGenerator::file_close(std::ofstream& f) const {
    // C linkage
    if (!this->cpp) {
      f << "#ifdef __cplusplus" << endl
        << "} /* extern \"C\" */" << endl
        << "#endif" << endl;
    }

    // Close file(s)
    f.close();
  }

  void CodeGenerator::generate_real_t(std::ostream &s) const {
    s << "#ifndef real_t" << endl
      << "#define real_t " << this->real_t << endl
      << "#endif /* real_t */" << endl << endl;
  }

  std::string CodeGenerator::generate(const std::string& prefix) const {
    // Throw an error if the prefix contains the filename, since since syntax
    // has changed
    casadi_assert_message(prefix.find(this->name + this->suffix)==string::npos,
       "The signature of CodeGenerator::generate has changed. "
       "Instead of providing the filename, only provide the prefix.");

    // Create c file
    ofstream s;
    string fullname = prefix + this->name + this->suffix;
    file_open(s, fullname);

    // Generate the actual function
    dump(s);

    // Mex entry point
    if (this->mex) generate_mex(s);

    // Main entry point
    if (this->main) generate_main(s);

    // Finalize file
    file_close(s);

    // Generate header
    if (this->with_header) {
      // Create a header file
      file_open(s, prefix + this->name + ".h");

      // Define the real_t type (typically double)
      generate_real_t(s);

      // Add declarations
      s << this->header.str();

      // Finalize file
      file_close(s);
    }
    return fullname;
  }

  void CodeGenerator::generate_mex(std::ostream &s) const {
    // Begin conditional compilation
    s << "#ifdef MATLAB_MEX_FILE" << endl;

    // Function prototype
    if (this->cpp) s << "extern \"C\"" << endl; // C linkage
    s << "void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {"
      << endl;

    // Create a buffer
    size_t buf_len = 0;
    for (int i=0; i<exposed_fname.size(); ++i) {
      buf_len = std::max(buf_len, exposed_fname[i].size());
    }
    s << "  char buf[" << (buf_len+1) << "];" << endl;

    // Read string argument
    s << "  int buf_ok = --argc >= 0 && !mxGetString(*argv++, buf, sizeof(buf));" << endl;

    // Create switch
    s << "  if (!buf_ok) {" << endl
      << "    /* name error */" << endl;
    for (int i=0; i<exposed_fname.size(); ++i) {
      s << "  } else if (strcmp(buf, \"" << exposed_fname[i] << "\")==0) {" << endl
        << "    return mex_" << exposed_fname[i] << "(resc, resv, argc, argv);" << endl;
    }
    s << "  }" << endl;

    // Error
    s << "  mexErrMsgTxt(\"First input should be a command string. Possible values:";
    for (int i=0; i<exposed_fname.size(); ++i) {
      s << " '" << exposed_fname[i] << "'";
    }
    s << "\");" << endl;

    // End conditional compilation and function
    s << "}" << endl
         << "#endif" << endl;
  }

  void CodeGenerator::generate_main(std::ostream &s) const {
    s << "int main(int argc, char* argv[]) {" << endl;

    // Create switch
    s << "  if (argc<2) {" << endl
      << "    /* name error */" << endl;
    for (int i=0; i<exposed_fname.size(); ++i) {
      s << "  } else if (strcmp(argv[1], \"" << exposed_fname[i] << "\")==0) {" << endl
        << "    return main_" << exposed_fname[i] << "(argc-2, argv+2);" << endl;
    }
    s << "  }" << endl;

    // Error
    s << "  fprintf(stderr, \"First input should be a command string. Possible values:";
    for (int i=0; i<exposed_fname.size(); ++i) {
      s << " '" << exposed_fname[i] << "'";
    }
    s << "\\n\");" << endl;

    // End main
    s << "  return 1;" << endl
      << "}" << endl;
  }

  void CodeGenerator::dump(std::ostream& s) const {
    // Consistency check
    casadi_assert(current_indent_ == 0);

    // Prefix internal symbols to avoid symbol collisions
    s << "#ifdef CODEGEN_PREFIX" << endl
         << "  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)" << endl
         << "  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID" << endl
         << "  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)" << endl
         << "#else /* CODEGEN_PREFIX */" << endl
         << "  #define CASADI_PREFIX(ID) " << this->name << "_ ## ID" << endl
         << "#endif /* CODEGEN_PREFIX */" << endl << endl;

    s << this->includes.str();
    s << endl;

    // Real type (usually double)
    generate_real_t(s);

    // Type conversion
    s << "#define to_double(x) "
      << (this->cpp ? "static_cast<double>(x)" : "(double) x") << endl
      << "#define to_int(x) "
      << (this->cpp ? "static_cast<int>(x)" : "(int) x") << endl
      << "#define CASADI_CAST(x,y) "
      << (this->cpp ? "static_cast<x>(y)" : "(x) y") << endl;

    // External function declarations
    if (!added_externals_.empty()) {
      s << "/* External functions */" << endl;
      for (auto&& i : added_externals_) {
        s << i << endl;
      }
      s << endl << endl;
    }

    // Pre-C99
    s << "/* Pre-c99 compatibility */" << endl
         << "#if __STDC_VERSION__ < 199901L" << endl
         << "real_t CASADI_PREFIX(fmin)(real_t x, real_t y) { return x<y ? x : y;}" << endl
         << "#define fmin(x,y) CASADI_PREFIX(fmin)(x,y)" << endl
         << "real_t CASADI_PREFIX(fmax)(real_t x, real_t y) { return x>y ? x : y;}" << endl
         << "#define fmax(x,y) CASADI_PREFIX(fmax)(x,y)" << endl
         << "#endif" << endl << endl;

    // Codegen auxiliary functions
    s << this->auxiliaries.str();

    // Print integer constants
    stringstream name;
    for (int i=0; i<integer_constants_.size(); ++i) {
      name.str(string());
      name << "CASADI_PREFIX(s" << i << ")";
      print_vector(s, name.str(), integer_constants_[i]);
      s << "#define s" << i << " CASADI_PREFIX(s" << i << ")" << endl;
    }

    // Print double constants
    for (int i=0; i<double_constants_.size(); ++i) {
      name.str(string());
      name << "CASADI_PREFIX(c" << i << ")";
      print_vector(s, name.str(), double_constants_[i]);
      s << "#define c" << i << " CASADI_PREFIX(c" << i << ")" << endl;
    }

    // Codegen body
    s << this->body.str();

    // End with new line
    s << endl;
  }

  std::string CodeGenerator::to_string(int n) {
    stringstream ss;
    ss << n;
    return ss.str();
  }

  std::string CodeGenerator::work(int n, int sz) const {
    if (n<0 || sz==0) {
      return "0";
    } else if (sz==1 && !this->codegen_scalars) {
      return "(&w" + to_string(n) + ")";
    } else {
      return "w" + to_string(n);
    }
  }

  std::string CodeGenerator::workel(int n) const {
    if (n<0) return "0";
    stringstream s;
    if (this->codegen_scalars) s << "*";
    s << "w" << n;
    return s.str();
  }

  std::string CodeGenerator::array(const std::string& type, const std::string& name, int len,
                                   const std::string& def) {
    std::stringstream s;
    s << type << " ";
    if (len==0) {
      s << "*" << name << " = 0";
    } else {
      s << name << "[" << len << "]";
      if (!def.empty()) s << " = " << def;
    }
    s << ";" << endl;
    return s.str();
  }

  void CodeGenerator::print_vector(std::ostream &s, const std::string& name, const vector<int>& v) {
    s << array("static const int", name, v.size(), initializer(v));
  }

  void CodeGenerator::print_vector(std::ostream &s, const std::string& name,
                                  const vector<double>& v) {
    s << array("static const real_t", name, v.size(), initializer(v));
  }

  void CodeGenerator::addInclude(const std::string& new_include, bool relative_path,
                                 const std::string& use_ifdef) {
    // Register the new element
    bool added = added_includes_.insert(new_include).second;

    // Quick return if it already exists
    if (!added) return;

    // Ifdef opening
    if (!use_ifdef.empty()) this->includes << "#ifdef " << use_ifdef << endl;

    // Print to the header section
    if (relative_path) {
      this->includes << "#include \"" << new_include << "\"" << endl;
    } else {
      this->includes << "#include <" << new_include << ">" << endl;
    }

    // Ifdef closing
    if (!use_ifdef.empty()) this->includes << "#endif" << endl;
  }

  bool CodeGenerator::simplifiedCall(const Function& f) {
    return f->simplifiedCall();
  }

  std::string CodeGenerator::
  operator()(const Function& f, const std::string& arg,
             const std::string& res, const std::string& iw,
             const std::string& w, const std::string& mem) const {
    return f->codegen_name(*this) + "(" + arg + ", " + res + ", "
      + iw + ", " + w + ", " + mem + ")";
  }

  std::string CodeGenerator::operator()(const Function& f,
                                        const std::string& arg, const std::string& res) const {
    return f->codegen_name(*this) + "(" + arg + ", " + res + ")";
  }

  void CodeGenerator::addExternal(const std::string& new_external) {
    added_externals_.insert(new_external);
  }

  int CodeGenerator::addSparsity(const Sparsity& sp) {
    // Get the current number of patterns before looking for it
    size_t num_patterns_before = added_sparsities_.size();

    // Get index of the pattern
    const void* h = static_cast<const void*>(sp.get());
    int& ind = added_sparsities_[h];

    // Generate it if it does not exist
    if (added_sparsities_.size() > num_patterns_before) {

      // Compact version of the sparsity pattern
      std::vector<int> sp_compact = sp.compress();

      // Codegen vector
      ind = getConstant(sp_compact, true);
    }

    return ind;
  }

  std::string CodeGenerator::sparsity(const Sparsity& sp) {
    return "s" + to_string(addSparsity(sp));
  }

  int CodeGenerator::get_sparsity(const Sparsity& sp) const {
    const void* h = static_cast<const void*>(sp.get());
    PointerMap::const_iterator it=added_sparsities_.find(h);
    casadi_assert(it!=added_sparsities_.end());
    return it->second;
  }

  size_t CodeGenerator::hash(const std::vector<double>& v) {
    // Calculate a hash value for the vector
    std::size_t seed=0;
    if (!v.empty()) {
      casadi_assert(sizeof(double) % sizeof(size_t)==0);
      const int int_len = v.size()*(sizeof(double)/sizeof(size_t));
      const size_t* int_v = reinterpret_cast<const size_t*>(&v.front());
      for (size_t i=0; i<int_len; ++i) {
        hash_combine(seed, int_v[i]);
      }
    }
    return seed;
  }

  size_t CodeGenerator::hash(const std::vector<int>& v) {
    size_t seed=0;
    hash_combine(seed, v);
    return seed;
  }

  int CodeGenerator::getConstant(const std::vector<double>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    auto eq = added_double_constants_.equal_range(h);
    for (auto i=eq.first; i!=eq.second; ++i) {
      if (equal(v, double_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      int ind = double_constants_.size();
      double_constants_.push_back(v);
      added_double_constants_.insert(make_pair(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  int CodeGenerator::getConstant(const std::vector<int>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    pair<multimap<size_t, size_t>::iterator, multimap<size_t, size_t>::iterator> eq =
      added_integer_constants_.equal_range(h);
    for (multimap<size_t, size_t>::iterator i=eq.first; i!=eq.second; ++i) {
      if (equal(v, integer_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      int ind = integer_constants_.size();
      integer_constants_.push_back(v);
      added_integer_constants_.insert(pair<size_t, size_t>(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  std::string CodeGenerator::constant(const std::vector<int>& v) {
    return "s" + to_string(getConstant(v, true));
  }

  std::string CodeGenerator::constant(const std::vector<double>& v) {
    return "c" + to_string(getConstant(v, true));
  }

  void CodeGenerator::addAuxiliary(Auxiliary f, const vector<string>& inst) {
    // Look for existing instantiations
    auto f_match = added_auxiliaries_.equal_range(f);
    // Look for duplicates
    for (auto it=f_match.first; it!=f_match.second; ++it) {
      if (it->second==inst) return;
    }
    added_auxiliaries_.insert(make_pair(f, inst));

    // Add the appropriate function
    switch (f) {
    case AUX_COPY:
      this->auxiliaries << sanitize_source(casadi_copy_str, inst);
      break;
    case AUX_SWAP:
      this->auxiliaries << sanitize_source(casadi_swap_str, inst);
      break;
    case AUX_SCAL:
      this->auxiliaries << sanitize_source(casadi_scal_str, inst);
      break;
    case AUX_AXPY:
      this->auxiliaries << sanitize_source(casadi_axpy_str, inst);
      break;
    case AUX_DOT:
      this->auxiliaries << sanitize_source(casadi_dot_str, inst);
      break;
    case AUX_BILIN:
      this->auxiliaries << sanitize_source(casadi_bilin_str, inst);
      break;
    case AUX_RANK1:
      this->auxiliaries << sanitize_source(casadi_rank1_str, inst);
      break;
    case AUX_IAMAX:
      this->auxiliaries << sanitize_source(casadi_iamax_str, inst);
      break;
    case AUX_INTERPN:
      addAuxiliary(AUX_INTERPN_WEIGHTS);
      addAuxiliary(AUX_INTERPN_INTERPOLATE);
      addAuxiliary(AUX_FLIP, {});
      addAuxiliary(AUX_FILL);
      addAuxiliary(AUX_FILL, {"int"});
      this->auxiliaries << sanitize_source(casadi_interpn_str, inst);
      break;
    case AUX_INTERPN_GRAD:
      addAuxiliary(AUX_INTERPN);
      this->auxiliaries << sanitize_source(casadi_interpn_grad_str, inst);
      break;
    case AUX_DE_BOOR:
      this->auxiliaries << sanitize_source(casadi_de_boor_str, inst);
      break;
    case AUX_ND_BOOR_EVAL:
      addAuxiliary(AUX_DE_BOOR);
      addAuxiliary(AUX_FILL);
      addAuxiliary(AUX_FILL, {"int"});
      addAuxiliary(AUX_LOW);
      this->auxiliaries << sanitize_source(casadi_nd_boor_eval_str, inst);
      break;
    case AUX_FLIP:
      this->auxiliaries << sanitize_source(casadi_flip_str, inst);
      break;
    case AUX_LOW:
      this->auxiliaries << sanitize_source(casadi_low_str, inst);
      break;
    case AUX_INTERPN_WEIGHTS:
      addAuxiliary(AUX_LOW);
      this->auxiliaries << sanitize_source(casadi_interpn_weights_str, inst);
      break;
    case AUX_INTERPN_INTERPOLATE:
      this->auxiliaries << sanitize_source(casadi_interpn_interpolate_str, inst);
      break;
    case AUX_NORM_1:
      this->auxiliaries << sanitize_source(casadi_norm_1_str, inst);
      break;
    case AUX_NORM_2:
      this->auxiliaries << sanitize_source(casadi_norm_2_str, inst);
      break;
    case AUX_NORM_INF:
      this->auxiliaries << sanitize_source(casadi_norm_inf_str, inst);
      break;
    case AUX_FILL:
      this->auxiliaries << sanitize_source(casadi_fill_str, inst);
      break;
    case AUX_MTIMES:
      this->auxiliaries << sanitize_source(casadi_mtimes_str, inst);
      break;
    case AUX_SQ:
      auxSq();
      break;
    case AUX_SIGN:
      auxSign();
      break;
    case AUX_PROJECT:
      this->auxiliaries << sanitize_source(casadi_project_str, inst);
      break;
    case AUX_DENSIFY:
      addAuxiliary(AUX_FILL);
      {
        std::vector<std::string> inst2 = inst;
        if (inst.size()==1) inst2.push_back(inst[0]);
        this->auxiliaries << sanitize_source(casadi_densify_str, inst2);
      }
      break;
    case AUX_TRANS:
      this->auxiliaries << sanitize_source(casadi_trans_str, inst);
      break;
    case AUX_TO_MEX:
      addAuxiliary(AUX_DENSIFY);
      this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                        << sanitize_source(casadi_to_mex_str, inst)
                        << "#endif" << endl << endl;
      break;
    case AUX_FROM_MEX:
      addAuxiliary(AUX_FILL);
      this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                        << sanitize_source(casadi_from_mex_str, inst)
                        << "#endif" << endl << endl;
      break;
    }
  }

  std::string CodeGenerator::to_mex(const Sparsity& sp, const std::string& arg) {
    addAuxiliary(AUX_TO_MEX);
    stringstream s;
    s << "to_mex(" << sparsity(sp) << ", " << arg << ");";
    return s.str();
  }

  std::string CodeGenerator::from_mex(std::string& arg,
                                      const std::string& res, std::size_t res_off,
                                      const Sparsity& sp_res, const std::string& w) {
    // Handle offset with recursion
    if (res_off!=0) return from_mex(arg, res+"+"+to_string(res_off), 0, sp_res, w);

    addAuxiliary(AUX_FROM_MEX);
    stringstream s;
    s << "from_mex(" << arg
      << ", " << res << ", " << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  void CodeGenerator::auxSq() {
    this->auxiliaries
      << "real_t CASADI_PREFIX(sq)(real_t x) "
      << "{ return x*x;}" << endl
      << "#define sq(x) CASADI_PREFIX(sq)(x)" << endl << endl;
  }

  void CodeGenerator::auxSign() {
    this->auxiliaries
      << "real_t CASADI_PREFIX(sign)(real_t x) "
      << "{ return x<0 ? -1 : x>0 ? 1 : x;}" << endl
      << "#define sign(x) CASADI_PREFIX(sign)(x)" << endl << endl;
  }

  std::string CodeGenerator::constant(double v) {
    stringstream s;
    if (isnan(v)) {
      s << "NAN";
    } else if (isinf(v)) {
      if (v<0) s << "-";
      s << "INFINITY";
    } else {
      int v_int(v);
      if (v_int==v) {
        // Print integer
        s << v_int << ".";
      } else {
        // Print real
        std::ios_base::fmtflags fmtfl = s.flags(); // get current format flags
        s << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1) << v;
        s.flags(fmtfl); // reset current format flags
      }
    }
    return s.str();
  }

  std::string CodeGenerator::initializer(const vector<double>& v) {
    stringstream s;
    s << "{";
    for (int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << constant(v[i]);
    }
    s << "}";
    return s.str();
  }

  std::string CodeGenerator::initializer(const vector<int>& v) {
    stringstream s;
    s << "{";
    for (int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << v[i];
    }
    s << "}";
    return s.str();
  }

  std::string CodeGenerator::copy(const std::string& arg,
                                  std::size_t n, const std::string& res) {
    stringstream s;
    // Perform operation
    addAuxiliary(AUX_COPY);
    s << "copy(" << arg << ", " << n << ", " << res << ");";
    return s.str();
  }

  std::string CodeGenerator::fill(const std::string& res,
                                  std::size_t n, const std::string& v) {
    stringstream s;
    // Perform operation
    addAuxiliary(AUX_FILL);
    s << "fill(" << res << ", " << n << ", " << v << ");";
    return s.str();
  }

  std::string CodeGenerator::dot(int n, const std::string& x,
                                        const std::string& y) {
    addAuxiliary(AUX_DOT);
    stringstream s;
    s << "dot(" << n << ", " << x << ", " << y << ")";
    return s.str();
  }

  std::string CodeGenerator::bilin(const std::string& A, const Sparsity& sp_A,
                                   const std::string& x, const std::string& y) {
    addAuxiliary(AUX_BILIN);
    stringstream s;
    s << "bilin(" << A << ", " << sparsity(sp_A) << ", " << x << ", " << y << ")";
    return s.str();
  }

  std::string CodeGenerator::rank1(const std::string& A, const Sparsity& sp_A,
                                   const std::string& alpha, const std::string& x,
                                   const std::string& y) {
    addAuxiliary(AUX_RANK1);
    stringstream s;
    s << "rank1(" << A << ", " << sparsity(sp_A) << ", "
      << alpha << ", " << x << ", " << y << ");";
    return s.str();
  }

  std::string CodeGenerator::interpn(int ndim, const std::string& grid, const std::string& offset,
                                   const std::string& values, const std::string& x,
                                   const std::string& lookup_mode,
                                   const std::string& iw, const std::string& w) {
    addAuxiliary(AUX_INTERPN);
    stringstream s;
    s << "interpn(" << ndim << ", " << grid << ", "  << offset << ", "
      << values << ", " << x << ", " << lookup_mode << ", " << iw << ", " << w << ");";
    return s.str();
  }

  std::string CodeGenerator::interpn_grad(const std::string& grad,
                                   int ndim, const std::string& grid, const std::string& offset,
                                   const std::string& values, const std::string& x,
                                   const std::string& lookup_mode,
                                   const std::string& iw, const std::string& w) {
    addAuxiliary(AUX_INTERPN_GRAD);
    stringstream s;
    s << "interpn_grad(" << grad << ", " << ndim << ", " << grid << ", " << offset << ", "
      << values << ", " << x << ", " << lookup_mode << ", " << iw << ", " << w << ");";
    return s.str();
  }

  std::string CodeGenerator::declare(std::string s) {
    // Add C linkage?
    if (this->cpp) {
      s = "extern \"C\" " + s;
    }

    // To header file
    if (this->with_header) {
      this->header << s << ";" << endl;
    }

    return s;
  }

  std::string
  CodeGenerator::project(const std::string& arg, const Sparsity& sp_arg,
                         const std::string& res, const Sparsity& sp_res,
                         const std::string& w) {
    // If sparsity match, simple copy
    if (sp_arg==sp_res) return copy(arg, sp_arg.nnz(), res);

    // Create call
    addAuxiliary(AUX_PROJECT);
    stringstream s;
    s << "project(" << arg << ", " << sparsity(sp_arg) << ", " << res << ", "
      << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  std::string CodeGenerator::printf(const std::string& str, const std::vector<std::string>& arg) {
    addInclude("stdio.h");
    stringstream s;
    s << "PRINTF(\"" << str << "\"";
    for (int i=0; i<arg.size(); ++i) s << ", " << arg[i];
    s << ");";
    return s.str();
  }

  std::string CodeGenerator::printf(const std::string& str, const std::string& arg1) {
    std::vector<std::string> arg;
    arg.push_back(arg1);
    return printf(str, arg);
  }

  std::string CodeGenerator::printf(const std::string& str, const std::string& arg1,
                                    const std::string& arg2) {
    std::vector<std::string> arg;
    arg.push_back(arg1);
    arg.push_back(arg2);
    return printf(str, arg);
  }

  std::string CodeGenerator::printf(const std::string& str, const std::string& arg1,
                                    const std::string& arg2, const std::string& arg3) {
    std::vector<std::string> arg;
    arg.push_back(arg1);
    arg.push_back(arg2);
    arg.push_back(arg3);
    return printf(str, arg);
  }

  std::string CodeGenerator::mtimes(const std::string& x, const Sparsity& sp_x,
                                    const std::string& y, const Sparsity& sp_y,
                                    const std::string& z, const Sparsity& sp_z,
                                    const std::string& w, bool tr) {
    addAuxiliary(AUX_MTIMES);
    stringstream s;
    s << "mtimes(" << x << ", " << sparsity(sp_x) << ", " << y << ", " << sparsity(sp_y) << ", "
      << z << ", " << sparsity(sp_z) << ", " << w << ", " <<  (tr ? "1" : "0") << ");";
    return s.str();
  }

  void CodeGenerator::print_formatted(const std::string& s) {
    // Quick return if empty
    if (s.empty()) return;

    // If new line, add indentation
    if (newline_) {
      int shift = s.front()=='}' ? -1 : 0;
      casadi_assert(current_indent_+shift>=0);
      this->buffer << string(indent_*(current_indent_+shift), ' ');
      newline_ = false;
    }

    // Print to body
    this->buffer << s;

    // Brackets change indentation for next row
    // NOTE(@jaeandersson): Should ignore strings, comments
    for (char c : s) {
      if (c=='{') {
        indent();
      } else if (c=='}') {
        unindent();
      }
    }
  }

  CodeGenerator& CodeGenerator::operator<<(const string& s) {
    // Loop over newline characters
    size_t off=0;
    while (true) {
      size_t pos = s.find('\n', off);
      if (pos==string::npos) {
        // No more newline characters
        print_formatted(s.substr(off));
        break;
      } else {
        // Ends with newline
        print_formatted(s.substr(off, pos-off));
        this->buffer << '\n';
        newline_ = true;
        off = pos+1;
      }
    }

    return *this;
  }

  void CodeGenerator::flush(std::ostream &s) {
    s << this->buffer.str();
    this->buffer.str(string());
  }

  void CodeGenerator::local(const std::string& name, const std::string& type,
                            const std::string& ref) {
    // Check if the variable already exists
    auto it = local_variables_.find(name);
    if (it==local_variables_.end()) {
      // Add it
      local_variables_[name] = make_pair(type, ref);
    } else {
      // Consistency check
      casadi_assert_message(it->second.first==type, "Type mismatch for " + name);
      casadi_assert_message(it->second.second==ref, "Type mismatch for " + name);
    }
  }

  void CodeGenerator::init_local(const string& name, const string& def) {
    bool inserted = local_default_.insert(make_pair(name, def)).second;
    casadi_assert_message(inserted, name + " already defined");
  }

  string CodeGenerator::
  sanitize_source(const std::string& src, const vector<string>& inst, bool shorthand) {
    // Create suffix if templates type are not all "real_t"
    string suffix;
    for (const string& s : inst) {
      if (s!="real_t") {
        for (const string& s : inst) suffix += "_" + s;
        break;
      }
    }

    // Construct map of name replacements
    std::map<std::string, std::string> rep;
    for (int i=0; i<inst.size(); ++i) rep["T" + to_string(i+1)] = inst[i];
    // Return object
    stringstream ret;
    // Number of template parameters
    size_t npar = 0;
    // Process C++ source
    string line, def, fname;
    istringstream stream(src);
    while (std::getline(stream, line)) {
      size_t n1, n2;

      // C++ template declaration
      if (line.find("template")==0) {
        npar = count(line.begin(), line.end(), ',') + 1;
        continue;
      }

      // Inline declaration
      if (line == "inline") {
        continue;
      }

      // Ignore C++ style comment
      if ((n1 = line.find("//")) != string::npos) line.erase(n1);

      // Remove trailing spaces
      if ((n1 = line.find_last_not_of(' ')) != string::npos) {
        line.erase(n1 + 1);
      } else {
        continue;
      }

      // Get function name (must be the first occurrence of CASADI_PREFIX)
      if (fname.empty()) {
        string s = "CASADI_PREFIX(";
        if ((n1 = line.find(s)) == string::npos) {
          casadi_error("Cannot find \"" + s + "\" in " + line);
        }
        n1 += s.size();
        s = ")(";
        if ((n2=line.find(s, n1)) == string::npos) {
          casadi_error("Cannot find \"" + s + "\" in " + line);
        }
        fname = line.substr(n1, n2-n1);
        casadi_assert(!fname.empty());
        for (char c : fname) {
          casadi_assert_message(isalnum(c) || c=='_', "Invalid filename: " + fname);
        }

        // Get argument list
        n1 = n2 + s.size();
        n2 = line.find(")", n1);

        casadi_assert(n2!=string::npos);
        string args = line.substr(n1, n2-n1) + ",";

        // Get argument list
        n1 = 0;
        while ((n2 = args.find(',', n1)) != string::npos) {
          n1 = args.rfind(' ', n2);
          s = args.substr(n1+1, n2-n1-1);
          def = def.empty() ? s : def + ", " + s;
          n1 = n2 + 1;

        }

        // Add suffix
        if (!suffix.empty()) {
          line.replace(line.find(fname), fname.size(), fname + suffix);
          fname += suffix;
        }

        // Finalize shorthand, e.g. #define fmin(x,y) CASADI_PREFIX(fmin)(x,y)
        def = "#define " + fname + "(" + def + ") CASADI_PREFIX(" + fname + ")(" + def + ")\n";
      }

      // Perform string replacements
      for (auto&& r : rep) {
        string::size_type n = 0;
        while ((n = line.find(r.first, n)) != string::npos) {
          line.replace(n, r.first.size(), r.second);
          n += r.second.size();
        }
      }

      // Append to return
      ret << line << "\n";
    }

    // Assert number of template parameters
    casadi_assert_message(npar==inst.size(),
      "Mismatching number of template parameters for " + fname);

    // Add shorthand
    if (shorthand) ret << def;

    // Trailing newline
    ret << "\n";
    return ret.str();
  }

} // namespace casadi
