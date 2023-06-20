/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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
#include "convexify.hpp"
#include <casadi_runtime_str.h>
#include <iomanip>

namespace casadi {

  CodeGenerator::CodeGenerator(const std::string& name, const Dict& opts) {
    // Default options
    this->verbose = true;
    this->verbose_runtime = false;
    this->mex = false;
    this->with_sfunction = false;
    this->cpp = false;
    this->main = false;
    this->casadi_real_type = "double";
    this->casadi_int_type = CASADI_INT_TYPE_STR;
    this->codegen_scalars = false;
    this->with_header = false;
    this->with_mem = false;
    this->with_export = true;
    this->with_import = false;
    this->include_math = true;
    this->infinity = "INFINITY";
    this->nan = "NAN";
    this->real_min = "";
    bool prefix_set = false;
    this->prefix = "";
    avoid_stack_ = false;
    indent_ = 2;

    // Read options
    for (auto&& e : opts) {
      if (e.first=="verbose") {
        this->verbose = e.second;
      } else if (e.first=="verbose_runtime") {
        this->verbose_runtime = e.second;
      } else if (e.first=="mex") {
        this->mex = e.second;
      } else if (e.first=="with_sfunction") {
        this->with_sfunction = e.second;
      } else if (e.first=="cpp") {
        this->cpp = e.second;
      } else if (e.first=="main") {
        this->main = e.second;
      } else if (e.first=="casadi_real") {
        this->casadi_real_type = e.second.to_string();
      } else if (e.first=="casadi_int") {
        this->casadi_int_type = e.second.to_string();
      } else if (e.first=="codegen_scalars") {
        this->codegen_scalars = e.second;
      } else if (e.first=="with_header") {
        this->with_header = e.second;
      } else if (e.first=="with_mem") {
        this->with_mem = e.second;
      } else if (e.first=="with_export") {
        this->with_export = e.second;
      } else if (e.first=="with_import") {
        this->with_import = e.second;
      } else if (e.first=="include_math") {
        this->include_math = e.second;
      } else if (e.first=="infinity") {
        this->infinity = e.second.to_string();
      } else if (e.first=="nan") {
        this->nan = e.second.to_string();
      } else if (e.first=="real_min") {
        this->real_min = e.second.to_string();
      } else if (e.first=="indent") {
        indent_ = e.second;
        casadi_assert_dev(indent_>=0);
      } else if (e.first=="avoid_stack") {
        avoid_stack_ = e.second;
      } else if (e.first=="prefix") {
        this->prefix = e.second.to_string();
        prefix_set = true;
      } else {
        casadi_error("Unrecognized option: " + str(e.first));
      }
    }

    // If real_min is not specified, make an educated guess
    if (this->real_min.empty()) {
      std::stringstream ss;
      ss << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
      if (casadi_real_type=="float") {
        ss << std::numeric_limits<float>::min();
        this->real_min = ss.str();
      } else if (casadi_real_type=="double") {
        ss << std::numeric_limits<double>::min();
        this->real_min = ss.str();
      } else {
        this->real_min = "<NOT SPECIFIED>";
      }
    }

    // Start at new line with no indentation
    newline_ = true;
    current_indent_ = 0;

    // Start off without the need for thread-local memory
    needs_mem_ = false;

    // Divide name into base and suffix (if any)
    std::string::size_type dotpos = name.rfind('.');
    if (dotpos==std::string::npos) {
      this->name = name;
      this->suffix = this->cpp ? ".cpp" : ".c";
    } else {
      this->name = name.substr(0, dotpos);
      this->suffix = name.substr(dotpos);
    }

    // Symbol prefix
    if (this->with_export) dll_export = "CASADI_SYMBOL_EXPORT ";
    if (this->with_import) dll_import = "CASADI_SYMBOL_IMPORT ";

    // Make sure that the base name is sane
    casadi_assert_dev(Function::check_name(this->name));

    // Includes needed
    if (this->include_math) add_include("math.h");
    if (this->main) add_include("stdio.h");
    if (this->verbose_runtime) add_auxiliary(AUX_PRINTF);

    // Mex and main need string.h
    if (this->mex || this->main) {
      add_include("string.h");
    }

    // Mex
    if (this->mex) {
      add_include("mex.h", false, "MATLAB_MEX_FILE");
    }

    // s-Function
    if (this->with_sfunction) {
      this->casadi_real_type = "real_T";
      this->casadi_int_type = "int_T";
      this->with_header = true;
      add_include("simstruc.h");
    }

    // Memory struct entry point
    if (this->with_mem) {
      this->header << "#include <casadi/mem.h>\n";
    }

    // Use name as default prefix
    if (!prefix_set) {
      this->prefix = this->name;
    }

  }

  void CodeGenerator::scope_enter() {
    local_variables_.clear();
    local_default_.clear();
  }

  void CodeGenerator::scope_exit() {
    // Order local variables
    std::map<std::string, std::set<std::pair<std::string, std::string>>> local_variables_by_type;
    for (auto&& e : local_variables_) {
      local_variables_by_type[e.second.first].insert(std::make_pair(e.first, e.second.second));
    }

    // Codegen local variables
    for (auto&& e : local_variables_by_type) {
      body << "  " << e.first;
      for (auto it=e.second.begin(); it!=e.second.end(); ++it) {
        body << (it==e.second.begin() ? " " : ", ") << it->second << it->first;
        // Insert definition, if any
        auto k=local_default_.find(it->first);
        if (k!=local_default_.end()) body << "=" << k->second;
      }
      body << ";\n";
    }
  }

  std::string CodeGenerator::add_dependency(const Function& f) {
    // Quick return if it already exists
    for (auto&& e : added_functions_) if (e.f==f) return e.codegen_name;

    // Give it a name
    std::string fname = shorthand("f" + str(added_functions_.size()));

    // Add to list of functions
    added_functions_.push_back({f, fname});

    // Generate declarations
    f->codegen_declarations(*this);

    // Print to file
    f->codegen(*this, fname);

    // Codegen reference count functions, if needed
    if (f->has_refcount_) {
      // Increase reference counter
      *this << "void " << fname << "_incref(void) {\n";
      f->codegen_incref(*this);
      *this << "}\n\n";

      // Decrease reference counter
      *this << "void " << fname << "_decref(void) {\n";
      f->codegen_decref(*this);
      *this << "}\n\n";
    }

    bool fun_needs_mem = !f->codegen_mem_type().empty();
    needs_mem_ |= fun_needs_mem;

    if (fun_needs_mem) {
      // Alloc memory
      *this << "int " << fname << "_alloc_mem(void) {\n";
      flush(this->body);
      scope_enter();
      f->codegen_alloc_mem(*this);
      scope_exit();
      *this << "}\n\n";

      // Initialize memory
      *this << "int " << fname << "_init_mem(int mem) {\n";
      flush(this->body);
      scope_enter();
      f->codegen_init_mem(*this);
      scope_exit();
      *this << "}\n\n";

      // Clear memory
      *this << "void " << fname << "_free_mem(int mem) {\n";
      flush(this->body);
      scope_enter();
      f->codegen_free_mem(*this);
      scope_exit();
      *this << "}\n\n";

      // Checkout
      *this << "int " << fname << "_checkout(void) {\n";
      flush(this->body);
      scope_enter();
      f->codegen_checkout(*this);
      scope_exit();
      *this << "}\n\n";

      // Clear memory
      *this << "void " << fname << "_release(int mem) {\n";
      flush(this->body);
      scope_enter();
      f->codegen_release(*this);
      scope_exit();
      *this << "}\n\n";

    }

    // Flush to body
    flush(this->body);

    return fname;
  }

    void CodeGenerator::add(const Function& f, bool with_jac_sparsity) {
    // Add if not already added
    std::string codegen_name = add_dependency(f);

    // Define function
    *this << declare(f->signature(f.name())) << "{\n"
          << "return " << codegen_name <<  "(arg, res, iw, w, mem);\n"
          << "}\n\n";

    // Generate meta information
    f->codegen_meta(*this);

    // Generate Jacobian sparsity information
    if (with_jac_sparsity) {
      // Generate/get Jacobian sparsity
      std::vector<Sparsity> jac = f.jac_sparsity();
      // Code generate the sparsity pattern
      add_io_sparsities("jac_" + f.name(), f->sparsity_in_, jac);

      // Flush buffers
      flush(this->body);
    }

    // Generate function specific code for Simulink sfunction
    if (this->with_sfunction) this->added_sfunctions.push_back( this->codegen_sfunction(f) );

    // Add to list of exposed symbols
    this->exposed_fname.push_back(f.name());
  }

  std::string CodeGenerator::dump() {
    std::stringstream s;
    dump(s);
    return s.str();
  }

  std::string CodeGenerator::casadi_version() {
    std::stringstream ss;
    ss << CASADI_MAJOR_VERSION << "." << CASADI_MINOR_VERSION << "." << CASADI_PATCH_VERSION;
    if (!CASADI_IS_RELEASE) ss << "+";
    return ss.str();
  }

  void CodeGenerator::file_open(std::ofstream& f, const std::string& name, bool cpp) {
    // Open a file for writing
    f.open(name);

    // Print header
    f << "/* This file was automatically generated by CasADi " << casadi_version() << ".\n"
      << " *  It consists of: \n"
      << " *   1) content generated by CasADi runtime: not copyrighted\n"
      << " *   2) template code copied from CasADi source: permissively licensed (MIT-0)\n"
      << " *   3) user code: owned by the user\n"
      << " *\n"
      << " */\n";

    // C linkage
    if (!cpp) {
      f << "#ifdef __cplusplus\n"
        << "extern \"C\" {\n"
        << "#endif\n\n";
    }
  }

  void CodeGenerator::file_close(std::ofstream& f, bool cpp) {
    // C linkage
    if (!cpp) {
      f << "#ifdef __cplusplus\n"
        << "} /* extern \"C\" */\n"
        << "#endif\n";
    }

    // Close file(s)
    f.close();
  }

  void CodeGenerator::generate_casadi_real(std::ostream &s) const {
    s << "#ifndef casadi_real\n"
      << "#define casadi_real " << this->casadi_real_type << std::endl
      << "#endif\n\n";
  }

  void  CodeGenerator::generate_export_symbol(std::ostream &s) const {
      s << "/* Symbol visibility in DLLs */\n"
      << "#ifndef CASADI_SYMBOL_EXPORT\n"
      << "  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)\n"
      << "    #if defined(STATIC_LINKED)\n"
      << "      #define CASADI_SYMBOL_EXPORT\n"
      << "    #else\n"
      << "      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)\n"
      << "    #endif\n"
      << "  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)\n"
      << "    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility (\"default\")))\n"
      << "  #else"  << std::endl
      << "    #define CASADI_SYMBOL_EXPORT\n"
      << "  #endif\n"
      << "#endif\n\n";
  }

  void  CodeGenerator::generate_import_symbol(std::ostream &s) const {
      s << "/* Symbol visibility in DLLs */\n"
      << "#ifndef CASADI_SYMBOL_IMPORT\n"
      << "  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)\n"
      << "    #if defined(STATIC_LINKED)\n"
      << "      #define CASADI_SYMBOL_IMPORT\n"
      << "    #else\n"
      << "      #define CASADI_SYMBOL_IMPORT __declspec(dllimport)\n"
      << "    #endif\n"
      << "  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)\n"
      << "    #define CASADI_SYMBOL_IMPORT __attribute__ ((visibility (\"default\")))\n"
      << "  #else"  << std::endl
      << "    #define CASADI_SYMBOL_IMPORT\n"
      << "  #endif\n"
      << "#endif\n\n";
  }

  void CodeGenerator::generate_casadi_int(std::ostream &s) const {
    s << "#ifndef casadi_int\n"
      << "#define casadi_int " << this->casadi_int_type << std::endl
      << "#endif\n\n";
  }

  std::string CodeGenerator::generate(const std::string& prefix) {
    // Throw an error if the prefix contains the filename, since since syntax
    // has changed
    casadi_assert(prefix.find(this->name + this->suffix)==std::string::npos,
       "The signature of CodeGenerator::generate has changed. "
       "Instead of providing the filename, only provide the prefix.");

    // Create c file
    std::ofstream s;
    std::string fullname = prefix + this->name + this->suffix;
    file_open(s, fullname, this->cpp);

    // Dump code to file
    dump(s);

    // Mex entry point
    if (this->mex) generate_mex(s);

    // Main entry point
    if (this->main) generate_main(s);

    // Finalize file
    file_close(s, this->cpp);

    // Generate s-function
    if (this->with_sfunction) {
      for (unsigned ii=0; ii<this->added_sfunctions.size(); ii++) {
        std::string sfunction_code = this->added_sfunctions.at(ii);
        std::string sfunction_name = this->exposed_fname.at(ii);
        generate_sfunction(sfunction_name, sfunction_code);
      }
    }

    // Generate header
    if (this->with_header) {
      // Create a header file
      file_open(s, prefix + this->name + ".h", this->cpp);

      // Define the casadi_real type (typically double)
      generate_casadi_real(s);

      // Define the casadi_int type
      generate_casadi_int(s);

      // Generate export symbol macros
      if (this->with_import) generate_import_symbol(s);

      // Add declarations
      s << this->header.str();

      // Finalize file
      file_close(s, this->cpp);
    }
    return fullname;
  }

  void CodeGenerator::generate_mex(std::ostream &s) const {
    // Begin conditional compilation
    s << "#ifdef MATLAB_MEX_FILE\n";

    // Function prototype
    if (this->cpp) s << "extern \"C\"\n"; // C linkage
    s << "void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {"
      << std::endl;

    // Create a buffer
    size_t buf_len = 0;
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      buf_len = std::max(buf_len, exposed_fname[i].size());
    }
    s << "  char buf[" << (buf_len+1) << "];\n";

    // Read std::string argument
    s << "  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));\n";

    // Create switch
    s << "  if (!buf_ok) {\n";
    // Allow stringless call when unambiguous
    if (exposed_fname.size()==1) {
      s << "    mex_" << exposed_fname[0] << "(resc, resv, argc, argv);\n"
        << "    return;\n";
    } else {
      s << "    /* name error */\n";
    }
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << "  } else if (strcmp(buf, \"" << exposed_fname[i] << "\")==0) {\n"
        << "    mex_" << exposed_fname[i] << "(resc, resv, argc-1, argv+1);\n"
        << "    return;\n";
    }
    s << "  }\n";

    // Error
    s << "  mexErrMsgTxt(\"First input should be a command string. Possible values:";
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << " '" << exposed_fname[i] << "'";
    }
    s << "\");\n";

    // End conditional compilation and function
    s << "}\n"
         << "#endif\n";
  }

  void CodeGenerator::generate_sfunction(const std::string& name,
      const std::string& sfunction) const {
    // Create c file
    std::ofstream f;
    f.open("sfun_"+ name + ".c");

    // Print header
    f << "// Must specify the S_FUNCTION_NAME as the name of the S-function\n"
      << "#define S_FUNCTION_NAME sfun_" << name << "\n"
      << "#define S_FUNCTION_LEVEL 2\n\n"
      << "// Need to include simstruc.h for the definition of the SimStruct and its\n"
      << "// associated macro definitions\n"
      << "#ifndef __SIMSTRUC__\n"
      << "#include \"simstruc.h\"\n"
      << "#endif\n\n"
      << "// Specific header file(s) required by the legacy code function\n"
      << "#include \"" << this->name << ".h\"\n\n\n";

    // Codegenerate s-function
    f << sfunction;

    // Close file(s)
    f.close();
  }

  std::string CodeGenerator::codegen_sfunction(const Function& f) const {
    std::stringstream g;
    // TODO(@jaeandersson): These helper functions really should be moved
    // to the runtime directory

    // Initialize function
    g << "/* Function: mdlInitializeSizes ===========================================\n"
      << "* Abstract:\n"
      << "*   The sizes information is used by Simulink to determine the S-function\n"
      << "*   blocks characteristics (number of inputs, outputs, states, etc.).\n"
      << "*/\n"
      << "static void mdlInitializeSizes(SimStruct *S)\n"
      << "{\n\n"
      << "  /* Declare auxilary variables */\n"
      << "  int_T ii;\n"
      << "  const int_T* sp;\n\n"
      << "  /* Set number of simulink s-function block parameters "
        "(the ones which appear by double click on simulink block) */\n"
      << "  ssSetNumSFcnParams(S, 0);\n\n"
      << "  /* Report if parameter mismatch occurs */\n"
      << "  if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;\n\n"
      << "  /* Specify the number of states for which a block detects "
        "zero crossings that occur between sample points */\n"
      << "  ssSetNumNonsampledZCs(S, 0);\n\n"
      << "  /* Set number of simulink input ports */\n"
      << "  if (!ssSetNumInputPorts(S, " << f->n_in_ << ")) return;\n\n"
      << "  /* Configure simulink input ports (inputs are assumed to be dense "
        "vectors or matrices) */\n"
      << "  for (ii=0; ii<" << f->n_in_ << "; ++ii) {\n"
      << "    sp = " << f.name() << "_sparsity_in(ii);\n"
      << "    if (sp[1]==1) {\n"
      << "      ssSetInputPortWidth(S, ii, sp[0]);\n"
      << "    }\n"
      << "    else {\n"
      << "      ssSetInputPortMatrixDimensions(S, ii, sp[0], sp[1]);\n"
      << "    }\n"
      << "    ssSetInputPortDirectFeedThrough(S, ii, 1);\n"
      << "  }\n\n"
      << "  /* Set number of simulink output ports */\n"
      << "  if (!ssSetNumOutputPorts(S, " << f->n_out_ << ")) return;\n\n"
      << "  /* Configure simulink output ports (dense or sparse vectors or matrices allowed) */\n"
      << "  for (ii=0; ii<" << f->n_out_ << "; ++ii) {\n"
      << "    sp = " << f.name() << "_sparsity_out(ii);\n"
      << "    if (sp[1]==1) {\n"
      << "      ssSetOutputPortWidth(S, ii, sp[0]);\n"
      << "    }\n"
      << "    else {\n"
      << "      ssSetOutputPortMatrixDimensions(S, ii, sp[0], sp[1]);\n"
      << "    }\n"
      << "  }\n"
      << "  ssSetOutputPortOutputExprInRTW(S, 0, 0);\n\n"
      << "  /* This S-function can be used in referenced model simulating in normal mode */\n"
      << "  ssSetModelReferenceNormalModeSupport(S, MDL_START_AND_MDL_PROCESS_PARAMS_OK);\n\n"
      << "  /* Set the number of sample time */\n"
      << "  ssSetNumSampleTimes(S, 1);\n\n"
      << "  /* Set the compliance with the SimState feature */\n"
      << "  ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);\n\n"
      << "  /**\n"
      << "  * All options have the form SS_OPTION_<name> and are documented in\n"
      << "    * matlabroot/simulink/include/simstruc.h. The options should be\n"
      << "    * bitwise ord together as in\n"
      << "    *    ssSetOptions(S, (SS_OPTION_name1 | SS_OPTION_name2))\n"
      << "    */\n"
      << "}\n\n\n";

    // Initialize sample times function
    g << "/* Function: mdlInitializeSampleTimes =====================================\n"
      << " * Abstract:\n"
      << " *   This function is used to specify the sample time(s) for your\n"
      << " *   S-function. You must register the same number of sample times as\n"
      << " *   specified in ssSetNumSampleTimes.\n"
      << " */\n"
      << "static void mdlInitializeSampleTimes(SimStruct *S)\n"
      << "{\n"
      << "    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);\n"
      << "    ssSetOffsetTime(S, 0, FIXED_IN_MINOR_STEP_OFFSET);\n"
      << "    #if defined(ssSetModelReferenceSampleTimeDefaultInheritance)\n"
      << "    ssSetModelReferenceSampleTimeDefaultInheritance(S);\n"
      << "    #endif\n"
      << "}\n\n\n";

    // Model output function
    g << "/* Function: mdlOutputs ===================================================\n"
      << " * Abstract:\n"
      << " *   In this function, you compute the outputs of your S-function\n"
      << " *   block. Generally outputs are placed in the output vector(s),\n"
      << " *   ssGetOutputPortSignal.\n"
      << " */\n"
      << "static void mdlOutputs(SimStruct *S, int_T tid)\n"
      << "{\n\n"
      << "  /* Declare auxilary variables */\n"
      << "  int_T ii, jj, row, col, nnz_col, ind_start_row_index, offset = 0, jj_total = 0;\n"
      << "  const int_T* sp;\n\n"
      << "  /* Allocate buffers for casadi input and output and simulink output */\n"
      << "  " + CodeGenerator::array("real_T", "w", f->sz_w()+f->nnz_out())
      << "  " + CodeGenerator::array("int_T", "iw", f->sz_iw())
      << "  const real_T* arg[" << f->sz_arg() <<"] = {0};\n"
      << "  real_T* res[" << f->sz_res() << "] = {0};\n"
      << "  real_T* y[" << f->n_out_ << "] = {0};\n\n"
      << "  /* Point inputs directly to casadi input buffer */\n"
      << "  for (ii=0; ii<" << f->n_in_ << ";++ii) {\n"
      << "    arg[ii] = *ssGetInputPortRealSignalPtrs( S, ii );\n"
      << "  }\n\n"
      << "  /* Point outputs to buffer */\n"
      << "  for (ii=0; ii<" << f->n_out_ << ";++ii) {\n"
      << "    y[ii] = ssGetOutputPortRealSignal( S, ii );\n"
      << "  }\n\n"
      << "  /* Point allocated working array to casadi output buffer */\n";
    for (casadi_int ii=0; ii<f->n_out_; ++ii) {
      g << "  res[" << ii << "] = w + offset;\n"
        << "  offset += " << f.nnz_out(ii) << ";\n";
    }
    g << "  \n"
      << "  /* Call CasADi function */\n"
      << "  " << f.name() << "( arg, res, iw, w+offset, 0 );\n\n"
      << "  /* Assign results to Simulink output array */\n"
      << "  for (ii=0; ii<" << f->n_out_ << "; ++ii){\n\n"
      << "    /* Get sparsity information of casadi function output "
        "(sp[0] - n_rows, sp[1] - n_cols, sp[2] - dense/sparse) */\n"
      << "    sp = " << f.name() << "_sparsity_out(ii);\n\n"
      << "    /* Check if output is dense (sp[2]=1) or sparse (sp[2]=0) */\n"
      << "    if (sp[2]==0) {\n"
      << "      jj_total = 0;\n"
      << "      ind_start_row_index = 2 + sp[1] + 1;\n\n"
      << "      /* Distribute nonzero elements column by column */\n"
      << "      for (col=0; col<sp[1]; col++) {\n\n"
      << "        /* The cumulative sum of nonzero elements after each column starts at index 2, "
        "after last entry of CCS array col_ptr; number of nonzero elements in current column is "
        "obtained by the difference of two consecutive values */\n"
      << "        nnz_col = sp[2+col+1] - sp[2+col];\n\n"
      << "        /* Distribute nonzero elements of current column to correct row position */\n"
      << "        for (jj=0; jj<nnz_col; jj++) {\n"
      << "          row = sp[ind_start_row_index+jj_total];\n"
      << "          y[ii][row + sp[0]*col] = res[ii][jj_total];\n"
      << "          jj_total++;\n"
      << "        }\n"
      << "      }\n"
      << "    }\n"
      << "    else {\n"
      << "      y[ii] = res[ii];\n"
      << "    }\n"
      << "  }\n"
      << "}\n\n\n";

    // Model terminate function
    g << "/* Function: mdlTerminate =================================================\n"
      << " * Abstract:\n"
      << " *   In this function, you should perform any actions that are necessary\n"
      << " *   at the termination of a simulation.\n"
      << " */\n"
      << "static void mdlTerminate(SimStruct *S)\n"
      << "{\n"
      << "}\n\n\n"
      << "/* Required S-function trailer */\n"
      << "#ifdef MATLAB_MEX_FILE\n"
      << "# include \"simulink.c\"\n"
      << "#else\n"
      << "# include \"cg_sfun.h\"\n"
      << "#endif";

    return g.str();
  }

  void CodeGenerator::generate_main(std::ostream &s) const {
    s << this->dll_export << "int main(int argc, char* argv[]) {\n";

    // Create switch
    s << "  if (argc<2) {\n"
      << "    /* name error */\n";
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << "  } else if (strcmp(argv[1], \"" << exposed_fname[i] << "\")==0) {\n"
        << "    return main_" << exposed_fname[i] << "(argc-2, argv+2);\n";
    }
    s << "  }\n";

    // Error
    s << "  fprintf(stderr, \"First input should be a command string. Possible values:";
    for (casadi_int i=0; i<exposed_fname.size(); ++i) {
      s << " '" << exposed_fname[i] << "'";
    }
    s << "\\n";
    s << "Note: you may use function.generate_input to create a command string.";
    s << "\\n\");\n";

    // End main
    s << "  return 1;\n"
      << "}\n";
  }

  void CodeGenerator::define_rom_double(const void* id, casadi_int size) {
    auto it = file_scope_double_.find(id);
    casadi_assert(it==file_scope_double_.end(), "Already defined.");
    shorthand("rd" + str(file_scope_double_.size()));
    file_scope_double_[id] = size;
  }

  std::string CodeGenerator::rom_double(const void* id) const {
    auto it = file_scope_double_.find(id);
    casadi_assert(it!=file_scope_double_.end(), "Not defined.");
    casadi_int size = std::distance(file_scope_double_.begin(), it);
    return "casadi_rd" + str(size);
  }

  void CodeGenerator::define_rom_integer(const void* id, casadi_int size) {
    auto it = file_scope_double_.find(id);
    casadi_assert(it==file_scope_double_.end(), "Already defined.");
    shorthand("ri" + str(file_scope_double_.size()));
    file_scope_double_[id] = size;
  }

  std::string CodeGenerator::rom_integer(const void* id) const {
    auto it = file_scope_double_.find(id);
    casadi_assert(it!=file_scope_double_.end(), "Not defined.");
    casadi_int size = std::distance(file_scope_double_.begin(), it);
    return "casadi_ri" + str(size);
  }

  void CodeGenerator::dump(std::ostream& s) {
    // Consistency check
    casadi_assert_dev(current_indent_ == 0);

    // Prefix internal symbols to avoid symbol collisions
    s << "/* How to prefix internal symbols */\n"
      << "#ifdef CASADI_CODEGEN_PREFIX\n"
      << "  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)\n"
      << "  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID\n"
      << "  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)\n"
      << "#else\n"
      << "  #define CASADI_PREFIX(ID) " << this->prefix << "_ ## ID\n"
      << "#endif\n\n";

    s << this->includes.str();
    s << std::endl;

    // Numeric types after includes: may depend on them. e.g. mex type
    // Real type (usually double)
    generate_casadi_real(s);

    // Integer type (usually long long)
    generate_casadi_int(s);

    if (needs_mem_) {
      s << "#ifndef CASADI_MAX_NUM_THREADS\n";
      s << "#define CASADI_MAX_NUM_THREADS 1\n";
      s << "#endif\n\n";
    }

    // casadi/mem after numeric types to define derived types
    // Memory struct entry point
    if (this->with_mem) {
      s << "#include <casadi/mem.h>\n" << std::endl;
    }

    // Macros
    if (!added_shorthands_.empty()) {
      s << "/* Add prefix to internal symbols */\n";
      for (auto&& i : added_shorthands_) {
        s << "#define " << "casadi_" << i <<  " CASADI_PREFIX(" << i <<  ")\n";
      }
      s << std::endl;
    }

    if (this->with_export) generate_export_symbol(s);

    // Check if inf/nan is needed
    for (const auto& d : double_constants_) {
      for (double e : d) {
        if (isinf(e)) add_auxiliary(AUX_INF);
        if (isnan(e)) add_auxiliary(AUX_NAN);
      }
    }

    // Codegen auxiliary functions
    s << this->auxiliaries.str();

    // Print integer constants
    if (!integer_constants_.empty()) {
      for (casadi_int i=0; i<integer_constants_.size(); ++i) {
        print_vector(s, "casadi_s" + str(i), integer_constants_[i]);
      }
      s << std::endl;
    }

    // Print double constants
    if (!double_constants_.empty()) {
      for (casadi_int i=0; i<double_constants_.size(); ++i) {
        print_vector(s, "casadi_c" + str(i), double_constants_[i]);
      }
      s << std::endl;
    }

    // Print file scope double work
    if (!file_scope_double_.empty()) {
      casadi_int i=0;
      for (const auto& it : file_scope_double_) {
        s << "static casadi_real casadi_rd" + str(i++) + "[" + str(it.second) + "];\n";
      }
      s << std::endl;
    }

    // Print file scope integer work
    if (!file_scope_integer_.empty()) {
      casadi_int i=0;
      for (const auto& it : file_scope_integer_) {
        s << "static casadi_real casadi_ri" + str(i++) + "[" + str(it.second) + "];\n";
      }
      s << std::endl;
    }

    // External function declarations
    if (!added_externals_.empty()) {
      s << "/* External functions */\n";
      for (auto&& i : added_externals_) {
        s << i << std::endl;
      }
      s << std::endl << std::endl;
    }

    // Codegen body
    s << this->body.str();

    // End with new line
    s << std::endl;
  }

  std::string CodeGenerator::work(casadi_int n, casadi_int sz) const {
    if (n<0 || sz==0) {
      return "0";
    } else if (sz==1 && !this->codegen_scalars) {
      return "(&w" + str(n) + ")";
    } else {
      return "w" + str(n);
    }
  }

  std::string CodeGenerator::workel(casadi_int n) const {
    if (n<0) return "0";
    std::stringstream s;
    if (this->codegen_scalars) s << "*";
    s << "w" << n;
    return s.str();
  }

  std::string CodeGenerator::array(const std::string& type, const std::string& name, casadi_int len,
                                   const std::string& def) {
    std::stringstream s;
    s << type << " ";
    if (len==0) {
      s << "*" << name << " = 0";
    } else {
      s << name << "[" << len << "]";
      if (!def.empty()) s << " = " << def;
    }
    s << ";\n";
    return s.str();
  }

  void CodeGenerator::print_vector(std::ostream &s, const std::string& name,
      const std::vector<casadi_int>& v) {
    s << array("static const casadi_int", name, v.size(), initializer(v));
  }

  void CodeGenerator::print_vector(std::ostream &s, const std::string& name,
                                  const std::vector<double>& v) {
    s << array("static const casadi_real", name, v.size(), initializer(v));
  }

  std::string CodeGenerator::print_op(casadi_int op, const std::string& a0) {
    switch (op) {
      case OP_FABS:
        add_auxiliary(AUX_FABS);
        return "casadi_fabs("+a0+")";
      case OP_SQ:
        add_auxiliary(AUX_SQ);
        return "casadi_sq("+a0+")";
      case OP_SIGN:
        add_auxiliary(AUX_SIGN);
        return "casadi_sign("+a0+")";
      case OP_LOG1P:
        add_auxiliary(AUX_LOG1P);
        return "casadi_log1p("+a0+")";
      case OP_EXPM1:
        add_auxiliary(AUX_EXPM1);
        return "casadi_expm1("+a0+")";
      default:
        return casadi_math<double>::print(op, a0);
    }
  }
  std::string CodeGenerator::print_op(casadi_int op, const std::string& a0, const std::string& a1) {
    switch (op) {
      case OP_FMIN:
        add_auxiliary(AUX_FMIN);
        return "casadi_fmin("+a0+","+a1+")";
      case OP_FMAX:
        add_auxiliary(AUX_FMAX);
        return "casadi_fmax("+a0+","+a1+")";
      case OP_HYPOT:
        add_auxiliary(AUX_HYPOT);
        return "casadi_hypot("+a0+","+a1+")";
      default:
        return casadi_math<double>::print(op, a0, a1);
    }
  }

  void CodeGenerator::add_include(const std::string& new_include, bool relative_path,
                                 const std::string& use_ifdef) {
    // Register the new element
    bool added = added_includes_.insert(new_include).second;

    // Quick return if it already exists
    if (!added) return;

    // Ifdef opening
    if (!use_ifdef.empty()) this->includes << "#ifdef " << use_ifdef << std::endl;

    // Print to the header section
    if (relative_path) {
      this->includes << "#include \"" << new_include << "\"\n";
    } else {
      this->includes << "#include <" << new_include << ">\n";
    }

    // Ifdef closing
    if (!use_ifdef.empty()) this->includes << "#endif\n";
  }

  std::string CodeGenerator::
  operator()(const Function& f, const std::string& arg,
             const std::string& res, const std::string& iw,
             const std::string& w) {
    std::string name = add_dependency(f);
    bool needs_mem = !f->codegen_mem_type().empty();
    if (needs_mem) {
      std::string mem = "mid";
      local("flag", "int");
      local(mem, "int");
      *this << mem << " = " << name << "_checkout();\n";
      *this << "if (" << mem << "<0) return 1;\n";
      *this << "flag = " + name + "(" + arg + ", " + res + ", "
              + iw + ", " + w + ", " << mem << ");\n";
      *this << name << "_release(" << mem << ");\n";
      return "flag";
    } else {
      return name + "(" + arg + ", " + res + ", "
              + iw + ", " + w + ", 0)";
    }
  }

  void CodeGenerator::add_external(const std::string& new_external) {
    added_externals_.insert(new_external);
  }

  std::string CodeGenerator::shorthand(const std::string& name) const {
    casadi_assert(added_shorthands_.count(name), "No such macro: " + name);
    return "casadi_" + name;
  }

  std::string CodeGenerator::shorthand(const std::string& name, bool allow_adding) {
    bool added = added_shorthands_.insert(name).second;
    if (!allow_adding) {
      casadi_assert(added, "Duplicate macro: " + name);
    }
    return "casadi_" + name;
  }

  casadi_int CodeGenerator::add_sparsity(const Sparsity& sp) {
    return get_constant(sp, true);
  }

  std::string CodeGenerator::sparsity(const Sparsity& sp) {
    return shorthand("s" + str(add_sparsity(sp)));
  }

  casadi_int CodeGenerator::get_sparsity(const Sparsity& sp) const {
    return const_cast<CodeGenerator&>(*this).get_constant(sp, false);
  }

  size_t CodeGenerator::hash(const std::vector<double>& v) {
    // Calculate a hash value for the vector
    std::size_t seed=0;
    if (!v.empty()) {
      casadi_assert_dev(sizeof(double) % sizeof(size_t)==0);
      const casadi_int int_len = v.size()*(sizeof(double)/sizeof(size_t));
      const size_t* int_v = reinterpret_cast<const size_t*>(&v.front());
      for (size_t i=0; i<int_len; ++i) {
        hash_combine(seed, int_v[i]);
      }
    }
    return seed;
  }

  size_t CodeGenerator::hash(const std::vector<casadi_int>& v) {
    size_t seed=0;
    hash_combine(seed, v);
    return seed;
  }

  casadi_int CodeGenerator::get_constant(const std::vector<double>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    auto eq = added_double_constants_.equal_range(h);
    for (auto i=eq.first; i!=eq.second; ++i) {
      if (equal(v, double_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      casadi_int ind = double_constants_.size();
      double_constants_.push_back(v);
      added_double_constants_.insert(std::make_pair(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  casadi_int CodeGenerator::get_constant(const std::vector<casadi_int>& v, bool allow_adding) {
    // Hash the vector
    size_t h = hash(v);

    // Try to locate it in already added constants
    std::pair<std::multimap<size_t, size_t>::iterator, std::multimap<size_t, size_t>::iterator> eq =
      added_integer_constants_.equal_range(h);
    for (std::multimap<size_t, size_t>::iterator i=eq.first; i!=eq.second; ++i) {
      if (equal(v, integer_constants_[i->second])) return i->second;
    }

    if (allow_adding) {
      // Add to constants
      casadi_int ind = integer_constants_.size();
      integer_constants_.push_back(v);
      added_integer_constants_.insert(std::pair<size_t, size_t>(h, ind));
      return ind;
    } else {
      casadi_error("Constant not found");
      return -1;
    }
  }

  std::string CodeGenerator::constant(const std::vector<casadi_int>& v) {
    return shorthand("s" + str(get_constant(v, true)));
  }

  void CodeGenerator::constant_copy(
      const std::string& name, const std::vector<casadi_int>& v, const std::string& type) {
    std::string ref = constant(v);
    if (!v.empty()) {
      local(name+"[" + str(v.size()) + "]", type);
    } else {
      local(name, type, "*");
    }
    if (!v.empty()) {
      local("i", type);
      (*this) << "for (i=0;i<" << v.size() << ";++i) " + name + "[i] = " + ref + "[i];\n";
    } else {
      init_local(name, "0");
    }
  }

  std::string CodeGenerator::constant(const std::vector<double>& v) {
    return shorthand("c" + str(get_constant(v, true)));
  }

  void CodeGenerator::add_auxiliary(Auxiliary f, const std::vector<std::string>& inst) {
    // Look for existing instantiations
    auto f_match = added_auxiliaries_.equal_range(f);
    // Look for duplicates
    for (auto it=f_match.first; it!=f_match.second; ++it) {
      if (it->second==inst) return;
    }
    added_auxiliaries_.insert(std::make_pair(f, inst));

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
      add_auxiliary(AUX_INTERPN_WEIGHTS);
      add_auxiliary(AUX_INTERPN_INTERPOLATE);
      add_auxiliary(AUX_FLIP, {});
      add_auxiliary(AUX_CLEAR);
      add_auxiliary(AUX_CLEAR, {"casadi_int"});
      this->auxiliaries << sanitize_source(casadi_interpn_str, inst);
      break;
    case AUX_INTERPN_GRAD:
      add_auxiliary(AUX_INTERPN);
      this->auxiliaries << sanitize_source(casadi_interpn_grad_str, inst);
      break;
    case AUX_DE_BOOR:
      this->auxiliaries << sanitize_source(casadi_de_boor_str, inst);
      break;
    case AUX_ND_BOOR_EVAL:
      add_auxiliary(AUX_DE_BOOR);
      add_auxiliary(AUX_FILL);
      add_auxiliary(AUX_FILL, {"casadi_int"});
      add_auxiliary(AUX_CLEAR);
      add_auxiliary(AUX_CLEAR, {"casadi_int"});
      add_auxiliary(AUX_LOW);
      this->auxiliaries << sanitize_source(casadi_nd_boor_eval_str, inst);
      break;
    case AUX_FLIP:
      this->auxiliaries << sanitize_source(casadi_flip_str, inst);
      break;
    case AUX_LOW:
      this->auxiliaries << sanitize_source(casadi_low_str, inst);
      break;
    case AUX_INTERPN_WEIGHTS:
      add_auxiliary(AUX_LOW);
      this->auxiliaries << sanitize_source(casadi_interpn_weights_str, inst);
      break;
    case AUX_INTERPN_INTERPOLATE:
      this->auxiliaries << sanitize_source(casadi_interpn_interpolate_str, inst);
      break;
    case AUX_NORM_1:
      this->auxiliaries << sanitize_source(casadi_norm_1_str, inst);
      break;
    case AUX_NORM_2:
      add_auxiliary(AUX_DOT);
      this->auxiliaries << sanitize_source(casadi_norm_2_str, inst);
      break;
    case AUX_NORM_INF:
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_norm_inf_str, inst);
      break;
    case AUX_VECTOR_FMAX:
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_vector_fmax_str, inst);
      break;
    case AUX_VECTOR_FMIN:
      add_auxiliary(AUX_FMIN);
      this->auxiliaries << sanitize_source(casadi_vector_fmin_str, inst);
      break;
    case AUX_MASKED_NORM_INF:
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_masked_norm_inf_str, inst);
      break;
    case AUX_CLIP_MIN:
      // add_auxiliary(AUX_CLIP_MIN);
      this->auxiliaries << sanitize_source(casadi_clip_min_str, inst);
      break;
    case AUX_CLIP_MAX:
      // add_auxiliary(AUX_CLIP_MAX);
      this->auxiliaries << sanitize_source(casadi_clip_max_str, inst);
      break;
    case AUX_CLEAR:
      this->auxiliaries << sanitize_source(casadi_clear_str, inst);
      break;
    case AUX_FILL:
      this->auxiliaries << sanitize_source(casadi_fill_str, inst);
      break;
    case AUX_MV:
      this->auxiliaries << sanitize_source(casadi_mv_str, inst);
      break;
    case AUX_MV_DENSE:
      this->auxiliaries << sanitize_source(casadi_mv_dense_str, inst);
      break;
    case AUX_MTIMES:
      this->auxiliaries << sanitize_source(casadi_mtimes_str, inst);
      break;
    case AUX_TRILSOLVE:
      this->auxiliaries << sanitize_source(casadi_trilsolve_str, inst);
      break;
    case AUX_TRIUSOLVE:
      this->auxiliaries << sanitize_source(casadi_triusolve_str, inst);
      break;
    case AUX_PROJECT:
      this->auxiliaries << sanitize_source(casadi_project_str, inst);
      break;
    case AUX_TRI_PROJECT:
      this->auxiliaries << sanitize_source(casadi_tri_project_str, inst);
      break;
    case AUX_DENSIFY:
      add_auxiliary(AUX_CLEAR);
      add_auxiliary(AUX_CAST);
      {
        std::vector<std::string> inst2 = inst;
        if (inst.size()==1) inst2.push_back(inst[0]);
        this->auxiliaries << sanitize_source(casadi_densify_str, inst2);
      }
      break;
    case AUX_SPARSIFY:
      add_auxiliary(AUX_CAST);
      {
        std::vector<std::string> inst2 = inst;
        if (inst.size()==1) inst2.push_back(inst[0]);
        this->auxiliaries << sanitize_source(casadi_sparsify_str, inst2);
      }
      break;
    case AUX_TRANS:
      this->auxiliaries << sanitize_source(casadi_trans_str, inst);
      break;
    case AUX_TO_MEX:
      add_auxiliary(AUX_TO_DOUBLE);
      this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                        << sanitize_source(casadi_to_mex_str, inst)
                        << "#endif\n\n";
      break;
    case AUX_FROM_MEX:
      add_auxiliary(AUX_FILL);
      this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                        << sanitize_source(casadi_from_mex_str, inst)
                        << "#endif\n\n";
      break;
    case AUX_FINITE_DIFF:
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_finite_diff_str, inst);
      break;
    case AUX_QR:
      add_auxiliary(AUX_IF_ELSE);
      add_auxiliary(AUX_SCAL);
      add_auxiliary(AUX_DOT);
      add_auxiliary(AUX_CLEAR);
      this->auxiliaries << sanitize_source(casadi_qr_str, inst);
      break;
    case AUX_LSQR:
      add_auxiliary(AUX_COPY);
      add_auxiliary(AUX_CLEAR);
      add_auxiliary(AUX_NORM_2);
      add_auxiliary(AUX_MV);
      add_auxiliary(AUX_FABS);
      add_auxiliary(AUX_SIGN);
      this->auxiliaries << sanitize_source(casadi_lsqr_str, inst);
      break;
    case AUX_QP:
      this->auxiliaries << sanitize_source(casadi_qp_str, inst);
      break;
    case AUX_QRQP:
      add_auxiliary(AUX_QP);
      add_auxiliary(AUX_COPY);
      add_auxiliary(AUX_QR);
      add_auxiliary(AUX_MAX);
      add_auxiliary(AUX_FMIN);
      add_auxiliary(AUX_FMAX);
      add_auxiliary(AUX_TRANS);
      add_auxiliary(AUX_AXPY);
      add_auxiliary(AUX_MV);
      add_auxiliary(AUX_BILIN);
      add_auxiliary(AUX_INF);
      add_auxiliary(AUX_REAL_MIN);
      add_auxiliary(AUX_CLEAR);
      add_include("stdarg.h");
      add_include("stdio.h");
      add_include("math.h");

      this->auxiliaries << sanitize_source(casadi_qrqp_str, inst);
      break;
    case AUX_NLP:
      this->auxiliaries << sanitize_source(casadi_nlp_str, inst);
      break;
    case AUX_SQPMETHOD:
      add_auxiliary(AUX_QP);
      add_auxiliary(AUX_NLP);
      this->auxiliaries << sanitize_source(casadi_sqpmethod_str, inst);
      break;
    case AUX_FEASIBLESQPMETHOD:
      add_auxiliary(AUX_QP);
      add_auxiliary(AUX_NLP);
      this->auxiliaries << sanitize_source(casadi_feasiblesqpmethod_str, inst);
      break;
    case AUX_LDL:
      this->auxiliaries << sanitize_source(casadi_ldl_str, inst);
      break;
    case AUX_NEWTON:
      add_auxiliary(AUX_COPY);
      add_auxiliary(AUX_AXPY);
      add_auxiliary(AUX_NORM_INF);
      add_auxiliary(AUX_QR);
      this->auxiliaries << sanitize_source(casadi_newton_str, inst);
      break;
    case AUX_MAX_VIOL:
      this->auxiliaries << sanitize_source(casadi_max_viol_str, inst);
      break;
    case AUX_SUM_VIOL:
      this->auxiliaries << sanitize_source(casadi_sum_viol_str, inst);
      break;
    case AUX_SUM:
      this->auxiliaries << sanitize_source(casadi_sum_str, inst);
      break;
    case AUX_VFMIN:
      add_auxiliary(AUX_FMIN);
      this->auxiliaries << sanitize_source(casadi_vfmin_str, inst);
      break;
    case AUX_VFMAX:
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_vfmax_str, inst);
      break;
    case AUX_REGULARIZE:
      add_auxiliary(AUX_FABS);
      add_auxiliary(AUX_FMIN);
      this->auxiliaries << sanitize_source(casadi_regularize_str, inst);
      break;
    case AUX_BOUNDS_CONSISTENCY:
      add_auxiliary(AUX_ISINF);
      add_auxiliary(AUX_FMAX);
      add_auxiliary(AUX_FMIN);
      this->auxiliaries << sanitize_source(casadi_bound_consistency_str, inst);
      break;
    case AUX_FILE_SLURP:
      add_include("stdio.h");
      this->auxiliaries << sanitize_source(casadi_file_slurp_str, inst);
      break;
    case AUX_CACHE:
      this->auxiliaries << sanitize_source(casadi_cache_str, inst);
      break;
    case AUX_CVX:
      add_auxiliary(AUX_CLEAR);
      add_auxiliary(AUX_FABS);
      add_auxiliary(AUX_COPY);
      add_auxiliary(AUX_DOT);
      add_auxiliary(AUX_AXPY);
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_cvx_str, inst);
      break;
    case AUX_CONVEXIFY:
      add_auxiliary(AUX_CVX);
      add_auxiliary(AUX_PROJECT);
      add_auxiliary(AUX_REGULARIZE);
      add_auxiliary(AUX_COPY);
      this->auxiliaries << sanitize_source(casadi_convexify_str, inst);
      break;
    case AUX_LOGSUMEXP:
      add_auxiliary(AUX_LOG1P);
      add_auxiliary(AUX_FMAX);
      this->auxiliaries << sanitize_source(casadi_logsumexp_str, inst);
      break;
    case AUX_SPARSITY:
      this->auxiliaries << sanitize_source(casadi_sparsity_str, inst);
      break;
    case AUX_TO_DOUBLE:
      this->auxiliaries << "#define casadi_to_double(x) "
                        << "(" << (this->cpp ? "static_cast<double>(x)" : "(double) x") << ")\n\n";
      break;
    case AUX_TO_INT:
      this->auxiliaries << "#define casadi_to_int(x) "
                        << "(" << (this->cpp ? "static_cast<casadi_int>(x)" : "(casadi_int) x")
                        << ")\n\n";
      break;
    case AUX_CAST:
      this->auxiliaries << "#define CASADI_CAST(x,y) "
                        << "(" << (this->cpp ? "static_cast<x>(y)" : "(x) y") << ")\n\n";
      break;
    case AUX_SQ:
      shorthand("sq");
      this->auxiliaries << "casadi_real casadi_sq(casadi_real x) { return x*x;}\n\n";
      break;
    case AUX_SIGN:
      shorthand("sign");
      this->auxiliaries << "casadi_real casadi_sign(casadi_real x) "
                        << "{ return x<0 ? -1 : x>0 ? 1 : x;}\n\n";
      break;
    case AUX_IF_ELSE:
      shorthand("if_else");
      this->auxiliaries << "casadi_real casadi_if_else"
                        << "(casadi_real c, casadi_real x, casadi_real y) "
                        << "{ return c!=0 ? x : y;}\n\n";
      break;
    case AUX_PRINTF:
      this->auxiliaries << "#ifndef CASADI_PRINTF\n";
      if (this->mex) {
        this->auxiliaries << "#ifdef MATLAB_MEX_FILE\n"
                          << "  #define CASADI_PRINTF mexPrintf\n"
                          << "#else\n"
                          << "  #define CASADI_PRINTF printf\n"
                          << "#endif\n";
      } else {
        add_include("stdio.h");
        this->auxiliaries << "#define CASADI_PRINTF printf\n";
        this->auxiliaries << "#ifndef CASADI_SNPRINTF\n";
        this->auxiliaries << "#define CASADI_SNPRINTF snprintf\n";
        this->auxiliaries << "#endif\n\n";
      }
      this->auxiliaries << "#endif\n\n";
      break;
    case AUX_FMIN:
      shorthand("fmin");
      this->auxiliaries << "casadi_real casadi_fmin(casadi_real x, casadi_real y) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return x<y ? x : y;\n"
                        << "#else\n"
                        << "  return fmin(x, y);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    case AUX_FMAX:
      shorthand("fmax");
      this->auxiliaries << "casadi_real casadi_fmax(casadi_real x, casadi_real y) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return x>y ? x : y;\n"
                        << "#else\n"
                        << "  return fmax(x, y);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    case AUX_FABS:
      shorthand("fabs");
      this->auxiliaries << "casadi_real casadi_fabs(casadi_real x) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return x>0 ? x : -x;\n"
                        << "#else\n"
                        << "  return fabs(x);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    case AUX_ISINF:
      shorthand("isinf");
      this->auxiliaries << "casadi_real casadi_isinf(casadi_real x) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return x== INFINITY || x==-INFINITY;\n"
                        << "#else\n"
                        << "  return isinf(x);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    case AUX_MIN:
      this->auxiliaries << "casadi_int casadi_min(casadi_int x, casadi_int y) {\n"
                        << "  return x>y ? y : x;\n"
                        << "}\n\n";
      break;
    case AUX_MAX:
      this->auxiliaries << "casadi_int casadi_max(casadi_int x, casadi_int y) {\n"
                        << "  return x>y ? x : y;\n"
                        << "}\n\n";
      break;
    case AUX_MMIN:
      add_auxiliary(AUX_VFMIN);
      add_auxiliary(AUX_INF);
      this->auxiliaries << sanitize_source(casadi_mmin_str, inst);
      break;
    case AUX_MMAX:
      add_auxiliary(AUX_VFMAX);
      add_auxiliary(AUX_INF);
      this->auxiliaries << sanitize_source(casadi_mmax_str, inst);
      break;
    case AUX_INF:
      this->auxiliaries << "#ifndef casadi_inf\n"
                        << "  #define casadi_inf " << this->infinity << "\n"
                        << "#endif\n\n";
      break;
    case AUX_NAN:
      this->auxiliaries << "#ifndef casadi_nan\n"
                        << "  #define casadi_nan " << this->nan << "\n"
                        << "#endif\n\n";
      break;
    case AUX_REAL_MIN:
      this->auxiliaries << "#ifndef casadi_real_min\n"
                        << "  #define casadi_real_min " << this->real_min << "\n"
                        << "#endif\n\n";
      break;
    case AUX_LOG1P:
      shorthand("log1p");
      this->auxiliaries << "casadi_real casadi_log1p(casadi_real x) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return log(1+x);\n"
                        << "#else\n"
                        << "  return log1p(x);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    case AUX_EXPM1:
      shorthand("expm1");
      this->auxiliaries << "casadi_real casadi_expm1(casadi_real x) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return exp(x)-1;\n"
                        << "#else\n"
                        << "  return expm1(x);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    case AUX_HYPOT:
      shorthand("hypot");
      this->auxiliaries << "casadi_real casadi_hypot(casadi_real x, casadi_real y) {\n"
                        << "/* Pre-c99 compatibility */\n"
                        << "#if __STDC_VERSION__ < 199901L\n"
                        << "  return sqrt(x*x+y*y);\n"
                        << "#else\n"
                        << "  return hypot(x, y);\n"
                        << "#endif\n"
                        << "}\n\n";
      break;
    }
  }

  std::string CodeGenerator::to_mex(const Sparsity& sp, const std::string& arg) {
    add_auxiliary(AUX_TO_MEX);
    std::stringstream s;
    s << "casadi_to_mex(" << sparsity(sp) << ", " << arg << ");";
    return s.str();
  }

  std::string CodeGenerator::from_mex(std::string& arg,
                                      const std::string& res, std::size_t res_off,
                                      const Sparsity& sp_res, const std::string& w) {
    // Handle offset with recursion
    if (res_off!=0) return from_mex(arg, res+"+"+str(res_off), 0, sp_res, w);

    add_auxiliary(AUX_FROM_MEX);
    std::stringstream s;
    s << "casadi_from_mex(" << arg
      << ", " << res << ", " << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  std::string CodeGenerator::fmu_helpers(const std::string& modelname) {
    // Process C++ source
    std::stringstream ret;
    std::string line;
    std::istringstream stream(casadi_fmu_str);
    while (std::getline(stream, line)) {
      // Replacements
      if (line.find("MODELNAME") != std::string::npos) {
        line = replace(line, "MODELNAME", modelname);
      }
      // Append to return
      ret << line << "\n";
    }
    return ret.str();
  }

  std::string CodeGenerator::constant(const std::string& v) {
    std::string ret = v;
    ret = replace(ret, "\\", "\\\\");
    ret = replace(ret, "\"", "\\\"");
    return "\"" + ret + "\"";
  }

  std::string CodeGenerator::constant(casadi_int v) {
    return constant(static_cast<double>(v));
  }
  std::string CodeGenerator::constant(double v) {
    std::stringstream s;
    if (isnan(v)) {
      add_auxiliary(AUX_NAN);
      s << "casadi_nan";
    } else if (isinf(v)) {
      add_auxiliary(AUX_INF);
      if (v<0) s << "-";
      s << "casadi_inf";
    } else {
      casadi_int v_int = static_cast<casadi_int>(v);
      if (static_cast<double>(v_int)==v) {
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

  std::string CodeGenerator::initializer(const std::vector<double>& v) {
    std::stringstream s;
    s << "{";
    for (casadi_int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << constant(v[i]);
    }
    s << "}";
    return s.str();
  }

  std::string CodeGenerator::initializer(const std::vector<casadi_int>& v) {
    std::stringstream s;
    s << "{";
    for (casadi_int i=0; i<v.size(); ++i) {
      if (i!=0) s << ", ";
      s << v[i];
    }
    s << "}";
    return s.str();
  }

  std::string CodeGenerator::copy(const std::string& arg,
                                  std::size_t n, const std::string& res) {
    std::stringstream s;
    // Perform operation
    add_auxiliary(AUX_COPY);
    s << "casadi_copy(" << arg << ", " << n << ", " << res << ");";
    return s.str();
  }

  void CodeGenerator::copy_check(const std::string& arg, size_t n, const std::string& res,
      bool check_lhs, bool check_rhs) {
    std::vector<std::string> checks;
    if (check_lhs) checks.push_back(arg);
    if (check_rhs) checks.push_back(res);
    if (!checks.empty()) *this << "if (" << join(checks, " && ") << ") ";
    *this << copy(arg, n, res) << "\n";
  }

  void CodeGenerator::copy_default(const std::string& arg, size_t n, const std::string& res,
      const std::string& def, bool check_rhs) {
    *this << "if (" << arg << ") {\n";
    if (check_rhs) *this << "if (" << res << ") ";
    *this << copy(arg, n, res) << "\n";
    *this << "} else {\n";
    if (check_rhs) *this << "if (" << res << ") ";
    *this << fill(res, n, def) << "\n";
    *this << "}\n";
  }

  std::string CodeGenerator::clear(const std::string& res, std::size_t n) {
    std::stringstream s;
    // Perform operation
    add_auxiliary(AUX_CLEAR);
    s << "casadi_clear(" << res << ", " << n << ");";
    return s.str();
  }

  std::string CodeGenerator::arg(casadi_int i) const {
    return "arg[" + str(i) + "]";
  }

  std::string CodeGenerator::res(casadi_int i) const {
    return "res[" + str(i) + "]";
  }

  std::string CodeGenerator::mem(const Function& f) {
    std::string name = f->codegen_name(*this, false);
    std::string mem_array = shorthand(name + "_mem");
    return mem_array+"[mem]";
  }

  std::string CodeGenerator::fill(const std::string& res,
                                  std::size_t n, const std::string& v) {
    if (v=="0") return clear(res, n);
    std::stringstream s;
    // Perform operation
    add_auxiliary(AUX_FILL);
    s << "casadi_fill(" << res << ", " << n << ", " << v << ");";
    return s.str();
  }

  std::string CodeGenerator::dot(casadi_int n, const std::string& x,
                                 const std::string& y) {
    add_auxiliary(AUX_DOT);
    std::stringstream s;
    s << "casadi_dot(" << n << ", " << x << ", " << y << ")";
    return s.str();
  }

  std::string CodeGenerator::bilin(const std::string& A, const Sparsity& sp_A,
                                   const std::string& x, const std::string& y) {
    add_auxiliary(AUX_BILIN);
    std::stringstream s;
    s << "casadi_bilin(" << A << ", " << sparsity(sp_A) << ", " << x << ", " << y << ")";
    return s.str();
  }

  std::string CodeGenerator::rank1(const std::string& A, const Sparsity& sp_A,
                                   const std::string& alpha, const std::string& x,
                                   const std::string& y) {
    add_auxiliary(AUX_RANK1);
    std::stringstream s;
    s << "casadi_rank1(" << A << ", " << sparsity(sp_A) << ", "
      << alpha << ", " << x << ", " << y << ");";
    return s.str();
  }

  std::string CodeGenerator::interpn(
      const std::string& res, casadi_int ndim, const std::string& grid,
      const std::string& offset,
      const std::string& values, const std::string& x,
      const std::string& lookup_mode, casadi_int m,
      const std::string& iw, const std::string& w) {
    add_auxiliary(AUX_INTERPN);
    std::stringstream s;
    s << "casadi_interpn(" << res << ", " << ndim << ", " << grid << ", "  << offset << ", "
      << values << ", " << x << ", " << lookup_mode << ", " << m << ", " << iw << ", " << w << ");";
    return s.str();
  }

  std::string CodeGenerator::interpn_grad(const std::string& grad,
      casadi_int ndim, const std::string& grid, const std::string& offset,
      const std::string& values, const std::string& x,
      const std::string& lookup_mode, casadi_int m,
      const std::string& iw, const std::string& w) {
    add_auxiliary(AUX_INTERPN_GRAD);
    std::stringstream s;
    s << "casadi_interpn_grad(" << grad << ", " << ndim << ", " << grid << ", " << offset << ", "
      << values << ", " << x << ", " << lookup_mode << "," << m << ", " << iw << ", " << w << ");";
    return s.str();
  }

  std::string CodeGenerator::trans(const std::string& x, const Sparsity& sp_x,
                                   const std::string& y, const Sparsity& sp_y,
                                   const std::string& iw) {
    add_auxiliary(CodeGenerator::AUX_TRANS);
    return "casadi_trans(" + x + "," + sparsity(sp_x) + ", "
            + y + ", " + sparsity(sp_y) + ", " + iw + ")";
  }

  std::string CodeGenerator::declare(std::string s) {
    // Add c linkage
    std::string cpp_prefix = this->cpp ? "extern \"C\" " : "";

    // To header file
    if (this->with_header) {
      this->header << cpp_prefix << this->dll_import << s << ";\n";
    }

    // Return name with declarations
    return cpp_prefix + this->dll_export + s;
  }

  std::string
  CodeGenerator::project(const std::string& arg, const Sparsity& sp_arg,
                         const std::string& res, const Sparsity& sp_res,
                         const std::string& w) {
    // If sparsity match, simple copy
    if (sp_arg==sp_res) return copy(arg, sp_arg.nnz(), res);

    // Create call
    add_auxiliary(AUX_PROJECT);
    std::stringstream s;
    s << "casadi_project(" << arg << ", " << sparsity(sp_arg) << ", " << res << ", "
      << sparsity(sp_res) << ", " << w << ");";
    return s.str();
  }

  std::string
  CodeGenerator::tri_project(const std::string& arg, const Sparsity& sp_arg,
                         const std::string& res, bool lower) {
    // Create call
    add_auxiliary(AUX_TRI_PROJECT);
    std::stringstream s;
    s << "casadi_tri_project(" << arg << ", " << sparsity(sp_arg) << ", ";
    s << res  << ", " << (lower ? 1: 0) << ");";
    return s.str();
  }

  std::string
  CodeGenerator::densify(const std::string& arg, const Sparsity& sp_arg,
                        const std::string& res, bool tr) {
    // Create call
    add_auxiliary(AUX_DENSIFY);
    std::stringstream s;
    s << "casadi_densify(" << arg << ", " << sparsity(sp_arg) << ", " << res << ", "
      << (tr ? 1 : 0) << ");";
    return s.str();
  }

  std::string
  CodeGenerator::sparsify(const std::string& arg, const std::string& res,
                          const Sparsity& sp_res, bool tr) {
    // Create call
    add_auxiliary(AUX_SPARSIFY);
    std::stringstream s;
    s << "casadi_sparsify(" << arg << ", " << res << ", "
      << sparsity(sp_res) << ", " << (tr ? 1 : 0) << ");";
    return s.str();
  }

  std::string CodeGenerator::printf(const std::string& str, const std::vector<std::string>& arg) {
    add_auxiliary(AUX_PRINTF);
    std::stringstream s;
    s << "CASADI_PRINTF(\"" << str << "\"";
    for (casadi_int i=0; i<arg.size(); ++i) s << ", " << arg[i];
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

  std::string CodeGenerator::axpy(casadi_int n, const std::string& a,
                                  const std::string& x, const std::string& y) {
    add_auxiliary(AUX_AXPY);
    return "casadi_axpy(" + str(n) + ", " + a + ", " + x + ", " + y + ");";
  }

  std::string CodeGenerator::clip_min(const std::string& x, casadi_int n,
                                  const std::string& min, const std::string& mask) {
    add_auxiliary(AUX_CLIP_MIN);
    return "casadi_clip_min(" + x + ", " + str(n) + ", " + min + ", " + mask + ");";
  }

  std::string CodeGenerator::clip_max(const std::string& x, casadi_int n,
                                  const std::string& min, const std::string& mask) {
    add_auxiliary(AUX_CLIP_MAX);
    return "casadi_clip_max(" + x + ", " + str(n) + ", " + min + ", " + mask + ");";
  }

  std::string CodeGenerator::vector_fmax(casadi_int n, const std::string& x,
                                  const std::string& y, const std::string& z) {
    add_auxiliary(AUX_VECTOR_FMAX);
    return "casadi_vector_fmax(" + str(n) + ", " + x + ", " + y + ", " + z + ");";
  }

  std::string CodeGenerator::vector_fmin(casadi_int n, const std::string& x,
                                  const std::string& y, const std::string& z) {
    add_auxiliary(AUX_VECTOR_FMIN);
    return "casadi_vector_fmin(" + str(n) + ", " + x + ", " + y + ", " + z + ");";
  }

  std::string CodeGenerator::masked_norm_inf(casadi_int n, const std::string& x,
                                  const std::string& mask) {
    add_auxiliary(AUX_MASKED_NORM_INF);
    return "casadi_masked_norm_inf(" + str(n) + ", " + x + ", " + mask + ")";
  }

  std::string CodeGenerator::scal(casadi_int n, const std::string& alpha, const std::string& x) {
    add_auxiliary(AUX_SCAL);
    return "casadi_scal(" + str(n) + ", " + alpha + ", " + x + ");";
  }

  std::string CodeGenerator::mv(const std::string& x, const Sparsity& sp_x,
                                const std::string& y, const std::string& z, bool tr) {
    add_auxiliary(AUX_MV);
    return "casadi_mv(" + x + ", " + sparsity(sp_x) + ", " + y + ", "
           + z + ", " +  (tr ? "1" : "0") + ");";
  }

  std::string CodeGenerator::mv(const std::string& x, casadi_int nrow_x, casadi_int ncol_x,
                                const std::string& y, const std::string& z, bool tr) {
    add_auxiliary(AUX_MV_DENSE);
    return "casadi_mv_dense(" + x + ", " + str(nrow_x) + ", " + str(ncol_x) + ", "
           + y + ", " + z + ", " +  (tr ? "1" : "0") + ");";
  }

  std::string CodeGenerator::mtimes(const std::string& x, const Sparsity& sp_x,
                                    const std::string& y, const Sparsity& sp_y,
                                    const std::string& z, const Sparsity& sp_z,
                                    const std::string& w, bool tr) {
    add_auxiliary(AUX_MTIMES);
    return "casadi_mtimes(" + x + ", " + sparsity(sp_x) + ", " + y + ", " + sparsity(sp_y) + ", "
      + z + ", " + sparsity(sp_z) + ", " + w + ", " +  (tr ? "1" : "0") + ");";
  }

  std::string CodeGenerator::trilsolve(const Sparsity& sp_x, const std::string& x,
      const std::string& y, bool tr, bool unity, casadi_int nrhs) {
    add_auxiliary(AUX_TRILSOLVE);
    return "casadi_trilsolve(" + sparsity(sp_x) + ", " + x + ", " + y + ", " + str(tr) + ", "
        + str(unity) + ", " + str(nrhs) + ");";
  }

  std::string CodeGenerator::triusolve(const Sparsity& sp_x, const std::string& x,
      const std::string& y, bool tr, bool unity, casadi_int nrhs) {
    add_auxiliary(AUX_TRIUSOLVE);
    return "casadi_triusolve(" + sparsity(sp_x) + ", " + x + ", " + y + ", " + str(tr) + ", "
        + str(unity) + ", " + str(nrhs) + ");";
  }


  std::string CodeGenerator::logsumexp(const std::string& A, casadi_int n) {
    add_auxiliary(AUX_LOGSUMEXP);
    std::stringstream s;
    s << "casadi_logsumexp(" << A << ", " << n << ");";
    return s.str();
  }

  void CodeGenerator::print_formatted(const std::string& s) {
    // Quick return if empty
    if (s.empty()) return;

    // If new line, add indentation
    if (newline_) {
      casadi_int shift = s.front()=='}' ? -1 : 0;
      casadi_assert_dev(current_indent_+shift>=0);
      this->buffer << std::string(indent_*(current_indent_+shift), ' ');
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

  CodeGenerator& CodeGenerator::operator<<(const std::string& s) {
    // Loop over newline characters
    size_t off=0;
    while (true) {
      size_t pos = s.find('\n', off);
      if (pos==std::string::npos) {
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
    this->buffer.str(std::string());
  }

  void CodeGenerator::local(const std::string& name, const std::string& type,
                            const std::string& ref) {
    // Check if the variable already exists
    auto it = local_variables_.find(name);
    if (it==local_variables_.end()) {
      // Add it
      local_variables_[name] = std::make_pair(type, ref);
    } else {
      // Consistency check
      casadi_assert(it->second.first==type, "Type mismatch for " + name);
      casadi_assert(it->second.second==ref, "Type mismatch for " + name);
    }
  }

  std::string CodeGenerator::sx_work(casadi_int i) {
    if (avoid_stack_) {
      return "w[" + str(i) + "]";
    } else {
      std::string name = "a"+str(i);

      // Make sure work vector element has been declared
      local(name, "casadi_real");

      return name;
    }
  }

  void CodeGenerator::init_local(const std::string& name, const std::string& def) {
    bool inserted = local_default_.insert(std::make_pair(name, def)).second;
    casadi_assert(inserted, name + " already defined");
  }

  std::string CodeGenerator::
  sanitize_source(const std::string& src,
                  const std::vector<std::string>& inst, bool add_shorthand) {
    // Create suffix if templates type are not all "casadi_real"
    std::string suffix;
    for (const std::string& s : inst) {
      if (s!="casadi_real") {
        for (const std::string& s : inst) suffix += "_" + s;
        break;
      }
    }

    // Construct map of name replacements
    std::vector<std::pair<std::string, std::string> > rep;
    for (casadi_int i=0; i<inst.size(); ++i) {
      rep.push_back(std::make_pair("T" + str(i+1), inst[i]));
    }

    // Return object
    std::stringstream ret;
    // Process C++ source
    std::string line;
    std::istringstream stream(src);
    while (std::getline(stream, line)) {
      size_t n1, n2;

      // C++ template declarations are ignored
      if (line.find("template")==0) continue;

      // Macro definitions are ignored
      if (line.find("#define")==0) continue;
      if (line.find("#undef")==0) continue;

      // Inline declaration
      if (line == "inline") continue;

      // If line starts with "// SYMBOL", add shorthand
      if (line.find("// SYMBOL") != std::string::npos) {
        n1 = line.find("\"");
        n2 = line.find("\"", n1+1);
        std::string sym = line.substr(n1+1, n2-n1-1);
        if (add_shorthand) shorthand(sym + suffix);
        if (!suffix.empty()) {
          rep.push_back(std::make_pair(sym, sym + suffix));
        }
        continue;
      }

      // If line starts with "// C-REPLACE", add to list of replacements
      if (line.find("// C-REPLACE") != std::string::npos) {
        // Get C++ string
        n1 = line.find("\"");
        n2 = line.find("\"", n1+1);
        std::string key = line.substr(n1+1, n2-n1-1);
        // Get C string
        n1 = line.find("\"", n2+1);
        n2 = line.find("\"", n1+1);
        std::string sub = line.substr(n1+1, n2-n1-1);
        // Add to replacements
        rep.push_back(std::make_pair(key, sub));
        continue;
      }

      // If line starts with "// C-VERBOSE", skip the next line
      if (!verbose_runtime && line.find("// C-VERBOSE") != std::string::npos) {
        // Ignore next line
        std::getline(stream, line);
        continue;
      }

      // Ignore other C++ style comment
      if ((n1 = line.find("//")) != std::string::npos) line.erase(n1);

      // Remove trailing spaces
      if ((n1 = line.find_last_not_of(' ')) != std::string::npos) {
        line.erase(n1 + 1);
      } else {
        continue;
      }

      // Perform string replacements
      for (auto&& it = rep.rbegin(); it!=rep.rend(); ++it) {
        line = replace(line, it->first, it->second);
      }

      // Append to return
      ret << line << "\n";
    }

    // Trailing newline
    ret << "\n";
    return ret.str();
  }

  void CodeGenerator::comment(const std::string& s) {
    if (verbose) {
      *this << "/* " << s << " */\n";
    }
  }

  void CodeGenerator::
  add_io_sparsities(const std::string& name,
                    const std::vector<Sparsity>& sp_in,
                    const std::vector<Sparsity>& sp_out) {
    // Insert element, quick return if it already exists
    if (!sparsity_meta.insert(name).second) return;

    // Input sparsities
    *this << declare("const casadi_int* " + name + "_sparsity_in(casadi_int i)") << " {\n"
      << "switch (i) {\n";
    for (casadi_int i=0; i<sp_in.size(); ++i) {
      *this << "case " << i << ": return " << sparsity(sp_in[i]) << ";\n";
    }
    *this << "default: return 0;\n}\n"
      << "}\n\n";

    // Output sparsities
    *this << declare("const casadi_int* " + name + "_sparsity_out(casadi_int i)") << " {\n"
      << "switch (i) {\n";
    for (casadi_int i=0; i<sp_out.size(); ++i) {
      *this << "case " << i << ": return " << sparsity(sp_out[i]) << ";\n";
    }
    *this << "default: return 0;\n}\n"
      << "}\n\n";
  }

  std::string CodeGenerator::
  qr(const std::string& sp, const std::string& A, const std::string& w,
      const std::string& sp_v, const std::string& v, const std::string& sp_r,
      const std::string& r, const std::string& beta, const std::string& prinv,
      const std::string& pc) {
    add_auxiliary(CodeGenerator::AUX_QR);
    return "casadi_qr(" + sp + ", " + A + ", " + w + ", "
           + sp_v + ", " + v + ", " + sp_r + ", " + r + ", "
           + beta + ", " + prinv + ", " + pc + ");";
  }

  std::string CodeGenerator::
  qr_solve(const std::string& x, casadi_int nrhs, bool tr,
      const std::string& sp_v, const std::string& v,
      const std::string& sp_r, const std::string& r,
      const std::string& beta, const std::string& prinv,
      const std::string& pc, const std::string& w) {
    add_auxiliary(CodeGenerator::AUX_QR);
    return "casadi_qr_solve(" + x + ", " + str(nrhs) + ", " + (tr ? "1" : "0") + ", "
           + sp_v + ", " + v + ", " + sp_r + ", " + r + ", "
           + beta + ", " + prinv + ", " + pc + ", " + w + ");";
  }

  std::string CodeGenerator::
  lsqr_solve(const std::string& A, const std::string&x,
             casadi_int nrhs, bool tr, const std::string& sp, const std::string& w) {
    add_auxiliary(CodeGenerator::AUX_LSQR);
    return "casadi_lsqr_solve(" + A + ", " + x + ", " + str(nrhs) + ", "
           + (tr ? "1" : "0") + ", " + sp + ", " + w + ");";
  }

  std::string CodeGenerator::
  ldl(const std::string& sp_a, const std::string& a,
      const std::string& sp_lt, const std::string& lt, const std::string& d,
      const std::string& p, const std::string& w) {
    add_auxiliary(CodeGenerator::AUX_LDL);
    return "casadi_ldl(" + sp_a + ", " + a + ", " + sp_lt + ", " + lt + ", "
           + d + ", " + p + ", " + w + ");";
  }

  std::string CodeGenerator::
  ldl_solve(const std::string& x, casadi_int nrhs,
    const std::string& sp_lt, const std::string& lt, const std::string& d,
    const std::string& p, const std::string& w) {
    add_auxiliary(CodeGenerator::AUX_LDL);
    return "casadi_ldl_solve(" + x + ", " + str(nrhs) + ", " + sp_lt + ", "
           + lt + ", " + d + ", " + p + ", " + w + ");";
  }

  std::string CodeGenerator::
  fmax(const std::string& x, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_FMAX);
    return "casadi_fmax(" + x + ", " + y + ");";
  }

  std::string CodeGenerator::
  fmin(const std::string& x, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_FMIN);
    return "casadi_fmin(" + x + ", " + y + ");";
  }

  std::string CodeGenerator::
  vfmax(const std::string& x, casadi_int n, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_VFMAX);
    return "casadi_vfmax(" + x + ", " + str(n) + ", " + y + ");";
  }

  std::string CodeGenerator::
  vfmin(const std::string& x, casadi_int n, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_VFMIN);
    return "casadi_vfmin(" + x + ", " + str(n) + ", " + y + ");";
  }

  std::string CodeGenerator::
  vfmax(const std::string& x, const std::string& n, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_VFMAX);
    return "casadi_vfmax(" + x + ", " + n + ", " + y + ");";
  }

  std::string CodeGenerator::
  vfmin(const std::string& x, const std::string& n, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_VFMIN);
    return "casadi_vfmin(" + x + ", " + n + ", " + y + ");";
  }

  std::string CodeGenerator::
  max(const std::string& x, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_MAX);
    return "casadi_max(" + x + ", " + y + ")";
  }

  std::string CodeGenerator::
  min(const std::string& x, const std::string& y) {
    add_auxiliary(CodeGenerator::AUX_MIN);
    return "casadi_min(" + x + ", " + y + ")";
  }

  std::string CodeGenerator::
  mmax(const std::string& x, casadi_int n, bool is_dense) {
    add_auxiliary(CodeGenerator::AUX_MMAX);
    return "casadi_mmax(" + x + ", " + str(n) + ", " + str(casadi_int(is_dense)) + ")";
  }

  std::string CodeGenerator::
  mmin(const std::string& x, casadi_int n, bool is_dense) {
    add_auxiliary(CodeGenerator::AUX_MMIN);
    return "casadi_mmin(" + x + ", " + str(n) + ", " + str(casadi_int(is_dense)) + ")";
  }

  std::string CodeGenerator::
  max_viol(casadi_int n, const std::string& x, const std::string& lb, const std::string& ub) {
    add_auxiliary(CodeGenerator::AUX_MAX_VIOL);
    return "casadi_max_viol(" + str(n) + ", " + x+ ", " + lb + ", " + ub + ")";
  }

  std::string CodeGenerator::
  sum_viol(casadi_int n, const std::string& x, const std::string& lb, const std::string& ub) {
    add_auxiliary(CodeGenerator::AUX_SUM_VIOL);
    return "casadi_sum_viol(" + str(n) + ", " + x+ ", " + lb + ", " + ub + ")";
  }

  std::string CodeGenerator::
  norm_inf(casadi_int n, const std::string& x) {
    add_auxiliary(CodeGenerator::AUX_NORM_INF);
    return "casadi_norm_inf(" + str(n) + ", " + x + ")";
  }

  std::string CodeGenerator::
  norm_2(casadi_int n, const std::string& x) {
    add_auxiliary(CodeGenerator::AUX_NORM_2);
    return "casadi_norm_2(" + str(n) + ", " + x + ")";
  }

  std::string CodeGenerator::
  lb_eig(const Sparsity& sp_h, const std::string& h) {
    add_auxiliary(CodeGenerator::AUX_REGULARIZE);
    return "casadi_lb_eig(" + sparsity(sp_h) + ", " + h + ")";
  }

  std::string CodeGenerator::
  regularize(const Sparsity& sp_h, const std::string& h, const std::string& reg) {
    add_auxiliary(CodeGenerator::AUX_REGULARIZE);
    return "casadi_regularize(" + sparsity(sp_h) + ", " + h + ", " + reg + ");";
  }

  std::string CodeGenerator::
  convexify_eval(const ConvexifyData &d,
    const std::string& Hin, const std::string& Hout, const std::string& iw, const std::string& w) {
      add_auxiliary(CodeGenerator::AUX_CONVEXIFY);
      return Convexify::generate(*this, d, Hin, Hout, iw, w);
  }

  std::string CodeGenerator::
  low(const std::string& x, const std::string& grid, casadi_int ng, casadi_int lookup_mode) {
    add_auxiliary(CodeGenerator::AUX_LOW);
    return "casadi_low(" + x + ", " + grid + ", " + str(ng) + ", " + str(lookup_mode) + ");";
  }

  std::string CodeGenerator::
  bound_consistency(casadi_int n, const std::string& x,
    const std::string& lam, const std::string& lbx, const std::string& ubx) {
      add_auxiliary(CodeGenerator::AUX_BOUNDS_CONSISTENCY);
      return "casadi_bound_consistency(" + str(n) + ", " + x + ", " + lam +
             ", " + lbx + ", " + ubx + ")";
    }

  std::string CodeGenerator::
  file_slurp(const std::string& fname, casadi_int n, const std::string& a) {
    add_auxiliary(CodeGenerator::AUX_FILE_SLURP);
    return "casadi_file_slurp(\"" + fname + "\", " + str(n) + ", " + a + ")";
  }

  std::string CodeGenerator::
  cache_check(const std::string& key, const std::string& cache, const std::string& loc,
        casadi_int stride, casadi_int sz, casadi_int key_sz, const std::string& val) {
    add_auxiliary(CodeGenerator::AUX_CACHE);
    return "cache_check(" + key + ", " + cache + ", " + loc + ", " +
    str(stride) + ", " + str(sz) + ", " + str(key_sz) + ", " + val + ")";
  }

  void CodeGenerator::sz_work(size_t& sz_arg, size_t& sz_res, size_t& sz_iw, size_t& sz_w) const {
    sz_arg = sz_res = sz_iw = sz_w = 0;
    for (auto&& f : added_functions_) {
      sz_arg = std::max(sz_arg, f.f.sz_arg());
      sz_res = std::max(sz_res, f.f.sz_res());
      sz_iw = std::max(sz_iw, f.f.sz_iw());
      sz_w = std::max(sz_w, f.f.sz_w());
    }
  }

} // namespace casadi
