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


#include "dump.hpp"
#include "casadi_os.hpp"
#include "dm.hpp"
#include "filesystem_impl.hpp"
#include <iomanip>
namespace casadi {

  void Dump::reset_dump_count() {
    dump_count_ = 0;
  }

  Dump::Dump(const MX& x, const std::string& base_filename,
             const std::string& dir, const std::string& format, bool verbose)
      : base_filename_(base_filename), dir_(dir), format_(format), verbose_(verbose) {
    casadi_assert_dev(x.nnz()>0);
    set_dep(x);
    set_sparsity(x.sparsity());
    finalize();
  }

  void ensure_directory_exists(const std::string& dir) {
    if (Filesystem::is_enabled()) {
      std::string effective_dir = Filesystem::ensure_trailing_slash(dir);
      casadi_assert(Filesystem::ensure_directory_exists(effective_dir),
        "Unable to create the required directory for '" + effective_dir + "'.");
    }
  }

  void Dump::finalize() {
    ensure_directory_exists(dir_);
  }

  std::string Dump::disp(const std::vector<std::string>& arg) const {
    return "dump(" + arg.at(0) + ", " + base_filename_ + ")";
  }

  void Dump::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
      const std::vector<bool>& unique) const {
    Dict opts;
    if (!dir_.empty()) opts["dir"] = dir_;
    if (!format_.empty()) opts["format"] = format_;
    if (verbose_) opts["verbose"] = true;
    res[0] = arg[0].dump(base_filename_, opts);
  }

  int Dump::eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const {
    if (arg[0]!=res[0]) {
      std::copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Dump::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    // Build filename with counter
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(6) << dump_count_++;
    std::string format = format_.empty() ? "mtx" : format_;
    std::string filename = base_filename_ + "." + ss.str() + "." + format;
    if (!dir_.empty()) filename = dir_ + filesep() + filename;
    // Dump to file
    if (verbose_) {
      uout() << "dump -> " << filename << std::endl;
    }
    DM::to_file(filename, sparsity(), arg[0], format);
    // Perform operation
    if (arg[0]!=res[0]) {
      std::copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Dump::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    if (arg[0]!=res[0]) {
      std::copy(arg[0], arg[0]+nnz(), res[0]);
    }
    return 0;
  }

  int Dump::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const {
    bvec_t *a = arg[0];
    bvec_t *r = res[0];
    casadi_int n = nnz();
    if (a != r) {
      for (casadi_int i=0; i<n; ++i) {
        *a++ |= *r;
        *r++ = 0;
      }
    }
    return 0;
  }

  void Dump::add_dependency(CodeGenerator& g, const Instance& inst,
      const Function& owner) const {
    ensure_directory_exists(g.dump_dir_prefix + dir_ + g.dump_dir_suffix);
  }

  void Dump::generate(CodeGenerator& g,
                          const std::vector<casadi_int>& arg,
                          const std::vector<casadi_int>& res,
                          const std::vector<bool>& arg_is_ref,
                          std::vector<bool>& res_is_ref,
                          bool prefer_inline) const {
    std::string format = format_.empty() ? "mtx" : format_;
    std::string effective_dir = g.dump_dir_prefix + dir_ + g.dump_dir_suffix;
    std::string prefix;
    if (!effective_dir.empty()) prefix = effective_dir + "/";
    // prefix + base_filename + "." + 6 digits + "." + format + null
    casadi_int buf_size = prefix.size() + base_filename_.size()
      + 1 + 6 + 1 + format.size() + 1;
    std::string a = g.work(arg[0], dep(0).nnz(), arg_is_ref[0]);
    // Block scope for dump variables
    g << "{\n";
    g << "static int dump_id = 0;\n";
    g << "char dump_fname[" << buf_size << "];\n";
    g << "FILE* dump_file;\n";
    g << "snprintf(dump_fname, " << buf_size << ", \"" << prefix << base_filename_
      << ".%06d." << format << "\", dump_id++);\n";
    if (verbose_) {
      g << g.printf("dump -> %s\\n", "dump_fname") << "\n";
    }
    g << "dump_file = fopen(dump_fname, \"w\");\n";
    g << "if (dump_file) {\n";
    g << g.to_file("dump_file", dep(0).sparsity(), a) << ";\n";
    g << "fclose(dump_file);\n";
    g << "}\n";
    g << "}\n";
    // Copy
    generate_copy(g, arg, res, arg_is_ref, res_is_ref, 0);
  }

  void Dump::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.version("Dump", 1);
    s.pack("Dump::base_filename", base_filename_);
    s.pack("Dump::dir", dir_);
    s.pack("Dump::format", format_);
    s.pack("Dump::verbose", verbose_);
  }

  Dump::Dump(DeserializingStream& s) : MXNode(s) {
    s.version("Dump", 1);
    s.unpack("Dump::base_filename", base_filename_);
    s.unpack("Dump::dir", dir_);
    s.unpack("Dump::format", format_);
    s.unpack("Dump::verbose", verbose_);
    finalize();
  }

} // namespace casadi
