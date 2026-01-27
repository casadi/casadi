
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


#include "mumps_interface.hpp"
#include "casadi/core/global_options.hpp"
#include <cstdio>

namespace casadi {

  static void dump_matrix_file(DMUMPS_STRUC_C* mumps_data, const char* filename, bool symmetric) {
    FILE* fh = fopen(filename, "w");
    if (symmetric) {
      fprintf(fh, "%%%%MatrixMarket matrix coordinate real symmetric\n");
    } else {
      fprintf(fh, "%%%%MatrixMarket matrix coordinate real general\n");
    }
    fprintf(fh, "%d %d %lld\n", mumps_data->n, mumps_data->n, (long long)mumps_data->nnz);
    for (long long i = 0; i < mumps_data->nnz; i++) {
      fprintf(fh, "%d %d %25.18e\n", mumps_data->irn[i], mumps_data->jcn[i], mumps_data->a[i]);
    }
    fclose(fh);
  }

  static void dump_rhs_file(const double* rhs, int n, int nrhs, const char* filename) {
    FILE* fh = fopen(filename, "w");
    fprintf(fh, "%%%%MatrixMarket matrix array real general\n");
    fprintf(fh, "%d %d\n", n, nrhs);
    // MatrixMarket array format is column-major
    for (int j = 0; j < nrhs; j++) {
      for (int i = 0; i < n; i++) {
        fprintf(fh, "%25.18e\n", rhs[j * n + i]);
      }
    }
    fclose(fh);
  }

  extern "C"
  int CASADI_LINSOL_MUMPS_EXPORT
  casadi_register_linsol_mumps(LinsolInternal::Plugin* plugin) {
    plugin->creator = MumpsInterface::creator;
    plugin->name = "mumps";
    plugin->doc = MumpsInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &MumpsInterface::options_;
    plugin->deserialize = &MumpsInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_MUMPS_EXPORT casadi_load_linsol_mumps() {
    LinsolInternal::registerPlugin(casadi_register_linsol_mumps);
  }

  MumpsInterface::MumpsInterface(const std::string& name, const Sparsity& sp)
    : LinsolInternal(name, sp) {
  }

  MumpsInterface::~MumpsInterface() {
    clear_mem();
  }

  const Options MumpsInterface::options_
  = {{&ProtoFunction::options_},
     {{"symmetric",
      {OT_BOOL,
       "Symmetric matrix"}},
      {"posdef",
       {OT_BOOL,
       "Positive definite"}},
      {"dump_mtx",
       {OT_BOOL,
       "Dump matrices to MatrixMarket files for debugging"}},
      {"dump_stats",
       {OT_BOOL,
       "Dump MUMPS statistics to log files"}},
      {"error_analysis",
       {OT_BOOL,
       "Enable MUMPS error analysis (ICNTL(11)=1)"}},
      {"print_level",
       {OT_INT,
       "MUMPS print level (0-4, default 0)"}}
     }
  };

  void MumpsInterface::init(const Dict& opts) {
    // Call the init method of the base class
    LinsolInternal::init(opts);

    // Default options
    symmetric_ = false;
    posdef_ = false;
    dump_mtx_ = false;
    dump_stats_ = false;
    error_analysis_ = false;
    print_level_ = 0;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="symmetric") {
        symmetric_ = op.second;
      } else if (op.first=="posdef") {
        posdef_ = op.second;
      } else if (op.first=="dump_mtx") {
        dump_mtx_ = op.second;
      } else if (op.first=="dump_stats") {
        dump_stats_ = op.second;
      } else if (op.first=="error_analysis") {
        error_analysis_ = op.second;
      } else if (op.first=="print_level") {
        print_level_ = op.second;
      }
    }

    // Options consistency
    if (posdef_ && !symmetric_) casadi_error("Inconsistent options");
  }

  int MumpsInterface::init_mem(void* mem) const {
    if (LinsolInternal::init_mem(mem)) return 1;
    auto m = static_cast<MumpsMemory*>(mem);
    // Free existing MUMPS instance, if any
    if (m->id) {
      // Terminate the instance of the package
      m->id->job = -2;
      dmumps_c(m->id);
      // Delete the memory structure
      delete m->id;
    }
    m->id = new DMUMPS_STRUC_C();

    // Initialize MUMPS
    m->id->job = -1;  // initializes an instance of the package
    m->id->par = 1;
    m->id->sym = symmetric_ ? posdef_ ? 2 : 1 : 0;
    m->id->comm_fortran = -987654;
    dmumps_c(m->id);

    // Sparsity pattern in MUMPS format
    casadi_int n = this->nrow();
    casadi_int nnz = symmetric_ ? sp_.nnz_upper() : this->nnz();
    m->nz.resize(nnz);
    m->irn.clear();
    m->jcn.clear();
    m->irn.reserve(nnz);
    m->jcn.reserve(nnz);
    const casadi_int* colind = this->colind();
    const casadi_int* row = this->row();
    for (casadi_int c = 0; c < n; ++c) {
      for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
        if (!symmetric_ || row[k] <= c) {
          m->irn.push_back(row[k] + 1);
          m->jcn.push_back(c + 1);
        }
      }
    }

    return 0;
  }

  int MumpsInterface::nfact(void* mem, const double* A) const {
    auto m = static_cast<MumpsMemory*>(mem);
    casadi_assert_dev(A!=nullptr);

    // Copy nonzero entries to m->nz
    auto nz_it = m->nz.begin();
    if (symmetric_) {
      // Get upper triangular entries
      casadi_int n = this->nrow();
      const casadi_int* colind = this->colind();
      const casadi_int* row = this->row();
      for (casadi_int c = 0; c < n; ++c) {
        for (casadi_int k = colind[c]; k < colind[c + 1]; ++k) {
          if (row[k] <= c) *nz_it++ = A[k];
        }
      }
    } else {
      // Copy all entries
      std::copy(A, A + this->nnz(), nz_it);
    }

    // Define problem
    m->id->n = this->nrow();
    m->id->nnz = m->nz.size();
    m->id->irn = get_ptr(m->irn);
    m->id->jcn = get_ptr(m->jcn);
    m->id->a = get_ptr(m->nz);

    // Output control
    if (print_level_ > 0) {
      m->id->icntl[1 - 1] = 6;  // error messages
      m->id->icntl[2 - 1] = 6;  // diagnostic messages
      m->id->icntl[3 - 1] = 6;  // global info
      m->id->icntl[4 - 1] = print_level_;  // print level
    } else {
      m->id->icntl[1 - 1] = -1;
      m->id->icntl[2 - 1] = -1;
      m->id->icntl[3 - 1] = -1;
      m->id->icntl[4 - 1] = 0;
    }

    // Error analysis
    if (error_analysis_) {
      m->id->icntl[11 - 1] = 1;
    }

    // Dump matrix before factorization
    if (dump_mtx_) {
      char buffer[64];
      sprintf(buffer, "mumps_mtx_it%06d.mtx", m->fact_counter);
      dump_matrix_file(m->id, buffer, symmetric_);
    }
    m->solve_counter = 0;
    m->fact_counter++;

    // Symbolic and numeric factorization
    m->id->job = 4;
    dmumps_c(m->id);

    // Dump statistics after factorization
    if (dump_stats_) {
      int it = m->fact_counter - 1;
      int sol = m->solve_counter;  // not incremented yet for factorization stats
      bool write_header = false;
      FILE* fh = fopen("mumps_stats.csv", "r");
      if (fh == NULL) {
        write_header = true;
      } else {
        fclose(fh);
      }
      fh = fopen("mumps_stats.csv", "a");
      if (write_header) {
        fprintf(fh, "iter,solve");
        for (int i = 0; i < 40; ++i) fprintf(fh, ",INFOG(%d)", i+1);
        for (int i = 0; i < 20; ++i) fprintf(fh, ",RINFOG(%d)", i+1);
        for (int i = 0; i < 40; ++i) fprintf(fh, ",ICNTL(%d)", i+1);
        for (int i = 0; i < 15; ++i) fprintf(fh, ",CNTL(%d)", i+1);
        fprintf(fh, "\n");
      }
      fprintf(fh, "%d,%d", it, sol);
      for (int i = 0; i < 40; ++i) fprintf(fh, ",%d", m->id->infog[i]);
      for (int i = 0; i < 20; ++i) fprintf(fh, ",%e", m->id->rinfog[i]);
      for (int i = 0; i < 40; ++i) fprintf(fh, ",%d", m->id->icntl[i]);
      for (int i = 0; i < 15; ++i) fprintf(fh, ",%e", m->id->cntl[i]);
      fprintf(fh, "\n");
      fclose(fh);
    }

    return 0;
  }

  int MumpsInterface::solve(void* mem, const double* A, double* x, casadi_int nrhs, bool tr) const {
    auto m = static_cast<MumpsMemory*>(mem);

    // Transpose or not?
    m->id->icntl[9 - 1] = tr ? 0 : 1;

    // Solve factorized linear system
    m->id->job = 3;

    // Dump all RHS before solve
    if (dump_mtx_) {
      char buffer[64];
      sprintf(buffer, "mumps_rhs_it%06d_%06d.mtx", m->fact_counter - 1, m->solve_counter++);
      dump_rhs_file(x, m->id->n, nrhs, buffer);
    }

    for (casadi_int i=0;i<nrhs;++i) {
      m->id->rhs = x;
      dmumps_c(m->id);
      x += m->id->n;
    }

    return 0;
  }

  MumpsMemory::MumpsMemory() {
    this->id = 0;
    this->fact_counter = 0;
    this->solve_counter = 0;
  }

  MumpsMemory::~MumpsMemory() {
    if (this->id) {
      // Terminate the instance of the package
      this->id->job = -2;
      dmumps_c(this->id);
      // Delete the structure
      delete this->id;
    }
  }

  MumpsInterface::MumpsInterface(DeserializingStream& s) : LinsolInternal(s) {
    int v = s.version("Mumps", 1, 2);
    s.unpack("MumpsInterface::symmetric", symmetric_);
    s.unpack("MumpsInterface::posdef", posdef_);
    if (v >= 2) {
      s.unpack("MumpsInterface::dump_mtx", dump_mtx_);
      s.unpack("MumpsInterface::dump_stats", dump_stats_);
      s.unpack("MumpsInterface::error_analysis", error_analysis_);
      s.unpack("MumpsInterface::print_level", print_level_);
    } else {
      dump_mtx_ = false;
      dump_stats_ = false;
      error_analysis_ = false;
      print_level_ = 0;
    }
  }

  void MumpsInterface::serialize_body(SerializingStream &s) const {
    LinsolInternal::serialize_body(s);
    s.version("Mumps", 2);
    s.pack("MumpsInterface::symmetric", symmetric_);
    s.pack("MumpsInterface::posdef", posdef_);
    s.pack("MumpsInterface::dump_mtx", dump_mtx_);
    s.pack("MumpsInterface::dump_stats", dump_stats_);
    s.pack("MumpsInterface::error_analysis", error_analysis_);
    s.pack("MumpsInterface::print_level", print_level_);
  }

} // namespace casadi
