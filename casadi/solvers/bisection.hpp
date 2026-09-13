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

#ifndef CASADI_BISECTION_HPP
#define CASADI_BISECTION_HPP

#include "casadi/core/rootfinder_impl.hpp"
#include <casadi/solvers/casadi_rootfinder_bisection_export.h>

/// \cond
namespace casadi {
  struct CASADI_ROOTFINDER_BISECTION_EXPORT BisectionMemory : public RootfinderMemory {
    int return_status;
    casadi_int iter;
    casadi_int search_iter;
    double f_mid;
    double bracket_width;
  };

  class CASADI_ROOTFINDER_BISECTION_EXPORT Bisection : public Rootfinder {
  public:
    explicit Bisection(const std::string &name, const Function &f);

    ~Bisection() override;

    const char *plugin_name() const override { return "bisection"; }

    std::string class_name() const override { return "Bisection"; }

    static Rootfinder *creator(const std::string &name, const Function &f) {
      return new Bisection(name, f);
    }

    static const Options options_;

    const Options &get_options() const override { return options_; }

    Dict get_stats(void *mem) const override;

    void init(const Dict &opts) override;

    void *alloc_mem() const override { return new BisectionMemory(); }

    int init_mem(void *mem) const override;

    void free_mem(void *mem) const override { delete static_cast<BisectionMemory *>(mem); }

    void set_work(void *mem, const double **&arg, double **&res,
                  casadi_int *&iw, double *&w) const override;

    int solve(void *mem) const override;

    static const std::string meta_doc;

    // void codegen_body(CodeGenerator &g) const override;

    // void codegen_declarations(CodeGenerator &g) const override;

    void serialize_body(SerializingStream &s) const override;

    static ProtoFunction *deserialize(DeserializingStream &s) { return new Bisection(s); }

  protected:
    explicit Bisection(DeserializingStream &s);

    casadi_int max_iter_;
    casadi_int max_search_;

    double abstol_;
    double abstol_step_;

    double search_step_;

    double lb_;
    double ub_;

    int finish(BisectionMemory *m, double x_sol, double f_sol, double width,
               int status, bool success, UnifiedReturnStatus urs) const;
    static std::string status_str(int status);
  };

} // namespace casadi

/// \endcond
#endif // CASADI_BISECTION_HPP
