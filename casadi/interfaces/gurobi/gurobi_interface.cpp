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


#include "gurobi_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/nlp_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_GUROBI_EXPORT
  casadi_register_conic_gurobi(Conic::Plugin* plugin) {
    plugin->creator = GurobiInterface::creator;
    plugin->name = "gurobi";
    plugin->doc = GurobiInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &GurobiInterface::options_;
    plugin->deserialize = &GurobiInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_GUROBI_EXPORT casadi_load_conic_gurobi() {
    Conic::registerPlugin(casadi_register_conic_gurobi);
  }

  GurobiInterface::GurobiInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
  }

  GurobiInterface::~GurobiInterface() {
    clear_mem();
  }

  const Options GurobiInterface::options_
  = {{&Conic::options_},
     {{"vtype",
       {OT_STRINGVECTOR,
        "Type of variables: [CONTINUOUS|binary|integer|semicont|semiint]"}},
      {"gurobi",
       {OT_DICT,
        "Options to be passed to gurobi."}},
      {"sos_groups",
       {OT_INTVECTORVECTOR,
        "Definition of SOS groups by indices."}},
      {"sos_weights",
       {OT_DOUBLEVECTORVECTOR,
        "Weights corresponding to SOS entries."}},
      {"sos_types",
       {OT_INTVECTOR,
        "Specify 1 or 2 for each SOS group."}}
     }
  };

  void GurobiInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    std::vector<std::string> vtype;

    std::vector< std::vector<casadi_int> > sos_groups;
    std::vector< std::vector<double> > sos_weights;
    std::vector<casadi_int> sos_types;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="vtype") {
        vtype = op.second;
      } else if (op.first=="gurobi") {
        opts_ = op.second;
      } else if (op.first=="sos_groups") {
        sos_groups = op.second.to_int_vector_vector();
      } else if (op.first=="sos_weights") {
        sos_weights = op.second.to_double_vector_vector();
      } else if (op.first=="sos_types") {
        sos_types = op.second.to_int_vector();
      }
    }

    // Validaty SOS constraints
    check_sos(nx_, sos_groups, sos_weights, sos_types);

    // Populate SOS structures
    flatten_nested_vector(sos_groups, sos_ind_, sos_beg_);
    if (!sos_weights.empty())
      flatten_nested_vector(sos_weights, sos_weights_);

    sos_types_ = to_int(sos_types);

    // Variable types
    if (!vtype.empty()) {
      casadi_assert(vtype.size()==nx_, "Option 'vtype' has wrong length");
      vtype_.resize(nx_);
      for (casadi_int i=0; i<nx_; ++i) {
        if (vtype[i]=="continuous") {
          vtype_[i] = GRB_CONTINUOUS;
        } else if (vtype[i]=="binary") {
          vtype_[i] = GRB_BINARY;
        } else if (vtype[i]=="integer") {
          vtype_[i] = GRB_INTEGER;
        } else if (vtype[i]=="semicont") {
          vtype_[i] = GRB_SEMICONT;
        } else if (vtype[i]=="semiint") {
          vtype_[i] = GRB_SEMIINT;
        } else {
          casadi_error("No such variable type: " + vtype[i]);
        }
      }
    }

    // Initialize SDP to SOCP memory
    sdp_to_socp_init(sdp_to_socp_mem_);

    // Temporary memory
    alloc_w(sdp_to_socp_mem_.indval_size, true); // val
    alloc_iw(sdp_to_socp_mem_.indval_size, true); // ind
    alloc_iw(nx_, true); // ind2
    alloc_iw(nx_, true); // vtypes
  }

  int GurobiInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<GurobiMemory*>(mem);

    // Load environment
    casadi_int flag = GRBloadenv(&m->env, nullptr); // no log file
    casadi_assert(!flag && m->env, "Failed to create GUROBI environment. Flag: " + str(flag)
      + ":" + GRBgeterrormsg(m->env));

    m->sos_weights = sos_weights_;
    m->sos_beg = sos_beg_;
    m->sos_ind = sos_ind_;
    m->sos_types = sos_types_;

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");
    return 0;
  }

  inline const char* return_status_string(casadi_int status) {
    switch (status) {
    case GRB_LOADED:
      return "LOADED";
    case GRB_OPTIMAL:
      return "OPTIMAL";
    case GRB_INFEASIBLE:
      return "INFEASIBLE";
    case GRB_INF_OR_UNBD:
      return "INF_OR_UNBD";
    case GRB_UNBOUNDED:
      return "UNBOUNDED";
    case GRB_CUTOFF:
      return "CUTOFF";
    case GRB_ITERATION_LIMIT:
      return "ITERATION_LIMIT";
    case GRB_NODE_LIMIT:
      return "NODE_LIMIT";
    case GRB_TIME_LIMIT:
      return "TIME_LIMIT";
    case GRB_SOLUTION_LIMIT:
      return "SOLUTION_LIMIT";
    case GRB_INTERRUPTED:
      return "INTERRUPTED";
    case GRB_NUMERIC:
      return "NUMERIC";
    case GRB_SUBOPTIMAL:
      return "SUBOPTIMAL";
    case GRB_INPROGRESS:
      return "INPROGRESS";
    }
    return "Unknown";
  }

  int GurobiInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<GurobiMemory*>(mem);
    const SDPToSOCPMem& sm = sdp_to_socp_mem_;

    // Statistics
    m->fstats.at("preprocessing").tic();

    // Problem has not been solved at this point
    m->return_status = -1;

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    // Inputs
    const double *h=arg[CONIC_H],
      *g=arg[CONIC_G],
      *a=arg[CONIC_A],
      *lba=arg[CONIC_LBA],
      *uba=arg[CONIC_UBA],
      *lbx=arg[CONIC_LBX],
      *ubx=arg[CONIC_UBX],
      *p=arg[CONIC_P],
      *q=arg[CONIC_Q],
      *x0=arg[CONIC_X0];
      //*lam_x0=arg[CONIC_LAM_X0];

    // Outputs
    double *x=res[CONIC_X],
      *cost=res[CONIC_COST],
      *lam_a=res[CONIC_LAM_A],
      *lam_x=res[CONIC_LAM_X];

    // Temporary memory
    double *val=w; w+=sm.indval_size;
    int *ind=reinterpret_cast<int*>(iw); iw+=sm.indval_size;
    int *ind2=reinterpret_cast<int*>(iw); iw+=nx_;
    char *vtypes=reinterpret_cast<char*>(iw); iw+=nx_;

    // Greate an empty model
    GRBmodel *model = nullptr;
    try {
      casadi_int flag = GRBnewmodel(m->env, &model, name_.c_str(), 0,
        nullptr, nullptr, nullptr, nullptr, nullptr);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      // Add variables
      for (casadi_int i=0; i<nx_; ++i) {
        // Get bounds
        double lb = lbx ? lbx[i] : 0., ub = ubx ? ubx[i] : 0.;
        if (isinf(lb)) lb = -GRB_INFINITY;
        if (isinf(ub)) ub =  GRB_INFINITY;

        // Get variable type
        char vtype;
        if (!vtype_.empty()) {
          // Explicitly set 'vtype' takes precedence
          vtype = vtype_.at(i);
        } else if (!discrete_.empty() && discrete_.at(i)) {
          // Variable marked as discrete (integer or binary)
          vtype = lb==0 && ub==1 ? GRB_BINARY : GRB_INTEGER;
        } else {
          // Continious variable
          vtype = GRB_CONTINUOUS;
        }
        vtypes[i] = vtype;

        // Pass to model
        flag = GRBaddvar(model, 0, nullptr, nullptr, g ? g[i] : 0., lb, ub, vtype, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      GRBupdatemodel(model);
      for (casadi_int i=0; i<nx_; ++i) {
        // If it is a discrete variable, we can pass the start guess
        if (vtypes[i] != GRB_CONTINUOUS) {
          flag = GRBsetdblattrelement(model, "Start", i, x0[i]);
          casadi_assert(!flag, GRBgeterrormsg(m->env));
        }
      }


      /*  Treat SOCP constraints */

      // Add helper variables for SOCP
      for (casadi_int i=0;i<sm.r.size()-1;++i) {
        for (casadi_int k=0;k<sm.r[i+1]-sm.r[i]-1;++k) {
          flag = GRBaddvar(model, 0, nullptr, nullptr, 0, -GRB_INFINITY, GRB_INFINITY,
                           GRB_CONTINUOUS, nullptr);
          casadi_assert(!flag, GRBgeterrormsg(m->env));
        }
        flag = GRBaddvar(model, 0, nullptr, nullptr, 0, 0, GRB_INFINITY, GRB_CONTINUOUS, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      flag = GRBupdatemodel(model);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      // Add quadratic terms
      const casadi_int *H_colind=H_.colind(), *H_row=H_.row();
      for (int i=0; i<nx_; ++i) {

        // Quadratic term nonzero indices
        casadi_int numqnz = H_colind[1]-H_colind[0];
        for (casadi_int k=0;k<numqnz;++k) ind[k]=H_row[k];
        H_colind++;
        H_row += numqnz;

        // Corresponding column
        casadi_fill(ind2, numqnz, i);

        // Quadratic term nonzeros
        if (h) {
          casadi_copy(h, numqnz, val);
          casadi_scal(numqnz, 0.5, val);
          h += numqnz;
        } else {
          casadi_clear(val, numqnz);
        }

        // Pass to model
        flag = GRBaddqpterms(model, numqnz, ind, ind2, val);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      // Add constraints
      const casadi_int *AT_colind=sm.AT.colind(), *AT_row=sm.AT.row();
      for (casadi_int i=0; i<na_; ++i) {
        // Get bounds
        double lb = lba ? lba[i] : 0., ub = uba ? uba[i] : 0.;

        casadi_int numnz = 0;
        // Loop over rows
        for (casadi_int k=AT_colind[i]; k<AT_colind[i+1]; ++k) {
          casadi_int j = AT_row[k];

          ind[numnz] = j;
          val[numnz] = a ? a[sm.A_mapping[k]]  : 0;

          numnz++;
        }

        // Pass to model
        if (isinf(lb)) {
          if (isinf(ub)) {
            // Neither upper or lower bounds, skip
          } else {
            // Only upper bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_LESS_EQUAL, ub, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          }
        } else {
          if (isinf(ub)) {
            // Only lower bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_GREATER_EQUAL, lb, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          } else if (lb==ub) {
            // Upper and lower bounds equal
            flag = GRBaddconstr(model, numnz, ind, val, GRB_EQUAL, lb, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          } else {
            // Both upper and lower bounds
            flag = GRBaddrangeconstr(model, numnz, ind, val, lb, ub, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
          }
        }
      }

      // Add SOS constraints when applicable
      if (!m->sos_ind.empty()) {
        flag = GRBaddsos(model, m->sos_beg.size()-1, m->sos_ind.size(),
            get_ptr(m->sos_types), get_ptr(m->sos_beg), get_ptr(m->sos_ind),
            get_ptr(m->sos_weights));
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      // SOCP helper constraints
      const Sparsity& sp = sm.map_Q.sparsity();
      const casadi_int* colind = sp.colind();
      const casadi_int* row = sp.row();
      const casadi_int* data = sm.map_Q.ptr();

      // Loop over columns
      for (casadi_int i=0; i<sp.size2(); ++i) {

        casadi_int numnz = 0;
        // Loop over rows
        for (casadi_int k=colind[i]; k<colind[i+1]; ++k) {
          casadi_int j = row[k];

          ind[numnz] = j;
          val[numnz] = (q && j<nx_) ? q[data[k]] : -1;

          numnz++;
        }

        // Get bound
        double bound = sm.map_P[i]==-1 ? 0 : -p[sm.map_P[i]];

        flag = GRBaddconstr(model, numnz, ind, val, GRB_EQUAL, bound, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      // Loop over blocks
      for (casadi_int i=0; i<sm.r.size()-1; ++i) {
        casadi_int block_size = sm.r[i+1]-sm.r[i];

        // Indicate x'x - y^2 <= 0
        for (casadi_int j=0;j<block_size;++j) {
          ind[j] = nx_ + sm.r[i] + j;
          val[j] = j<block_size-1 ? 1 : -1;
        }

        flag = GRBaddqconstr(model, 0, nullptr, nullptr,
          block_size, ind, ind, val,
          GRB_LESS_EQUAL, 0, nullptr);
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      flag = 0;
      for (auto && op : opts_) {
        int ret = GRBgetparamtype(m->env, op.first.c_str());
        switch (ret) {
          case -1:
            casadi_error("Parameter '" + op.first + "' unknown to Gurobi.");
          case 1:
            {
              flag = GRBsetintparam(GRBgetenv(model), op.first.c_str(), op.second);
              break;
            }
          case 2:
              flag = GRBsetdblparam(GRBgetenv(model), op.first.c_str(), op.second);
              break;
          case 3:
            {
              std::string s = op.second;
              flag = GRBsetstrparam(GRBgetenv(model), op.first.c_str(), s.c_str());
              break;
            }
          default:
            casadi_error("Not implememented : " + str(ret));
        }
        casadi_assert(!flag, GRBgeterrormsg(m->env));
      }

      m->fstats.at("preprocessing").toc();
      m->fstats.at("solver").tic();

      // Solve the optimization problem
      flag = GRBoptimize(model);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      m->fstats.at("solver").toc();
      m->fstats.at("postprocessing").tic();

      int optimstatus;
      flag = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      if (verbose_) uout() << "return status: " << return_status_string(optimstatus) <<
        " (" << optimstatus <<")" << std::endl;

      m->return_status = optimstatus;
      m->success = optimstatus==GRB_OPTIMAL;
      if (optimstatus==GRB_ITERATION_LIMIT || optimstatus==GRB_TIME_LIMIT
          || optimstatus==GRB_NODE_LIMIT || optimstatus==GRB_SOLUTION_LIMIT)
        m->unified_return_status = SOLVER_RET_LIMITED;

      // Get the objective value, if requested
      if (cost) {
        flag = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, cost);
        if (flag) cost[0] = casadi::nan;
      }

      // Get the optimal solution, if requested
      if (x) {
        flag = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nx_, x);
        if (flag) fill_n(x, nx_, casadi::nan);
      }

      if (lam_x) fill_n(lam_x, nx_, casadi::nan);
      if (lam_a) fill_n(lam_a, na_, casadi::nan);

      // Free memory
      GRBfreemodel(model);
      m->fstats.at("postprocessing").toc();

    } catch (...) {
      // Free memory
      if (model) GRBfreemodel(model);
      throw;
    }

    return 0;
  }

  Dict GurobiInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<GurobiMemory*>(mem);
    stats["return_status"] = return_status_string(m->return_status);
    return stats;
  }

  GurobiMemory::GurobiMemory() {
    this->env = nullptr;
  }

  GurobiMemory::~GurobiMemory() {
    if (this->env) GRBfreeenv(this->env);
  }

  GurobiInterface::GurobiInterface(DeserializingStream& s) : Conic(s) {
    s.version("GurobiInterface", 1);
    s.unpack("GurobiInterface::vtype", vtype_);
    s.unpack("GurobiInterface::opts", opts_);
    s.unpack("GurobiInterface::sos_weights", sos_weights_);
    s.unpack("GurobiInterface::sos_beg", sos_beg_);
    s.unpack("GurobiInterface::sos_ind", sos_ind_);
    s.unpack("GurobiInterface::sos_types", sos_types_);
    Conic::deserialize(s, sdp_to_socp_mem_);
  }

  void GurobiInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);
    s.version("GurobiInterface", 1);
    s.pack("GurobiInterface::vtype", vtype_);
    s.pack("GurobiInterface::opts", opts_);
    s.pack("GurobiInterface::sos_weights", sos_weights_);
    s.pack("GurobiInterface::sos_beg", sos_beg_);
    s.pack("GurobiInterface::sos_ind", sos_ind_);
    s.pack("GurobiInterface::sos_types", sos_types_);
    Conic::serialize(s, sdp_to_socp_mem_);
  }

} // namespace casadi
