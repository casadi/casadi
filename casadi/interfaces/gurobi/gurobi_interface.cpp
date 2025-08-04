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


#include "gurobi_interface.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/nlp_tools.hpp"

namespace casadi {

  // Helper functions for cleaner C API usage
  // Error handling macro
  #define GUROBI_CALL_WARN(func, msg) do { \
      int error = (func); \
      if (error) { \
          casadi_warning(msg ": Gurobi error " + std::to_string(error)); \
          return; \
      } \
  } while(0)

  // Convert string sense to Gurobi sense
  char sense_to_gurobi(const std::string& sense) {
    if (sense == "<=") return GRB_LESS_EQUAL;
    if (sense == ">=") return GRB_GREATER_EQUAL;
    if (sense == "=") return GRB_EQUAL;
    return GRB_LESS_EQUAL; // default
  }

  // RAII wrapper for callback data access
  class CallbackDataHelper {
  private:
    void* cbdata_;
    int where_;

  public:
    CallbackDataHelper(void* cbdata, int where) : cbdata_(cbdata), where_(where) {}

    bool getDouble(int what, double& value) {
      int error = GRBcbget(cbdata_, where_, what, &value);
      return error == 0;
    }

    bool getInt(int what, int& value) {
      int error = GRBcbget(cbdata_, where_, what, &value);
      return error == 0;
    }

    bool getSolution(std::vector<double>& solution) {
      int error = GRBcbget(cbdata_, where_, GRB_CB_MIPSOL_SOL, solution.data());
      return error == 0;
    }

    bool addLazyConstraint(const std::vector<int>& ind,
                                          const std::vector<double>& val,
                                          char sense, double rhs) const {
      if (where_ != GRB_CB_MIPSOL) {
        casadi_warning("Can only add lazy constraints in MIPSOL context");
        return 1;
      }

      int error = GRBcblazy(cbdata_, static_cast<int>(ind.size()), ind.data(), val.data(), sense, rhs);
      if (error != 0) {
        casadi_warning("GRBcblazy failed with code: " + std::to_string(error));
      }
      return error == 0;
    }
  };

  // Static C callback function
  static int gurobi_mipsol_callback(GRBmodel *model, void *cbdata, int where, void *usrdata) {
    try {
      if (where == GRB_CB_MIPSOL && usrdata) {
        GurobiInterface* interface = static_cast<GurobiInterface*>(usrdata);
        interface->handle_mipsol_callback(model, cbdata, where);
      }
    } catch (const std::exception& e) {
      // Don't let exceptions escape to C code
      casadi_warning("Exception in Gurobi callback: " + std::string(e.what()));
    }
    return 0; // Continue optimization
  }

  extern "C"
  int CASADI_CONIC_GUROBI_EXPORT
  casadi_register_conic_gurobi(Conic::Plugin* plugin) {
    plugin->creator = GurobiInterface::creator;
    plugin->name = "gurobi";
    plugin->doc = GurobiInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &GurobiInterface::options_;
    plugin->deserialize = &GurobiInterface::deserialize;
    #ifdef GUROBI_ADAPTOR
      char buffer[400];
      int ret = gurobi_adaptor_load(buffer, sizeof(buffer));
      if (ret!=0) {
        casadi_warning("Failed to load Gurobi adaptor: " + std::string(buffer) + ".");
        return 1;
      }
    #endif
    return 0;
  }

  extern "C"
  void CASADI_CONIC_GUROBI_EXPORT casadi_load_conic_gurobi() {
    Conic::registerPlugin(casadi_register_conic_gurobi);
  }

  GurobiInterface::GurobiInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st),
      enable_mipsol_callback_(false),
      mipsol_callback_(Function())
  {

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
        "Specify 1 or 2 for each SOS group."}},
      {"mipsol_callback",
        {OT_FUNCTION,
          "User callback function for MIPSOL events. "
          "Input: dict with solution data. Output: dict with lazy constraints."}},
      {"enable_mipsol_callback",
        {OT_BOOL,
          "Enable MIPSOL callbacks for dynamic constraint generation"}}
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
        else if (op.first == "mipsol_callback") {
            try {
                // Attempt to convert to function
                mipsol_callback_ = op.second;
                casadi_message("Successfully obtained callback function");

                // Verify signature
                if (mipsol_callback_.n_in() != 5 || mipsol_callback_.n_out() != 3) {
                    casadi_error("Callback function has wrong signature. Expected 5 inputs, 3 outputs");
                }
            } catch (const std::exception& e) {
                casadi_error("Failed to get callback function: " + std::string(e.what()));
            }
        }
        else if (op.first == "enable_mipsol_callback") {
            enable_mipsol_callback_ = op.second;
        }
      }
    // Final validation
    if (enable_mipsol_callback_) {
      if (mipsol_callback_.is_null()) {
          casadi_error("mipsol_callback must be provided when enable_mipsol_callback is true");
      }
      casadi_message("Callback setup complete");
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

    m->pool_sol_nr = 0;

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
          // Continuous variable
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

      std::vector<char> constraint_type(na_); // For each a entry: 0 absent, 1 linear
      casadi_int npi = 0;

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

        constraint_type[i] = 1;
        // Pass to model
        if (isinf(lb)) {
          if (isinf(ub)) {
            constraint_type[i] = 0;
            // Neither upper or lower bounds, skip
          } else {
            // Only upper bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_LESS_EQUAL, ub, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
            npi++;
          }
        } else {
          if (isinf(ub)) {
            // Only lower bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_GREATER_EQUAL, lb, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
            npi++;
          } else if (lb==ub) {
            // Upper and lower bounds equal
            flag = GRBaddconstr(model, numnz, ind, val, GRB_EQUAL, lb, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
            npi++;
          } else {
            // Both upper and lower bounds
            flag = GRBaddrangeconstr(model, numnz, ind, val, lb, ub, nullptr);
            casadi_assert(!flag, GRBgeterrormsg(m->env));
            npi++;
          }
        }
      }
      std::vector<double> pi(npi);

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

      // Enable lazy constraints and callback if needed
      if (enable_mipsol_callback_) {
        // Enable lazy constraints
        int error = GRBsetintparam(GRBgetenv(model), "LazyConstraints", 1);
        if (error) {
          casadi_warning("Failed to enable LazyConstraints parameter");
          return 1;
        }
        // Set the callback function
        error = GRBsetcallbackfunc(model, gurobi_mipsol_callback,
                                  const_cast<GurobiInterface*>(this));
        if (error) {
          casadi_warning("Failed to set callback function");
          return 1;
        }
      }

      m->fstats.at("preprocessing").toc();
      m->fstats.at("solver").tic();

      // Solve the optimization problem
      flag = GRBoptimize(model);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      m->fstats.at("solver").toc();
      m->fstats.at("postprocessing").tic();

      int optimstatus;
      flag = GRBgetintattr(model, "Status", &optimstatus);
      casadi_assert(!flag, GRBgeterrormsg(m->env));

      if (verbose_) uout() << "return status: " << return_status_string(optimstatus) <<
        " (" << optimstatus <<")" << std::endl;

      m->return_status = optimstatus;
      m->d_qp.success = optimstatus==GRB_OPTIMAL;
      if (optimstatus==GRB_ITERATION_LIMIT || optimstatus==GRB_TIME_LIMIT
          || optimstatus==GRB_NODE_LIMIT || optimstatus==GRB_SOLUTION_LIMIT)
        m->d_qp.unified_return_status = SOLVER_RET_LIMITED;

      // Get the objective value, if requested
      if (cost) {
        flag = GRBgetdblattr(model, "ObjVal", cost);
        if (flag) cost[0] = casadi::nan;
      }

      // Get the optimal solution, if requested
      if (x) {
        flag = GRBgetdblattrarray(model, "X", 0, nx_, x);
        if (flag) std::fill_n(x, nx_, casadi::nan);
      }
      if (lam_x) {
        flag = GRBgetdblattrarray(model, "RC", 0, nx_, lam_x);
        if (!flag) casadi_scal(nx_, -1.0, lam_x);
        if (flag) std::fill_n(lam_x, nx_, casadi::nan);
      }
      if (lam_a) {
        flag = GRBgetdblattrarray(model, "Pi", 0, npi, get_ptr(pi));
        if (flag) {
          std::fill_n(lam_a, na_, casadi::nan);
        } else {
          const double * p = get_ptr(pi);
          for (casadi_int i=0;i<na_;++i) {
            if (constraint_type[i]==0) {
              lam_a[i] = 0;
            } else if (constraint_type[i]==1) {
              lam_a[i] = -(*p++);
            }
          }
        }
      }

      // Get solutions from solution pool
      flag = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &(m->pool_sol_nr));
      if (!flag && m->pool_sol_nr > 0) {
        m->pool_obj_vals = std::vector<double>(m->pool_sol_nr, casadi::nan);
        for (int idx = 0; idx < m->pool_sol_nr; ++idx) {
          std::vector<double> x_pool(nx_);

          flag = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_SOLUTIONNUMBER, idx);
          if (!flag) flag = GRBgetdblattr(model, GRB_DBL_ATTR_POOLOBJVAL, &(m->pool_obj_vals[idx]));
          if (!flag) flag = GRBgetdblattrarray(model, GRB_DBL_ATTR_XN, 0, nx_, get_ptr(x_pool));

          if (flag) {
            m->pool_obj_vals[idx] = casadi::nan;
            std::fill(x_pool.begin(), x_pool.end(), casadi::nan);
          }

          m->pool_solutions.push_back(x_pool);
        }
      }

      if (enable_mipsol_callback_) {
      // Clear callback to avoid dangling pointer
      GRBsetcallbackfunc(model, nullptr, nullptr);
      }

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
    stats["pool_sol_nr"] = m->pool_sol_nr;
    stats["pool_obj_val"] = m->pool_obj_vals;
    stats["pool_solutions"] = m->pool_solutions;
    return stats;
  }

  GurobiMemory::GurobiMemory() {
    this->env = nullptr;
  }

  GurobiMemory::~GurobiMemory() {
    if (this->env) GRBfreeenv(this->env);
  }

  GurobiInterface::GurobiInterface(DeserializingStream& s) : Conic(s) {
    s.version("GurobiInterface", 2);
    s.unpack("GurobiInterface::enable_mipsol_callback", enable_mipsol_callback_);
    s.unpack("GurobiInterface::mipsol_callback", mipsol_callback_);
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
    s.version("GurobiInterface", 2);
    s.pack("GurobiInterface::enable_mipsol_callback", enable_mipsol_callback_);
    s.pack("GurobiInterface::mipsol_callback", mipsol_callback_);
    s.pack("GurobiInterface::vtype", vtype_);
    s.pack("GurobiInterface::opts", opts_);
    s.pack("GurobiInterface::sos_weights", sos_weights_);
    s.pack("GurobiInterface::sos_beg", sos_beg_);
    s.pack("GurobiInterface::sos_ind", sos_ind_);
    s.pack("GurobiInterface::sos_types", sos_types_);
    Conic::serialize(s, sdp_to_socp_mem_);
  }

void GurobiInterface::handle_mipsol_callback(GRBmodel *model, void *cbdata, int where) {
  try {
    if (!enable_mipsol_callback_ || mipsol_callback_.is_null()) {
      casadi_warning("Callback triggered but not properly initialized");
      return;
    }

    CallbackDataHelper helper(cbdata, where);

    // Get the number of variables
    int numvars;
    int error = GRBgetintattr(model, "NumVars", &numvars);
    if (error) {
      casadi_warning("Failed to get number of variables in callback");
      return;
    }

    // Get the integer solution
    std::vector<double> x_vals(numvars);
    if (!helper.getSolution(x_vals)) {
      casadi_warning("Failed to get solution in MIPSOL callback");
      return;
    }

    // Convert to CasADi DM (only use the variables we care about)
    DM x_solution = DM::zeros(nx_, 1);
    for (casadi_int i = 0; i < nx_ && i < numvars; ++i) {
      x_solution(i) = x_vals[i];
    }

    // Get additional solution information
    double obj_val, obj_best, obj_bound;
    int sol_count;

    if (!helper.getDouble(GRB_CB_MIPSOL_OBJ, obj_val) ||
        !helper.getDouble(GRB_CB_MIPSOL_OBJBST, obj_best) ||
        !helper.getDouble(GRB_CB_MIPSOL_OBJBND, obj_bound) ||
        !helper.getInt(GRB_CB_MIPSOL_SOLCNT, sol_count)) {
      casadi_warning("Failed to get callback information");
      return;
    }

    // Prepare input data for callback
    std::vector<double> input_data;

    // Add x_solution data
    const double* x_ptr = x_solution.ptr();
    input_data.insert(input_data.end(), x_ptr, x_ptr + nx_);

    // Add scalar values
    input_data.push_back(obj_val);
    input_data.push_back(obj_best);
    input_data.push_back(obj_bound);
    input_data.push_back(static_cast<double>(sol_count));

    // Prepare output buffers
    double flag = 0.0;                      // First output: boolean (as double)
    std::vector<double> a_vec(nx_, 0.0);    // Second output: vector A
    double b_val = 0.0;                     // Third output: scalar b

    // Get memory requirements for the callback
    casadi_int sz_arg = mipsol_callback_.sz_arg();
    casadi_int sz_res = mipsol_callback_.sz_res();
    casadi_int sz_iw = mipsol_callback_.sz_iw();
    casadi_int sz_w = mipsol_callback_.sz_w();

    // Allocate work arrays
    std::vector<casadi_int> iw(sz_iw);
    std::vector<double> w(sz_w);

    // Set up argument pointers
    std::vector<const double*> arg(sz_arg);
    std::vector<double*> res(sz_res);

    // Fill input arguments
    if (sz_arg > 0) {
      size_t offset = 0;
      if (sz_arg > 0) {
        arg[0] = &input_data[offset];  // x
        offset += nx_;
      }
      for (casadi_int i = 1; i < sz_arg && offset < input_data.size(); ++i) {
        arg[i] = &input_data[offset];  // Scalars
        offset += 1;
      }
    }

    // Set up output pointers
    if (sz_res >= 1) res[0] = &flag;
    if (sz_res >= 2) res[1] = a_vec.data();
    if (sz_res >= 3) res[2] = &b_val;

    // Call the callback function
    int callback_result = mipsol_callback_(arg.data(), res.data(),
                                           iw.data(), w.data(), 0);

    if (callback_result != 0) {
      casadi_warning("Callback function returned error code: " + std::to_string(callback_result));
      return;
    }

    // Process results if the callback requests constraint
    if (flag > 0.5) {
      DM A = DM::zeros(1, nx_);
      for (casadi_int i = 0; i < nx_; ++i) {
        A(0, i) = a_vec[i];
      }
      DM b = DM(b_val);

      std::vector<DM> callback_outputs = {flag, A, b};
      process_lazy_constraints(model, cbdata, callback_outputs);
    }

  } catch (const std::exception& e) {
    casadi_warning("Error in handle_mipsol_callback: " + std::string(e.what()));
  }
}

void GurobiInterface::process_lazy_constraints(GRBmodel *model, void *cbdata,
                                               const std::vector<DM>& callback_result) {
  try {
    if (callback_result.size() < 3) {
      casadi_warning("Callback result does not contain enough data");
      return;
    }

    // Expected order: [add_lazy_constraint_flag, A_lazy, b_lazy, (optional) sense]
    double add_lazy_flag = static_cast<double>(callback_result[0]);
    if (add_lazy_flag < 0.5) return;

    DM A_lazy = callback_result[1];
    DM b_lazy = callback_result[2];
    std::string sense = "<=";

    CallbackDataHelper helper(cbdata, GRB_CB_MIPSOL);

    for (casadi_int i = 0; i < A_lazy.size1(); ++i) {
      std::vector<int> cind;
      std::vector<double> cval;

      for (casadi_int j = 0; j < A_lazy.size2() && j < nx_; ++j) {
        double coeff = static_cast<double>(A_lazy(i, j));
        if (std::abs(coeff) > 1e-12) {
          cind.push_back(static_cast<int>(j));
          cval.push_back(coeff);
        }
      }

      if (!cind.empty()) {
        double rhs = static_cast<double>(b_lazy(i));
        char gurobi_sense = sense_to_gurobi(sense);

        if (!helper.addLazyConstraint(cind, cval, gurobi_sense, rhs)) {
          casadi_warning("Failed to add lazy constraint " + std::to_string(i));
        }
      }
    }

  } catch (const std::exception& e) {
    casadi_warning("Error processing lazy constraints: " + std::string(e.what()));
  }
}

} // namespace casadi
