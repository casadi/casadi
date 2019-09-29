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

#include "cbc_interface.hpp"

#include "CbcSOS.hpp"
#include "casadi/core/nlp_tools.hpp"

namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_CONIC_CBC_EXPORT
  casadi_register_conic_cbc(Conic::Plugin* plugin) {
    plugin->creator = CbcInterface::creator;
    plugin->name = "cbc";
    plugin->doc = CbcInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &CbcInterface::options_;
    plugin->deserialize = &CbcInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_CBC_EXPORT casadi_load_conic_cbc() {
    Conic::registerPlugin(casadi_register_conic_cbc);
  }


  CbcInterface::CbcInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {

    hot_start_ = false;
  }

  const Options CbcInterface::options_
  = {{&Conic::options_},
     {{"cbc",
       {OT_DICT,
        "Options to be passed to CBC."
        "Three sets of options are supported. "
        "The first can be found in OsiSolverParameters.hpp. "
        "The second can be found in CbcModel.hpp. "
        "The third are options that can be passed to CbcMain1."
        }},
      {"sos_groups",
       {OT_INTVECTORVECTOR,
        "Definition of SOS groups by indices."}},
      {"sos_weights",
       {OT_DOUBLEVECTORVECTOR,
        "Weights corresponding to SOS entries."}},
      {"sos_types",
       {OT_INTVECTOR,
        "Specify 1 or 2 for each SOS group."}},
      {"hot_start",
       {OT_BOOL,
        "Hot start with x0 [Default false]."}},
     }
   };

  inline std::string return_status_string(int status) {
    switch (status) {
    case -1:
      return "before branchAndBound";
    case 0:
      return "finished";
    case 1:
      return "stopped - on maxnodes, maxsols, maxtime";
    case 2:
      return "difficulties so run was abandoned";
    case 5:
      return "stopped by event handler";
    default:
      return "unknown";
    }
  }

  std::map<std::string, CbcModel::CbcIntParam> CbcInterface::param_map_int =  {
    {"MaxNumNode", CbcModel::CbcMaxNumNode},
    {"MaxNumSol", CbcModel::CbcMaxNumSol},
    {"FathomDiscipline", CbcModel::CbcFathomDiscipline},
    {"Printing", CbcModel::CbcPrinting},
    {"NumberBranches", CbcModel::CbcNumberBranches},
  };

  std::map<std::string, CbcModel::CbcDblParam> CbcInterface::param_map_double =  {
    {"IntegerTolerance", CbcModel::CbcIntegerTolerance},
    {"InfeasibilityWeight", CbcModel::CbcInfeasibilityWeight},
    {"CutoffIncrement", CbcModel::CbcCutoffIncrement},
    {"AllowableGap", CbcModel::CbcAllowableGap},
    {"AllowableFractionGap", CbcModel::CbcAllowableFractionGap},
    {"MaximumSeconds", CbcModel::CbcMaximumSeconds},
    {"CurrentCutoff", CbcModel::CbcCurrentCutoff},
    {"OptimizationDirection", CbcModel::CbcOptimizationDirection},
    {"CurrentObjectiveValue", CbcModel::CbcCurrentObjectiveValue},
    {"CurrentMinimizationObjectiveValue", CbcModel::CbcCurrentMinimizationObjectiveValue},
    {"StartSeconds", CbcModel::CbcStartSeconds},
    {"HeuristicGap", CbcModel::CbcHeuristicGap},
    {"HeuristicFractionGap", CbcModel::CbcHeuristicFractionGap},
    {"SmallestChange", CbcModel::CbcSmallestChange},
    {"SumChange", CbcModel::CbcSumChange},
    {"LargestChange", CbcModel::CbcLargestChange},
    {"SmallChange", CbcModel::CbcSmallChange}
  };

  std::map<std::string, OsiIntParam> CbcInterface::osi_param_map_int =  {
    {"MaxNumIteration", OsiMaxNumIteration},
    {"MaxNumIterationHotStart", OsiMaxNumIterationHotStart},
    {"NameDiscipline", OsiNameDiscipline}
  };

  std::map<std::string, OsiDblParam> CbcInterface::osi_param_map_double =  {
    {"DualObjectiveLimit", OsiDualObjectiveLimit},
    {"PrimalObjectiveLimit", OsiPrimalObjectiveLimit},
    {"DualTolerance", OsiDualTolerance},
    {"PrimalTolerance", OsiPrimalTolerance},
    {"ObjOffset", OsiObjOffset}
  };

  inline std::string return_secondary_status_string(int status) {
    switch (status) {
    case -1:
      return "unset";
    case 0:
      return "search completed with solution";
    case 1:
      return "linear relaxation not feasible (or worse than cutoff)";
    case 2:
      return "stopped on gap";
    case 3:
      return "stopped on nodes";
    case 4:
      return "stopped on time";
    case 5:
      return "stopped on user event";
    case 6:
      return "stopped on solutions";
    case CbcEventHandler::CbcEvent::node:
      return "node";
    case CbcEventHandler::CbcEvent::treeStatus:
      return "treeStatus";
    case CbcEventHandler::CbcEvent::solution:
      return "solution";
    case CbcEventHandler::CbcEvent::heuristicSolution:
      return "heuristicSolution";
    case CbcEventHandler::CbcEvent::beforeSolution1:
      return "beforeSolution1";
    case CbcEventHandler::CbcEvent::beforeSolution2:
      return "beforeSolution2";
    case CbcEventHandler::CbcEvent::afterHeuristic:
      return "afterHeuristic";
    case CbcEventHandler::CbcEvent::smallBranchAndBound:
      return "smallBranchAndBound";
    case CbcEventHandler::CbcEvent::heuristicPass:
      return "heuristicPass";
    case CbcEventHandler::CbcEvent::convertToCuts:
      return "convertToCuts";
    case CbcEventHandler::CbcEvent::endSearch:
      return "endSearch";
    default:
      return "unknown";
    }
  }

  class CasadiHandler : public CoinMessageHandler {
    public:
      virtual int print() ;
  };

  int CasadiHandler::print() {
    uout() << messageBuffer() << std::endl;
    return 0;
  }

  void CbcInterface::copy_cbc_results(const CbcModel& model, double** res) const {
    // Primal solution
    const double* x = model.getColSolution();
    casadi_copy(x, nx_, res[CONIC_X]);

    // Dual solution (x)
    const double* minus_lam_x = model.getReducedCost();
    if (res[CONIC_LAM_X]) {
      casadi_copy(minus_lam_x, nx_, res[CONIC_LAM_X]);
      casadi_scal(nx_, -1., res[CONIC_LAM_X]);
    }

    // Dual solution (A)
    const double* minus_lam_a = model.getRowPrice();
    if (res[CONIC_LAM_A]) {
      casadi_copy(minus_lam_a, na_, res[CONIC_LAM_A]);
      casadi_scal(na_, -1., res[CONIC_LAM_A]);
    }

    // Optimal cost
    double f = model.getObjValue();
    if (res[CONIC_COST]) *res[CONIC_COST] = f;
  }

  void CbcInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Conic::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="sos_groups") {
        sos_groups_ = to_int(op.second.to_int_vector_vector());
      } else if (op.first=="sos_weights") {
        sos_weights_ = op.second.to_double_vector_vector();
      } else if (op.first=="sos_types") {
        sos_types_ = op.second.to_int_vector();
      } else if (op.first=="hot_start") {
        hot_start_ = op.second;
      } else if (op.first=="cbc") {
        opts_ = op.second;
      }
    }

    // Validaty SOS constraints
    check_sos(nx_, sos_groups_, sos_weights_, sos_types_);

    // Default options
    casadi_assert(H_.nnz()==0, "Not an LP");

    // Allocate work vectors
    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nnz_in(CONIC_H), true); // H
    alloc_w(nnz_in(CONIC_A), true); // A
  }

  int CbcInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    if (!mem) return 1;
    auto m = static_cast<CbcMemory*>(mem);

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");

    m->colind.resize(A_.size2()+1);
    m->row.resize(A_.nnz());

    return 0;
  }

  int CbcInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<CbcMemory*>(mem);

    // Problem has not been solved at this point
    m->return_status = -1;
    m->secondary_return_status = -1;

    m->fstats.at("preprocessing").tic();

    // Get inputs
    double* g=w; w += nx_;
    casadi_copy(arg[CONIC_G], nx_, g);
    double* lbx=w; w += nx_;
    casadi_copy(arg[CONIC_LBX], nx_, lbx);
    double* ubx=w; w += nx_;
    casadi_copy(arg[CONIC_UBX], nx_, ubx);
    double* lba=w; w += na_;
    casadi_copy(arg[CONIC_LBA], na_, lba);
    double* uba=w; w += na_;
    casadi_copy(arg[CONIC_UBA], na_, uba);
    double* H=w; w += nnz_in(CONIC_H);
    casadi_copy(arg[CONIC_H], nnz_in(CONIC_H), H);
    double* A=w; w += nnz_in(CONIC_A);
    casadi_copy(arg[CONIC_A], nnz_in(CONIC_A), A);

    copy_vector(A_.colind(), m->colind);
    copy_vector(A_.row(), m->row);

    // Create osi_model
    OsiClpSolverInterface osi_model;

    osi_model.loadProblem(A_.size2(), A_.size1(), get_ptr(m->colind), get_ptr(m->row), A,
                      lbx, ubx, g, lba, uba);

    // Pass information on discreteness
    if (!discrete_.empty()) {
      for (casadi_int i=0; i<A_.size2();++i) {
        if (discrete_[i]) osi_model.setInteger(i);
      }
    }

    CbcModel model(osi_model);

    if (hot_start_) {
      model.setBestSolution(arg[CONIC_X0], nx_, COIN_DBL_MAX, true);

      // We store the result here already, because when CbcMain1 cannot do
      // better than setBestSolution() it will return a bogus result.
      copy_cbc_results(model, res);
    }

    // Construct SOS constraints
    std::vector<CbcSOS> sos_objects;
    for (casadi_int i=0;i<sos_groups_.size();++i) {
      const std::vector<int>& sos_group = sos_groups_[i];
      sos_objects.emplace_back(&model, sos_group.size(), get_ptr(sos_group),
                    sos_weights_.empty() ? nullptr : get_ptr(sos_weights_[i]), i, sos_types_[i]);
    }
    std::vector<CbcObject*> sos_objects_ptr;
    for (casadi_int i=0;i<sos_groups_.size();++i) {
      sos_objects_ptr.push_back(&sos_objects[i]);
    }
    if (!sos_objects.empty()) {
      model.addObjects(sos_objects.size(), get_ptr(sos_objects_ptr));
    }

    // Reset options
    CbcMain0(model);

    std::vector<std::string> main1_options(1, "CbcInterface");

    // Read Osi options
    for (auto&& op : opts_) {
      {
        // Check for double params
        auto it = param_map_double.find(op.first);
        if (it!=param_map_double.end()) {
          casadi_assert(model.setDblParam(it->second, op.second.to_double()),
            "Error setting option '" + op.first + "'.");
          continue;
        }
      }
      {
        // Check for integer params
        auto it = param_map_int.find(op.first);
        if (it!=param_map_int.end()) {
          casadi_assert(model.setIntParam(it->second, op.second.to_int()),
            "Error setting option '" + op.first + "'.");
          continue;
        }
      }
      {
        // Check for double params
        auto it = osi_param_map_double.find(op.first);
        if (it!=osi_param_map_double.end()) {
          casadi_assert(model.solver()->setDblParam(it->second, op.second.to_double()),
            "Error setting option '" + op.first + "'.");
          continue;
        }
      }
      {
        // Check for integer params
        auto it = osi_param_map_int.find(op.first);
        if (it!=osi_param_map_int.end()) {
          casadi_assert(model.solver()->setIntParam(it->second, op.second.to_int()),
            "Error setting option '" + op.first + "'.");
          continue;
        }
      }
      if (op.first=="startalg") {
        std::string startalg = op.second.to_string();
        main1_options.push_back("-" + op.second.to_string());
      } else {
        main1_options.push_back("-" + op.first);
        main1_options.push_back(str(op.second));
      }
    }

    main1_options.push_back("-solve");
    main1_options.push_back("-quit");

    std::vector<const char*> main_options_char;
    for (const auto& s : main1_options) main_options_char.push_back(s.c_str());

    CasadiHandler ch;
    model.passInMessageHandler(&ch);

    m->fstats.at("preprocessing").toc();
    m->fstats.at("solver").tic();

    CbcMain1(main_options_char.size(), get_ptr(main_options_char), model);

    m->fstats.at("solver").toc();
    m->fstats.at("postprocessing").tic();

    if (hot_start_ && model.status() == 0 &&
      model.isProvenOptimal() && model.secondaryStatus() == 1) {
      // Solution found by setBestSolution is best and only correct one.
    } else {
      copy_cbc_results(model, res);
    }

    m->fstats.at("postprocessing").toc();

    m->return_status = model.status();
    m->success = m->return_status==0 && model.isProvenOptimal() && model.secondaryStatus() <= 1;
    m->secondary_return_status = model.secondaryStatus();
    m->iter_count = model.getIterationCount();
    m->node_count = model.getNodeCount();
    if (m->return_status==1) m->unified_return_status = SOLVER_RET_LIMITED;

    if (verbose_) casadi_message("CBC return status: " + return_status_string(m->return_status));
    if (verbose_) casadi_message(
      "CBC secondary return status: " + return_secondary_status_string(m->secondary_return_status));

    return 0;
  }

  CbcInterface::~CbcInterface() {
    clear_mem();
  }

  CbcMemory::CbcMemory() {
  }

  CbcMemory::~CbcMemory() {
  }


  Dict CbcInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<CbcMemory*>(mem);
    stats["return_status"] = return_status_string(m->return_status);
    stats["secondary_return_status"] = return_secondary_status_string(m->secondary_return_status);
    stats["iter_count"] = m->iter_count;
    stats["node_count"] = m->node_count;
    return stats;
  }

  CbcInterface::CbcInterface(DeserializingStream& s) : Conic(s) {
    s.version("CbcInterface", 1);
    s.unpack("CbcInterface::opts", opts_);
    s.unpack("CbcInterface::sos_groups", sos_groups_);
    s.unpack("CbcInterface::sos_weights", sos_weights_);
    s.unpack("CbcInterface::sos_types", sos_types_);
    s.unpack("CbcInterface::hot_start", hot_start_);
  }

  void CbcInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("CbcInterface", 1);
    s.pack("CbcInterface::opts", opts_);
    s.pack("CbcInterface::sos_groups", sos_groups_);
    s.pack("CbcInterface::sos_weights", sos_weights_);
    s.pack("CbcInterface::sos_types", sos_types_);
    s.pack("CbcInterface::hot_start", hot_start_);
  }

} // end namespace casadi
