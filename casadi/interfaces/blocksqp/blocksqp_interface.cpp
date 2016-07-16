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


#include "blocksqp_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_BLOCKSQP_EXPORT
  casadi_register_nlpsol_blocksqp(Nlpsol::Plugin* plugin) {
    plugin->creator = BlocksqpInterface::creator;
    plugin->name = "blocksqp";
    plugin->doc = BlocksqpInterface::meta_doc.c_str();
    plugin->version = 30;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_BLOCKSQP_EXPORT casadi_load_nlpsol_blocksqp() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_blocksqp);
  }

  BlocksqpInterface::BlocksqpInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }


  BlocksqpInterface::~BlocksqpInterface() {
    clear_memory();
  }

  Options BlocksqpInterface::options_
  = {{&Nlpsol::options_},
     {{"blocksqp",
       {OT_DICT,
        "Options to be passed to BLOCKSQP"}},
      {"contype",
       {OT_INTVECTOR,
        "Type of constraint"}}
     }
  };

  void BlocksqpInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="blocksqp") {
        //opts_ = op.second;
      } else if (op.first=="contype") {
        //contype_ = op.second;
      }
    }
  }

  void BlocksqpInterface::init_memory(void* mem) const {
    Nlpsol::init_memory(mem);
    auto m = static_cast<BlocksqpMemory*>(mem);
  }

  void BlocksqpInterface::set_work(void* mem, const double**& arg, double**& res,
                                 int*& iw, double*& w) const {
    auto m = static_cast<BlocksqpMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);
    casadi_error("Not implemented");
  }

  void BlocksqpInterface::solve(void* mem) const {
    auto m = static_cast<BlocksqpMemory*>(mem);
    casadi_error("Not implemented");
  }

} // namespace casadi
