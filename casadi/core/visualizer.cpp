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


#include "function.hpp"
#include "casadi_meta.hpp"
#include "function_internal.hpp"
#include "sx_function.hpp"
#include "sx.hpp"
#include "serializer.hpp"
#include "mx_function.hpp"
#include "filesystem_impl.hpp"
#include <casadi/core/resource_casadi_viz.hpp>

#include <algorithm>
#include <iomanip>
#include <map>
#include <sstream>

namespace casadi {
namespace {

// Escape JSON and HTML script delimiters without interpreting expression labels.
std::string graph_string(const std::string& value) {
  std::ostringstream s;
  s << '"';
  for (unsigned char c : value) {
    if (c == '"' || c == '\\') {
      s << '\\' << c;
    } else if (c < 32 || c == '<' || c == '>' || c == '&') {
      s << "\\u" << std::hex << std::setfill('0') << std::setw(4)
        << static_cast<unsigned int>(c) << std::dec;
    } else {
      s << c;
    }
  }
  s << '"';
  return s.str();
}

void graph_sparsity(std::ostream& s, const Sparsity& sp) {
  s << "{\"shape\":[" << sp.size1() << "," << sp.size2() << "],\"colind\":[";
  for (casadi_int i = 0; i <= sp.size2(); ++i) {
    if (i) s << ",";
    s << sp.colind()[i];
  }
  s << "],\"row\":[";
  for (casadi_int i = 0; i < sp.nnz(); ++i) {
    if (i) s << ",";
    s << sp.row()[i];
  }
  s << "]}";
}

struct GraphNode {
  casadi_int id, op, io = -1, offset = 0;
  bool binary, ordered;
  std::string kind, display, symbol;
  std::vector<std::string> input_names, output_names, constants;
  std::vector<Sparsity> inputs, outputs;
};
struct GraphEdge { casadi_int from, output, to, input; };
struct GraphModel {
  std::vector<GraphNode> nodes;
  std::vector<GraphEdge> edges;
};

std::string graph_model(Function function, const std::string& direction,
    std::vector<Function>& functions, std::map<const FunctionInternal*, casadi_int>& indices,
    bool include_functions, GraphModel& model) {
  auto sx = dynamic_cast<const SXFunction*>(function.get());
  auto mx = dynamic_cast<const MXFunction*>(function.get());
  const auto sx_inputs = sx ? function.sx_in() : std::vector<SX>{};
  const auto mx_inputs = mx ? function.mx_in() : std::vector<MX>{};
  std::ostringstream data, edges;
  data << "{\"version\":1,\"name\":" << graph_string(function.name())
       << ",\"type\":" << graph_string(function.class_name())
       << ",\"direction\":" << graph_string(direction) << ",\"inputs\":[";
  for (casadi_int i = 0; i < function.n_in(); ++i) {
    if (i) data << ",";
    data << "{\"name\":" << graph_string(function.name_in(i)) << ",\"sparsity\":";
    graph_sparsity(data, function.sparsity_in(i));
    data << "}";
  }
  data << "],\"outputs\":[";
  for (casadi_int i = 0; i < function.n_out(); ++i) {
    if (i) data << ",";
    data << "{\"name\":" << graph_string(function.name_out(i)) << ",\"sparsity\":";
    graph_sparsity(data, function.sparsity_out(i));
    data << "}";
  }
  data << "],\"nodes\":[";

  // Work slots are reused: resolve dependencies before replacing their producers.
  std::map<casadi_int, std::pair<casadi_int, casadi_int>> producer;
  bool first_edge = true;
  for (casadi_int k = 0; k < function.n_instructions(); ++k) {
    casadi_int op = function.instruction_id(k);
    auto arg = function.instruction_input(k), res = function.instruction_output(k);
    std::string label = casadi_math<double>::name(op);
    std::string expression = sx ? sx->print(sx->algorithm_.at(k))
                               : mx->print(mx->algorithm_.at(k));
    GraphNode node;
    node.id = k; node.op = op;
    std::string kind = "operation", io_metadata, call_metadata;
    if (op == OP_INPUT || op == OP_OUTPUT) {
      bool input = op == OP_INPUT;
      casadi_int io = (input ? arg : res).at(0);
      casadi_int offset = sx ? (input ? arg : res).at(1)
                            : mx->algorithm_.at(k).data->offset();
      node.io = io; node.offset = offset;
      io_metadata = ",\"io_index\":" + str(io) + ",\"io_offset\":" + str(offset);
      if (input) {
        std::string symbol = sx ? str(sx_inputs.at(io).nonzeros().at(offset))
                                : str(mx_inputs.at(io));
        node.symbol = symbol;
        io_metadata += ",\"symbol\":" + graph_string(symbol);
      }
      label = (input ? function.name_in(io) : function.name_out(io));
      if ((sx && (input ? function.nnz_in(io) : function.nnz_out(io)) != 1) || offset != 0) {
        label += "[" + str(offset) + "]";
      }
      kind = input ? "input" : "output";
    } else if (op == OP_PARAMETER) {
      kind = "symbol";
      label = sx ? str(sx->free_vars_.at(sx->algorithm_.at(k).i1))
                 : mx->algorithm_.at(k).data.name();
    } else if (op == OP_CONST) {
      kind = "constant";
    } else if (op == OP_CALL) {
      kind = "call";
      label = (sx ? sx->call_.el.at(sx->algorithm_.at(k).i1).f.name()
                             : mx->algorithm_.at(k).data.which_function().name());
    }
    if (op == OP_INPUT) arg.clear();
    if (op == OP_OUTPUT) res.clear();
    std::vector<std::string> input_names, output_names, constants;
    for (casadi_int i = 0; i < arg.size(); ++i) input_names.push_back("arg" + str(i));
    for (casadi_int i = 0; i < res.size(); ++i) output_names.push_back("out" + str(i));
    if (op == OP_CALL) {
      const Function& f = sx ? sx->call_.el.at(sx->algorithm_.at(k).i1).f
                             : mx->algorithm_.at(k).data.which_function();
      call_metadata = ",\"callee_type\":" + graph_string(f.class_name());
      Dict info = f.info();
      auto source = info.find("model_path");
      if (source != info.end()) {
        std::string path = source->second.to_string();
        if (!path.empty()) {
          call_metadata += ",\"model_path\":" + graph_string(path);
          label += " [" + path.substr(path.find_last_of("/\\") + 1) + "]";
        }
      }
      if (include_functions && (f.is_a("SXFunction") || f.is_a("MXFunction"))) {
        auto found = indices.find(f.get());
        if (found == indices.end()) {
          found = indices.emplace(f.get(), functions.size()).first;
          functions.push_back(f);
        }
        call_metadata += ",\"callee\":" + str(found->second);
      }
      input_names.clear(); output_names.clear();
      for (casadi_int side = 0; side < 2; ++side) {
        auto& names = side == 0 ? input_names : output_names;
        for (casadi_int i = 0; i < (side == 0 ? f.n_in() : f.n_out()); ++i) {
          std::string name = side == 0 ? f.name_in(i) : f.name_out(i);
          casadi_int count = side == 0 ? f.nnz_in(i) : f.nnz_out(i);
          if (sx) {
            for (casadi_int j = 0; j < count; ++j) {
              names.push_back(name + (count == 1 ? "" : "[" + str(j) + "]"));
            }
          } else {
            names.push_back(name);
          }
        }
      }
    } else if (op == OP_MTIMES) {
      input_names = {"accumulator", "left", "right"};
    } else if (op == OP_GETNONZEROS) {
      input_names = {"source"};
    } else if (op == OP_SETNONZEROS || op == OP_ADDNONZEROS) {
      input_names = {"base", "values"};
    }
    std::string display = label, formula = label;
    if (kind == "operation") {
      if (mx) {
        formula = MX::print_operator(mx->algorithm_.at(k).data, input_names);
      } else {
        formula = casadi_math<double>::pre(op);
        for (casadi_int i = 0; i < input_names.size(); ++i) {
          if (i) formula += casadi_math<double>::sep(op);
          formula += input_names[i];
        }
        formula += casadi_math<double>::post(op);
      }
      switch (op) {
        case OP_ADD: display = "+"; break;
        case OP_SUB: display = "-"; break;
        case OP_MUL: display = "*"; break;
        case OP_DIV: display = "/"; break;
        case OP_NEG: display = "-"; break;
        case OP_SQ: display = "(.)^2"; break;
        case OP_POW: case OP_CONSTPOW: display = "pow"; break;
        case OP_TWICE: display = "2*(.)"; break;
        case OP_INV: display = "1/(.)"; break;
        default: break;
      }
    } else if (op == OP_CONST) {
      std::vector<double> values = sx ? std::vector<double> {sx->algorithm_.at(k).d}
                                      : static_cast<DM>(mx->algorithm_.at(k).data).nonzeros();
      for (double value : values) {
        std::ostringstream number;
        number << std::setprecision(17) << value;
        constants.push_back(number.str());
      }
      display = constants.size() == 1 ? constants.front() : "constant";
      formula = display;
    }
    node.kind = kind; node.display = display;
    node.binary = casadi_math<double>::is_binary(op);
    node.ordered = arg.size() > 1 && !operation_checker<CommChecker>(op);
    node.input_names = input_names; node.output_names = output_names;
    node.constants = constants;
    if (k) data << ",";
    data << "{\"id\":" << k << ",\"op\":" << op << ",\"label\":"
         << graph_string(label) << ",\"kind\":" << graph_string(kind)
         << ",\"expression\":" << graph_string(expression)
         << ",\"display\":" << graph_string(display)
         << ",\"formula\":" << graph_string(formula)
         << ",\"binary\":" << (casadi_math<double>::is_binary(op) ? "true" : "false")
         << ",\"ordered\":" << (arg.size() > 1 && !operation_checker<CommChecker>(op)
                                         ? "true" : "false");
    for (const auto& entry : std::map<std::string, std::vector<std::string>>{
        {"input_names", input_names}, {"output_names", output_names}, {"constants", constants}}) {
      data << "," << graph_string(entry.first) << ":[";
      for (casadi_int i = 0; i < entry.second.size(); ++i) {
        if (i) data << ",";
        data << graph_string(entry.second[i]);
      }
      data << "]";
    }
    if (mx && (op == OP_GETNONZEROS || op == OP_SETNONZEROS || op == OP_ADDNONZEROS)) {
      data << ",\"mapping_kind\":" << graph_string(op == OP_GETNONZEROS ? "extract"
        : op == OP_SETNONZEROS ? "assign" : "add") << ",\"mapping\":[";
      const auto mapping = mx->algorithm_.at(k).data.mapping().nonzeros();
      for (casadi_int i = 0; i < mapping.size(); ++i) {
        if (i) data << ",";
        data << mapping[i];
      }
      data << "]";
    }
    data << io_metadata << call_metadata << ",\"inputs\":[";
    for (casadi_int i = 0; i < arg.size(); ++i) {
      if (i) data << ",";
      node.inputs.push_back(sx ? Sparsity::scalar() : mx->algorithm_.at(k).data->dep(i).sparsity());
      graph_sparsity(data, node.inputs.back());
      if (arg[i] < 0) continue;
      auto p = producer.find(arg[i]);
      casadi_assert(p != producer.end(), "Missing graph producer for instruction " + str(k));
      model.edges.push_back({p->second.first, p->second.second, k, i});
      if (!first_edge) edges << ",";
      first_edge = false;
      edges << "{\"from\":" << p->second.first << ",\"output\":" << p->second.second
            << ",\"to\":" << k << ",\"input\":" << i << "}";
    }
    data << "],\"outputs\":[";
    for (casadi_int i = 0; i < res.size(); ++i) {
      if (i) data << ",";
      node.outputs.push_back(sx ? Sparsity::scalar() : mx->algorithm_.at(k).data->sparsity(i));
      graph_sparsity(data, node.outputs.back());
      if (res[i] >= 0) producer[res[i]] = {k, i};
    }
    data << "]}";
    model.nodes.push_back(std::move(node));
  }
  data << "],\"edges\":[" << edges.str() << "]}";

  return data.str();
}

// DOT quoted strings and HTML labels have different escaping rules from JSON.
std::string dot_quote(const std::string& value) {
  std::string out = "\"";
  for (char c : value) {
    if (c == '\n') {
      out += "\\n";
    } else if (c == '\r') {
      out += "\\r";
    } else {
      if (c == '\\' || c == '"') out += '\\';
      out += c;
    }
  }
  return out + "\"";
}
std::string dot_html(const std::string& value) {
  std::string out;
  for (char c : value) {
    switch (c) {
      case '&': out += "&amp;"; break;
      case '<': out += "&lt;"; break;
      case '>': out += "&gt;"; break;
      case '"': out += "&quot;"; break;
      default: out += c;
    }
  }
  return out;
}
std::string graph_size(const Sparsity& sp) {
  return str(sp.size1()) + "-by-" + str(sp.size2());
}
void dot_matrix(std::ostream& s, const std::string& id, const Sparsity& sp,
    const std::string& title, const std::vector<std::string>* constants = nullptr,
    const std::string& color = "#666666", bool show_sizes = true) {
  s << id << " [shape=plain, fontcolor=\"#666666\", fontsize=10, label=<"
    << "<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"3\">";
  if (show_sizes || !title.empty()) s << "<TR><TD COLSPAN=\""
    << std::max<casadi_int>(1, std::min<casadi_int>(8, sp.size2()))
    << "\"><FONT COLOR=\"" << color << "\">"
    << dot_html(show_sizes ? (title.empty() ? graph_size(sp) : title + " : " + graph_size(sp))
      : title)
    << "</FONT></TD></TR>";
  std::map<std::pair<casadi_int, casadi_int>, casadi_int> lookup;
  for (casadi_int c = 0; c < std::min<casadi_int>(8, sp.size2()); ++c) {
    for (casadi_int k = sp.colind()[c]; k < sp.colind()[c+1]; ++k) {
      if (sp.row()[k] < 8) lookup[{sp.row()[k], c}] = k;
    }
  }
  for (casadi_int r = 0; r < std::min<casadi_int>(8, sp.size1()); ++r) {
    s << "<TR>";
    for (casadi_int c = 0; c < std::min<casadi_int>(8, sp.size2()); ++c) {
      auto nz = lookup.find({r, c});
      s << "<TD BGCOLOR=\"white\"";
      if (nz != lookup.end()) s << " PORT=\"nz" << nz->second << "\"";
      s << ">";
      if (constants) s << (nz == lookup.end() ? "." : dot_html(constants->at(nz->second)));
      else
        s << "<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLSPACING=\"0\" CELLPADDING=\"0\">"
        << "<TR><TD FIXEDSIZE=\"TRUE\" WIDTH=\"9\" HEIGHT=\"9\" BGCOLOR=\""
        << (nz == lookup.end() ? "white" : "#111111") << "\"></TD></TR></TABLE>";
      s << "</TD>";
    }
    if (!sp.size2()) s << "<TD>empty</TD>";
    s << "</TR>";
  }
  if (!sp.size1()) s << "<TR><TD>empty</TD></TR>";
  if (sp.size1() > 8 || sp.size2() > 8) s << "<TR><TD>...</TD></TR>";
  s << "</TABLE>>];\n";
}
bool dot_record(const GraphNode& n) {
  return (n.ordered && !n.binary) || n.kind == "call";
}
bool dot_has_table(const GraphNode& n, const Sparsity& sp, bool show_contents) {
  if (!show_contents || sp.numel() == 1) return false;
  if (n.kind == "call" || n.inputs.empty() || n.op == OP_GETNONZEROS
      || n.op == OP_SETNONZEROS || n.op == OP_ADDNONZEROS) return true;
  for (const auto& input : n.inputs) if (input != sp) return true;
  return false;
}
void dot_ports(std::ostream& s, const GraphNode& n,
    const std::vector<std::string>& names, const std::string& prefix, bool show_sizes) {
  if (names.empty()) return;
  s << "<TR><TD><TABLE BORDER=\"0\" CELLBORDER=\"1\" COLOR=\"#d26666\"><TR>";
  for (casadi_int i = 0; i < names.size(); ++i) {
    s << "<TD PORT=\"" << prefix << i << "\"><FONT POINT-SIZE=\"10\">"
      << dot_html(names[i]);
    const auto& sp = prefix == "in" ? n.inputs.at(i) : n.outputs.at(i);
    if (show_sizes && sp.numel() != 1) s << " " << graph_size(sp);
    s << "</FONT></TD>";
  }
  s << "</TR></TABLE></TD></TR>";
}
void graph_dot(std::ostream& s, const Function& f, const GraphModel& model,
    const std::string& direction, const std::string& view, bool show_contents, bool show_sizes) {
  bool function_view = view == "function";
  std::map<std::pair<casadi_int, casadi_int>, std::string> sources;
  std::map<casadi_int, std::string> targets;
  s << "digraph G {\ngraph [rankdir=" << (function_view ? "TB" : direction)
    << ", bgcolor=\"white\", pad=0.4, nodesep=0.55, ranksep=0.65];\n"
    << "node [shape=ellipse, style=filled, color=\"#b00000\", fillcolor=\"#b00000\", "
    << "fontcolor=white, fontname=Helvetica, fontsize=14, margin=\"0.15,0.08\", "
    << "width=0.5, height=0.5];\nedge [color=\"#34658b\", penwidth=1.5, arrowsize=0.65];\n";
  for (const auto& n : model.nodes) {
    if (n.kind == "output" || (function_view && n.kind == "input")) continue;
    const std::string id = "n" + str(n.id);
    std::string display = !function_view && n.kind == "input" ? n.symbol : n.display;
    std::string color = n.kind == "input" || n.kind == "symbol" ? "#34658b"
      : n.kind == "constant" ? "#38754d" : "#b00000";
    s << id << " [color=" << dot_quote(color) << ", fillcolor=" << dot_quote(color);
    if (dot_record(n)) {
      s << ", shape=plain, label=<<TABLE BORDER=\"0\" CELLBORDER=\"0\" "
        << "CELLSPACING=\"0\" CELLPADDING=\"6\" BGCOLOR=\"#b00000\">";
      dot_ports(s, n, n.input_names, "in", show_sizes);
      s << "<TR><TD>" << dot_html(display) << "</TD></TR>";
      if (n.kind == "call") dot_ports(s, n, n.output_names, "out", show_sizes);
      s << "</TABLE>>";
    } else {
      const auto& ports = n.outputs.empty() ? n.inputs : n.outputs;
      bool table = false;
      for (const auto& sp : n.outputs) table = table || dot_has_table(n, sp, show_contents);
      if (show_sizes && !table && !ports.empty() && ports.front().numel() != 1) {
        display += "\n" + graph_size(ports.front());
      }
      s << ", label=" << dot_quote(display);
    }
    s << "];\n";
    for (casadi_int i = 0; i < n.outputs.size(); ++i) {
      const auto& sp = n.outputs[i];
      if (!dot_has_table(n, sp, show_contents)) continue;
      std::string value = "v" + str(n.id) + "_" + str(i);
      dot_matrix(s, value, sp, n.kind == "call" ? n.output_names.at(i) : "",
        n.kind == "constant" ? &n.constants : nullptr, "#666666", show_sizes);
      s << id << (n.kind == "call" ? ":out" + str(i) : "") << " -> " << value << ";\n";
      sources[{n.id, i}] = value;
    }
  }
  for (bool input : {true, false}) {
    std::vector<std::string> row;
    for (casadi_int i = 0; i < (input ? f.n_in() : f.n_out()); ++i) {
      const auto& sp = input ? f.sparsity_in(i) : f.sparsity_out(i);
      if (!function_view && (input || !show_contents || !f.is_a("SXFunction")
          || sp.numel() == 1)) continue;
      std::string kind = input ? "input" : "output", id = "b" + kind + "_" + str(i);
      std::string name = function_view ? (input ? f.name_in(i) : f.name_out(i)) : "";
      if (show_contents && sp.numel() != 1) dot_matrix(s, id, sp, name, nullptr,
        function_view ? (input ? "#34658b" : "#b00000") : "#666666", show_sizes);
      else
        s << id << " [label=" << dot_quote(name
        + (show_sizes && sp.numel() != 1 ? "\n" + graph_size(sp) : "")) << ", color=\""
        << (input ? "#34658b" : "#b00000") << "\", fillcolor=\""
        << (input ? "#34658b" : "#b00000") << "\"];\n";
      for (const auto& n : model.nodes) {
        if (n.kind != kind || n.io != i) continue;
        std::string endpoint = id;
        if (show_contents && f.is_a("SXFunction") && sp.numel() != 1 && n.offset < sp.nnz()) {
          const auto col = std::upper_bound(sp.colind(), sp.colind()+sp.size2()+1, n.offset)
            - sp.colind() - 1;
          if (col < 8 && sp.row()[n.offset] < 8) endpoint += ":nz" + str(n.offset);
        }
        if (input) sources[{n.id, 0}] = endpoint;
        else
          targets[n.id] = endpoint;
      }
      row.push_back(id);
    }
    if (function_view && !row.empty()) {
      s << "{rank=" << (input ? "source" : "sink") << ";";
      for (const auto& id : row) s << id << ";";
      s << "}\n";
      for (casadi_int i = 1; i < row.size(); ++i) {
        s << row[i-1] << " -> " << row[i] << " [style=invis, weight=100];\n";
      }
    }
  }
  for (const auto& e : model.edges) {
    const auto& from = model.nodes.at(e.from);
    const auto& to = model.nodes.at(e.to);
    auto source = sources.find({e.from, e.output});
    auto target = targets.find(e.to);
    if (to.kind == "output" && target == targets.end()) continue;
    s << (source != sources.end() ? source->second : "n" + str(e.from)
      + (from.kind == "call" ? ":out" + str(e.output) : "")) << " -> ";
    if (target != targets.end()) s << target->second;
    else
      s << "n" << e.to << (to.binary ? (e.input == 0 ? ":w" : ":e")
      : dot_record(to) ? ":in" + str(e.input) : "");
    s << ";\n";
  }
  s << "}\n";
}

struct GraphOptions {
  std::string direction = "TB", viz_js, view = "function";
  std::string viewer_url = "https://unpkg.com/@casadi/casadi-viz@"
    + std::string(CasadiMeta::version()).substr(
      0, std::string(CasadiMeta::version()).find_last_of('.'))
    + "/dist/index.js";
  bool include_functions = true, show_matrix_contents = true, show_matrix_sizes = true;
  explicit GraphOptions(const Dict& opts) {
    for (const auto& opt : opts) {
      if (opt.first == "direction") {
        direction = opt.second.to_string();
      } else if (opt.first == "view") {
        view = opt.second.to_string();
      } else if (opt.first == "include_functions") {
        include_functions = opt.second.to_bool();
      } else if (opt.first == "show_matrix_contents") {
        show_matrix_contents = opt.second.to_bool();
      } else if (opt.first == "show_matrix_sizes") {
        show_matrix_sizes = opt.second.to_bool();
      } else if (opt.first == "viewer_url") {
        viewer_url = opt.second.to_string();
      } else if (opt.first == "viz_js") {
        viz_js = opt.second.to_string();
      } else {
        casadi_error("Unknown export_graph option: " + opt.first);
      }
    }
    casadi_assert(direction == "LR" || direction == "RL" || direction == "TB"
      || direction == "BT", "Invalid export_graph direction: " + direction);

    casadi_assert(view == "function" || view == "expression",
      "Invalid export_graph view: " + view);

  }
};

std::string graph_bundle(const Function& f, const GraphOptions& opts, bool dot = false) {
  casadi_assert(f.is_a("SXFunction") || f.is_a("MXFunction"),
    "export_graph requires an SXFunction or MXFunction");
  std::vector<Function> functions{f};
  std::map<const FunctionInternal*, casadi_int> indices{{f.get(), 0}};
  std::vector<std::string> models;
  for (casadi_int i = 0; i < functions.size(); ++i) {
    GraphModel model;
    models.push_back(graph_model(functions[i], opts.direction, functions, indices,
      opts.include_functions, model));
    if (dot) {
      std::ostringstream output;
      graph_dot(output, f, model, opts.direction, opts.view,
        opts.show_matrix_contents, opts.show_matrix_sizes);
      return output.str();
    }
  }
  std::ostringstream data;
  data << models.front().substr(0, models.front().size()-1)
       << ",\"format\":\"casadi_viz\",\"view\":"
       << graph_string(opts.view) << ",\"casadi_version\":"
       << graph_string(CasadiMeta::version()) << ",\"functions\":[";
  for (casadi_int i = 1; i < models.size(); ++i) {
    if (i > 1) data << ",";
    data << models[i];
  }
  data << "]}";

  return data.str();
}

void write_graph(const std::string& fname, const std::string& data,
    const GraphOptions& options) {
  if (fname.size() < 5 || fname.substr(fname.size()-5) != ".html") {
    auto output = Filesystem::ofstream_ptr(fname);
    *output << data << "\n";
    output->flush();
    casadi_assert(output->good(), "Failed to write graph to '" + fname + "'");
    return;
  }

  std::ostringstream config;
  config << "{\"viewer_url\":" << graph_string(options.viewer_url) << ",\"runtime\":";
  if (options.viz_js.empty()) {
    config << "null";
  } else {
    auto input = Filesystem::ifstream_ptr(options.viz_js);
    std::ostringstream contents;
    contents << input->rdbuf();
    config << "{\"source\":" << graph_string(contents.str()) << "}";
  }

  config << "}";

  const std::string html = resource_casadi_viz;
  const std::string graph_marker = "@@GRAPH_DATA@@", viz_marker = "@@VIZ_CONFIG@@";
  size_t graph_pos = html.find(graph_marker), viz_pos = html.find(viz_marker);
  casadi_assert_dev(graph_pos != std::string::npos && viz_pos > graph_pos);
  auto output = Filesystem::ofstream_ptr(fname);
  *output << html.substr(0, graph_pos) << data
          << html.substr(graph_pos + graph_marker.size(), viz_pos-graph_pos-graph_marker.size())
          << config.str() << html.substr(viz_pos + viz_marker.size());
  output->flush();
  casadi_assert(output->good(), "Failed to write graph to '" + fname + "'");
}

std::string graph_extension(const std::string& fname) {
  const auto dot = fname.find_last_of('.');
  const std::string extension = dot == std::string::npos ? "" : fname.substr(dot);
  casadi_assert(extension == ".html" || extension == ".dot" || extension == ".casadi_viz",
    "Unsupported export_graph extension '" + extension
    + "'. Supported extensions: .html, .dot, .casadi_viz");
  return extension;
}

template<typename MatType>
std::string expression_bundle(const std::vector<MatType>& expressions, const Dict& opts) {
  Dict defaults = opts;
  if (!defaults.count("view")) defaults["view"] = "expression";
  const GraphOptions options(defaults);
  StringSerializer serializer;
  serializer.pack(expressions);
  std::ostringstream data;
  data << "{\"format\":\"casadi_viz\",\"version\":1,\"source\":"
       << graph_string(serializer.encode()) << ",\"view\":" << graph_string(options.view)
       << ",\"direction\":" << graph_string(options.direction)
       << ",\"casadi_version\":" << graph_string(CasadiMeta::version())
       << ",\"include_functions\":" << (options.include_functions ? "true" : "false")
       << ",\"show_matrix_contents\":" << (options.show_matrix_contents ? "true" : "false")
       << ",\"show_matrix_sizes\":" << (options.show_matrix_sizes ? "true" : "false") << "}";
  return data.str();
}

template<typename MatType>
void export_expressions(const std::vector<MatType>& expressions, const std::string& fname,
    const Dict& opts) {
  if (graph_extension(fname) == ".dot") {
    Dict defaults = opts;
    if (!defaults.count("view")) defaults["view"] = "expression";
    Function("expression", symvar(veccat(expressions)), expressions).export_graph(fname, defaults);
  } else {
    write_graph(fname, expression_bundle(expressions, opts), GraphOptions(opts));
  }
}

} // namespace

std::string Function::export_graph(const Dict& opts) const {
  return graph_bundle(*this, GraphOptions(opts));
}

void Function::export_graph(const std::string& fname, const Dict& opts) const {
  const std::string extension = graph_extension(fname);
  const GraphOptions options(opts);
  write_graph(fname, graph_bundle(*this, options, extension == ".dot"), options);
}

void export_graph(const SX& expression, const std::string& fname, const Dict& opts) {
  export_expressions(std::vector<SX>{expression}, fname, opts);
}

void export_graph(const std::vector<SX>& expressions, const std::string& fname, const Dict& opts) {
  export_expressions(expressions, fname, opts);
}

void export_graph(const MX& expression, const std::string& fname, const Dict& opts) {
  export_expressions(std::vector<MX>{expression}, fname, opts);
}

void export_graph(const std::vector<MX>& expressions, const std::string& fname, const Dict& opts) {
  export_expressions(expressions, fname, opts);
}

std::string export_graph(const std::vector<SX>& expressions, const Dict& opts) {
  return expression_bundle(expressions, opts);
}

std::string export_graph(const SX& expression, const Dict& opts) {
  return export_graph(std::vector<SX>{expression}, opts);
}

std::string export_graph(const std::vector<MX>& expressions, const Dict& opts) {
  return expression_bundle(expressions, opts);
}

std::string export_graph(const MX& expression, const Dict& opts) {
  return export_graph(std::vector<MX>{expression}, opts);
}

} // namespace casadi
