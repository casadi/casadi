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


#include "onnx_translator.hpp"

const std::string casadi::OnnxTranslator::meta_doc =
  "ONNX translator for importing and exporting ONNX computational graphs.\n\n"
  "Supports:\n"
  "  - Loading ONNX models from file\n"
  "  - Converting CasADi Functions to ONNX\n"
  "  - Setting dimensions for symbolic inputs\n"
  "  - Creating CasADi Functions from ONNX graphs\n"
  "  - Saving ONNX models to file\n"
  "  - Control flow with subgraphs (If operator)\n"
  "  - Outer scope variable access in subgraphs\n\n"
  "Supported Operations:\n"
  "  Binary: Add, Sub, Mul, Div, Pow\n"
  "  Trigonometric: Sin, Cos, Tan\n"
  "  Inverse Trigonometric: Asin, Acos, Atan\n"
  "  Hyperbolic: Sinh, Cosh, Tanh\n"
  "  Inverse Hyperbolic: Asinh, Acosh, Atanh\n"
  "  Exponential/Logarithmic: Exp, Log, Sqrt\n"
  "  Rounding: Ceil, Floor\n"
  "  Other: Abs, Sign, Neg, Erf\n"
  "  Utilities: Identity, Constant\n"
  "  Tensor: MatMul, Split, Slice (partial), GatherElements (partial)\n"
  "  Control Flow: If, Loop, Scan\n\n"
  "Control Flow Operations:\n"
  "  If: Conditional execution with then_branch and else_branch subgraphs.\n"
  "      Subgraphs can access variables from outer scope automatically.\n"
  "      Maps to CasADi Function::if_else for efficient evaluation.\n"
  "      Supports nested If operators.\n\n"
  "  Loop: Iterative execution with state (like while/for loops).\n"
  "        Body receives iteration_num, condition, and loop-carried dependencies.\n"
  "        Body outputs condition, updated dependencies, and scan outputs.\n"
  "        Maps to CasADi Function::mapaccum for stateful iteration.\n"
  "        Currently executes fixed iterations (conditional termination TODO).\n\n"
  "  Scan: Functional iteration with NO outer scope access.\n"
  "        All data must be passed explicitly as inputs.\n"
  "        Maps to CasADi Function::map for parallel/serial execution.\n"
  "        Simpler and more restrictive than Loop.\n\n"
  "Example (If operator):\n"
  "  x = MX.sym('x')\n"
  "  condition = x > 0\n"
  "  f_then = Function('then', [x], [x * 2])\n"
  "  f_else = Function('else', [x], [x / 2])\n"
  "  f_if = Function.if_else('conditional', f_then, f_else)\n"
  "  y = f_if(condition, [x])[0]\n"
  "  f = Function('example', [x], [y])\n"
  "  # Export to ONNX\n"
  "  t = translator('onnx')\n"
  "  t.load(f)\n"
  "  t.save('example.onnx')\n";
