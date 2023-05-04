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


#include <casadi/casadi.hpp>

using namespace casadi;
/**
 *  Solve an NLP using codegen  
 *  Part 2: solution
 *  Joel Andersson, KU Leuven 2013
 */

int main(){
    
  // Create an NLP solver passing derivative information
  Function solver = nlpsol("solver", "ipopt", "nlp.casadi");

  // Bounds and initial guess
  std::map<std::string, DM> arg = {{"lbx", -inf},
                                   {"ubx",  inf},
                                   {"lbg",    0},
                                   {"ubg",  inf},
                                   {"x0",     0}};

  // Solve the NLP
  auto res = solver(arg);

  // Print solution
  std::cout << "-----" << std::endl;
  std::cout << "objective at solution = " << res.at("f") << std::endl;
  std::cout << "primal solution = " << res.at("x") << std::endl;
  std::cout << "dual solution (x) = " << res.at("lam_x") << std::endl;
  std::cout << "dual solution (g) = " << res.at("lam_g") << std::endl;
  
  return 0;
}

