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


#include <casadi/casadi.hpp>

using namespace casadi;
/**
 * How to use Callback
 * Joel Andersson
 */

 class MyCallback : public Callback {
 public:
   // Constructor
   MyCallback(double d) : d(d) {
     construct("f");
   }

   // Destructor
   ~MyCallback() override { std::cout << "MyCallback is destroyed here." << std::endl; };

   // Initialize the object
   void init() override {
     std::cout << "initializing object" << std::endl;
   }

   // Number of inputs and outputs
   casadi_int get_n_in() override { return 1;}
   casadi_int get_n_out() override { return 1;}

   // Evaluate numerically
   std::vector<DM> eval(const std::vector<DM>& arg) const override {
     DM x = arg.at(0);
     DM f = sin(d*x);
     return {f};
   }

 private:
   // Data members
   double d;
 };

 int main() {
   // Create a callvack
   MyCallback cb(0.5);
   // Create a function object
   Function f = cb;
   std::vector<DM> arg = {2};
   std::vector<DM> res = f(arg);
   std::cout << res << std::endl;

   // Single reference
   std::cout << "Let's overwrite f here." << std::endl;
   f = Function();
   std::cout << "Done." << std::endl;

   // Function reference
   f = cb;
   Function g = f;
   f = Function();
   std::cout << "Let's overwrite g here." << std::endl;
   g = Function();
   std::cout << "Done." << std::endl;

   // Embedding in a graph
   f = cb;
   MX x = MX::sym("x");
   g = Function("g",{x},{f(x)});
   f = Function();
   std::cout << "Let's overwrite g here." << std::endl;
   g = Function();
   std::cout << "Done." << std::endl;

   return 0;
 }
