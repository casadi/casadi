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
 private:
   // Data members
   double d;
   // Private constructor
   MyCallback(double d) : d(d) {}
 public:
   // Creator function, creates an owning reference
   static Function create(const std::string& name, double d,
                          const Dict& opts=Dict()) {
     return Callback::create(name, new MyCallback(d), opts);
   }

   // Initialize the object
   virtual void init() {
     std::cout << "initializing object" << std::endl;
   }

   // Number of inputs and outputs
   virtual int get_n_in() { return 1;}
   virtual int get_n_out() { return 1;}

   // Evaluate numerically
   virtual std::vector<DM> eval(const std::vector<DM>& arg) {
     DM x = arg.at(0);
     DM f = sin(d*x);
     return {f};
   }
 };

 int main() {
   Function f = MyCallback::create("f", 0.5);
   std::vector<DM> arg = {2};
   std::vector<DM> res = f(arg);
   std::cout << res << std::endl;
   return 0;
 }
