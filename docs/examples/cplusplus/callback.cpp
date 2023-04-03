/*
 *    MIT No Attribution
 *
 *    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
   // Create a callback
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
