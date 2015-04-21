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


/**
 *  Example program demonstrating the usage of C-code generated from CasADi
 *  Generated code is encapsulated in self-contained code with the entry points
 *  defined below. 
 *  Note that other usage, e.g. accessing the internal data structures in the 
 *  generated files is not recommended and subject to change.
 *
 *  We show how the generated code can be used from C (or C++), without requiring
 *  any CasADi classes as well as how to use it from C++ using CasADi's ExternalFunction.
 *
 *  Joel Andersson, 2013-2015
 */

#include <stdio.h>
#include <dlfcn.h>

// Usage from C
int usage_c(){
  printf("---\n");
  printf("Usage from C.\n");
  printf("\n");

  /* Signature of the entry points */
  typedef int (*nargPtr)(int *n_in, int *n_out);
  typedef int (*sparsityPtr)(int n_in, int *n_col, int *n_row, int **colind, int **row);
  typedef int (*workPtr)(int *ni, int *nr);
  typedef int (*evalPtr)(const double* const* arg, double* const* res, int* iw, double* w);

  /* Handle to the dll */
  void* handle;

  /* Load the dll */
  handle = dlopen("./f.so", RTLD_LAZY);  
  if(handle==0){
    printf("Cannot open f.so, error %s\n", dlerror()); 
    return 1;
  }

  /* Reset error */
  dlerror(); 

  /* Function for initialization and getting the number of inputs and outputs */
  nargPtr narg = (nargPtr)dlsym(handle, "f_narg");
  if(dlerror()){
    printf("Failed to retrieve \"narg\" function.\n");
    return 1;
  }

  /* Function for getting the sparsities of the inputs and outptus */
  sparsityPtr sparsity = (sparsityPtr)dlsym(handle, "f_sparsity");
  if(dlerror()){
    printf("Failed to retrieve \"sparsity\" function.\n");
    return 1;
  }

  /* Function for getting the size of the required work vectors */
  workPtr work = (workPtr)dlsym(handle, "f_work");
  if(dlerror()){
    printf("Failed to retrieve \"work\" function.\n");
    return 1;
  }

  /* Function for numerical evaluation */
  evalPtr eval = (evalPtr)dlsym(handle, "f");
  if(dlerror()){
    printf("Failed to retrieve \"f\" function.\n");
    return 1;
  }

  /* Initialize and get the number of inputs and outputs */
  int n_in=-1, n_out=-1;
  narg(&n_in, &n_out);
  printf("n_in = %d, n_out = %d\n", n_in, n_out);

  /* Get sparsities */
  int ind;
  for(ind=0; ind<n_in + n_out; ++ind){
    if(ind<n_in){
      printf("Input %d\n",ind);
    } else {
      printf("Output %d\n",ind-n_in);
    }

    int nrow,ncol;
    int *colind, *row;
    sparsity(ind,&nrow,&ncol,&colind,&row);

    printf("  Dimension: %d-by-%d\n",nrow,ncol);
    printf("  Nonzeros: {");
    int rr,cc,el;
    for(cc=0; cc<ncol; ++cc){ /* loop over rows */
      for(el=colind[cc]; el<colind[cc+1]; ++el){ /* loop over nonzeros */
        /* Separate the entries */
        if(el!=0) printf(", ");

        /* Get the column */
        rr = row[el]; 
        
        /* Print the nonzero */
        printf("{%d,%d}",rr,cc);
      }
    }

    printf("}\n\n");
  }

  /* Allocate work vectors */
  int n_iw=-1, n_w=-1;
  work(&n_iw, &n_w);
  const int n=256;
  if (n_iw>n || n_w>n) {
    printf("Too small work vectors allocated. Required int array of length %d"
           "(%d allocated) and double array of length %d (%d allocated)\n",
           n_iw, n, n_w, n);
    return 1;
  }
  int iw[n];
  double w[n];

  /* Evaluate the function */
  const double x_val[] = {1,2,3,4};
  const double y_val = 5;
  double res0;
  double res1[4];
  const double *all_inputs[] = {x_val,&y_val};
  double *all_outputs[] = {&res0,res1};
  eval(all_inputs,all_outputs, iw, w);
  
  printf("result (0): %g\n",res0);
  printf("result (1): [%g,%g;%g,%g]\n",res1[0],res1[1],res1[2],res1[3]);
  
  /* Free the handle */
  dlclose(handle);

  return 0;
}


// C++ (and CasADi) from here on
#include <casadi/casadi.hpp>
using namespace casadi;
using namespace std;

void usage_cplusplus(){
  cout << "---" << endl;
  cout << "Usage from C++" << endl;
  cout << endl;

  // Use CasADi's "ExternalFunction" to load the compiled function
  ExternalFunction ff("f");
  ff.init();

  // Use like any other CasADi function (note that derivatives are not supported)
  double x_val[] = {1,2,3,4};
  ff.setInput(x_val,0);
  double y_val = 5;
  ff.setInput(y_val,1);

  ff.evaluate();

  cout << "result (0): " << ff.output(0) << endl;
  cout << "result (1): " << ff.output(1) << endl;
}

int main(){
    
  // Variables
  SX x = SX::sym("x",2,2);
  SX y = SX::sym("y",1); 

  // Simple function
  vector<SX> f_in;
  f_in.push_back(x);
  f_in.push_back(y);
  vector<SX> f_out;
  f_out.push_back(sqrt(y)-1);
  f_out.push_back(sin(x)-y);
  SXFunction f(f_in,f_out);
  f.init();

  // Generate C-code
  f.generate("f");

  // Compile the C-code to a shared library
  string compile_command = "gcc -fPIC -shared -O3 f.c -o f.so";
  int flag = system(compile_command.c_str());
  casadi_assert_message(flag==0, "Compilation failed");

  // Usage from C
  flag = usage_c();
  casadi_assert_message(flag==0, "Example failed");

  // Usage from C++
  usage_cplusplus();

  return 0;
}

