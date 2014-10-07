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


#include "core/mx/mx.hpp"
#include "core/std_vector_tools.hpp"
#include "core/function/sx_function.hpp"
#include "core/sx/sx_tools.hpp"
#include "core/function/function.hpp"
#include <ctime>

using namespace std;
using namespace casadi;

/** Determinant example from ADOL-C */

#include <cstdio>
#include <iostream>

/*----------------------------------------------------------------------------
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     detexam.cpp
 Revision: $Id: detexam.cpp 91 2010-02-24 07:56:58Z awalther $
 Contents: modified computation of determinants

 Copyright (c) Andrea Walther, Andreas Griewank, Andreas Kowarz, 
               Hristo Mitev, Sebastian Schlenkrich, Jean Utke, Olaf Vogel 
  
 This file is part of ADOL-C. This software is provided as open source.
 Any use, reproduction, or distribution of the software constitutes 
 recipient's acceptance of the terms of the accompanying license file.
 
---------------------------------------------------------------------------*/

/****************************************************************************/
/*                                                                 INCLUDES */
//#include <adolc.h>
//#include <../examples/additional_examples/clock/myclock.h>


double myclock(){
  return clock() / double(CLOCKS_PER_SEC);
}


/****************************************************************************/
/*                                                           DOUBLE ROUTINE */
int n,it;
double** PA;
double pdet( int k, int m ) {
    if (m == 0)
        return 1.0 ;
    else {
        double* pt = PA[k-1];
        double t = 0;
        int p = 1;
        int s;
        if (k%2)
            s = 1;
        else
            s = -1;
        for (int i=0; i<n; i++) {
            int p1 = 2*p;
            if (m%p1 >= p) {
                if (m == p) {
                    if (s>0)
                        t += *pt;
                    else
                        t -= *pt;
                } else {
                    if (s>0)
                        t += *pt*pdet(k-1, m-p);
                    else
                        t -= *pt*pdet(k-1, m-p);
                }
                s = -s;
            }
            ++pt;
            p = p1;
        }
        return t;
    }
}

/****************************************************************************/
/*                                                          SX ROUTINE */
SX** A;
SX zero = 0;
SX det( int k, int m ) {
    if (m == 0)
        return 1.0;
    else {
        SX* pt = A[k-1];
        SX t = zero;
        int p = 1;
        int s;
        if (k%2)
            s = 1;
        else
            s = -1;
        for (int i=0; i<n; i++) {
            int p1 = 2*p;
            if (m%p1 >= p) {
                if (m == p) {
                    if (s>0)
                        t += *pt;
                    else
                        t -= *pt;
                } else {
                    if (s>0)
                        t += *pt*det(k-1, m-p);
                    else
                        t -= *pt*det(k-1, m-p);
                }
                s = -s;
            }
            ++pt;
            p = p1;
        }
        return t;
    }
}

/****************************************************************************/
/*                                                             MAIN PROGRAM */
int main() {
    int i, j;
    int tag = 1;
    fprintf(stdout,"COMPUTATION OF DETERMINANTS Type 1 (ADOL-C Example)\n\n");
    fprintf(stdout,"order of matrix = ? \n");
    scanf("%d",&n);
    A  = new SX*[n];
    PA = new double*[n];
    int n2 = n*n;
    double* a = new double[n2];

    /*--------------------------------------------------------------------------*/
    /* Preparation */
    double diag = 0;
    int m = 1;
    double* pa = a;
    for (i=0; i<n; i++) {
        m *= 2;
        PA[i] = new double[n];
        double* ppt = PA[i];
        for (j=0; j<n; j++) {
            *ppt++ = j/(1.0+i);
            *pa++  = j/(1.0+i);
        }
        diag += PA[i][i];   // val corrected to value 2/23/91
        PA[i][i] += 1.0;
        a[i*n+i] += 1.0;
    }
    diag += 1;

    /*--------------------------------------------------------------------------*/
    double t0 = myclock();                               /* building graph */

    for (i=0; i<n; i++) {
        A[i] = new SX[n];
        SX* pt = A[i];
    }
vector<SX> symA = ssym("A",n*n).data();
for (int i=0; i<n*n; i++)
  A[i%n][i/n] = symA[i];

    SX deter;
    deter = det(n,m-1);

    double t1 = myclock();


// making function
SXFunction fcn(symA,deter);
fcn.init();
vector<double>& inp = fcn.input().data();
    for (i=0; i<n; i++) {
        double* ppt = PA[i];
        for (j=0; j<n; j++)
            inp[i+n*j] = *ppt++;
    }

    double t2 = myclock();

// evaluate
fcn.evaluate();


// get result
    double detout = 0.0;
    detout = fcn.output().front();

    double t3 = myclock();


// evaluate 100 times
for(int k=0; k<100; ++k)
  fcn.evaluate();

  double t4 = myclock();


// evaluate 100 times
fcn.adjSeed()[0] = 1;
for(int k=0; k<100; ++k)
  fcn.evaluate(0,1);

  double t5 = myclock();

// evaluate 100 times
for(int k=0; k<100; ++k)
  fcn.evaluate(1,0);

  double t6 = myclock();

//    deter >>= detout;
//    trace_off();
    fprintf(stdout,"\n %f =? %f should be the same \n",detout,diag);

if (abs(detout-diag)>1e-8)
  return -1;    

//ld " << t1-t0 << endl;
cout << "time for constructing graph " << t1-t0 << endl;
cout << "time for making function " << t2-t1 << endl;
cout << "graph+function+1 evaluation " << t3-t0 << endl;
cout << "1 evaluation " << (t4-t3)/100 << endl;
cout << "1 reverse " << (t5-t4)/100 << endl;
cout << "1 forward " << (t6-t5)/100 << endl;

    return 0;
}
