/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "casadi/sx/sx_matrix.hpp"
#include "casadi/mx/mx.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/fx.hpp"
#include <ctime>

using namespace std;
using namespace CasADi;

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
vector<SX> symA(n*n);
make_symbolic(symA,"A");
for (int i=0; i<n*n; i++)
  A[i%n][i/n] = symA[i];

    SX deter;
    deter = det(n,m-1);

    double t1 = myclock();


// making function
SXFunction fcn(symA,deter);
fcn.setOption("ad_order",1);
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
    detout = fcn.output().data()[0];

    double t3 = myclock();


// evaluate 100 times
for(int k=0; k<100; ++k)
  fcn.evaluate();

  double t4 = myclock();


// evaluate 100 times
fcn.input().dataA()[0] = 1;
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
//ld " << t1-t0 << endl;
cout << "time for constructing graph " << t1-t0 << endl;
cout << "time for making function " << t2-t1 << endl;
cout << "graph+function+1 evaluation " << t3-t0 << endl;
cout << "1 evaluation " << (t4-t3)/100 << endl;
cout << "1 reverse " << (t5-t4)/100 << endl;
cout << "1 forward " << (t6-t5)/100 << endl;



#if 0

    /*--------------------------------------------------------------------------*/
    int tape_stats[STAT_SIZE];

    tapestats(tag,tape_stats);

    fprintf(stdout,"\n    independents            %d\n",tape_stats[NUM_INDEPENDENTS]);
    fprintf(stdout,"    dependents              %d\n",tape_stats[NUM_DEPENDENTS]);
    fprintf(stdout,"    operations              %d\n",tape_stats[NUM_OPERATIONS]);
    fprintf(stdout,"    operations buffer size  %d\n",tape_stats[OP_BUFFER_SIZE]);
    fprintf(stdout,"    locations buffer size   %d\n",tape_stats[LOC_BUFFER_SIZE]);
    fprintf(stdout,"    constants buffer size   %d\n",tape_stats[VAL_BUFFER_SIZE]);
    fprintf(stdout,"    maxlive                 %d\n",tape_stats[NUM_MAX_LIVES]);
    fprintf(stdout,"    valstack size           %d\n\n",tape_stats[TAY_STACK_SIZE]);

    /*--------------------------------------------------------------------------*/
    int itu = 8-n;
    itu = itu*itu*itu*itu;
    itu = itu > 0 ? itu : 1;
    double raus;

    /*--------------------------------------------------------------------------*/
    double t10 = myclock();                             /* 1. time (original) */
    for (it = 0; it < itu; it++)
        raus = pdet(n,m-1);
    double t11 = myclock();
    double rtu = itu/(t11-t10);

    double* B = new double[n2];
    double* detaut = new double[1];

    /*--------------------------------------------------------------------------*/
    double t40 = myclock();                      /* 4. time (forward no keep) */
    for (it = 0; it < itu; it++)
        forward(tag,1,n2,0,a,detaut);
    double t41 = myclock();

    /*--------------------------------------------------------------------------*/
    double t20 = myclock();                         /* 2. time (forward+keep) */
    for(it = 0; it < itu; it++)
        forward(tag,1,n2,1,a,detaut);
    double t21 = myclock();
    // fprintf(stdout,"\n %f =? %f should be the same \n",detout,*detaut);

    double u[1];
    u[0] = 1.0;

    /*--------------------------------------------------------------------------*/
    double t30 = myclock();                              /* 3. time (reverse) */
    for (it = 0; it < itu; it++)
        reverse(tag,1,n2,0,u,B);
    double t31 = myclock();

    /*--------------------------------------------------------------------------*/
    /* output of results */
    // optional generation of tape_doc.tex
    // tape_doc(tag,1,n2,a,detaut);
    fprintf(stdout,"\n first base? :   \n");
    for (i=0; i<n; i++) {
        adouble sum = 0;
        adouble* pt;
        pt = A[i];
        for (j=0; j<n; j++)
            sum += (*pt++)*B[j];
        fprintf(stdout,"%E ",sum.value());
    }
    fprintf(stdout,"\n\n times for ");
    fprintf(stdout,"\n tracing          : \t%E",(t01-t00)*rtu);
    fprintf(stdout," units \t%E    seconds",(t01-t00));
    fprintf(stdout,"\n forward (no keep): \t%E",(t41-t40)*rtu/itu);
    fprintf(stdout," units \t%E    seconds",(t41-t40)/itu);
    fprintf(stdout,"\n forward + keep   : \t%E",(t21-t20)*rtu/itu);
    fprintf(stdout," units \t%E    seconds",(t21-t20)/itu);
    fprintf(stdout,"\n reverse          : \t%E",(t31-t30)*rtu/itu);
    fprintf(stdout," units \t%E    seconds\n",(t31-t30)/itu);

#endif
    return 1;
}


























#if 0
int n;
SX **A;                        // A is an n x n matrix
SX zero = 0;

SX det(int k, int m)           // k <= n is the order of the submatrix
{ if (m == 0)                       // its column indices
        return 1.0;
    else                              // are encoded in m
    {
        SX *pt = A[k-1];
        SX   t = zero;
        int p = 1;
        int s;
        if (k%2)
            s = 1;
        else
            s = -1;
        for(int i=0; i<n; i++) {
            int p1 = 2*p;
            if (m%p1 >= p) {
                if (m == p) {
                    if (s>0)
                        t += *pt;
                    else
                        t -= *pt;
                } else {
                    if (s>0)
                        t += *pt*det(k-1, m-p); // recursive call to det
                    else
                        t -= *pt*det(k-1, m-p); // recursive call to det
                }
                s = -s;
            }
            ++pt;
            p = p1;
        }
        return t;
    }
}

int main() {
  clock_t time1 = clock();

    int i,j, m = 1;
    int tag = 1;
    int keep = 1;

    cout << "COMPUTATION OF DETERMINANTS (ADOL-C Documented Example)\n\n";
    n = 5;

//    cin >> n;

    A = new SX*[n];
    SX ad;

    SXMatrix Amat("A",n,n);
    double Aval[n][n];

    double detout = 0.0, diag = 1.0;// here keep the intermediates for
    for (i=0; i<n; i++)             // the subsequent call to reverse
    { m *= 2;
        A[i] = new SX[n];
        for (j=0; j<n; j++){
            Aval[i][j] = j/(1.0+i);      // make all elements of A independent
	    A[i][j] = Amat(i,j);
	}
        diag += Aval[i][i];       // value(adouble) converts to double
        A[i][i] += 1.0;
    }
    ad = det(n,m-1);                // actual function call.

   clock_t time2 = clock();
   cout << "graph built after " << (time2-time1)/double(CLOCKS_PER_SEC) << " seconds." << endl;

  // Create a function
  SXFunction ffcn(Amat,ad);
  ffcn->initAD(1);

   clock_t time3 = clock();
   cout << "graph sorted after " << (time3-time1)/double(CLOCKS_PER_SEC) << " seconds." << endl;

  // Give a value
  vector<double>& ip = ffcn->getInput();
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      ip[j+i*n] = Aval[i][j];

  // Evaluate
  for(int i=0; i<100; ++i){
    ffcn->evaluate();
  }

   clock_t time4 = clock();
   cout << "function evaluated 100 times after " << (time4-time1)/double(CLOCKS_PER_SEC) << " seconds." << endl;

return 0;


  // Get output
  detout = ffcn->getOutput()[0];

    printf("\n %f - %f = %f  (should be 0)\n",detout,diag,detout-diag);

    double u[1];
    u[0] = 1.0;
    double* B = new double[n*n];

    // set a backward seed
    ffcn->getOutputSeed()[0] = u[0];

    // Evaluate backwards
    ffcn->evaluateAdj();

    // Get the result
    vector<double>& ider = ffcn->getInputSeed();
    for(int i=0; i<n*n; ++i)
      B[i] = ider[i];

    cout << " \n first base? : ";
    for (i=0; i<n; i++) {
        double sum = 0;
        for (j=0; j<n; j++)             // the matrix A times the first n
            sum += Aval[i][j]*B[j];          // components of the gradient B
        cout << sum << " ";      // must be a Cartesian basis vector
    }
    cout << "\n";

    return 1;

}

#endif
