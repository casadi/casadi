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
/*#define CASADI_NDEBUG*/
#include <casadi/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>


#include "nlp.h"
#undef f_OPNUM
#include "opcode.hd"
#include "asl_reader.hpp"
#define asl ((ASL_fg*)cur_ASL)
#include "assert.h"


using namespace CasADi;

// Variables
std::vector<SX> vars;
std::vector<SX> objs;
std::vector<SX> cons;
std::vector<SX> intermediates;



extern "C"{

typedef struct maxinfo {
  int ndv;
  int nvt;
} maxinfo;

maxinfo bmax, cmax, omax;

list **com_ifset;
int ndv, ndvmax, npd, nv1, nv2, nvt, nvtmax;
int *cvmap;
char *intsk;
extern efunc *r_op[];
#define f_OPNUM r_op[79]
#define f_OPHOL r_op[80]
#define f_OPVARVAL      r_op[81]

expr_nx *nums;

int *dvtfree, *rownos;
int ndvfree, ndvtmax, nvtfree;

int new_vt(){
  int ret = nvtfree < ndvtmax ? dvtfree[nvtfree++] : nvt++;
  if(ret>=intermediates.size()){
    intermediates.resize(ret+1);
  }
  return ret;
}

#define NDVTGULP 1000
void vt_free(int k){
  if (ndvfree >= nvtfree){
    int i = ndvtmax;
    int j = ndvtmax += NDVTGULP;
    dvtfree = (int *)Realloc(dvtfree, j*sizeof(int));
    while(i > nvtfree)
      dvtfree[--j] = dvtfree[--i];
    nvtfree += NDVTGULP;
  }
  dvtfree[--nvtfree] = k;
}
void Lset(linpart *L, int nlin){
  linpart *Le;
  double t;
  register dLR *d;
  
  for(Le = L + nlin; L < Le; L++) {
    t = L->fac;
    d = (dLR *)(&L->fac);
    if(t == 1.){
      d->kind = dLR_one;
    } else if (t == -1.) {
      d->kind = dLR_negone;
    } else {
      d->kind = dLR_VP;
      *(d->o.vp = (real *)mem(sizeof(real))) = t;
    }
  }
}

char *f_OPVARVAL1(expr *e, char *buf){
  int k = e->a;
  casadi_assert(k>=0);
  casadi_assert(k<nv1);
  return buf;
}

char *f_OPNUM1(register expr *e, char *rv){
  return rv;
}
SX get_expression(expr *e){
  int k = (int)(e->op);
  e->op = r_op[k];
  int k1 = op_type[k];
  double t;
  int i;
/*  if(e->a>0){
    return intermediates.at(e->a);
  } else {*/
    switch(k1){
      case 9: // number
        t = ((expr_nx *)e)->v;
        return t;
      case 10: // variable value
        i = (expr_v *)e - var_e;
        if(i<vars.size()){
          return vars.at(i);
        } else {
          casadi_assert_message(0,"Dependent variables not supported");
 
          // The dependent expressions are available:
          int i = comb;
          int j = combc;
          cexp *c;
          for(c = cexps + i; i < j; i++, c++){
            linpart *L = c->L;
            expr *e = c->e;
            int nlin = c->nlin;
          }
        }
    }
/*  }*/
  casadi_assert(0);
}

int ewalk(expr *e){
  int operation = (int)(e->op);
  e->op = r_op[operation];
  int k1 = op_type[operation];

  int i, j, k;
  expr *e1, **ep, **epe;
  expr_if *eif;
  expr_va *eva;
  expr_f *ef;
  argpair *ap, *ape;
  de *d;
  derp *dp;
  dLR *LR, **LRp;
  efunc *op;
  double t;
  SX x,y,f,cx,cy;
  switch(k1){
    case 11: /* OP2POW */
      e->dL = 0;
      /* no break */
      case 1: /* unary */
        j = ewalk(e->L.e);
        i = new_vt();
        
        if (j > 0){
          vt_free(j);
        }
        e->a = i;
        
        if(j==0){
          x = get_expression(e->L.e);
        } else {
          x = intermediates.at(j);
        }
        switch(operation){
          case FLOOR: f = floor(x); break;
          case CEIL: f = ceil(x); break;
          case ABS: f = abs(x); break;
          case OPUMINUS: f = -x; break;
          case OPNOT: f = !x; break;
          case OP_tanh: f = tanh(x); break;
          case OP_tan: f = tan(x); break;
          case OP_sqrt: f = sqrt(x); break;
          case OP_sinh: f = sinh(x); break;
          case OP_sin: f = sin(x); break;
/*          case OP_log10: f = log10(x); break;*/
          case OP_log: f = log(x); break;
          case OP_exp: f = exp(x); break;
          case OP_cosh: f = cosh(x); break;
          case OP_cos: f = cos(x); break;
/*          case OP_atanh: f = atanh(x); break;*/
          case OP_atan: f = atan(x); break;
/*          case OP_asinh: f = asinh(x); break;*/
          case OP_asin: f = asin(x); break;
/*          case OP_acosh: f = acosh(x); break;*/
          case OP_acos: f = acos(x); break;
          case OPCOUNT: casadi_assert_message(0,"Unary operation OPCOUNT not implemented. What is it supposed to do?"); break;
          case OPNUMBEROF: casadi_assert_message(0,"Unary operation OPNUMBEROF not implemented. What is it supposed to do?"); break;
          case OPNUMBEROFs: casadi_assert_message(0,"Unary operation OPNUMBEROFs not implemented. What is it supposed to do?"); break;
          case OPALLDIFF: casadi_assert_message(0,"Unary operation OPALLDIFF not implemented. What is it supposed to do?"); break;
          case OP2POW: f = x*x; break;
          default:
            casadi_assert_message(0,"unknown: " << operation);
        };
        intermediates.at(i) = f;
        
        return i;
        
      case 2: /* binary */
        k1 = 1;
/*        casadi_assert(0);*/
        switch(operation) { // WHY THIS????
          case OPPLUS:
          case OPMINUS:
          case OPREM:
            e->dL = 1.;
            break;
          case OPMULT:
            k1 = 3;
            /* no break */
            default:
              e->dL = 0.;
        }
        k = ewalk(e->L.e);
        j = ewalk(e->R.e);
        i = new_vt();
        if (j > 0){
          vt_free(j);
        }
        if (k > 0){
          vt_free(k);
        }
        e->a = i;
        if(k==0){
          x = get_expression(e->L.e);
        } else {
          x = intermediates.at(k);
        }
        
        if(j==0){
          y = get_expression(e->R.e);
        } else {
          y = intermediates.at(j);
        }
        
        cx = e->dL;
        if(e->dL==0) cx = 1; // NOTE: Sometimes 0, why???
        cy = e->dR;
        if(e->dR==0) cy = 1; // NOTE: Sometimes 0, why???
        
        switch(operation){
          case OPPLUS: f = cx*x + cy*y; break;
/*          case OPMINUS: f = cx*x - cy*y; break;*/
          case OPMINUS: f = cx*x + cy*y; break; // NOTE: changed - -> +
          case OPMULT: f = (cx*x) * (cy*y); break;
          case OPDIV: f = (cx*x) / (cy*y); break;
/*          case OPREM: f = cx*x ? cy*y; break;*/
          case OPPOW: f = pow(cx*x,cy*y); break;
/*          case OPLESS: f = cx*x ? cy*y; break;*/
          case OPOR: f = (cx*x) || (cy*y); break;
          case OPAND: f = (cx*x) && (cy*y); break;
          case LT: f = (cx*x) < (cy*y); break;
          case LE: f = (cx*x) <= (cy*y); break;
          case EQ: f = (cx*x) == (cy*y); break;
          case GE: f = (cx*x) >= (cy*y); break;
          case GT: f = (cx*x) > (cy*y); break;
          case NE: f = (cx*x) != (cy*y); break;
/*          case OP_atan2: f = cx*x ? cy*y; break;*/
/*          case OPintDIV: f = cx*x ? cy*y; break;*/
/*          case OPprecision: f = cx*x ? cy*y; break;*/
/*          case OPround: f = cx*x ? cy*y; break;*/
/*          case OPtrunc: f = cx*x ? cy*y; break;*/
/*          case OPATLEAST: f = cx*x ? cy*y; break;
          case OPATMOST: f = cx*x ? cy*y; break;
          case OPEXACTLY: f = cx*x ? cy*y; break;
          case OPNOTATLEAST: f = cx*x ? cy*y; break;
          case OPNOTATMOST: f = cx*x ? cy*y; break;
          case OPNOTEXACTLY: f = cx*x ? cy*y; break;*/
/*          case OP_IFF: f = cx*x ? cy*y; break;*/
          case OP1POW: f = pow(cx*x,cy*y); break;
          case OPCPOW: f = pow(cx*x,cy*y); break;
          default:
            casadi_assert_message(0,"unknown: " << operation);
        };
        
        if(false){
          using namespace std;
          cout << "operation = " << operation << endl;
          cout << "x = " << x << endl;
          cout << "y = " << y << endl;
          cout << "cx = " << cx << endl;
          cout << "cy = " << cy << endl;
          cout << "f = " << f << endl;
        }
        
        
        intermediates.at(i) = f;
        return i;
        
      case 3: /* vararg (min, max) */
        casadi_assert(0);
        
      case 4: /* piece-wise linear */
        casadi_assert(0);
        
      case 5: /* if */
        casadi_assert(0);
        
      case 6: /* sumlist */
        ep = e->L.ep;
        epe = e->R.ep;
        i = ewalk(e1 = *ep++);
        if(i==0){
          x = get_expression(e1);
        } else {
          x = intermediates.at(i);
        }
        f = x;
        
        j = ewalk(e1 = *ep++);
        if(j==0){
          x = get_expression(e1);
        } else {
          x = intermediates.at(j);
        }
        f += x;
        if (i > 0) {
          if (j > 0)
            vt_free(j);
        } else if (!(i = j)){
          i = new_vt();
        } 
        do {
          if ((j = ewalk(e1 = *ep++)) > 0){
            vt_free(j);
          }
          if(j==0){
            x = get_expression(e1);
          } else {
            x = intermediates.at(j);
          }
          f += x;
        } while(ep < epe);
        
        e->a = i;
        intermediates.at(i) = f;
        
        return i;
        
      case 7: /* function call */
        ef = (expr_f *)e;
        i = k = 0;
        casadi_assert(0);
        for(ap = ef->ap, ape = ef->ape; ap < ape; ap++){
          if (j = ewalk(ap->e)){
            if (j < 0) {
              i = j;
            } else {
              k++;
            }
          }
        }
        if(k){
          for(ap = ef->ap; ap < ape; ap++) {
            op = ap->e->op;
            if(op != f_OPNUM && op != f_OPHOL && op != f_OPVARVAL && (j = ap->e->a) > 0){
              vt_free(j);
            }
          }
        }
        if(!i){
          return e->a = new_vt();
        }
        return e->a = i;
              
      case 8: /* Hollerith */
        break;
        
      case 9: /* number */
        t = ((expr_n *)e)->v;
        ((expr_nx *)e)->v = t;
        ((expr_nx *)e)->next = nums;
        nums = (expr_nx *)e;
        break;
        
      case 10: /* variable value */
        i = (expr_v *)e - var_e;
        
        if(intermediates.empty()){
          intermediates.resize(1);
        }
/*        intermediates.at(0) = vars.at(i);*/
        
        
        break;
        /*DEBUG*/default:
        /*DEBUG*/ fprintf(Stderr, "bad opnumber %d in ewalk\n", operation);
        /*DEBUG*/ exit(1);
  }
  return 0;
}

void max_save(maxinfo *m){
  if (nvtmax < nvt)
    nvtmax = nvt;
  m->nvt = nvtmax - 1;
  m->ndv = ndvmax;
}

void zset(register int *z){
  register int *ze;
  if(z){
    for(ze = z + 2**z; z < ze; z += 2){
      var_e[z[1]].a = z[2];
    }
  }
}

void ndvtreset(int *z){
  nvt = 1;
  nvtfree = ndvtmax;
  ndv = 0;
  ndvfree = 0;
  zset(z);
}

void comwalk(int i, int n){
  if (i >= n) return;
  zset(zaC[i]);
  list **ifset = com_ifset + i;
  cexp *c = cexps + i;
  cexp *ce;
  for(ce = cexps + n; c < ce; c++) {
    Lset(c->L, c->nlin);
    ndvtreset(zaC[i]);
    nvtfree = ndvtmax;
    ewalk(c->e);
    ifset++;
    if (nvtmax < nvt){
      nvtmax = nvt;
    }
  }
}

void cde_walk(cde *d, int n, maxinfo *m, int *c1st, int **z, bool obj){
  int ii=0;
  
  int j = *c1st++;
  cde *De = d + n;
  for(; De > d; d++, ii++){
    int i = j;
    int k = *c1st++;
    ndvtreset(*z++);
    while(j < k){
      cexp1 *c1 = cexps1 + j++;
      Lset(c1->L, c1->nlin);
      ewalk(c1->e);
    }
    nvtfree = ndvtmax;
    int kk = ewalk(d->e);
    
    // What are we calculating?
    SX& v = obj ? objs.at(ii) : cons.at(ii);
    if(kk==0){
      v = get_expression(d->e);
      
    } else {
      v = intermediates.at(kk);
    }
    
    while(i < k){
      int i1 = i++ + ncom0;
      if (!cvmap[i1]){
        cvmap[i1] = nvt++;
      }
    }
    if (nvtmax < nvt){
      nvtmax = nvt;
    }
  }
  max_save(m);
}

void output(){
  // Objective
  static char rv[] = "rv";
  casadi_assert(n_obj==1);
  ograd *og;
  for(og = Ograd[0]; og; og = og->next){
    if (og->coef){
      break;
    }
  }
  
  char *eval;
  int *c1 = o_cexp1st;
  int i=0;
  for(; i < n_obj; i++, c1++){
    ograd *og = Ograd[i];
    for(; og; og = og->next){
      if(og->coef){
        break;
      }
    }
    if (og) {
      SX f = og->coef * vars.at(og->varno) + objs.at(0);
      while(og = og->next){
        if (og->coef){
          f += og->coef * vars.at(og->varno);
        }
      }
      objs.at(0) = f;
    }
  }
  
  // Constraints
  if(n_con){
    int i;
    int *c1 = c_cexp1st;
    for(i = 0; i < n_con; i++, c1++) {
      char *s = const_cast<char*>("aa");
      
      // Add linear ?
      cgrad *cg;
      for(cg = Cgrad[i]; cg; cg = cg->next){
        if (cg->coef) {
          SX f = cg->coef * vars.at(cg->varno);
          while(cg = cg->next){
            if (cg->coef){
              f += cg->coef * vars.at(cg->varno);
            }
          }
          cons.at(i) += f;
          break;
        }
      }
    }
  }
}

void get_rownos(){
  int i = n_con, i1, j, j1;
  cgrad *cg;
  memset((char *)rownos, 0, nzc*sizeof(int));
  while(i1 = i){
    for(cg = Cgrad[--i]; cg; cg = cg->next){
      rownos[cg->goff] = i1;
    }
  }
  
  for(i = 0; i <= n_var; i++){
    A_colstarts[i]++;
  }
  
  if(intsk){
    i1 = j = 0;
    --intsk;
    for(i = 1; i <= n_var; i++){
      j1 = A_colstarts[i] - 1;
      if (!intsk[i]){
        for(; j < j1; j++){
          rownos[i1++] = rownos[j];
        }
      }
      
      A_colstarts[i] = i1 + 1;
      j = j1;
    }
    
    ++intsk;
    nzc = i1;
    if(nlvbi && (nlvci < nlvc || nlvoi < nlvo) || nlvci && nlvoi < nlvo){
      for(i = 0; i < n_var; i++){
        A_colstarts[i] = A_colstarts[i+1];
      }
      i = n_con;
      while(--i >= 0){
        for(cg = Cgrad[i]; cg; cg = cg->next){
          if (!intsk[cg->varno]){
            cg->goff = --A_colstarts[cg->varno] - 1;
          }
        }
      }
    }
  }
}
    
int main(int argc, char **argv){
  ASL_alloc(ASL_read_fg);
  g_fmt_decpt = 1;
  want_derivs = 0;
  return_nofile = 1;
  casadi_assert(argc==2);
  char *progname = argv[1];
  
  fint L = strlen(progname);
  FILE *nl = jacdim0(progname, L);

  int i;
  for(i = 0; i < N_OPS; i++){
    r_ops[i] = (efunc *)i;
  }
  
  nv1 = c_vars > o_vars ? c_vars : o_vars;
  int ncom = (i = comb + comc + como) + comc1 + como1;
  nv2 = nv1 + ncom;
        
  c_cexp1st = (int *)Malloc((n_con + n_obj + 2)*sizeof(int));
  o_cexp1st = c_cexp1st + n_con + 1;
  zac = (int **)Malloc((n_con + n_obj + i)*sizeof(int*));
  zao = zac + n_con;
  zaC = zao + n_obj;
  
  if (n_con){
    rownos = (int *)Malloc((nzc + nv1 + 1)*sizeof(int));
    A_colstarts = rownos + nzc;
  }
    
  LUv = (real *)Malloc((3*nv1+2*n_con)*sizeof(real));
  LUrhs = LUv + 2*nv1;
  X0 = LUrhs + 2*n_con;
      
  size_expr_n = sizeof(expr_nx);
  fg_read(nl,0);
        
  list **c_ifset = (list **)Malloc((n_con + n_obj + ncom0)*sizeof(list *));
  list **o_ifset = c_ifset + n_con;
  com_ifset = o_ifset + n_obj;
  op_type[OP1POW] = 2;
  op_type[OP2POW] = 11;
  op_type[OPCPOW] = 2; // NOTE: Joel: I added this
  
  dvtfree = (int *)Malloc(NDVTGULP*sizeof(int));
  ndvtmax = nvtfree = NDVTGULP;
  
  if (n_con) get_rownos();
  cvmap = (int *)Malloc(ncom*sizeof(int));
  npd = 0;
  for(i = 0; i < ncom0; i++){
    cvmap[i] = -(++npd);
  }
  if (ncom1){
    memset((char *)&cvmap[ncom0], 0, ncom1*sizeof(int));
  }
  
  ndv = 0;
  
  memset((char *)adjoints, 0, amax*sizeof(real));
  
  
  // Allocate CasADi variables
  vars = ssym("x",n_var).data();
  
  // Allocate CasADi objective functions
  objs.resize(n_obj);
  
  // Allocate CasADi constraint functions
  cons.resize(n_con);
  
  
  
  
  comwalk(0,comb);
  max_save(&bmax);
  comwalk(comb, combc);
  cde_walk(con_de, n_con, &cmax, c_cexp1st, zac,false);
  nvtmax = bmax.nvt + 1;
  ndvmax = bmax.ndv;
  comwalk(combc, ncom0);
  cde_walk(obj_de, n_obj, &omax, o_cexp1st, zao,true);

  int nv = nv1 + ncom;
  for(i = 0; i < nv; i++){
    var_e[i].op = (efunc *)f_OPVARVAL1;
  }
  
  for(expr_nx *enx = nums; enx; enx = enx->next){
    enx->op = f_OPNUM1;
  }
  
  output();

  if(false){
    std::cout << "objs = " << std::endl;
    for(std::vector<SX>::const_iterator it=objs.begin(); it!=objs.end(); ++it){
      std::cout << *it << std::endl;
    }
    
    std::cout << "cons = " << std::endl;
    for(std::vector<SX>::const_iterator it=cons.begin(); it!=cons.end(); ++it){
      std::cout << *it << std::endl;
    }
  }
  
  using namespace std;
  cout << "solving problem" << endl;

  // Get the variable bounds
  vector<double> x_lb(nv1), x_ub(nv1);
  for(int i=0; i<nv1; ++i){
    x_lb[i] = LUv[2*i];
    x_ub[i] = LUv[2*i+1];
  }
  
  // Get the constraint bounds
  vector<double> g_lb(n_con), g_ub(n_con);
  for(int i=0; i<n_con; ++i){
    g_lb[i] = LUrhs[2*i];
    g_ub[i] = LUrhs[2*i+1];
  }

  // Get the starting guess
  vector<double> x_guess(X0,X0+nv1);
    
  // NLP expressions
  SXMatrix x = vars;
  SXMatrix f = objs;
  SXMatrix g = cons;
  
  // NLP functions
  SXFunction ffcn(x,f);
  SXFunction gfcn(x,g);
  
  // NLP solver
  IpoptSolver nlp_solver(ffcn,gfcn);
  
  // Set options
/*  nlp_solver.setOption("max_iter",10);*/
/*  nlp_solver.setOption("verbose",true);*/
//  nlp_solver.setOption("linear_solver","ma57");
  nlp_solver.setOption("generate_hessian",true);
/*  nlp_solver.setOption("hessian_approximation","limited-memory");*/
  
  // Initialize NLP solver
  nlp_solver.init();
  
  // Pass the bounds and initial guess
  nlp_solver.setInput(x_lb,NLP_LBX);
  nlp_solver.setInput(x_ub,NLP_UBX);
  nlp_solver.setInput(g_lb,NLP_LBG);
  nlp_solver.setInput(g_ub,NLP_UBG);
  nlp_solver.setInput(x_guess,NLP_X_INIT);
  
  // Solve NLP
  nlp_solver.solve();
  
  return 0;
}
    
}
