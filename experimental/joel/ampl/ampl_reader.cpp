/****************************************************************
Copyright (C) 1997, 2001 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/
/*#define CASADI_NDEBUG*/
#include <casadi/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>


#include "nlp.h"
#undef f_OPNUM
#include "opcode.hd"
#include "ampl_reader.hpp"
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
  int ncond;
  int ndv;
  int nvt;
  int needT1;
  /* !! plus goo for function evaluations */
} maxinfo;

maxinfo bmax, cmax, omax;

list **c_ifset, **com_ifset, **o_ifset;
int branches = 0;
int condlevel, fwalk, iflevel, maxa, ncond, ndv, ndvmax, npd, nplterm, nv1, nv2, nvt, nvtmax;
int *cvmap;
fint nJ;
expr *Plterms;
char *intsk;
extern efunc *r_op[];
#define f_OPNUM r_op[79]
#define f_OPHOL r_op[80]
#define f_OPVARVAL      r_op[81]

extern char *progname;
char *T = "t1";
expr_nx *nums;

int *dvtfree, *rownos;
int ndvfree, ndvtmax, needT1, nvtfree;

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

char *fpval(real r){
  static char buf[32];
  if(r >= Infinity){
    return "1.7e308";
  } else if (r <= negInfinity){
    return "-1.7e308";
  }
  g_fmt(buf, r);
  return buf;
}


void binop(char *a, char *b, char *op, char *c){
  if(a){
    if(a == b){
      printf("\t%s %s= %s;\n", a, op, c);
    } else {
      printf("\t%s = %s %s %s;\n", a, b, op, c);
    }
  }
}

void Lset(linpart *L, int nlin){
  linpart *Le;
  double t;
  register dLR *d;
  
  for(Le = L + nlin; L < Le; L++) {
    t = L->fac;
    d = Make_dLR(&L->fac);
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
  sprintf(buf, "x[%d]", k);
  return buf;
}

// char *f_OPVARVAL1(expr *e, char *buf){
//   int k = e->a;
//   casadi_assert(k>=0);
//   if(k >= nv1) {
//     k = cvmap[(expr_v *)e - var_e - nv1];
//     k += 1;
//   }
//   sprintf(buf, "x[%d]", k);
//   return buf;
// }





char *f_OPNUM1(register expr *e, char *rv){
  g_fmt(rv,((expr_nx *)e)->v);
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
        return vars.at(std::min<int>(i,vars.size()-1));
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
        
/*        printf("variable value = %d, vars.size() = %d, intermediates.empty() = %d\n",i,vars.size(),intermediates.empty());*/
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
  m->ncond = ncond;
  m->needT1 = needT1;
}

void max_restore(maxinfo *m){
  nvtmax = m->nvt + 1;
  ndvmax = m->ndv;
  ncond = m->ncond;
  needT1 = m->needT1;
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

static char * vprod(real t, int k){
  static char buf[64];
  int i;
  
  if (t == 1.){
    sprintf(buf, "x[%d]", k);
  } else if (t == -1.) {
    buf[0] = '-';
    sprintf(buf+1, "x[%d]", k);
  } else {
    i = g_fmt(buf, t);
    buf[i++] = '*';
    sprintf(buf+i, "x[%d]", k);
  }
  return buf;
}

char *rv_output(char *rv, char *eval, ograd *og){
  char *s;
  for(; og; og = og->next){
    if(og->coef){
      break;
    }
  }
  if (og) {
    
    SX f = og->coef * vars.at(og->varno);
    
    s = vprod(og->coef, og->varno);
    if (strcmp(eval, "0.")){
      binop(rv, eval, "+", s);
      f += objs.at(0);
    } else {
      printf("\t%s = %s;\n", rv, s);
    }
    while(og = og->next){
      if (og->coef){
        binop(rv, rv, "+",vprod(og->coef, og->varno));
        f += og->coef * vars.at(og->varno);
      }
    }
    objs.at(0) = f;
    return rv;
  }
  return eval;
}

char *con_linadd(int i, char *s){
  cgrad *cg;
  char *s1;
  
  for(cg = Cgrad[i]; cg; cg = cg->next){
    if (cg->coef) {
      SX f = cg->coef * vars.at(cg->varno);
      s1 = vprod(cg->coef, cg->varno);
      if (strcmp(s,"0.")){
        binop(T, s, "+", s1);
      } else {
        printf("\t%s = %s;\n", T, s1);
      }
      s = T;
      while(cg = cg->next){
        if (cg->coef){
          f += cg->coef * vars.at(cg->varno);
          binop(s, s, "+",vprod(cg->coef, cg->varno));
        }
      }
      cons.at(i) += f;
      break;
    }
  }
  return s;
}

int nzcgrad(){
  cgrad *cg;
  int i;
  for(i = 0; i < n_con; i++)
    for(cg = Cgrad[i]; cg; cg = cg->next)
      if (cg->coef)
        return 1;
      return 0;
}

char *cv_name(linpart *L, char *buf){
  expr_v ev;
  expr *ep = (expr *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
  return f_OPVARVAL1(ep, buf);
}

void com_out(expr *e, linpart *L, int nlin, int k0){
  char buf[32], bufg[32], res[32], vn[32];
  printf("\n%s\t/*** defined variable %d ***/\n\n", "", k0+1);
  int j = cvmap[k0];
  casadi_assert(j>=0);
  efuncb *eb = (efuncb *)e->op;
  if (eb != (efuncb *)OPNUM && eb != (efuncb *)OPVARVAL){
    e->a = j;
  }
  j--;
  char *s = "v[%d]";
  sprintf(res, s, j);
  s = (*(efuncb *)e->op)(e, buf);
  if(!L){
    if(strcmp(res,s)){
      printf("\t%s = %s;\n", res, s);
    }
    return;
  }
  int asg = !strcmp(s, "0.");
  
  dLR *d;
  int op;
  double t;
  int k;
  linpart *Le;
  for(Le = L + nlin; L < Le; L++, asg = 0) {
    d = dLRp(L->fac);
    op = '+';
    switch(k = d->kind){
      case dLR_negone:
        op = '-';
      case dLR_one:
        break;
      case dLR_VP:
        t = *d->o.vp;
        if (t < 0. && !asg) {
          t = -t;
          op = '-';
        }
        g_fmt(bufg, t);
        break;
      default:
        /*DEBUG*/ 
        fprintf(Stderr,"Bad d->kind = %d in com_walk\n", d->kind);
        exit(14);
    }
    if(asg){
      printf(op == '-' ? "\t%s = -" : "\t%s = ", res);
    } else if (res != s) {
      printf("\t%s = %s %c ", res, s, op);
    } else {
      printf("\t%s %c= ", res, op);
    }
    if (k == dLR_VP){
      printf("%s*", bufg);
    }
    printf("%s%s", cv_name(L,vn), ";\n");
    s = res;
  }
}

char *putout(expr *e, int i, int j, char *what, int k, int *z){
  static char buf[32];
  cexp1 *c1;
  
  ++k;
  ndvtreset(z);
  if (i < j) {
    if(what){
      printf("\n\n%s\t/*** defined variables for %s %d ***/\n","", what, k);
    }
    
    for(c1 = cexps1 + i; i < j; i++, c1++){
      com_out(c1->e, c1->L, c1->nlin, i + ncom0);
    }
  }
  printf(what ? "\n%s  /***  %s %d  ***/\n\n" : "\n%s  /***  objective ***/\n\n", "", what, k);
  return (*(efuncb *)e->op)(e, buf);
}

void obj_output(){
  static char rv[] = "rv";
  casadi_assert(n_obj==1);
  printf(" real feval0_(long *nobj, real *x){\n");
  ograd *og;
  for(og = Ograd[0]; og; og = og->next){
    if (og->coef){
      break;
    }
  }
/*  casadi_assert(og==0);*/
  printf(" /* Work vector */\n");
  printf(" real v[%d];\n", omax.nvt);
  if(omax.ncond){
    printf(" static int cond[%d];\n", omax.ncond);
  }
  if (omax.needT1){
    printf(" real t1, t2;\n");
  }
  
  char *eval, *s;
  int *c1 = o_cexp1st;
  int i=0;
  for(; i < n_obj; i++, c1++){
    eval = putout(obj_de[i].e, c1[0], c1[1],n_obj > 1 ? const_cast<char*>("objective") : NULL, i, zao[i]);
    s = rv_output(rv, eval, Ograd[i]);
    printf("\n\treturn %s;\n", s);
  }
  
  printf("}\n");
  branches = 0;
}
  
void con_output(){
  printf("\n void ceval0_(real *x, real *c)\n{");

  if(n_con){
    printf("\n\treal v[%d];\n", cmax.nvt);
    printf("\tstatic int cond[%d];\n", cmax.ncond);
    printf("\treal t1, t2;\n");
      
    int i;
    int *c1 = c_cexp1st;
    for(i = 0; i < n_con; i++, c1++) {
      char *s = putout(con_de[i].e, c1[0], c1[1], "constraint", i, zac[i]);
      printf("\tc[%d] = %s;\n", i, con_linadd(i,s));
    }
  }
  printf("}\n");
  branches = 0;
}
      
void output(){
  printf("#include \"math.h\"\n#define real double\n");
  printf(" real boundc_[1+%d+%d] /* Infinity, variable bounds, constraint bounds */ = {1.7e308",2*nv1, 2*n_con);
  real *b = LUv;
  real *be = b + 2*(n_con + nv1);
  while(b < be){
    printf(",%s", fpval(*b++));
  }
  
  printf("};\n\n");
        
  obj_output();
  con_output();
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
  {
    using namespace std;
    cout << "erfinv(0.4) = " << erfinv(0.4) << endl;
    cout << "erfinv(0.8) = " << erfinv(0.8) << endl;
    cout << "erfinv(-0.8) = " << erfinv(-0.8) << endl;
    cout << "erfinv(-1.4) = " << erfinv(-1.4) << endl;
    cout << "erfinv(2.4) = " << erfinv(2.4) << endl;
    cout << "erfinv(1.0) = " << erfinv(1.0) << endl;
    cout << "erfinv(-1.0) = " << erfinv(-1.0) << endl;
    return 0;
  }
  
  ASL_alloc(ASL_read_fg);
  g_fmt_decpt = 1;
  want_derivs = 0;
  return_nofile = 1;
  // progname = "/home/janderss/src/ampl/netlib.org/ampl/solvers/examples/cork.nl";
  // progname = "/home/janderss/Desktop/ampl_test/corkscrw_small.nl";
  // progname = "/home/janderss/Desktop/ampl_test/corkscrw.nl";
  // progname = "/home/janderss/Desktop/ampl_test/clnlbeam.nl";
  // progname = "/home/janderss/Desktop/ampl_test/clnlbeam_small.nl";
  // progname = "/home/janderss/Desktop/ampl_test/optmass.nl";
  // progname = "/home/janderss/Desktop/ampl_test/clnlbeam_small.nl";
  // progname = "/home/janderss/Desktop/ampl_test/svanberg.nl";
  // progname = "/home/janderss/Desktop/ampl_test/svanberg_small.nl";

//   progname = "/home/janderss/dev/casadi/trunk/misc/cuter/cuter_nl/aircrfta.nl";
//   progname = "/home/janderss/dev/casadi/trunk/misc/cuter/cuter_nl/allinitc.nl";
//   progname = "/home/janderss/dev/casadi/trunk/misc/cuter/cuter_nl/allinitc.nl";
//   progname = "/home/janderss/dev/casadi/trunk/misc/cuter/cuter_nl/alsotame.nl";
  casadi_assert(argc==2);
  progname = argv[1];
  
  
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
        
  c_ifset = (list **)Malloc((n_con + n_obj + ncom0)*sizeof(list *));
  o_ifset = c_ifset + n_con;
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
  
  using namespace std;
  cout << "ncom0 = " << ncom0 << endl;
  cout << "ncom1 = " << ncom1 << endl;
  return 0;
  
  ndv = ncond = 0;
  
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
  max_restore(&bmax);
  comwalk(combc, ncom0);
  cde_walk(obj_de, n_obj, &omax, o_cexp1st, zao,true);

  int nv = nv1 + ncom;
  for(i = 0; i < nv; i++){
    var_e[i].op = (efunc *)f_OPVARVAL1;
  }
  
  expr_nx *enx;
  for(enx = nums; enx; enx = enx->next){
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