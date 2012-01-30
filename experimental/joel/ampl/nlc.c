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

#include "nlp.h"
#undef f_OPNUM
#include "opcode.hd"
#include "nlc.h"
#define asl ((ASL_fg*)cur_ASL)
#include "assert.h"

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

real *bounds;
extern char *progname;
char *T = "t1";
expr_nx *nums;

int *dvtfree, *rownos;
int ndvfree, ndvtmax, needT1, nvtfree;

#define NDVTGULP 1000

void dvtfree_inc(){
  int i, j;
  i = ndvtmax;
  j = ndvtmax += NDVTGULP;
  dvtfree = (int *)Realloc(dvtfree, j*sizeof(int));
  while(i > nvtfree)
    dvtfree[--j] = dvtfree[--i];
  nvtfree += NDVTGULP;
}

int new_vt(){
  return nvtfree < ndvtmax ? dvtfree[nvtfree++] : nvt++;
}

void vt_free(int k){
  if (ndvfree >= nvtfree){
    dvtfree_inc();
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

static int ewalk(expr*);

char *f_OPVARVAL1(expr *e, char *buf){
  int k = e->a;
  assert(k>=0);
  assert(k<nv1);
  sprintf(buf, "x[%d]", k);
  return buf;
}

char *f_OPNUM1(register expr *e, char *rv){
  g_fmt(rv,((expr_nx *)e)->v);
  return rv;
}

static int ewalk(expr *e){
  static int achk[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 10, 11 };
  int k = (int)(e->op);
  e->op = r_op[k];
  int k1 = op_type[k];

  int i, j;
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
        return e->a = i;
        
      case 2: /* binary */
        k1 = 1;
        switch(k) {
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
        i = ewalk(e->L.e);
        j = ewalk(e->R.e);
        k = i;
        i = new_vt();
        if (j > 0){
          vt_free(j);
        }
        if (k > 0){
          vt_free(k);
        }
        return e->a = i;
        
      case 3: /* vararg (min, max) */
        assert(0);
        
      case 4: /* piece-wise linear */
        assert(0);
        
      case 5: /* if */
        assert(0);
        
      case 6: /* sumlist */
        ep = e->L.ep;
        epe = e->R.ep;
        i = ewalk(*ep++);
        j = ewalk(*ep++);
        if (i > 0) {
          if (j > 0)
            vt_free(j);
        } else if (!(i = j)){
          i = new_vt();
        } 
        do {
          if ((j = ewalk(*ep++)) > 0){
            vt_free(j);
          }
        } while(ep < epe);
        return e->a = i;
        
      case 7: /* function call */
        ef = (expr_f *)e;
        i = k = 0;
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
        break;
        /*DEBUG*/default:
        /*DEBUG*/ fprintf(Stderr, "bad opnumber %d in ewalk\n", k);
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

void cde_walk(cde *d, int n, maxinfo *m, int *c1st, int **z){
  int kk;
  printf("c1st = {");
  for(kk=0; kk<=n; ++kk){
    printf("%d,",c1st[kk]);
  }
  printf("}\n");
  
  int j = *c1st++;
  cde *De = d + n;
  for(; De > d; d++){
    int i = j;
    int k = *c1st++;
    ndvtreset(*z++);
    while(j < k){
      cexp1 *c1 = cexps1 + j++;
      Lset(c1->L, c1->nlin);
      ewalk(c1->e);
    }
    nvtfree = ndvtmax;
    ewalk(d->e);
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
    s = vprod(og->coef, og->varno);
    if (strcmp(eval, "0.")){
      binop(rv, eval, "+", s);
    } else {
      printf("\t%s = %s;\n", rv, s);
    }
    while(og = og->next){
      if (og->coef){
        binop(rv, rv, "+",vprod(og->coef, og->varno));
      }
    }
    return rv;
  }
  return eval;
}

char *con_linadd(int i, char *s){
  cgrad *cg;
  char *s1;
  
  for(cg = Cgrad[i]; cg; cg = cg->next){
    if (cg->coef) {
      s1 = vprod(cg->coef, cg->varno);
      if (strcmp(s,"0.")){
        binop(T, s, "+", s1);
      } else {
        printf("\t%s = %s;\n", T, s1);
      }
      s = T;
      while(cg = cg->next)
        if (cg->coef)
          binop(s, s, "+",
                vprod(cg->coef, cg->varno));
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
  assert(j>=0);
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
  assert(n_obj==1);
  printf(" real feval0_(long *nobj, real *x){\n");
  ograd *og;
  for(og = Ograd[0]; og; og = og->next){
    if (og->coef){
      break;
    }
  }
  assert(og==0);
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
    eval = putout(obj_de[i].e, c1[0], c1[1],n_obj > 1 ? "objective" : NULL, i, zao[i]);
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
  int i, j;
  plterm *p;
  expr *e;
  real *b, *be;
  dLR *LR;
  char buf[32], *x0;
  cexp *c;
  printf("#include \"math.h\"\n#define real double\n");
    
  printf(" real boundc_[1+%d+%d] /* Infinity, variable bounds, constraint bounds */ = {1.7e308",2*nv1, 2*n_con);
  b = bounds;
  be = b + 2*(n_con + nv1);
  while(b < be){
    printf(",q%s", fpval(*b++));
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
  ASL_alloc(ASL_read_fg);
  g_fmt_decpt = 1;
  want_derivs = 0;
  return_nofile = 1;
  progname = "../examples/cork.nl";
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
    
  LUv = bounds = (real *)Malloc((3*nv1+2*n_con)*sizeof(real));
  LUrhs = LUv + 2*nv1;
  X0 = LUrhs + 2*n_con;
      
  size_expr_n = sizeof(expr_nx);
  fg_read(nl,0);
        
  c_ifset = (list **)Malloc((n_con + n_obj + ncom0)*sizeof(list *));
  o_ifset = c_ifset + n_con;
  com_ifset = o_ifset + n_obj;
  op_type[OP1POW] = 2;
  op_type[OP2POW] = 11;
  
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
  
  ndv = ncond = 0;
  
  memset((char *)adjoints, 0, amax*sizeof(real));
  
  comwalk(0,comb);
  max_save(&bmax);
  comwalk(comb, combc);
  cde_walk(con_de, n_con, &cmax, c_cexp1st, zac);
  max_restore(&bmax);
  comwalk(combc, ncom0);
  cde_walk(obj_de, n_obj, &omax, o_cexp1st, zao);

  int nv = nv1 + ncom;
  for(i = 0; i < nv; i++){
    var_e[i].op = (efunc *)f_OPVARVAL1;
  }
  
  expr_nx *enx;
  for(enx = nums; enx; enx = enx->next){
    enx->op = f_OPNUM1;
  }
  
  output();
  return 0;
}
    
