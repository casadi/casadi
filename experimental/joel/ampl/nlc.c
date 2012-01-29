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

static maxinfo bmax, cmax, omax;

static int cwant = 2, derkind = 2, djoff = -1, owant = 2;
static int needx0check, nvseen, output_time, *vseen;
static list **c_ifset, **com_ifset, **o_ifset;
int branches = 0;
static int condlevel, fwalk, iflevel, maxa, ncond, ndv, ndvmax,
npd, nplterm, nv1, nv2, nvt, nvtmax;
static int *cvmap;
static fint nJ;
static expr *Plterms;
static char *intsk;
efunc *r_ops[N_OPS];
extern efunc *r_op[];
#define f_OPPLUS        r_op[0]
#define f_OPMINUS       r_op[1]
#define f_OPMULT        r_op[2]
#define f_OPDIV r_op[3]
#define f_OPREM r_op[4]
#define f_OPPOW r_op[5]
#define f_OPLESS        r_op[6]
#define f_MINLIST       r_op[11]
#define f_MAXLIST       r_op[12]
#define f_FLOOR r_op[13]
#define f_CEIL  r_op[14]
#define f_ABS   r_op[15]
#define f_OPUMINUS      r_op[16]
#define f_OPOR  r_op[20]
#define f_OPAND r_op[21]
#define f_LT    r_op[22]
#define f_LE    r_op[23]
#define f_EQ    r_op[24]
#define f_GE    r_op[28]
#define f_GT    r_op[29]
#define f_NE    r_op[30]
#define f_OPNOT r_op[34]
#define f_OPIFnl        r_op[35]
#define f_OP_tanh       r_op[37]
#define f_OP_tan        r_op[38]
#define f_OP_sqrt       r_op[39]
#define f_OP_sinh       r_op[40]
#define f_OP_sin        r_op[41]
#define f_OP_log10      r_op[42]
#define f_OP_log        r_op[43]
#define f_OP_exp        r_op[44]
#define f_OP_cosh       r_op[45]
#define f_OP_cos        r_op[46]
#define f_OP_atanh      r_op[47]
#define f_OP_atan2      r_op[48]
#define f_OP_atan       r_op[49]
#define f_OP_asinh      r_op[50]
#define f_OP_asin       r_op[51]
#define f_OP_acosh      r_op[52]
#define f_OP_acos       r_op[53]
#define f_OPSUMLIST     r_op[54]
#define f_OPintDIV      r_op[55]
#define f_OPprecision   r_op[56]
#define f_OPround       r_op[57]
#define f_OPtrunc       r_op[58]
#define f_OPCOUNT       r_op[59]
#define f_OPNUMBEROF    r_op[60]
#define f_OPNUMBEROFs   r_op[61]
#define f_OPATLEAST     r_op[62]
#define f_OPATMOST      r_op[63]
#define f_OPPLTERM      r_op[64]
#define f_OPIFSYM       r_op[65]
#define f_OPEXACTLY     r_op[66]
#define f_OPNOTATLEAST  r_op[67]
#define f_OPNOTATMOST   r_op[68]
#define f_OPNOTEXACTLY  r_op[69]
#define f_ANDLIST       r_op[70]
#define f_ORLIST        r_op[71]
#define f_OPIMPELSE     r_op[72]
#define f_OP_IFF        r_op[73]
#define f_OPALLDIFF     r_op[74]
#define f_OP1POW        r_op[75]
#define f_OP2POW        r_op[76]
#define f_OPCPOW        r_op[77]
#define f_OPFUNCALL     r_op[78]
#define f_OPNUM r_op[79]
#define f_OPHOL r_op[80]
#define f_OPVARVAL      r_op[81]

static real *bounds;
extern char *progname;
char *opEQ = "==", *opGE = ">=", *opGT = ">", *opLE = "<=", *opLT = "<", *opNE = "!=";
char *opAND = "&&", *opOR = "||";
char *T = "t1", *T1 = "t2";	/* temporary variable names */
static char seen[N_OPS], *declare[N_OPS];
static real negone = -1.0, one = 1.0;
extern real edagread_one;	/* &edagread_one = special derp->c value */

static expr_nx *nums;
static ograd **stog;

static Adjoint A0;
static int *acount, *dvtfree, *rownos;
static int densejac, krflag, ndvfree, ndvtmax, needT1, nvtfree,
pdts, stor_grad;

typedef struct
pdl {
  struct pdl *next;
  int i, ts;
} pdl;

static pdl *pdlbusy, *pdlfree;

static int pdlsave(pdl **opdb){
  condlevel++;
  if (iflevel++){
    *opdb = pdlbusy;
  } else {
    pdlfree = 0;
    *opdb = 0;
  }
  pdlbusy = 0;
  return ++pdts;
}

static void pdlreset(){
  pdl *p, *pnext;
  for(p = pdlbusy; p; p = pnext) {
    pnext = p->next;
    p->next = pdlfree;
    pdlfree = p;
  }
  pdlbusy = 0;
}

static void pdlrestore(pdl *opdb, int mts){
  pdl *p, *pnext;
  condlevel--;
  if (--iflevel){
    pnext = pdlbusy;
    while(p = pnext){
      pnext = p->next;
      p->next = opdb;
      opdb = p;
    }
    pnext = pdlfree;
    while((p = pnext) && p->ts >= mts){
      pnext = p->next;
      p->next = opdb;
      opdb = p;
    }
    pdlfree = pnext;
    pdlbusy = opdb;
  }
}

#define NDVTGULP 1000

static void dvtfree_inc(){
  int i, j;
  i = ndvtmax;
  j = ndvtmax += NDVTGULP;
  dvtfree = (int *)Realloc(dvtfree, j*sizeof(int));
  while(i > nvtfree)
    dvtfree[--j] = dvtfree[--i];
  nvtfree += NDVTGULP;
}

static int new_dv(){ 
  return ndvfree ? dvtfree[--ndvfree] : ndv++;
}

static void dv_free(int k){
  if (ndvfree >= nvtfree)
    dvtfree_inc();
  dvtfree[ndvfree++] = k;
}

int new_vt(){
  return nvtfree < ndvtmax ? dvtfree[nvtfree++] : nvt++;
}

static void vt_free(int k){
  if (ndvfree >= nvtfree)
    dvtfree_inc();
  dvtfree[--nvtfree] = k;
}

static list *list_freelist;

static list *new_list(list *nxt){
  list *rv;
  if (rv = list_freelist)
    list_freelist = rv->next;
  else
    rv = (list *)mem(sizeof(list));
  rv->next = nxt;
  return rv;
}

static char *fpval(real r){
  static char buf[32];
  if(r >= Infinity){
    return "1.7e308";
  } else if (r <= negInfinity){
    return "-1.7e308";
  }
  g_fmt(buf, r);
  return buf;
}

void assign(char *a, char *b){
  if(a)
    printf("\t%s = %s;\n", a, b);
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

void Goto(int n){ printf("\tgoto L%d;\n", n); }

void ifgo(char *a, char *op, char *b, int n){
  printf("\tif (%s %s %s) goto L%d;\n", a, op, b, n); 
}

void label(int n){ printf(" L%d:\n", n); }

void ifstart(char *a, char *op, char *b){ 
  printf("\tif (%s %s %s) {\n", a, op, b); 
}

void elsestart(){ 
  printf("\t} else {\n"); 
}

void elseif(char *a, char *b, char *c){ 
  printf("\t} else if (%s %s %s) {\n", a, b, c);
}

void endif(){ printf("\t}\n"); }

void endswitch(int lbl){ printf("\t}\n", lbl); }

void domain(char *s, char *x){
  printf("\tdomain_(\"%s\", &%s, %dL);\n", s, x, strlen(s));
}

void zerdiv(char *s){ 
  printf("\tzerdiv_(&%s);", s); 
}

char *call1(char *what, char *a){
  static char buf[48];
  sprintf(buf, "%s(%s)", what, a);
  return buf;
}

char *call2(char *what, char *a, char *b){
  static char buf[64];
  sprintf(buf, "%s(%s, %s)", what, a, b);
  return buf;
}

char *num(int x){
  static char buf[16];
  sprintf(buf, "%d", x);
  return buf;
}

static char op_type[] = {
  #include "op_type.hd"
};

static void Lset(linpart *L, int nlin){
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

static void lwalk(expr **ep, int neg){
  int i, j;
  expr *e = *ep;
  top:
  switch(Intcast e->op) {
    case OPNOT:
      neg = 1 - neg;
      e = *ep = e->L.e;
      goto top;
    case OPOR:
      e->op = neg ? f_OPAND : f_OPOR;
      goto andor;
    case OPAND:
      e->op = neg ? f_OPOR : f_OPAND;
      andor:
      lwalk(&e->L.e, neg);
      ep = &e->R.e;
      e = *ep;
      goto top;
    case LT:
      e->op = neg ? f_GE : f_LT;
      goto compare;
    case LE:
      e->op = neg ? f_GT : f_LE;
      goto compare;
    case EQ:
      e->op = neg ? f_NE : f_EQ;
      goto compare;
    case GE:
      e->op = neg ? f_LT : f_GE;
      goto compare;
    case GT:
      e->op = neg ? f_LE : f_GT;
      goto compare;
    case NE:
      e->op = neg ? f_EQ : f_NE;
      compare:
      i = ewalk(e->L.e);
      j = ewalk(e->R.e);
      if (i > 0)
        vt_free(i);
      if (j > 0)
        vt_free(j);
  }
}

static char *f_OPVARVAL1(expr *e, char *buf){
  int k;
  Adjoint *A;
  char *fmt;
  
  if ((k = e->a) < nv1) {
    if (k < 0) {
      fmt = "pd[%d]";
      k = (-1) - k;
    }
    else {
      A = Adjp(&adjoints[k]);
      if (!A->seen) {
        A->seen = 1;
        vseen[nvseen++] = e->a;
      }
      fmt = "x[%d]";
    }
  }
  else {
    k = cvmap[(expr_v *)e - var_e - nv1];
    if (k < 0) {
      fmt = "pd[%d]";
      k = (-1) - k;
    }
    else {
      fmt = "v[%d]";
      k += (-1);
    }
  }
  sprintf(buf, fmt, k);
  return buf;
}

char *f_OPNUM1(register expr *e, char *rv){
  g_fmt(rv,((expr_nx *)e)->v);
  return rv;
}

static int ewalk(expr *e){
  int i, j, k, k1, mts;
  expr *e1, **ep, **epe;
  expr_if *eif;
  expr_va *eva;
  expr_f *ef;
  argpair *ap, *ape;
  de *d;
  derp *dp;
  dLR *LR, **LRp;
  efunc *op;
  Adjoint *A;
  double t;
  pdl *opdb;
  static int achk[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 10, 11 };
  
  k = Intcast e->op;
  e->op = r_op[k];
  achk[k1 = op_type[k]];
  switch(k1) {
    case 11: /* OP2POW */
      e->dL = 0;
      /* no break */
      case 1: /* unary */
        j = ewalk(e->L.e);
        i = new_vt();
        if (j > 0)
          vt_free(j);
        seen[k] = 1;
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
        if (j > 0)
          vt_free(j);
        if (k > 0)
          vt_free(k);
        return e->a = i;
        
            case 3: /* vararg (min, max) */
              eva = (expr_va *)e;
              d = eva->L.d;
              condlevel++;
              if (!(i = ewalk(d->e)))
                i = new_vt();
              while(e1 = (++d)->e) {
                pdlreset();
                if ((j = ewalk(e1)) > 0)
                  vt_free(j);
              }
              condlevel--;
              if (eva->a != nv1) {
                eva->next = (expr_va *)ncond++;
                /* arrange to find this expr_va */
                /* when walking derps */
                dp = eva->R.D;
                LRp = Make_dLRp(dp->c.rp);
                *LRp = LR = (dLR *)mem(sizeof(dLR));
                LR->kind = dLR_VARARG;
                LR->o.eva = eva;
              }
              return e->a = i;
              
            case 4: /* piece-wise linear */
              LR = Make_dLR(&e->dR);
              LR->kind = nplterm++;
              LR->o.ep = Plterms;
              Plterms = e;
              return e->a = new_vt();
              
            case 5: /* if */
              eif = (expr_if *)e;
              eif->next = (expr_if*) -1;
              lwalk(&eif->e, 0);
              mts = pdlsave(&opdb);
              if ((i = ewalk(eif->T)) > 0)
                vt_free(i);
              pdlreset();
              if ((i = ewalk(eif->F)) > 0)
                vt_free(i);
              pdlrestore(opdb, mts);
              if (eif->a != nv1) {
                eif->next = (expr_if *)ncond++;
                /* arrange to find this expr_if */
                /* when walking derps */
                dp = eif->D;
                LRp = Make_dLRp(dp->c.rp);
                *LRp = LR = (dLR *)mem(sizeof(dLR));
                LR->kind = dLR_IF;
                LR->o.eif = eif;
              }
              return e->a = new_vt();
              
            case 6: /* sumlist */
              ep = e->L.ep;
              epe = e->R.ep;
              i = ewalk(*ep++);
              j = ewalk(*ep++);
              if (i > 0) {
                if (j > 0)
                  vt_free(j);
              }
              else if (!(i = j))
                i = new_vt();
              do {
                if ((j = ewalk(*ep++)) > 0)
                  vt_free(j);
              }
              while(ep < epe);
              return e->a = i;
              
            case 7: /* function call */
              ef = (expr_f *)e;
              i = k = 0;
              for(ap = ef->ap, ape = ef->ape; ap < ape; ap++)
                if (j = ewalk(ap->e))
                  if (j < 0) {
                    i = j;
                  }
                  else
                    k++;
                  if (k)
                    for(ap = ef->ap; ap < ape; ap++) {
                      op = ap->e->op;
                      if (op != f_OPNUM
                        && op != f_OPHOL
                        && op != f_OPVARVAL
                        && (j = ap->e->a) > 0)
                        vt_free(j);
                    }
                    if (!i)
                      return e->a = new_vt();
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
              A = Adjp(&adjoints[e->a]);
              A->ifset = condlevel ? 1 : 0;
              if (!A->seen) {
                A->seen = 1;
                vseen[nvseen++] = i;
              }
              break;
              /*DEBUG*/default:
              /*DEBUG*/ fprintf(Stderr, "bad opnumber %d in ewalk\n", k);
              /*DEBUG*/ exit(1);
  }
  return 0;
}

static void max_save(maxinfo *m){
  if (nvtmax < nvt)
    nvtmax = nvt;
  m->nvt = nvtmax - 1;
  m->ndv = ndvmax;
  m->ncond = ncond;
  m->needT1 = needT1;
}

static void max_restore(maxinfo *m){
  nvtmax = m->nvt + 1;
  ndvmax = m->ndv;
  ncond = m->ncond;
  needT1 = m->needT1;
}

static void vreset(list **ifset){
  Adjoint *A;
  list *lastif = 0;
  int k;
  
  while(nvseen) {
    k = vseen[--nvseen];
    A = Adjp(&adjoints[var_e[k].a]);
    if (A->ifset)
      (lastif = new_list(lastif))->item.i = k;
    *A = A0;
  }
  *ifset = lastif;
}

static void ewalkvt(expr *e, list **ifset){
  nvtfree = ndvtmax;
  ewalk(e);
  vreset(ifset);
}

static void zset(register int *z){
  register int *ze;
  if (z)
    for(ze = z + 2**z; z < ze; z += 2)
      var_e[z[1]].a = z[2];
}

static void ndvtreset(int *z){
  nvt = 1;
  nvtfree = ndvtmax;
  ndv = 0;
  ndvfree = 0;
  zset(z);
}

static void comwalk(int i, int n){
  cexp *c, *ce;
  list **ifset;
  
  if (i >= n)
    return;
  zset(zaC[i]);
  ifset = com_ifset + i;
  c = cexps + i;
  for(ce = cexps + n; c < ce; c++) {
    Lset(c->L, c->nlin);
    ndvtreset(zaC[i]);
    ewalkvt(c->e, ifset++);
    if (nvtmax < nvt)
      nvtmax = nvt;
  }
}

static void cde_walk(cde *d, int n, maxinfo *m, list **ifset, int *c1st, int **z) {
  cde *De;
  int i, i1, j, k;
  cexp1 *c1;
  
  j = *c1st++;
  for(De = d + n; d < De; d++, ifset++) {
    i = j;
    k = *c1st++;
    ndvtreset(*z++);
    while(j < k) {
      c1 = cexps1 + j++;
      Lset(c1->L, c1->nlin);
      ewalk(c1->e);
    }
    ewalkvt(d->e, ifset);
    while(i < k)
      if (!cvmap[i1 = i++ + ncom0])
        cvmap[i1] = nvt++;
      if (nvtmax < nvt)
        nvtmax = nvt;
  }
  max_save(m);
}

static list **dstored;

static void dstore(Adjoint *a, derp *d){
  a->stored = 1;
  if (dstored){
    (*dstored = new_list(*dstored))->item.D = d;
  }
}

static char * vprod(real t, int k){
  static char buf[64];
  int i;
  
  if (t == 1.)
    sprintf(buf, "x[%d]", k);
  else if (t == -1.) {
    buf[0] = '-';
    sprintf(buf+1, "x[%d]", k);
  }
  else {
    i = g_fmt(buf, t);
    buf[i++] = '*';
    sprintf(buf+i, "x[%d]", k);
  }
  return buf;
}

static char *rv_output(char *rv, char *eval, ograd *og){
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

static char *con_linadd(int i, char *s){
  cgrad *cg;
  char *s1;
  
  for(cg = Cgrad[i]; cg; cg = cg->next)
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
    return s;
}

static int nzcgrad(){
  cgrad *cg;
  int i;
  for(i = 0; i < n_con; i++)
    for(cg = Cgrad[i]; cg; cg = cg->next)
      if (cg->coef)
        return 1;
      return 0;
}

static char *cv_name(linpart *L, char *buf){
  expr_v ev;
  expr *ep = (expr *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
  return f_OPVARVAL1(ep, buf);
}

static void com_out(expr *e, linpart *L, int nlin, int k0){
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
    if (k == dLR_VP)
      printf("%s*", bufg);
    printf("%s%s", cv_name(L,vn), ";\n");
    s = res;
  }
}

static char *putout(expr *e, int i, int j, char *what, int k, int *z){
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
  assert(krflag==0);
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
  int i;
  for(i = 0; i < n_obj; i++, c1++){
    eval = putout(obj_de[i].e, c1[0], c1[1],n_obj > 1 ? "objective" : NULL, i, zao[i]);
    s = rv_output(rv, eval, Ograd[i]);
    printf("\n\treturn %s;\n", s);
  }
  
  printf("}\n");
  branches = 0;
}
  
void con_output(){
  assert(krflag==0);

  int *c1, i;
  char *s;
  printf("\n void ceval0_(real *x, real *c)\n{");
  
  if(!n_con){
    printf("}\n");
    return;
  }
  
  if(cmax.nvt){
    printf("\n\treal v[%d];\n", cmax.nvt);
  }
  if(cmax.ncond)
    printf("\tstatic int cond[%d];\n", cmax.ncond);
  
  if (cmax.needT1)
    printf("\treal t1, t2;\n");
  else if(nzcgrad())
    printf("\treal t1;\n");
  
  c1 = c_cexp1st;
  for(i = 0; i < n_con; i++, c1++) {
    s = putout(con_de[i].e, c1[0], c1[1], "constraint", i, zac[i]);
    printf("\tc[%d] = %s;\n", i, con_linadd(i,s));
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
  
  for(i = j = 0; i < N_OPS; i++)
    if (seen[i] && declare[i])
      printf(j++ ? ", %s_()" : " extern real %s_()",
      declare[i]);
    if (j)
      printf(";\n");
    for(e = Plterms; e; e = LR->o.ep) {
      LR = dLRp(e->dR);
      p = e->L.p;
      i = 2*p->n - 1;
      printf(" real bs%d[%d] = {\n", LR->kind, i);
      for(b = p->bs, be = b + i; b < be; b++) {
        g_fmt(buf, *b);
        printf("\t%s%s\n", buf,
                b + 1 < be ? "," : "};\n");
      }
    }
    
    printf(" real boundc_[1+%d+%d] /* Infinity, variable bounds, constraint bounds */ = {1.7e308",2*nv1, 2*n_con);
    b = bounds;
    be = b + 2*(n_con + nv1);
    while(b < be)
      printf(",q%s", fpval(*b++));
    printf("};\n\n");
        
    obj_output();
    con_output();
}

static void get_rownos(){
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
  cwant = owant = derkind = 1;
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
        
  needx0check = comb > 0 || derkind & 2;
  c_ifset = (list **)Malloc((n_con + n_obj + ncom0)*sizeof(list *));
  o_ifset = c_ifset + n_con;
  com_ifset = o_ifset + n_obj;
  op_type[OP1POW] = 2;
  op_type[OP2POW] = 11;
      
  stog = (ograd **)Malloc(nv1*sizeof(ograd *));

  
  declare[OP_asinh] = "asinh";
  declare[OP_asinh+1] = "asinhd";
  declare[OP_acosh] = "acosh";
  declare[OP_acosh+1] = "acoshd";
  declare[OPPLTERM] = "plterm";
    
  dvtfree = (int *)Malloc(NDVTGULP*sizeof(int));
  ndvtmax = nvtfree = NDVTGULP;
  
  if (n_con) get_rownos();
  
  vseen = (int *)Malloc((nv2 + ncom)*sizeof(int));
  cvmap = vseen + nv2;
  npd = 0;
  for(i = 0; i < ncom0; i++)
    cvmap[i] = -(++npd);
  if (ncom1)
    memset((char *)&cvmap[ncom0], 0, ncom1*sizeof(int));
  
  ndv = ncond = 0;
  
  memset((char *)adjoints, 0, amax*sizeof(real));
  if ((i = amax - nv1) > 0) {
    acount = (int *)Malloc(i*sizeof(int));
    memset((char *)acount, 0, i*sizeof(int));
  }
  
  comwalk(0,comb);
  max_save(&bmax);
  comwalk(comb, combc);
  cde_walk(con_de, n_con, &cmax, c_ifset, c_cexp1st, zac);
  max_restore(&bmax);
  comwalk(combc, ncom0);
  cde_walk(obj_de, n_obj, &omax, o_ifset, o_cexp1st, zao);
  
  int nv = nv1 + ncom;
  for(i = 0; i < nv; i++){
    var_e[i].op = (efunc *)f_OPVARVAL1;
  }
  
  expr_nx *enx;
  for(enx = nums; enx; enx = enx->next){
    enx->op = f_OPNUM1;
  }
  
  output_time = 1;
  output();
  return 0;
}
    
char *e_val(expr *e, char *buf){
  assert(e->a >= 0);
  sprintf(buf, "v[%d]", (-1) + e->a);
  return buf;
}
