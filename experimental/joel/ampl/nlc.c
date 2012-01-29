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

#define Egulp 400

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
#include "r_ops1.hd"

static real *bounds;
extern char *progname;
char *Half = "0.5", *One = "1.", *Negone = "-1.", *Zero = "0.";
char *opEQ = "==", *opGE = ">=", *opGT = ">", *opLE = "<=",
*opLT = "<", *opNE = "!=";
char *opAND = "&&", *opOR = "||";
char *T = "t1", *T1 = "t2";	/* temporary variable names */
char *offlfmt1 = "&%s", *offlfmt2 = "&%s, &%s";
static char *assign_fmt = "\t%s = %s;\n";
static char *binop_fmt = "\t%s = %s %s %s;\n";
static char *goto_fmt = "\tgoto L%d;\n";
static char *ifgo_fmt = "\tif (%s %s %s) goto L%d;\n";
static char *label_fmt = " L%d:\n";
static char *ifstart_fmt = "\tif (%s %s %s) {\n";
static char *elseif_fmt = "\t} else if (%s %s %s) {\n";
static char *else_fmt = "\t} else {\n";
static char *endif_fmt = "\t}\n";
static char *case_fmt = "   case %d:\n";
char *cond_fmt = "cond[%d]";
static char *break_fmt = "\tbreak;\n";
static char *endswitch_fmt = "\t}\n";
static char *zerdiv_fmt = "\tzerdiv_(&%s);";
static char *call_fmt = "\t%s(&%s);\n";
static char *dv_fmt = "dv[%d]";
static char *Fortstar = "";
char *pd_fmt = "pd[%d]";
static char *tv_fmt = "v[%d]";
static char *x_fmt = "x[%d]";
static char seen[N_OPS], *declare[N_OPS];
static char *grad_fmt = "g[%d]", *jac_fmt = "J[%d]";
static char *g0_fmt = "\tg[%d] = %s;\n", *j0_fmt = "\tJ[%d] = %s;\n";
static char *eos = ";\n";
static char *star = "";
static char *Void = "void";
static char *xcheck = "\tfint wantfg = *needfg;\n\
if (xcheck(x) && wantfg == 2)\n\t\twantfg = 3;\n";
static char *xcheck0 = "\txcheck(x);\n";
static char *xcheckdcl = "(real *x)";
static char *xkind = "\tif (!(xkind & %d)) {\n\t\txkind |= %d;\n";
static real negone = -1.0, one = 1.0;
extern real edagread_one;	/* &edagread_one = special derp->c value */

static expr_nx *nums;
static ograd **stog;

static Adjoint *dwalk_one;
static real *rdwalk_one;
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

static int new_pd(){
  pdl *p;
  if (!iflevel){
    return npd++;
  }
  if(p = pdlfree){
    pdlfree = p->next;
  } else {
    p = (pdl *)mem(sizeof(pdl));
    p->i = npd++;
  }
  p->next = pdlbusy;
  pdlbusy = p;
  p->ts = pdts;
  return p->i;
}

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

static int new_vt(int deriv){
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
    printf(assign_fmt, a, b);
}

void binop(char *a, char *b, char *op, char *c){
  if(a){
    if(a == b){
      printf("\t%s %s= %s;\n", a, op, c);
    } else {
      printf(binop_fmt, a, b, op, c);
    }
  }
}

void Goto(int n){ printf(goto_fmt, n); }

void ifgo(char *a, char *op, char *b, int n){ 
  printf(ifgo_fmt, a, op, b, n); 
}

void label(int n){ printf(label_fmt, n); }

void ifstart(char *a, char *op, char *b){ 
  printf(ifstart_fmt, a, op, b); 
}

void elsestart(){ 
  printf(else_fmt); 
}

void elseif(char *a, char *b, char *c){ 
  printf(elseif_fmt, a, b, c);
}

void endif(){ printf(endif_fmt); }

void endswitch(int lbl){ printf(endswitch_fmt, lbl); }

void domain(char *s, char *x){
  printf("\tdomain_(\"%s\", &%s, %dL);\n", s, x, strlen(s));
}

void zerdiv(char *s){ 
  printf(zerdiv_fmt, s); 
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

void introuble(char *who, char *x){
  printf("\tif (errno) in_trouble(\"%s\",%s);\n",who, x);
}

void introuble2(char *who, char *x, char *y){
  printf("\tif (errno) in_trouble2(\"%s\",%s,%s);\n",who, x, y);
}

char *num(int x){
  static char buf[16];
  sprintf(buf, "%d", x);
  return buf;
}

static char op_type[] = {
  #include "op_type.hd"
};

static void dvset(register real *dp, int deriv){
  register dLR *d;
  real t;
  
  t = *dp;
  d = Make_dLR(dp);
  if (!deriv) {
    d->kind = dLR_UNUSED;
    return;
  }
  if (t == -1.) {
    d->kind = dLR_negone;
    d->o.vp = &negone;
  } else if(t == 1.){
    d->kind = dLR_one;
    d->o.vp = &one;
  } else {
    d->kind = dLR_PD;
    d->o.i = new_pd();
  }
}

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

static int ewalk(expr*, int);

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
      i = ewalk(e->L.e,0);
      j = ewalk(e->R.e,0);
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
      fmt = pd_fmt;
      k = (-1) - k;
    }
    else {
      A = Adjp(&adjoints[k]);
      if (!A->seen) {
        A->seen = 1;
        vseen[nvseen++] = e->a;
      }
      fmt = x_fmt;
    }
  }
  else {
    k = cvmap[(expr_v *)e - var_e - nv1];
    if (k < 0) {
      fmt = pd_fmt;
      k = (-1) - k;
    }
    else {
      fmt = tv_fmt;
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

static int ewalk(expr *e, int deriv){
  assert(deriv==0);
  int deriv1, i, j, k, k1, lderiv, mts, rderiv;
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
  deriv1 = deriv & 1;
  if (achk[k1 = op_type[k]] && e->a == nv1)
    deriv = deriv1 = 0;
  switch(k1) {
    case 11: /* OP2POW */
      e->dL = 0;
      /* no break */
      case 1: /* unary */
        j = ewalk(e->L.e, deriv1);
        i = new_vt(deriv);
        if (j > 0)
          vt_free(j);
        seen[k] = 1;
        reti:
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
        if (lderiv = rderiv = deriv1) {
          lderiv = achk[op_type[Intcast e->L.e->op]]
          && e->L.e->a != nv1 ? k1 : 0;
          rderiv = achk[op_type[Intcast e->R.e->op]]
          && e->R.e->a != nv1 ? k1 : 0;
        }
        i = ewalk(e->L.e,lderiv);
        j = ewalk(e->R.e,rderiv);
        k = i;
        i = new_vt(0);
        if (j > 0)
          vt_free(j);
        if (k > 0)
          vt_free(k);
        goto reti;
        
            case 3: /* vararg (min, max) */
              eva = (expr_va *)e;
              d = eva->L.d;
              condlevel++;
              if (!(i = ewalk(d->e,deriv)))
                i = new_vt(deriv);
              while(e1 = (++d)->e) {
                pdlreset();
                if ((j = ewalk(e1,deriv1)) > 0)
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
              goto reti;
              
            case 4: /* piece-wise linear */
              LR = Make_dLR(&e->dR);
              LR->kind = nplterm++;
              LR->o.ep = Plterms;
              Plterms = e;
              retnew:
              return e->a = new_vt(deriv);
              
            case 5: /* if */
              eif = (expr_if *)e;
              eif->next = (expr_if*) -1;
              lwalk(&eif->e, 0);
              mts = pdlsave(&opdb);
              if ((i = ewalk(eif->T,deriv1)) > 0)
                vt_free(i);
              pdlreset();
              if ((i = ewalk(eif->F,deriv1)) > 0)
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
              goto retnew;
              
            case 6: /* sumlist */
              ep = e->L.ep;
              epe = e->R.ep;
              i = ewalk(*ep++, deriv);
              if (i < 0)
                deriv = deriv1;
              j = ewalk(*ep++, deriv);
              if (i > 0) {
                if (j > 0)
                  vt_free(j);
              }
              else if (!(i = j))
                i = new_vt(deriv);
              do {
                if ((j = ewalk(*ep++, deriv1)) > 0)
                  vt_free(j);
              }
              while(ep < epe);
              goto reti;
              
            case 7: /* function call */
              ef = (expr_f *)e;
              i = k = 0;
              for(ap = ef->ap, ape = ef->ape; ap < ape; ap++)
                if (j = ewalk(ap->e,deriv))
                  if (j < 0) {
                    i = j;
                    deriv = deriv1;
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
                      goto retnew;
                    goto reti;
                    
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
  ewalk(e,0);
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

static void dwalk0(derp *, derp*);
static void dwalk1(derp *, derp*);

static void dwalk0(derp *d, derp *d0){
  int i;
  
  if (d == d0)
    return;
  do {
    if ((i = d->b.rp - adjoints - nv1) > 0)
      acount[i]++;
  }
  while((d = d->next) != d0);
}

static void cond_magic(dLR *LR){
  expr_va *eva;
  expr_if *eif;
  derp *D, *D0, Dsave;
  de *d;
  int iftype;
  
  if (LR->kind == dLR_IF) {
    eif = LR->o.eif;
    D = eif->D;
    iftype = 1;
  }
  else {
    eva = LR->o.eva;
    D = eva->R.D;
    iftype = 0;
  }
  Dsave = *D;
  D->c.rp = &edagread_one;
  if (iftype) {
    D->a.rp = eif->Tv.rp;
    D->next = eif->dT;
    dwalk1(D, eif->d0);
    D->a.rp = eif->Fv.rp;
    D->next = eif->dF;
    dwalk1(D, eif->d0);
  }
  else {
    D0 = eva->d0;
    for(d = eva->L.d; d->e; d++) {
      D->a.rp = d->dv.rp;
      D->next = d->d;
      dwalk1(D, D0);
    }
    *D = Dsave;
  }
  *D = Dsave;
}

static void dwalk1(derp *d,  derp *d0){
  Adjoint *a, *b;
  dLR *c;
  int i, j, j0, neg;
  
  if (d == d0)
    return;
  dwalk0(d, d0);
  do {
    j = -1;
    if ((j0 = d->b.rp - adjoints - nv1) > 0)
      j = --acount[j0];
    b = Adjp(d->b.rp);
    if (!b->storage)	/* defined var used only in if */
      continue;
    i = d->a.rp - adjoints;
    if (maxa < i)
      maxa = i;
    c = dLRp(*d->c.rp);
    if (c->kind >= dLR_VARARG) {
      if (!output_time)
        cond_magic(c);
      goto jfree;
    }
    if (i < nv1 && !fwalk)
      continue;
    a = Adjp(d->a.rp);
    switch(a->storage) {
      case STOR_VP:
        a->storage = STOR_DV;
        a->o.i = new_dv();
      default:
        goto jfree;
      case 0:
        break;
    }
    if (b->storage == STOR_IMPLICIT) {
      a->neg = b->neg;
      switch(c->kind) {
        case dLR_negone:
          a->neg = 1 - b->neg;
          /* no break */
          case dLR_one:
            a->storage = STOR_IMPLICIT;
            break;
          case dLR_VP:
            a->storage = STOR_VP;
            a->o.vp = c->o.vp;
            break;
          case dLR_PD:
            a->storage = STOR_PD;
            a->o.i = c->o.i;
            break;
          case dLR_VARVAL:
            a->storage = STOR_VARVAL;
            a->o.i = c->o.i;
            break;
          default:
            /*DEBUG*/ fprintf(Stderr,
            "\nBad c->kind = %d with STOR_IMPLICIT\n",
                              c->kind);
                              /*DEBUG*/ exit(10);
      }
      jfree:
      if (!j && b->storage == STOR_DV)
        dv_free(b->o.i);
      continue;
    }
    neg = b->neg;
    switch(c->kind) {
      case dLR_negone:
        neg = 1 - neg;
        /* no break */
        case dLR_one:
          if (j || fwalk && b->storage == STOR_PD)
            break;
          a->storage = b->storage;
          a->o = b->o;
          a->neg = neg;
          continue;
    }
    a->storage = STOR_DV;
    a->o.i = !j && b->storage == STOR_DV ? b->o.i : new_dv();
  }
  while((d = d->next) != d0);
  if (b && b->storage == STOR_DV)
    dv_free(b->o.i);
}

static void areset(cgrad *cg){
  Adjoint *A;
  int k;
  while(nvseen) {
    A = Adjp(&adjoints[var_e[vseen[--nvseen]].a]);
    *A = A0;
  }
  
  for(k = nv1; k <= maxa; k++) {
    A = Adjp(&adjoints[k]);
    *A = A0;
  }
  maxa = 0;
  for(; cg; cg = cg->next) {
    A = Adjp(&adjoints[cg->varno]);
    *A = A0;
  }
}

static void dwalk(register derp *d, cgrad *cg, int *z, list *L){
  register real *rp;
  register derp *d1;
  int i;
  Adjoint *A;
  static cgrad *cg0;
  
  areset(cg0);
  cg0 = cg;
  ndv = 0;
  ndvfree = 0;
  zset(z);
  for(; L; L = L->next) {
    i = L->item.i;
    A = Adjp(&adjoints[var_e[i].a]);
    A->seen = 1;
    vseen[nvseen++] = i;
    A->storage = STOR_DEFV;
    A->o.i = ndv++;
  }
  if (d1 = d) {
    dwalk_one->storage = STOR_IMPLICIT;
    rp = d->b.rp;
    do {
      if (d1->b.rp == rp)
        d1->b.rp = rdwalk_one;
    }
    while(d1 = d1->next);
    dwalk1(d, 0);
  }
  if (ndvmax < ndv)
    ndvmax = ndv;
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
      ewalk(c1->e, 0);
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
  if (dstored)
    (*dstored = new_list(*dstored))->item.D = d;
}

static char * vprod(real t, int k){
  static char buf[64];
  int i;
  
  if (t == 1.)
    sprintf(buf, x_fmt, k);
  else if (t == -1.) {
    buf[0] = '-';
    sprintf(buf+1, x_fmt, k);
  }
  else {
    i = g_fmt(buf, t);
    buf[i++] = '*';
    sprintf(buf+i, x_fmt, k);
  }
  return buf;
}

static void plcommon(int deriv){
  register plterm *p;
  register expr *e;
  register dLR *LR;
  
  for(e = Plterms; e; e = LR->o.ep) {
    p = e->L.p;
    LR = dLRp(e->dR);
    printf("\tcommon /pltcom/ bs%d\n\
    double precision bs%d(%d)\n", LR->kind, LR->kind, 2*p->n - 1);
  }
  printf(
  "\tcommon /xkindc/ xkind\n\tinteger xkind\n\tsave /xkindc/\n");
  
  if (npd)
    printf(
    "\tcommon /pdcomn/ pd\n\tdouble precision pd(%d)\n\tsave /pdcomn/\n",
          npd);
          if (!deriv && needx0check)
            printf("\tlogical xchk\n");
}

static char *rv_output(char *rv, char *eval, ograd *og){
  char *s;
  
  for(; og; og = og->next)
    if (og->coef)
      break;
    if (og) {
      s = vprod(og->coef, og->varno);
      if (strcmp(eval, Zero))
        binop(rv, eval, "+", s);
      else
        assign(rv, s);
      while(og = og->next)
        if (og->coef)
          binop(rv, rv, "+",
          vprod(og->coef, og->varno));
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
      if (strcmp(s,Zero))
        binop(T, s, "+", s1);
      else
        assign(T, s1);
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
  expr *ep;
  
  ep = (expr *)((char *)L->v.rp - ((char *)&ev.v - (char *)&ev));
  return f_OPVARVAL1(ep, buf);
}

static void com_out(expr *e, linpart *L, int nlin, int k0){
  char buf[32], bufg[32], res[32], vn[32];
  char *s;
  double t;
  linpart *Le;
  int asg, j, k, op;
  efuncb *eb;
  dLR *d;
  int Fortran1 = -1;
  
  printf("\n%s\t/*** defined variable %d ***/\n\n", star, k0+1);
  if ((j = cvmap[k0]) < 0) {
    j = Fortran1 - j;
    s = pd_fmt;
  }
  else {
    eb = (efuncb *)e->op;
    if (eb != (efuncb *)OPNUM && eb != (efuncb *)OPVARVAL)
      e->a = j;
    j += (-1);
    s = tv_fmt;
  }
  sprintf(res, s, j);
  s = callb(e,buf);
  if (!L) {
    if (strcmp(res,s))
      assign(res, s);
    return;
  }
  asg = !strcmp(s, Zero);
  for(Le = L + nlin; L < Le; L++, asg = 0) {
    d = dLRp(L->fac);
    op = '+';
    switch(k = d->kind) {
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
      default:/*DEBUG*/ fprintf(Stderr,
        "Bad d->kind = %d in com_walk\n", d->kind);
        /*DEBUG*/ exit(14);
    }
    if (asg)
      printf(op == '-' ? "\t%s = -" : "\t%s = ", res);
    else if (res != s)
      printf("\t%s = %s %c ", res, s, op);
    else
      printf("\t%s %c= ", res, op);
    if (k == dLR_VP)
      printf("%s*", bufg);
    printf("%s%s", cv_name(L,vn), eos);
    s = res;
  }
}

static char *putout(expr *e, int i, int j, char *what, int k, int *z){
  static char buf[32];
  cexp1 *c1;
  
  ++k;
  ndvtreset(z);
  if (i < j) {
    if (what)
      printf("\n\n%s\t/*** defined variables for %s %d ***/\n",
      star, what, k);
    for(c1 = cexps1 + i; i < j; i++, c1++)
      com_out(c1->e, c1->L, c1->nlin, i + ncom0);
  }
  printf(what ? "\n%s  /***  %s %d  ***/\n\n"
  : "\n%s  /***  objective ***/\n\n", star, what, k);
  return callb(e,buf);
}

static void funnel_set(funnel *f, char *gj0_fmt){
  assert(0);
}

void obj_output(int deriv){
  int *c1, i;
  ograd *og;
  static char rv[] = "rv";
  static char *header[4] = {
    "0_(fint *nobj, real *x)",
    "_(fint *nobj, fint *needfg, real *x, real *g)",
    "0_(nobj, x)\n\tfint *nobj;\n real *x;",
    "_(nobj, needfg, x, g)\n fint *nobj, *needfg;\n real *x, *g;"
  };
  char *eval, *s;
  
  printf(" real\nfeval%s\n", header[deriv + krflag]);
  
  if (!n_obj) {
    /*{*/ printf("{ return 0.; }\n");
    return;
  }
  printf("{");
  for(og = 0, i = 0; i < n_obj; i++) {
    for(og = Ograd[i]; og; og = og->next)
      if (og->coef)
        goto break2;
  }
  break2:
  if (omax.nvt) {
    printf("\n\treal v[%d]", omax.nvt);
    if (omax.ndv)
      printf(", dv[%d]", omax.ndv);
    printf("%s;\n", og ? ", rv" : "");
  }
  else if (og)
    printf("\n\treal rv;\n");
  
  if (omax.ncond)
    printf("\tstatic int cond[%d];\n", omax.ncond);
  
  if (omax.needT1)
    printf("\treal t1, t2;\n");
  
  if ((f_b || f_o) && deriv) {
    printf("\tif (wantfg & 2) {\n");
    if (f_b)
      printf("\t\tfunnelb(x);\n");
    if (f_o)
      funnel_set(f_o, g0_fmt);
    printf("\t\t}\n");
  }
  
  printf(n_obj > 1 ? "\n\tswitch(*nobj) {\n" : "\n");
  c1 = o_cexp1st;
  for(i = 0; i < n_obj; i++, c1++) {
    if (n_obj > 1)
      printf("\n  case %d:\n", i);
    eval = putout(obj_de[i].e, c1[0], c1[1],
                  n_obj > 1 ? "objective" : NULL, i, zao[i]);
                  s = rv_output(rv, eval, Ograd[i]);
                  printf("\n\treturn %s;\n", s);
      }
      if (n_obj > 1)
        printf("\n\t}\n");
      printf("}\n");
      branches = 0;
  }
  
void con_output(int deriv){
    assert(deriv==0);

  int *c1, i;
  char *s;
  static char *header[4] = {
    "0_(real *x, real *c)",
    "_(fint *needfg, real *x, real *c, real *J)",
    "0_(x, c)\n real *x, *c;",
    "_(needfg, x, c, J)\n fint *needfg;\n real *x, *c, *J;"
  };
  printf("\n void\nceval%s\n{", header[krflag]);
  
  if (!n_con) {
    printf("}\n");
    return;
  } /*}*/
  
  if (cmax.nvt) {
    printf("\n\treal v[%d]", cmax.nvt);
    printf(";\n", cmax.ndv);
  }
  if (cmax.ncond)
    printf("\tstatic int cond[%d];\n", cmax.ncond);
  
  if (cmax.needT1)
    printf("\treal t1, t2;\n");
  else if (nzcgrad())
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
  static char *fhead[2] = { "(real *x)", "(x) real *x;" };
  
  printf("#include \"math.h\"\n#include \"errno.h\"\n\
  #ifndef fint\n\
  #ifndef Long\n\
  #include \"arith.h\"	/* for Long */\n\
  #ifndef Long\n\
  #define Long long\n\
  #endif\n\
  #endif\n\
  #define fint Long\n\
  #endif\n\
  #ifndef real\n\
  #define real double\n\
  #endif\n");
  if (!krflag)
    printf("#ifdef __cplusplus\nextern \"C\" {\n#endif\n");
  if (krflag < 2)
    printf(" %s\n %s\n %s\n %s\n %s\n %s\n %s\n %s\n",
    "real acosh_(real *);",
    "real asinh_(real *);",
    "real acoshd_(real *, real *);",
    "real asinhd_(real *, real *);",
    "void in_trouble(char *, real);",
    "void in_trouble2(char *, real, real);",
    "void domain_(char *, real *, fint);",
    "void zerdiv_(real *);");
  printf(" fint auxcom_[1] = { %d /* nlc */ };\n", nlc);
  printf(" fint funcom_[%d] = {\n\
  %d /* nvar */,\n\
  %d /* nobj */,\n\
  %d /* ncon */,\n\
  %d /* nzc */,\n\
  %d /* densejac */",
  n_con ? nzc + nv1 + n_obj + 6 : n_obj + 5,
  nv1, n_obj, n_con, 0 ? nv1*n_con : nzc, 0);
  
  if (n_obj) {
    printf(",\n\n\t/* objtype (0 = minimize, 1 = maximize) */\n");
    for(i = 0; i < n_obj; i++)
      printf("%s\n\t%d", i ? "," : "", objtype[i]);
  }
  
  if (n_con) {
    printf(",\n\n\t/* colstarts */\n");
    for(i = 0; i <= n_var; i++)
      printf("\t%d,\n", A_colstarts[i]);
    printf("\n\t/* rownos */\n\t1");
    for(i = 1; i < nzc; i++)
      printf(",\n\t%d", rownos[i]);
  }
  printf(" };\n\n");
  
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
    
    printf(" real boundc_[1+%d+%d] /* Infinity, variable bounds, constraint bounds */ = {\n\t\t1.7e308",2*nv1, 2*n_con);
    b = bounds;
    be = b + 2*(n_con + nv1);
    while(b < be)
      printf(",\n\t\t%s", fpval(*b++));
    printf("};\n\n");
    
    printf(" real x0comn_[%d] = {\n", nv1);
    for(i = 0; i < nv1; i++)
      printf("%s\t\t%s", i ? ",\n" : "", fpval(X0[i]));
    printf(" };\n\n");
    
    if (npd)
      printf(" static real pd[%d];\n", npd);
            
      if (f_b && derkind & 2) {
        printf("\n static void\nfunnelb%s\n{\n", fhead[krflag>>1]);
        if ((i = bmax.ndv) > 0)
          printf("\treal dv[%d];\n", i);
        funnel_set(f_b, "Botch<%d> = %s;");
        printf("\t}\n");
      }
            
      if (needx0check) {
        printf("static real old_x[%d];\nstatic int xkind = -1;\n\n\
        static int\nxcheck%s\n{\n\treal", nv1, xcheckdcl);
        if (comb > 0) {
          printf(" *x0 = x,");
          x0 = "x0";
        } else {
          x0 = "x";
        }
        printf(" *x1 = old_x, *xe = x + %d;\n", nv1);
        if ((i = bmax.nvt) > 0)
          printf("\treal v[%d];\n", i);
        printf("\terrno = 0;\n\
        if (xkind >= 0) {\n\t\twhile(*%s++ == *x1++)\n\
          \tif (%s == xe)\n\t\t\t\treturn 0;\n\t\t--%s, --x1;\n\t\t}\n\
          do *x1++ = *%s++;\n\t\twhile(%s < xe);\n\txkind = 0;\n",
            x0,x0,x0,x0,x0);
          for(i = 0, c = cexps; i < comb; c++, i++)
            com_out(c->e, c->L, c->nlin, i);
          printf("\treturn 1;\n\t}\n");
  }
  for(i = 1; i < 3; i++) {
    if (owant & i)
      obj_output(i-1);
    if (cwant & i)
      con_output(i-1);
  }
  if (!krflag)
    printf("#ifdef __cplusplus\n\t}\n#endif\n");
}

static void cant(char *s1, char *s2){
  fprintf(Stderr, "Can't open %s", s1);
      if (s2)
        fprintf(Stderr, " or %s.nl", s2);
      fputc('\n', Stderr);
      exit(1);
}
    
static void get_rownos(){
  int i = n_con, i1, j, j1;
  cgrad *cg;
  memset((char *)rownos, 0, nzc*sizeof(int));
  while(i1 = i)
    for(cg = Cgrad[--i]; cg; cg = cg->next)
      rownos[cg->goff] = i1;
    for(i = 0; i <= n_var; i++)
      A_colstarts[i]++;
    if (intsk) {
      i1 = j = 0;
      --intsk;
      for(i = 1; i <= n_var; i++) {
        j1 = A_colstarts[i] - 1;
        if (!intsk[i])
          for(; j < j1; j++)
            rownos[i1++] = rownos[j];
          A_colstarts[i] = i1 + 1;
        j = j1;
      }
      ++intsk;
      nzc = i1;
      if (nlvbi && (nlvci < nlvc || nlvoi < nlvo)
        || nlvci && nlvoi < nlvo) {
        for(i = 0; i < n_var; i++)
          A_colstarts[i] = A_colstarts[i+1];
        i = n_con;
      while(--i >= 0)
        for(cg = Cgrad[i]; cg; cg = cg->next)
          if (!intsk[cg->varno])
            cg->goff = --A_colstarts[cg->varno] - 1;
      }
    }
}
    
static void nlvzap(int i, int j){
      memset(intsk + i - j, 1, j);
}
    
int main(int argc, char **argv){
  ASL_alloc(ASL_read_fg);
  progname = "../examples/cork.nl";
  g_fmt_decpt = 1;
  want_derivs = 0;
  cwant = owant = derkind = 1;
  return_nofile = 1;
  fint L;
  FILE *nl = jacdim0(progname, L = strlen(progname));

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
  
  dvset(&edagread_one, 1);
  #ifdef X64_bit_pointers
  {Adjoint *Ap = (Adjoint *)Malloc(amax*sizeof(Adjoint));
  memset((char *)Ap, 0, amax*sizeof(Adjoint));
  for(i = 0; i < amax; i++)
    *(Adjoint **)&adjoints[i] = Ap++;
  }
  #else
  memset((char *)adjoints, 0, amax*sizeof(real));
  #endif
  rdwalk_one = &adjoints[nv1];
  dwalk_one = Adjp(rdwalk_one);
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
  for(i = 0; i < nv; i++)
    var_e[i].op = (efunc *)f_OPVARVAL1;
  
  expr_nx *enx;
  for(enx = nums; enx; enx = enx->next)
    enx->op = f_OPNUM1;
  
  if (n_obj || n_con) {
    output_time = 1;
    output();
  }
  return 0;
}
    
char *e_val(expr *e, char *buf){
  int i;
  if ((i = e->a) >= 0)
    sprintf(buf, tv_fmt, (-1) + i);
  else
    sprintf(buf, pd_fmt, (-1) - i);
  return buf;
}
