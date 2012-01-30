/****************************************************************
Copyright (C) 1997 Lucent Technologies
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
#include "nlc.h"
#define asl cur_ASL
#include "assert.h"

#ifdef __cplusplus
extern "C" real log(real);
#else
extern real log ANSI((real));
#endif

char *opEQ = "==", *opGE = ">=", *opGT = ">", *opLE = "<=", *opLT = "<", *opNE = "!=";
char *opAND = "&&", *opOR = "||";

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

static char *pdval(dLR *d, char *buf){
  sprintf(buf, "pd[%d]", d->o.i + Fortran);
  return buf;
}

static char *commute(expr *e, char *rv, char *op){
  char buf1[32], buf2[32];
  expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, buf1);
  e1 = e->R.e;
  (*(efuncb *)e1->op)(e1, buf2);
  e_val(e, rv);
  if (!strcmp(rv,buf1)){
    printf("\t%s %s= %s;\n", rv, op, buf2);
  } else if (!strcmp(rv,buf2)){
    printf("\t%s %s= %s;\n", rv, op, buf1);
  } else {
    printf("\t%s = %s %s %s;\n", rv, buf1, op, buf2);
  }
  return rv;
}

static char *f_OPPLUS(expr *e, char *buf){ return commute(e, buf, "+"); }

static char *f_OPSUMLIST(register expr *e, char *buf){
  char buf1[32], buf2[32];

  expr **ep = e->L.ep;
  expr **epe = e->R.ep;
  expr *e1 = *ep++;
  (*(efuncb *)e1->op)(e1, buf1);
  e1 = *ep++;
  (*(efuncb *)e1->op)(e1, buf2);
  char *rv = e_val(e, buf);
  if (!strcmp(rv,buf2)){
    printf("\t%s += %s;\n", rv, buf1);
  } else if (!strcmp(rv,buf1)){
    printf("\t%s += %s;\n", rv, buf2);
  } else {
    printf("\t%s = %s + %s;\n", rv, buf1, buf2);
  } while(ep < epe) {
    e1 = *ep++;
    printf("\t%s += %s;\n", rv,(*(efuncb *)e1->op)(e1, buf1));
  }
  return rv;
}

void assign(char *a, char *b){
  if(a)
    printf("\t%s = %s;\n", a, b);
}

static char *binop1(char *a, char *b, char *op, char *c){
  if (!strcmp(a,b)){
    printf("\t%s %s= %s;\n", a, op, c);
  } else {
    printf("\t%s = %s %s %s;\n", a, b, op, c);
  }
  return a;
}

static char *noncommute(expr *e, char *rv, char *op){
  char buf1[32], buf2[32];
  expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, buf1);
  e1 = e->R.e;
  (*(efuncb *)e1->op)(e1, buf2);
  e_val(e, rv);
  return binop1(rv, buf1, op, buf2);
}

static char *f_OPMINUS(expr *e, char *buf){ return noncommute(e, buf, "-"); }

static char *f_OPMULT(expr *e, char *buf){ return commute(e, buf, "*"); }

static char *f_OPDIV(expr *e, char *rv){
  char buf[32], s1[32], s2[32];
  static char *fmt[] = { "\t%s = -%s %s %s;\n", "\t%s = -%s %s %s\n" };

  expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, s1);
  e1 = e->R.e;
  (*(efuncb *)e1->op)(e1, s2);
  if (e1->op != (efunc *)f_OPNUM1) {
    ifstart(s2, opEQ, "0.");
    printf("\tzerdiv_(&%s);", s2); 
    printf("\t}\n");
  }
  binop1(e_val(e,rv), s1, "/", s2);
  return rv;
}
static void opstart(register expr *e, char *who, char *val, char *L){
  register expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, L);
  printf("\t%s = %s;\n", e_val(e,val), call1(who, L));
}

static void op2start(register expr *e, char *who, char *val, char *L, char *R){
  register expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, L);
  e1 = e->R.e;
  (*(efuncb *)e1->op)(e1, R);
  printf("\t%s = %s;\n", e_val(e,val), call2(who, L, R));
}

static void foppowstart(register expr *e, char *val, char *L, char *R){
  op2start(e, "pow", val, L, R);
}

static char *f_OPREM(expr *e, char *rv){
  char L[32], R[32];
  op2start(e, "fmod", rv, L, R);
  return rv;
}

char *f_OPCPOW(expr *e, char *rv){
  char L[32], R[32];
  foppowstart(e, rv, L, R);
  return rv;
}

char *f_OPPOW(register expr *e, char *rv){
  char L[32], R[32];
  foppowstart(e, rv, L, R);
  return rv;
}

char *f_OP1POW(expr *e, char *rv){
  char L[32], R[32];
  foppowstart(e, rv, L, R);
  return rv;
}

char *f_OP2POW(expr *e, char *rv){
  char L[32];
  register expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, L);
  binop1(e_val(e,rv), L, "*", L);
  return rv;
}

static char *f_OPLESS(register expr *e, char *rv){
  char L[32], R[32];
  register expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, L);
  e1 = e->R.e;
  binop1(e_val(e,rv), L, "-",(*(efuncb *)e1->op)(e1, R));
  ifstart(rv, opLT, "0.");
  printf("\t%s = %s;\n", rv, "0.");
  printf("\t}\n");
  return rv;
}

static char *minmax(expr_va *e, char *rv, char *cmp){
  char buf[32], cbuf[16];
  de *d = e->L.d;
  expr *e1 = d->e;
  (*(efuncb *)e1->op)(e1, buf);
  if (strcmp(e_val((expr*)e,rv), buf)){
    printf("\t%s = %s;\n", rv, buf);
  }
  for(d++; e1 = d->e; d++) {
    ifstart(rv, cmp,(*(efuncb *)e1->op)(e1, buf));
    printf("\t%s = %s;\n", rv, buf);
    printf("\t}\n");
  }
  return rv;
}

static char *f_MINLIST(expr *e, char *rv){ return minmax((expr_va *)e, rv, opGT); }

static char *f_MAXLIST(expr *e, char *rv){ return minmax((expr_va *)e, rv, opLT); }

static char *f_FLOOR(register expr *e, char *rv){
  char buf[32];
  register expr *e1 = e->L.e;
  printf("\t%s = %s;\n", rv, call1("floor",(*(efuncb *)e1->op)(e1, buf)));
  return rv;
}

static char *f_CEIL(register expr *e, char *rv){
  char buf[32];
  register expr *e1 = e->L.e;
  printf("\t%s = %s;\n", rv, call1("ceil", (*(efuncb *)e1->op)(e1,buf)));
  return rv;
}

static char *f_ABS(expr *e, char *rv){
  char *pp, *pn, *sp, *sn;
  char buf[32], L[32];

  expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1,L+1);
  if (L[1] == '-') {
    sn = L + 1;
    sp = L + 2;
    pn = "1.";
    pp = "-1.";
  } else {
    L[0] = '-';
    sn = L;
    sp = L + 1;
    pn = "-1.";
    pp = "1.";
  }
  ifstart(sp, opLT, "0.");
  printf("\t%s = %s;\n", e_val(e,rv), sn);
  elsestart();
  printf("\t%s = %s;\n", rv, sp);
  printf("\t}\n");
  return rv;
}

static char *f_OPUMINUS(register expr *e, char *rv){
  char buf[32];
  register char *s;

  register expr *e1 = e->L.e;
  (*(efuncb *)e1->op)(e1,buf+1);
  if (buf[1] == '-'){
    s = buf+2;
  } else {
    buf[0] = '-';
    s = buf;
  }
  printf("\t%s = %s;\n", e_val(e,rv), s);
  return rv;
}

static char *f_OP_tanh(register expr *e, char *rv){
  char L[32];
  opstart(e, "tanh", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_tan(register expr *e, char *rv){
  char L[32];
  opstart(e, "tan", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_sqrt(register expr *e, char *rv){
  char L[32], pd[32];
  opstart(e, "sqrt", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_sinh(register expr *e, char *rv){
  char L[32], pd[32];
  opstart(e, "sinh", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_sin(register expr *e, char *rv){
  char L[32];
  opstart(e, "sin", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_log10(register expr *e, char *rv){
  char L[32];
  static char Le10[24];

  if (!Le10[0])
    sprintf(Le10, "%.17g", 1./log(10.));
  opstart(e, "log10", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_log(register expr *e, char *rv){
  char L[32];
  opstart(e, "log", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_exp(register expr *e, char *rv){
  char L[32];
  opstart(e, "exp", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_cosh(register expr *e, char *rv){
  char L[32];
  opstart(e, "cosh", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_cos(register expr *e, char *rv){
  char L[32];
  opstart(e, "cos", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_atanh(register expr *e, char *rv){
  register expr *e1;
  char L[32], buf[64];

  e1 = e->L.e;
  (*(efuncb *)e1->op)(e1,L);
  ifstart(L, opLE, "-1.");
  printf("\tdomain_(\"%s\", &%s, %dL);\n", "atanh", L, strlen("atanh"));
  printf("\t}\n");
  ifstart(L, opGE, "1.");
  printf("\tdomain_(\"%s\", &%s, %dL);\n", "atanh", L, strlen("atanh"));
  printf("\t}\n");
  sprintf(buf, "%s * log((%s + %s) / (%s - %s))",
          "0.5", "1.", L, "1.", L);
  printf("\t%s = %s;\n", e_val(e,rv), buf);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_atan2(register expr *e, char *rv){
  char buf[64], L[32], R[32];
  op2start(e, "atan2", rv, L, R);
  dLR *Ld = dLRp(e->dL);
  dLR *Rd = dLRp(e->dR);
  return rv;
}

static char *f_OP_atan(register expr *e, char *rv){
  char L[32], buf[96];
  opstart(e, "atan", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *offline(expr *e, char *who, char *whod, char *whof, char *whofd, char *rv){
  char L[32], buf[64];
  register expr *e1 = e->L.e;
  dLR *Ld = dLRp(e->dL);
  (*(efuncb *)e1->op)(e1,L);
  sprintf(buf, "&%s", L);
  printf("\t%s = %s;\n", rv, call1(who, buf));
  return rv;
}

static char *f_OP_asinh(register expr *e, char *rv){ 
  return offline(e, "asinh_", "asinhd_", "asinh", "asinhd", rv); 
}

static char *f_OP_asin(register expr *e, char *rv){
  char L[32];
  opstart(e, "asin", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

static char *f_OP_acosh(register expr *e, char *rv){ 
  return offline(e, "acosh_", "acoshd_", "acosh", "acoshd", rv); 
}

static char *f_OP_acos(register expr *e, char *rv){
  char L[32];
  opstart(e, "acos", rv, L);
  dLR *Ld = dLRp(e->dL);
  return rv;
}

typedef void lfunc(expr *, int, int);

static char *f_OPIFnl(expr *e, char *rv){
  register expr_if *eif = (expr_if *)e;
  register lfunc *op;
  int elselabel, endlabel, k;
  char cbuf[16], val[32];

  elselabel = ++branches;
  endlabel = ++branches;
  e_val(e, rv);
  e = eif->e;
  op = (lfunc *)e->op;
  op(e, 0, elselabel);
  k = -1;
  e = eif->T;
  (*(efuncb *)e->op)(e,e_val(e,val));
  if (strcmp(rv,val))
          printf("\t%s = %s;\n", rv, val);
  Goto(endlabel);
  label(elselabel);
  if (k != -1)
          printf("\t%s = %s;\n", cbuf, "0");
  e = eif->F;
  (*(efuncb *)e->op)(e,val);
  if (strcmp(rv,val))
          printf("\t%s = %s;\n", rv, val);
  label(endlabel);
  return rv;
}

static void vf_OPOR(register expr *e, int TL, int FL){
  register lfunc *op;
  register expr *e1;
  int mylbl;
  if (TL) {
    e1 = e->L.e;
    op = (lfunc *)e1->op;
    op(e1, TL, 0);
    e1 = e->R.e;
    op = (lfunc *)e1->op;
    op(e1, TL, 0);
  } else {
    mylbl = ++branches;
    e1 = e->L.e;
    op = (lfunc *)e1->op;
    op(e1, mylbl, 0);
    e1 = e->R.e;
    op = (lfunc *)e1->op;
    op(e1, 0, FL);
    label(mylbl);
  }
}

static void vf_OPAND(expr *e, int TL, int FL){
  register lfunc *op;
  register expr *e1;
  int mylbl;

  if (TL) {
    mylbl = ++branches;
    e1 = e->L.e;
    op = (lfunc *)e1->op;
    op(e1, 0, mylbl);
    e1 = e->R.e;
    op = (lfunc *)e1->op;
    op(e1, TL, 0);
    label(mylbl);
  } else {
    e1 = e->L.e;
    op = (lfunc *)e1->op;
    op(e1, 0, FL);
    e1 = e->R.e;
    op = (lfunc *)e1->op;
    op(e1, 0, FL);
  }
}

static void compare(register expr *e, char *cmp, int lbl){
  register expr *e1;
  char L[32], R[32];

  e1 = e->L.e;
  (*(efuncb *)e1->op)(e1, L);
  e1 = e->R.e;
  ifgo(L, cmp, (*(efuncb *)e1->op)(e1,R), lbl);
}

static void vf_LT(expr *e, int TL, int FL){
  if (TL){
    compare(e, opLT, TL);
  } else {
    compare(e, opGE, FL);
  }
}
#define f_LT (efuncb *)vf_LT

static void vf_LE(expr *e, int TL, int FL){
  if (TL) {
    compare(e, opLE, TL);
  } else {
    compare(e, opGT, FL);
  }
}

static void vf_EQ(expr *e, int TL, int FL){
  if (TL){
    compare(e, opEQ, TL);
  } else {
    compare(e, opNE, FL);
  }
}

static void vf_GE(expr *e, int TL, int FL){
  if (TL){
    compare(e, opGE, TL);
  } else {
    compare(e, opLT, FL);
  }
}

static void vf_GT(expr *e, int TL, int FL){
  if (TL){
    compare(e, opGT, TL);
  } else {
    compare(e, opLE, FL);
  }
}

static void vf_NE(expr *e, int TL, int FL){
  if (TL)
    compare(e, opNE, TL);
  else
    compare(e, opEQ, FL);
}

static void vf_OPNOT(register expr *e, int TL, int FL){
  e = e->L.e;
  register lfunc *op = (lfunc *)e->op;
  op(e, FL, TL);
}

#define f_OPOR (efuncb *)vf_OPOR
#define f_OPAND (efuncb *)vf_OPAND
#define f_LT (efuncb *)vf_LT
#define f_LE (efuncb *)vf_LE
#define f_EQ (efuncb *)vf_EQ
#define f_NE (efuncb *)vf_NE
#define f_GE (efuncb *)vf_GE
#define f_GT (efuncb *)vf_GT
#define f_OPNOT (efuncb *)vf_OPNOT
static char Bug[] = "$$Bug$$";

static char *f_OPIFSYM(expr *e, char *rv){ fprintf(Stderr,"OPHOL not yet implemented\n",e,rv); return Bug; }

static char * f_OPintDIV(register expr *e, char *rv){ 
  return offline(e, "intdiv_", Bug, "intdiv", Bug, rv); 
}

static char *f_OPprecision(register expr *e, char *rv){ 
  return offline(e, "precis_", Bug, "precision", Bug, rv); 
}

static char *f_OPround(register expr *e, char *rv){ 
  return offline(e, "round_", Bug, "round", Bug, rv); 
}

static char *f_OPtrunc(register expr *e, char *rv){ 
  return offline(e, "trunc_", Bug, "trunc", Bug, rv); 
}

static char *f_OPHOL(register expr *e, char *rv){ 
  fprintf(Stderr,"OPHOL not yet implemented\n",e,rv); 
  return Bug;
}

static char *f_OPFUNCALL(register expr *e, char *rv){ 
  fprintf(Stderr,"OPFUNCALL not yet implemented\n",e,rv); 
  return Bug; 
}

static char *f_OPPLTERM(register expr *e, char *rv){ 
  fprintf(Stderr,"OPPLTERM not yet implemented\n",e,rv); 
  return Bug; 
}

#define f_OPVARVAL (efuncb *)OPVARVAL
#undef f_OPNUM
#define f_OPNUM (efuncb *)OPNUM

efuncb *r_op[] = {
  f_OPPLUS,f_OPMINUS,f_OPMULT,f_OPDIV,f_OPREM,f_OPPOW,f_OPLESS,0,0,0,0,f_MINLIST,f_MAXLIST,
  f_FLOOR,f_CEIL,f_ABS,f_OPUMINUS,0,0,0,f_OPOR,f_OPAND,f_LT,f_LE,f_EQ,0,0,0,f_GE,f_GT,
  f_NE,0,0,0,f_OPNOT,f_OPIFnl,0,f_OP_tanh,f_OP_tan,f_OP_sqrt,f_OP_sinh,f_OP_sin,f_OP_log10,f_OP_log,f_OP_exp,
  f_OP_cosh,f_OP_cos,f_OP_atanh,f_OP_atan2,f_OP_atan,f_OP_asinh,f_OP_asin,f_OP_acosh,f_OP_acos,f_OPSUMLIST,
  f_OPintDIV,f_OPprecision,f_OPround,f_OPtrunc,0,0,0,0,0,f_OPPLTERM,f_OPIFSYM,0,0,0,0,0,0,0,0,0,f_OP1POW,
  f_OP2POW,f_OPCPOW,f_OPFUNCALL,f_OPNUM,f_OPHOL,f_OPVARVAL
};

char *e_val(expr *e, char *buf){
  assert(e->a >= 0);
  sprintf(buf, "v[%d]", (-1) + e->a);
  return buf;
}
