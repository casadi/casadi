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

#define callb(a,b) (*(efuncb *)a->op)(a, b)

#ifdef __cplusplus
	extern "C" real log(real);
#else
	extern real log ANSI((real));
#endif

static char *pdval(dLR *d, char *buf){
  sprintf(buf, pd_fmt, d->o.i + Fortran);
  return buf;
}

static char *commute(expr *e, char *rv, char *op){
	char buf1[32], buf2[32];
	expr *e1;
	e1 = e->L.e;
	callb(e1, buf1);
	e1 = e->R.e;
	callb(e1, buf2);
	e_val(e, rv);
	if (Fortran)
		printf("\t%s = %s %s %s\n", rv, buf1, op, buf2);
	else {
		if (!strcmp(rv,buf1))
			printf("\t%s %s= %s;\n", rv, op, buf2);
		else if (!strcmp(rv,buf2))
			printf("\t%s %s= %s;\n", rv, op, buf1);
		else
			printf("\t%s = %s %s %s;\n", rv, buf1, op, buf2);
		}
	return rv;
	}

 static char *f_OPPLUS(expr *e, char *buf){ return commute(e, buf, "+"); }

 static char *f_OPSUMLIST(register expr *e, char *buf){
	expr *e1, **ep, **epe;
	char *rv;
	char buf1[32], buf2[32];

	ep = e->L.ep;
	epe = e->R.ep;
	e1 = *ep++;
	callb(e1, buf1);
	e1 = *ep++;
	callb(e1, buf2);
	rv = e_val(e, buf);
	if (Fortran) {
		printf("\t%s = %s + %s\n", rv, buf1, buf2);
		while(ep < epe) {
			e1 = *ep++;
			printf("\t%s = %s + %s\n", rv, rv,
				callb(e1, buf1));
			}
		}
	else {
		if (!strcmp(rv,buf2))
			printf("\t%s += %s;\n", rv, buf1);
		else if (!strcmp(rv,buf1))
			printf("\t%s += %s;\n", rv, buf2);
		else
			printf("\t%s = %s + %s;\n", rv, buf1, buf2);
		while(ep < epe) {
			e1 = *ep++;
			printf("\t%s += %s;\n", rv, callb(e1,buf1));
			}
		}
	return rv;
	}

 static char *binop1(char *a, char *b, char *op, char *c){
	if (Fortran)
		printf("\t%s = %s %s %s\n", a, b, op, c);
	else {
		if (!strcmp(a,b))
			printf("\t%s %s= %s;\n", a, op, c);
		else
			printf("\t%s = %s %s %s;\n", a, b, op, c);
		}
	return a;
	}

 static char *noncommute(expr *e, char *rv, char *op){
	char buf1[32], buf2[32];
	expr *e1;
	e1 = e->L.e;
	callb(e1, buf1);
	e1 = e->R.e;
	callb(e1, buf2);
	e_val(e, rv);
	return binop1(rv, buf1, op, buf2);
	}

 static char *f_OPMINUS(expr *e, char *buf){ return noncommute(e, buf, "-"); }

 static char *f_OPMULT(expr *e, char *buf){ return commute(e, buf, "*"); }

 static char *f_OPDIV(expr *e, char *rv){
	expr *e1;
	char buf[32], s1[32], s2[32];
	dLR *Ld, *Rd;
	static char *fmt[] = { "\t%s = -%s %s %s;\n", "\t%s = -%s %s %s\n" };

	e1 = e->L.e;
	callb(e1, s1);
	e1 = e->R.e;
	callb(e1, s2);
	if (e1->op != (efunc *)f_OPNUM1) {
		ifstart(s2, opEQ, Zero);
		zerdiv(s2);
		endif();
		}
	binop1(e_val(e,rv), s1, "/", s2);
	if (want_derivs) {
		Ld = dLRp(e->dL);
		Rd = dLRp(e->dR);
		if (Ld->kind) {
			pdval(Ld,buf);
			if (Rd->kind) {
				binop(buf, One, "/", s2);
				printf(fmt[Fortran], pdval(Rd,s1), rv,
					"*", buf);
				}
			else
				binop(buf, One, "/", s2);
			}
		else if (Rd->kind)
			printf(fmt[Fortran], pdval(Rd,s1),
				rv, "/", s2);
		}
	return rv;
	}
 static void opstart(register expr *e, char *who, char *val, char *L){
	register expr *e1;
	e1 = e->L.e;
	callb(e1,L);
	assign(e_val(e,val), call1(who, L));
	introuble(who, L);
	}

 static void op2start(register expr *e, char *who, char *val, char *L, char *R){
	register expr *e1;
	e1 = e->L.e;
	callb(e1,L);
	e1 = e->R.e;
	callb(e1,R);
	assign(e_val(e,val), call2(who, L, R));
	introuble2(who, L, R);
	}

 static void foppowstart(register expr *e, char *val, char *L, char *R){
	register expr *e1;

	if (Fortran) {
		e1 = e->L.e;
		callb(e1,L);
		e1 = e->R.e;
		binop(e_val(e,val), L, "**", callb(e1,R));
		}
	else
		op2start(e, "pow", val, L, R);
	}

 static char *f_OPREM(expr *e, char *rv){
	char L[32], R[32];
	op2start(e, "fmod", rv, L, R);
	return rv;
	}

 char *f_OPCPOW(expr *e, char *rv){
	char L[32], R[32], logL[32], pd[32];
	real Lv;
	dLR *Rd;

	foppowstart(e, rv, L, R);
	if (!want_derivs)
		return rv;
	Rd = dLRp(e->dR);
	if (Rd->kind != dLR_UNUSED) {
		if ((Lv = e->L.en->v) <= 0) {
			fprintf(Stderr, "%s: can't differentiate pow(%g,%s)\n",
				progname, Lv, R);
			exit(1);
			}
		g_fmt(logL, log(Lv));
		binop(pdval(Rd,pd), logL, "*", rv);
		}
	return rv;
	}

 char *f_OPPOW(register expr *e, char *rv){
	char buf[96], L[32], R[32], pdL[32], pdR[32];
	dLR *Ld, *Rd;

	foppowstart(e, rv, L, R);
	if (!want_derivs)
		return rv;
	Ld = dLRp(e->dL);
	Rd = dLRp(e->dR);
	if (Ld->kind || Rd->kind) {
		ifstart(L, opGT, Zero);
		sprintf(buf, "(%s/%s) * %s", R, L, rv);
		if (Ld->kind)
			assign(pdval(Ld,pdL), buf);
		if (Rd->kind)
			binop(pdval(Rd,pdR), call1("log",L), "*", rv);
		elseif(R, opGT, One);
		if (Ld->kind)
			assign(pdL, Zero);
		if (Rd->kind)
			assign(pdR, Zero);
		elseif(R, opEQ, One);
		if (Ld->kind)
			assign(pdL, One);
		if (Rd->kind)
			assign(pdR, Zero);
		elsestart();
		introuble2("pow'", L, R);
		endif();
		}
	return rv;
	}

 char *f_OP1POW(expr *e, char *rv)	/* f_OPPOW for R = numeric constant */
{
	char buf[96], buf1[32], L[32], R[32];
	dLR *Ld;
	real Rv;
	Long L1;
	int isint;
	char *op;

	foppowstart(e, rv, L, R);
	if (!want_derivs)
		return rv;
	Ld = dLRp(e->dL);
	if (Ld->kind) {
		Rv = ((expr_nx *)e->R.e)->v;
		if (Rv >= -2147483647. && Rv <= 2147483647.
		 && (L1 = (Long)Rv, (real)L1 == Rv)) {
			isint = L1 == 3 ? 2 : 1;
			op = opNE;
			}
		else {
			isint = 0;
			op = opGT;
			}
		ifstart(L, op, Zero);
		if (isint == 2)
			sprintf(buf, "%s*(%s*%s)", R, L, L);
		else
			sprintf(buf, "%s*(%s/%s)", R, rv, L);
		assign(pdval(Ld,buf1), buf);
		if (isint)
			elsestart();
		else
			elseif(R, opGT, One);
		assign(buf1, Zero);
		if (!isint) {
			elsestart();
			introuble2("pow'", L, R);
			}
		endif();
		}
	return rv;
	}

 char *f_OP2POW(expr *e, char *rv)	/* f_OPPOW for R = 2 */{
	register expr *e1;
	char buf[32], L[32];
	dLR *Ld;

	e1 = e->L.e;
	callb(e1, L);
	binop1(e_val(e,rv), L, "*", L);
	if (want_derivs) {
		Ld = dLRp(e->dL);
		binop1(pdval(Ld,buf), L, "+", L);
		}
	return rv;
	}

static char *f_OPLESS(register expr *e, char *rv){
	register expr *e1;
	char L[32], R[32];
	dLR *Ld, *Rd;

	e1 = e->L.e;
	callb(e1, L);
	e1 = e->R.e;
	binop1(e_val(e,rv), L, "-", callb(e1,R));
	ifstart(rv, opLT, Zero);
	assign(rv, Zero);
	if (want_derivs) {
		Ld = dLRp(e->dL);
		Rd = dLRp(e->dR);
		if (Ld->kind)
			assign(pdval(Ld,L), Zero);
		if (Rd->kind)
			assign(pdval(Rd,R), Zero);
		elsestart();
		if (Ld->kind)
			assign(L, One);
		if (Rd->kind)
			assign(R, Negone);
		}
	endif();
	return rv;
	}

static char *minmax(expr_va *e, char *rv, char *cmp){
	char buf[32], cbuf[16];
	de *d;
	expr *e1;

	d = e->L.d;
	e1 = d->e;
	callb(e1,buf);
	if (strcmp(e_val((expr*)e,rv), buf))
		assign(rv, buf);
	if (want_derivs && e->R.D) {
		sprintf(cbuf, cond_fmt, (int)e->next);
		assign(cbuf, Fortran ? "1" : "0");
		}
	for(d++; e1 = d->e; d++) {
		ifstart(rv, cmp, callb(e1,buf));
		if (d->d && want_derivs)
			assign(cbuf, num((d - e->L.d) + Fortran));
		assign(rv, buf);
		endif();
		}
	return rv;
	}

 static char *f_MINLIST(expr *e, char *rv){ return minmax((expr_va *)e, rv, opGT); }

 static char *f_MAXLIST(expr *e, char *rv){ return minmax((expr_va *)e, rv, opLT); }

 static char *f_FLOOR(register expr *e, char *rv){
	register expr *e1;
	char buf[32];

	e1 = e->L.e;
	assign(rv, call1("floor", callb(e1,buf)));
	return rv;
	}

 static char *f_CEIL(register expr *e, char *rv){
	register expr *e1;
	char buf[32];

	e1 = e->L.e;
	assign(rv, call1("ceil", callb(e1,buf)));
	return rv;
	}

 static char *f_ABS(expr *e, char *rv){
	expr *e1;
	char *pp, *pn, *sp, *sn;
	dLR *Ld;
	char buf[32], L[32];

	e1 = e->L.e;
	callb(e1,L+1);
	if (L[1] == '-') {
		sn = L + 1;
		sp = L + 2;
		pn = One;
		pp = Negone;
		}
	else {
		L[0] = '-';
		sn = L;
		sp = L + 1;
		pn = Negone;
		pp = One;
		}
	ifstart(sp, opLT, Zero);
	assign(e_val(e,rv), sn);
	if (want_derivs) {
		Ld = dLRp(e->dL);
		assign(pdval(Ld,buf), pn);
		}
	elsestart();
	assign(rv, sp);
	if (want_derivs)
		assign(buf, pp);
	endif();
	return rv;
	}

 static char *f_OPUMINUS(register expr *e, char *rv){
	char buf[32];
	register expr *e1;
	register char *s;

	e1 = e->L.e;
	callb(e1,buf+1);
	if (buf[1] == '-')
		s = buf+2;
	else {
		buf[0] = '-';
		s = buf;
		}
	assign(e_val(e,rv), s);
	return rv;
	}

 static char *f_OP_tanh(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "tanh", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		assign(pdval(Ld,pd), call1("cosh",L));
		introuble("tanh'", L);
		binop(pd, pd, "*", pd);
		}
	return rv;
	}

 static char *f_OP_tan(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "tan", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		assign(pdval(Ld,pd), call1("cos",L));
		introuble("tan'", L);
		binop(pd, pd, "*", pd);
		}
	return rv;
	}

 static char *f_OP_sqrt(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "sqrt", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		ifstart(rv, opLE, Zero);
		domain("sqrt'", rv);
		endif();
		binop(pdval(Ld,pd), Half, "/", rv);
		}
	return rv;
	}

 static char *f_OP_sinh(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "sinh", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		assign(pdval(Ld,pd), call1("cosh",L));
		introuble("sinh'", L);
		}
	return rv;
	}

 static char *f_OP_sin(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "sin", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		assign(pdval(Ld,pd), call1("cos",L));
		introuble("sin'", L);
		}
	return rv;
	}

 static char *f_OP_log10(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;
	static char Le10[24];

	if (!Le10[0])
		sprintf(Le10, "%.17g", 1./log(10.));
	opstart(e, "log10", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind)
		binop(pdval(Ld,pd), Le10, "/", L);
	return rv;
	}

 static char *f_OP_log(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "log", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind)
		binop(pdval(Ld,pd), One, "/", L);
	return rv;
	}

 static char *f_OP_exp(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "exp", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind)
		assign(pdval(Ld,pd), rv);
	return rv;
	}

 static char *f_OP_cosh(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "cosh", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		assign(pdval(Ld,pd), call1("sinh",L));
		introuble("cosh'", L);
		}
	return rv;
	}

 static char *f_OP_cos(register expr *e, char *rv){
	char L[32], pd[32];
	dLR *Ld;

	opstart(e, "cos", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind)
		assign(pdval(Ld,pd), call1("-sin",L));
	return rv;
	}

 static char *f_OP_atanh(register expr *e, char *rv){
	register expr *e1;
	char L[32], buf[64], pd[32];
	dLR *Ld;

	e1 = e->L.e;
	callb(e1,L);
	ifstart(L, opLE, Negone);
	domain("atanh", L);
	endif();
	ifstart(L, opGE, One);
	domain("atanh", L);
	endif();
	sprintf(buf, "%s * log((%s + %s) / (%s - %s))",
		Half, One, L, One, L);
	assign(e_val(e,rv), buf);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		sprintf(buf, "%s / (%s - %s*%s)", One, One, L, L);
		assign(pdval(Ld,pd), buf);
		}
	return rv;
	}

 static char *f_OP_atan2(register expr *e, char *rv){
	char buf[64], L[32], R[32];
	dLR *Ld, *Rd;

	op2start(e, "atan2", rv, L, R);
	Ld = dLRp(e->dL);
	Rd = dLRp(e->dR);
	if (want_derivs && (Ld->kind || Rd->kind)) {
		binop(T, L, "/", R);
		sprintf(buf, "%s / (%s + %s*%s)", One, One, T, T);
		assign(T1, buf);
		binop(T1, T1, "/", R);
		if (Ld->kind)
			assign(pdval(Ld,L), T1);
		if (Rd->kind) {
			sprintf(buf, "-%s*%s", T, T1);
			assign(pdval(Rd,L), buf);
			}
		}
	return rv;
	}

 static char *f_OP_atan(register expr *e, char *rv){
	char L[32], buf[96], pd[32];
	dLR *Ld;

	opstart(e, "atan", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		sprintf(buf, "%s / (%s + %s*%s)", One, One, L, L);
		assign(pdval(Ld,pd), buf);
		}
	return rv;
	}

 static char *offline(expr *e, char *who, char *whod, char *whof, char *whofd, char *rv){
	register expr *e1;
	char pd[32], L[32], buf[64];
	dLR *Ld;

	if (Fortran) {
		who = whof;
		whod = whofd;
		}
	e1 = e->L.e;
	Ld = dLRp(e->dL);
	callb(e1,L);
	if (want_derivs && Ld->kind) {
		sprintf(buf, offlfmt2, L, pdval(Ld,pd));
		who = whod;
		}
	else
		sprintf(buf, offlfmt1, L);
	assign(rv, call1(who, buf));
	return rv;
	}

 static char *f_OP_asinh(register expr *e, char *rv){ return offline(e, "asinh_", "asinhd_", "asinh", "asinhd", rv); }

 static char *f_OP_asin(register expr *e, char *rv){
	char buf[64], L[32], pd[32];
	dLR *Ld;

	opstart(e, "asin", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		sprintf(buf, "%s / sqrt(%s - %s*%s)", One, One, L, L);
		assign(pdval(Ld,pd), buf);
		}
	return rv;
	}

 static char *f_OP_acosh(register expr *e, char *rv){ return offline(e, "acosh_", "acoshd_", "acosh", "acoshd", rv); }

 static char *f_OP_acos(register expr *e, char *rv){
	char buf[64], L[32], pd[32];
	dLR *Ld;

	opstart(e, "acos", rv, L);
	Ld = dLRp(e->dL);
	if (want_derivs && Ld->kind) {
		sprintf(buf, "-%s / sqrt(%s - %s*%s)", One, One, L, L);
		assign(pdval(Ld,pd), buf);
		}
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
	if (want_derivs && (k = (int)eif->next) != -1) {
		sprintf(cbuf, cond_fmt, k);
		assign(cbuf, "1");
		}
	e = eif->T;
	callb(e,e_val(e,val));
	if (strcmp(rv,val))
		assign(rv, val);
	Goto(endlabel);
	label(elselabel);
	if (k != -1)
		assign(cbuf, "0");
	e = eif->F;
	callb(e,val);
	if (strcmp(rv,val))
		assign(rv, val);
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
		}
	else {
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
		}
	else {
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
	callb(e1, L);
	e1 = e->R.e;
	ifgo(L, cmp, callb(e1,R), lbl);
	}

 static void vf_LT(expr *e, int TL, int FL){
	if (TL)
		compare(e, opLT, TL);
	else
		compare(e, opGE, FL);
	}
#define f_LT (efuncb *)vf_LT

 static void vf_LE(expr *e, int TL, int FL){
	if (TL)
		compare(e, opLE, TL);
	else
		compare(e, opGT, FL);
	}

 static void vf_EQ(expr *e, int TL, int FL){
	if (TL)
		compare(e, opEQ, TL);
	else
		compare(e, opNE, FL);
	}

 static void vf_GE(expr *e, int TL, int FL){
	if (TL)
		compare(e, opGE, TL);
	else
		compare(e, opLT, FL);
	}

 static void vf_GT(expr *e, int TL, int FL){
	if (TL)
		compare(e, opGT, TL);
	else
		compare(e, opLE, FL);
	}

 static void vf_NE(expr *e, int TL, int FL){
	if (TL)
		compare(e, opNE, TL);
	else
		compare(e, opEQ, FL);
	}

 static void vf_OPNOT(register expr *e, int TL, int FL){
	register lfunc *op;
	e = e->L.e;
	op = (lfunc *)e->op;
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

#include "opcode.hd"	/* for OPVARVAL and OPNUM */

#define f_OPVARVAL (efuncb *)OPVARVAL
#undef f_OPNUM
#define f_OPNUM (efuncb *)OPNUM

 efuncb *
r_op[] = {
#include "c_op.hd"
};
