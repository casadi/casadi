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
#ifndef ASL_READER_HPP
#define ASL_READER_HPP

extern "C"{

typedef struct v_i {
  union {
    real v;
    struct v_i *next;
  } u;
  int i;
} v_i;

typedef union vpi {
  int	i;
  real	*vp;
  cgrad	*cg;
  ograd	*og;
  v_i	*vi;
} vpi;

typedef struct dLR {
  int kind;
  union {
    double *vp;
    expr_if *eif;
    expr_va *eva;
    expr *ep;
    int i;
  } o;
} dLR;

#define dLR_UNUSED	0
#define dLR_PD		1
#define dLR_one		2
#define dLR_negone	3
#define dLR_VP		4
#define dLR_VARVAL	5
#define dLR_VARARG	6
#define dLR_IF		7

typedef char *efuncb(expr *, char *);

typedef struct expr_nx {
  efuncb *op;
  struct expr_nx *next;
  real v;
} expr_nx;

#define OPPLUS          0
#define OPMINUS         1
#define OPMULT          2
#define OPDIV           3
#define OPREM           4
#define OPPOW           5
#define OPLESS          6
#define MINLIST         11
#define MAXLIST         12
#define FLOOR           13
#define CEIL            14
#define ABS             15
#define OPUMINUS        16
#define OPOR            20
#define OPAND           21
#define LT              22
#define LE              23
#define EQ              24
#define GE              28
#define GT              29
#define NE              30
#define OPNOT           34
#define OPIFnl          35
#define OP_tanh         37
#define OP_tan          38
#define OP_sqrt         39
#define OP_sinh         40
#define OP_sin          41
#define OP_log10        42
#define OP_log          43
#define OP_exp          44
#define OP_cosh         45
#define OP_cos          46
#define OP_atanh        47
#define OP_atan2        48
#define OP_atan         49
#define OP_asinh        50
#define OP_asin         51
#define OP_acosh        52
#define OP_acos         53
#define OPSUMLIST       54
#define OPintDIV        55
#define OPprecision     56
#define OPround         57
#define OPtrunc         58
#define OPCOUNT         59
#define OPNUMBEROF      60
#define OPNUMBEROFs     61
#define OPATLEAST       62
#define OPATMOST        63
#define OPPLTERM        64
#define OPIFSYM         65
#define OPEXACTLY       66
#define OPNOTATLEAST    67
#define OPNOTATMOST     68
#define OPNOTEXACTLY    69
#define ANDLIST         70
#define ORLIST          71
#define OPIMPELSE       72
#define OP_IFF          73
#define OPALLDIFF       74
#define OP1POW          75
#define OP2POW          76
#define OPCPOW          77
#define OPFUNCALL       78
#define OPNUM           79
#define OPHOL           80
#define OPVARVAL        81
#define N_OPS           82


static char op_type[] = {
        2 /* OPPLUS */,
        2 /* OPMINUS */,
        2 /* OPMULT */,
        2 /* OPDIV */,
        2 /* OPREM */,
        2 /* OPPOW */,
        2 /* OPLESS */,
        0,
        0,
        0,
        0,
        3 /* MINLIST */,
        3 /* MAXLIST */,
        1 /* FLOOR */,
        1 /* CEIL */,
        1 /* ABS */,
        1 /* OPUMINUS */,
        0,
        0,
        0,
        2 /* OPOR */,
        2 /* OPAND */,
        2 /* LT */,
        2 /* LE */,
        2 /* EQ */,
        0,
        0,
        0,
        2 /* GE */,
        2 /* GT */,
        2 /* NE */,
        0,
        0,
        0,
        1 /* OPNOT */,
        5 /* OPIFnl */,
        0,
        1 /* OP_tanh */,
        1 /* OP_tan */,
        1 /* OP_sqrt */,
        1 /* OP_sinh */,
        1 /* OP_sin */,
        1 /* OP_log10 */,
        1 /* OP_log */,
        1 /* OP_exp */,
        1 /* OP_cosh */,
        1 /* OP_cos */,
        1 /* OP_atanh */,
        2 /* OP_atan2 */,
        1 /* OP_atan */,
        1 /* OP_asinh */,
        1 /* OP_asin */,
        1 /* OP_acosh */,
        1 /* OP_acos */,
        6 /* OPSUMLIST */,
        2 /* OPintDIV */,
        2 /* OPprecision */,
        2 /* OPround */,
        2 /* OPtrunc */,
        11 /* OPCOUNT */,
        11 /* OPNUMBEROF */,
        11 /* OPNUMBEROFs */,
        2 /* OPATLEAST */,
        2 /* OPATMOST */,
        4 /* OPPLTERM */,
        5 /* OPIFSYM */,
        2 /* OPEXACTLY */,
        2 /* OPNOTATLEAST */,
        2 /* OPNOTATMOST */,
        2 /* OPNOTEXACTLY */,
        6 /* ANDLIST */,
        6 /* ORLIST */,
        5 /* OPIMPELSE */,
        2 /* OP_IFF */,
        11 /* OPALLDIFF */,
        1 /* OP1POW */,
        1 /* f_OP2POW */,
        1 /* f_OPCPOW */,
        7 /* OPFUNCALL */,
        9 /* OPNUM */,
        8 /* OPHOL */,
        10 /* OPVARVAL */
};


}
#endif // ASL_READER_HPP
