//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file OpCode.h
 * \brief Declare the OpCodes used in Minotaur.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAUROPCODE_H
#define MINOTAUROPCODE_H

namespace Minotaur {
typedef enum {
  OpAbs,
  OpAcos,
  OpAcosh,
  OpAsin,
  OpAsinh,
  OpAtan,
  OpAtanh,
  OpCeil,
  OpCos,
  OpCosh,
  OpCPow, // k^x where k is constant.
  OpDiv,
  OpExp,
  OpFloor,
  OpInt,
  OpIntDiv,
  OpLog,
  OpLog10,
  OpMinus,
  OpMult,
  OpNone,
  OpNum,
  OpPlus,
  OpPow,  // y^x 
  OpPowK, // x^k, where k is constant. 
  OpRound,
  OpSin,
  OpSinh,
  OpSqr, // x^2, same as OP2POW in AMPL.
  OpSqrt,
  OpSumList,
  OpTan,
  OpTanh,
  OpUMinus,
  OpVar
} OpCode;
}

#endif

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
