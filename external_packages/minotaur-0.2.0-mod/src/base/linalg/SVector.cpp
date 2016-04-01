// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 -- 2010 The MINOTAUR Team.
// 

// /**
// \file SMatrix.h
// \brief Define sparse vectors
// \author Ashutosh Mahajan, Argonne National Laboratory
//
// This file contains definitions for methods of sparse vectors.
// */


#include "MinotaurConfig.h"
#include "SVector.h"

using namespace Minotaur::LinearAlgebra;

//template <class T> SVector<T>::Svector(UInt len, UInt nnz) 
//  : Vector<T>(len), nnz_(nnz), max_(max)
//{
//  idx_ = new UInt[max_];
//  dat_ = new T[max_];
//}
//
//
//SVector<T>::Svector(UInt len, UInt nnz, UInt max)
//  : Vector(len), nnz_(nnz), max_(max)
//{
//  idx_ = new UInt[max_];
//  dat_ = new T[max_];
//}
//
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
