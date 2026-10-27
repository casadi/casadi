package org.casadi;

/** Implemented by the matrix types DM, SX and MX: what they have in common, for code over any of them. */
public interface GenericMatrix {
  long size1();
  long size2();
  long numel();
  long nnz();
  String dim();
  Sparsity sparsity();
  boolean is_scalar();
  boolean is_dense();
  boolean is_empty();
  boolean is_column();
  boolean is_row();
  boolean is_square();
  boolean is_vector();
  GenericMatrix T();
}
