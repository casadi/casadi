// NOLINT(legal/copyright)
// SYMBOL "sp_nnz"
inline
casadi_int casadi_sp_nnz(const casadi_int* sp) {
  return sp[2+sp[1]];
}
