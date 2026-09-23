//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// CBOR (RFC 8949) encoding into a caller-provided byte array, and decoding from one.
// Each casadi_cbor_write_X writes one item at p and returns the position after it;
// casadi_cbor_X_size gives the byte count up front, so space can be reserved first.
// Each casadi_cbor_read_X reads one item at p, before end: the position after it, or null if
// the item is not an X or does not end before end. Nothing is written then.

// SYMBOL "cbor_little_endian"
inline
int casadi_cbor_little_endian(void) {
  // Is the host little-endian? Then the bytes of a double need reversing to be big-endian
  union { unsigned short s; unsigned char c[2]; } e;
  e.s = 1;
  return e.c[0];
}

// SYMBOL "cbor_head_size"
inline
casadi_int casadi_cbor_head_size(casadi_int v) {
  // v >= 0: inline below 24, else in the next 1, 2, 4 or 8 bytes. Bytes are counted with
  // shifts by 8 only, so no constant or shift exceeds casadi_int, whatever its width
  casadi_int k;
  if (v < 24) return 1;
  for (k=0; v; v >>= 8) k++;
  return 1 + (k==1 ? 1 : k==2 ? 2 : k<=4 ? 4 : 8);
}

// SYMBOL "cbor_write_head"
inline
unsigned char* casadi_cbor_write_head(unsigned char* p, casadi_int major, casadi_int v) {
  // Major type in the top 3 bits; v inline if < 24, else in the next 1, 2, 4 or 8 bytes
  casadi_int n, i;
  n = casadi_cbor_head_size(v);
  if (n==1) {
    *p++ = (unsigned char) (major*32 + v);
    return p;
  }
  *p++ = (unsigned char) (major*32 + (n==2 ? 24 : n==3 ? 25 : n==5 ? 26 : 27));
  for (i=n-2; i>=0; --i) *p++ = (unsigned char) ((v >> (8*i)) & 0xff);
  return p;
}

// SYMBOL "cbor_int_size"
inline
casadi_int casadi_cbor_int_size(casadi_int v) {
  return casadi_cbor_head_size(v<0 ? -1-v : v);
}

// SYMBOL "cbor_write_int"
inline
unsigned char* casadi_cbor_write_int(unsigned char* p, casadi_int v) {
  return v<0 ? casadi_cbor_write_head(p, 1, -1-v) : casadi_cbor_write_head(p, 0, v);
}

// SYMBOL "cbor_write_text"
inline
unsigned char* casadi_cbor_write_text(unsigned char* p, const char* s, casadi_int n) {
  casadi_int i;
  p = casadi_cbor_write_head(p, 3, n);
  for (i=0; i<n; ++i) *p++ = (unsigned char) s[i];
  return p;
}

// SYMBOL "cbor_write_array"
inline
unsigned char* casadi_cbor_write_array(unsigned char* p, casadi_int n) {
  // Header only: the n items follow
  return casadi_cbor_write_head(p, 4, n);
}

// SYMBOL "cbor_write_real"
template<typename T1>
unsigned char* casadi_cbor_write_real(unsigned char* p, T1 x) {
  // float64, big-endian: always 9 bytes
  union { double d; unsigned char c[8]; } u;
  casadi_int i;
  int swap = casadi_cbor_little_endian();
  u.d = x;
  *p++ = 0xfb;
  for (i=0; i<8; ++i) *p++ = u.c[swap ? 7-i : i];
  return p;
}

// SYMBOL "cbor_write_bool"
inline
unsigned char* casadi_cbor_write_bool(unsigned char* p, int v) {
  *p++ = v ? 0xf5 : 0xf4;
  return p;
}

// SYMBOL "cbor_write_null"
inline
unsigned char* casadi_cbor_write_null(unsigned char* p) {
  *p++ = 0xf6;
  return p;
}

// SYMBOL "cbor_nullable_size"
inline
casadi_int casadi_cbor_nullable_size(casadi_int v) {
  return v<0 ? 1 : casadi_cbor_head_size(v);
}

// SYMBOL "cbor_write_nullable"
inline
unsigned char* casadi_cbor_write_nullable(unsigned char* p, casadi_int v) {
  // v >= 0 as an unsigned int, or null for none (v < 0)
  return v<0 ? casadi_cbor_write_null(p) : casadi_cbor_write_head(p, 0, v);
}

// SYMBOL "cbor_reserve"
inline
unsigned char* casadi_cbor_reserve(unsigned char* p, casadi_int cap, casadi_int* needed,
    casadi_int len, casadi_int* pos) {
  // len bytes for the next items of a CBOR sequence (RFC 8742) in p[0, cap), of which
  // *needed are taken: their offset *pos, and where to write them, or null if they do not fit.
  // The first that do not fit get a break (0xff) instead, which ends the sequence, and leave
  // *needed at cap+1 for good, so it cannot overflow. Concurrent callers must hold a lock
  *pos = *needed;
  if (*pos <= cap) *needed = len <= cap - *pos ? *pos + len : cap + 1;
  if (len <= cap - *pos) return p + *pos;
  if (*pos < cap) p[*pos] = 0xff;
  return 0;
}

// SYMBOL "cbor_read_head"
inline
const unsigned char* casadi_cbor_read_head(const unsigned char* p, const unsigned char* end,
    casadi_int* major, casadi_int* v) {
  // Major type and argument, as casadi_cbor_write_head writes them; null at a break or beyond
  // casadi_int
  casadi_int n, i, a;
  if (p>=end) return 0;
  *major = *p >> 5;
  n = *p++ & 31;
  if (n<24) {
    *v = n;
    return p;
  }
  // Argument in the next 1, 2, 4 or 8 bytes; no indefinite lengths
  if (n>27) return 0;
  n = 1 << (n-24);
  if (end-p<n) return 0;
  a = 0;
  for (i=0; i<n; ++i) {
    if (a >> (8*sizeof(casadi_int)-9)) return 0;
    a = (a << 8) | *p++;
  }
  *v = a;
  return p;
}

// SYMBOL "cbor_read_int"
inline
const unsigned char* casadi_cbor_read_int(const unsigned char* p, const unsigned char* end,
    casadi_int* v) {
  casadi_int major, a;
  p = casadi_cbor_read_head(p, end, &major, &a);
  if (!p || major>1) return 0;
  *v = major==0 ? a : -1-a;
  return p;
}

// SYMBOL "cbor_read_text"
inline
const unsigned char* casadi_cbor_read_text(const unsigned char* p, const unsigned char* end,
    const unsigned char** s, casadi_int* n) {
  // The n bytes at *s, not null-terminated
  casadi_int major, a;
  p = casadi_cbor_read_head(p, end, &major, &a);
  if (!p || major!=3 || a>end-p) return 0;
  *s = p;
  *n = a;
  return p + a;
}

// SYMBOL "cbor_read_array"
inline
const unsigned char* casadi_cbor_read_array(const unsigned char* p, const unsigned char* end,
    casadi_int* n) {
  // Header only: the n items follow
  casadi_int major, a;
  p = casadi_cbor_read_head(p, end, &major, &a);
  if (!p || major!=4) return 0;
  *n = a;
  return p;
}

// SYMBOL "cbor_read_real"
inline
const unsigned char* casadi_cbor_read_real(const unsigned char* p, const unsigned char* end,
    double* x) {
  // float64, big-endian, as casadi_cbor_write_real writes it
  union { double d; unsigned char c[8]; } u;
  casadi_int i;
  int swap;
  if (end-p<9 || *p!=0xfb) return 0;
  swap = casadi_cbor_little_endian();
  for (i=0; i<8; ++i) u.c[swap ? 7-i : i] = p[1+i];
  *x = u.d;
  return p+9;
}

// SYMBOL "cbor_read_bool"
inline
const unsigned char* casadi_cbor_read_bool(const unsigned char* p, const unsigned char* end,
    int* v) {
  if (p>=end || (*p!=0xf4 && *p!=0xf5)) return 0;
  *v = *p==0xf5;
  return p+1;
}

// SYMBOL "cbor_read_null"
inline
const unsigned char* casadi_cbor_read_null(const unsigned char* p, const unsigned char* end) {
  if (p>=end || *p!=0xf6) return 0;
  return p+1;
}

// SYMBOL "cbor_read_nullable"
inline
const unsigned char* casadi_cbor_read_nullable(const unsigned char* p, const unsigned char* end,
    casadi_int* v) {
  // As casadi_cbor_write_nullable writes it: -1 for null
  const unsigned char* q = casadi_cbor_read_null(p, end);
  if (q) {
    *v = -1;
    return q;
  }
  p = casadi_cbor_read_int(p, end, v);
  return p && *v>=0 ? p : 0;
}

// SYMBOL "cbor_skip"
inline
const unsigned char* casadi_cbor_skip(const unsigned char* p, const unsigned char* end) {
  // Past the item at p, of any kind casadi_cbor.hpp writes; null if it is none of them
  casadi_int major, a, i;
  if (p<end && *p==0xfb) return end-p<9 ? 0 : p+9;
  if (p<end && *p>=0xf4 && *p<=0xf6) return p+1;
  p = casadi_cbor_read_head(p, end, &major, &a);
  if (!p) return 0;
  if (major==0 || major==1) return p;
  if (major==3) return a>end-p ? 0 : p+a;
  if (major!=4) return 0;
  for (i=0; i<a && p; ++i) p = casadi_cbor_skip(p, end);
  return p;
}
