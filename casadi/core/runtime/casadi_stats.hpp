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

// Stats stream: CBOR sequence (RFC 8742) of records; nodes are referenced by stream offset
//   BEGIN_CALL       [1, parent|null, id, mem|null]   (node)
//   END_CALL         [2, call, flag]
//   SET              [3, node, key, value]
//   APPEND           [4, node, key, value]
//   BEGIN_SECTION    [5, parent, name]                (node)
//   END_SCOPE        [6, node]                        (section or iteration)
//   BEGIN_ITERATION  [7, parent, index]               (node)
//   DECLARE_FIELDS   [8, call, [names]]
//   SET_FIELD        [9, iteration, field, value]

// FILTER-MACROS OFF
enum casadi_stats_kind {
  CASADI_STATS_BEGIN_CALL = 1,
  CASADI_STATS_END_CALL = 2,
  CASADI_STATS_SET = 3,
  CASADI_STATS_APPEND = 4,
  CASADI_STATS_BEGIN_SECTION = 5,
  CASADI_STATS_END_SCOPE = 6,
  CASADI_STATS_BEGIN_ITERATION = 7,
  CASADI_STATS_DECLARE_FIELDS = 8,
  CASADI_STATS_SET_FIELD = 9
};
// FILTER-MACROS ON

// SYMBOL "stats_strlen"
inline
casadi_int casadi_stats_strlen(const char* s) {
  casadi_int n = 0;
  while (s[n]) n++;
  return n;
}

// SYMBOL "stats_reserve_record"
inline
unsigned char* casadi_stats_reserve_record(struct casadi_stats_sink* s, casadi_int n,
    enum casadi_stats_kind kind, casadi_int ref, casadi_int len, casadi_int* pos) {
  // Write [n, kind, ref, ...]; returns where the remaining len bytes go
  unsigned char* p = casadi_stats_reserve(s, 2 + casadi_cbor_nullable_size(ref) + len, pos);
  if (!p) return 0;
  p = casadi_cbor_write_array(p, n);
  p = casadi_cbor_write_int(p, kind);
  return casadi_cbor_write_nullable(p, ref);
}

// SYMBOL "stats_read_head"
inline
const unsigned char* casadi_stats_read_head(const unsigned char* p, const unsigned char* end,
    casadi_int* n, casadi_int* kind, casadi_int* ref) {
  // Returns where the other n-2 items start; null at a break or cut
  p = casadi_cbor_read_array(p, end, n);
  if (!p || *n<2) return 0;
  p = casadi_cbor_read_int(p, end, kind);
  if (!p) return 0;
  return casadi_cbor_read_nullable(p, end, ref);
}

// SYMBOL "stats_begin_call"
inline
casadi_int casadi_stats_begin_call(struct casadi_stats_sink* s, const char* id, casadi_int mem,
    casadi_int parent) {
  casadi_int pos, n;
  unsigned char* p;
  if (!s) return -1;
  n = casadi_stats_strlen(id);
  p = casadi_stats_reserve_record(s, 4, CASADI_STATS_BEGIN_CALL, parent,
    casadi_cbor_head_size(n) + n + casadi_cbor_nullable_size(mem), &pos);
  if (p) {
    p = casadi_cbor_write_text(p, id, n);
    casadi_cbor_write_nullable(p, mem);
    casadi_stats_commit(s, pos);
  }
  return pos;
}

// SYMBOL "stats_end_call"
inline
void casadi_stats_end_call(struct casadi_stats_sink* s, casadi_int call, casadi_int flag) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_record(s, 3, CASADI_STATS_END_CALL, call,
    casadi_cbor_int_size(flag), &pos);
  if (p) {
    casadi_cbor_write_int(p, flag);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_begin_section"
inline
casadi_int casadi_stats_begin_section(struct casadi_stats_sink* s, casadi_int parent,
    const char* name) {
  casadi_int pos, n;
  unsigned char* p;
  if (!s) return -1;
  n = casadi_stats_strlen(name);
  p = casadi_stats_reserve_record(s, 3, CASADI_STATS_BEGIN_SECTION, parent,
    casadi_cbor_head_size(n) + n, &pos);
  if (p) {
    casadi_cbor_write_text(p, name, n);
    casadi_stats_commit(s, pos);
  }
  return pos;
}

// SYMBOL "stats_end_scope"
inline
void casadi_stats_end_scope(struct casadi_stats_sink* s, casadi_int node) {
  casadi_int pos;
  if (!s || node<0) return;
  if (casadi_stats_reserve_record(s, 2, CASADI_STATS_END_SCOPE, node, 0, &pos)) {
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_begin_iteration"
inline
casadi_int casadi_stats_begin_iteration(struct casadi_stats_sink* s, casadi_int parent,
    casadi_int index) {
  // Close with casadi_stats_end_scope
  casadi_int pos;
  unsigned char* p;
  if (!s) return -1;
  p = casadi_stats_reserve_record(s, 3, CASADI_STATS_BEGIN_ITERATION, parent,
    casadi_cbor_int_size(index), &pos);
  if (p) {
    casadi_cbor_write_int(p, index);
    casadi_stats_commit(s, pos);
  }
  return pos;
}

// SYMBOL "stats_reserve_value"
inline
unsigned char* casadi_stats_reserve_value(struct casadi_stats_sink* s, casadi_int ref,
    enum casadi_stats_kind op, const char* key, casadi_int nval, casadi_int* pos) {
  // key null for SET_FIELD
  casadi_int nkey;
  unsigned char* p;
  nkey = key ? casadi_stats_strlen(key) : 0;
  p = casadi_stats_reserve_record(s, 4, op, ref,
    key ? casadi_cbor_head_size(nkey) + nkey + nval : nval, pos);
  if (!p || !key) return p;
  return casadi_cbor_write_text(p, key, nkey);
}

// SYMBOL "stats_reserve_field"
inline
unsigned char* casadi_stats_reserve_field(struct casadi_stats_sink* s, casadi_int it,
    casadi_int field, casadi_int nval, casadi_int* pos) {
  unsigned char* p = casadi_stats_reserve_value(s, it, CASADI_STATS_SET_FIELD, 0,
    casadi_cbor_int_size(field) + nval, pos);
  return p ? casadi_cbor_write_int(p, field) : 0;
}

// SYMBOL "stats_put_int"
inline
void casadi_stats_put_int(struct casadi_stats_sink* s, casadi_int ref,
    enum casadi_stats_kind op, const char* key, casadi_int v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_value(s, ref, op, key, casadi_cbor_int_size(v), &pos);
  if (p) {
    casadi_cbor_write_int(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_put_bool"
inline
void casadi_stats_put_bool(struct casadi_stats_sink* s, casadi_int ref,
    enum casadi_stats_kind op, const char* key, int v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_value(s, ref, op, key, 1, &pos);
  if (p) {
    casadi_cbor_write_bool(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_put_real"
template<typename T1>
void casadi_stats_put_real(struct casadi_stats_sink* s, casadi_int ref,
    enum casadi_stats_kind op, const char* key, T1 v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_value(s, ref, op, key, 9, &pos);
  if (p) {
    casadi_cbor_write_real(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_put_text"
inline
void casadi_stats_put_text(struct casadi_stats_sink* s, casadi_int ref,
    enum casadi_stats_kind op, const char* key, const char* v) {
  casadi_int n, pos;
  unsigned char* p;
  if (!s) return;
  n = casadi_stats_strlen(v);
  p = casadi_stats_reserve_value(s, ref, op, key, casadi_cbor_head_size(n) + n, &pos);
  if (p) {
    casadi_cbor_write_text(p, v, n);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_put_reals"
template<typename T1>
void casadi_stats_put_reals(struct casadi_stats_sink* s, casadi_int ref,
    enum casadi_stats_kind op, const char* key, const T1* v, casadi_int n) {
  casadi_int i, pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_value(s, ref, op, key, casadi_cbor_head_size(n) + 9*n, &pos);
  if (p) {
    p = casadi_cbor_write_array(p, n);
    for (i=0; i<n; ++i) p = casadi_cbor_write_real(p, v[i]);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_texts_size"
inline
casadi_int casadi_stats_texts_size(const char** v, casadi_int n) {
  casadi_int i, len;
  len = casadi_cbor_head_size(n);
  for (i=0; i<n; ++i) {
    len += casadi_cbor_head_size(casadi_stats_strlen(v[i])) + casadi_stats_strlen(v[i]);
  }
  return len;
}

// SYMBOL "stats_write_texts"
inline
unsigned char* casadi_stats_write_texts(unsigned char* p, const char** v, casadi_int n) {
  casadi_int i;
  p = casadi_cbor_write_array(p, n);
  for (i=0; i<n; ++i) p = casadi_cbor_write_text(p, v[i], casadi_stats_strlen(v[i]));
  return p;
}

// SYMBOL "stats_put_texts"
inline
void casadi_stats_put_texts(struct casadi_stats_sink* s, casadi_int ref,
    enum casadi_stats_kind op, const char* key, const char** v, casadi_int n) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_value(s, ref, op, key, casadi_stats_texts_size(v, n), &pos);
  if (p) {
    casadi_stats_write_texts(p, v, n);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_declare_fields"
inline
void casadi_stats_declare_fields(struct casadi_stats_sink* s, casadi_int call,
    const char** names, casadi_int n) {
  // Once, before the first iteration
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_record(s, 3, CASADI_STATS_DECLARE_FIELDS, call,
    casadi_stats_texts_size(names, n), &pos);
  if (p) {
    casadi_stats_write_texts(p, names, n);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_field_int"
inline
void casadi_stats_set_field_int(struct casadi_stats_sink* s, casadi_int it, casadi_int field,
    casadi_int v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_field(s, it, field, casadi_cbor_int_size(v), &pos);
  if (p) {
    casadi_cbor_write_int(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_field_bool"
inline
void casadi_stats_set_field_bool(struct casadi_stats_sink* s, casadi_int it, casadi_int field,
    int v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_field(s, it, field, 1, &pos);
  if (p) {
    casadi_cbor_write_bool(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_field_real"
template<typename T1>
void casadi_stats_set_field_real(struct casadi_stats_sink* s, casadi_int it, casadi_int field,
    T1 v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_field(s, it, field, 9, &pos);
  if (p) {
    casadi_cbor_write_real(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_field_text"
inline
void casadi_stats_set_field_text(struct casadi_stats_sink* s, casadi_int it, casadi_int field,
    const char* v) {
  casadi_int n, pos;
  unsigned char* p;
  if (!s) return;
  n = casadi_stats_strlen(v);
  p = casadi_stats_reserve_field(s, it, field, casadi_cbor_head_size(n) + n, &pos);
  if (p) {
    casadi_cbor_write_text(p, v, n);
    casadi_stats_commit(s, pos);
  }
}

// Reading, as StatsRecorder::get_stat

// SYMBOL "stats_equals"
inline
int casadi_stats_equals(const unsigned char* p, casadi_int n, const char* s) {
  casadi_int i;
  for (i=0; i<n; ++i) {
    if (!s[i] || p[i]!=(unsigned char) s[i]) return 0;
  }
  return !s[n];
}

// SYMBOL "stats_id_is"
inline
int casadi_stats_id_is(const unsigned char* id, casadi_int n, const char* fname) {
  // id: "<symbol>:<name>"
  casadi_int k;
  for (k=0; k<n && id[k]!=':'; ++k) {}
  k = k<n ? k+1 : 0;
  return casadi_stats_equals(id+k, n-k, fname);
}

// SYMBOL "stats_find"
inline
const unsigned char* casadi_stats_find(const unsigned char* s, casadi_int n, const char* fname,
    const char* key) {
  // Value of key in the last call of fname, or null
  const unsigned char *p, *q, *t, *val, *v, *end;
  casadi_int len, kind, ref, nt, i, call;
  int named;
  end = s + n;
  v = t = 0;
  nt = 0;
  call = -1;
  for (p=s; p<end; p=q) {
    q = casadi_stats_read_head(p, end, &len, &kind, &ref);
    if (!q) break;
    named = kind==CASADI_STATS_BEGIN_CALL || (kind==CASADI_STATS_SET && call>=0 && ref==call);
    if (named) q = casadi_cbor_read_text(q, end, &t, &nt);
    val = q;
    for (i=named ? 3 : 2; i<len && q; ++i) q = casadi_cbor_skip(q, end);
    if (!q) break;
    if (kind==CASADI_STATS_BEGIN_CALL) {
      if (casadi_stats_id_is(t, nt, fname)) {
        call = p - s;
        v = 0;
      }
    } else if (named && casadi_stats_equals(t, nt, key)) {
      v = val;
    }
  }
  return v;
}

// SYMBOL "stats_get_int"
inline
int casadi_stats_get_int(const unsigned char* s, casadi_int n, const char* fname,
    const char* key, casadi_int* v) {
  const unsigned char* p;
  int b;
  p = casadi_stats_find(s, n, fname, key);
  if (!p) return 1;
  if (casadi_cbor_read_int(p, s+n, v)) return 0;
  if (!casadi_cbor_read_bool(p, s+n, &b)) return 1;
  *v = b;
  return 0;
}

// C-REPLACE "static_cast<double>" "(double) "
// SYMBOL "stats_get_real"
inline
int casadi_stats_get_real(const unsigned char* s, casadi_int n, const char* fname,
    const char* key, double* v) {
  const unsigned char* p;
  casadi_int i;
  int b;
  p = casadi_stats_find(s, n, fname, key);
  if (!p) return 1;
  if (casadi_cbor_read_real(p, s+n, v)) return 0;
  if (casadi_cbor_read_int(p, s+n, &i)) {
    *v = static_cast<double>(i);
    return 0;
  }
  if (!casadi_cbor_read_bool(p, s+n, &b)) return 1;
  *v = b;
  return 0;
}

// C-REPLACE "static_cast<char>" "(char) "
// SYMBOL "stats_get_text"
inline
int casadi_stats_get_text(const unsigned char* s, casadi_int n, const char* fname,
    const char* key, char* dst, casadi_int cap) {
  // Truncated to cap, null-terminated
  const unsigned char *p, *t;
  casadi_int i, nt;
  p = casadi_stats_find(s, n, fname, key);
  if (!p || !casadi_cbor_read_text(p, s+n, &t, &nt)) return 1;
  if (cap<=0) return 0;
  for (i=0; i<nt && i<cap-1; ++i) dst[i] = static_cast<char>(t[i]);
  dst[i] = 0;
  return 0;
}

// SYMBOL "get_stat_int"
// EXPORT
inline
int casadi_get_stat_int(const struct casadi_stats_sink* s, const char* fname, const char* key,
    casadi_int* v) {
  // 0 if found as int or bool
  casadi_int n;
  const unsigned char* p = casadi_stats_data(s, &n);
  return casadi_stats_get_int(p, n, fname, key, v);
}

// SYMBOL "get_stat_real"
// EXPORT
inline
int casadi_get_stat_real(const struct casadi_stats_sink* s, const char* fname, const char* key,
    double* v) {
  // 0 if found as real, int or bool
  casadi_int n;
  const unsigned char* p = casadi_stats_data(s, &n);
  return casadi_stats_get_real(p, n, fname, key, v);
}

// SYMBOL "get_stat_text"
// EXPORT
inline
int casadi_get_stat_text(const struct casadi_stats_sink* s, const char* fname, const char* key,
    char* dst, casadi_int cap) {
  // 0 if found as text; truncated to cap, null-terminated
  casadi_int n;
  const unsigned char* p = casadi_stats_data(s, &n);
  return casadi_stats_get_text(p, n, fname, key, dst, cap);
}
