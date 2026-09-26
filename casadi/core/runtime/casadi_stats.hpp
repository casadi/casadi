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
//   BEGIN_CALL       [1, parent|null, id, mem|null, end]   (node)
//   END_CALL         [2, call, flag]
//   SET              [3, node, key, value]
//   BEGIN_SECTION    [4, parent, name, end]                (node)
//   END_SCOPE        [5, node]                             (section or iteration)
//   BEGIN_ITERATION  [6, parent, index, end]               (node)
//   DECLARE_FIELDS   [7, call, [names]]
//   SET_FIELD        [8, iteration, field, value]
// Nothing refers to a node after its END record. end: offset of that record, as a 4-byte
// uint (0x1a); 0 until casadi_stats_index patches it in place

// FILTER-MACROS OFF
enum casadi_stats_kind {
  CASADI_STATS_BEGIN_CALL = 1,
  CASADI_STATS_END_CALL = 2,
  CASADI_STATS_SET = 3,
  CASADI_STATS_BEGIN_SECTION = 4,
  CASADI_STATS_END_SCOPE = 5,
  CASADI_STATS_BEGIN_ITERATION = 6,
  CASADI_STATS_DECLARE_FIELDS = 7,
  CASADI_STATS_SET_FIELD = 8
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
  // Returns where the other n-2 items start; null at a break or cut (outputs zero)
  *n = *kind = *ref = 0;
  p = casadi_cbor_read_array(p, end, n);
  if (!p || *n<2) return 0;
  p = casadi_cbor_read_int(p, end, kind);
  if (!p) return 0;
  return casadi_cbor_read_nullable(p, end, ref);
}

// SYMBOL "stats_write_end"
inline
unsigned char* casadi_stats_write_end(unsigned char* p, casadi_int e) {
  // Fixed width, to patch in place; shifts by 8 only: safe for any casadi_int width
  casadi_int i;
  *p++ = 0x1a;
  for (i=3; i>=0; --i) {
    p[i] = (unsigned char) (e & 0xff);
    e >>= 8;
  }
  return p + 4;
}

// SYMBOL "stats_begin_call"
inline
casadi_int casadi_stats_begin_call(struct casadi_stats_sink* s, const char* id, casadi_int mem,
    casadi_int parent) {
  casadi_int pos, n;
  unsigned char* p;
  if (!s) return -1;
  n = casadi_stats_strlen(id);
  p = casadi_stats_reserve_record(s, 5, CASADI_STATS_BEGIN_CALL, parent,
    casadi_cbor_head_size(n) + n + casadi_cbor_nullable_size(mem) + 5, &pos);
  if (p) {
    p = casadi_cbor_write_text(p, id, n);
    p = casadi_cbor_write_nullable(p, mem);
    casadi_stats_write_end(p, 0);
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
  p = casadi_stats_reserve_record(s, 4, CASADI_STATS_BEGIN_SECTION, parent,
    casadi_cbor_head_size(n) + n + 5, &pos);
  if (p) {
    p = casadi_cbor_write_text(p, name, n);
    casadi_stats_write_end(p, 0);
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
  p = casadi_stats_reserve_record(s, 4, CASADI_STATS_BEGIN_ITERATION, parent,
    casadi_cbor_int_size(index) + 5, &pos);
  if (p) {
    p = casadi_cbor_write_int(p, index);
    casadi_stats_write_end(p, 0);
    casadi_stats_commit(s, pos);
  }
  return pos;
}

// SYMBOL "stats_reserve_set"
inline
unsigned char* casadi_stats_reserve_set(struct casadi_stats_sink* s, casadi_int node,
    const char* key, casadi_int nval, casadi_int* pos) {
  casadi_int nkey;
  unsigned char* p;
  nkey = casadi_stats_strlen(key);
  p = casadi_stats_reserve_record(s, 4, CASADI_STATS_SET, node,
    casadi_cbor_head_size(nkey) + nkey + nval, pos);
  return p ? casadi_cbor_write_text(p, key, nkey) : 0;
}

// SYMBOL "stats_reserve_field"
inline
unsigned char* casadi_stats_reserve_field(struct casadi_stats_sink* s, casadi_int it,
    casadi_int field, casadi_int nval, casadi_int* pos) {
  unsigned char* p = casadi_stats_reserve_record(s, 4, CASADI_STATS_SET_FIELD, it,
    casadi_cbor_int_size(field) + nval, pos);
  return p ? casadi_cbor_write_int(p, field) : 0;
}

// SYMBOL "stats_set_int"
inline
void casadi_stats_set_int(struct casadi_stats_sink* s, casadi_int node,
    const char* key, casadi_int v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_set(s, node, key, casadi_cbor_int_size(v), &pos);
  if (p) {
    casadi_cbor_write_int(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_bool"
inline
void casadi_stats_set_bool(struct casadi_stats_sink* s, casadi_int node,
    const char* key, int v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_set(s, node, key, 1, &pos);
  if (p) {
    casadi_cbor_write_bool(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_real"
template<typename T1>
void casadi_stats_set_real(struct casadi_stats_sink* s, casadi_int node,
    const char* key, T1 v) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_set(s, node, key, 9, &pos);
  if (p) {
    casadi_cbor_write_real(p, v);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_text"
inline
void casadi_stats_set_text(struct casadi_stats_sink* s, casadi_int node,
    const char* key, const char* v) {
  casadi_int n, pos;
  unsigned char* p;
  if (!s) return;
  n = casadi_stats_strlen(v);
  p = casadi_stats_reserve_set(s, node, key, casadi_cbor_head_size(n) + n, &pos);
  if (p) {
    casadi_cbor_write_text(p, v, n);
    casadi_stats_commit(s, pos);
  }
}

// SYMBOL "stats_set_reals"
template<typename T1>
void casadi_stats_set_reals(struct casadi_stats_sink* s, casadi_int node,
    const char* key, const T1* v, casadi_int n) {
  casadi_int i, pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_set(s, node, key, casadi_cbor_head_size(n) + 9*n, &pos);
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

// SYMBOL "stats_set_texts"
inline
void casadi_stats_set_texts(struct casadi_stats_sink* s, casadi_int node,
    const char* key, const char** v, casadi_int n) {
  casadi_int pos;
  unsigned char* p;
  if (!s) return;
  p = casadi_stats_reserve_set(s, node, key, casadi_stats_texts_size(v, n), &pos);
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

// Reading: a node is the stream offset of its BEGIN record, -1 the root

// SYMBOL "stats_is"
inline
int casadi_stats_is(const char* a, const char* b) {
  while (*a && *a==*b) {
    a++;
    b++;
  }
  return *a==*b;
}

// SYMBOL "stats_equals"
inline
int casadi_stats_equals(const unsigned char* p, casadi_int n, const char* s) {
  casadi_int i;
  for (i=0; i<n; ++i) {
    if (!s[i] || p[i]!=(unsigned char) s[i]) return 0;
  }
  return !s[n];
}

// SYMBOL "stats_name"
inline
const unsigned char* casadi_stats_name(const unsigned char* id, casadi_int* n) {
  // id: "<symbol>:<name>"; the name follows the first ':', if any
  casadi_int k;
  for (k=0; k<*n && id[k]!=':'; ++k) {}
  if (k==*n) return id;
  *n -= k + 1;
  return id + k + 1;
}

// SYMBOL "stats_matches"
inline
int casadi_stats_matches(const unsigned char* id, casadi_int n, const char* pattern) {
  // A pattern with ':' is the full id, otherwise the name; null or empty: any
  casadi_int k;
  if (!pattern || !*pattern) return 1;
  for (k=0; pattern[k]; ++k) {
    if (pattern[k]==':') return casadi_stats_equals(id, n, pattern);
  }
  id = casadi_stats_name(id, &n);
  return casadi_stats_equals(id, n, pattern);
}

// SYMBOL "stats_ends"
inline
int casadi_stats_ends(casadi_int kind, casadi_int ref, casadi_int node) {
  // Does the record end node? Nothing refers to a node after its end
  return ref==node && (kind==CASADI_STATS_END_CALL || kind==CASADI_STATS_END_SCOPE);
}

// SYMBOL "stats_record"
inline
const unsigned char* casadi_stats_record(const unsigned char* p, const unsigned char* end,
    casadi_int* kind, casadi_int* ref, const unsigned char** items) {
  // End of the record at p, null at a break or cut; *items: after the reference
  casadi_int n, i;
  p = casadi_stats_read_head(p, end, &n, kind, ref);
  *items = p;
  for (i=2; i<n && p; ++i) p = casadi_cbor_skip(p, end);
  return p;
}

// SYMBOL "stats_slot"
inline
casadi_int casadi_stats_slot(const unsigned char* s, casadi_int n, casadi_int node) {
  // Offset of the end field of node, -1 if not a node
  const unsigned char *p, *end;
  casadi_int len, kind, ref, i;
  if (node<0 || node>=n) return -1;
  end = s + n;
  p = casadi_stats_read_head(s+node, end, &len, &kind, &ref);
  if (!p || (kind!=CASADI_STATS_BEGIN_CALL && kind!=CASADI_STATS_BEGIN_SECTION
      && kind!=CASADI_STATS_BEGIN_ITERATION)) return -1;
  for (i=3; i<len && p; ++i) p = casadi_cbor_skip(p, end);
  return p && end-p>=5 && *p==0x1a ? p - s : -1;
}

// SYMBOL "stats_read_end"
inline
casadi_int casadi_stats_read_end(const unsigned char* s, casadi_int n, casadi_int k) {
  // Value of the end field at k if below n, else -1
  casadi_int e, i;
  e = 0;
  for (i=1; i<=4; ++i) {
    if (e > (n >> 8)) return -1;
    e = (e << 8) | s[k+i];
  }
  return e<n ? e : -1;
}

// SYMBOL "stats_end"
inline
casadi_int casadi_stats_end(const unsigned char* s, casadi_int n, casadi_int node) {
  // Offset of the END record of node once indexed, else 0
  casadi_int k, e, len, kind, ref;
  k = casadi_stats_slot(s, n, node);
  if (k<0) return 0;
  e = casadi_stats_read_end(s, n, k);
  if (e<=node || !casadi_stats_read_head(s+e, s+n, &len, &kind, &ref)) return 0;
  return casadi_stats_ends(kind, ref, node) ? e : 0;
}

// SYMBOL "stats_over"
inline
const unsigned char* casadi_stats_over(const unsigned char* s, casadi_int n, casadi_int pos,
    casadi_int kind, const unsigned char* q) {
  // Where to continue after the record at pos, ending at q: at the END record of an indexed
  // node, skipping its subtree
  casadi_int e;
  if (kind!=CASADI_STATS_BEGIN_CALL && kind!=CASADI_STATS_BEGIN_SECTION
      && kind!=CASADI_STATS_BEGIN_ITERATION) return q;
  e = casadi_stats_end(s, n, pos);
  return e ? s + e : q;
}

// SYMBOL "stats_start"
inline
const unsigned char* casadi_stats_start(const unsigned char* s, casadi_int n, casadi_int from,
    int skip) {
  // Past the record at from, or over its subtree (skip); the stream start for -1
  casadi_int kind, ref;
  const unsigned char *items, *q;
  if (from<0) return s;
  if (from>=n) return 0;
  q = casadi_stats_record(s+from, s+n, &kind, &ref, &items);
  return q && skip ? casadi_stats_over(s, n, from, kind, q) : q;
}

// SYMBOL "stats_below"
inline
int casadi_stats_below(const unsigned char* s, casadi_int n, casadi_int node, casadi_int scope) {
  // Is node strictly below scope? References point backwards, which ends the walk
  casadi_int len, kind, ref;
  while (node>=0) {
    if (!casadi_stats_read_head(s+node, s+n, &len, &kind, &ref) || ref>=node) return 0;
    if (ref==scope) return 1;
    node = ref;
  }
  return 0;
}

// SYMBOL "stats_next_named"
inline
casadi_int casadi_stats_next_named(const unsigned char* s, casadi_int n,
    enum casadi_stats_kind begin, casadi_int scope, casadi_int cursor, const char* pattern,
    int direct) {
  // Next call or section after cursor matching pattern: a child of scope (direct), skipping
  // indexed subtrees, or any descendant
  const unsigned char *p, *q, *items, *t, *end;
  casadi_int kind, ref, nt;
  end = s + n;
  for (p=casadi_stats_start(s, n, cursor>=0 ? cursor : scope, direct && cursor>=0); p; p=q) {
    q = casadi_stats_record(p, end, &kind, &ref, &items);
    if (!q || casadi_stats_ends(kind, ref, scope)) break;
    if (kind==begin && casadi_cbor_read_text(items, end, &t, &nt)
        && casadi_stats_matches(t, nt, pattern)
        && (direct ? ref==scope : casadi_stats_below(s, n, p-s, scope))) return p - s;
    if (direct) q = casadi_stats_over(s, n, p-s, kind, q);
  }
  return -1;
}

// SYMBOL "stats_find_function"
// EXPORT
inline
casadi_int casadi_stats_find_function(const struct casadi_stats_sink* s, casadi_int scope,
    casadi_int cursor, const char* pattern) {
  // Next call after cursor (-1: first) at any depth below scope (-1: root); pattern with ':'
  // matches the full id, otherwise the function name, empty or null any; -1 if none
  casadi_int n;
  const unsigned char* p = casadi_stats_data(s, &n);
  return casadi_stats_next_named(p, n, CASADI_STATS_BEGIN_CALL, scope, cursor, pattern, 0);
}

// SYMBOL "stats_select_function"
// EXPORT
inline
casadi_int casadi_stats_select_function(const struct casadi_stats_sink* s, casadi_int parent,
    casadi_int cursor, const char* pattern) {
  // As find_function, among the direct children of parent
  casadi_int n;
  const unsigned char* p = casadi_stats_data(s, &n);
  return casadi_stats_next_named(p, n, CASADI_STATS_BEGIN_CALL, parent, cursor, pattern, 1);
}

// SYMBOL "stats_select_iteration"
// EXPORT
inline
casadi_int casadi_stats_select_iteration(const struct casadi_stats_sink* s, casadi_int parent,
    casadi_int cursor, casadi_int index) {
  // Next iteration of parent after cursor (-1: first) with index (-1: any), or -1
  const unsigned char *b, *p, *q, *items, *end;
  casadi_int n, kind, ref, i;
  b = casadi_stats_data(s, &n);
  end = b + n;
  for (p=casadi_stats_start(b, n, cursor>=0 ? cursor : parent, cursor>=0); p; p=q) {
    q = casadi_stats_record(p, end, &kind, &ref, &items);
    if (!q || casadi_stats_ends(kind, ref, parent)) break;
    if (kind==CASADI_STATS_BEGIN_ITERATION && ref==parent
        && casadi_cbor_read_int(items, end, &i) && (index<0 || i==index)) return p - b;
    q = casadi_stats_over(b, n, p-b, kind, q);
  }
  return -1;
}

// SYMBOL "stats_select_last_iteration"
// EXPORT
inline
casadi_int casadi_stats_select_last_iteration(const struct casadi_stats_sink* s,
    casadi_int parent) {
  // Iteration of parent recorded last, or -1
  casadi_int r = -1, it = -1;
  while ((it = casadi_stats_select_iteration(s, parent, it, -1))>=0) r = it;
  return r;
}

// SYMBOL "stats_select_section"
// EXPORT
inline
casadi_int casadi_stats_select_section(const struct casadi_stats_sink* s, casadi_int parent,
    casadi_int cursor, const char* name) {
  // Next section of parent after cursor (-1: first) named name (null or empty: any), or -1
  casadi_int n;
  const unsigned char* p = casadi_stats_data(s, &n);
  return casadi_stats_next_named(p, n, CASADI_STATS_BEGIN_SECTION, parent, cursor, name, 1);
}

// SYMBOL "stats_field"
inline
casadi_int casadi_stats_field(const unsigned char* s, casadi_int n, casadi_int call,
    const char* key) {
  // Index of key among the iteration fields call declares, or -1
  const unsigned char *p, *q, *items, *t, *end;
  casadi_int kind, ref, m, j, nt;
  end = s + n;
  for (p=casadi_stats_start(s, n, call, 0); p; p=q) {
    q = casadi_stats_record(p, end, &kind, &ref, &items);
    if (!q || casadi_stats_ends(kind, ref, call)) break;
    if (kind==CASADI_STATS_DECLARE_FIELDS && ref==call) {
      items = casadi_cbor_read_array(items, end, &m);
      for (j=0; j<m && items; ++j) {
        items = casadi_cbor_read_text(items, end, &t, &nt);
        if (items && casadi_stats_equals(t, nt, key)) return j;
      }
      break;
    }
    q = casadi_stats_over(s, n, p-s, kind, q);
  }
  return -1;
}

// SYMBOL "stats_lookup"
inline
const unsigned char* casadi_stats_lookup(const unsigned char* s, casadi_int n, casadi_int node,
    const char* key) {
  // Item holding key of node, or null. Keys: the last SET, a field of an iteration,
  // and the reserved id, name, mem and flag (call), name (section) and index (iteration);
  // name is the whole id
  const unsigned char *p, *q, *items, *t, *v, *end;
  casadi_int kind, ref, nt, f, field, e;
  int flag;
  end = s + n;
  if (node<0 || node>=n) return 0;
  q = casadi_stats_record(s+node, end, &kind, &ref, &items);
  if (!q) return 0;
  if (kind==CASADI_STATS_BEGIN_CALL && casadi_stats_is(key, "id")) return items;
  if ((kind==CASADI_STATS_BEGIN_CALL || kind==CASADI_STATS_BEGIN_SECTION)
      && casadi_stats_is(key, "name")) return items;
  if (kind==CASADI_STATS_BEGIN_CALL && casadi_stats_is(key, "mem")) {
    return casadi_cbor_skip(items, end);
  }
  if (kind==CASADI_STATS_BEGIN_ITERATION && casadi_stats_is(key, "index")) return items;
  flag = kind==CASADI_STATS_BEGIN_CALL && casadi_stats_is(key, "flag");
  field = kind==CASADI_STATS_BEGIN_ITERATION ? casadi_stats_field(s, n, ref, key) : -1;
  // Indexed: the flag is in the END record
  e = flag ? casadi_stats_end(s, n, node) : 0;
  if (e) return casadi_stats_record(s+e, end, &kind, &ref, &items) ? items : 0;
  v = 0;
  for (p=q; p; p=q) {
    q = casadi_stats_record(p, end, &kind, &ref, &items);
    if (!q) break;
    if (ref==node) {
      if (flag && kind==CASADI_STATS_END_CALL) return items;
      if (casadi_stats_ends(kind, ref, node)) break;
      if (!flag && kind==CASADI_STATS_SET) {
        items = casadi_cbor_read_text(items, end, &t, &nt);
        if (items && casadi_stats_equals(t, nt, key)) v = items;
      } else if (!flag && kind==CASADI_STATS_SET_FIELD && field>=0) {
        items = casadi_cbor_read_int(items, end, &f);
        if (items && f==field) v = items;
      }
    }
    q = casadi_stats_over(s, n, p-s, kind, q);
  }
  return v;
}

// SYMBOL "stats_contiguous"
inline
int casadi_stats_contiguous(const unsigned char* s, casadi_int n, casadi_int node,
    casadi_int e) {
  // Do the records from node up to e all belong to its subtree? Indexed ones are skipped
  const unsigned char *p, *q, *items, *end;
  casadi_int kind, ref;
  end = s + n;
  for (p=casadi_stats_start(s, n, node, 0); p && p<s+e; p=q) {
    q = casadi_stats_record(p, end, &kind, &ref, &items);
    if (!q || (ref!=node && !casadi_stats_below(s, n, ref, node))) return 0;
    q = casadi_stats_over(s, n, p-s, kind, q);
  }
  return 1;
}

// SYMBOL "stats_index"
// EXPORT
inline
void casadi_stats_index(struct casadi_stats_sink* s) {
  // After the last write: store where each node ends, which speeds up queries. Nodes whose
  // records interleave with others (parallel evaluation) are not indexed
  unsigned char* b;
  const unsigned char *p, *q, *items, *end;
  casadi_int n, kind, ref, k, e;
  casadi_stats_data(s, &n);
  b = s->p;
  end = b + n;
  for (p=b; p; p=q) {
    q = casadi_stats_record(p, end, &kind, &ref, &items);
    if (!q) break;
    e = p - b;
    if ((k = casadi_stats_slot(b, n, e))>=0) {
      // Children end first: reset what an earlier index stored
      casadi_stats_write_end(b+k, 0);
    } else if ((kind==CASADI_STATS_END_CALL || kind==CASADI_STATS_END_SCOPE) && ref<e
        && ((e >> 16) >> 16)==0 && (k = casadi_stats_slot(b, n, ref))>=0
        && casadi_stats_contiguous(b, n, ref, e)) {
      casadi_stats_write_end(b+k, e);
    }
  }
}

// SYMBOL "stats_copy_record"
inline
unsigned char* casadi_stats_copy_record(const unsigned char* b, casadi_int n,
    const unsigned char* p, const unsigned char* q, unsigned char* w, const unsigned char* lim) {
  // Copy the record at p (ending at q) to w, its reference renumbered through the end field of
  // the node it refers to; null if it does not fit before lim
  const unsigned char* items;
  casadi_int m, kind, ref;
  items = casadi_stats_read_head(p, q, &m, &kind, &ref);
  if (ref>=0) ref = casadi_stats_read_end(b, n, casadi_stats_slot(b, n, ref));
  if (lim - w < 2 + casadi_cbor_nullable_size(ref) + (q - items)) return 0;
  w = casadi_cbor_write_array(w, m);
  w = casadi_cbor_write_int(w, kind);
  w = casadi_cbor_write_nullable(w, ref);
  while (items<q) *w++ = *items++;
  return w;
}

// SYMBOL "stats_deinterleave"
// EXPORT
inline
int casadi_stats_deinterleave(struct casadi_stats_sink* s) {
  // After the last write: reorder depth first, the records of each node contiguous and in
  // stream order, then index. Writes to the unwritten part of the buffer first; returns 1,
  // leaving the stream as it was, if that does not fit or the stream is truncated
  unsigned char *b, *o, *w;
  const unsigned char *p, *q, *items, *end;
  casadi_int n, node, kind, ref, len;
  if (casadi_stats_truncated(s)) return 1;
  casadi_stats_data(s, &n);
  b = s->p;
  end = b + n;
  o = w = b + n;
  // Depth first without a stack: back to the parent past the node's BEGIN record
  node = -1;
  for (p=b; w; ) {
    q = casadi_stats_record(p, end, &kind, &ref, &items);
    if (!q || casadi_stats_ends(kind, ref, node)) {
      if (q) w = casadi_stats_copy_record(b, n, p, q, w, b + s->cap);
      if (node<0 || !w) break;
      p = casadi_stats_start(b, n, node, 0);
      casadi_stats_read_head(b+node, end, &len, &kind, &node);
    } else if (ref==node) {
      if (kind==CASADI_STATS_BEGIN_CALL || kind==CASADI_STATS_BEGIN_SECTION
          || kind==CASADI_STATS_BEGIN_ITERATION) {
        // The end field maps the node to its new offset until indexed
        if ((((w - o) >> 16) >> 16)!=0) break;
        casadi_stats_write_end(b + casadi_stats_slot(b, n, p-b), w - o);
        node = p - b;
      }
      w = casadi_stats_copy_record(b, n, p, q, w, b + s->cap);
      p = q;
    } else {
      p = q;
    }
  }
  if (w && node<0) {
    for (p=o; p<w; ) *b++ = *p++;
    s->needed = w - o;
  }
  casadi_stats_index(s);
  return w && node<0 ? 0 : 1;
}

// C-REPLACE "static_cast<double>" "(double) "
// SYMBOL "stats_read_number"
inline
const unsigned char* casadi_stats_read_number(const unsigned char* p, const unsigned char* end,
    double* v) {
  // A real, int or bool as double
  const unsigned char* q;
  casadi_int i;
  int b;
  if ((q = casadi_cbor_read_real(p, end, v))) return q;
  if ((q = casadi_cbor_read_int(p, end, &i))) {
    *v = static_cast<double>(i);
    return q;
  }
  if ((q = casadi_cbor_read_bool(p, end, &b))) *v = b;
  return q;
}

// SYMBOL "stats_get_int"
// EXPORT
inline
int casadi_stats_get_int(const struct casadi_stats_sink* s, casadi_int node, const char* key,
    casadi_int* v) {
  // 0 if found as int or bool
  const unsigned char *b, *p;
  casadi_int n;
  int t;
  b = casadi_stats_data(s, &n);
  p = casadi_stats_lookup(b, n, node, key);
  if (!p) return 1;
  if (casadi_cbor_read_int(p, b+n, v)) return 0;
  if (!casadi_cbor_read_bool(p, b+n, &t)) return 1;
  *v = t;
  return 0;
}

// SYMBOL "stats_get_real"
// EXPORT
inline
int casadi_stats_get_real(const struct casadi_stats_sink* s, casadi_int node, const char* key,
    double* v) {
  // 0 if found as real, int or bool
  const unsigned char *b, *p;
  casadi_int n;
  b = casadi_stats_data(s, &n);
  p = casadi_stats_lookup(b, n, node, key);
  return p && casadi_stats_read_number(p, b+n, v) ? 0 : 1;
}

// C-REPLACE "static_cast<char>" "(char) "
// SYMBOL "stats_get_text"
// EXPORT
inline
int casadi_stats_get_text(const struct casadi_stats_sink* s, casadi_int node, const char* key,
    char* v, casadi_int cap, casadi_int* len) {
  // 0 if found as text; at most cap-1 chars in v, null-terminated if cap>0; *len: full
  // length without the null, truncated if >= cap (len may be null)
  const unsigned char *b, *p, *t;
  casadi_int n, nt, i;
  b = casadi_stats_data(s, &n);
  p = casadi_stats_lookup(b, n, node, key);
  if (!p || !casadi_cbor_read_text(p, b+n, &t, &nt)) return 1;
  if (casadi_stats_is(key, "name")) t = casadi_stats_name(t, &nt);
  if (cap>0) {
    for (i=0; i<nt && i<cap-1; ++i) v[i] = static_cast<char>(t[i]);
    v[i] = 0;
  }
  if (len) *len = nt;
  return 0;
}

// SYMBOL "stats_get_reals"
// EXPORT
inline
int casadi_stats_get_reals(const struct casadi_stats_sink* s, casadi_int node, const char* key,
    double* v, casadi_int cap, casadi_int* len) {
  // 0 if found as an array of reals, ints or bools; first cap in v, *len: full length
  // (len may be null)
  const unsigned char *b, *p;
  casadi_int n, m, i;
  double x;
  b = casadi_stats_data(s, &n);
  p = casadi_stats_lookup(b, n, node, key);
  if (!p || !(p = casadi_cbor_read_array(p, b+n, &m))) return 1;
  for (i=0; i<m; ++i) {
    p = casadi_stats_read_number(p, b+n, &x);
    if (!p) return 1;
    if (i<cap) v[i] = x;
  }
  if (len) *len = m;
  return 0;
}
