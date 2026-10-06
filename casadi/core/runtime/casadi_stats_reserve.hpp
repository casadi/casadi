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

// SYMBOL "stats_window"
inline
unsigned char* casadi_stats_window(struct casadi_stats_sink* s, casadi_int len, casadi_int pos) {
  // Null if it does not fit; the first misfit leaves a 0xff break to end the stream
  if (len <= s->cap - pos) return s->p + pos;
  if (pos < s->cap) s->p[pos] = 0xff;
  return 0;
}

// SYMBOL "stats_bump"
inline
unsigned char* casadi_stats_bump(struct casadi_stats_sink* s, casadi_int len, casadi_int* pos) {
  // Default reserve: bump needed (atomic or under mutex when thread-safe)
  // C-THREAD-SAFE
#if defined(CASADI_ATOMIC_FETCH_ADD)
  // C-THREAD-SAFE
  *pos = CASADI_ATOMIC_FETCH_ADD(&s->needed, len);
  // C-THREAD-SAFE
#elif defined(CASADI_THREAD_TYPE)
  // C-THREAD-SAFE
  CASADI_MUTEX_LOCK(&casadi_stats_mutex);
  // C-THREAD-SAFE
  *pos = s->needed; s->needed += len;
  // C-THREAD-SAFE
  CASADI_MUTEX_UNLOCK(&casadi_stats_mutex);
  // C-THREAD-SAFE
#else
  *pos = s->needed; s->needed += len;
  // C-THREAD-SAFE
#endif
  return casadi_stats_window(s, len, *pos);
}

// SYMBOL "stats_reserve"
// EXPORT
inline
unsigned char* casadi_stats_reserve(struct casadi_stats_sink* s, casadi_int len, casadi_int* pos) {
  // Where to write len bytes, at stream offset *pos; null if they do not fit
  // C-FUNCTION-POINTERS
  if (s->reserve) return s->reserve(s, len, pos);
  // Cannot call the owner's callbacks: drop the record
  // C-NO-FUNCTION-POINTERS
  if (s->reserve || s->commit) *pos = -1;
  // C-NO-FUNCTION-POINTERS
  if (s->reserve || s->commit) return 0;
  return casadi_stats_bump(s, len, pos);
}

// SYMBOL "stats_commit"
inline
void casadi_stats_commit(struct casadi_stats_sink* s, casadi_int pos) {
  // C-FUNCTION-POINTERS
  if (s->commit) s->commit(s, pos);
  // C-NO-FUNCTION-POINTERS
  (void) s; (void) pos;
}

// SYMBOL "stats_init_sink"
// EXPORT
inline
void casadi_stats_init_sink(struct casadi_stats_sink* s, unsigned char* p, casadi_int cap) {
  // Libraries called into reserve through this file's bump
  s->p = p;
  s->cap = cap;
  s->needed = 0;
  s->reserve = 0;
  s->commit = 0;
  s->data = 0;
  // C-FUNCTION-POINTERS
  s->reserve = casadi_stats_bump;
}

// SYMBOL "stats_clear"
// EXPORT
inline
void casadi_stats_clear(struct casadi_stats_sink* s) {
  s->needed = 0;
}

// SYMBOL "stats_nbytes"
// EXPORT
inline
casadi_int casadi_stats_nbytes(const struct casadi_stats_sink* s) {
  // Including dropped records
  return s->needed;
}

// SYMBOL "stats_truncated"
// EXPORT
inline
int casadi_stats_truncated(const struct casadi_stats_sink* s) {
  return s->needed > s->cap;
}

// SYMBOL "stats_data"
// EXPORT
inline
const unsigned char* casadi_stats_data(const struct casadi_stats_sink* s, casadi_int* n) {
  // Recorded bytes held in memory
  *n = s->needed < s->cap ? s->needed : s->cap;
  return s->p;
}
