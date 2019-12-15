// NOLINT(legal/copyright)
// SYMBOL "cache_check"
template<typename T1>
int cache_check(const T1* key, T1* cache, int* loc, casadi_int stride, casadi_int sz, casadi_int key_sz, T1** val) { // NOLINT(whitespace/line_length)
    char match;
    int i, c;
    casadi_int k;
    T1 *lkey;
    // Walk through cache locations
    for (i=0;i<sz;++i) {
      match = 1;
      if (loc[i]<0) { // Never filled
        loc[i] = i;
        *val = cache + i*stride;
        break;
      } else {
        *val = cache + loc[i]*stride;

        // Check for a key hit
        lkey = *val;
        for (k=0;k<key_sz;++k) {
          if (lkey[k]!=key[k]) {
            match = 0;
            break;
          }
        }
      }

      if (match) {
        // Move current location to front
        c = loc[i];
        for (k=i;k>0;--k) loc[k] = loc[k-1];
        loc[0] = c;
        // Indicate cache hit
        return 1;
      }
    }

    // Indicate cache miss
    return 0;
}
