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
