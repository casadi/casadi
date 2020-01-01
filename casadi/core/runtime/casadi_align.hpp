// NOLINT(legal/copyright)

// SYMBOL "align"
template<typename T1>
T1* casadi_align(T1* a, int p) {
    // e.g p = 64= 1000000
    p--; // 0111111
    // C-REPLACE "reinterpret_cast<uintptr_t>" "(uintptr_t)"
    uintptr_t r = reinterpret_cast<uintptr_t>(a);
    if (r & p) { // any trailing 6 bits set
    // C-REPLACE "reinterpret_cast<T1*>" "(T1*)"
      return reinterpret_cast<T1*>((r | p) + 1);
    } else {     // none set
      return a;
    }
}
