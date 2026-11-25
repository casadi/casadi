using System;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java VectortoolsTest (test/python/vectortools.py)</summary>
public class VectortoolsTest : CasadiTestCase {
  public void test_complement() {
    message("complement");
    long[] w = { 2, 1, 4, 6 };
    Check.Throws<Exception>(() => ca.complement(w, 3));
    Check.Throws<Exception>(() => ca.complement(w, 6));
    var wc = ca.complement(w, 8);
    checkarray(wc.Select(x => (double)x).ToArray(), ca.vertcat(0, 3, 5, 7), "complement");
  }
}

}
