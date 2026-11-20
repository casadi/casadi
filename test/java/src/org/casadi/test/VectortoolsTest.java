package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertThrows;

import org.casadi.DM;
import org.casadi.ca;
import org.junit.jupiter.api.Test;

/** Port of test/python/vectortools.py */
public class VectortoolsTest extends CasadiTestCase {

  @Test
  public void test_complement() {
    message("complement");
    long[] w = {2, 1, 4, 6};
    assertThrows(RuntimeException.class, () -> ca.complement(w, 3));
    assertThrows(RuntimeException.class, () -> ca.complement(w, 6));
    long[] wc = ca.complement(w, 8);
    double[] wcd = new double[wc.length];
    for (int i = 0; i < wc.length; i++) wcd[i] = wc[i];
    checkarray(new DM(wcd), dm(0, 3, 5, 7), "complement");
  }
}
