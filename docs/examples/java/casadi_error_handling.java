//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
import org.casadi.*;

/** \brief Demonstration of how to work with the error handling in CasADi
 *
 * In Java, CasADi errors (casadi_assert/casadi_error in C++) arrive as
 * org.casadi.CasadiException, an unchecked RuntimeException.
 */

public class casadi_error_handling {
  static boolean bad_test() {
    try {
      // This will fail
      if (!false) throw new CasadiException("Assertion \"false\" failed:\nFailing assert");

      // Returns true, but the code won't reach this place
      return true;
    } catch (RuntimeException e) {
      throw new CasadiException("bad_test3 failed:\n" + e.getMessage());
    }
  }

  public static void main(String[] args) {
    // Warning
    System.err.println("CasADi warning: This function will fail.");

    // Developer error (uncaught, it terminates the program like in C++; caught here to show the message)
    try {
      if (!bad_test()) throw new CasadiException("");
    } catch (CasadiException e) {
      System.out.println("terminate called after throwing a CasadiException");
      System.out.println("  what():  " + e.getMessage());
    }
  }
}
