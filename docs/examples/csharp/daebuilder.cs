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
using System;
using CasADi;

// Example on how to use the DaeBuilder class
// Joel Andersson, UW Madison 2017

public static class daebuilder {
  public static void Main() {

    // Start with an empty DaeBuilder instance
    var dae = ca.DaeBuilder("rocket");

    // Add input expressions
    MX a = dae.add("a", "parameter", "tunable");
    MX b = dae.add("b", "parameter", "tunable");
    MX u = dae.add("u", "input");
    MX h = dae.add("h");
    MX v = dae.add("v");
    MX m = dae.add("m");

    // Constants
    double g = 9.81; // gravity

    // Set ODE right-hand-side
    dae.eq(dae.der(h), v);
    dae.eq(dae.der(v), (u - a * ca.pow(v, 2)) / m - g);
    dae.eq(dae.der(m), -b * ca.pow(u, 2));

    // Specify initial conditions
    dae.set_start("h", 0);
    dae.set_start("v", 0);
    dae.set_start("m", 1);

    // Add meta information
    dae.set_unit("h", "m");
    dae.set_unit("v", "m/s");
    dae.set_unit("m", "kg");

    // Print DAE
    dae.disp(true);
  }

}
