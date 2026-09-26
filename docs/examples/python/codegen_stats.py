#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# Solver statistics from generated code, written to a caller-owned buffer (no I/O, no malloc)
import casadi as ca
import os
import subprocess

x = ca.MX.sym("x", 2)
p = ca.MX.sym("p")
nlp = {"x": x, "p": p, "f": (1-x[0])**2 + p*(x[1]-x[0]**2)**2, "g": x[0]+x[1]}
solver = ca.nlpsol("solver", "ipopt", nlp, {"ipopt.print_level": 0, "print_time": False})

# Option "stats" adds solver_with_stats(..., struct casadi_stats_sink* sink, casadi_int parent)
solver.generate("solver.c", {"stats": True, "with_header": True})

main = r"""
#include <stdio.h>
#include "solver.h"
static unsigned char buffer[1 << 16];
int main(void) {
  casadi_real x0[2] = {-1, 1}, p = 100, ubg = 1, x[2];
  const casadi_real* arg[64] = {x0, &p, 0, 0, 0, &ubg, 0, 0};
  casadi_real* res[64] = {x};
  casadi_int iw[64];
  casadi_real w[256];
  struct casadi_stats_sink stats;
  casadi_int n, iter_count;
  const unsigned char* data;
  FILE* f;
  int mem, flag;
  stats = solver_stats_make_sink(buffer, sizeof(buffer));
  mem = solver_checkout();
  flag = solver_with_stats(arg, res, iw, w, mem, &stats, -1);   /* -1: root call */
  printf("flag %d, x = [%g, %g], %lld stats bytes%s\n", flag, x[0], x[1],
         (long long) solver_stats_nbytes(&stats),
         solver_stats_truncated(&stats) ? " (truncated)" : "");
  data = solver_stats_data(&stats, &n);
  f = fopen("stats.cbor", "wb");
  fwrite(data, 1, n, f);
  fclose(f);
  if (!solver_get_stat_int(&stats, "solver", "iter_count", &iter_count)) {
    printf("%lld iterations\n", (long long) iter_count);
  }
  solver_release(mem);
  return flag;
}
"""
with open("main.c", "w") as f:
  f.write(main)

inc = ca.GlobalOptions.getCasadiIncludePath()
lib = ca.GlobalOptions.getCasadiPath()
subprocess.run(["gcc", "main.c", "solver.c", "-I" + inc, "-I" + inc + "/coin-or", "-L" + lib,
                "-Wl,-rpath," + lib, "-lipopt", "-lm", "-o", "main"], check=True)
# Windows: no rpath
env = dict(os.environ, PATH=lib + os.pathsep + os.environ.get("PATH", ""))
subprocess.run([os.path.join(".", "main")], check=True, env=env)

stats = ca.StatsRecorder.load("stats.cbor")
print(stats.to_json()[:300], "...")
print("iterations:", stats.get_stat("solver", "iter_count"))
print("objective per iteration:", stats.get_stat("solver", "iterations")["obj"])

# Same stream from the virtual machine
vm = ca.StatsRecorder()
solver(vm, x0=[-1, 1], p=100, ubg=1)
print("VM:", vm.get_stat("solver", "return_status"))
