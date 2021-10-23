import functools
from os.path import join, dirname, abspath
from typing import List, cast
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from difflib import unified_diff
# import tikzplotlib
from util.loader import load_raw_data, convert_data

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 15,
    "lines.linewidth": 1,
})

machine = "XPS-15-9560"
this_folder = dirname(__file__) if "__file__" in globals() else abspath("")
root_folder = dirname(dirname(this_folder))
out_folder = join(root_folder, "test", "testresults", machine)
get_test_result_folder = lambda testname: join(out_folder, testname, "CUTEst")
fig_out_folder = "/tmp/fig"

names = {
    # "gaapga-29-5": "Guarded Anderson accelerated PGA", 
    # "panoc-2nd-29-5-upd-ls-crit": "Structured PANOC (FD, L-BFGS, impr.)",
    # "lbfgsbpp-29-5": "L-BFGS-B",
    # "panoc-2nd-29-5-ad": "Structured PANOC (AD, L-BFGS, impr.)",
    # "panoc-2nd-29-5-baseline": "Structured PANOC (FD, L-BFGS, orig.)",
    # "panoc-29-5-baseline": "PANOC (orig.)",
    # "panoc-29-5-cbfgs": "PANOC (CBFGS, orig.)",
    # "panoc-29-5-upd-ls-crit": "PANOC (impr.)",
    # "pga-29-5": "Projected gradient",
    # "panoc-2nd-29-5-single-penalty": "Structured PANOC (single penalty factor)",
    # "panoc-2nd-4-6-no-backtrack": "Structured PANOC (no ALM backtracking)",
    "strucpanoc-20-10-baseline": "Structured PANOC (baseline)",
    "strucpanoc-21-10-baseline": "Structured PANOC (baseline new)",
    "strucpanoc-20-10-hessheur": "Structured PANOC (Hessian heuristic /1)",
    "strucpanoc-20-10-hessheur10": "Structured PANOC (Hessian heuristic /10)",
    "strucpanoc-20-10-hessheur50": "Structured PANOC (Hessian heuristic /50)",
    "strucpanoc-21-10-hessheur10": "Structured PANOC (Hessian heuristic new /10)",
    "strucpanoc-21-10-hessheur10-n":"Structured PANOC (Hessian heuristic new2 /10)",
    "strucpanoc-21-10-hessheur10-n2":"Structured PANOC (Hessian heuristic new2-2 /10)",
    "strucpanoc-21-10-hessheur10-n2-5":"Structured PANOC (Hessian heuristic new2-2 /10x5)",
    "strucpanoc-21-10-hessheur2-n2-5":"Structured PANOC (Hessian heuristic new2-2 /2x5)",
    "strucpanoc-21-10-hessheur20-n2-5":"Structured PANOC (Hessian heuristic new2-2 /00x5)",
}
labelnames = {
    "time": "Time",
    "inner iterations": "Inner iterations",
    "inner convergence failures": "Inner convergence failures",
    "f evaluations": "Objective evaluations",
    "grad_f evaluations": "Objective gradient evaluations",
    "outer iterations": "Outer interations",
    "L-BFGS rejected": "L-BFGS update rejections",
    "average τ": "Average line search parameter $\\tau$",
    "avg time per it": "Average time per iteration",
    "‖x‖": "Norm of decision variables $\\|x\\|$",
    "‖y‖": "Norm of Lagrange multipliers $\\|y\\|$",
    "‖Σ‖": "Norm of penalty factors $\\|\\Sigma\\|$",
}
labelnoglog = {
    "inner convergence failures",
    "average τ",
}
labelhigherbetter = {
    "average τ"
}

sel = [
    # ["panoc-29-5-baseline", "panoc-29-5-upd-ls-crit",         ["time", "inner iterations", "avg time per it", "average τ"]],
    # ["panoc-2nd-29-5-baseline", "panoc-2nd-29-5-upd-ls-crit", ["time", "inner iterations", "avg time per it", "average τ"]],
    # ["panoc-2nd-29-5-ad", "panoc-2nd-29-5-upd-ls-crit",       ["time", "inner iterations", "avg time per it", "grad_f evaluations"]],
    # ["panoc-29-5-upd-ls-crit", "panoc-2nd-29-5-upd-ls-crit",  ["time", "inner iterations", "avg time per it", "average τ"]],
    # ["panoc-29-5-baseline", "panoc-29-5-cbfgs",               ["time", "inner iterations", "L-BFGS rejected"]],
    # ["pga-29-5", "panoc-2nd-29-5-upd-ls-crit",                ["time", "inner iterations", "avg time per it", "outer iterations"]],
    # ["gaapga-29-5", "panoc-2nd-29-5-upd-ls-crit",             ["time", "inner iterations", "avg time per it", "outer iterations"]],
    # ["lbfgsbpp-29-5", "panoc-2nd-29-5-upd-ls-crit",           ["time", "inner iterations", "avg time per it", "f evaluations", "‖x‖", "‖y‖", "‖Σ‖"]],
    # ["panoc-2nd-29-5-single-penalty", "panoc-2nd-29-5-upd-ls-crit",  ["time", "inner iterations", "‖Σ‖", "outer iterations"]],
    # ["panoc-2nd-4-6-no-backtrack", "panoc-2nd-29-5-upd-ls-crit",  ["time", "inner iterations"]],
    ["strucpanoc-21-10-hessheur10-n2-5", "strucpanoc-21-10-hessheur20-n2-5", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    ["strucpanoc-21-10-baseline", "strucpanoc-21-10-hessheur2-n2-5", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    ["strucpanoc-21-10-baseline", "strucpanoc-21-10-hessheur10-n2-5", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    ["strucpanoc-21-10-hessheur10-n2", "strucpanoc-21-10-hessheur10-n2-5", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-baseline", "strucpanoc-21-10-hessheur10-n2", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-baseline", "strucpanoc-21-10-hessheur10-n", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-baseline", "strucpanoc-21-10-hessheur10", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-baseline", "strucpanoc-20-10-hessheur", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-baseline", "strucpanoc-20-10-hessheur10", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-hessheur", "strucpanoc-20-10-hessheur10", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-hessheur10", "strucpanoc-20-10-hessheur50", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
    # ["strucpanoc-20-10-hessheur10", "strucpanoc-21-10-hessheur10", ["time", "inner iterations", "avg time per it", "average τ", "f evaluations", "grad_f evaluations"]],
]

outf = open('/tmp/input.tex', 'w')

for s in sel:
    testnames = cast(List[str], s[0:2])
    dispnames = list(map(lambda n: names[n], testnames))
    folders = list(map(get_test_result_folder, testnames))
    raw_data = list(map(load_raw_data, folders))
    dfs = list(map(convert_data, raw_data))
    paramfiles = list(map(lambda f: join(f, "parameters.yaml"), folders))
    params = []
    for pf in paramfiles:
        with open(pf) as f:
            params.append(f.readlines())
    sys.stdout.writelines(
        unified_diff(params[0], params[1], fromfile=testnames[0], tofile=testnames[1])
    )
    print("\n")

    both_converged = (dfs[0]["status"] == "Converged") & (dfs[1]["status"] == "Converged")


    def df_stats(df):
        conv = df["status"].value_counts()["Converged"]
        tot = df["status"].count()
        tot_time = df["time"].sum()
        conv_time = df.where(df["status"] == "Converged")["time"].sum()
        print(f"Converged:      {conv}/{tot} = {100*conv/tot:.02f}%")
        print(f"Total time:     {tot_time:.03f}s")
        print(f"Converged time: {conv_time:.03f}s")


    print(f"{dispnames[0]}\n---\n")
    df_stats(dfs[0])
    print("\n")
    print(f"{dispnames[1]}\n---\n")
    df_stats(dfs[1])
    print("\n")

    statusses = dfs[0][["status"]].join(dfs[1][["status"]], lsuffix=" 0", rsuffix=" 1")
    not_conv_to_conv = (statusses["status 0"] != "Converged") & (
        statusses["status 1"] == "Converged"
    )
    conv_to_not_conv = (statusses["status 0"] == "Converged") & (
        statusses["status 1"] != "Converged"
    )

    print(
        f"{len(dfs[0][not_conv_to_conv])} tests that didn't converge before do converge after the change"
    )
    print(
        f"{len(dfs[0][conv_to_not_conv])} tests that converged before no longer converge after the change"
    )

    print()

    print("The following tests went from not converging to converging")
    print(dfs[0][not_conv_to_conv])
    print(dfs[1][not_conv_to_conv])

    print()

    print("The following tests went from converging to no longer converging")
    print(dfs[0][conv_to_not_conv])
    print(dfs[1][conv_to_not_conv])

    print()

    print(f'[0] #L-BFGS rej == 0: {np.sum(dfs[0][both_converged]["L-BFGS rejected"] == 0)}')
    print(f'[1] #L-BFGS rej == 0: {np.sum(dfs[1][both_converged]["L-BFGS rejected"] == 0)}')

    out_name = "-vs-".join(testnames)

    def compare_plot(ax, t0, t1, log=True, regression=False, lowerbetter=True):
        plotfun = ax.loglog if log else ax.plot
        lower = t1 < t0
        higher = t1 > t0
        nochange = np.logical_not(np.logical_or(lower, higher))
        neutralcolor = "b"
        lowercolor = "g" if lowerbetter else "r"
        highercolor = "r" if lowerbetter else "g"
        drawLines = []
        for tt0, tt1 in zip(t0, t1):
            if tt1 != tt0:
                drawLines.append(np.array([tt0, tt0]))
                drawLines.append(np.array([tt0, tt1]))
        plotfun(
            *drawLines,
            color='lightgray',
            linewidth=0.5,
        )
        # drawLines = []
        # for tt0, tt1 in zip(t0, t1):
        #     if tt1 > tt0:   
        #         drawLines.append(np.array([tt0, tt0]))
        #         drawLines.append(np.array([tt0, tt1]))
        # plotfun(
        #     *drawLines,
        #     color='lightgray',
        #     linewidth=0.5,
        # )
        plotfun(
            t0[nochange],
            t1[nochange],
            ".",
            color=neutralcolor,
            markersize="4",
            label=f"Identical ({np.sum(nochange)})",
        )
        plotfun(
            t0[lower],
            t1[lower],
            ".",
            color=lowercolor,
            markersize="4",
            label=f"Lower ({np.sum(lower)})",
        )
        plotfun(
            t0[higher],
            t1[higher],
            ".",
            color=highercolor,
            markersize="4",
            label=f"Higher ({np.sum(higher)})",
        )
        m, M = min(np.min(t0[np.isfinite(t0)]), np.min(t1[np.isfinite(t1)])), max(np.max(t0[np.isfinite(t0)]), np.max(t1[np.isfinite(t1)]))
        plotfun([m, M], [m, M], "k:", linewidth=1)
        if regression:
            if log:
                a, b = np.polyfit(np.log10(t0), np.log10(t1), deg=1)
                reg_eval = lambda x: 10 ** (np.log10(x) * a + b)
            else:
                a, b = np.polyfit(t0, t1, deg=1)
                reg_eval = lambda x: x * a + b
            print(f"{a=}, {b=}")
            plotfun([m, M], [reg_eval(m), reg_eval(M)], "r:", linewidth=1)

    for prop in s[2]:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        compare_plot(
            ax,
            dfs[0][both_converged][prop],
            dfs[1][both_converged][prop],
            log=prop not in labelnoglog,
            lowerbetter=prop not in labelhigherbetter,
        )
        ax.set_xlabel(f'{dispnames[0]}', fontsize=12) # ({dfs[0]["solver"][0]})')
        ax.set_ylabel(f'{dispnames[1]}', fontsize=12) # ({dfs[1]["solver"][0]})')
        ax.set_title(labelnames[prop])
        ax.axis("equal")
        ax.legend(loc='upper left', fontsize=10)
        plt.tight_layout()
        plt.savefig(join(fig_out_folder, out_name + '-' + prop + ".pdf"))

        # tikzplotlib.clean_figure()
        fname = out_name + '-' + prop + ".tex"
        outf.write(f'\\input{{fig/{fname}}}\n')
        # tikzplotlib.save(join(fig_out_folder, fname), 
        #                 axis_width='0.5\linewidth')
outf.close()
plt.show()
