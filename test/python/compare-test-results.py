import functools
from os.path import join, dirname, abspath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
from difflib import unified_diff
from util.loader import load_raw_data, convert_data

machine = "XPS-15-9500"
this_folder = dirname(__file__) if "__file__" in globals() else abspath("")
root_folder = dirname(dirname(this_folder))
out_folder = join(root_folder, "test", "testresults", machine)
get_test_result_folder = lambda testname: join(out_folder, testname, "CUTEst")

testnames = ["panoc-29-5-baseline", "panoc-29-5-cbfgs"]
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


print(f"{testnames[0]}\n---\n")
df_stats(dfs[0])
print("\n")
print(f"{testnames[1]}\n---\n")
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
with pd.ExcelWriter(join(out_folder, out_name + ".xlsx")) as writer:
    dfs[0].to_excel(writer, sheet_name=testnames[0])
    dfs[1].to_excel(writer, sheet_name=testnames[1])


def compare_plot(ax, t0, t1, log=True, regression=False, lowerbetter=True):
    plotfun = ax.loglog if log else ax.plot
    lower = t1 < t0
    higher = t1 > t0
    nochange = np.logical_not(np.logical_or(lower, higher))
    neutralcolor = "b"
    lowercolor = "g" if lowerbetter else "r"
    highercolor = "r" if lowerbetter else "g"
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
    m, M = min(np.min(t0), np.min(t1)), max(np.max(t0), np.max(t1))
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


cmp_prop = [
    [
        "inner iterations",
        "outer iterations",
        ("inner convergence failures", {"log": False}),
        # "f",
    ],
    ["time", "f evaluations", "grad_f evaluations"],
    [
        ("‖Σ‖"),
        # ("fraction τ=1 accepted", {"log": False, "lowerbetter": False}),
        ("L-BFGS rejected"),
        ("average τ", {"log": False, "lowerbetter": False}),
    ],
]

if functools.reduce(lambda a, b: a * b, map(lambda df: df["count τ"].sum(), dfs)) == 0:
    cmp_prop.pop()

nr = len(cmp_prop)
nc = max(map(len, cmp_prop))

fig, axs = plt.subplots(nr, nc, squeeze=False, figsize=(20 * nc / 3, 20 * nr / 3))
for r, rprop in enumerate(cmp_prop):
    for c, prop in enumerate(rprop):
        prop, opt = prop if type(prop) is tuple else (prop, {})
        ax = axs[r, c]
        compare_plot(
            ax,
            dfs[0][both_converged][prop],
            dfs[1][both_converged][prop],
            **opt,
        )
        if r == nr - 1:
            ax.set_xlabel(f'{testnames[0]} ({dfs[0]["solver"][0]})')
        if c == 0:
            ax.set_ylabel(f'{testnames[1]} ({dfs[1]["solver"][0]})')
        ax.set_title(prop)
        ax.axis("equal")
        ax.legend()
plt.tight_layout()
plt.savefig(join(out_folder, out_name + ".pdf"))


fig, axs = plt.subplots(2, 1, squeeze=False, figsize=(20 * 1 / 3, 20 * 2 / 3), sharex='all', sharey='all')
axs[0][0].loglog(
    dfs[0][both_converged]["n"],
    dfs[0][both_converged]["time"],
    ".",
)
axs[0][0].set_title(f'{testnames[0]}\n({dfs[0]["solver"][0]})')
axs[0][0].set_ylabel('Runtime')
axs[1][0].loglog(
    dfs[1][both_converged]["n"],
    dfs[1][both_converged]["time"],
    ".",
)
axs[1][0].set_title(f'{testnames[1]}\n({dfs[1]["solver"][0]})')
axs[1][0].set_ylabel('Runtime')
axs[1][0].set_xlabel('Problem size $n$')

plt.show()
