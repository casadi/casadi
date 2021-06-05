from os.path import join, basename, splitext, exists
from glob import glob
import yaml
import pandas as pd
import numpy as np
from datetime import timedelta


def load_raw_data(folder):
    if not exists(folder):
        raise FileNotFoundError(folder)
    output_files = glob(join(folder, "*.yaml"))
    raw_data = {}
    for filename in output_files:
        if basename(filename) == "parameters.yaml":
            continue
        with open(filename, "r") as f:
            all_content = yaml.safe_load_all(f)
            content = next(all_content)
            name = splitext(basename(filename))[0]
            raw_data[name] = content
    return raw_data


def convert_data(raw_data):
    data = []
    for name, content in raw_data.items():
        element = {
            "name": content.get("name", name),
            "solver": content.get("solver", "<unknown>"),
            "status": content["status"],
            # 'time': timedelta(seconds=content['elapsed time']),
            "time": float(content["elapsed time"]),
            "inner iterations": content["inner"]["iterations"],
            "outer iterations": content["outer iterations"],
            "inner convergence failures": content["inner convergence failures"],
            "initial penalty reduced": content.get("initial penalty reduced", -1),
            "penalty reduced": content.get("penalty reduced", -1),
            "f": float(content["f"]),
            "ε": float(content["ε"]),
            "δ": float(content["δ"]),
            "f evaluations": int(content["counters"]["f"]),
            "grad_f evaluations": int(content["counters"]["grad_f"]),
            "g evaluations": int(content["counters"]["g"]),
            "grad_g evaluations": int(content["counters"]["grad_g_prod"]),
            "linesearch failures": int(content["inner"].get("linesearch_failures", 0)),
            "newton failures": int(content["inner"].get("newton_failures", 0)),
            "L-BFGS failures": int(content["inner"].get("lbfgs_failures", 0)),
            "L-BFGS rejected": int(content["inner"].get("lbfgs_rejected", 0)),
            "τ=1 accepted": int(content["inner"].get("τ_1_accepted", 0)),
            "accelerated_steps_accepted": int(
                content["inner"].get("accelerated_steps_accepted", 0)
            ),
            "sum τ": float(content["inner"].get("sum_τ", 0)),
            "count τ": int(content["inner"].get("count_τ", 0)),
            "‖Σ‖": float(content["‖Σ‖"]),
            "‖x‖": float(content["‖x‖"]),
            "‖y‖": float(content["‖y‖"]),
            "n": int(content.get("n", -1)),
            "m": int(content.get("m", -1)),
            "box constr x": int(content.get("box constraints x", -1)),
        }
        with np.errstate(divide="ignore", invalid="ignore"):
            element["average τ"] = np.divide(
                element["sum τ"],
                element["count τ"],
            )
            element["fraction τ=1 accepted"] = np.divide(
                element["τ=1 accepted"],
                element["count τ"],
            )
            element["avg time per it"] = np.divide(
                element["time"],
                element["inner iterations"],
            )
        data.append(element)
    df = pd.DataFrame(data)
    df.set_index("name", inplace=True)
    df.sort_index(inplace=True)
    return df
