from os.path import join, basename, splitext, exists
from glob import glob
import yaml
import pandas as pd

def load_raw_data(folder):
    if not exists(folder):
        raise FileNotFoundError(folder)
    output_files = glob(join(folder, '*.yaml'))
    raw_data = {}
    for filename in output_files:
        if basename(filename) == 'parameters.yaml':
            continue
        with open(filename, 'r') as f:
            all_content = yaml.safe_load_all(f)
            content = next(all_content)
            name = splitext(basename(filename))[0]
            raw_data[name] = content
    return raw_data

def convert_data(raw_data):
    data = []
    for name, content in raw_data.items():
        element = {
            'name': content.get('name', name),
            'status': content['status'],
            # 'time': timedelta(seconds=content['elapsed time']),
            'time': float(content['elapsed time']),
            'inner iterations': content['inner iterations'],
            'outer iterations': content['outer iterations'],
            'inner convergence failures': content['inner convergence failures'],
            'initial penalty reduced': content.get('initial penalty reduced', -1),
            'penalty reduced': content.get('penalty reduced', -1),
            'f': float(content['f']),
            'ε': float(content['ε']),
            'δ': float(content['δ']),
            'f evaluations': content['counters']['f'],
            'grad_f evaluations': content['counters']['grad_f'],
            'g evaluations': content['counters']['g'],
            'grad_g evaluations': content['counters']['grad_g'],
            'linesearch failures': content['linesearch failures'],
            'L-BFGS failures': content['L-BFGS failures'],
            'L-BFGS rejected': content['L-BFGS rejected'],
            '‖Σ‖': float(content['‖Σ‖']),
            '‖x‖': float(content['‖x‖']),
            '‖y‖': float(content['‖y‖']),
            'n': int(content.get('n', -1)),
            'm': int(content.get('m', -1)),
            'box constr x': int(content.get('box constraints x', -1)),
        }
        data.append(element)
    df = pd.DataFrame(data)
    df.set_index('name', inplace=True)
    df.sort_index(inplace=True)
    return df
