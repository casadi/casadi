import casadi

import pathlib
from collections import defaultdict

wiki = pathlib.Path("../../casadi.wiki/")

import re

entries = defaultdict(dict)

def parse_all(a,parent=[]):
    if hasattr(a,"__dict__"):
        for k,v in a.__dict__.items():
            if k.startswith("_") or k in ["os","np","numpy","re","inspect","types"]: continue
            if isinstance(v.__doc__,str):
                if "casadi/casadi/wiki/L_" in v.__doc__:
                    initial_indent = len(v.__doc__)-len(v.__doc__.lstrip())
                    target_indent = 6
                    assert initial_indent>=target_indent
                    for m in re.findall("casadi/casadi/wiki/L_(\w+)",v.__doc__):
                        drepr = "\n".join(["> " + e[initial_indent-target_indent:] for e in v.__doc__.split("\n")])
                        entries[m][v] = (v.__qualname__,drepr)

                                   
            parse_all(v,parent+[k])
            
        
parse_all(casadi)

for L,parts in entries.items():
  md = wiki / f"L_{L}.md"
  qualnames = [qualname for k,(qualname,descr) in parts.items()]
  qualnames = ",".join("`"+e+"`" for e in qualnames)
  descrs = [descr for k,(qualname,descr) in parts.items()]
  descrs = list(dict.fromkeys(descrs))
  drepr = "\n".join(descrs)
  with open(md,"w") as outp:
    outp.write(f"""
# Standard documentation for {qualnames}:
{drepr}

# Extra documentation

_To edit, see [writing tips](Extra-doc-tips)_.

""")

    if len(parts)>1:
        print(L,[qualname for k,(qualname,descr) in parts.items()])

             
