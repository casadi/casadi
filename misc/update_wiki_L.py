import casadi

import pathlib

wiki = pathlib.Path("../../casadi.wiki/")

import re

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
                        md = wiki / f"L_{m}.md"
                        
                        with open(md,"w") as outp:
                            drepr = "\n".join(["> " + e[initial_indent-target_indent:] for e in v.__doc__.split("\n")])
                            
                            outp.write(f"""
# Standard documentation for `{v.__qualname__}`:
{drepr}

# Extra documentation

_To edit, see [writing tips](Extra-doc-tips)_.

""")
                        
            parse_all(v,parent+[k])
        
        
parse_all(casadi)
