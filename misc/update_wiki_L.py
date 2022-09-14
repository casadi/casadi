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
                    for m in re.findall("casadi/casadi/wiki/L_(\w+)",v.__doc__):
                        md = wiki / f"L_{m}.md"
                        with open(md,"w") as outp:
                            outp.write(f"""
# Standard documentation for `{v.__qualname__}`:
```
{v.__doc__}
```
# Extra documentation

""")
                        
            parse_all(v,parent+[k])
        
        
parse_all(casadi)
