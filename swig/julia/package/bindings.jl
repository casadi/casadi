# Binary locations of the registered package: the JLLs (casadi CI substitutes its build tree).
using CasADi_jll
using CasADiSWIG_jll
_plugin_dir() = dirname(CasADi_jll.libcasadi_path)
const _BODY = CasADiSWIG_jll.casadi_native_body
