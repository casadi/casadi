# src/CasADiNative.jl of the registered CasADiNative package (github.com/casadi/CasADiNative.jl).
# bindings.jl locates the binaries; the ergonomics layer is spliced in as this module's body.
module CasADiNative

include("bindings.jl")
import LinearAlgebra
import OpenBLAS32_jll

function __init__()
    # Ipopt_jll needs LP64 BLAS behind libblastrampoline (same as Ipopt.jl).
    if !Base.any(lib -> lib.interface == :lp64, LinearAlgebra.BLAS.lbt_get_config().loaded_libs)
        LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
    end
    # casadi finds its plugins via CASADI_PLUGIN_SEARCH_PATH; prepend the plugin dir.
    libdir = _plugin_dir()
    sep = Sys.iswindows() ? ';' : ':'
    ENV["CASADI_PLUGIN_SEARCH_PATH"] = haskey(ENV, "CASADI_PLUGIN_SEARCH_PATH") ?
        string(libdir, sep, ENV["CASADI_PLUGIN_SEARCH_PATH"]) : libdir
    _layer_init()
end

include(_BODY)

end # module CasADiNative
