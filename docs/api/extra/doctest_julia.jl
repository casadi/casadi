# Run Julia doctest blocks (Documenter style: `julia> expr` + expected output) against CasADiNative.
#
# Usage: julia docs/api/extra/doctest_julia.jl <blocks.txt>
# The file holds blocks separated by lines `#### <name>`; each block is the text between
# \doctest{julia} and \enddoctest. Needs CASADI_JL pointing at a dir holding CasADiNative.jl
# (default: the worktree build dir swig/julia/CasADiNative). Exit code = number of failing blocks.

const _jl_dir = get(ENV, "CASADI_JL", joinpath(@__DIR__, "..", "..", "..", "build", "swig", "julia", "CasADiNative"))
include(joinpath(_jl_dir, "CasADiNative.jl"))

function parse_blocks(path)
    blocks = Tuple{String,String}[]
    name, buf = "", String[]
    for line in eachline(path)
        if startswith(line, "#### ")
            isempty(name) || push!(blocks, (name, join(buf, "\n")))
            name, buf = strip(line[6:end]), String[]
        else
            push!(buf, line)
        end
    end
    isempty(name) || push!(blocks, (name, join(buf, "\n")))
    blocks
end

# Split a block into (source, expected) pairs. Continuation lines are indented by 7 spaces.
function parse_examples(text)
    out = Tuple{String,String}[]
    lines = split(text, "\n")
    i = 1
    while i <= length(lines)
        l = lines[i]
        if startswith(l, "julia> ")
            src = [l[8:end]]
            i += 1
            while i <= length(lines) && startswith(lines[i], "       ")
                push!(src, lines[i][8:end]); i += 1
            end
            exp = String[]
            while i <= length(lines) && !startswith(lines[i], "julia> ") && !isempty(strip(lines[i]))
                push!(exp, lines[i]); i += 1
            end
            push!(out, (join(src, "\n"), join(exp, "\n")))
        else
            i += 1
        end
    end
    out
end

norm(s) = join(filter(!isempty, [strip(l) for l in split(strip(s), "\n")]), "\n")

function matches(expected, got)
    e, g = norm(expected), norm(got)
    e == g && return true
    occursin("...", e) || return false
    pat = join(map(x -> replace(x, r"[.*+?^\${}()|\[\]\\]" => s -> "\\" * s), split(e, "...")), ".*")
    occursin(Regex("^" * replace(pat, "\n" => "\\n") * "\$", "s"), g)
end

function run_block(name, text)
    mod = Module(Symbol("Doctest_", replace(name, r"\W" => "_")))
    Core.eval(mod, :(import $(nameof(@__MODULE__)).CasADiNative as ca))
    Core.eval(mod, :(const CasADiNative = $(CasADiNative)))
    ok = true
    for (src, expected) in parse_examples(text)
        path, io = mktemp()
        got = ""
        try
            result = redirect_stdout(io) do
                Core.eval(mod, Meta.parseall(src))
            end
            close(io)
            got = read(path, String)
            if !endswith(rstrip(src), ";") && result !== nothing
                rio = IOBuffer(); show(rio, MIME("text/plain"), result)
                got *= String(take!(rio))
            end
        catch err
            isopen(io) && close(io)
            got = read(path, String) * "ERROR: " * sprint(showerror, err)
        end
        rm(path, force=true)
        if !matches(expected, got)
            ok = false
            println("**********************************************************************")
            println("$name\nFailed example:\n    $(replace(src, "\n" => "\n    "))\nExpected:\n    $(replace(expected, "\n" => "\n    "))\nGot:\n    $(replace(rstrip(got), "\n" => "\n    "))")
        end
    end
    ok
end

blocks = parse_blocks(ARGS[1])
nfail = count(!, [run_block(n, t) for (n, t) in blocks])
println("doctest_julia: $(length(blocks)) blocks, $nfail failures")
exit(nfail)
