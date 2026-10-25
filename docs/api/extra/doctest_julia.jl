# Run Julia doctest blocks (Documenter style: `julia> expr` + expected output) against CasADiNative.
#
# Usage: julia docs/api/extra/doctest_julia.jl <blocks.txt>
# The file holds blocks separated by lines `#### <name>`; each block is the text between
# \doctest and \enddoctest of a \lang{julia} group. Needs CASADI_JL pointing at a dir holding CasADiNative.jl
# (default: the worktree build dir swig/julia/CasADiNative). Exit code = number of failing blocks.
# Each block runs in a fresh module where `ca` is CasADiNative (`import CasADiNative as ca`).

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

# Dict iteration order differs between Julia versions: compare a Dict's entry lines sorted
function sortdicts(s)
    lines = split(s, "\n")
    out, i = String[], 1
    while i <= length(lines)
        push!(out, lines[i])
        if occursin(r"Dict\{.*\} with \d+ entr(y|ies):$", lines[i])
            j = i + 1
            while j <= length(lines) && startswith(lines[j], "  ")
                j += 1
            end
            append!(out, sort(lines[i+1:j-1]))
            i = j
        else
            i += 1
        end
    end
    join(out, "\n")
end

norm(s) = join(filter(!isempty, [strip(l) for l in split(sortdicts(strip(s)), "\n")]), "\n")

# Numbers as in doctest_headers.py: half a unit in the last written digit, at least 1e-9 relative
const NUMBER = r"(?<![\w.@])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?(?![\w.])"
numbers(s) = [m for m in eachmatch(NUMBER, s) if !occursin(r"^[-+]?0\d", m.match)]
function tolerance(lit)
    m = match(r"^[-+]?\d*(?:\.(\d*))?(?:[eE]([-+]?\d+))?$", lit)
    tol = m[1] === nothing && m[2] === nothing ? 0.0 :
        0.5 * 10.0^(parse(Int, something(m[2], "0")) - length(something(m[1], "")))
    x = abs(parse(Float64, lit))
    max(tol, 1e-9 * (x == 0 ? 1.0 : x))
end
function as_expected(expected, got)
    w, g = numbers(expected), numbers(got)
    length(w) == length(g) || return got
    out, pos = IOBuffer(), 1
    for (a, b) in zip(w, g)
        if abs(parse(Float64, a.match) - parse(Float64, b.match)) <= tolerance(a.match) * (1 + 1e-9)
            print(out, got[pos:prevind(got, b.offset)], a.match)
            pos = b.offset + ncodeunits(b.match)
        end
    end
    print(out, got[pos:end])
    String(take!(out))
end

function matches(expected, got)
    e = norm(expected)
    g = as_expected(e, norm(got))
    e == g && return true
    occursin("...", e) || return false
    pat = join(map(x -> replace(x, r"[.*+?^\${}()|\[\]\\]" => s -> "\\" * s), split(e, "...")), ".*")
    occursin(Regex("^" * replace(pat, "\n" => "\\n") * "\$", "s"), g)
end

# A fresh module per block, with CasADiNative imported as `ca`
function block_module(name)
    mod = Module(Symbol("Doctest_", replace(name, r"\W" => "_")))
    Core.eval(mod, :(import $(nameof(@__MODULE__)).CasADiNative as ca))
    Core.eval(mod, :(const CasADiNative = $(CasADiNative)))
    mod
end

# Printed output plus, unless the source ends with `;`, the REPL display of the result
function eval_example(mod, src)
    path, io = mktemp()
    got = ""
    try
        result = redirect_stdout(io) do
            r = Core.eval(mod, Meta.parseall(src))
            Libc.flush_cstdio()
            r
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
    got
end

function run_block(name, text)
    mod = block_module(name)
    ok = true
    for (src, expected) in parse_examples(text)
        got = eval_example(mod, src)
        if !matches(expected, got)
            ok = false
            println("**********************************************************************")
            println("$name\nFailed example:\n    $(replace(src, "\n" => "\n    "))\nExpected:\n    $(replace(expected, "\n" => "\n    "))\nGot:\n    $(replace(rstrip(got), "\n" => "\n    "))")
        end
    end
    ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    blocks = parse_blocks(ARGS[1])
    nex = sum(length(parse_examples(t)) for (_, t) in blocks; init=0)
    nfail = count(!, [run_block(n, t) for (n, t) in blocks])
    println("julia: $(length(blocks)) blocks, $nex examples, $nfail failures")
    exit(nfail)
end
