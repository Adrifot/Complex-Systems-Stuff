module BoolNets

using Random
using JSON3

export BoolNet, getbit, nextstate,  getbasins, 
        randomnet, erdos_renyi, erdos_renyi_k, watts_strogatz, 
        plot_topology, plot_statediagram, save_graph, analyze,
        getattractors, savenet, loadnet, validnet,
        kin, kout, avg_kin, avg_kout

"""
    BoolNet(inputs, tables)

A Boolean network with N genes, numbered 1...N.

# Fields
- `inputs::Vector{Vector{Int}}`: the input dependencies of each gene.
    `inputs[i]` contains the genes that affect gene `i`. The **first**
    listed input is the most significant bit of the boolean table position.
- `tables::Vector{Vector{Int}}`: truth tables for each gene.
    `tables[i]` contains `2^K` output bits for gene `i` ordered
    by position in the truth table.

# Example
For the following Boolean Network:\n 
`G = {{A, B, C}, {A' = !A or !B, B' = A and C, C' = !A and !B}}`
The `BoolNet` struct will contain:\n
- `inputs`: `[1,2], [1,3], [1,2]`
- `tables`: `[[1,1,1,0], [0,0,0,1], [1,0,0,0]]`
"""
struct BoolNet
    inputs::Vector{Vector{Int}}
    tables::Vector{Vector{Int}}
end


"""
    getbit(num::Int, pos::Int)

Get bit `pos` of `num`. Bit 0 is the least significant.
"""
getbit(num, pos) = (num>>pos)%2


"""
    nextstate(state::Int, net::BoolNet) -> Int

Return the next state of the boolean network `net` given current state `state`.
"""
function nextstate(state::Int, net::BoolNet)
    N = length(net.inputs)
    res = 0
    for gene in 1:N 
        pos = 1
        pow = length(net.inputs[gene]) -1
        for i in net.inputs[gene]
            pos += getbit(state, N - i) * 1<<pow
            pow -= 1
        end
        res += net.tables[gene][pos] * 2^(N - gene)
    end
    return res
end

"""
    getbasins(N, f) -> Dict{Int, Int}
    getbasins(net) -> Dict{Int, Int}

Get the basins of each attractor for a given boolean network `net`. Returns a dictionary with `AttractorID => BasinSize` entries.
Runs in `O(2^N)` time. Consumes `O(2^N)` space. Mind the size of the boolean network.
"""
function getbasins(N::Int, f)
    N > 25 && throw(ArgumentError("N is too big! This function will try to assign 2^$N Ints! Try again with an N ≤ 25"))
    v = zeros(Int, 1<<N)
    path = Int[]

    for s0 in 0:(1<<N)-1
        v[s0+1] != 0 && continue
        empty!(path)
        s = s0
        label = 0

        while true
            val = v[s+1]
            if val > 0
                label = val 
                break
            elseif val < 0
                label = s + 1
                break
            else
                push!(path, s)
                v[s+1] = -1
                s = f(s)
            end
        end

        for t in path
            v[t+1] = label
        end
    end

    basins = Dict{Int, Int}()
    for l in v
        basins[l-1] = get(basins, l-1, 0) + 1
    end
    return basins
end

getbasins(net::BoolNet) = getbasins(length(net.inputs), s -> nextstate(s, net)) 

"""
    randomnet(N, K, rng) -> BoolNet

Generate a random Boolean network using the NK model
with `N` genes where each gene receives input from `K` distinct genes.
"""
function randomnet(N::Int, K::Int, rng::AbstractRNG=Random.default_rng())
    1 ≤ K ≤ N || throw(ArgumentError("K must be between 1 and N"))

    inputs = Vector{Vector{Int}}(undef, N)
    tables = Vector{Vector{Int}}(undef, N)

    num_entries = 1 << K # = 2^K

    for i in 1:N 
        # Randomly choose K distinct genes and populate the truth table
        inputs[i] = randperm(rng, N)[1:K]
        tables[i] = rand(rng, 0:1, num_entries)
    end

    return BoolNet(inputs, tables)
end


"""
    erdos_renyi(N, p, rng; force_selfloops=false) -> BoolNet
    erdos_renyi(N, m, rng; force_selfloops=false) -> BoolNet

Generate an Erdos-Renyi-like boolean network of `N` nodes with link probability `p`.
OR 
Generate an Erdos-Renyi-like boolean network of `N` nodes with `m` edges.

Self-loops will be created for each isolated node if `force_selfloops` is set to `true`.
Naturally occuring self-loops are allowed.
"""
function erdos_renyi(N::Int, p::Float64, rng::AbstractRNG=Random.default_rng(); force_selfloops::Bool=false)
    0.0 ≤ p ≤ 1.0 || throw(ArgumentError("p is a probability and should be between 0 and 1."))

    inputs = Vector{Vector{Int}}(undef, N)
    tables = Vector{Vector{Int}}(undef, N)

    for node in 1:N
        gene_inputs = [i for i in 1:N if rand(rng) < p]
        (force_selfloops && isempty(gene_inputs)) && (gene_inputs = [node]) # if no connections, make self-loop
        inputs[node] = gene_inputs
        K = length(gene_inputs)
        tables[node] = rand(rng, 0:1, 1<<K)
    end

    return BoolNet(inputs, tables)
end

function erdos_renyi(N::Int, m::Int, rng::AbstractRNG=Random.default_rng(); force_selfloops::Bool=false)
    0 ≤ m ≤ N^2 || throw(ArgumentError("m is too big. With the provided N = $N, maximum value is m = $(N^2)"))

    inputs = Vector{Vector{Int}}(undef, N)
    tables = Vector{Vector{Int}}(undef, N)

    inputs = [Int[] for _ in 1:N]
    tables = Vector{Vector{Int}}(undef, N)

    # pick m distinct edges out of the N^2 possible ones
    for idx in Random.randperm(rng, N^2)[1:m]
        q, r = divrem(idx-1, N)
        node = q + 1
        input = r + 1
        push!(inputs[node], input)
    end

    for node in 1:N
        (force_selfloops && isempty(inputs[node])) && (push!(inputs[node], node)) # if no connections, make self-loop
        sort!(inputs[node])
        K = length(inputs[node])
        tables[node] = rand(rng, 0:1, 1<<K)
    end
    return BoolNet(inputs, tables)
end


"""
    erdos_renyi_k(N, k, rng) -> BoolNet

Generate an Erdos-Renyi-like boolean network with a desired average in-degree `k`
"""
function erdos_renyi_k(N::Int, k::Real, rng::AbstractRNG=Random.default_rng())
    p = k / N
    return erdos_renyi(N, p, rng)
end


"""
    watts_strogatz(N, k, p, rng) -> BoolNet

Generate a Watts-Strogatz small-world random boolean network. 
Start from a circular 1D lattice of `N` nodes, each connected to its `k` neighbors.
Then rewire each link with a probability of `p`.
"""
function watts_strogatz(N::Int, k::Int, p::Float64, rng::AbstractRNG=Random.default_rng())
    iseven(k) || throw(ArgumentError("Initial degree k must be even."))
    k < N || throw(ArgumentError("`k` has to be smaller than `N` in order to avoid duplicate links."))
    0.0 ≤ p ≤ 1.0 || throw(ArgumentError("p is a probability; it should be between 0 and 1."))

    # Initial ring lattice
    inputs = [Int[] for _ in 1:N]
    halfk = k ÷ 2

    for i in 1:N
        for offset in -halfk:halfk
            if offset != 0
                target = mod1(i + offset, N)
                push!(inputs[i], target)
            end
        end
    end

    # Rewire the edges
    for i in 1:N
        for idx in 1:length(inputs[i])
            if rand(rng) < p
                available = setdiff(1:N, inputs[i])
                if !isempty(available) inputs[i][idx] = rand(rng, available) end
            end
        end
    end

    # Generate truth tables
    tables = [rand(rng, 0:1, 1 << length(inputs[i])) for i in 1:N]

    return BoolNet(inputs, tables)
end



"""
    getattractors(net) -> Vector{Tuple{Vector{Int}, Int}}

Return `(cycle_states, basin_size)` for every attractor of `net`, largest basin first.
"""
function getattractors(net::BoolNet)
    N = length(net.inputs)
    f = s -> nextstate(s, net)
    result = Tuple{Vector{Int}, Int}[]
    for (label, basin) in getbasins(N, f)
        start = label      
        cycle = [start]
        s = f(start)
        while s != start
            push!(cycle, s)
            s = f(s)
        end
        push!(result, (cycle, basin))
    end
    sort!(result, by = x -> -x[2])
    return result
end

fmt_state(s::Int, N::Int) = string(s, base=2, pad=N)


"""
    plot_topology(net::BoolNet) -> String

DOT source for the gene-interaction graph. Edges point from input genes to target genes.
"""
function plot_topology(net::BoolNet)
    N = length(net.inputs)
    io = IOBuffer()
    println(io, "digraph BoolNetTopology {")
    println(io, "  node [shape=circle, style=filled, fillcolor=\"#E3F2FD\", fontname=\"Helvetica\"];")
    println(io, "  edge [color=\"#555555\", arrowhead=vee];")
    for gene in 1:N, input_gene in net.inputs[gene]
        println(io, "  $input_gene -> $gene;")
    end
    println(io, "}")
    return String(take!(io))
end


"""
    plot_statediagram(net::BoolNet; maxN=10) -> String

DOT source for the state-transition diagram over all 2^N states.
Attractor states are highlighted in red. Refuses nets with `N > maxN`.
"""
function plot_statediagram(net::BoolNet; maxN::Int=10)
    N = length(net.inputs)
    N <= maxN || throw(ArgumentError("State diagram has 2^$N nodes; raise maxN if you really want it."))

    onattractor = Set(Iterators.flatten(first.(getattractors(net))))

    io = IOBuffer()
    println(io, "digraph StateSpace {")
    println(io, "  rankdir=LR;")
    println(io, "  node [fontname=\"Courier\", shape=box, style=rounded];")

    for s in 0:(1 << N) - 1
        if s in onattractor
            println(io, "  \"$(fmt_state(s, N))\" [style=\"rounded,filled\", fillcolor=\"#FFCDD2\", color=\"#D32F2F\", penwidth=2.0];")
        else
            println(io, "  \"$(fmt_state(s, N))\" [style=\"rounded,filled\", fillcolor=\"#F5F5F5\", color=\"#BDBDBD\"];")
        end
    end

    for s in 0:(1 << N) - 1
        ns = nextstate(s, net)
        color = s in onattractor ? "#D32F2F" : "#757575"
        println(io, "  \"$(fmt_state(s, N))\" -> \"$(fmt_state(ns, N))\" [color=\"$color\"];")
    end

    println(io, "}")
    return String(take!(io))
end


function save_graph(dot_str::String, filename::String="graph.png")
    Sys.which("dot") === nothing &&
        error("Graphviz `dot` not found. Install it (e.g. `sudo apt install graphviz`).")
    ext = lstrip(splitext(filename)[2], '.')
    dir = dirname(filename)
    isempty(dir) || mkpath(dir)

    open(`dot -T$ext -o $filename`, "w", stdout) do io
        print(io, dot_str)
    end
    return filename
end


"""
    analyze(net, name; outdir="img", maxN=10)

Print the attractors of `net` and save its topology (and, if `N ≤ maxN`,
its state diagram) as PNGs in `outdir`. Returns the attractor list.
"""
function analyze(net::BoolNet, name::String; outdir::String="img", maxN::Int=10)
    N = length(net.inputs)
    attractors = getattractors(net)

    println("$(length(attractors)) attractor(s):")
    for (cycle, basin) in attractors
        println("  length $(length(cycle)), basin $basin: ", join(fmt_state.(cycle, N), " → "))
    end

    save_graph(plot_topology(net), joinpath(outdir, "$(name)_topology.png"))
    if N <= maxN
        save_graph(plot_statediagram(net; maxN), joinpath(outdir, "$(name)_state_diagram.png"))
    else
        @info "Skipping state diagram: N = $N > maxN = $maxN"
    end
    return attractors
end

"""
    validnet(net::BoolNet) -> Bool

Structural check: one table per gene, table length `2^K`, inputs are genes in `1:N`,
table entries are 0 or 1.
"""
function validnet(net::BoolNet)
    N = length(net.inputs)
    length(net.tables) == N || return false
    for i in 1:N
        length(net.tables[i]) == 2^length(net.inputs[i]) || return false
        all(g -> 1 <= g <= N, net.inputs[i]) || return false
        all(b -> b == 0 || b == 1, net.tables[i]) || return false
    end
    return true
end

"""
    savenet(path, net::BoolNet)

Write `net` as JSON. Genes are numbered 1..N in the file.
"""
function savenet(path::AbstractString, net::BoolNet)
    validnet(net) || throw(ArgumentError("refusing to save an invalid BoolNet"))
    doc = (format = "boolnet", version = 1, n = length(net.inputs),
           inputs = net.inputs, tables = net.tables)
    open(path, "w") do io
        JSON3.pretty(io, JSON3.write(doc))
    end
    return path
end


"""
    loadnet(path) -> BoolNet

Read a JSON file written by `savenet`. Throws `ArgumentError` for anything malformed.
"""
function loadnet(path::AbstractString)
    obj = JSON3.read(read(path, String))
    for key in (:format, :version, :n, :inputs, :tables)
        haskey(obj, key) || throw(ArgumentError("missing field `$key`"))
    end
    obj.format == "boolnet" || throw(ArgumentError("not a boolnet file"))
    obj.version == 1 || throw(ArgumentError("unsupported version $(obj.version)"))

    inputs = [Vector{Int}(x) for x in obj.inputs]
    tables = [Vector{Int}(x) for x in obj.tables]
    obj.n == length(inputs) || throw(ArgumentError("`n` disagrees with `inputs`"))

    net = BoolNet(inputs, tables)
    validnet(net) || throw(ArgumentError("file describes an invalid network"))
    return net
end


function kin(net::BoolNet)
    kins = [length(ins) for ins in net.inputs]
    return kins
end

function avg_kin(net::BoolNet) 
    kins = kin(net)
    return sum(kins) / length(kins) 
end

function kout(net::BoolNet)
    N = length(net.inputs)
    kouts = zeros(Int, N)
    for input in net.inputs
        for gene in input
            kouts[gene] += 1
        end
    end
    return kouts
end

function avg_kout(net::BoolNet)
    kouts = kout(net)
    return sum(kouts) / length(kouts)
end

end # module

# using .BoolNets

