module BoolNets

using Random

export BoolNet, getbit, nextstate, findcycle, getbasins, randomnet


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
`B = {{A, B, C}, {A' = !A or !B, B' = A and C, C' = !A and !B}}`
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
getbit(num, pos) = num>>pos %2


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
            pos += getbit(state, N - i) * 2^pow
            pow -= 1
        end
        res += net.tables[gene][pos] * 2^(N - gene)
    end
    return res
end


"""
    findcycle(state, statefunc) -> (transient length, cycle length, cycle states)

Follow the trajectory from `state` via `statefunc` until a cycle is found.
The cycle is returned rotated so it always starts with its smallest state.
"""
function findcycle(state::Int, statefunc)
    step = 0
    seen = Dict{Int, Int}()
    states = Vector{Int}()
    while true
         if haskey(seen, state)
            t0 = seen[state]
            cyclelength = step - t0
            cyclestates = states[t0+1:end]
            cycle = circshift(cyclestates, 1-argmin(cyclestates))
            return (t0, cyclelength, cycle)
        else
            seen[state] = step 
            push!(states, state)
            state = statefunc(state)
            step += 1
        end 
    end
end


"""
    getbasins(net::BoolNet) -> Dict{attractor => basin size}
    getbasins(N::Int, f) -> Dict{attractor => basin size}

Return a `Dict{Vector{Int}, {Int}}` containing attractors and their basin sized of `net` 
"""
function getbasins(N::Int, f)
    basins = Dict{Vector{Int}, Int}()
    for k in 0:2^N-1
        _, _, c = findcycle(k, f)
        basins[c] = get(basins, c, 0) + 1
    end
    return basins
end

function getbasins(net::BoolNet)
    N = length(net.inputs)
    return getbasins(N, s->nextstate(s, net))
end

"""
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

end # module