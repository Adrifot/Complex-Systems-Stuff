using Test

include("BoolNets.jl")
using .BoolNets

# ---------------------------------------------------------------------------
# Hand-written oracles and fixtures (N = 3, gene A = bit 2, B = bit 1, C = bit 0)
# ---------------------------------------------------------------------------

# K=1 network: A' = C, B' = A, C' = !B
f1(s::Int) = getbit(s, 0)*4 + getbit(s, 2)*2 + (getbit(s, 1) ⊻ 1)

# K=2 network: A' = A OR B, B' = A AND !C, C' = A AND B
f2(s::Int) = (getbit(s, 2) | getbit(s, 1))*4 +
             (getbit(s, 2) & (getbit(s, 0) ⊻ 1))*2 +
             (getbit(s, 2) & getbit(s, 1))

netk1 = BoolNet(
    [[3], [1], [2]],
    [[0,1], [0,1], [1,0]]
)

netk2 = BoolNet(
    [[1,2], [1,3], [1,2]],
    [[0,1,1,1], [0,0,1,0], [0,0,0,1]]
)

# ---------------------------------------------------------------------------
# Helpers for checks that hold for ANY network
# ---------------------------------------------------------------------------

# Structural sanity: sizes match, inputs are valid genes, table entries are bits.
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

# A cycle is closed if each state's successor is the next listed state
# (and the last one leads back to the first).
function isclosedcycle(cycle, f)
    n = length(cycle)
    return all(f(cycle[i]) == cycle[mod1(i + 1, n)] for i in 1:n)
end

# Run the invariant checks on a network.
function checkinvariants(net::BoolNet)
    N = length(net.inputs)
    f = s -> nextstate(s, net)
    basins = getbasins(net)
    @test validnet(net)
    @test sum(values(basins)) == 2^N              # every state lands somewhere
    @test all(isclosedcycle(c, f) for c in keys(basins))
end

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@testset "f1 cycles" begin
    @test findcycle(0, f1) == (0, 6, [0, 1, 5, 7, 6, 2])
    @test findcycle(6, f1) == (0, 6, [0, 1, 5, 7, 6, 2])
    @test findcycle(3, f1) == (0, 2, [3, 4])
    @test findcycle(4, f1) == (0, 2, [3, 4])
end

@testset "f2 cycles" begin
    @test findcycle(0, f2) == (0, 1, [0])
    @test findcycle(1, f2) == (1, 1, [0])
    @test findcycle(4, f2) == (0, 4, [4, 6, 7, 5])
    @test findcycle(7, f2) == (0, 4, [4, 6, 7, 5])
    @test findcycle(2, f2) == (1, 4, [4, 6, 7, 5])
end

@testset "attractors & basins (hand-written functions)" begin
    @test getbasins(3, f1) == Dict([3, 4] => 2, [0, 1, 5, 7, 6, 2] => 6)
    @test getbasins(3, f2) == Dict([0] => 2, [4, 6, 7, 5] => 6)
end

@testset "BoolNet agrees with the hand-written oracles" begin
    for s in 0:7
        @test nextstate(s, netk1) == f1(s)
        @test nextstate(s, netk2) == f2(s)
    end
    @test getbasins(netk1) == getbasins(3, f1)
    @test getbasins(netk2) == getbasins(3, f2)
end

@testset "invariants on the toy networks" begin
    checkinvariants(netk1)
    checkinvariants(netk2)
end

@testset "validnet rejects malformed networks" begin
    # one gene, but an input that is not a gene and a table of the wrong size
    @test !validnet(BoolNet([[1, 2]], [[0, 1]]))
    # table entry that is not a bit
    @test !validnet(BoolNet([[1]], [[0, 2]]))
    # fewer tables than genes
    @test !validnet(BoolNet([[1], [1]], [[0, 1]]))
end

@testset "edge cases outside N = 3" begin
    # N=1, gene negates itself: 0 -> 1 -> 0
    flip = BoolNet([[1]], [[1, 0]])
    @test getbasins(flip) == Dict([0, 1] => 2)
    @test findcycle(0, s -> nextstate(s, flip)) == (0, 2, [0, 1])
    checkinvariants(flip)

    # N=1, gene copies itself: two fixed points
    copy1 = BoolNet([[1]], [[0, 1]])
    @test getbasins(copy1) == Dict([0] => 1, [1] => 1)
    checkinvariants(copy1)

    # N=1, K=0: a constant gene that is always ON (empty input list)
    const_on = BoolNet([Int[]], [[1]])
    @test nextstate(0, const_on) == 1
    @test nextstate(1, const_on) == 1
    @test findcycle(0, s -> nextstate(s, const_on)) == (1, 1, [1])
    @test getbasins(const_on) == Dict([1] => 2)
    checkinvariants(const_on)
end

@testset "N = 4 with a K = 3 gene" begin
    # gene 1 = majority(gene 2, gene 3, gene 4); genes 2, 3, 4 copy their predecessor.
    # State bits are g1 g2 g3 g4, with gene 1 as the most significant bit.
    majnet = BoolNet(
        [[2, 3, 4], [1], [2], [3]],
        [[0,0,0,1,0,1,1,1], [0,1], [0,1], [0,1]]
    )
    # hand-computed transitions
    @test nextstate(7, majnet)  == 11   # 0111 -> 1011
    @test nextstate(1, majnet)  == 0    # 0001 -> 0000
    @test nextstate(8, majnet)  == 4    # 1000 -> 0100
    @test nextstate(5, majnet)  == 10   # 0101 -> 1010
    checkinvariants(majnet)
end