using Test
using Random
include("./BoolNets.jl")
using .BoolNets

# structural checks every generated network must pass
function check_structure(net::BoolNet, N::Int)
    @test length(net.inputs) == N
    @test length(net.tables) == N
    for i in 1:N
        ins = net.inputs[i]
        @test all(1 .<= ins .<= N)                       # inputs are valid genes
        @test length(unique(ins)) == length(ins)         # no duplicate inputs
        @test length(net.tables[i]) == 2^length(ins)     # table size matches K
        @test all(t -> t in (0, 1), net.tables[i])       # bits only
    end
end

meandeg(net) = sum(length, net.inputs) / length(net.inputs)

@testset "nextstate" begin
    # docstring example, checked against the boolean formulas directly
    net = BoolNet([[1,2],[1,3],[1,2]],
                  [[1,1,1,0],[0,0,0,1],[1,0,0,0]])
    for s in 0:7
        A, B, C = getbit(s,2), getbit(s,1), getbit(s,0)
        expected = (1-(A&B))*4 + (A&C)*2 + (1-(A|B))*1   # !A|!B, A&C, !A&!B
        @test nextstate(s, net) == expected
    end
end

@testset "findcycle / getbasins" begin
    @test findcycle(0, s -> mod(s+1, 4)) == (0, 4, [0,1,2,3])
    @test findcycle(0, s -> s == 3 ? 2 : s+1) == (2, 2, [2,3])

    net = BoolNet([[1,2],[1,3],[1,2]],
                  [[1,1,1,0],[0,0,0,1],[1,0,0,0]])
    @test sum(values(getbasins(net))) == 2^3

    rnet = randomnet(8, 2, MersenneTwister(1))
    @test sum(values(getbasins(rnet))) == 2^8
end

@testset "randomnet" begin
    net = randomnet(20, 3, MersenneTwister(1))
    check_structure(net, 20)
    @test all(length.(net.inputs) .== 3)
    @test_throws ArgumentError randomnet(5, 6)
    @test_throws ArgumentError randomnet(5, 0)
end

@testset "erdos_renyi (N, p)" begin
    check_structure(erdos_renyi(30, 0.1, MersenneTwister(1)), 30)

    # mean in-degree should be near p*N (bias from the self-input fix is negligible here)
    N, p = 200, 0.05
    rng = MersenneTwister(2)
    avg = sum(meandeg(erdos_renyi(N, p, rng)) for _ in 1:20) / 20
    @test isapprox(avg, p*N; atol=0.5)

    # extremes
    @test all(length.(erdos_renyi(10, 1.0, MersenneTwister(1)).inputs) .== 10)
    @test all(length.(erdos_renyi(10, 0.0, MersenneTwister(1)).inputs) .== 1)  # all self-loop fallback

    @test_throws ArgumentError erdos_renyi(10, 1.5)
    @test_throws ArgumentError erdos_renyi(10, -0.1)
end

@testset "erdos_renyi (N, m)" begin
    N, m = 20, 300
    net = erdos_renyi(N, m, MersenneTwister(1))
    check_structure(net, N)
    total = sum(length, net.inputs)
    @test m <= total <= m + N          # self-input fix can only add edges

    # small m: lots of empty nodes get a self-loop
    net0 = erdos_renyi(10, 0, MersenneTwister(1))
    check_structure(net0, 10)
    @test all(net0.inputs[i] == [i] for i in 1:10)

    @test_throws ArgumentError erdos_renyi(10, 101)   # > N^2
    @test_throws ArgumentError erdos_renyi(10, -1)
end

@testset "erdos_renyi_k" begin
    N, k = 200, 5.0
    rng = MersenneTwister(3)
    avg = sum(meandeg(erdos_renyi_k(N, k, rng)) for _ in 1:20) / 20
    @test isapprox(avg, k; atol=0.3)
end

@testset "watts_strogatz" begin
    N, k = 30, 4

    # p = 0 gives the plain ring lattice
    ring = watts_strogatz(N, k, 0.0, MersenneTwister(1))
    check_structure(ring, N)
    for i in 1:N
        @test Set(ring.inputs[i]) == Set(mod1(i+o, N) for o in (-2,-1,1,2))
    end

    # rewiring never changes in-degree
    net = watts_strogatz(N, k, 0.3, MersenneTwister(1))
    check_structure(net, N)
    @test all(length.(net.inputs) .== k)

    # p = 1 should actually change the wiring
    full = watts_strogatz(N, k, 1.0, MersenneTwister(1))
    @test any(Set(full.inputs[i]) != Set(ring.inputs[i]) for i in 1:N)

    @test_throws ArgumentError watts_strogatz(N, 3, 0.1)
    @test_throws ArgumentError watts_strogatz(N, k, 1.5)
end

@testset "reproducibility" begin
    a = erdos_renyi(15, 0.2, MersenneTwister(42))
    b = erdos_renyi(15, 0.2, MersenneTwister(42))
    @test a.inputs == b.inputs && a.tables == b.tables

    a = erdos_renyi(15, 40, MersenneTwister(42))
    b = erdos_renyi(15, 40, MersenneTwister(42))
    @test a.inputs == b.inputs && a.tables == b.tables
end