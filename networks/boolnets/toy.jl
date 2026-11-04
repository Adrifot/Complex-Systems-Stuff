using Test

function getbit(num::Int, pos::Int)
    num>>pos % 2
end

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

function getbasins(N::Int, f)
    basins = Dict{Vector{Int}, Int}()
    for k in 0:2^N-1
        _, _, c = findcycle(k, f)
        basins[c] = get(basins, c, 0) + 1
    end
    return basins
end

function f1(s::Int) 
    getbit(s, 0)*4 + getbit(s, 2)*2 + (getbit(s, 1) ⊻ 1)
end

function f2(s::Int)
    (getbit(s, 2) | getbit(s, 1))*4 + (getbit(s, 2) & (getbit(s, 0) ⊻ 1))*2 + (getbit(s, 2) & getbit(s, 1))
end

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
