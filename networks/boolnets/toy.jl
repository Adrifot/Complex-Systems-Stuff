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

function getattractors(N::Int, f)
    attractors = Set{Vector{Int}}()
    for k in 0:2^N-1
        _, _, att = findcycle(k, f)
        push!(attractors, att)
    end
    return attractors
end

function f1(s::Int) 
    getbit(s, 0)*4 + getbit(s, 2)*2 + (getbit(s, 1) ⊻ 1)
end

function f2(s::Int)
    (getbit(s, 2) | getbit(s, 1))*4 + (getbit(s, 2) & (getbit(s, 0) ⊻ 1))*2 + (getbit(s, 2) & getbit(s, 1))
end