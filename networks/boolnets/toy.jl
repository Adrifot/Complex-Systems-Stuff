function getbit(num::Int, pos::Int)
    num>>pos % 2
end

function findcycle(state::Int, statefunc)
    step = 0
    seen = Dict{Int, Int}()
    while true
        if haskey(seen, state)
            t0 = seen[state]
            cl = step - t0
            return (t0, cl)
        else
            seen[state] = step 
            state = statefunc(state)
            step += 1
        end 
    end
end