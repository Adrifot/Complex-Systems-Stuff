struct State
    A::Bool
    B::Bool
    C::Bool
end

function nextstate(curr::State)
    State(curr.C, curr.A, !curr.B)
end

function nextstate(curr::Int)
    (curr%2)*4 + ((curr>>2)%2)*2 + ((curr>>1)%2)⊻1
end

function showstate(s::State)
    print(s.A |> Int, s.B |> Int, s.C |> Int)
end

function iteratestates(s::State, k::Int)
    if k == 0 return
    else 
        showstate(s)
        print("→")
        iteratestates(nextstate(s), k-1)
    end
end

start = State(false, false, false)

function state2int(s::State)
    return Int(s.A) * 4 + Int(s.B) * 2 + Int(s.C) 
end

function findcycle(start::State)
    s = state2int(start)
    step = 0
    seen = Dict{Int, Int}()
    while true
        if haskey(seen, s)
            t0 = seen[s]
            cl = step - t0
            println("Cycle of length $cl found after a transient length of $t0")
            return (t0, cl)
        else
            seen[s] = step 
            s = nextstate(s)
            step += 1
        end 
    end
end