import Cartesian.@forcartesian

#
# RateMatrix handles converting between tuples
# describing the state of a species and the index into
# the rate matrix corresponding to that species
#
immutable RateMatrix
    data::Matrix
    mask::BitMatrix
    smax::(Int...)

    # constructor
    function RateMatrix(data::Matrix,smax::(Int...))
        # data must be square
        if size(data,1) != size(data,2)
            error("RateMatrix must be square")
        end
        # dimensions suggested by smax tuple must
        # match actual dimensionality of data matrix
        if reduce(*,smax) != size(data,1)
            error("RateMatrix smax tuple is incompatible with data matrix.")
        end

        # construct temporary lookup table
        lookup = Array((Int...),reduce(*,smax))
        i = 1
        @forcartesian c smax begin
            lookup[i] = tuple(reverse(c)...)
            i += 1
        end

        # from the lookup table,construct the mask
        mask = falses(size(data))
        @forcartesian c size(data) begin
            tups = (lookup[c[1]],lookup[c[2]])
            diffs = 0
            for i in 1:length(tups[1])
                if tups[1][i] != tups[2][i]
                    diffs += 1
                end
            end
            if diffs ==1
                mask[c...] = true
            end
        end

        new(data,mask,smax)
    end
end


function getindex(m::RateMatrix,state1::(Int...),state2::(Int...))
    index1 = 0
    for (i,v) in enumerate(state1)
        index1 += v*m.jump[i]
    end
    index2 = 0
    for (i,v) in enumerate(state2)
        index2 += v*m.jump[i]
    end
    return m.data[index1+1,index2+1]
end

# macro for delegating undefined methods to the first field of a type
# see https://github.com/JuliaLang/julia/pull/3292
macro delegate(source, targets)
    typename = esc(source.args[1])
    fieldname = esc(Expr(:quote, source.args[2].args[1]))
    funcnames = targets.args
    n = length(funcnames)
    fdefs = Array(Any, n)
    for i in 1:n
        funcname = esc(funcnames[i])
        fdefs[i] = quote
                     ($funcname)(a::($typename), args...) =
                       ($funcname)(a.($fieldname), args...)
                   end
    end
    return Expr(:block, fdefs...)
end


importall Base
@delegate RateMatrix.data [ndims size]


# TipState

expandstate(states,smax) = expandstate(1,states,smax)

# expand a vector of states to account for missing data
function expandstate(c::Int,states::Vector{(Int...)},smax::(Int...))
    # base case
    if c == length(states[1])+1  # this seems jank
        return states
    end
    tmp = (Int...)[]
    for i in 1:length(states)
        if states[i][c] == -1
            for j in 0:smax[c]-1
                push!(tmp, tuple(setindex!([states[i]...],j,c)...))
            end
        end
    end
    # base case: if no missing data; don't change anything
    expandstate(c+1,length(tmp)==0 ? states : tmp ,smax)
end


function state2int(state::(Int...),maxstates::(Int...))
    index = 0::Int
    for (i,v) in enumerate(state)
        index += v * reduce(*,maxstates[i+1:end])
    end
    return index+1  # +1 becase julia is 1 indexed
end


type TipState
    states::Vector{(Int...)}
    smax::(Int...)
    is::Vector{Int}

    # constructor from a state tuple
    function TipState(state::(Int...), smax::(Int...))

        # tuples must be same length
        if length(state) != length(smax)
            error("TipState cannot be constructed.")
        end

        # ensure state is valid
        for i in 1:length(state)
            if state[1] >= smax[i]
                error("TipState cannot be constructed.")
            end
        end

        # compute the integer vector representation of this state
        # remember that -1 means missing data, first expand these
        states = expandstate([state],smax)

        #convert each state to an integer
        is = [state2int(state,smax) for state in states]

        new(states,smax,is)
    end

    # TODO: constructor from an integer
end
