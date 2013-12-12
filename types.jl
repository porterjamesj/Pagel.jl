import Cartesian.@forcartesian

# stub type so we can dispatch on model
immutable Model{M} end

# Transition represents a change in a state
immutable Transition
    allowed::Bool # is this an allowable transition
    trait::Int  # the trait index in the state tuple
    # the traits we are going from and to
    from::Int
    to::Int
end

Base.show(io::IO,t::Transition) =
    print(io,t.allowed ? (t.trait, t.from, t.to) : false)

isallowed(t::Transition) = t.allowed

# RateMatrix
#
# a Matrix along with a Transitions array that describes
# which indices are valid transitions; and, if so which
# transition they represent
#
type RateMatrix
    data::Matrix
    transitions::Matrix{Transition}
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
            lookup[i] = tuple(c...)
            i += 1
        end
        lookup = sort(lookup)  # TODO: This is janky as all hell

        # from the lookup table,construct the transitions
        transitions = Array(Transition,size(data))
        mask = falses(size(data))
        @forcartesian c size(data) begin
            tups = (lookup[c[1]],lookup[c[2]])
            diffs = 0
            loc = 0
            for i in 1:length(tups[1])
                if tups[1][i] != tups[2][i]
                    diffs += 1
                    loc = i
                end
            end
            if diffs == 1
                transitions[c...] = Transition(true,loc,
                                               tups[1][loc]-1,
                                               tups[2][loc]-1)
            else
                transitions[c...] = Transition(false,0,0,0)
            end
        end

        new(data,transitions,smax)
    end
end

function RateMatrix(smax::(Int...))
    s = reduce(*,smax)
    RateMatrix(zeros(s,s),smax)
end

#
# construct from a preexisting vector of data using the computed mask
# model indicates the model type, since the length of the rates vector
# should be different for each
#
# model should be either dependant (i.e. identical transitions are
# allowed to have different rates depending on the values of other
# traits) or independant (i.e. the invariant that identical
# transitions have the same rate is enforced) model. To illustrate the
# difference, the transitions:
#
# (0,0,0) => (0,0,1) and (0,1,0) => (0,1,1)
#
# must have identical rates in an independant model, but must have
# different rates in a dependant model
#
#
# The rule we use for mapping the rates vector onto the RateMatrix
# is as follows:
# If the model is dependant, it's relatively simple. We figure out
# which indices of the RateMatrix correspond to allowed transitions,
# call these valid indices. Iterate over the rates vector and copy
# rates[i] to RateMatrix[valid[i]]
#
# If the model is independant, on the other hand, things are a bit
# more complicated. In this case. we walk over the valid indices of
# the RateMatrix, checking if we have seen the corrsponding transition
# before. If not, associate this transition with the current rate
# index. If so look at which rate index was associated
# with the first instance of this transition, and copy the corresponding
# rate to this location in the RateMatrix
#


function RateMatrix(model::Model{:independant},
                    rates::Vector,
                    smax::(Int...))
    r = RateMatrix(smax)
    valid = find(isallowed,r.transitions)

    if length(rates) != length(unique(r.transitions[valid]))
        error("Rate vector not compatible with possible states")
    end

    firsts = Dict{Transition,Int}()
    j = 1
    for i in 1:length(valid)
        t = r.transitions[valid[i]]
        if haskey(firsts,t)
            r.data[valid[i]] = rates[firsts[t]]
        else
            r.data[valid[i]] = rates[j]
            firsts[t] = j
            j += 1
        end
    end

    # fill in the determined indices along the diagonal
    for i in 1:size(r.data,1)
        r.data[i,i] = -sum(r.data[i,:])
    end

    return r
end

function RateMatrix(model::Model{:dependant},
                    rates::Vector,
                    smax::(Int...))
    r = RateMatrix(smax)
    valid = find(isallowed,r.transitions)

    if length(rates) != length(valid)
        error("Rate vector not compatible with possible states")
    end

    for i in 1:length(valid)
        r.data[valid[i]] = rates[i]
    end

    # fill in the determined indices along the diagonal
    for i in 1:size(r.data,1)
        r.data[i,i] = -sum(r.data[i,:])
    end

    return r
end


# TipState

# expand a vector of states to account for missing data
function expandstate(c::Int,states::Vector{(Int...)},smax::(Int...))
    # base case
    if c == length(states[1])+1  # this seems jank (because of the indexing)
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

expandstate(states,smax) = expandstate(1,states,smax)

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


# define getindex with a tip state
Base.getindex(m::Matrix,t::TipState) = m[t.is]
Base.getindex(m::Matrix,t1::TipState,t2::TipState) = m[t1.is,t2.is]
