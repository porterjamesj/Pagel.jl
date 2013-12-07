import Cartesian.@forcartesian

#
# RateMatrix handles converting between tuples
# describing the state of a species and the index into
# the rate matrix corresponding to that species
#
immutable RateMatrix
    data::Matrix
    states::(Int...)
    jump::(Int...)
    lookup::Vector

    # constructor
    function RateMatrix(data::Matrix,states::(Int...))
        # data must be square
        if size(data,1) != size(data,2)
            error("RateMatrix must be square")
        end
        # dimensions suggested by states tuple must
        # match actual dimensionality of data matrix
        if reduce(*,states) != size(data,1)
            error("RateMatrix states tuple is incompatible with data matrix.")
        end

        # construct jump table
        tmp = Array(Int,length(states))
        for i in 1:length(states)
            tmp[i] = reduce(*,states[i+1:end])
        end
        jump = tuple(tmp...)

        # construct the lookup table
        lookup = Array((Int...),reduce(*,states))
        i = 1
        @forcartesian c states begin
            lookup[i] = tuple(reverse(c)...)
            i += 1
        end

        new(data,states,jump,lookup)
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
