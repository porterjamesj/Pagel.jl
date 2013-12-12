module Pagel

include("newick.jl")
include("types.jl")

#
# Convert a whitespace separated file with state information into a
# states dictionary, which maps from tip labels to TipStates
#
function statedict(filepath::String)
    data = readdlm(filepath,'\t')
    dict = Dict{String,TipState}()
    maxstate = tuple(map(x->int(x+1),reducedim(max,data[:,2:end],[1],-1))...)
    for i in 1:size(data,1)
        row = [datum=="-"?-1:datum for datum in data[i,:]]
        dict[row[1]] = TipState(tuple(map(int,row[2:end])...), maxstate)
    end
    return dict
end


# compute the transition probabilites in a given interval
function P(Q::Matrix,t::Real)
    return expm(Q*t)
end

#
# recursively compute the likelihood of a single node.
#
# Q is the rate matrix, states maps from tip labels to
# an array of character states (e.g. (0,1) (1,1) etc.)
# i is an integer index into a dimension of Q
#
function likelihood(node::PhyloNode, i::Int, Q::Matrix, states::Dict)
    if istip(node)
        return 1  # base case
    else
        # array to hold results from each child
        res = Float64[]
        for child in node.children
            # compute allowed states for the child
            if istip(child)
                childindex = state2ind(states[child.label])
            else
                # the child is an internal node and allowed any state
                childindex = 1:4
            end
            push!(res,sum([P(Q,child.length)[i,c]*likelihood(child, c, Q, states)
                           for c in childindex]))
        end
        return reduce(*,res)
    end
end

likelihood(node::PhyloNode,
           Q::Matrix,
           states::Dict) = sum([likelihood(node,i,Q,states) for i in 1:4])



end # module Pagel
