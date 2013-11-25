module Pagel

using Phylogenetics

# wrap Phylogenetics.jl in my own API

type PhyloNode
    label::String
    length::Float64
    children::Vector{PhyloNode}
end

getroot(p::Phylogeny) = Phylogenetics.getroot(p.edge[:,2],p.edge[:,1])

# produce a nested array of arrays of the same shape as getkids,
# but which contains edge lengths instead
function getlengths(phy::Phylo)
    N = length(phy.tipLabel)
    lengths = [Float64[] for i in 1:N+phy.Nnode]
    for i in 1:N+phy.Nnode-1
        push!(lengths[phy.edge[i,1]],phy.edgeLength[i])
    end
    return lengths
end

# turn a Phylo into a tree of PhyloNodes
function wrap(p::Phylo)

    # given the index of a node and it's length, wrap it
    function wrapnode(n::Int,l::Float64)
        return PhyloNode(labels[n],l,
                         convert(Vector{PhyloNode},
                                 [wrapnode(n,l) for (n,l) in
                                  zip(kids[n],lengths[n])]))
    end

    kids = getkids(p)
    lengths = getlengths(p)
    root = getroot(p)
    labels = [p.tipLabel, p.nodeLabel]
    return wrapnode(root,0.0)
end

# convenience function to read in a tree and wrap it
function wrap(filepath::String)
    return wrap(readtree(filepath)[1])
end


istip(p::PhyloNode) = p.children == []


# Given a vector of rates, construct a rate matrix
#
# if model is dependant rates should be of length 8, as follows:
#
# q12   q13   q21   q24   q31   q34   q42   q43
#
# if model is independant, rates should be of length 4, as follows:
#
# 1->3, 3->1, 1->2, 2->1
#
# this is equivalent to the forward and backward rates in the first
# character, followed by the forward and backward rates in the second
# character
#
# this uses the memory already allocated in R
# R should be 4 by 4
function gen_rate_matrix!(r::Vector,
                          R::Matrix,
                          model::Symbol)
    if model == :independant
        # extend the rates, mirroring across antidiagonal
        r = [ r[4], r[2], r[3], r[2], r[1], r[4], r[1], r[3]]
    end
    # each column corresonds to a column in the rate matrix
    # copy rates over
    R[2,1] = r[3]
    R[3,1] = r[5]
    R[1,2] = r[1]
    R[4,2] = r[7]
    R[1,3] = r[2]
    R[4,3] = r[8]
    R[2,4] = r[4]
    R[3,4] = r[6]
    # compute no-change rates
    R[1,1] = -(R[1,2] + R[1,3])
    R[2,2] = -(R[2,1] + R[2,4])
    R[3,3] = -(R[3,1] + R[3,4])
    R[4,4] = -(R[4,2] + R[4,3])
    return R
end

gen_rate_matrix(r::Vector, model::Symbol) = gen_rate_matrix!(r,zeros(4,4),model)


# convert a state tuple to an index into the
# rate matrix
function state2ind(state::(Int,Int))
    if state == (0,0)
        return 1
    elseif state == (0,1)
        return 2
    elseif state == (1,0)
        return 3
    elseif state == (1,1)
        return 4
    else
        error("state is invald")
        # TODO: this is super janky, need to implement
        # it for missing data (-1) as well.
    end
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

#
# Convert a whitespace separated file with state information into a
# states dictionary. A -1 indicates missing data
#
# NB: It would probably be better to use NA for missing data,
# but I'm not sure I want to have all of DataFrames as a dependancy
#
function statedict(filepath::String)
    data = readdlm(filepath,'\t')
    dict = Dict{String,(Int,Int)}()
    for i in 1:size(data,1)
        row = [datum=="-"?-1:datum for datum in data[i,:]]
        dict[row[1]] = (int(row[2]),int(row[3]))
    end
    return dict
end

end # module Pagel
