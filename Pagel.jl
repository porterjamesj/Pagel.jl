# wrap Phylogenetics.jl in my own API

type PhyloNode
    label::String
    length::Float64
    children::Vector{PhyloNode}
end

getroot(p::Phylogeny) = getroot(p.edge[:,2],p.edge[:,1])

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


istip(p::PhyloNode) = p.children == []


# Given a vector of rates, construct a rate matrix
#
# if model is dependant rates should be of length 8, as follows:
#
# 2->1, 3->1, 1->2, 4->2, 1->3, 4->3, 2->4, 3->4
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
    @show rates = reshape(r,(2,4))
    # copy rates over
    R[2,1], R[3,1] = rates[:,1]
    R[1,2], R[4,2] = rates[:,2]
    R[1,3], R[4,3] = rates[:,3]
    R[2,4], R[3,4] = rates[:,4]
    # compute no-change rates
    R[1,1] = -(R[2,1] + R[3,1])
    R[2,2] = -(R[1,2] + R[4,2])
    R[3,3] = -(R[1,3] + R[4,3])
    R[4,4] = -(R[2,4] + R[3,4])
    return R
end

gen_rate_matrix(r:Vector, model::Symbol) = gen_rate_matrix!(r,zeros(4,4),model)


# Perform a postorder reduction (Felsenstein's algorithm)
# of a phylogram
function postorder(f::Function, r::Function, node::PhyloNode)
    if istip(node)
        return f(node)
    else
        pos = [postorder(f,r,kid) for kid in node.children]
        push!(pos,f(node))
        return reduce(r,pos)
    end
end

# given a vector of rates, a tree, and a dict of character states,
# compute the likelihood
#
function likelihood(rates::Vector,tree::Phylogram,states::Dict)
end
