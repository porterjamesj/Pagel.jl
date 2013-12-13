# this is necessary for me on Mavericks;
# on <= 10.8 it may break things, I am not sure
@osx_only dlopen("libstdc++",RTLD_GLOBAL)
using NLopt

# use NLOpt to maximize the likelihood of the given model
# for the given tree and states.
#
# the user can optionally specify an algorithm
function maxlikelihood(model::Symbol,tree::PhyloNode,states::StateDict;
                       algorithm=:LN_NELDERMEAD,ftol_rel=1e-6)
    # first we need to determine the dimensionality of the problem
    # TODO: this is wonky; there's a better way
    r = RateMatrix(states.smax)
    trans = r.transitions[find(isallowed,r.transitions)]
    if model == :dependant
        dims = length(trans)
    else # algorithm == :independant
        dims = trans |> unique |> length
    end

    opt = Opt(algorithm, dims)

    # the function we will optimize
    # grad is necessary for the NLOpt interface,
    # we don't use it since we can't take derivatives
    function likelihood(params::Vector,grad::Vector)
        Q = RateMatrix(Model{model}(),params,states.smax)
        log(Pagel.likelihood(tree,Q,states))
    end

    # we want to maximize the log likelihood
    max_objective!(opt, likelihood)

    # rates should be greater than zero
    lower_bounds!(opt::Opt, zeros(dims))
    # upper_bounds!(opt::Opt, fill(300,dims))

    # set stopping criteria
    ftol_rel!(opt, ftol_rel)

    # initial search point; as long as the dimensionality of the problem
    # I use the same starting points for the search as Mesquite
    inits = Array(Float64,dims)
    for i in 1:dims
        inits[i] = randbool() ? 10 : 0.01
    end
    optimize(opt,inits)
end

# interface that lets you pass filenames
function maxlikelihood(model::Symbol,
                       treefile::String,statesfile::String;
                       algorithm=:LN_NELDERMEAD,ftol_rel=1e-6)
    tree = readchomp(treefile) |> parsenewick
    states = states = statedict(statesfile)
    maxlikelihood(model, tree,states,
                  algorithm=algorithm,
                  ftol_rel=ftol_rel)
end
