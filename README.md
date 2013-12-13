# Pagel.jl

A Julia implementation of the method described in:

Pagel, Mark. "Detecting correlated evolution on phylogenies: a general method for the comparative analysis of discrete characters." _Proceedings of the Royal Society of London. Series B: Biological Sciences_ 255.1342 (1994): 37-45.

With some extensions that, in principle, allow you to model arbitrary numbers
of characters and states.

#### WARNING

This is a research prototype; you would be a fool to attempt to use it.

## Use

You'll need a recent version of Julia (>=0.2.0).

    Pkg.clone("https://github.com/porterjamesj/Pagel.jl.git")
    using Pagel

The module exports a single function, `maxlikelihood`, whoose interface is:

    maxlikelihood(model,treefile,statesfile)

Where model is one of `:dependant` or `:independant`, and `treefile` and `statefile`
are your data. `treefile` should just have a single newick string in it (not a NEXUS
file), states should be pretty much exactly what you would use for BayesTraits.
