# test a trivial example
tree = "$testdir/simple.tree" |> open |> readchomp |> Pagel.parsenewick
states = Pagel.statedict("$testdir/simple.txt")

# dependant
Q = Pagel.RateMatrix(Pagel.Model{:dependant}(),
                     [113.463768,
                      156.880815,
                      2.057176,
                      0.000000,
                      0.000000,
                      2.023122,
                      209.392213,
                      81.301830],
                     states.smax)

@test_approx_eq log(Pagel.likelihood(tree,Q,states)) -0.017479829591899487

#independant
Q = Pagel.RateMatrix(Pagel.Model{:independant}(),
                     [6.758346, 13.969125, 6.758346, 13.969126],
                     states.smax)

@test_approx_eq log(Pagel.likelihood(tree,Q,states)) -1.3862943611217102

# weird parameter set: [56.424218003678135,143.81537923337802,0.0,0.0005864443949624926,0.0,0.0017145650558412413,0.3503027852542958,0.18811228041188838]
