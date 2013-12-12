# test a trivial example
tree = "$testdir/simple.tree" |> open |> readchomp |> Pagel.parsenewick
states = Pagel.statedict("$testdir/simple.txt")
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
