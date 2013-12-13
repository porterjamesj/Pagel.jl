using Pagel

tree = "test.tree" |> open |> readchomp |> Pagel.parsenewick
states = Pagel.statedict("test.states")
Q = Pagel.RateMatrix(Pagel.Model{:independant}(),[1,2,3,4],states.smax)

# Qd = Pagel.gen_rate_matrix([0.020778,0.037929,0.009609,0.040891],:independant)
# Qd = Pagel.gen_rate_matrix([0.020778,0.037929,0.009609,0.06],:independant)
