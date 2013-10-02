intrahosts
=========

source("SIR.R")
reed_frost_multi_demo()


SIR.R
-------------------------------------------------
reed_frost(I0,N,q,niters)
-------------------------------------------------
# runs original Reed Frost model (chain binomial)

I0      initial number of infected I(0)
N       size of population
q       avoidance probability
niters  number of iterations

--------------------------------------------------
reed_frost_multi(I0, N, NHosts, q0, q1)
--------------------------------------------------
# runs Reed Frost model across multiple hosts (chain multinomial)

I0      initial number of infected (in Host 1)
N       size of population in each host
NHosts  number of hosts
q0      avoidance probability within host
q1      avoidance probability between hosts

Returns I, niters by NHosts matrix where the rows is the viral population in the hosts over time.

---------------------------------------------------
reed_frost_multi.2(I0, N NHosts, q0, q1)
---------------------------------------------------
# runs Reed Frost model across multiple hosts (chain multinomial)
# same as above but uses a C routine for main loop

---------------------------------------------------
reed_frost_multi_demo(N, NHosts)
----------------------------------------------------
# runs the demo, calling reed_frost_multi.2

----------------------------------------------------
plot_reed_frost(rf,psfile)
-----------------------------------------------------
# plotting function
