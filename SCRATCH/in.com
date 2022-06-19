dt      nstep   noutt     noutf   
0.005   23       1        10         
xl      yl       zl 
1.0     8.0     32.0 
t0      gridalfa   gridbeta   gridgamma
-55.0d0   35.0d0     1.60d0     60.0d0
de    rhos   beta_e   tau    eta     mm        pso     psoeq  phieq  eq_l   asym (n/m)
0.0    0.0      2.       0.0  0.0003     1      1.0d-8     0.8    1.0    1.      0.
istart  omegax   omega   alfa1
 0       0.45    0.66     2.0
kbc_p (Poisson bc)    kbc_h (Helmoltz bc)   0=DR, 1=NS   
 1                       1
front_l    front_r    delta_ps               
 -8.         -0.4      0.1
irestart_zero
 0
