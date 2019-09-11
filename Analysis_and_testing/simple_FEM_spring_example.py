import numpy as np
import matplotlib.pyplot as plt
'''
this is supposed to show nummerical instability in pure neumann oundary conditions finite elements
'''


d=[]
e=[]
s=[] # strains at first point
for epsilon_f in np.linspace(0.999,1.001,1000):
    k1=1
    k2=2
    f1=1
    f2=0.5
    f3=-1.5 *epsilon_f
    F=np.array([f1,f2,f3])
    epsilon_k=0.00001  ## small error  # such an error realisitc
    K=np.array([[k1,-k1,0],[-k1,k1+k2+epsilon_k,-k2],[0,-k2,k2]]) # is exactly singular without that epsilon

    e.append(epsilon_f)
    ds=np.linalg.solve(K,F)
    strain =np.array([ds[0]-ds[1],ds[1]-ds[2]])
    s.append(strain[0])
    d.append(np.abs(np.sum(ds)))# displacement of total system, should be close to zero, if force balance holds

   # K_cond = np.array([[k1, -k1, 0], [-k1, k1 + k2 + epsilon_k, -k2], [0, -k2, k2],[1, 1, 1]])  # K with condition, that sum of all displacements must e zero
   # F_cond=np.array([f1,f2,f3,0])
   # np.linalg.lstsq(K_cond, F_cond)


    ## here just as diffrence between displaccemntes
plt.figure()
plt.plot(e,d)
plt.plot(e,s)


## conclusion: unstable deformation, but stable strain--> also observed when usiong solids py