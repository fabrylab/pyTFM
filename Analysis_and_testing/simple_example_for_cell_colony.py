import random
import numpy as np
import sympy as symp
import itertools


##
#alternative parametrization
l=[[1,1,0,0,0],[-1,0,-1,-1,0],[0,-1,+1,0,1],[0,0,0,1,-1],[0,-1,1,1,0],[-1,0,-1,0,-1]]


ar1,ar2,ar3,ar4=[1,-1,-1,1]
A=np.array(list(itertools.combinations(l,6))[0])
F=np.array([[ar1],[ar2],[ar3],[ar4],[ar3+ar4],[ar2+ar4]])
F=F.flatten().astype(int)
#F=A*x


np.linalg.lstsq(A,F)

#x1=np.array([1,0,0,0,1])
x2=np.linalg.lstsq(A[:5],F[:5])[0]
#np.matmul(A,x1.T)
np.matmul(A,x2.T)



### gibt so keine eindeutige lösung!!!!!!!