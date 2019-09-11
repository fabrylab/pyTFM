import numpy as np
import matplotlib.pyplot as plt

f1=np.zeros((1000))
f1[550:600]=1
f1_ind=np.linspace(-500,500,1000)




def ft(ft_f1_ind):
    return np.sum(f1*np.cos(ft_f1_ind*f1_ind))

ft_f1=[]
for i in np.linspace(-2*np.pi,2*np.pi,1000):
    ft_f1.append(ft(i))
ft_f1=np.array(ft_f1)
plt.figure()
plt.plot(f1_ind,f1)
plt.figure()
plt.plot(np.linspace(-2*np.pi,2*np.pi,1000),ft_f1)

