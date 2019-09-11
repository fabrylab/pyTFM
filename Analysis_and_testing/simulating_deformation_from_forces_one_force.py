import numpy as np
import matplotlib.pyplot as plt

def get_xy_for_quiver(u):
    xs=np.zeros(np.shape(u))
    for i in range(np.shape(u)[0]):
        xs[i,:]=np.arange(0,np.shape(u)[1],1)
    ys = np.zeros(np.shape(u))
    for j in range(np.shape(u)[1]):  ## is inverted in other skript...
        ys[:, j] = np.arange(0,np.shape(u)[0], 1)
    return xs, ys



# single force

sigma=0.49#poison ratio
young=25536# youngsmodulus

def_abs= np.zeros((100,100))
def_x=np.zeros((100,100))
def_y=np.zeros((100,100))
F_pos=(50.5,50.5)
F=np.array([10000,10000]) ##? unit




A=(1+sigma)/(np.pi*young)

dist_test=np.zeros(np.shape(def_abs))
for i in range(np.shape(def_abs)[0]):
    for j in range(np.shape(def_abs)[1]):

        x=i-F_pos[0]
        y=j-F_pos[1]
        r=np.sqrt(x**2+y**2)
        K=(A/(r**3))*np.array([[(1-sigma)*(r**2)+sigma*(x**2), sigma*x*y],
                   [sigma*x*y, (1-sigma)*(r**2)+sigma*(y**2)]])

        defo=np.matmul(K,np.transpose(F))

        def_abs[i,j] = np.sqrt(defo[0]**2+defo[1]**2)
        def_x[i, j] =defo[0]
        def_y[i, j] =defo[1]


def_x_pos,def_y_pos=get_xy_for_quiver(def_x)# positions for quiver


plt.close("all")
plt.figure();plt.imshow(def_abs)
plt.colorbar()

scale=1
plt.arrow(F_pos[0],F_pos[1],F[0]*scale,F[1]*scale,head_width=3,facecolor="red",edgecolor="red")
plt.quiver(def_x_pos,def_y_pos,def_x,-def_y)



