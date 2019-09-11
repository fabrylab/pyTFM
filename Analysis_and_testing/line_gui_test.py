

from __future__ import print_function
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.image import AxesImage
import numpy as np
from numpy.random import rand
import os
import copy

plt.close("all")
folder1="/media/user/GINA1-BK/traktion_force_microscopy/new_data_from_magdalena26_02_19/WTshift/12"
tx=np.load(os.path.join(folder1,"tx.npy"))
ty=np.load(os.path.join(folder1,"ty.npy"))

fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.set_title('click on points, rectangles or text', picker=True)
ax1.set_ylabel('ylabel', picker=True, bbox=dict(facecolor='red'))
im= ax1.imshow(tx)  #
event_out=0
pos1=[]
pos2=[]
line=0


def on_press(event):
    global pos1,pos2, event_out,line
    event_out=event
    # button
    print('you pressed', event.button, event.xdata, event.ydata)
    if event.button==1: # right mouse button
        if len(pos1)==0:
            pos1=(event.xdata, event.ydata)#
        else:
            pos2=(event.xdata, event.ydata)
            line,=ax1.plot(pos1,pos2,"-",color="red",linewidth="3")
        print(pos1,pos2)
    if event.button==3: # left mouse button
        pos1=[]
        pos2=[]
        line.remove()
    fig.canvas.draw_idle()



   # thisline = event.artist
  #  xdata = thisline.get_xdata()
   # ydata = thisline.get_ydata()
   # ind = event.ind
  #  print('onpick1 line:', zip(np.take(xdata, ind), np.take(ydata, ind)))

fig.canvas.mpl_connect('button_press_event', on_press)