
import clickpoints
from pyTFM.database_functions import *

# renaming new convention:
import os
import copy
import shutil
rename_dict = {"x-1_y-1":0,"x-1_y0":1,"x0_y-1":2,"x0_y0":3,"x0_y1":4,"x1_y0":5,"x1_y1":6,"x-1_y1":7,"x1_y-1":8}
folder = os.getcwd()
files = os.listdir(folder)
for f in files:
    match = [k for k in rename_dict.keys() if k in f]
    if len(match)==1:
        f_new = copy.deepcopy(f)
        f_new.replace(match[0])
        shutil.move(f,f_new)




def test_n(phi, st):
    n = np.array([np.cos(phi), np.sin(phi)])
    n_f = np.abs(np.dot(np.matmul(st, n), n))
    s_f = np.sqrt(np.linalg.norm(np.matmul(st, n))**2-n_f**2)
    return n_f, s_f


st = np.array([[5000,0],[0,0.5]])
r = np.linspace(0,2*np.pi,1000)
nvs = [test_n(i,st)[0] for i in r]
svs = [test_n(i,st)[1] for i in r]
print(np.mean(nvs))
print(np.mean(svs))

plt.figure()
plt.plot(r,nvs)
plt.figure()
plt.plot(r,svs)



db=clickpoints.DataFile("/home/user/Software/tracktion_force_microscopy"
                        "/tracktion_force_microscopy/test_images_database_setup/database.cdb")

folders = {"folder1_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/",
                "folder2_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/",
                "folder3_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/",
                "folder_out_txt": "/home/user/Software/tracktion_force_microscopy/tracktion_force_microscopy/test_images_database_setup/"}
search_keys = {"after": "\d{1,4}after", "before": "\d{1,4}before",
                    "cells": "\d{1,4}bf_before",
                    "frames": "(\d{1,4})"}
db.getPath(id=1)
db._AddOption(key="key",value="vlaue")