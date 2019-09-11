import numpy as np
import pickle
from tqdm import tqdm
import sympy as symp
Ax=np.load("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/Ax.npy")
Ax=Ax.astype(int)
id=[]
for s in tqdm(range(int(Ax.shape[0]/100)-1),total=int(Ax.shape[0]/100)-1):
    ids=symp.Matrix(Ax[s*100:(s+1)*100,:]).rref()
    ids=np.array(ids[1])+s*100
    id.append(ids)
with open("/media/user/GINA1-BK/data_traktion_force_microscopy/temp_data_for_scripts/Ax_rref.npy", "wb") as f:
    pickle.dump(id, f)
