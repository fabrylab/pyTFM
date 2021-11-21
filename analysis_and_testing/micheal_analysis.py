from PIL import Image
import numpy as np
import os
import matplotlib.pyplot as plt
from pyTFM.frame_shift_correction import correct_stage_drift_stack
from scipy.ndimage.interpolation import shift
from pyTFM.TFM_functions import calculate_deformation, ffttc_traction
from pyTFM.plotting import show_quiver
from tqdm import tqdm

# loading masks and images
ims = Image.open("/home/andreas/Desktop/MichealTFM/ForceImage_TimeDependent.tif")
images = []
for i in range(ims.n_frames):
    ims.seek(i)
    images.append(np.array(ims))
images = np.array(images)

ms = Image.open("/home/andreas/Desktop/MichealTFM/Masks_TimeDependent.tif")
masks = []
for i in range(ms.n_frames):
    ms.seek(i)
    masks.append(np.array(ms))

#ref = np.array(Image.open("/home/andreas/Desktop/MichealTFM/Null_SingleCell.tif"))
# correting the stage drift

# calculate drifts with respect to reference image
ref = images[0] # reference doesnt look much like the images in the stack
images, ref = correct_stage_drift_stack(images, ref)


# deformation calculation
deformations = []
for i, im in tqdm(enumerate(images)):
    defx, defy, _, _ = calculate_deformation(ref, im)
    deformations.append((defx, defy))

# tractions
tractions = []
for (defx, defy) in deformations:
    tx, ty = ffttc_traction(defx, defy, pixelsize1=1, pixelsize2=2, young=1, sigma=0.5, filter="gaussian", fs=1)
    tractions.append((tx, ty))


# display of the fields
plt.ioff() # turns off pop up plots
folder = "/home/andreas/Desktop/MichealTFM/corrected"
os.makedirs(folder, exist_ok=True)
resolution = 200
figsize = (ref.shape[1]/resolution, ref.shape[0]/resolution)
cbar_style = "clickpoints"
plot_style = "clickpoints"
for i in range(len(images)):
    Image.fromarray(images[i]).save(os.path.join(folder, "%03d.png" % i))
    if i > 0:
        fig, axis = show_quiver(deformations[i][0], deformations[i][1], cbar_str="deformations",
                                figsize=figsize, resolution=resolution, cbar_style=cbar_style, plot_style=plot_style)
        fig.savefig(os.path.join(folder, "def_%03d.png" % i ), facecolor=fig.get_facecolor(),
                    edgecolor=fig.get_facecolor(), dpi=resolution)
        fig, axis = show_quiver(tractions[i][0], tractions[i][1], cbar_str="tractions",
                                figsize=figsize, resolution=resolution, cbar_style=cbar_style, plot_style=plot_style)
        fig.savefig(os.path.join(folder, "tractions_%03d.png" % i ), facecolor=fig.get_facecolor(),
                    edgecolor=fig.get_facecolor(), dpi=resolution)
plt.ion()

## make a clickpoints database
import clickpoints
db = clickpoints.DataFile(os.path.join(folder, "db.cdb"),"w")

db.setLayer(layer_name="images")
db.setLayer(layer_name="deformations")
db.setLayer(layer_name="tractions")


for im in os.listdir(folder):
    if im.endswith(".png"):
        frame = int(im[-7:-4])
        layer = "images"
        if "def" in im:
            layer = "deformations"
        if "tractions" in im:
            layer = "tractions"
        db.setImage(sort_index=frame, filename=os.path.join(folder, im), layer=layer)
db.db.close()


# analysis



# traction and deformation analysis:

