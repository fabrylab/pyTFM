


## generating gifs of shift corrected "before" and "after" images.
# the programm will search trough all sub folders in the given folder,
#output is saved in the "folder" as single gifs

# beware: sometimes files on usb drive are saved as hidden. check  with ls commadn

# data structure:
#....

import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import register_translation
from scipy.ndimage.interpolation import shift
from PIL import Image
import os
import re

import imageio



def normalizing(img):
    img = img - np.percentile(img, 1)  # 1 Percentile
    img = img / np.percentile(img, 99.99)  # norm to 99 Percentile
    img[img < 0] = 0.0
    img[img > 1] = 1.0
    return img



def cut_images_gifs(folder,selector=""):

    files_dict = {}
    for subdir, dirs, files in os.walk(folder):
        bf_check = sum(["before" in x or "Before" in x for x in dirs] + ["after" in x or "After" in x for x in dirs])
        if not selector in subdir:
            continue
        if bf_check == 2 and not subdir in files_dict:  # also cheks if key already exists
            experiment = os.path.split(subdir)[1]
            files_dict[experiment] = []
            print(experiment)
        if not bf_check == 2 and "before" not in subdir and "after" not in subdir and "Before" not in subdir and "After" not in subdir :  #
            continue


        for file in files:
            if re.match("\d{0,3}_{0,1}fluo.*",file):
                files_dict[experiment].append(os.path.join(subdir,file))
        print(subdir)





    for keys,values in files_dict.items():
        numbers = [re.match("(\d{0,3})_{0,1}fluo.*", os.path.split(file)[1]).group(1) for file in values]

        for number in np.unique(numbers):
            paire_img=[file for file in values if re.match(number+"_{0,1}fluo.*", os.path.split(file)[1])]
            print(paire_img)
            file_before=paire_img[0]
            file_after =paire_img[1]

            img_b = plt.imread(file_before)
            img_a = plt.imread(file_after)

            shift_values = register_translation(img_b, img_a, upsample_factor=100)

            shift_y = shift_values[0][0]
            shift_x = shift_values[0][1]

            # using interpolation to shift subpixel precision
            img_shift_b = shift(img_b, shift=(-shift_y, -shift_x), order=5)

            if shift_x < 0:
                img_shift_b = img_shift_b[:, int(np.ceil(-shift_x)):]
                img_shift_a = img_a[:, int(np.ceil(-shift_x)):]
            else:
                img_shift_b = img_shift_b[:, :-int(np.ceil(shift_x))]
                img_shift_a = img_a[:, :-int(np.ceil(shift_x))]

            if shift_y < 0:
                img_shift_b = img_shift_b[int(np.ceil(-shift_y)):, :]
                img_shift_a = img_shift_a[int(np.ceil(-shift_y)):, :]
            else:
                img_shift_b = img_shift_b[:-int(np.ceil(shift_y)), :]
                img_shift_a = img_shift_a[:-int(np.ceil(shift_y)), :]

            img_shift_a = np.array(img_shift_a, dtype=np.float)
            img_shift_b = np.array(img_shift_b, dtype=np.float)
            plt.close("all")

            img_shift_a_norm = normalizing(img_shift_a)
            img_shift_b_norm = normalizing(img_shift_b)

            images=[img_shift_a_norm * 255,img_shift_b_norm * 255]
            print(os.path.join(folder, keys)+number+".gif")
            imageio.mimsave(os.path.join(folder, keys)+number+".gif", images)


if __name__ == '__main__':
    folder="/home/user/Desktop/20190216/SNU"
    folder = "/home/user/Desktop/20190216/SNU/"








