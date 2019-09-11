


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


def croping_after_shift(image, shift_x, shift_y):
    if shift_x < 0:
        image = image[:, int(np.ceil(-shift_x)):]
    else:
        image = image[:, :-int(np.ceil(shift_x))]
    if shift_y < 0:
        image = image[int(np.ceil(-shift_y)):, :]
    else:
        image = image[:-int(np.ceil(shift_y)), :]
    return np.array(image, dtype=np.float)




def createFolder(directory):
    '''
    function to create directories, if they dont already exist
    '''
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)




folder="/home/user/Desktop/20190216/SNU"

def cut_images(folder,selector=""):
    files_dict = {}
    for subdir, dirs, files in os.walk(folder):
        bf_check = sum(["before" in x or "Before" in x for x in dirs] + ["after" in x or "After" in x for x in dirs])
        if not selector in subdir:
            continue
        if bf_check == 2 and not subdir in files_dict:  # also cheks if key already exists
            experiment = os.path.split(subdir)[1]
            files_dict[experiment] = []
            print(experiment)
        if "before" not in subdir and "after" not in subdir and "Before" not in subdir and "After" not in subdir :  #
            continue


        for file in files:
            if re.match("\d{0,3}_{0,1}fluo.*",file,flags=re.I):  # finding flou images of beads
                files_dict[experiment].append(os.path.join(subdir,file))
            if re.match("\d{0,3}_{0,1}BF_before.*", file,flags=re.I):
                files_dict[experiment].append(os.path.join(subdir, file))






    for keys,values in files_dict.items():
        new_folder =os.path.join(folder, keys + "shift")
        createFolder(new_folder)

        numbers=[re.match("(\d{0,3})_{0,1}[fluo|BF_before].*", os.path.split(file)[1]).group(1) for file in values]
        for number in np.unique(numbers):

            print(keys,number)
            file_before=[image for image in values if "before" in image and re.match(".*"+number+"_{0,1}fluo.*", image,flags=re.I)][0]


            file_after =[image for image in values if "after" in image and re.match(".*"+number+"_{0,1}fluo.*", image,flags=re.I)][0]


            file_bf_before=[image for image in values if "before" in image and re.match(".*"+number+"_{0,1}BF_before.*", image,flags=re.I)][0]

            img_b = plt.imread(file_before)
            img_a = plt.imread(file_after)
            imb_b_BF=plt.imread(file_bf_before)

            shift_values = register_translation(img_b, img_a, upsample_factor=100)
            shift_y = shift_values[0][0]
            shift_x = shift_values[0][1]

            # using interpolation to shift subpixel precision
            img_shift_b = shift(img_b, shift=(-shift_y, -shift_x), order=5)
            img_shift_bf = shift(imb_b_BF, shift=(-shift_y, -shift_x), order=5)





            b = normalizing(croping_after_shift(img_shift_b, shift_x, shift_y))
            a = normalizing(croping_after_shift(img_a, shift_x, shift_y))
            bf = normalizing(croping_after_shift(img_shift_bf, shift_x, shift_y))

            b_save = Image.fromarray(b * 255)
            a_save = Image.fromarray(a * 255)
            bf_save = Image.fromarray(bf * 255)


            b_save.save(os.path.join(new_folder, number+"after_shift.tif"))
            a_save.save(os.path.join(new_folder,  number+"before_shift.tif"))
            bf_save.save(os.path.join(new_folder,  number+"bf_before_shift.tif"))




if __name__ == '__main__':
    folder="/home/user/Desktop/20190216/SNU"
    cut_images(folder)









'''
folder="/media/user/GINA1-BK/traktion_force_microscopy/7_8_PERCENT/"
file_before="/media/user/MAGDALENA/TFM/20190216/Images for analysis/SNU/7_85_percent/before/01fluo.tif"
file_after="/media/user/MAGDALENA/TFM/20190216/Images for analysis/SNU/7_85_percent/after/01fluo.tif"
file_before_bf="/media/user/MAGDALENA/TFM/20190216/Images for analysis/SNU/7_85_percent/before_BF.tif"
img_b = plt.imread(file_before)
img_a = plt.imread(file_after)
img_bf= plt.imread(file_before_bf)
shift_values = register_translation(img_b, img_a, upsample_factor=100)

shift_y = shift_values[0][0]
shift_x = shift_values[0][1]
img_shift_b = shift(img_b, shift=(-shift_y, -shift_x), order=5)
img_shift_bf = shift(img_bf, shift=(-shift_y, -shift_x), order=5)

b=normalizing(croping_after_shift(img_shift_b,shift_x,shift_y))
a=normalizing(croping_after_shift(img_a,shift_x,shift_y))
bf=normalizing(croping_after_shift(img_shift_bf,shift_x,shift_y))


'''





