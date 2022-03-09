import time, os, sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import skimage.io
from skimage import morphology
from skimage.filters import threshold_otsu
from scipy.ndimage import gaussian_filter
#sys.path.insert(0,"../")
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from skimage.measure import label
import functions_common

### read parameters

parameters = functions_common.read_parameters(base_folder = "base/", local_folder = "local/")
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())

experiments = "all"
injection_time_dpi = parameters["dpi"]
filter_experiments = parameters["filter_experiments"]
skip_experiments = parameters["skip_experiments"]

### setup paths to data

data_path = parameters["data_folder"]
folder_3d_data = "/03_Preprocessed_Data/02_3D/"
use_gpu = parameters["use_gpu"]
output_folder = parameters["output_folder"]

# filter key file

key_file = key_file[key_file["dpi"]==injection_time_dpi]
key_file = key_file[~key_file["short_name"].isin(skip_experiments)]

# iterate over key file, every rowlinks to one image time series and contains corresponding meta information

for index, row in key_file.iterrows():
    short_name = str(row["short_name"])
    file_path = data_path + folder_3d_data + short_name + ".tif"
    fish_id = short_name.split("_")[0] + "_" + short_name.split("_")[1] + "_" + short_name.split("_")[3]

    if not os.path.exists(file_path):
        print(file_path)
        print("file path does not exist!")
        continue

    print("open image file:", file_path)

    img = skimage.io.imread(file_path)
    print("Image loaded.")
    img_tumor = img[:, :, :, :, 0]
    print(img_tumor.shape)
    # img_macrophages = img[:, :, :, :, 1]
    # img_vessels = img[:, :, :, :, 2]
    del img

    # Otsu thresholding of tumor channel
    #print("Applying Gaussian filter...")
    img_tumor_blurred = np.zeros(img_tumor.shape)
    img_tumor_thresh = np.zeros(img_tumor.shape)
    img_labeled_tumor = np.zeros(img_tumor.shape)

    segmentation_path = output_folder + "/tumor_segmentation/"
    if not os.path.exists(segmentation_path):
        os.makedirs(segmentation_path)

    print("Otsu thresholding...")
    for frame in range(img_tumor.shape[0]):
        print(str(frame))
        img_tumor_blurred[frame] = gaussian_filter(img_tumor[frame], sigma=parameters['sigma_tumor'])
        for slice in range(img_tumor.shape[1]):
            img_tumor_thresh[frame, slice, :, :] = np.where(img_tumor_blurred[frame, slice, :, :] > threshold_otsu(img_tumor_blurred[frame, slice, :, :]), True, False)
        #img_tumor_thresh[frame] = np.where(img_tumor_blurred[frame] > threshold_otsu(img_tumor_blurred[frame]), True, False)
        img_labeled_tumor[frame] = label(img_tumor_thresh[frame].astype(int))

        #print("random pixel value:", img_tumor_thresh[frame, 10, 100, 100])
        #print("random pixel value:", img_tumor_thresh[frame, 10, 100, 100].astype(bool))
        #print("random pixel value:", np.invert(img_tumor_thresh[frame, 10, 100, 100].astype(bool)))

        fig, axes = plt.subplots(1, 4, figsize=(16, 12))
        axes[0].imshow(np.amax(img_tumor[frame], axis=0))
        axes[1].imshow(np.amax(img_tumor_blurred[frame], axis=0))
        axes[2].imshow(np.amax(img_labeled_tumor[frame], axis=0))
        axes[3].imshow(np.amax(np.multiply(img_tumor[frame], np.invert(img_tumor_thresh[frame].astype(bool))), axis=0))
        #axes[3].imshow(np.amax(np.invert(img_tumor_thresh[frame].astype(int)), axis=0))
        plt.savefig(segmentation_path + short_name + "_" + str(frame) + ".png", format="png", bbox_inches="tight", dpi=150)
        plt.close()

    del img_tumor
    del img_tumor_blurred
    del img_tumor_thresh
    del img_labeled_tumor
