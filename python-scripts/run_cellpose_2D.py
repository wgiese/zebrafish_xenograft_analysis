import pandas as pd
import yaml
import numpy as np
import scipy.ndimage
import time, os, sys
import skimage.io
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
from cellpose import utils
from matplotlib import cm
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops, regionprops_table
import math
from scipy.ndimage import gaussian_filter
import sys
from cellpose import models, io, plot
import functions_common

### read parameters

parameters = functions_common.read_parameters()
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())


data_path = parameters["data_folder"]
folder_2d_data = "/03_Preprocessed_Data/01_2D/"
use_gpu = parameters["use_gpu"]
output_folder = data_path + "/cellpose_segmentation/"
experiments = "annotated"

for index, row in key_file.iterrows():

    short_name = str(row["short_name"]) 
    file_path = data_path + folder_2d_data + short_name + ".tif"
    
    if not os.path.exists(file_path):
        continue

    if (experiments == "annotated") and (not row["macrophages_annotated"]):
        print("No annotation for %s" % short_name)
        continue

    print("open image file %s" % file_path)
        
    img = skimage.io.imread(file_path)
    
    print(img.shape)
   
    #for time in range(img.shape[0]):
    #    print(time)
    #    macrophage_img = img[time,:,:,parameters["channel_macrophages"]]

    time = 0
    macrophage_img = img[time,:,:,parameters["channel_macrophages"]]
    
    channels = [0,0]
    model = models.Cellpose(gpu=use_gpu, model_type='cyto')
    masks, flows, styles, diams = model.eval(macrophage_img, diameter=25, channels=channels)

    fig, ax = plt.subplots(figsize=(15,15))
    ax.imshow(macrophage_img, cm.binary)
    ax.imshow(masks, cm.Set3, alpha = 0.5)

    annotations_file = parameters["data_folder"] + "04_Processed_Data/01_Annotated_Macrophages/" + short_name + '.csv'
    if os.path.exists(annotations_file):
        annotated_positions = pd.read_csv(annotations_file, sep = ";")
        print("annotated postions file exists ...")
        print(annotated_positions.head())
        annotated_df = annotated_positions[annotated_positions["time_point"] == time]
        ax.plot(annotated_df['X'], annotated_df['Y'], 'rx', markersize = 15)
    
    plt.savefig(output_folder + short_name + "-%s.png" % time)
    
           
    
#io.masks_flows_to_seg(img , masks, flows, diams, output_filepath , channels)

# load result
#cellpose_seg = np.load(output_filepath + "_seg.npy", allow_pickle=True)

#mask = cellpose_seg.item()['masks']
#for key in cellpose_seg.item():
#    print(key)
#
#print(mask.shape)

