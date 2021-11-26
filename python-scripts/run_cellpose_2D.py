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
#experiments = "all"

for index, row in key_file.iterrows():

    short_name = str(row["short_name"]) 
    file_path = data_path + folder_2d_data + short_name + ".tif"
    
    if not os.path.exists(file_path):
        continue

    if (experiments == "annotated") and (not row["macrophages_annotated"]):
        print("No annotation available for %s" % short_name)
        continue

    print("open image file %s" % file_path)
        
    img = skimage.io.imread(file_path)
    
    print(img.shape)

    coordinates_2D = pd.DataFrame()
    index = 0  

    for time in range(img.shape[0]):
        print(time)

        macrophage_img = img[time,:,:,parameters["channel_macrophages"]]
    
        channels = [0,0]
        model = models.Cellpose(gpu=use_gpu, model_type='cyto')
        masks, flows, styles, diams = model.eval(macrophage_img, diameter=25, channels=channels)
        masks = skimage.segmentation.clear_border(masks)


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
        else:
            print("Annotation can not be loaded, file does not exist.")
        
        plt.savefig(output_folder + short_name + "-%s-cellpose.png" % time)

       
        if os.path.exists(annotations_file):
            fig, ax = plt.subplots(figsize=(15,15))
            ax.imshow(macrophage_img, cm.binary)
            plt.savefig(output_folder + short_name + "-%s-annotations.png" % time)
            annotated_positions = pd.read_csv(annotations_file, sep = ";")
            print("annotated postions file exists ...")
            print(annotated_positions.head())
            annotated_df = annotated_positions[annotated_positions["time_point"] == time]
            ax.plot(annotated_df['X'], annotated_df['Y'], 'rx', markersize = 15)
            plt.savefig(output_folder + short_name + "-%s-annotations.png" % time)
        else:
            print("Annotation can not be loaded, file does not exist.")
        

        for label in range(1,np.max(masks)-1):
            single_cell_mask = np.where(masks == label, 1, 0)
            regions = skimage.measure.regionprops(single_cell_mask, intensity_image = macrophage_img)
            
            for props in regions:
                x_cell, y_cell = props.centroid
                #minor_axis_length = props.minor_axis_length
                #major_axis_length = props.major_axis_length
                #eccentricity = props.eccentricity
                mean_intensity = props.mean_intensity
                max_intensity = props.max_intensity
                min_intensity = props.min_intensity
                area = props.area
                perimeter = props.perimeter  
            
            coordinates_2D.at[index,"time_point"] = time
            coordinates_2D.at[index,"number"] = label
            coordinates_2D.at[index,"Area"] = area
            coordinates_2D.at[index,"Mean"] = mean_intensity
            coordinates_2D.at[index,"Min"] = min_intensity
            coordinates_2D.at[index,"Max"] = max_intensity
            coordinates_2D.at[index,"X"] = x_cell
            coordinates_2D.at[index,"Y"] = y_cell

            index +=1

        coordinates_2D.to_csv(output_folder + short_name + ".csv", sep=";")

