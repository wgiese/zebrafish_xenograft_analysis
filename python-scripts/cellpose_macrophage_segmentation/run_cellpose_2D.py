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
import scipy.ndimage as ndi
import sys
from cellpose import models, io, plot
sys.path.insert(0,"../")
import functions_common

### read parameters

parameters = functions_common.read_parameters(base_folder = "../base/", local_folder = "../local/")
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())


data_path = parameters["data_folder"]
folder_2d_data = "/03_Preprocessed_Data/01_2D/"
use_gpu = parameters["use_gpu"]
output_folder = parameters["cp_output_path"] #data_path + "/cellpose_segmentation/"
#experiments = "annotated"
experiments = "all"
injection_time_dpi = parameters["dpi"]
filter_experiments = parameters["filter_experiments"]
skip_experiments = parameters["skip_experiments"]

key_file = key_file[key_file["dpi"]==injection_time_dpi]
key_file = key_file[~key_file["short_name"].isin(skip_experiments)]

if filter_experiments == "annotated":
    print("only use annotated (2D coordinates) macrophages")
    key_file = key_file[key_file["macrophages_annotated"]==1.0]
    print(key_file)


analysis_summary = pd.DataFrame()
index_summary = 0

for index, row in key_file.iterrows():

    print(row['short_name'])
    
    #if row["dpi"] != injection_time_dpi:
    #    print("Skip, injection time not as specified in the parameter file")
    #    continue
	
    #for key in key_file.columns:
    #    analysis_summary.at[index_summary, key] = row[key]
    
    analysis_summary.at[index_summary, "short_name"] = row["short_name"]
    analysis_summary.at[index_summary, "dpi"] = row["dpi"]
    analysis_summary.at[index_summary, "properties_saved"] = "no"
	
    analysis_summary.to_csv(output_folder + "analysis_summary.csv", index = False)
    #if index < 11:
    #    continue

    short_name = str(row["short_name"]) 
    file_path = data_path + folder_2d_data + short_name + ".tif"
    fish_id = short_name.split("_")[0] + "_" + short_name.split("_")[1] + "_" + short_name.split("_")[3]


    if not os.path.exists(file_path):
        print(file_path)
        print("file path does not exist!")
        index_summary += 1
        continue

    #if (experiments == "annotated") and (row["macrophages_annotated"] == 0):
    #    print("No annotation available for %s" % short_name)
    #    continue

    print("open image file %s" % file_path)
        
    img = skimage.io.imread(file_path)
    
    print(img.shape)

    coordinates_2D = pd.DataFrame()
    index = 0  

    for time in range(img.shape[0]):
        print(time)

        macrophage_img = img[time,:,:,parameters["channel_macrophages"]]
    
        channels = [0,0]

        if parameters["cp_model_path"] == "None":
            model = models.Cellpose(gpu=use_gpu, model_type='cyto')
        else:
            model = models.CellposeModel(gpu=use_gpu, pretrained_model = parameters["cp_model_path"]) 
        
        
        #masks, flows, styles, diams = model.eval(macrophage_img, diameter=25, channels=channels)
        
        if parameters["diameter"] == "None":
            masks, flows, styles = model.eval(macrophage_img, diameter=None, channels=channels)
        else:
            masks, flows, styles = model.eval(macrophage_img, diameter=parameters["diameter"], channels=channels, flow_threshold = 0.4)

        
        #print(styles)

        masks = skimage.segmentation.clear_border(masks)
        outline_list= np.array(utils.outlines_list(masks))
        #print(outline_list)
        outlines = np.zeros((macrophage_img.shape[0],macrophage_img.shape[1]))
        for mask_id, outline_coords in enumerate(outline_list):
            outlines[tuple(outline_coords.T)] = mask_id + 1

        width = 2
        outlines_d = ndi.morphology.binary_dilation(outlines.astype(bool), iterations = width)
        outlines_ = np.where(outlines_d == True, 30, 0).T

        fig, ax = plt.subplots(figsize=(15,15))
        ax.imshow(macrophage_img, cm.binary)
        ax.imshow(np.ma.masked_where(masks == 0, masks), cm.Set3, alpha = 0.5)
        ax.imshow(np.ma.masked_where(outlines_ == 0, outlines_),  plt.cm.Reds, vmin=0, vmax=100, alpha = 0.5)
        
       
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
        
        fig, ax = plt.subplots(figsize=(15,15))
        ax.imshow(macrophage_img, cm.binary)
        ax.imshow(np.ma.masked_where(masks == 0, masks), cm.Set3, alpha = 0.5)
        ax.imshow(np.ma.masked_where(outlines_ == 0, outlines_),  plt.cm.Reds, vmin=0, vmax=100, alpha = 0.5)
        
        for mask_id in np.unique(masks):
    
            if mask_id == 0:
                continue
                
            #print(mask_id)
            single_cell_mask = np.where(masks ==mask_id, 1, 0)
            regions = skimage.measure.regionprops(single_cell_mask)
            
            
            for props in regions:
                x_cell, y_cell = props.centroid
                # note, the values of orientation from props are in [-pi/2,pi/2] with zero along the y-axis
                orientation = np.pi/2.0 - props.orientation
                minor_axis_length = props.minor_axis_length
                major_axis_length = props.major_axis_length
                eccentricity = props.eccentricity
                area = props.area
                perimeter = props.perimeter
                
            
            
            
            x0 = x_cell
            y0 = y_cell

            x1_major = x0 + math.sin(orientation) * 0.5 * major_axis_length
            y1_major = y0 + math.cos(orientation) * 0.5 * major_axis_length
            x2_major = x0 - math.sin(orientation) * 0.5 * major_axis_length
            y2_major = y0 - math.cos(orientation) * 0.5 * major_axis_length

            x1_minor = x0 + math.cos(orientation) * 0.5 * minor_axis_length
            y1_minor = y0 - math.sin(orientation) * 0.5 * minor_axis_length
            x2_minor = x0 - math.cos(orientation) * 0.5 * minor_axis_length
            y2_minor = y0 + math.sin(orientation) * 0.5 * minor_axis_length

            ax.plot((y1_major, y2_major), (x1_major, x2_major), '--k', linewidth=2.5)
            ax.plot((y1_minor, y2_minor), (x1_minor, x2_minor), '--k', linewidth=2.5)
            
            ax.text( y_cell,x_cell, str(mask_id), color = "red", fontsize=20)


        #for label in range(1,np.max(masks)-1):
            #single_cell_mask = np.where(masks == label, 1, 0)
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

            coordinates_2D.at[index,"short_name"] = short_name
            coordinates_2D.at[index,"fish_id"] = fish_id
            coordinates_2D.at[index,"time_point"] = time
            coordinates_2D.at[index,"number"] = mask_id
            coordinates_2D.at[index,"Area"] = area
            coordinates_2D.at[index,"Mean"] = mean_intensity
            coordinates_2D.at[index,"Min"] = min_intensity
            coordinates_2D.at[index,"Max"] = max_intensity
            coordinates_2D.at[index,"X"] = x_cell
            coordinates_2D.at[index,"Y"] = y_cell
            coordinates_2D.at[index,"dt_min"] = row["dt_min"]
            coordinates_2D.at[index,"time_in_min"] = row["dt_min"]*time
            coordinates_2D.at[index,"minor_axis_length"] = minor_axis_length
            coordinates_2D.at[index,"major_axis_length"] = major_axis_length
            coordinates_2D.at[index,"major_minor_ratio"] = major_axis_length/minor_axis_length
            coordinates_2D.at[index,"perimeter"] = perimeter
            coordinates_2D.at[index,"eccentricity"] = eccentricity
            coordinates_2D.at[index,"cancer_cells"] = row["cancer_cells"]
            
            index +=1

        plt.savefig(output_folder + short_name + "-%s-cell_properties.png" % time)
        coordinates_2D.to_csv(output_folder + short_name + ".csv", sep=";", index = False)

    analysis_summary.at[index_summary, "properties_saved"] = "yes"
    analysis_summary.to_csv(output_folder + "analysis_summary.csv", index = False)
    index_summary += 1
    

