import pandas as pd
# import yaml
import numpy as np
#import scipy.ndimage
import time, os, sys
import skimage.io
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
from cellpose import utils
from matplotlib import cm
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops, regionprops_table
from skimage import morphology
import math
#from scipy.ndimage import gaussian_filter
import scipy.ndimage as ndi
import sys
from cellpose import models, io, plot
sys.path.insert(0,"../")
import functions_common
import datetime
import argparse
import yaml

### read parameters

parser = argparse.ArgumentParser(description='Run cellpose macrophage segmentation.')
parser.add_argument('param', type=str, help='Path to the parameter file.')

args = parser.parse_args()
print("-------")
print("reading parameters from: ", args.param)
print("-------")

parameter_file  = args.param

#parameters = functions_common.read_parameters(base_folder = "../base/", local_folder = "../local/")
parameters = functions_common.read_parameters(base_params = "../base/parameters.yml", local_params = parameter_file)
print("#"*5,"parameters","#"*5)
print(parameters)

key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())

#experiments = "annotated"
experiments = "all"
injection_time_dpi = parameters["dpi"]
filter_experiments = parameters["filter_experiments"]
skip_experiments = parameters["skip_experiments"]

print("injection time:")
print(injection_time_dpi)

### setup paths to data

data_path = parameters["data_folder"]
folder_2d_data = parameters['image_file_path_2D']#"/03_Preprocessed_Data/01_2D/"
use_gpu = parameters["use_gpu"]
output_folder = parameters["cp_output_path"]

with open(output_folder + "/parameters.yml", 'w') as outfile:
    yaml.dump(parameters, outfile)


# filter key file
if (not isinstance(injection_time_dpi, list)): 
  key_file = key_file[key_file["dpi"]==injection_time_dpi]
key_file = key_file[~key_file["short_name"].isin(skip_experiments)]
key_file = key_file.dropna(subset=["short_name"])


if filter_experiments == "annotated":
    print("only use annotated (2D coordinates) macrophages")
    key_file = key_file[key_file["macrophages_annotated"]==1.0]
    print(key_file)

# create summary file that contains meta information and the possibility to resume

analysis_summary_path = output_folder + "analysis_summary.csv"
index_summary = 0

if parameters["resume"]:
    analysis_summary = pd.read_csv(analysis_summary_path)
    analysis_summary = analysis_summary[analysis_summary["properties_saved"]=="yes"]
    index_summary = analysis_summary.shape[0]
    set_types = {"short_name" : "object", "dpi" : "int8", 
                "properties_saved" : "object", "cellpose_diameter": "float32",
                "cellpose_model" : "object"}
    analysis_summary = analysis_summary.astype(set_types)
    print(analysis_summary.dtypes)
else:
    analysis_summary = pd.DataFrame()


# iterate over key file, every row links to one image time series and contains corresponding meta information

for index, row in key_file.iterrows():

    #if np.isnan(row["short_name"]):
    #    continue

    print("short name (of experiment):", row['short_name'])

    if parameters["resume"]:
        if row["short_name"] in list(analysis_summary["short_name"]):
            print("%s  already analysed -> skip")
            continue
    
    now = datetime.datetime.now()
    date_time = now.strftime("%Y-%m-%d %H:%M:%S")
    
    analysis_summary.at[index_summary, "short_name"] = row["short_name"]
    analysis_summary.at[index_summary, "dpi"] = row["dpi"]
    analysis_summary.at[index_summary, "properties_saved"] = "no"
    analysis_summary.at[index_summary, "cellpose_diameter"] = parameters["diameter"]
    analysis_summary.at[index_summary, "cellpose_model"] = str(parameters["cp_model_path"])
    analysis_summary.at[index_summary, "date_time"] = date_time
    analysis_summary.at[index_summary, "t_start"] = row['t_start']
    analysis_summary.at[index_summary, "t_end"] = row['t_end']
	
    analysis_summary.to_csv(output_folder + "analysis_summary.csv", index = False)

    short_name = str(row["short_name"]) 
    file_path = data_path + folder_2d_data + short_name + parameters["file_ext"]
    fish_id = short_name.split("_")[0] + "_" + short_name.split("_")[1] + "_" + short_name.split("_")[3]


    if not os.path.exists(file_path):
        print(file_path)
        print("file path does not exist!")
        index_summary += 1
        continue

    print("open image file:", file_path)
        
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
        
        if parameters["diameter"] == "None":
            masks, flows, styles = model.eval(macrophage_img, diameter=None, channels=channels)
        else:
            masks, flows, styles = model.eval(macrophage_img, diameter=parameters["diameter"], channels=channels, flow_threshold = parameters["flow_threshold"])


        masks = skimage.segmentation.clear_border(masks)
        masks = morphology.remove_small_objects(masks, parameters["min_macrophage_area"], connectivity=2)

        number_labels = np.max(masks)

        print("Found %s macrophages" % number_labels)

        
        fig, ax = plt.subplots(figsize=(15,15))
        ax.imshow(macrophage_img, cm.binary)
        ax.imshow(np.ma.masked_where(masks == 0, masks), cm.Set3, alpha = 0.5)
        
        if (number_labels > 1):
        
            outline_list= np.array(utils.outlines_list(masks))
            outlines = np.zeros((macrophage_img.shape[0],macrophage_img.shape[1]))
            for mask_id, outline_coords in enumerate(outline_list):
                
                if (outline_coords.T.shape[1] < 10):
                    print(np.count_nonzero(masks == mask_id + 1))
                    print(outline_coords.T.shape)
                    print(outline_coords.T.dtype)
                outlines[tuple(outline_coords.T)] = mask_id + 1

            width = 2
            outlines_d = ndi.morphology.binary_dilation(outlines.astype(bool), iterations = width)
            outlines_ = np.where(outlines_d == True, 30, 0).T

            ax.imshow(np.ma.masked_where(outlines_ == 0, outlines_),  plt.cm.Reds, vmin=0, vmax=100, alpha = 0.5)
        
       
        annotations_file = parameters["data_folder"] + "04_Processed_Data/01_Annotated_Macrophages/" + short_name + '.csv'
        if os.path.exists(annotations_file):
            annotated_positions = pd.read_csv(annotations_file, sep = ";")
            print("annotated positions file exists ...")
            print(annotated_positions.head())
            annotated_df = annotated_positions[annotated_positions["time_point"] == time]
            ax.plot(annotated_df['X'], annotated_df['Y'], 'go', markersize = 15)
        else:
            print("Annotation can not be loaded, file does not exist.")
                
        out_filepath = output_folder + short_name + "-%s-cellpose.png" % time
        print("Saving to ", out_filepath)
        plt.savefig(out_filepath)
        plt.close()
       
        #if os.path.exists(annotations_file):
        #    fig, ax = plt.subplots(figsize=(15,15))
        #    ax.imshow(macrophage_img, cm.binary)
        #    plt.savefig(output_folder + short_name + "-%s-annotations.png" % time)
        #    annotated_positions = pd.read_csv(annotations_file, sep = ";")
        #    print("annotated postions file exists ...")
        #    print(annotated_positions.head())
        #    annotated_df = annotated_positions[annotated_positions["time_point"] == time]
        #    ax.plot(annotated_df['X'], annotated_df['Y'], 'rx', markersize = 15)
        #    #plt.savefig(output_folder + short_name + "-%s-annotations.png" % time)
        #    plt.close()
        #else:
        #    print("Annotation can not be loaded, file does not exist.")
        
        if parameters["plot_properties"]:
            fig, ax = plt.subplots(figsize=(15,15))
            ax.imshow(macrophage_img, cm.binary)
            ax.imshow(np.ma.masked_where(masks == 0, masks), cm.Set3, alpha = 0.5)
            ax.imshow(np.ma.masked_where(outlines_ == 0, outlines_),  plt.cm.Reds, vmin=0, vmax=100, alpha = 0.5)
        
            if os.path.exists(annotations_file):
                annotated_positions = pd.read_csv(annotations_file, sep = ";")
                print("annotated positions file exists ...")
                print(annotated_positions.head())
                annotated_df = annotated_positions[annotated_positions["time_point"] == time]
                ax.plot(annotated_df['X'], annotated_df['Y'], 'go', markersize = 15)
            else:
                print("Annotation can not be loaded, file does not exist.")
         
        cell_eccentricity_mask = np.zeros((macrophage_img.shape[0], macrophage_img.shape[1]))
        cell_circularity_mask = np.zeros((macrophage_img.shape[0], macrophage_img.shape[1]))

        for mask_id in np.unique(masks):
    
            if mask_id == 0:
                continue
                
            #print(mask_id)
            single_cell_mask = np.where(masks ==mask_id, 1, 0)
            regions = skimage.measure.regionprops(single_cell_mask)
            
            
            if parameters["plot_properties"]:
                
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
                minor_axis_length = props.minor_axis_length
                major_axis_length = props.major_axis_length
                eccentricity = props.eccentricity
                mean_intensity = props.mean_intensity
                max_intensity = props.max_intensity
                min_intensity = props.min_intensity
                area = props.area
                perimeter = props.perimeter  

            cell_eccentricity_mask += single_cell_mask * eccentricity
            cell_circularity_mask += single_cell_mask * 4.0*np.pi*area/(perimeter*perimeter)
 
            coordinates_2D.at[index,"short_name"] = short_name
            coordinates_2D.at[index,"fish_id"] = fish_id
            coordinates_2D.at[index,"time_point"] = time
            coordinates_2D.at[index,"number"] = mask_id
            coordinates_2D.at[index,"Area"] = area
            coordinates_2D.at[index,"Mean"] = mean_intensity
            coordinates_2D.at[index,"Min"] = min_intensity
            coordinates_2D.at[index,"Max"] = max_intensity
            #coordinates_2D.at[index,"X"] = x_cell
            #coordinates_2D.at[index,"Y"] = y_cell
            coordinates_2D.at[index,"X"] = y_cell
            coordinates_2D.at[index,"Y"] = x_cell
            coordinates_2D.at[index,"area_mum2"] = area*row['PixelSizeX']*row['PixelSizeY']
            coordinates_2D.at[index,"dt_min"] = row["dt_min"]
            coordinates_2D.at[index,"time_in_min"] = row["dt_min"]*time
            coordinates_2D.at[index,"minor_axis_length_px"] = minor_axis_length
            coordinates_2D.at[index,"minor_axis_length_mum"] = minor_axis_length*row['PixelSizeX']
            coordinates_2D.at[index,"major_axis_length_px"] = major_axis_length
            coordinates_2D.at[index,"major_axis_length_mum"] = major_axis_length*row['PixelSizeX']
            coordinates_2D.at[index,"major_minor_ratio"] = major_axis_length/minor_axis_length
            coordinates_2D.at[index,"perimeter_px"] = perimeter
            coordinates_2D.at[index,"perimeter_mum"] = perimeter*row['PixelSizeX']
            coordinates_2D.at[index,"eccentricity"] = eccentricity
            coordinates_2D.at[index,"circularity"] = 4.0*np.pi*area/(perimeter*perimeter)
            coordinates_2D.at[index,"cancer_cells"] = row["cancer_cells"]
            coordinates_2D.at[index, "t_start"] = row['t_start']
            coordinates_2D.at[index, "t_end"] = row['t_end']
	           
            index +=1

        
        if parameters["plot_properties"]:
            plt.savefig(output_folder + short_name + "-%s-cell_properties.png" % time)
        coordinates_2D.to_csv(output_folder + short_name + ".csv", sep=";", index = False)
        plt.close()

        if parameters["plot_eccentricity"]:
            fig, ax = plt.subplots(figsize=(15,15))
            ax.imshow(macrophage_img, cm.binary)
            #print("Max eccentricity: ", np.max(cell_eccentricity_mask))
            cax = ax.imshow(np.ma.masked_where(masks == 0, cell_eccentricity_mask), cmap=plt.cm.bwr, vmin=0.0, vmax=1.0, alpha = 0.5)
            color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)
            color_bar.set_label("eccentricity")
            propereties_ = coordinates_2D[coordinates_2D['time_point'] == time]
            for k, row_prop in propereties_.iterrows():
                ax.text( row_prop['X'], row_prop['Y'], str(np.round(row_prop["eccentricity"],2)), color = "black", fontsize=10)
            plt.savefig(output_folder + short_name + "-%s-cell_eccentricity.png" % time)
            plt.close()
        if parameters["plot_circularity"]:
            fig, ax = plt.subplots(figsize=(15,15))
            ax.imshow(macrophage_img, cm.binary)
            cax = ax.imshow(np.ma.masked_where(masks == 0, cell_circularity_mask), cmap=plt.cm.bwr, vmin=0.0, vmax=1.0, alpha = 0.5)
            color_bar = fig.colorbar(cax, ax=ax, shrink=0.3)
            color_bar.set_label("circularity")
            propereties_ = coordinates_2D[coordinates_2D['time_point'] == time]
            for k, row_prop in propereties_.iterrows():
                ax.text( row_prop['X'], row_prop['Y'], str(np.round(row_prop["circularity"],2)), color = "black", fontsize=10)
            plt.savefig(output_folder + short_name + "-%s-cell_circularity.png" % time)
            plt.close()

        

    analysis_summary.at[index_summary, "properties_saved"] = "yes"
    analysis_summary.to_csv(output_folder + "analysis_summary.csv", index = False)
    index_summary += 1
    

