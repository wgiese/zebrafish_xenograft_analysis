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
sys.path.insert(0,"../")
import functions_common
import datetime

### read parameters

parameters = functions_common.read_parameters(base_folder = "../base/", local_folder = "../local/")
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())

### filter experiments

#experiments = "annotated"
experiments = "all"
injection_time_dpi = parameters["dpi"]
filter_experiments = parameters["filter_experiments"]
skip_experiments = parameters["skip_experiments"]

print("injection time:")
print(injection_time_dpi)

data_path = parameters["data_folder"]
folder_3d_data = "/03_Preprocessed_Data/02_3D/"
use_gpu = parameters["use_gpu"]
output_folder = parameters["cp_3D_output_path"]

# filter key file
if (not isinstance(injection_time_dpi, list)):
  key_file = key_file[key_file["dpi"]==injection_time_dpi]
key_file = key_file[~key_file["short_name"].isin(skip_experiments)]
key_file = key_file.dropna(subset=["short_name"])


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

    analysis_summary.to_csv(output_folder + "analysis_summary.csv", index = False)

    short_name = str(row["short_name"])
    file_path = data_path + folder_3d_data + short_name + ".tif"
    fish_id = short_name.split("_")[0] + "_" + short_name.split("_")[1] + "_" + short_name.split("_")[3]
    
    if not os.path.exists(file_path):
        print(file_path)
        print("file path does not exist!")
        index_summary += 1
        continue

    print("open image file:", file_path)

    img = skimage.io.imread(file_path)

    print("image dimensions: ")
    print(img.shape)

    coordinates_3D = pd.DataFrame()
    index = 0

    for time in range(img.shape[0]):
        print(time)

        # TODO: check channel order needed for cellpose, in the documentation it says "segment list of images x, or 4D array - Z x nchan x Y x X"

        macrophage_img = img[time,:,:,:,parameters["channel_macrophages"]]

        channels = [0,0]

        if parameters["cp_model_path"] == "None":
            model = models.Cellpose(gpu=use_gpu, model_type='cyto')
        else:
            model = models.CellposeModel(gpu=use_gpu, pretrained_model = parameters["cp_model_path"]) 

        if parameters["diameter"] == "None":
            masks, flows, styles = model.eval(macrophage_img, diameter=None, channels=channels)
        else:
            masks, flows, styles = model.eval(macrophage_img, diameter=parameters["diameter"], channels=channels, anisotropy = parameters["anisotropy"], z_axis = parameters["z_axis"], flow_threshold = parameters["flow_threshold"], stitch_threshold = parameters["stitch_threshold"])

### old code

#img = skimage.io.imread(filename)
#output_path = parameters['output_folder']
#output_filename = parameters["output_filename"]
#output_filepath = output_path + output_filename

# extract channels

#print(img.shape)

#img = img.reshape((img.shape[1],img.shape[2],img.shape[0]))

#print(img.shape)

#use_gpu = parameters["use_gpu"]
#model = models.Cellpose(gpu=use_gpu, model_type='cyto')

# choice for greyscale
#channels = [0,0]

#masks, flows, styles, diams = model.eval(img, diameter=40, channels=channels, anisotropy = 2.5, do_3D = True)

#io.masks_flows_to_seg(img , masks, flows, diams, output_filepath , channels)

# load result
#cellpose_seg = np.load(output_filepath + "_seg.npy", allow_pickle=True)

#mask = cellpose_seg.item()['masks']
#for key in cellpose_seg.item():
#    print(key)

#print(mask.shape)

'''
im_junction = img[:,:,parameters["channel_junction"]]
im_nucleus = img[:,:,parameters["channel_nucleus"]]
im_golgi = img[:,:,parameters["channel_golgi"]]

# use junction and nucleus channel for cellpose segmentation 
im_seg = np.array([im_junction, im_nucleus])


model = models.Cellpose(gpu=True, model_type='cyto')

channels = [1,2]

masks, flows, styles, diams = model.eval(im_seg, diameter=100, channels=channels)

io.masks_flows_to_seg(im_seg , masks, flows, diams, output_filename, channels)

#read cellpose segmentation or feature extraction

cellpose_seg = np.load(output_filename + "_seg.npy", allow_pickle=True)
mask = cellpose_seg.item()['masks']
for key in cellpose_seg.item():
    print(key)

img_golgi_blur = gaussian_filter(img[:,:,0], sigma= 3)
img_nuclei_blur = gaussian_filter(img[:,:,1], sigma=3)

nuclei_mask = np.where(img_nuclei_blur > threshold_otsu(img_nuclei_blur), True, False)
golgi_mask = np.where(img_golgi_blur > threshold_otsu(img_golgi_blur), True, False)

nuclei_label = nuclei_mask*mask
golgi_label = golgi_mask*mask
print(np.max(mask))

# feature extraction

single_cell_props = pd.DataFrame()
counter = 0

for label in range(1,np.max(mask)-1):
    
    single_cell_mask = np.where(mask ==label, 1, 0)
    single_nucleus_mask = np.where(nuclei_label==label, 1, 0)
    single_golgi_mask = np.where(golgi_label ==label, 1, 0)    
    
    area_cell_px2 = len(single_cell_mask[single_cell_mask==1])
    area_golgi_px2 = len(single_golgi_mask[single_golgi_mask==1])
    area_nucleus_px2 = len(single_nucleus_mask[single_nucleus_mask==1])    
    
    if (area_nucleus_px2) < 10 or (area_golgi_px2 < 10):
        continue
                
    #print(len(single_cell_mask[single_cell_mask==1]))
    #print(len(single_nucleus_mask[single_nucleus_mask==1]))
    #print(len(single_golgi_mask[single_golgi_mask==1]))
      
    regions = regionprops(single_cell_mask)
    for props in regions:
        x_cell, y_cell = props.centroid
        orientation = props.orientation
        minor_axis_length = props.minor_axis_length
        major_axis_length = props.major_axis_length
                
    single_cell_props.at[counter, "label"] = label
    single_cell_props.at[counter, "X_cell"] = x_cell
    single_cell_props.at[counter, "Y_cell"] = y_cell
    single_cell_props.at[counter, "shape_orientation"] = x_cell
    single_cell_props.at[counter, "major_axis_length"] = major_axis_length
    single_cell_props.at[counter, "minor_axis_length"] = minor_axis_length
    
    regions = regionprops(single_nucleus_mask)
    for props in regions:
        x_nucleus, y_nucleus = props.centroid
    
    single_cell_props.at[counter, "X_nuc"] = x_nucleus
    single_cell_props.at[counter, "Y_nuc"] = y_nucleus
    
    regions = regionprops(single_golgi_mask)
    for props in regions:
        x_golgi, y_golgi = props.centroid
    
    single_cell_props.at[counter, "X_golgi"] = x_golgi
    single_cell_props.at[counter, "Y_golgi"] = y_golgi
    
    distance2 = (x_golgi - x_nucleus)**2
    distance2 += (y_golgi - y_nucleus)**2
    
    single_cell_props.at[counter, "distance"] = np.sqrt(distance2)
    
    vec_x = x_golgi - x_nucleus
    vec_y = y_golgi - y_nucleus
    angle_rad_ = np.arctan2(vec_x, vec_y)
    
    angle_rad = angle_rad_
    
    if (angle_rad_ < 0.0):
        angle_rad = 2.0*np.pi + angle_rad_
    
    single_cell_props.at[counter, "vec_X"] = vec_x
    single_cell_props.at[counter, "vec_Y"] = vec_y
    single_cell_props.at[counter, "angle_rad"] = angle_rad
    single_cell_props.at[counter, "angle_deg"] = 180.0*angle_rad/np.pi   
    
    counter += 1
    
single_cell_props.to_csv(output_filepath + "_.csv") 

# visualiztion of nuclei-golgi vectors

fig, ax = plt.subplots()

nuclei_mask_ = np.where(nuclei_mask == True, 70, 0)
golgi_mask_ = np.where(golgi_mask == True, 20, 0)

ax.imshow(im_junction, cmap=plt.cm.gray, alpha = 1.0)
ax.imshow(np.ma.masked_where(nuclei_mask_ == 0, nuclei_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha=1.0)

ax.imshow(np.ma.masked_where(golgi_mask_ == 0, golgi_mask_),  plt.cm.gist_rainbow, vmin=0, vmax=100, alpha=1.0)

for index, row in single_cell_props.iterrows():
    ax.plot( row["Y_nuc"], row["X_nuc"], '.g', markersize=1)
    ax.plot( row["Y_golgi"], row["X_golgi"], '.m', markersize=1)
    ax.arrow(row["Y_nuc"], row["X_nuc"], row["Y_golgi"]- row["Y_nuc"],row["X_golgi"]- row["X_nuc"], color = 'red', width = 2)

ax.set_xlim(0,im_junction.shape[0])
ax.set_ylim(0,im_junction.shape[1])
plt.savefig(output_filepath + "_nuclei_golgi_vector.pdf")
plt.savefig(output_filepath + "_nuclei_golgi_vector.png")

# visualiztion of nuclei-golgi vectors

fig, ax = plt.subplots()
ax.imshow(im_junction, cmap=plt.cm.gray)

regions = regionprops(mask)
for props in regions:
    y0, x0 = props.centroid
    orientation = props.orientation
    x1 = x0 + math.cos(orientation) * 0.5 * props.minor_axis_length
    y1 = y0 - math.sin(orientation) * 0.5 * props.minor_axis_length
    x2 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
    y2 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length

    ax.plot((x0, x1), (y0, y1), '--r', linewidth=0.5)
    ax.plot((x0, x2), (y0, y2), '--r', linewidth=0.5)
    ax.plot(x0, y0, '.b', markersize=5)

ax.set_xlim(0,im_junction.shape[0])
ax.set_ylim(0,im_junction.shape[1])
plt.savefig(output_filepath + "_cellshape_orientation.pdf")
plt.savefig(output_filepath + "_cellshape_orientation.png")

'''
