import h5py
import numpy as np
from matplotlib import pyplot as plt
import scipy.ndimage
import os
import pandas as pd
import imageio
from skimage import io
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage import morphology
from skimage import restoration
import skimage.filters as skifi
import tifffile as tif
from scipy.ndimage import gaussian_filter
import matplotlib
import yaml


def read_parameters():
    parameters = dict()

    with open("../base/parameters.yml") as file:
        parameters = yaml.load(file, Loader=yaml.FullLoader)
    with open("../local/parameters.yml") as file:
        parameters_local = yaml.load(file, Loader=yaml.FullLoader)

    # overwrite global parameters with local setting
    for key in parameters_local:
        parameters[key] = parameters_local[key]

    return parameters


def read_key_file(parameters):
    file_path = parameters['data_folder'] + "00_Extra/" + parameters["key_file_name"]

    key_file = pd.read_excel(file_path)

    return key_file


def background_substraction(parameters, image):

    shape_z = 1
    shape_x = 5
    shape_y = 5
    intensity = 255
    radius = 10

    print("Starting background substraction ...")

    if parameters["substract_background"] == "ellipsoid_kernel": 
        
        print("Apply background substraction with ellipsoid kernel ((%s, %s, %s), %s)" % (shape_z, shape_x, shape_y, intensity) )

        background = restoration.rolling_ball(
            image,
            kernel=restoration.ellipsoid_kernel((shape_z, shape_x, shape_y), intensity)
        )
    else:
        print("Apply background substraction using rolling_ball method with radius: %s" % radius )
        background = restoration.rolling_ball(image, radius = radius)
        

    return background 


def thresholding_3D(parameters, image_blurred):
    '''
        input:
            parameters  - dict() with all meta parameters for analysis
            image - numpy array containing image data 2D or 3D
        output:
            image_labeled - numpy array 2D/3D with labeled macrophages
            thresh - threshold that was applied to obtain the mask

        TODO: - split into thresholding and labelling
              - move to common if used in 2D nd 3D
    '''

    # use Gaussian filter to smooth images
    # print(" Gaussian filtering...")
    # image_blurred = gaussian_filter(image, sigma=parameters['sigma'])

    # thresholding images
    thresh = 10
    if parameters["thresholding_method"] == "isodata":
    	thresh = skifi.threshold_isodata(image_blurred)
    if parameters["thresholding_method"] == "li":
    	thresh = skifi.threshold_li(image_blurred)
    if parameters["thresholding_method"] == "mean":
    	thresh = skifi.threshold_mean(image_blurred)
    if parameters["thresholding_method"] == "minimum":
    	thresh = skifi.threshold_minimum(image_blurred)
    if parameters["thresholding_method"] == "niblack":
    	thresh = skifi.threshold_niblack(image_blurred)
    if parameters["thresholding_method"] == "otsu":
    	thresh = skifi.threshold_otsu(image_blurred)
    	print("Using thresholding method otsu with value %s" % thresh)
    if parameters["thresholding_method"] == "sauvola":
    	thresh = skifi.threshold_sauvola(image_blurred, window_size = parameters["window_size_thresholding"])
    	print("Using thresholding method sauvola")
    if parameters["thresholding_method"] == "yen":
        thresh = skifi.threshold_yen(image_blurred)
    if parameters["thresholding_method"] == "triangle":
    	thresh = skifi.threshold_yen(image_blurred)
    if parameters["thresholding_method"] == "local":
    	thresh = skifi.threshold_local(image_blurred)

    #image_mask = np.where(image_blurred > thresh, True, False)
    image_mask = image_blurred > thresh
    del image_blurred

    # remove artefacts connected to border
    #print(" Removing small objects...")
    #image_mask = clear_border(image_mask)
    #image_mask = morphology.remove_small_objects(image_mask, parameters["macrophages_small_objects"]) #, connectivity=2)
       
    return image_mask, thresh

def get_macrophage_pairs(parameters, macrophage_pos_algorithm, macrophage_pos_annotated):

    macrophage_pairs_ = pd.DataFrame()
    macrophage_pairs = pd.DataFrame()
   
    counter = 1

    for index_1, row_1 in macrophage_pos_algorithm.iterrows():
        min_distance2 = np.Infinity
        label_annotated = 0
        label_algorithm = 0
        for index_2, row_2 in macrophage_pos_annotated.iterrows():
            distance2 = (row_1["x_centroid"] - row_2["Y"])**2
            distance2 += (row_1["y_centroid"] - row_2["X"])**2
                    
            if distance2 < min_distance2:
                min_distance2 = distance2
                label_algorithm = row_1["macrophage_label"]
                label_annotated = row_2["number"]
                x_ann = row_2["X"]
                y_ann = row_2["Y"]
         
        macrophage_pairs_.at[counter, "label_algorithm"] = label_algorithm
        macrophage_pairs_.at[counter, "label_annotated"] = label_annotated
        macrophage_pairs_.at[counter, "distance"] = np.sqrt(min_distance2) 
        macrophage_pairs_.at[counter, "x_alg"] = row_1["y_centroid"] 
        macrophage_pairs_.at[counter, "x_ann"] = x_ann 
        macrophage_pairs_.at[counter, "y_alg"] = row_1["x_centroid"] 
        macrophage_pairs_.at[counter, "y_ann"] = y_ann 
        macrophage_pairs_.at[counter, "z_alg"] = row_1["z_centroid"] 
    
        counter += 1
    
    counter = 1
    
    for index_1, row_1 in macrophage_pos_annotated.iterrows():
        min_distance2 = np.Infinity
        label_annotated = 0
        label_algorithm = 0
        for index_2, row_2 in macrophage_pairs_.iterrows():
            distance2 = (row_2["x_alg"] - row_1["X"])**2
            distance2 += (row_2["y_alg"] - row_1["Y"])**2
                    
            if distance2 < min_distance2:
                min_distance2 = distance2
                label_annotated = row_1["number"]
                label_algorithm = row_2["label_algorithm"]
                
        
        check_df =  macrophage_pairs_[macrophage_pairs_["label_algorithm"] == label_algorithm]
        
        #print(check_df)

        if check_df["label_annotated"].iloc[0] == label_annotated:
            
            #if check_df["distance"]
            macrophage_pairs.at[counter, "label_algorithm"] = label_algorithm
            macrophage_pairs.at[counter, "label_annotated"] = label_annotated
            macrophage_pairs.at[counter, "distance"] = np.sqrt(min_distance2) 

            counter += 1
    
    return macrophage_pairs, macrophage_pairs_
