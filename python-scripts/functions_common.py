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

 
    print("Type of threshold")
    print(type(thresh))
    if (type(thresh) == int or np.int64) or (type(thresh) == float):
        thresh_ = thresh
    else:
        thresh_ = "local"

    #image_mask = np.where(image_blurred > thresh, True, False)
    image_mask = image_blurred > thresh
    del image_blurred

    # remove artefacts connected to border
    # print(" Removing small objects...")
    image_mask = clear_border(image_mask)
    image_mask = morphology.remove_small_objects(image_mask, parameters["macrophages_small_objects"], connectivity=2)
       
    return image_mask, thresh_
