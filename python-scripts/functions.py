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
from skimage.filters import threshold_otsu
from skimage.external import tifffile as tif
from scipy.ndimage import gaussian_filter
import matplotlib
import yaml


font = {'family': 'normal', 'weight' : 'bold', 'size' : 22}
matplotlib.rc('font', **font)

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
    file_path = parameters['data_folder'] + "00_Extra/Details_imaging_data.xls"

    key_file = pd.read_excel(file_path)

    return key_file

def otsu_thresholding_3D(parameters, image):

    image_labeled = np.zeros(image.shape, dtype=bool)

    for tp in range(image.shape[0]):
        # use Gaussian filter to smooth images
        image_labeled[tp] = gaussian_filter(image[tp], sigma=parameters['sigma'])

        # Otsu thresholding images
        image_labeled[tp] = np.where(image_labeled[tp] > threshold_otsu(image_labeled[tp]), True, False)

        # remove artefacts connected to border
        image_labeled[tp] = clear_border(image_labeled[tp])
        image_labeled[tp] = morphology.remove_small_objects(image_labeled[tp], 2500, connectivity=2)

        # label image regions
        image_labeled[tp] = label(image_labeled[tp])

    return image_labeled

def loop_over_key_file(parameters, key_file):

    for file in key_file["New name"].unique():
        file_path = parameters["data_folder"] + "03_Processed_Data/3D/" + file
        if os.path.exists(file_path):
            movie = np.array(io.imread(file_path))
            im_tumor = movie[:, :, :, :, 0]
            im_macrophages = movie[:, :, :, :, 1]
            im_vessels = movie[:, :, :, :, 2]

            # Otsu thresholding
            mask_tumor = otsu_thresholding_3D(parameters, im_tumor)
            mask_macrophages = otsu_thresholding_3D(parameters, im_macrophages)
            mask_vessels = otsu_thresholding_3D(parameters, im_vessels)

            # save 2D mask projections to file
            mask_2D_tumor = np.zeros([:im_tumor.shape[0], :im_tumor.shape[2], :im_tumor.shape[3]])
            mask_2D_macrophages =
            mask_2D_vessels =

        else:
            print("The file '%s' does not exist.", %file_path)

    return None

#def calculate_tumor_macrophage_distance_3D(parameters, key_file):
