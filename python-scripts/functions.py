import h5py
import numpy as np
from matplotlib import pyplot as plt
import scipy.ndimage
import os
import pandas as pd
import imageio
from skimage import io
from skimage.filters import threshold_otsu
from scipy.ndimage import gaussian_filter
import matplotlib
import yaml

font = {'family': 'normal', 'weight' : 'bold', 'size' : 22}
matplotlib.rc('font', **font)

def read_parameters():
    parameters = dict()

    with open("base/parameters.yml") as file:
        parameters = yaml.load(file, Loader=yaml.FullLoader)
    with open("local/parameters.yml") as file:
        parameters_local = yaml.load(file, Loader=yaml.FullLoader)

    # overwrite global parameters with local setting
    for key in parameters_local:
        parameters[key] = parameters_local[key]

    return parameters


def read_key_file(parameters):
    file_path = parameters['data_folder'] + "00_Extra/Details_imaging_data.xls"

    key_file = pd.read_excel(file_path)

    return key_file

#def read_key_file(parameters):
 #   key_file = pd.DataFrame()
  #  return key_file



def calculate_tumor_macrophage_distance_3D(parameters, key_file):
