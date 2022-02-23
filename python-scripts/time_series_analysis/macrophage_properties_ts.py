import time, os, sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
sys.path.insert(0,"../")
import functions_common

# read parameters

parameters = functions_common.read_parameters(base_folder = "../base/", local_folder = "../local/")
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())

# plotting parameters
plt.rcParams.update({'font.size': 16})

# setup paths to data

macrophages_path_1dpi = parameters["macrophage_input_path_1dpi"]
summary_key_file_1dpi = pd.read_csv(macrophages_path_1dpi + 'analysis_summary.csv')

macrophages_path_5dpi = parameters["macrophage_input_path_5dpi"]
summary_key_file_5dpi = pd.read_csv(macrophages_path_5dpi + 'analysis_summary.csv')

# read macrophage properties from segmentation results 

def get_macrophage_properties_df(summary_key_file, macrophages_path):
    macrophage_props = pd.DataFrame()
    for index, row in summary_key_file.iterrows():
        single_movie_fn = row["short_name"] + ".csv"
        short_name = single_movie_fn.split(".")[0]
        single_movie_props = pd.read_csv(macrophages_path + single_movie_fn, sep=";")
    if macrophage_props.shape[0]>1:
        macrophage_props = pd.concat([macrophage_properties, single_movie_props], ignore_index = True)
    else:
        macrophage_props =  single_movie_props.copy()

    return macrophage_props


def time_frame_to_min(macrophage_properties_,dpi):
    
    dt_min = macrophage_properties_["dt_min"].iloc[0]
    macrophage_properties = macrophage_properties_.copy()
    
    macrophage_properties["time_in_min"] = macrophage_properties["time_point"]*dt_min + dpi*24.0*60.0
    macrophage_properties["time_in_h"] = macrophage_properties["time_in_min"]/60.0
    macrophage_properties["dpi"] = dpi

    return macrophage_properties
    

macrophage_properties_1dpi = get_macrophage_properties_df(summary_key_file_1dpi, macrophages_path_1dpi)
macrophage_properties_1dpi = time_frame_to_min(macrophage_properties_1dpi, dpi=1)
print(macrophage_properties_1dpi)

macrophage_properties_5dpi = get_macrophage_properties_df(summary_key_file_5dpi, macrophages_path_5dpi)
macrophage_properties_5dpi = time_frame_to_min(macrophage_properties_5dpi, dpi=5)
print(macrophage_properties_1dpi) 

