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
import tifffile as tif
from scipy.ndimage import gaussian_filter
import matplotlib
import yaml

font = {'family': 'normal', 'weight': 'bold', 'size': 22}
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
    file_path = parameters['data_folder'] + "00_Extra/" + parameters["key_file_name"]

    key_file = pd.read_excel(file_path)

    return key_file

#######################################################################################################################
#                                                 Analysis functions                                                  #
#######################################################################################################################

def get_dimension_folder(parameters):
    if parameters['dimension'] == '2D':
        return '01_2D/'
    elif parameters['dimension'] == '3D':
        return '02_3D/'
    else:
        print("Dimensions can only be '2D' or '3D'.")


def get_segmentation_file_path(parameters):
    dim_folder = get_dimension_folder(parameters)
    if parameters['segmentation_method'] == 'ilastik':
        return parameters['data_folder'] + '04_Processed_Data/03_Ilastik/' + dim_folder + '03_Probability_Maps/'

    elif parameters['segmentation_method'] == 'Thresholding':
        return parameters['data_folder'] + '04_Processed_Data/02_Threshold/' + dim_folder

    elif parameters['segmentation_method'] == 'Otsu':
        return parameters['data_folder'] + '04_Processed_Data/04_Otsu/' + dim_folder

    else:
        print('Segmentation method ' + parameters['segmentation_method'] + 'is not specified!')


def load_mask(parameters, mask_name):
    file_path = get_segmentation_file_path(parameters)
    if os.path.exists(file_path + mask_name):
        return np.array(io.imread(file_path + mask_name))
    else:
        print(file_path + mask_name + ' does not exist!')


def calculate_mean_distance_of_macrophages(distance_mask, macrophage_mask):
    distance_map = scipy.ndimage.morphology.distance_transform_edt(np.invert(distance_mask)).flatten()
    # print(distance_map.shape)
    return np.mean(distance_map[macrophage_mask.flatten()])


def plot_distance_map(distance_mask, out_path):
    """
    Function for plotting the distance map to check it by eye. Implemented because ilastik distance maps didn't work
    and for Thresholding the plots looked like for just taking the macrophage volume.
    :param distance_mask:
    :param out_path:
    :return:
    """
    distance_map = scipy.ndimage.morphology.distance_transform_edt(np.invert(distance_mask))
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.imshow(distance_map)
    plt.savefig(out_path, format="png", bbox_inches="tight", dpi=100)
    plt.close()


def get_volume(frame):
    return np.count_nonzero(frame > 0)


def get_volume_percentage(frame):
    return 100.0*len(frame[frame == True])/len(frame.flatten())


def do_analysis_on_all_files(parameters, key_file):
    dim_folder = get_dimension_folder(parameters)
    seg_method = parameters['segmentation_method']

    filenames = []
    all_tumor_volumes = []
    all_macrophage_volumes = []
    all_macrophage_volumes_labeled = []
    all_macrophage_number = []
    all_mean_macrophage_to_tumor_distances = []
    all_mean_macrophage_to_tumor_distances_labeled = []
    all_vessel_volumes = []
    all_mean_macrophage_to_vessel_distances = []

    all_macrophage_percentages = []
    all_tumor_percentages = []
    all_vessel_percentages = []

    for file in key_file["New name"].unique():
        filenames.append(file)

        # only load mask if at least one channel is good enough
        if (key_file[key_file['New name'] == file][seg_method + '_C1'] == 1).any() or (key_file[key_file['New name'] == file][seg_method + '_C2'] == 1).any() or (key_file[key_file['New name'] == file][seg_method + '_C3'] == 1).any():
            path_for_checking_frames = parameters['output_folder'] + '05_Data_Analysis/' + dim_folder + '03_Check_masks/' + seg_method + '/' + file + '/'
            if not os.path.exists(path_for_checking_frames):
                os.makedirs(path_for_checking_frames)

            # load mask or ilastik file
            print('Loading mask ' + file + '...')
            if seg_method == 'ilastik':
                file_path = get_segmentation_file_path(parameters)
                tumor_file = file_path + file + '_C1_Probabilities.h5'
                macrophage_file = file_path + file + '_C2_Probabilities.h5'
                vessel_file = file_path + file + '_C3_Probabilities.h5'

                with h5py.File(tumor_file, "r") as f:
                    a_group_key = list(f.keys())[0]
                    # Get the data
                    tumor_h5 = np.array(f[a_group_key])[:, :, :, :, 0]
                #tumor_data = tumor_h5[:, :, :, 0]

                with h5py.File(macrophage_file, "r") as f:
                    a_group_key = list(f.keys())[0]
                    # Get the data
                    macrophage_h5 = np.array(f[a_group_key])[:, :, :, :, 0]
                # macrophage_data = macrophage_h5[:, :, :, 0]

                with h5py.File(vessel_file, "r") as f:
                    a_group_key = list(f.keys())[0]
                    # Get the data
                    vessel_h5 = np.array(f[a_group_key])[:, :, :, :, 0]

                # only use probabilities higher than 0.5
                tumor_mask = np.where(np.array(tumor_h5) > 0.5, True, False)
                macrophage_mask = np.where(np.array(macrophage_h5) > 0.9, True, False)
                vessel_mask = np.where(np.array(vessel_h5) > 0.5, True, False)

            else:
                mask = load_mask(parameters, file + '.tif')
                tumor_mask = np.where(np.array(mask[:, :, :, :, 0]) > 0, True, False)
                macrophage_mask = np.where(np.array(mask[:, :, :, :, 1]) > 0, True, False)
                vessel_mask = np.where(np.array(mask[:, :, :, :, 2]) > 0, True, False)

            ##### Analysis functions

            # check if tumor mask looks okay
            if (key_file[key_file['New name'] == file][seg_method + '_C1'] == 1).any():
                print('    Get tumor volumes...')
                tumor_volumes = []
                tumor_percentages = []
                for frame in range(tumor_mask.shape[0]):
                    # save that frame for checking mask by eye...
                    fig, ax = plt.subplots(figsize=(15, 15))
                    ax.imshow(tumor_mask[frame])
                    plt.savefig(path_for_checking_frames + 'tumor_frame' + str(frame) + '.png', format="png", bbox_inches="tight", dpi=100)
                    plt.close()

                    plot_distance_map(tumor_mask[frame], path_for_checking_frames + 'tumor_distance_map_frame' + str(frame) + '.png')

                    tumor_volumes.append(get_volume(tumor_mask[frame]))
                    tumor_percentages.append(get_volume_percentage(tumor_mask[frame]))
                all_tumor_volumes.append(tumor_volumes)
                all_tumor_percentages.append(tumor_percentages)
            else:
                all_tumor_volumes.append(None)
                all_tumor_percentages.append(None)

            # check if macrophage mask looks okay
            if (key_file[key_file['New name'] == file][seg_method + '_C2'] == 1).any():
                labeled_macrophage_mask = get_labeled_macrophage_image(parameters, macrophage_mask)
                print('    Get macrophage volumes and number...')
                macrophage_volumes = []
                macrophage_volumes_labeled = []
                macrophage_number = []
                macrophage_percentages = []
                for frame in range(macrophage_mask.shape[0]):
                    fig, ax = plt.subplots(figsize=(15, 15))
                    ax.imshow(macrophage_mask[frame, :, :])
                    plt.savefig(path_for_checking_frames + 'macrophage_cutoff90percent_frame' + str(frame) + '.png', format="png",
                                 bbox_inches="tight", dpi=100)
                    plt.close()

                    fig, ax = plt.subplots(figsize=(15, 15))
                    ax.imshow(labeled_macrophage_mask[frame, :, :])
                    plt.savefig(path_for_checking_frames + 'macrophage_labeled_cutoff90percent_frame' + str(frame) + '.png',
                                format="png",
                                bbox_inches="tight", dpi=100)
                    plt.close()

                    macrophage_volumes.append(get_volume(macrophage_mask[frame]))
                    macrophage_volumes_labeled.append(get_volume(labeled_macrophage_mask[frame]))
                    macrophage_number.append(get_macrophage_number(labeled_macrophage_mask[frame]))
                    macrophage_percentages.append(get_volume_percentage(macrophage_mask[frame]))
                all_macrophage_volumes.append(macrophage_volumes)
                all_macrophage_volumes_labeled.append(macrophage_volumes_labeled)
                all_macrophage_number.append(macrophage_number)
                all_macrophage_percentages.append(macrophage_percentages)
            else:
                all_macrophage_volumes.append(None)
                all_macrophage_volumes_labeled.append(None)
                all_macrophage_number.append(None)
                all_macrophage_percentages.append(None)

            # check if vessel mask looks okay
            if (key_file[key_file['New name'] == file][seg_method + '_C3'] == 1).any():
                print('    Get vessel volumes...')
                vessel_volumes = []
                vessel_percentages = []
                for frame in range(vessel_mask.shape[0]):
                    fig, ax = plt.subplots(figsize=(15, 15))
                    ax.imshow(vessel_mask[frame, :, :])
                    plt.savefig(path_for_checking_frames + 'vessel_frame' + str(frame) + '.png', format="png",
                                bbox_inches="tight", dpi=100)
                    plt.close()

                    plot_distance_map(vessel_mask[frame], path_for_checking_frames + 'vessel_distance_map_frame' + str(frame) + '.png')

                    vessel_volumes.append(get_volume(vessel_mask[frame]))
                    vessel_percentages.append(get_volume_percentage(vessel_mask[frame]))
                all_vessel_volumes.append(vessel_volumes)
                all_vessel_percentages.append(vessel_percentages)
            else:
                all_vessel_volumes.append(None)
                all_vessel_percentages(None)

            # if tumor and macrophage mask look okay, go for the distance transform :)
            if (key_file[key_file['New name'] == file][seg_method + '_C1'] == 1).any() and \
                    (key_file[key_file['New name'] == file][seg_method + '_C2'] == 1).any():
                print('    Get macrophage to tumor distances...')
                mean_macrophage_to_tumor_distances = []
                mean_macrophage_to_tumor_distances_labeled = []
                for frame in range(macrophage_mask.shape[0]):
                    print(frame)
                    # print(macrophage_mask[frame].shape)
                    mean_macrophage_to_tumor_distances.append(calculate_mean_distance_of_macrophages(tumor_mask[frame],
                                                                                                  macrophage_mask[frame]))
                    # TO DO: I also tried to get the distances with the labeled image but somehow it didn't work
                    macrophage_labeled_boolean = np.where(np.array(labeled_macrophage_mask[frame]) > 0, True, False)
                    mean_macrophage_to_tumor_distances_labeled.append(calculate_mean_distance_of_macrophages(tumor_mask[frame], macrophage_labeled_boolean))
                all_mean_macrophage_to_tumor_distances.append(mean_macrophage_to_tumor_distances)
                all_mean_macrophage_to_tumor_distances_labeled.append(mean_macrophage_to_tumor_distances_labeled)

            else:
                all_mean_macrophage_to_tumor_distances.append(None)
                all_mean_macrophage_to_tumor_distances_labeled.append(None)

            # if vessel and macrophage masks look okay, go fo the distance transform :)
            if (key_file[key_file['New name'] == file][seg_method + '_C2'] == 1).any() and (key_file[key_file['New name'] == file][seg_method + '_C3'] == 1).any():
                print('    Get macrophage to vessel distances...')
                mean_macrophage_to_vessel_distances = []
                for frame in range(macrophage_mask.shape[0]):
                    print(frame)
                    mean_macrophage_to_vessel_distances.append(calculate_mean_distance_of_macrophages(vessel_mask[frame], macrophage_mask[frame]))
                all_mean_macrophage_to_vessel_distances.append(mean_macrophage_to_vessel_distances)
            else:
                all_mean_macrophage_to_vessel_distances.append(None)

        else:
            print('Looks like the masks of ' + file + ' sucked for the ' + parameters['segmentation_method'] + ' segmentation :(')
            all_tumor_volumes.append(None)
            all_macrophage_volumes.append(None)
            all_macrophage_volumes_labeled.append(None)
            all_macrophage_number.append(None)
            all_mean_macrophage_to_tumor_distances.append(None)
            all_mean_macrophage_to_tumor_distances_labeled.append(None)
            all_vessel_volumes.append(None)
            all_mean_macrophage_to_vessel_distances.append(None)

            all_tumor_percentages.append(None)
            all_macrophage_percentages.append(None)
            all_vessel_percentages.append(None)

    df = pd.DataFrame({
        'New name': filenames,
        'tumor_volume': all_tumor_volumes,
        'tumor_volume_percentage': all_tumor_percentages,
        'macrophage_volume': all_macrophage_volumes,
        'macrophage_volume_labeled': all_macrophage_volumes_labeled,
        'macrophage_number': all_macrophage_number,
        'macrophage_volume_percentage': all_macrophage_percentages,
        'vessel_volume': all_vessel_volumes,
        'vessel_volume_percentage': all_vessel_percentages,
        'mean_macrophage_to_tumor_distances': all_mean_macrophage_to_tumor_distances,
        'mean_macrophage_to_vessel_distances': all_mean_macrophage_to_vessel_distances,
        'mean_macrophage_to_tumor_distances_labeled': all_mean_macrophage_to_tumor_distances_labeled
    })

    path_for_saving_df = parameters['output_folder'] + '05_Data_Analysis/' + dim_folder + '01_Dataframes/' + seg_method + '/'
    if not os.path.exists(path_for_saving_df):
        os.makedirs(path_for_saving_df)
    df.to_pickle(path_for_saving_df + 'analysis_dataframe.pkl')

    return df
