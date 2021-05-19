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
#                                            OTSU THRESHOLDING FUNCTIONS                                              #
#######################################################################################################################

def otsu_thresholding_3D(parameters, image):
    # use Gaussian filter to smooth images
    # print(" Gaussian filtering...")
    image_labeled = gaussian_filter(image, sigma=parameters['sigma'])

    # Otsu thresholding images
    # print(" Otsu thresholding...")
    image_labeled = np.where(image_labeled > threshold_otsu(image_labeled), True, False)

    # remove artefacts connected to border
    # print(" Removing small objects...")
    image_labeled = clear_border(image_labeled)
    image_labeled = morphology.remove_small_objects(image_labeled, 2500, connectivity=2)

    # label image regions
    image_labeled = label(image_labeled)

    return image_labeled


def otsu_thresholding_3D_allFrames(parameters, movie):
    im_tumor = movie[:, :, :, :, 0]
    im_macrophages = movie[:, :, :, :, 1]
    im_vessels = movie[:, :, :, :, 2]

    mask_tumor = np.zeros(im_tumor.shape)
    mask_macrophages = np.zeros(im_tumor.shape)
    mask_vessels = np.zeros(im_tumor.shape)
    # mask = np.zeros(movie.shape)

    print("Total number of time points: " + str(movie.shape[0]))

    for tp in range(movie.shape[0]):
        print("Time point: " + str(tp))
        # mask[tp] = otsu_thresholding_3D(parameters, movie[tp])
        mask_tumor[tp] = otsu_thresholding_3D(parameters, im_tumor[tp])
        mask_macrophages[tp] = otsu_thresholding_3D(parameters, im_macrophages[tp])
        mask_vessels[tp] = otsu_thresholding_3D(parameters, im_vessels[tp])

    return mask_tumor, mask_macrophages, mask_vessels
    # return mask


def otsu_thresholding_2D_allFrames(parameters, movie2D):
    im_tumor = movie2D[:, :, :, 0]
    im_macrophages = movie2D[:, :, :, 1]
    im_vessels = movie2D[:, :, :, 2]

    mask_tumor = np.zeros(im_tumor.shape)
    mask_macrophages = np.zeros(im_tumor.shape)
    mask_vessels = np.zeros(im_tumor.shape)
    # mask = np.zeros(movie.shape)

    print("Total number of time points: " + str(movie2D.shape[0]))

    for tp in range(movie2D.shape[0]):
        print("Time point: " + str(tp))
        # mask[tp] = otsu_thresholding_3D(parameters, movie[tp])
        mask_tumor[tp] = otsu_thresholding_3D(parameters, im_tumor[tp])
        mask_macrophages[tp] = otsu_thresholding_3D(parameters, im_macrophages[tp])
        mask_vessels[tp] = otsu_thresholding_3D(parameters, im_vessels[tp])

    return mask_tumor, mask_macrophages, mask_vessels
    # return mask


def otsu_loop_over_key_file(parameters, key_file):
    """
    NEEDS IMPROVEMENT! Otsu masks are not saved yet. Careful with 2D/3D. It doesn't read that information from the
    parameters file yet.
    :param parameters:
    :param key_file:
    :return:
    """
    for file in key_file["New name"].unique():
        # file_path = parameters["data_folder"] + "03_Processed_Data/3D/" + file + '.tif'
        file_path = parameters["data_folder"] + "05_BGsubtracted/02_3D/" + file + '.tif'
        if os.path.exists(file_path):
            print(file_path)
            output_path = parameters["output_folder"] + "05_BGsubtracted/02_3D/Otsu/sigma" + str(
                parameters['sigma']) + '/' + file + '/'
            if not os.path.exists(output_path):
                os.makedirs(output_path)

                print("Loading data...")
                movie = np.array(io.imread(file_path))

                # im_tumor = movie[:, :, :, :, 0]
                # im_macrophages = movie[:, :, :, :, 1]
                # im_vessels = movie[:, :, :, :, 2]

                # Otsu thresholding
                print(" Otsu thresholding...")
                [mask_tumor, mask_macrophages, mask_vessels] = otsu_thresholding_3D_allFrames(parameters, movie)
                # [mask_2D_tumor, mask_2D_macrophages, mask_2D_vessels] = otsu_thresholding_2D_allFrames(parameters, movie)
                # print("     Tumor")
                # mask_tumor = otsu_thresholding_3D_allFrames(parameters, im_tumor)
                # print("     Macrophages")
                # mask_macrophages = otsu_thresholding_3D_allFrames(parameters, im_macrophages)
                # print("     Vessels")
                # mask_vessels = otsu_thresholding_3D_allFrames(parameters, im_vessels)

                # save 2D mask projections to file
                mask_2D_tumor = np.zeros((mask_tumor.shape[0], mask_tumor.shape[2], mask_tumor.shape[3]))
                mask_2D_macrophages = np.zeros((mask_tumor.shape[0], mask_tumor.shape[2], mask_tumor.shape[3]))
                mask_2D_vessels = np.zeros((mask_tumor.shape[0], mask_tumor.shape[2], mask_tumor.shape[3]))

                for tp in range(mask_tumor.shape[0]):
                    mask_2D_tumor[tp] = np.amax(mask_tumor[tp], axis=0)
                    mask_2D_macrophages[tp] = np.amax(mask_macrophages[tp], axis=0)
                    mask_2D_vessels[tp] = np.amax(mask_vessels[tp], axis=0)

                # numRows = mask_2D_tumor.shape[0]
                # save image of comparison between mask and original to assess how well segmentation works
                for tp in range(mask_2D_tumor.shape[0]):
                    fig, axes = plt.subplots(1, 2, figsize=(16, 12))
                    axes[0].imshow(np.amax(movie[tp, :, :, :, 0], axis=0))
                    # axes[0].imshow(movie[tp, :, :, 0])
                    axes[1].imshow(mask_2D_tumor[tp])
                    plt.savefig(output_path + "/tumor_" + str(tp) + ".png", format="png", bbox_inches="tight", dpi=150)
                    plt.close()

                    fig, axes = plt.subplots(1, 2, figsize=(16, 12))
                    axes[0].imshow(np.amax(movie[tp, :, :, :, 1], axis=0))
                    # axes[0].imshow(movie[tp, :, :, 1])
                    axes[1].imshow(mask_2D_macrophages[tp])
                    plt.savefig(output_path + "/macrophages_" + str(tp) + ".png", format="png", bbox_inches="tight",
                                dpi=150)
                    plt.close()

                    fig, axes = plt.subplots(1, 2, figsize=(16, 12))
                    axes[0].imshow(np.amax(movie[tp, :, :, :, 2], axis=0))
                    # axes[0].imshow(movie[tp, :, :, 2])
                    axes[1].imshow(mask_2D_vessels[tp])
                    plt.savefig(output_path + "/vessels_" + str(tp) + ".png", format="png", bbox_inches="tight",
                                dpi=150)
                    plt.close()

                # del movie
                # del mask_macrophages, mask_tumor, mask_vessels

            else:
                print("Otsu thresholding already done for " + file + " and sigma " + str(parameters['sigma']))

        else:
            print("The file '%s' does not exist." % file_path)

    return None


#######################################################################################################################
#######################################################################################################################


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


def calculate_tumor_macrophage_distance(tumor_mask, macrophage_mask):
    tumor_distance = scipy.ndimage.morphology.distance_transform_edt(np.invert(tumor_mask[:, :, :]))
    return np.mean(tumor_distance[macrophage_mask.flatten()])


def get_tumor_volume(tumor_mask):
    return np.count_nonzero(tumor_mask > 0)


def get_labeled_macrophage_image(parameters, macrophage_mask):
    all_frames_labeled = np.zeros(macrophage_mask.shape)

    for frame in range(macrophage_mask.shape[0]):
        # remove artifacts connected to image border
        cleared = clear_border(macrophage_mask[frame])
        cleared = morphology.remove_small_objects(cleared, parameters['macrophages_small_objects'], connectivity=2)

        # label image regions
        all_frames_labeled[frame] = label(cleared)

    return all_frames_labeled


def get_macrophage_volume(macrophage_mask):
    return np.count_nonzero(macrophage_mask > 0)


def get_macrophage_number(labeled_macrophage_mask):
    return np.max(labeled_macrophage_mask)


def do_analysis_on_all_files(parameters, key_file):
    filenames = []
    all_tumor_volumes = []
    all_macrophage_volumes = []
    all_macrophage_volumes_labeled = []
    all_macrophage_number = []
    all_mean_macrophage_to_tumor_distances = []
    all_mean_macrophage_to_tumor_distances_labeled = []

    for file in key_file["New name"].unique():
        filenames.append(file)
        # check if mask is good enough
        seg_method = parameters['segmentation_method']

        # only load mask if at least one channel is good enough
        if (key_file[key_file['New name'] == file][seg_method + '_C1'] == 1).any() or (key_file[key_file['New name'] == file][seg_method + '_C2'] == 1).any() or (key_file[key_file['New name'] == file][seg_method + '_C3'] == 1).any():
            # load mask or ilastik file
            print('Loading mask ' + file + '...')
            if seg_method == 'ilastik':
                tumor_file = load_mask(parameters, file + '_C1.tif')
                macrophage_file = load_mask(parameters, file + '_C2.tif')
                vessel_file = load_mask(parameters, file + '_C3.tif')

                with h5py.File(tumor_file, "r") as f:
                    a_group_key = list(f.keys())[0]
                    # Get the data
                    tumor_h5 = np.array(f[a_group_key])[:, :, :, 0]
                #tumor_data = tumor_h5[:, :, :, 0]

                with h5py.File(macrophage_file, "r") as f:
                    a_group_key = list(f.keys())[0]
                    # Get the data
                    macrophage_h5 = np.array(f[a_group_key])[:, :, :, 0]
                # macrophage_data = macrophage_h5[:, :, :, 0]

                with h5py.File(vessel_file, "r") as f:
                    a_group_key = list(f.keys())[0]
                    # Get the data
                    vessel_h5 = np.array(f[a_group_key])[:, :, :, 0]

                # only use probabilities higher than 0.5
                tumor_mask = np.where(tumor_h5 > 0.5, True, False)
                macrophage_mask = np.where(macrophage_h5 > 0.5, True, False)
                vessel_mask = np.where(vessel_h5 > 0.5, True, False)

            else:
                mask = load_mask(parameters, file + '.tif')
                tumor_mask = mask[:, :, :, 0]
                macrophage_mask = mask[:, :, :, 1]
                vessel_mask = mask[:, :, :, 2]

            ##### Analysis functions

            # check if tumor mask looks okay
            if (key_file[key_file['New name'] == file][seg_method + '_C1'] == 1).any():
                print('    Get tumor volumes...')
                tumor_volumes = []
                for frame in range(tumor_mask.shape[0]):
                    tumor_volumes.append(get_tumor_volume(tumor_mask[frame]))
                all_tumor_volumes.append(tumor_volumes)
            else:
                all_tumor_volumes.append(None)

            # check if macrophage mask looks okay
            if (key_file[key_file['New name'] == file][seg_method + '_C2'] == 1).any():
                labeled_macrophage_mask = get_labeled_macrophage_image(parameters, macrophage_mask)
                print('    Get macrophage volumes and number...')
                macrophage_volumes = []
                macrophage_volumes_labeled = []
                macrophage_number = []
                for frame in range(macrophage_mask.shape[0]):
                    macrophage_volumes.append(get_macrophage_volume(macrophage_mask[frame]))
                    macrophage_volumes_labeled.append(get_macrophage_volume(labeled_macrophage_mask[frame]))
                    macrophage_number.append(get_macrophage_number(labeled_macrophage_mask[frame]))
                all_macrophage_volumes.append(macrophage_volumes)
                all_macrophage_volumes_labeled.append(macrophage_volumes_labeled)
                all_macrophage_number.append(macrophage_number)
            else:
                all_macrophage_volumes.append(None)
                all_macrophage_volumes_labeled.append(None)
                all_macrophage_number.append(None)

            # if both of them look okay, go for the distance transform :)
            if (key_file[key_file['New name'] == file][seg_method + '_C1'] == 1).any() and (key_file[key_file['New name'] == file][seg_method + '_C2'] == 1).any():
                print('    Get macrophage to tumor distances...')
                mean_macrophage_to_tumor_distances = []
                mean_macrophage_to_tumor_distances_labeled = []
                for frame in range(macrophage_mask.shape[0]):
                    print(frame)
                    mean_macrophage_to_tumor_distances.append(calculate_tumor_macrophage_distance(tumor_mask[frame], macrophage_mask[frame]))
                    mean_macrophage_to_tumor_distances_labeled.append(
                        calculate_tumor_macrophage_distance(tumor_mask[frame], labeled_macrophage_mask[frame]))
                all_mean_macrophage_to_tumor_distances.append(mean_macrophage_to_tumor_distances)
                all_mean_macrophage_to_tumor_distances_labeled.append(mean_macrophage_to_tumor_distances_labeled)

            else:
                all_mean_macrophage_to_tumor_distances.append(None)
                all_mean_macrophage_to_tumor_distances_labeled.append(None)

        else:
            print('Looks like the masks of ' + file + ' sucked for the ' + parameters['segmentation_method'] + ' segmentation :(')
            filenames.append(None)
            all_tumor_volumes.append(None)
            all_macrophage_volumes.append(None)
            all_macrophage_volumes_labeled.append(None)
            all_macrophage_number.append(None)
            all_mean_macrophage_to_tumor_distances.append(None)
            all_mean_macrophage_to_tumor_distances_labeled.append(None)

    df = pd.DataFrame({
        'file_name': filenames,
        'tumor_volume': all_tumor_volumes,
        'macrophage_volume': all_macrophage_volumes,
        'macrophage_volume_labeled': all_macrophage_volumes_labeled,
        'macrophage_number': all_macrophage_number,
        'mean_macrophage_to_tumor_distances': all_mean_macrophage_to_tumor_distances,
        'mean_macrophage_to_tumor_distances': all_mean_macrophage_to_tumor_distances_labeled
    })

    dim_folder = get_dimension_folder(parameters)
    path_for_saving_df = parameters['output_folder'] + '05_Data_Analysis/' + dim_folder + '01_Dataframes/' + seg_method + '/'
    if not os.path.exists(path_for_saving_df):
        os.makedirs(path_for_saving_df)
    df.to_pickle(path_for_saving_df + 'analysis_dataframe.pkl')

    return df
