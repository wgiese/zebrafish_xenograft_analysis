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


def tumor_volume(tumor_mask):
    return np.count_nonzero(tumor_mask > 0)


def macrophage_volume(macrophage_mask):
    return np.count_nonzero(macrophage_mask > 0)


def macrophage_number(parameters, macrophage_mask):
    # remove artifacts connected to image border
    cleared = clear_border(macrophage_mask)
    cleared = morphology.remove_small_objects(cleared, parameters['macrophages_small_objects'], connectivity=2)

    # label image regions
    label_image = label(cleared)

    return np.max(label_image)


def do_analysis_on_all_files(parameters, key_file):
    filenames = []
    all_tumor_volumes = []
    all_macrophage_volumes = []
    all_macrophage_number = []
    all_mean_macrophage_to_tumor_distances = []

    for file in key_file["New name"].unique():
        filenames.append(file)
        # check if mask is good enough
        seg_method = parameters['segmentation_method']
        if key_file[file][seg_method + '_C1'] == 1 or key_file[file][seg_method + '_C2'] == 1 or key_file[file][seg_method + '_C3'] == 1:
            # load mask or ilastik file
            print('Loading mask...')
            if seg_method == 'ilastik':
                tumor_file = load_mask(parameters, file + '_C1')
                macrophage_file = load_mask(parameters, file + '_C2')
                vessel_file = load_mask(parameters, file + '_C3')

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
                mask = load_mask(parameters, file)
                tumor_mask = mask[:, :, :, 0]
                macrophage_mask = mask[:, :, :, 1]
                vessel_mask = mask[:, :, :, 2]

            ##### Analysis functions

        else:
            print('Looks like the masks of ' + file + 'sucked for the ' + parameters['segmentation_method'] + ' segmentation :(')
            all_tumor_volumes.append(None)
            all_macrophage_volumes.append(None)
            all_macrophage_number.append(None)
            all_mean_macrophage_to_tumor_distances.append(None)

    return
