import h5py
import numpy as np
from pyevtk.hl import imageToVTK
from matplotlib import pyplot as plt
from matplotlib import cm
import seaborn as sns
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
from scipy import ndimage as ndi
from skimage.segmentation import watershed
from skimage.feature import peak_local_max


font = {'family': 'normal', 'weight': 'bold', 'size': 22}
matplotlib.rc('font', **font)


import functions_common


#######################################################################################################################
#                                                 Analysis functions                                                  #
#######################################################################################################################

def otsu_thresholding_3D(parameters, image):
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
    image_blurred = gaussian_filter(image, sigma=parameters['sigma'])

    # Otsu thresholding images
    # print(" Otsu thresholding...")
    thresh = threshold_otsu(image_blurred)

    #image_mask = np.where(image_blurred > thresh, True, False)
    image_mask = image_blurred > thresh
    del image_blurred

    # remove artefacts connected to border
    # print(" Removing small objects...")
    image_mask = clear_border(image_mask)
    image_mask = morphology.remove_small_objects(image_mask, parameters["macrophages_small_objects"], connectivity=2)

    # label image regions
    image_labeled = label(image_mask)
    del image_mask

    return image_labeled, thresh


def otsu_thresholding_3D_allFrames(parameters, movie):
    '''
        input:
            parameters  - dict() with all meta parameters for analysis
            movie - numpy array containing 3D + time image data with dimensions in the following order time, x , y, z 
        output:
            tumor_labels - numpy array 3D + time with tumor labels
            macrophage_labels - numpy array 3D + time with macrophage labels
            vessel_labels - numpy array 3D + time with vessel labels

        TODO: - split into thresholding and labelling
              - move to common if used in 2D nd 3D
              - think about doing channels separately due to high memory demand in 3D
    '''

    im_tumor = movie[:, :, :, :, 0]
    im_macrophages = movie[:, :, :, :, 1]
    im_vessels = movie[:, :, :, :, 2]

    tumor_labels = np.zeros(im_tumor.shape)
    macrophage_labels = np.zeros(im_tumor.shape)
    vessel_labels = np.zeros(im_tumor.shape)
    # mask = np.zeros(movie.shape)

    print("Total number of time points: " + str(movie.shape[0]))

    for tp in range(movie.shape[0]):
        print("Time point: " + str(tp))
        # mask[tp] = otsu_thresholding_3D(parameters, movie[tp])
        tumor_labels[tp], thresh = otsu_thresholding_3D(parameters, im_tumor[tp])
        macrophages_labels[tp], thresh = otsu_thresholding_3D(parameters, im_macrophages[tp])
        vessel_labels[tp], thresh = otsu_thresholding_3D(parameters, im_vessels[tp])

    return tumor_labels, macrophages_labels, vessel_labels
    # return mask


def otsu_thresholding_2D_allFrames(parameters, movie2D):
    im_tumor = movie2D[:, :, :, 0]
    im_macrophages = movie2D[:, :, :, 1]
    im_vessels = movie2D[:, :, :, 2]

    mask_tumor = np.zeros(im_tumor.shape)
    mask_macrophages = np.zeros(im_tumor.shape)
    mask_vessels = np.zeros(im_tumor.shape)
    # mask = np.zeros(movie.shape)
    
    #TODO: Check if gaussian smoothing is helpful here

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
    for filename in key_file["short_name"].unique():
        if pd.isna(filename):
            continue
        print(filename)
        file_path = parameters["data_folder"] + "03_Processed_Data/3D/" + filename + '.tif'
        #file_path = parameters["data_folder"] + "05_BGsubtracted/02_3D/" + file + '.tif'

        if os.path.exists(file_path):
            print(file_path)
            #output_path = parameters["output_folder"] + "05_BGsubtracted/02_3D/Otsu/sigma" + str(
            #    parameters['sigma']) + '/' + file + '/'
            #if not os.path.exists(output_path):
            #    os.makedirs(output_path)

            print("Loading data...")
            movie = np.array(io.imread(file_path))

            # Otsu thresholding
            print(" Otsu thresholding...")
            [mask_tumor, mask_macrophages, mask_vessels] = otsu_thresholding_3D_allFrames(parameters, movie)
            
            del movie
            del mask_macrophages, mask_tumor, mask_vessels


            #if parameters["export_masks"] == "vtk":
    		#time_stamp = str(frame).zfill(3)
    		#filename = output_folder + filename_short + "-" + time_stamp
    
    		#imageToVTK(filename, cellData = {"tumor": im_tumor_smoothed, 
               #                      "macrophages" : im_macrophages_smoothed, 
               #                      "vessel": im_vessel_smoothed})
            	
            
            
                # [mask_2D_tumor, mask_2D_macrophages, mask_2D_vessels] = otsu_thresholding_2D_allFrames(parameters, movie)
                # print("     Tumor")
                # mask_tumor = otsu_thresholding_3D_allFrames(parameters, im_tumor)
                # print("     Macrophages")
                # mask_macrophages = otsu_thresholding_3D_allFrames(parameters, im_macrophages)
                # print("     Vessels")
                # mask_vessels = otsu_thresholding_3D_allFrames(parameters, im_vessels)

                # save 2D mask projections to file
            #    mask_2D_tumor = np.zeros((mask_tumor.shape[0], mask_tumor.shape[2], mask_tumor.shape[3]))
            #    mask_2D_macrophages = np.zeros((mask_tumor.shape[0], mask_tumor.shape[2], mask_tumor.shape[3]))
            #    mask_2D_vessels = np.zeros((mask_tumor.shape[0], mask_tumor.shape[2], mask_tumor.shape[3]))

              
                # numRows = mask_2D_tumor.shape[0]
                # save image of comparison between mask and original to assess how well segmentation works
                
            '''
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
            '''

                # del movie
                # del mask_macrophages, mask_tumor, mask_vessels

            #else:
            #    print("Otsu thresholding already done for " + file + " and sigma " + str(parameters['sigma']))

        else:
            print("The file '%s' does not exist." % file_path)

    return None

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


def get_macrophage_properties(parameters, key_file, experiment = "all", vtk_out = False, npy_out = False, proj_2D = True, hist = True):
    '''
        Input:
            parameters  - dict() with all meta parameters for analysis
            key_file - pandas data frame with relative filepaths to images and necessary meta info


    '''
    properties_df = pd.DataFrame()
    index_counter = 0

    print("Experiment(s): %s" % experiment)

    for index, row in key_file.iterrows():
        
        filename = row["short_name"]
        #print(filename)
        #print(experiment)

        if pd.isna(filename):
            continue
        if experiment == "annotated":
            if row["macrophages_annotated"] == 0:
                continue

        properties_df.at[index_counter, "short_name"] = filename
        properties_df.at[index_counter, "status"] = "not processed"

 

        if (experiment in ["all", "annotated", filename]) or (experiment in filename):
            file_path = parameters["data_folder"] + parameters["image_file_path"] + filename + '.tif'
            #file_path = parameters["data_folder"] + "02_Primary_Data/" + filename + '.tif'

            
            if np.isnan(row["PixelSizeX"]):
                properties_df.at[index_counter, "status"] = "pixel size is missing"
                index_counter += 1
                continue

            volume_conv_px_to_mum3 = row["PixelSizeX"]*row["PixelSizeY"]*row["PixelSizeZ"]
            print("Processing experiment: %s" % filename)

            if os.path.exists(file_path):
                print(file_path)

                print("Loading data...")
                movie = np.array(io.imread(file_path))
                movie_macrophages = movie[:, :, :, :, parameters["channel_macrophages"]]
                if parameters["substract_tumor_from_macrophages"]:
                    movie_tumor = movie[:, :, :, :, parameters["channel_tumor"]]

                del movie
                
                print("Shape of the macrophage movie: ")
                print(movie_macrophages.shape)                

                #mask_macrophages = np.zeros(movie_macrophages.shape)
                #labeled_macrophages = np.zeros(movie_macrophages.shape)

                for tp in range(movie_macrophages.shape[0]):
                    print("Time point: " + str(tp))
                    if parameters["substract_background"] in ["rolling_ball","ellipsoid_kernel"]:
                        macrophages_bg = functions_common.background_substraction(parameters, movie_macrophages[tp], row["PixelSizeZ"], row["PixelSizeX"], row["PixelSizeY"])
                    else:
                        macrophages_bg = movie_macrophages[tp]


                    #macrophages_blurred = gaussian_filter(macrophages_bg, sigma=parameters['sigma'])
                    sigma_z = row["PixelSizeZ"]*parameters['sigma_macrophage']
                    sigma_x = row["PixelSizeX"]*parameters['sigma_macrophage']
                    sigma_y = row["PixelSizeY"]*parameters['sigma_macrophage']
                    macrophages_blurred = gaussian_filter(macrophages_bg, sigma=[sigma_z,sigma_x,sigma_y])
                    macrophages_thresh, threshold = functions_common.thresholding_3D(parameters, macrophages_blurred)
                    
                    

                    #del macrophages_bg

                    if parameters["substract_tumor_from_macrophages"]:
                        print("Substracting tumor from macrophages ...")
                        tumor_blurred = gaussian_filter(movie_tumor[tp], sigma=1)
                        tumor_thresh, threshold_ = functions_common.thresholding_3D(parameters, tumor_blurred)
                        macrophages_thresh[tumor_thresh] = False   

                    #labeled_macrophages[tp], threshold = otsu_thresholding_3D(parameters, movie_macrophages[tp])
                    
                    if parameters["segmentation_method"] == "Watershed":
                        print("Apply distance transform ...")
                        distance = ndi.distance_transform_edt(macrophages_thresh)
                        coords = peak_local_max(distance, footprint=np.ones((3, 3, 3)), labels=macrophages_thresh)
                        mask = np.zeros(distance.shape, dtype=bool)
                        mask[tuple(coords.T)] = True
                        markers, num_labels = ndi.label(mask)
                        print("Apply watershed ...")
                        labeled_macrophages = watershed(-distance, markers, mask=macrophages_thresh)
                    else:
                        #image_mask = clear_border(image_mask)
                        #image_mask = morphology.remove_small_objects(image_mask, parameters["macrophages_small_objects"], connectivity=2)
                        print("Apply skimage labelling to threshold mask ...")
                        macrophages_thresh = clear_border(macrophages_thresh)
                        macrophages_thresh = morphology.remove_small_objects(macrophages_thresh, parameters["macrophages_small_objects"]/volume_conv_px_to_mum3) #, connectivity=2)
                        labeled_macrophages = label(macrophages_thresh)
                        num_labels = np.max(labeled_macrophages)

                    print("Number of labels: %s" % num_labels)                    

                    print("Extract coordinates")
                    x_centroids = []
                    y_centroids = []
                    for label_id in range(1, num_labels + 1):
                        macrophage_prop = regionprops(np.where(labeled_macrophages == label_id, 1, 0))
                        macrophage_volume_px = len(labeled_macrophages[labeled_macrophages == label_id].flatten())                       
                        macrophage_volume_mum3 = macrophage_volume_px*volume_conv_px_to_mum3                       

                        if (macrophage_volume_mum3 < parameters["macrophages_small_objects"]):
                            continue

                        for props in macrophage_prop:
                            z_centroid, x_centroid, y_centroid = props.centroid
                        x_centroids.append(x_centroid)
                        y_centroids.append(y_centroid)                    
    
                        properties_df.at[index_counter, "short_name"] = filename
                        properties_df.at[index_counter, "time_frame"] = tp
                        if parameters["thresholding_method"]  in ["niblack","sauvola"]:
                            properties_df.at[index_counter, "threshold"] = "local"
                        else:
                            properties_df.at[index_counter, "threshold"] = threshold
                        properties_df.at[index_counter, "macrophage_label"] = label_id
                        properties_df.at[index_counter, "macrophage_volume_px"] = macrophage_volume_px
                        properties_df.at[index_counter, "macrophage_volume_mum3"] = macrophage_volume_mum3
                        properties_df.at[index_counter, "x_centroid"] = x_centroid
                        properties_df.at[index_counter, "y_centroid"] = y_centroid
                        properties_df.at[index_counter, "z_centroid"] = z_centroid
                        properties_df.at[index_counter, "status"] = "OK"
                     
                        index_counter += 1
                    
                    time_stamp = "-" + str(tp).zfill(3)
                    
                    # debug output
                    if hist:                        
                        if parameters["thresholding_method"]  in ["niblack","sauvola"]:
                            fig, ax = plt.subplots(1,3, figsize=(45,15))
                        else:
                            fig, ax = plt.subplots(1,2, figsize=(30,15))

                        n_bins = np.max(movie_macrophages[tp].flatten())
                        sns.histplot(movie_macrophages[tp].flatten(), bins = n_bins + 1,ax = ax[0])
                        ax[0].set_xlim(0, n_bins + 1)
                        n_bins = np.max(macrophages_blurred.flatten())
                        sns.histplot(macrophages_blurred.flatten(), bins = n_bins + 1,ax = ax[1])
                        ax[1].set_xlim(0, n_bins + 1)
                        ax[0].set_ylim(0, len(movie_macrophages[tp][movie_macrophages[tp] > 0].flatten()))
                        ax[1].set_ylim(0, len(macrophages_blurred[macrophages_blurred > 0].flatten()))
                        #ax[0].set_ylim(0, 0.1*len(movie_macrophages[tp].flatten()))
                        #ax[1].set_ylim(0, 0.1*len(movie_macrophages[tp].flatten()))
                        ax[0].set_title("histogram original")
                        ax[1].set_title("histogram of blurred image")
                        if parameters["thresholding_method"]  in ["niblack","sauvola"]:
                            mean_2D_proj = np.zeros((labeled_macrophages.shape[1],labeled_macrophages.shape[2])) 
                            for x in range(labeled_macrophages.shape[1]):
                                for y in range(labeled_macrophages.shape[2]):
                                    mean_2D_proj[x,y] = np.mean(threshold[:,x,y])
                            ax[2].imshow(mean_2D_proj[:,:])
                            ax[2].set_title("thresholding filter")
                        else:
                            print("plot histogram with threshold threshold : %s" % threshold) 
                            ax[1].axvline(threshold, color = "r", linewidth = 3)
                            ax[1].set_title("histogram of blurred image with threshold at %s" % threshold)
                        plt.savefig(parameters["output_folder"] + filename + "-histogram" + time_stamp + ".pdf")
                        plt.savefig(parameters["output_folder"] + filename + "-histogram" + time_stamp + ".png")
                        plt.close()
                    if vtk_out:
                        imageToVTK(parameters["output_folder"] + filename + time_stamp, cellData = {"macrophages" : labeled_macrophages} )
                    
                    if proj_2D:
                        label_2D_proj = np.zeros((labeled_macrophages.shape[1],labeled_macrophages.shape[2])) 
                        sum_bg_2D_proj = np.zeros((labeled_macrophages.shape[1],labeled_macrophages.shape[2])) 
                        sum_2D_proj = np.zeros((labeled_macrophages.shape[1],labeled_macrophages.shape[2])) 
                        sum_blurred_2D_proj = np.zeros((labeled_macrophages.shape[1],labeled_macrophages.shape[2])) 
                        plot_df = properties_df[properties_df["time_frame"] == tp]
                        plot_df = plot_df[plot_df["short_name"] == filename]
                        for x in range(labeled_macrophages.shape[1]):
                            for y in range(labeled_macrophages.shape[2]):
                                label_2D_proj[x,y] = np.max(labeled_macrophages[:,x,y])
                                sum_2D_proj[x,y] = np.sum(movie_macrophages[tp][:,x,y])
                                sum_blurred_2D_proj[x,y] = np.sum(macrophages_blurred[:,x,y])
                                sum_bg_2D_proj[x,y] = np.sum(macrophages_bg[:,x,y])
                        
                        fig, ax = plt.subplots(1,4, figsize=(60,15))
                        p = ax[0].imshow(sum_2D_proj[:,:], cm.binary)
                        #cbar = fig.colorbar(p, ax = ax[0])
                        ax[1].imshow(sum_bg_2D_proj[:,:], cm.binary)
                        ax[2].imshow(sum_blurred_2D_proj[:,:], cm.binary)
                        ax[3].imshow(label_2D_proj[:,:], cm.Set3)
                        ax[0].set_title("sum projection - original")
                        ax[1].set_title("sum projection - background removed")
                        ax[2].set_title("sum projection - gaussian filter")
                        ax[3].set_title("labels from segmentation")
                        
                        
                        legend = 0

                        annotated_positions_file = parameters["data_folder"] + "04_Processed_Data/01_Annotated_Macrophages/" + filename + '.csv'
                        if os.path.exists(annotated_positions_file):
                            annotated_positions = pd.read_csv(annotated_positions_file, sep = ";")
                            print("annotated postions file exists ...")
                            print(annotated_positions.head())
                            plot_annotated_df = annotated_positions[annotated_positions["time_point"] == tp]
                            ax[3].plot(plot_annotated_df['X'], plot_annotated_df['Y'], 'go', markersize = 15)
                            if legend == 0:
                                ax[2].plot(plot_annotated_df['X'], plot_annotated_df['Y'], 'go', markersize = 15, label='macrophages by annotation')
                                legend = 1
                            ax[2].plot(plot_annotated_df['X'], plot_annotated_df['Y'], 'go', markersize = 15)

                        else: 
                            print("annotated postions file does not exists ...")
                        if len(plot_df.index) > 1:
                            legend = 0
                            for ind, row_plt in plot_df.iterrows():
                                #if row_plt["macrophage_volume"] < 100000:
                                #ax[1].plot(row_plt['y_centroid'], row_plt['x_centroid'], 'rX', markersize = 15)
                                if legend == 0:
                                    ax[2].plot(row_plt['y_centroid'], row_plt['x_centroid'], 'rX', markersize = 15, label='macrophages by algorithm')
                                    legend = 1
                                ax[2].plot(row_plt['y_centroid'], row_plt['x_centroid'], 'rX', markersize = 15)
                                #else:
                                #    ax[0].plot(row_plt['y_centroid'], row_plt['x_centroid'], 'rx', markersize = 15)
                                #    ax[1].plot(row_plt['y_centroid'], row_plt['x_centroid'], 'rx', markersize = 15)
                            #ax.plot(y_centroids, x_centroids, 'rx', markersize = 15)
                            ax[2].legend()
                        plt.tight_layout()
                        plt.savefig(parameters["output_folder"] + filename + time_stamp + ".pdf")
                        plt.savefig(parameters["output_folder"] + filename + time_stamp + ".png")
                        plt.close()
                    properties_df.to_csv(parameters["output_folder"] + "macrophage_properties_" + experiment + ".csv")

                    if npy_out:
                        np.save(parameters["output_folder"] + filename + time_stamp + ".npy", labeled_macrophages, allow_pickle=False)
            else:
                print("Path %s does not exist" % file_path)
                properties_df.at[index_counter, "status"] = "file not found"
                index_counter += 1

    return properties_df 

def get_tumor_macrophage_point_distances(parameters, key_file, macrophage_properties, experiment = "all", proj_2D = True):
    
    distances_df = pd.DataFrame()

    index_counter = 0

    for index, row in key_file.iterrows():

        filename = row["short_name"]
        print(filename)
        print(experiment)
        print(macrophage_properties["short_name"].unique())

        if pd.isna(filename):
            continue
        if not (filename in macrophage_properties["short_name"].unique()):
            continue
        
        print("read info")

        if (experiment == "all" or experiment == filename):

            file_path = parameters["data_folder"] + "03_Preprocessed_Data/02_3D/" + filename + '.tif'
            
            print("Row:")
            print(row)

            print("Loading data...")
            movie = np.array(io.imread(file_path))
            movie_tumor = movie[:, :, :, :, parameters["channel_tumor"]]
            
            single_macrophage_properties = macrophage_properties[macrophage_properties["short_name"]==filename]

            time_points = int(np.max(single_macrophage_properties["time_frame"]) + 1)

            for tp in range(time_points):

                sigma_z = row["PixelSizeZ"]*parameters['sigma_tumor']
                sigma_x = row["PixelSizeX"]*parameters['sigma_tumor']
                sigma_y = row["PixelSizeY"]*parameters['sigma_tumor']
                tumor_blurred = gaussian_filter(movie_tumor[tp], sigma=[sigma_z,sigma_x,sigma_y])
                tumor_thresh, threshold = functions_common.thresholding_3D(parameters, tumor_blurred)
               
                print("Compute distance transform ... ")
                dx = row["PixelSizeX"] 
                dy = row["PixelSizeY"] 
                dz = row["PixelSizeZ"] 
                tumor_distances = scipy.ndimage.morphology.distance_transform_edt(np.invert(tumor_thresh), sampling = [dz,dx,dy])
    
                mean_2D_proj = np.zeros((tumor_thresh.shape[1],tumor_thresh.shape[2]))
                sum_2D_proj = np.zeros((tumor_thresh.shape[1],tumor_thresh.shape[2]))
                orig_2D_proj = np.zeros((movie_tumor[tp].shape[1],  movie_tumor[tp].shape[2]))
                plot_df = single_macrophage_properties[single_macrophage_properties["time_frame"] == tp]
                plot_df = plot_df[plot_df["short_name"] == filename]
                for x in range(tumor_thresh.shape[1]):
                    for y in range(tumor_thresh.shape[2]):
                        mean_2D_proj[x,y] = np.mean(tumor_distances[:,x,y])
                        sum_2D_proj[x,y] = np.sum(tumor_thresh[:,x,y])
                        orig_2D_proj[x,y] = np.sum(movie_tumor[tp][:,x,y])

                fig, ax = plt.subplots(1,3, figsize=(45,15))
                ax[0].imshow(orig_2D_proj[:,:])
                ax[1].imshow(sum_2D_proj[:,:])
                ax[2].imshow(mean_2D_proj[:,:])
                ax[0].set_title("sum projection of tumor channel")
                ax[1].set_title("threshold image of tumor channel")
                ax[2].set_title("mean along z-axis of distance transform")
                
                print(single_macrophage_properties.head())
                for ind, row_plt in plot_df.iterrows():
                    #if row_plt["macrophage_volume"] < 100000:

                    x = int(row_plt['x_centroid'])
                    y = int(row_plt['y_centroid'])
                    z = int(row_plt['z_centroid'])
                    ax[2].plot(y, x, 'rx', markersize = 15)
                    ax[2].text(y,x,str(round(tumor_distances[z,x,y],2)))

                    ax[1].plot(y, x, 'rX', markersize = 15)
                    for col in plot_df.columns:
                        print(col)
                        print(row_plt[col])
                        distances_df.at[index_counter, col] = row_plt[col]
                        #distances_df.at[index_counter, "short_name"] = filename
                    distances_df.at[index_counter, "tumor_distance"] = tumor_distances[z,x,y]
                    print("Distance:")
                    print(tumor_distances[z,x,y])
                        
                        #distances_df.at[index_counter, "time_frame"] = tp
                        #distances_df.at[index_counter, "short_name"] = filename
                    index_counter += 1

                    
                time_stamp = "-" + str(tp).zfill(3)
                plt.savefig(parameters["output_folder"] + filename + "-tumor-distance" + time_stamp + ".pdf")
                plt.savefig(parameters["output_folder"] + filename + "-tumor-distance" + time_stamp + ".png")

                distances_df.to_csv(parameters["output_folder"] + "tumor_distances_" + experiment + ".csv")


    return distances_df
    
def get_tumor_volume(parameters, key_file, experiment = "all", proj_2D = True):
    """
    Function to extract the tumor volumes from all movies, plot them and save them in a dataframe.
    :param parameters:
    :param key_file:
    :param macrophage_properties:
    :param experiment:
    :param proj_2D:
    :return:
    """
    tumor_volume_df = pd.DataFrame()

    index_counter = 0

    for index, row in key_file.iterrows():

        filename = row["short_name"]
        print(filename)
        print(experiment)

        if pd.isna(filename):
            continue

        print("read info")

        if (experiment == "all" or experiment == filename):
            tumor_volume_df.at[index_counter, "short_name"] = filename
            tumor_volume_df.at[index_counter, "dpi"] = row["dpi"]
            tumor_volume_df.at[index_counter, "cancer_cells"] = row["cancer_cells"]
            time_interval_in_minutes = row["Ti (min"]

            file_path = parameters["data_folder"] + "03_Preprocessed_Data/02_3D/" + filename + '.tif'

            print("Row:")
            print(row)

            print("Loading data...")
            movie = np.array(io.imread(file_path))
            movie_tumor = movie[:, :, :, :, parameters["channel_tumor"]]

            del movie

            time_points = movie_tumor.shape[0]

            tumor_volumes = []
            tumor_volume_percentages = []

            for tp in range(time_points):
                tumor_blurred = gaussian_filter(movie_tumor[tp], sigma=parameters['sigma'])
                tumor_thresh, threshold = functions_common.thresholding_3D(parameters, tumor_blurred)
                tumor_volumes.append(get_volume(tumor_thresh))
                tumor_volume_percentages.append(get_volume_percentage(tumor_thresh))

            tumor_volume_df.at[index_counter, "tumor_volumes"] = tumor_volumes
            tumor_volume_df.at[index_counter, "tumor_volume_percentages"] = tumor_volume_percentages
            tumor_volume_df.at[index_counter, "time"] = range(0, len(tumor_volumes)*time_interval_in_minutes, time_interval_in_minutes)

            index_counter += 1


    tumor_volume_df.to_csv(parameters["output_folder"] + "tumor_volumes_" + experiment + ".csv")

    #tumor_volume_df_1dpi = tumor_volume_df[tumor_volume_df["dpi"==1]]
    #tumor_volume_df_5dpi = tumor_volume_df[tumor_volume_df["dpi"==5]]

    if experiment == 'all':
        # plot tumor volumes
        plt.figure()
        #sns.lineplot(data=tumor_volume_df_1dpi, x="time", y="tumor_volumes", hue="cancer_cells")
        sns.relplot(data=tumor_volume_df, x="time", y="tumor_volumes", hue="cancer_cells", col="dpi", col_wrap=2)

        plt.xlabel('time [min]', fontsize=14)
        plt.ylabel("Tumor volume (squ-pixel)", fontsize=14)
        plt.tick_params(labelsize=12)

        plt.savefig(parameters["output_folder"] + "tumor_volumes.png", format="png", bbox_inches="tight", dpi=150)
        plt.savefig(parameters["output_folder"] + "tumor_volumes.pdf", format="pdf", bbox_inches="tight", dpi=150)
        plt.close()

        # plot tumor volume percentages
        plt.figure()
        sns.relplot(data=tumor_volume_df, x="time", y="tumor_volume_percentages", hue="cancer_cells", col="dpi", col_wrap=2)

        plt.xlabel('time [min]', fontsize=14)
        plt.ylabel("Tumor volume (squ-pixel)", fontsize=14)
        plt.tick_params(labelsize=12)

        plt.savefig(parameters["output_folder"] + "tumor_volume_percentages.png", format="png", bbox_inches="tight", dpi=150)
        plt.savefig(parameters["output_folder"] + "tumor_volume_percentages.pdf", format="pdf", bbox_inches="tight", dpi=150)
        plt.close()

    return tumor_volume_df


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
