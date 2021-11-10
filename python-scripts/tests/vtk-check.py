from pyevtk.hl import imageToVTK
import numpy as np
import pandas as pd
import skimage.morphology

data_3d = pd.read_csv("test.csv", sep = ";")

print(data_3d)

#remove negative/ to small z coordinates

annotated_macrophages = np.zeros((np.max(data_3d['# time_point'])+1,1024,1024,256))

data_3d = data_3d[data_3d['Z'] > 0.0]

print(data_3d)

for frame in data_3d['# time_point'].unique():
    print(frame)
    data_3d_ = data_3d[data_3d['# time_point'] == frame]
    for ind, row in data_3d_.iterrows():
        #data_3d_macrophage = data_3d_[data_3d_['number'] == number]
        x = row['X']
        y = row['Y']
        z = row['Z']
        annotated_macrophages[frame,x,y,z] = 1
        annotated_macrophages[frame,x+1,y,z] = 1
        annotated_macrophages[frame,x-1,y,z] = 1
        annotated_macrophages[frame,x+1,y+1,z] = 1
        annotated_macrophages[frame,x+1,y-1,z] = 1
        annotated_macrophages[frame,x+1,y+1,z+1] = 1
        annotated_macrophages[frame,x+1,y+1,z-1] = 1
        annotated_macrophages[frame,x,y+1,z+1] = 1
        annotated_macrophages[frame,x,y+1,z-1] = 1

    #annotated_macrophages = skimage.morphology.binary_dilation(annotated_macrophages)

    time_stamp = str(frame).zfill(3)
    filename = "2d-3d-macrophage-" + time_stamp
   
    imageToVTK(filename, cellData = {"annotated_macrophage": annotated_macrophages[frame,:,:,:]})
    

#time_stamp = str(frame).zfill(3)
#filename = output_folder + filename_short + "-" + time_stamp
   
#imageToVTK(filename, cellData = {"tumor": im_tumor_smoothed, 
#                      "macrophages" : im_macrophages_smoothed, 
#                      "vessel": im_vessel_smoothed})
