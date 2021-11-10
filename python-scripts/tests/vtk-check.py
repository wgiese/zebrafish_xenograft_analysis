from pyevtk.hl import imageToVTK
import numpy as np
import pandas as pd
import skimage.morphology
import skimage.draw

data_3d = pd.read_csv("test.csv", sep = ";")

print(data_3d)

#remove negative/ to small z coordinates

annotated_macrophages = np.zeros((np.max(data_3d['# time_point'])+1,1024,1024,256))

data_3d = data_3d[data_3d['Z'] > 0.0]

print(data_3d)

radius = 10

for frame in data_3d['# time_point'].unique():
    print("Frame: %s" % frame)
    data_3d_ = data_3d[data_3d['# time_point'] == frame]
    
    for ind, row in data_3d_.iterrows():
        
        print("Macrophage number: %s" % row['number'])

        x_m = row['X']
        y_m = row['Y']
        z_m = row['Z']
        
        rr, cc = skimage.draw.disk((x_m, y_m), radius, shape=(1024,1024))
        annotated_macrophages[frame,rr,cc,z_m] = 1
    
        for z in range(z_m - radius +1, z_m + radius):
            if (z > 0) and (z < 256):
                radius_z = np.sqrt(radius**2 - (z-z_m)**2)
                
                rr, cc = skimage.draw.disk((x_m, y_m), radius_z, shape=(1024,1024))
                annotated_macrophages[frame,rr,cc,z] = 1



        #ellip = skimage.draw.ellipsoid(2.0, 2.0, 2.0, spacing=(x_m, y_m, z_m), levelset=False)
        #print(ellip.shape)
    
    #for x in range(1024):
    #	for y in range(1024):
    #        for z in range(np.max(data_3d['# time_point'])+1):
    #            for ind, row in data_3d_.iterrows():
    #                #print("MAcrophage number: %s" % row['number'])
                    #data_3d_macrophage = data_3d_[data_3d_['number'] == number]
    #                x_m = row['X']
    #                y_m = row['Y']
    #                z_m = row['Z']
    #                dist2 = (x-x_m)**2 + (y-y_m)**2 + (z-z_m)**2
    #                if dist2 < radius*radius:
    #                    annotated_macrophages[frame,x,y,z] = 1

        #annotated_macrophages[frame,x+1,y,z] = 1
        #annotated_macrophages[frame,x-1,y,z] = 1
        #annotated_macrophages[frame,x+1,y+1,z] = 1
        #annotated_macrophages[frame,x+1,y-1,z] = 1
        #annotated_macrophages[frame,x+1,y+1,z+1] = 1
        #annotated_macrophages[frame,x+1,y+1,z-1] = 1
        #annotated_macrophages[frame,x,y+1,z+1] = 1
        #annotated_macrophages[frame,x,y+1,z-1] = 1

    #annotated_macrophages = skimage.morphology.binary_dilation(annotated_macrophages)

    time_stamp = str(frame).zfill(3)
    filename = "2d-3d-macrophage-" + time_stamp
   
    imageToVTK(filename, cellData = {"annotated_macrophage": annotated_macrophages[frame,:,:,:]})
    

#time_stamp = str(frame).zfill(3)
#filename = output_folder + filename_short + "-" + time_stamp
   
#imageToVTK(filename, cellData = {"tumor": im_tumor_smoothed, 
#                      "macrophages" : im_macrophages_smoothed, 
#                      "vessel": im_vessel_smoothed})
