import functions_common
import functions_3D
import pandas as pd

parameters = functions_common.read_parameters()
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())
print(key_file.keys())


#macrophage_properties_df = functions_3D.get_macrophage_properties(parameters, key_file, experiment ="200804_LBT070_1dpi_Pos001", npy_out = False, vtk_out = True)#, 2D_projection = True)
#macrophage_properties_df = functions_3D.get_macrophage_properties(parameters, key_file, experiment ="all", npy_out = False, vtk_out = False)#, 2D_projection = True)
#macrophage_properties_df = functions_3D.get_macrophage_properties(parameters, key_file, experiment ="200811_LBT123_1dpi_Pos005", npy_out = False, vtk_out = False)#, 2D_projection = True)

macrophage_properties_df = pd.read_csv("/home/wgiese/Data/zebrafish_xenograft/output_data/threshold_sauvola/macrophage_properties_200811_LBT123_1dpi_Pos005.csv", index_col=0)

print(macrophage_properties_df)

distances_df = functions_3D.get_tumor_macrophage_point_distances(parameters, key_file, macrophage_properties_df, experiment = "all")

print(distances_df)



# Do Otsu Thresholding if specified in the parameters file
#if parameters["otsu_thresholding"]:
#    functions_3D.otsu_loop_over_key_file(parameters, key_file)




# Do analysis


