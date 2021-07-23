import functions_common
import functions_3D

parameters = functions_common.read_parameters()
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())
print(key_file.keys())


macrophage_properties_df = functions_3D.get_macrophage_properties(parameters, key_file, experiment ="200804_LBT070_1dpi_Pos001", npy_out = True, vtk_out = True)

print(macrophage_properties_df)



# Do Otsu Thresholding if specified in the parameters file
#if parameters["otsu_thresholding"]:
#    functions_3D.otsu_loop_over_key_file(parameters, key_file)




# Do analysis


