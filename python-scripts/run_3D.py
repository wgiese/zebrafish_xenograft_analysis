import functions_3D

parameters = functions_2D.read_parameters()
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_2D.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())
print(key_file.keys())

# Do Otsu Thresholding if specified in the parameters file
if parameters["otsu_thresholding"]:
    functions_2D.otsu_loop_over_key_file(parameters, key_file)

# Do analysis


