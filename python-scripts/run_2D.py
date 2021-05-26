import functions
import os
import pandas as pd

parameters = functions.read_parameters()
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())
print(key_file.keys())

# Do Otsu Thresholding if specified in the parameters file
if parameters["otsu_thresholding"]:
    functions.otsu_loop_over_key_file(parameters, key_file)

# Do analysis, if dataframe doesn't exist yet
file_path_df = parameters['output_folder'] + '05_Data_Analysis/01_2D/01_Dataframes/' + parameters['segmentation_method'] + '/analysis_dataframe.pkl'
if not os.path.isfile(file_path_df):
    print('Dataframe needs to be obtained from the data.')
    dataframe = pd.merge(key_file, functions.do_analysis_on_all_files(parameters, key_file))
else:
    print('Dataframe exists already and is loaded for plotting...')
    dataframe = pd.merge(key_file, pd.read_pickle(file_path_df))

# Plot results
functions.plot_all(parameters, dataframe)