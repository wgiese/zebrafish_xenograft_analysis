import pandas as pd
import yaml

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
    
    file_path = parameters['data_folder'] + "00_Extra/Details_imaging_data.xlsx"

    key_file = pd.read_excel(file_path)
        
    return key_file

parameters = read_parameters()
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())




