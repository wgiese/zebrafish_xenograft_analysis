import functions

parameters = functions.read_parameters()
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())
