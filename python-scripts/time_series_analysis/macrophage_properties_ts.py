import time, os, sys
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
sys.path.insert(0,"../")
import functions_common
import argparse

### read parameters

parser = argparse.ArgumentParser(description='Run cellpose macrophage segmentation.')
parser.add_argument('param', type=str, help='Path to the parameter file.')

args = parser.parse_args()
print("-------")
print("reading parameters from: ", args.param)
print("-------")

parameter_file  = args.param

parameters = functions_common.read_parameters(base_params = "../base/parameters.yml", local_params = parameter_file)
print("#"*5,"parameters","#"*5)
print(parameters)
key_file = functions_common.read_key_file(parameters)
print("#"*5,"key_file","#"*5)
print(key_file.head())

# plotting parameters
plt.rcParams.update({'font.size': 16})

# setup paths to data

macrophages_path_1dpi = parameters["macrophage_input_path_1dpi"]
summary_key_file_1dpi = pd.read_csv(macrophages_path_1dpi + 'analysis_summary.csv')

macrophages_path_5dpi = parameters["macrophage_input_path_5dpi"]
summary_key_file_5dpi = pd.read_csv(macrophages_path_5dpi + 'analysis_summary.csv')

start_time_points = parameters["start_time"]
end_time_points = parameters["end_time"]


# read macrophage properties from segmentation results 

def get_macrophage_properties_df(summary_key_file, macrophages_path):
    macrophage_props = pd.DataFrame()
    for index, row in summary_key_file.iterrows():
        single_movie_fn = row["short_name"] + ".csv"
        short_name = single_movie_fn.split(".")[0]
        path = macrophages_path + single_movie_fn
        if os.path.isfile(path):
            single_movie_props = pd.read_csv(macrophages_path + single_movie_fn, sep=";")
        else:
            continue
        # TODO: extract observed time only
        # note that key file entries t_end and t_start are assuming numbering of time frames from 1
        single_movie_props = single_movie_props[single_movie_props['time_point']<= row['t_end']-1]
        single_movie_props = single_movie_props[single_movie_props['time_point']>= row['t_start']-1]

        if "1dpi" in short_name:
            single_movie_props['fish_id'] = short_name.split("1dpi_")[0] + short_name.split("1dpi_")[1] 
        else:
            single_movie_props['fish_id'] = short_name.split("5dpi_")[0] + short_name.split("5dpi_")[1] 


        if macrophage_props.shape[0]>1:
            macrophage_props = pd.concat([macrophage_props, single_movie_props], ignore_index = True)
        else:
            macrophage_props =  single_movie_props.copy()

    return macrophage_props


def time_frame_to_min(macrophage_properties_,dpi):
    
    #dt_min = macrophage_properties_["dt_min"].iloc[0]
    macrophage_properties = macrophage_properties_.copy()
    # macrophage_properties["circularity"] = 4.0*macrophage_properties["Area"]/(np.pi*macrophage_properties["perimeter_px"]) 
   
    macrophage_properties["time_in_min"] = macrophage_properties["time_point"]*macrophage_properties_["dt_min"] + start_time_points[dpi]
    macrophage_properties["time_in_h"] = macrophage_properties["time_in_min"]/60.0
    macrophage_properties["dpi"] = dpi

    return macrophage_properties
    

macrophage_properties_1dpi = get_macrophage_properties_df(summary_key_file_1dpi, macrophages_path_1dpi)
macrophage_properties_1dpi = time_frame_to_min(macrophage_properties_1dpi, dpi='1dpi')
print(macrophage_properties_1dpi)

macrophage_properties_5dpi = get_macrophage_properties_df(summary_key_file_5dpi, macrophages_path_5dpi)
macrophage_properties_5dpi = time_frame_to_min(macrophage_properties_5dpi, dpi='5dpi')
print(macrophage_properties_5dpi) 

# TODO: save 1dpi and 5dpi with all time points
macrophage_properties_1dpi.to_csv("macrophage_props_1dpi.csv", index = False)

macrophage_properties_5dpi.to_csv("macrophage_props_5dpi.csv", index = False)

macrophage_properties = pd.concat([macrophage_properties_1dpi,macrophage_properties_5dpi], ignore_index = True)

# TODO: transform to hourly
start_1dpi = start_time_points['1dpi']
end_5dpi = end_time_points['5dpi']
obs_time_points = np.arange(start_1dpi,end_5dpi,60)
print("Target observation time points:")
print(obs_time_points)
print("Time points in the data set:")
print(macrophage_properties["time_in_min"].unique())



macrophage_properties = macrophage_properties[macrophage_properties["time_in_min"].isin(obs_time_points)]


print("Data types in macrophage properties data frame:")
print(macrophage_properties.dtypes)

new_types = {   'short_name' : 'object', 
                'fish_id' : 'object',
                'time_point' : 'int16',
                'number' : 'int16',
                'Area' : 'int16', 
                'Mean' : 'float16',
                'Min' : 'int16',
                'Max'  : 'int16',
                'X' : 'float16',
                'Y' : 'float16',
                'area_mum2' : 'float16',
                'dt_min' : 'float16',
                'time_in_min' : 'float16',
                'minor_axis_length_px' : 'float16',
                'minor_axis_length_mum' : 'float16',
                'major_axis_length_px' : 'float16',
                'major_axis_length_mum' : 'float16',
                'major_minor_ratio' : 'float16',
                'perimeter_px' : 'float16',
                'perimeter_mum' : 'float16',
                'eccentricity' : 'float16',
                'circularity' : 'float16',
                'cancer_cells' : 'object',
                't_start' : 'int16',
                't_end' : 'int16',
                'time_in_h' : 'float16',
                'dpi' : 'object' }

macrophage_properties = macrophage_properties.astype(new_types)
print("New data types in macrophage properties data frame (reduce space):")
print(macrophage_properties.dtypes)

macrophage_properties.to_csv("macrophage_props_1dpi_and_5dpi.csv", index = False)

# count macrophages per time point and image and store in data frame

macrophage_count = pd.DataFrame()
index = 0
for short_name in macrophage_properties["short_name"].unique():
    
    single_movie = macrophage_properties[macrophage_properties["short_name"]== short_name]
    
    for time_point in single_movie["time_point"].unique():
        single_frame = single_movie[single_movie["time_point"] == time_point]
        macrophage_count.at[index, "short_name"] = single_frame["short_name"].iloc[0]
        macrophage_count.at[index, "fish_id"] = single_frame["fish_id"].iloc[0]
        macrophage_count.at[index, "time_point"] = time_point
        macrophage_count.at[index, "macrophage_count"] = len(single_frame["number"])
        macrophage_count.at[index, "time_in_min"] = single_frame["time_in_min"].iloc[0]
        macrophage_count.at[index, "time_in_h"] = single_frame["time_in_h"].iloc[0]
        macrophage_count.at[index, "cancer_cells"] = single_frame["cancer_cells"].iloc[0]
        macrophage_count.at[index, "dpi"] = single_frame["dpi"].iloc[0]
        index +=1

macrophage_count.to_csv("macrophage_count_1dpi_and_5dpi.csv", index = False)

# prepare axis limits

start_1dpi = int(start_time_points['1dpi'])
end_1dpi = int(end_time_points['1dpi'])
delta_1dpi = end_1dpi - start_1dpi

start_5dpi = int(start_time_points['5dpi'])
end_5dpi = int(end_time_points['5dpi'])
delta_5dpi = end_5dpi - start_5dpi

ratio = delta_5dpi/delta_1dpi

max_count = macrophage_count["macrophage_count"].max()




if parameters["sample_full_hours"]:
    macrophage_properties = macrophage_properties[macrophage_properties["time_in_min"].isin(obs_time_points)]
    macrophage_count = macrophage_count[macrophage_count["time_in_min"].isin(obs_time_points)]

macrophage_count.to_csv("macrophage_count.csv", index =False)

# plot each fish_id as single line, color code by cancer cell line

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(30,10), gridspec_kw={'width_ratios': [delta_1dpi , delta_5dpi]})

sns.color_palette("Set3")
ax1.set_xlim(start_1dpi, end_1dpi)
ax2.set_xlim(start_5dpi, end_5dpi)

feature = "macrophage_count"

sns.scatterplot(x="time_in_min", y=feature, data=macrophage_count, ax=ax1, color="k", size =5, legend = False)
sns.scatterplot(x="time_in_min", y=feature, data=macrophage_count, ax=ax2, color="k", size =5,legend = False)

sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", data=macrophage_count, units="fish_id",
             linewidth = 1.0, palette = "colorblind",  estimator = None, ax=ax1, ci=95, color="c")
sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", data=macrophage_count, units="fish_id", 
             linewidth = 1.0, palette = "colorblind", estimator = None, ax=ax2, ci=95, color="c")

ax1.set_xlabel("time post injection [h]",  horizontalalignment='right', position=(0.75,25))
ax1.set_ylabel(feature)

ax2.set_xlabel("")
ax2.set_ylabel("")


ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_yticks([])

d = .02  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1 - d, 1 + d), (-2.0*d, +2.0*d), **kwargs)  # top-right diagonal

d_ = d/ratio 

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((- d_, + d_), (- 2.0*d, + 2.0*d), **kwargs)

xticks = np.arange(start_1dpi,end_1dpi,60)
xtick_labels = [i + 24  for i in range(len(xticks))] 
ax1.set_xticks(xticks)
ax1.set_xticklabels(xtick_labels)

xticks = np.arange(start_5dpi + 60,end_5dpi,60)
xtick_labels = [i + 5*24 + 1  for i in range(len(xticks))] 
ax2.set_xticks(xticks)
ax2.set_xticklabels(xtick_labels)

ax1.set_ylim(0,max_count + 10)
ax2.set_ylim(0,max_count + 10)

plt.savefig("macrophage_count_time_course.pdf")
plt.savefig("macrophage_count_time_course_all.png")




# plot comparison of selected conditions

filter_conditions = parameters["filter_conditions"]
macrophage_count_comparison = macrophage_count[macrophage_count["cancer_cells"].isin(filter_conditions)]
macrophage_comparison = macrophage_properties[macrophage_properties["cancer_cells"].isin(filter_conditions)]



fig, (ax1,ax2) = plt.subplots(1,2,figsize=(30,10), gridspec_kw={'width_ratios': [delta_1dpi , delta_5dpi]})

sns.color_palette("Set3")
ax1.set_xlim(start_1dpi, end_1dpi)
ax2.set_xlim(start_5dpi, end_5dpi)

feature = "macrophage_count"

#sns.lineplot(x="time_in_min", y=feature, data=count_macrophages, ax=ax1, ci=95, color="c")
#sns.lineplot(x="time_in_min", y=feature, data=count_macrophages, ax=ax2, ci=95, color="c")
sns.scatterplot(x="time_in_min", y=feature, data=macrophage_count_comparison , ax=ax1, 
                color="k", size =5, legend = False)
sns.scatterplot(x="time_in_min", y=feature, data=macrophage_count_comparison , ax=ax2, 
                color="k", size =5,legend = False)

sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", data=macrophage_count_comparison, 
             color="k", ax=ax1, ci=95)
sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", data=macrophage_count_comparison, 
             color="k", ax=ax2, ci=95)

sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", linewidth = 1.0, palette = "colorblind",
             data=macrophage_count_comparison, units="fish_id", estimator = None, ax=ax1, ci=95, color="c")
sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", linewidth = 1.0, palette = "colorblind", 
             data=macrophage_count_comparison, units="fish_id", estimator = None, ax=ax2, ci=95, color="c")

ax1.set_xlabel("time post injection [h]",  horizontalalignment='right', position=(0.75,25))
ax1.set_ylabel(feature)

ax2.set_xlabel("")
ax2.set_ylabel("")


ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_yticks([])

d = .02  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1 - d, 1 + d), (-2.0*d, +2.0*d), **kwargs)  # top-right diagonal

d_ = d/ratio 

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((- d_, + d_), (- 2.0*d, + 2.0*d), **kwargs)

xticks = np.arange(start_1dpi,end_1dpi,60)
xtick_labels = [i + 24  for i in range(len(xticks))] 
ax1.set_xticks(xticks)
ax1.set_xticklabels(xtick_labels)

xticks = np.arange(start_5dpi + 60,end_5dpi,60)
xtick_labels = [i + 5*24 + 1  for i in range(len(xticks))] 
ax2.set_xticks(xticks)
ax2.set_xticklabels(xtick_labels)

ax1.set_ylim(0,max_count + 10)
ax2.set_ylim(0,max_count + 10)

plt.savefig("macrophage_count_time_course_comparison.pdf")
plt.savefig("macrophage_count_time_course_comparison.png")




# plot comparison of selected feature 

feature = "Area"

max_feature = np.percentile(macrophage_properties[feature],95)
min_feature = np.percentile(macrophage_properties[feature],5)

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(30,10), gridspec_kw={'width_ratios': [delta_1dpi , delta_5dpi]})

sns.color_palette("Set3")
ax1.set_xlim(start_1dpi, end_1dpi)
ax2.set_xlim(start_5dpi, end_5dpi)

sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", data=macrophage_comparison, 
             color="k", ax=ax1, ci=90)
sns.lineplot(x="time_in_min", y=feature, hue="cancer_cells", data=macrophage_comparison, 
             color="k", ax=ax2, ci=90)

ax1.set_xlabel("time post injection [h]",  horizontalalignment='right', position=(0.75,25))
ax1.set_ylabel(feature)

ax2.set_xlabel("")
ax2.set_ylabel("")


ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_yticks([])

d = .02  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1 - d, 1 + d), (-2.0*d, +2.0*d), **kwargs)  # top-right diagonal

d_ = d/ratio 

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((- d_, + d_), (- 2.0*d, + 2.0*d), **kwargs)

xticks = np.arange(start_1dpi,end_1dpi,60)
xtick_labels = [i + 24  for i in range(len(xticks))] 
ax1.set_xticks(xticks)
ax1.set_xticklabels(xtick_labels)

xticks = np.arange(start_5dpi + 60,end_5dpi,60)
xtick_labels = [i + 5*24 + 1  for i in range(len(xticks))] 
ax2.set_xticks(xticks)
ax2.set_xticklabels(xtick_labels)

ax1.set_ylim(0.9*min_feature,1.1*max_feature)
ax2.set_ylim(0.9*min_feature,1.1*max_feature)

plt.savefig("macrophage_%s_time_course_comparison.pdf" % feature)
plt.savefig("macrophage_%s_time_course_comparison.png" % feature)

