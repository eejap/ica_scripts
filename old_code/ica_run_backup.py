#!/usr/bin/env python
print("importing modules") 
import os
from sklearn.decomposition import FastICA, PCA 
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from osgeo import gdal
import rasterio
from scipy.interpolate import interp2d
from numpy.ma import masked_array
import glob
import h5py
from datetime import datetime, timedelta 
import geopandas as gpd
from scipy.stats import linregress
from matplotlib.dates import date2num, num2date
import pandas as pd
from shapely.geometry import Point, Polygon
from rasterio.transform import from_origin
import netCDF4 as nc
import argparse

# this is Andrew Watson's library of functions, see 
# https://github.com/Active-Tectonics-Leeds/interseismic_practical
import sys
import interseis_lib as lib

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process a single frame.')
parser.add_argument('--frame', type=str, help='Frame name')
args = parser.parse_args()

# Use the provided frame name or default to 'frame1' (or any default frame name)
frame = args.frame
#frame="028A_05385_191813"
print("loading files")

# load files for frames
cumh5_dir = "/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/cumh5"
mask_dir = "/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/mask"
EQA_dir = "/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/EQA.dem_par"
subs_poly_path = 'vU_merge_161023_noBabRas_WGS84_fillednosmooth_wmean_rad2_dist18_ltemin10_polygonised_dissolved.shp'

ncomponents = 4

output_nc_path = "/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/inelastic_components/{}_comp/{}_inelastic_component_{}.nc".format(ncomponents, frame, ncomponents)
output_tif_path = "/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/inelastic_components/{}_comp/{}_inelastic_component_{}.tif".format(ncomponents, frame, ncomponents)

frames_data = []
EQA_par_pattern = os.path.join(EQA_dir,f"{frame}_GEOCml*GACOSmask_EQA.dem_par")
EQA_par_file = glob.glob(EQA_par_pattern)
cumh5_pattern = os.path.join(cumh5_dir,f"{frame}_GEOCml*GACOSmask_cum.h5")
cumh5_file = glob.glob(cumh5_pattern)
mask_pattern = os.path.join(mask_dir,f"{frame}_GEOCml*GACOSmask_coh_03_mask.geo.tif") 
mask_file = glob.glob(mask_pattern)

with h5py.File(cumh5_file[0], 'r') as file:
    imdates = file['imdates']
    imdates = imdates[:]
    vel = file['vel']
    vel = vel[:]
    cum = file['cum']
    cum = cum[:]
    
dates=[]
for date_value in imdates:
    date_string = str(date_value) # Convert int32 to string
    year = int(date_string[:4])
    month = int(date_string[4:6])
    day = int(date_string[6:])
    real_date = datetime(year, month, day)
    dates.append(real_date)
    width = int(lib.get_par(EQA_par_file[0],'width'))
    length = int(lib.get_par(EQA_par_file[0],'nlines'))
    
# get corner positions
corner_lat = float(lib.get_par(EQA_par_file[0], 'corner_lat'))
corner_lon = float(lib.get_par(EQA_par_file[0], 'corner_lon'))

# get post spacing (distance between velocity measurements)
post_lat = float(lib.get_par(EQA_par_file[0],'post_lat'))
post_lon = float(lib.get_par(EQA_par_file[0],'post_lon'))

# calculate grid spacings
lat = corner_lat + post_lat*np.arange(1,length+1) - post_lat/2
lon = corner_lon + post_lon*np.arange(1,width+1) - post_lon/2
frames_data.append({
    'frame': frame,
    'EQA_file': EQA_par_file,
    'cumh5_file': cumh5_file,
    'mask_file': mask_file,
    'imdates': imdates,
    'vel': vel,
    'cum': cum,
    'dates': dates,
    'width': width,
    'length': length,
    'corner_lat': corner_lat,
    'corner_lon': corner_lon,
    'post_lat': post_lat,
    'post_lon': post_lon,
    'lat': lat,
    'lon': lon
        })

# Create GeoDataFrame
frames_gdf = gpd.GeoDataFrame(frames_data, columns=['frame', 'EQA_file', 'cumh5_file', 'mask_file', 'imdates', 'vel', 'cum', 'dates', 'width', 'length', 'corner_lat', 'corner_lon', 'post_lat', 'post_lon','lat', 'lon']) 

print("ica prep")
# cum shape is (t, lat, lon) we want to make an array of (pixels, time) we want to reshape cum (202, 268, 327) 
# into (202,(268*327))
frames_gdf["cum_"] = ""
for index, row in frames_gdf.iterrows():
    frame = row['frame']
    cum_data = row['cum']
    # Check if 'cum' data is not empty
    if not np.isnan(cum_data).all():
        cum_shape = cum_data.shape
        n_pix = cum_shape[1] * cum_shape[2]
        cum_ = np.reshape(cum_data, (cum_shape[0], n_pix))
        frames_gdf.at[index, 'cum_'] = cum_
        
# we have cum of shape 202, 268, 327 and mask tif. We want to mask cum with mask tif i.e. make pixels in cum nan 
# where mask_tif = 0. NaNs are where no data e.g. outside of frame, low coh...
frames_gdf["cum_masked"] = ""
for index, row in frames_gdf.iterrows():
    frame = row['frame']
    cum_data = row['cum_']
    mask = row['mask_file']
    with rasterio.open(mask[0]) as tif:
        # Read the raster data
        mask_tif = tif.read(1)
        # Reshape the mask to 1D
        mask_1d = mask_tif.flatten()
        # Tile the first row along the rows to match the shape of asc_cum_
        mask_2d = np.tile(mask_1d, (cum_data.shape[0],1))
        
        # Apply the mask to every element in the 3D array
        cum_masked = cum_data * mask_2d
        frames_gdf.at[index, 'cum_masked'] = cum_masked
        
# Find columns (pixels) containing zeros
frames_gdf["cum_no_nans_zeros"] = ""
frames_gdf["non_zero_ind"] = ""
frames_gdf["non_nan_ind"] = ""

for index, row in frames_gdf.iterrows():
    frame = row['frame']
    cum_masked = row['cum_masked']
    zero_pixels = np.all(cum_masked == 0, axis=0)
    # Create a new data array without columns (pixels) containing all zeros
    cum_no_zeros = cum_masked[:, ~zero_pixels]
    # Print how many NaNs there are
    nan_indices = np.argwhere(np.isnan(cum_no_zeros))
    # find columns containing NaNs and zeroes
    nan_pixels = np.any(np.isnan(cum_no_zeros), axis=0)
    # create a new data array without nan columns (pixels)
    cum_no_zeros_no_nans = cum_no_zeros[:, ~nan_pixels]
    zero_ind = np.argwhere(zero_pixels).flatten()
    non_zero_ind = np.argwhere(~zero_pixels).flatten()
    nans = np.any(np.isnan(cum_masked), axis=0)
    nan_ind = np.argwhere(nans).flatten()
    non_nan_ind = np.argwhere(~nans).flatten()
    frames_gdf.at[index, 'cum_no_nans_zeros'] = cum_no_zeros_no_nans
    frames_gdf.at[index, 'non_zero_ind'] = non_zero_ind
    frames_gdf.at[index, 'non_nan_ind'] = non_nan_ind

print("run ica")
frames_gdf["S_ft"] = "" 
frames_gdf["restored_signals"] = ""

# attempt ICA
for index, row in frames_gdf.iterrows():
    frame = row['frame']
    data = row['cum_no_nans_zeros']
    # set up the transformer (ica). In MATLAB you do the whitening first then the transforming, here you do it in 
    # one.
    ica = FastICA(n_components=ncomponents, whiten="unit-variance")
    
    # fit the transformeter to the data array
    try:
        S_ft = ica.fit_transform(data) # fit model and recover signals
        S_t = ica.transform(data) # recover sources from x using unmixing matrix
    except ConvergenceWarning as e:
        print("ICA did not converge for frame:", frame)
        raise e
    
    ## S_ft and S_t results should be identical as ica.transform uses mixing matrix calculated by 
    ## ica.fit_transform
    frames_gdf.at[index, 'S_ft'] = S_ft
    # Take each signal and restore with outer product
    restored_signals_outer = []
    for j in range(ncomponents):
        S_j = np.copy(S_ft)
        signal = S_j[:,j]
        mixing = ica.mixing_[:,j]
        restored_signal_j = np.outer(signal, mixing)
        
        # Append the restored signal to the list
        restored_signals_outer.append(restored_signal_j)
        
    frames_gdf.at[index, 'restored_signals'] = restored_signals_outer
    
# Find the common indices
for index, row in frames_gdf.iterrows():
    frame = row['frame']
    non_nan_ind = row['non_nan_ind']
    non_zero_ind = row['non_zero_ind']
    restored_signals_outer = row['restored_signals']
    cum = row['cum']
    lon_plot = row['lon']
    lat_plot = row['lat']
    common_indices = np.intersect1d(non_nan_ind, non_zero_ind)
    
    # do some reshaping
    for restored_signal in restored_signals_outer:
        # Create a new matrix with NaNs
        cum_with_nans = np.full((cum.shape[2] * cum.shape[1],), np.nan)
         # Assign values from the restored signal to non-NaN positions
        cum_with_nans[common_indices] = restored_signal[-1]
        cum_with_nans_reshaped = cum_with_nans.reshape((cum.shape[1], cum.shape[2]))
print("checking for linearity")

# keep only subsiding polygon pixels from cum_nozeros_no_nans convert polygons shape into a geopandas dataframe
frames_gdf["subsiding_restored_signals"] = ""
frames_gdf["restored_signals_3d"] = ""
gdf_polygons = gpd.read_file(subs_poly_path)
for index, row in frames_gdf.iterrows():
    frame = row['frame']
    non_nan_ind = row['non_nan_ind']
    non_zero_ind = row['non_zero_ind']
    restored_signals_outer = row['restored_signals']
    cum = row['cum']
    lon_plot = row['lon']
    lat_plot = row['lat']
    common_indices = np.intersect1d(non_nan_ind, non_zero_ind)
    cum_with_nans = np.full((cum.shape[0], cum.shape[1], cum.shape[2]), np.nan)
    
    # Create meshgrid of lon and lat
    lon, lat = np.meshgrid(lon_plot, lat_plot)
    
    # Flatten lon and lat
    lon_1d = lon.flatten()
    lat_1d = lat.flatten()
    
    # Create a GeoDataFrame for the flattened lon and lat
    geometry = [Point(lon, lat) for lon, lat in zip(lon_1d, lat_1d)]
    gdf_points = gpd.GeoDataFrame(geometry=geometry, columns=['geometry'])
    gdf_points.crs = "EPSG:4326"
    gdf_points['Latitude'] = lat_1d
    gdf_points['Longitude'] = lon_1d
    
    # Perform Spatial Join
    joined = gpd.sjoin(gdf_points, gdf_polygons, predicate='within')
    extracted_coordinates = joined[['Latitude', 'Longitude']]
    extracted_indices = joined.index
    
    # Create a boolean mask for lon and lat arrays
    mask = np.isin(np.arange(lon.size), extracted_indices)
    subsiding_restored = []
    three_d_list = []
    
    for restored_signal in restored_signals_outer:
        # Convert common_indices to 3D indices i.e. flat indices into a tuple for 3d array
        indices_3d = np.unravel_index(common_indices, cum.shape[1:])
        # Create a copy of cum_with_nans to work with for each restored_signal
        cum_with_nans_copy = cum_with_nans.copy()
        
        # need to make sure a different time series assigned to timestep
        for j in range(restored_signal.shape[0]):
            
            # Assign values from the restored signal to non-NaN positions
            cum_with_nans_copy[j, indices_3d[0], indices_3d[1]] = restored_signal[j,:]
        
        # reshape cum_with_nans into time x pixels
        cum_with_nans_pix = cum_with_nans_copy.reshape(cum.shape[0], cum.shape[1] * cum.shape[2])
        
        #reshape cum_with_nans_pix into 3d to save signals
        restored_signal_3d = cum_with_nans_pix.reshape(cum.shape[0], cum.shape[1], cum.shape[2])
        three_d_list.append(restored_signal_3d)
        
        # mask cum_with_nans_pix with extracted indices
        masked_cum_with_nans_pix = cum_with_nans_pix * np.where(mask == 0, np.nan, 1)
        subsiding_restored.append(masked_cum_with_nans_pix)
        
    frames_gdf.at[index, 'subsiding_restored_signals'] = subsiding_restored
    frames_gdf.at[index, "restored_signals_3d"] = three_d_list

# remove nans and zeros from each subsiding restored signal
frames_gdf["subsiding_no_nans"] = ""
for index, row in frames_gdf.iterrows():
    frame = row['frame']
    subsiding_restored = row['subsiding_restored_signals']
    no_nans_store = []
    for restored_signal in subsiding_restored:
        # Print how many NaNs there are
        nan_indices = np.argwhere(np.isnan(restored_signal))
        # find columns containing NaNs and zeroes
        nan_pixels = np.any(np.isnan(restored_signal), axis=0)
        # create a new data array without nan columns (pixels)
        subsiding_no_nans = restored_signal[:, ~nan_pixels]
        no_nans_store.append(subsiding_no_nans)
    frames_gdf.at[index, 'subsiding_no_nans'] = no_nans_store

# find most linear!!
frames_gdf["inelastic_restored_signal_3d"] = ""
for index, row in frames_gdf.iterrows():
    frame = row['frame']
    restored_signals = row['subsiding_no_nans']
    dates = row['dates']
    restored_signals_3d = row['restored_signals_3d']
    restored_signals_full = row['restored_signals']
    corner_lat = row['corner_lat']
    corner_lon = row['corner_lon']
    post_lon = row['post_lon']
    post_lat = row['post_lat']
    width = row['width']
    height = row['length']
    lat = row['lat']
    lon = row['lon']
    nc_data = row['imdates']
    
 # Only proceed if there are enough coherent pixels
    if restored_signals[0].shape[1] >= 20:
        # List to store R-squared values for each trend
        mean_gradient = []
        median_r_squared = []
        median_dominant_period = []
        mean_second_derivative = []

        # Convert dates to numerical values
        num_dates = date2num(dates)

        for signal in restored_signals:
            r_squared_values = []
            gradient = []
            second_deriv = []
            dominant_periods = []

            # Loop through each trend at each pixel
            for i in range(signal.shape[1]):
                # Perform linear regression and calculate R-squared
                slope, _, r_value, _, _ = linregress(num_dates, signal[:, i])

                # Store the gradient at each pixel
                gradient.append(slope)

                # Append R-squared value to the list
                r_squared_values.append(r_value ** 2)

                # Calculate the second derivative using numpy.gradient at each pixel
                first_derivative = np.gradient(signal[:, i], num_dates)
                second_deriv.append(np.gradient(first_derivative, num_dates))

                # Perform FFT and find frequencies
                fft_result = np.fft.fft(signal[:, i])

                # Calculate the corresponding frequency using FFT frequencies
                freqs = np.fft.fftfreq(len(signal[:, i]))
                dominant_frequency = np.abs(freqs[np.argmax(np.abs(fft_result))])
                dominant_period = 1 / dominant_frequency
                dominant_periods.append(dominant_period)

            # take mean of gradients
            mean_gradient.append(np.mean(gradient))
            # take median of r_squared per IC
            median_r_squared.append(np.median(r_squared_values))
            # find median period of signal
            median_dominant_period.append(np.median(dominant_periods))
            # take mean of second derivative per IC
            mean_second_derivative.append(np.mean(second_deriv))
    
        print('median_r_squared:', median_r_squared)
        print('mean_gradient:', mean_gradient)
        print('median dominant period:', median_dominant_period, 'timesteps - do not use this data yet')
        print('mean_second_derivative:', mean_second_derivative, 'mm/yr^2')
    
        # If mean gradient is negative, find the index of the trend with the maximum R-squared value
        negative_indices = [i for i, val in enumerate(mean_gradient) if val < 0]
        if negative_indices:
            median_r_squared_negative_gradients = [median_r_squared[i] for i in negative_indices]
            max_r_squared_negative_grad = np.max(median_r_squared_negative_gradients)
            max_index_in_median_r_squared = median_r_squared.index(max_r_squared_negative_grad)

            # Choose the signal
            inelastic_signal = restored_signals_3d[max_index_in_median_r_squared]
            inelastic_signal_subsiding = restored_signals[max_index_in_median_r_squared]

            # Save as NetCDF
            output_nc_dir = f"/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/inelastic_components/{ncomponents}_comp/"
            if not os.path.exists(output_nc_dir):
                os.makedirs(output_nc_dir)

            output_nc_path = os.path.join(output_nc_dir, f"{frame}_inelastic_component_{ncomponents}.nc")
            
            print('Saving inelastic nc')
            with nc.Dataset(output_nc_path, 'w') as file:
                # Create dimensions
                file.createDimension('dates', inelastic_signal.shape[0])
                file.createDimension('latitude', inelastic_signal.shape[1])
                file.createDimension('longitude', inelastic_signal.shape[2])

                # Create variables
                time_var = file.createVariable('dates', 'f4', ('dates',))
                lat_var = file.createVariable('latitude', 'f4', ('latitude',))
                lon_var = file.createVariable('longitude', 'f4', ('longitude',))
                data_var = file.createVariable('data', 'f4', ('dates', 'latitude', 'longitude'))
                period_var = file.createVariable('period', 'f4')
                second_derivative_var = file.createVariable('second_derivative', 'f4')

                # Add data to variables
                time_var[:] = nc_data
                lat_var[:] = lat
                lon_var[:] = lon
                data_var[:] = inelastic_signal
                period_var[:] = median_dominant_period[max_index_in_median_r_squared]
                second_derivative_var[:] = mean_second_derivative[max_index_in_median_r_squared]
            
            print('Saving inelastic geotiff')
            # Save as GeoTIFF
            output_tif_dir = f"/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/inelastic_components/{ncomponents}_comp/"
            if not os.path.exists(output_tif_dir):
                os.makedirs(output_tif_dir)

            output_tif_path = os.path.join(output_tif_dir, f"{frame}_inelastic_component_{ncomponents}.tif")
            flipped_data = np.flipud(inelastic_signal[-1, :, :])  # Assuming this is the correct data to flip
            # Create a transformation for the GeoTIFF
            transform = from_origin(corner_lon, corner_lat, post_lon, post_lat)

            with rasterio.open(output_tif_path, 'w', driver='GTiff', height=height, width=width, count=1,
                               dtype='float32', crs='EPSG:4326', transform=transform) as dst:
                # Write the data to the GeoTIFF
                dst.write(flipped_data, 1)

        # Save other signals as NetCDF
        output_nc_dir = f"/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/other_components/{ncomponents}_comp/{frame}"
        if not os.path.exists(output_nc_dir):
            os.makedirs(output_nc_dir)

        output_nc_path = os.path.join(output_nc_dir, f"other_component.nc")

        print('Saving other component nc file')
        with nc.Dataset(output_nc_path, 'w') as file:
            # Create dimensions
            file.createDimension('dates', inelastic_signal.shape[0])
            file.createDimension('latitude', inelastic_signal.shape[1])
            file.createDimension('longitude', inelastic_signal.shape[2])

            # Create variables
            time_var = file.createVariable('dates', 'f4', ('dates',))
            lat_var = file.createVariable('latitude', 'f4', ('latitude'))
            lon_var = file.createVariable('longitude', 'f4', ('longitude'))
            
            # Create variables for each component
            data_vars = []
            period_vars = []
            second_derivative_vars = []

            for comp_index in range(ncomponents - 1):
                data_var = file.createVariable(f'data_{comp_index}', 'f4', ('dates', 'latitude', 'longitude'))
                period_var = file.createVariable(f'period_{comp_index}', 'f4')
                second_derivative_var = file.createVariable(f'second_derivative_{comp_index}', 'f4')
    
                data_vars.append(data_var)
                period_vars.append(period_var)
                second_derivative_vars.append(second_derivative_var)

            # Add data, period, and second derivative to variables
            comp_index=0
            for i, signal_3d in enumerate(restored_signals_3d):
                if i != max_index_in_median_r_squared:
                    data_vars[comp_index][:] = signal_3d
                    period_vars[comp_index][:] = median_dominant_period[i]
                    second_derivative_vars[comp_index][:] = mean_second_derivative[i]
                    comp_index += 1

            time_var[:] = nc_data
            lat_var[:] = lat
            lon_var[:] = lon

            print('Saving other component geotiffs')
            # Save as GeoTIFF
            output_tif_dir = f"/gws/nopw/j04/nceo_geohazards_vol1/projects/COMET/eejap002/ica_data/all_iran/other_components/{ncomponents}_comp/{frame}"
            # Check if the directory exists
            if os.path.exists(output_tif_dir):
                # Iterate over the files in the directory
                for filename in os.listdir(output_tif_dir):
                # Check if the file is a GeoTIFF
                    if filename.endswith(".tif"):
                        # Construct the full path to the file
                        file_path = os.path.join(output_tif_dir, filename)
                        # Delete the file
                        os.remove(file_path)
            else:
                print("Directory does not exist.")
                
            for i, signal_3d in enumerate(restored_signals_3d):
                if i != max_index_in_median_r_squared:
                    output_tif_path = os.path.join(output_tif_dir, f"other_component_{i}.tif")
                    flipped_final_disp_other = np.flipud(signal_3d[-1, :, :])  # Assuming this is the correct data to flip

                    with rasterio.open(output_tif_path, 'w', driver='GTiff', height=height, width=width, count=1,
                                       dtype='float32', crs='EPSG:4326', transform=transform) as dst:
                        # Write the data to the GeoTIFF
                        dst.write(flipped_final_disp_other, 1)
