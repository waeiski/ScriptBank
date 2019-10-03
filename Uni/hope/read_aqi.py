# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 16:38:00 2019

@author: tuomvais
"""


import boto3
import botocore
import os

# replace with your bucket name
BUCKET_NAME = 'eu-central-1' 

OBJECT_NAME = 'Enfuser2/Finland/pks'
FILE_NAME = 'allPollutants_2019_09_11T15.zip'

# replace with your object key
file = 'Enfuser2/Finland/pks/allPollutants_2019_09_11T15.zip'

s3 = boto3.resource('s3')

#s3 = boto3.client('s3')

s3.download_file(BUCKET_NAME, OBJECT_NAME, FILE_NAME)

try:
    s3.Bucket(BUCKET_NAME).download_file(file, os.path.basename(file))
except botocore.exceptions.ClientError as e:
    if e.response['Error']['Code'] == "404":
        print("The object does not exist.")
    else:
        raise

###############################################################################
import rioxarray
import xarray

# import AQI mesh as xarray
data = xarray.open_dataset('allPollutants_2019-09-11T15.nc')

# Fetch AQI as xarray, AQI.data has shape (time, lat, lon)
AQI = data['AQI']

# save AQI to raster
AQI.rio.to_raster('allPollutants_raster_test.tif')


###############################################################################
from netCDF4 import Dataset

# read mesh file in
file = Dataset(r'C:\LocalData\tuomvais\Lataukset\ENFUSER_example\allPollutants_2019-09-11T15.nc', 'r+', format='NETCDF4')

# retrieve AQI as numpy masked array (returns in channels first format)
aqi = file.variables['AQI'][:]

# retrieve lat lon
lat = file.variables['latitude'][:]
lon = file.variables['longitude'][:]

