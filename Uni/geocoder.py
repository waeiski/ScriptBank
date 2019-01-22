'''
This is a geocoding script based on Peter Gleeson's An A-Z of useful Python tricks.

USAGE
=====
    The script eats CSV files with addresses per line, geocodes them using OSM's nominatim
    and returns a pickled GeoDataframe in WGS84.
    Run it with the following command:

    python geocoder.py -i my_addresses.csv -o geocoded_addresses.pkl

ARGUMENTS
=========
    -i/--input: Path to the csv file with addresses
    -o/--output: Path to output pickled geodataframe

RETURNS
=======
    Pickle containing the address, latitude, longitude and shapely geometry

'''
from geopy import Nominatim
import pandas as pd
import geopandas as gpd
import argparse
from shapely.geometry import Point
import pycrs

# set up argument parser
ap = argparse.ArgumentParser()

# define arguments
ap.add_argument("-i", "--input", required=True,
                help="Path to the CSV file containing addresses and/or place "
                "names per line. No headers!")

ap.add_argument("-o","--output", required=True,
                help="Path to the output pickled geodataframe.")

# parse arguments
args = vars(ap.parse_args())

# read addresses in
df = pd.read_csv(args['input'], header=None, encoding='utf-8', names=['address'])
print('\n [INFO] - Addresses are read in')

# create list out of addresses
adr_df = gpd.GeoDataFrame()
adr_df['geometry'] = None
print('\n [INFO] - Empty GeoDataFrame created')

# geocode addresses
print('\n [INFO] - Geocoding with nominatim')
for i, row in df.iterrows():
    place = row['address']
    location = Nominatim(user_agent='myapp').geocode(place)
    lat = location.latitude
    lon = location.longitude
    adr_df.at[i, 'address'] = location.address
    adr_df.at[i, 'latitude'] = lat
    adr_df.at[i, 'longitude'] = lon
print('\n [INFO] - Geocoding done')

# Create geometries
# NOTE: Must be outside previous for loop to work
print('\n [INFO] - Creating geometries')
for i, row in adr_df.iterrows():
    lat = row['latitude']
    lon = row['longitude']
    adr_df.at[i, 'geometry'] = Point(lat, lon)
print('\n [INFO] - Geometries created')

# Defining CRS to WGS84
epsg_code = 4326
adr_df.crs = pycrs.parser.from_epsg_code(epsg_code).to_proj4()

# save output
print('\n [INFO] - Defined projection to WGS84. Saving output...')
adr_df.to_pickle(args['output'])
print('\n [INFO] - Done!')