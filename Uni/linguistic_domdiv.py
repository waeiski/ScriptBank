# -*- coding: utf-8 -*-

"""
This script calculates different dominance and diversity values to a grid based
on topics of social media posts.

Usage:
    Execute the script from the command line using the following command:

    python3 topic_linguistic.py -g grid.gpkg -p posts.gpkg -t 'topic' -ge 'geometry' -o output.shp

Arguments:
    -g/--grid: Path to the file with grid (square, hexagon etc.)
    -p/--points: Path to the file with social media posts with topics.
    -t/--topic: Name of column containing topics as strings.
    -ge/--geometry: Column name for geometry of spatial data.
    -o/--output: Path to output shapefile or geopackage.

Returns:
    Pickle containing topic dominance, diversity per grid cell.
"""

import argparse
import pandas as pd
import geopandas as gpd
import numpy as np
from collections import Counter
import skbio.diversity.alpha as sk

# Set up the argument parser
ap = argparse.ArgumentParser()

# Define arguments
ap.add_argument("-g", "--grid", required=True,
                help="Path to the grid file. Must contain 'id' column with "
                "unique values for each polygon in grid.")

ap.add_argument("-p", "--points", required=True,
                help="Path to file containing the point geometries. File "
                "should be geopackage, shapefile or similar")

ap.add_argument("-ft", "--fthresh", required=False, type=float,
                help="fastText threshold for including data into the plot. "
                     "The value must be in range [0..1].")

ap.add_argument("-l", "--language", required=True,
                help="Column name for language. Must be strings.")

ap.add_argument("-ge", "--geometry", required=True,
                help="Column name for geometry data, e.g. 'geometry'.")

ap.add_argument("-o", "--output", required=True,
                help="Path to output pickle")

# Parse arguments
args = vars(ap.parse_args())

# Read in the spatial files
print('[INFO] - Loading grid file...')
grid = gpd.read_file(args['grid'])
print('[INFO] - Grid file loaded')
print('[INFO] - Loading points file...')
points = gpd.read_file(args['points'])
print('[INFO] - Points loaded')

# Check if the layers have the same coordinate system
print('[INFO] - Checking coordinate systems for spatial joining')
if grid.crs == points.crs:
    print('[INFO] - Coordinate systems match, continuing...')
else:
    print('[INFO] - Coordinate systems do not match! Make sure the data are in'
          ' the same CRS and try again.')
    print('[INFO] - Exiting....')
    exit

# Check if a confidence thresholds has been provided for fastText
if args['fthresh']:

    # Assign threshold to variable
    ft = args['fthresh']

# If thresholds have been defined, drop the predictions below the threshold
if args['fthresh']:

    # Filter posts based on fastText prediction confidence
    points = points.loc[points['probability'].apply(lambda x: float(x) >= ft)]

# Simplifying points geodataframe
points = points[[args['language'], args['geometry']]]

# Spatial join
print('[INFO] - Spatially joining points to grid')
joined = gpd.sjoin(grid, points, how='inner', op='contains')

# Dissolving
print('[INFO] - Dissolving results')
dissolved = joined.dissolve(by='id', aggfunc=lambda x: list(x))

# Getting output length
dis_len = len(dissolved)

# Counting language dominance, menhinick diversity and simpson index
print('[INFO] - Calculating variables..')
for i, row in dissolved.iterrows():
    print("[INFO] - Calculating grid cell {}/{}...".format(i, dis_len))
    lang_counts = list(Counter(row[args['language']]).values()) # occurence counts
    lang_counts = np.asarray(lang_counts) # cast as numpy array for skbio
    dissolved.at[i, 'dominance'] = sk.dominance(lang_counts)
    dissolved.at[i, 'menhinick'] = sk.menhinick(lang_counts)
    dissolved.at[i, 'simpson'] = sk.simpson(lang_counts)
    dissolved.at[i, 'berger'] = sk.berger_parker_d(lang_counts)
    dissolved.at[i, 'singles'] = sk.singles(lang_counts)
    dissolved.at[i, 'shannon'] = np.exp(sk.shannon(lang_counts, base=np.e))
    dissolved.at[i, 'unique'] = sk.observed_otus(lang_counts)

# Select columns for output
cols = ['geometry', 'dominance', 'menhinick', 'simpson', 'berger', 'singles',
        'shannon', 'unique']
output = dissolved[cols]

# Save the output to pickle
print('[INFO] - Saving to shapefile')
output.to_file(args['output'], encoding='utf-8')

# Print status
print("[INFO] - ... Done.")
