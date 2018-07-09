# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 09:38:12 2018

DESCRIPTION
===========
This script reads in GeoPackages with lemmatized texts from Vantaa great
regions and creates wordclouds corresponding to the regions from columns
containing lemmatized texts. The script has been written to extract texts
specific to Yleiskaava 2020 -maptionnaire questions. The script user (YOU)
has to decide which question to use.

USAGE
=====
Check that filepaths are valid for your system.
Check that columns are valid for your files.
Tinker with visualization options if you want different outputs.
Run the script line by line using your preferences. Do not run the whole script
in one go.

REQUIREMENTS
============
Python 3.5+
GeoPandas
Pandas
WordCloud
Imageio

@author: Tuomas Väisänen, Vantaan kaupunki 2018
"""

import geopandas as gpd
import pandas as pd
from wordcloud import WordCloud
from collections import defaultdict
from imageio import imread

# Initializing the filepaths
av = r'C:\GIS\Vantaa\points\Lems\Lems_Aviapolis.gpkg'
hk = r'C:\GIS\Vantaa\points\Lems\Lems_Hakunila.gpkg'
kv = r'C:\GIS\Vantaa\points\Lems\Lems_Kivisto.gpkg'
kk = r'C:\GIS\Vantaa\points\Lems\Lems_Koivukyla.gpkg'
ko = r'C:\GIS\Vantaa\points\Lems\Lems_Korso.gpkg'
my = r'C:\GIS\Vantaa\points\Lems\Lems_Myyrmaki.gpkg'
tk = r'C:\GIS\Vantaa\points\Lems\Lems_Tikkurila.gpkg'

vp = r'C:\GIS\Vantaa\points\Virkistys_Lem.gpkg'

# Reading the files in as GeoDataFrames
aviapolis = gpd.read_file(av)
hakunila = gpd.read_file(hk)
kivisto = gpd.read_file(kv)
koivukyla = gpd.read_file(kk)
korso = gpd.read_file(ko)
myyrmaki = gpd.read_file(my)
tikkurila = gpd.read_file(tk)

# For Virkistys-themed answers only, do not use otherwise
virkistys =gpd.read_file(vp)

# Selecting the wanted rows from GeoDataFrames, select ONLY ONE group of the 
# following.
aviapolis = aviapolis.loc[aviapolis['buttonID'] == 'Omaleimainen alue']
hakunila = hakunila.loc[hakunila['buttonID'] == 'Omaleimainen alue']
kivisto = kivisto.loc[kivisto['buttonID'] == 'Omaleimainen alue']
koivukyla = koivukyla.loc[koivukyla['buttonID'] == 'Omaleimainen alue']
korso = korso.loc[korso['buttonID'] == 'Omaleimainen alue']
myyrmaki = myyrmaki.loc[myyrmaki['buttonID'] == 'Omaleimainen alue']
tikkurila = tikkurila.loc[tikkurila['buttonID'] == 'Omaleimainen alue']

aviapolis = aviapolis.loc[aviapolis['buttonID'] == 'Avo: alueen muutostoive']
hakunila = hakunila.loc[hakunila['buttonID'] == 'Avo: alueen muutostoive']
kivisto = kivisto.loc[kivisto['buttonID'] == 'Avo: alueen muutostoive']
koivukyla = koivukyla.loc[koivukyla['buttonID'] == 'Avo: alueen muutostoive']
korso = korso.loc[korso['buttonID'] == 'Avo: alueen muutostoive']
myyrmaki = myyrmaki.loc[myyrmaki['buttonID'] == 'Avo: alueen muutostoive']
tikkurila = tikkurila.loc[tikkurila['buttonID'] == 'Avo: alueen muutostoive']

aviapolis = aviapolis.loc[aviapolis['buttonID'] == 'Muu idea: avoin palaute kenttä']
hakunila = hakunila.loc[hakunila['buttonID'] == 'Muu idea: avoin palaute kenttä']
kivisto = kivisto.loc[kivisto['buttonID'] == 'Muu idea: avoin palaute kenttä']
koivukyla = koivukyla.loc[koivukyla['buttonID'] == 'Muu idea: avoin palaute kenttä']
korso = korso.loc[korso['buttonID'] == 'Muu idea: avoin palaute kenttä']
myyrmaki = myyrmaki.loc[myyrmaki['buttonID'] == 'Muu idea: avoin palaute kenttä']
tikkurila = tikkurila.loc[tikkurila['buttonID'] == 'Muu idea: avoin palaute kenttä']

# These are for virkistys
aviapolis = virkistys.loc[virkistys['nimi'] == 'Aviapolis']
hakunila = virkistys.loc[virkistys['nimi'] == 'Hakunila']
kivisto = virkistys.loc[virkistys['nimi'] == 'Kivistö']
koivukyla = virkistys.loc[virkistys['nimi'] == 'Koivukylä']
korso = virkistys.loc[virkistys['nimi'] == 'Korso']
myyrmaki = virkistys.loc[virkistys['nimi'] == 'Myyrmäki']
tikkurila = virkistys.loc[virkistys['nimi'] == 'Tikkurila']

# Reading in Vantaa great region masks
aviMask = imread(r'C:\GIS\Vantaa\wcmask\Maskit\Aviapolis-maski.png')
hakMask = imread(r'C:\GIS\Vantaa\wcmask\Maskit\Hakunila-maski.png')
kivMask = imread(r'C:\GIS\Vantaa\wcmask\Maskit\Kivisto-maski.png')
koiMask = imread(r'C:\GIS\Vantaa\wcmask\Maskit\Koivukyla-maski.png')
korMask = imread(r'C:\GIS\Vantaa\wcmask\Maskit\Korso-maski.png')
myrMask = imread(r'C:\GIS\Vantaa\wcmask\Maskit\Myyrmaki-maski.png')
tikMask = imread(r'C:\GIS\Vantaa\wcmask\Maskit\Tikkurila-maski.png')

# Reading in finnish stopwords, listing and setting them
swfp = r'C:\GIS\Vantaa\stopwords\stopwords_suomi.txt'
sWords = pd.read_csv(swfp, names=('a'), header=None)
stopwords = sWords['a'].values.tolist()
stopwords.extend(['alue','voida'])
stoplist = set(stopwords)

# Creating lists for WordCloud creation
alueet = [aviapolis, hakunila, kivisto, koivukyla, korso, myyrmaki, tikkurila]
nimet = ['aviapolis', 'hakunila', 'kivisto', 'koivukyla', 'korso', 'myyrmaki', 'tikkurila']
maskit = [aviMask, hakMask, kivMask, koiMask, korMask, myrMask, tikMask]

# Creating WordClouds for all Vantaa large regions
for i in range(len(alueet)):
    docs = alueet[i]['Fin_lemmas'].tolist()
    texts = [[word for word in document.lower().split() if word not in stoplist] for document in docs]
    texts = [[word.replace('vannaskoski', 'vantaankoski') for word in entry] for entry in texts]
    texts = [[word.replace('vannas', 'vantaa') for word in entry] for entry in texts]
    frequency = defaultdict(int)
    for text in texts:
        for token in text:
            frequency[token] += 1
    texts = [[token for token in text if frequency[token] > 1] for text in texts]
    wc = WordCloud(background_color='white', colormap='summer', max_words=100, width=1000, height=800, mask=maskit[i])
    wc.generate_from_frequencies(frequency)
    wc.to_file(r'C:\GIS\Vantaa\sanapilvet\{}_virkistys_wordcloud.png'.format(nimet[i]))
