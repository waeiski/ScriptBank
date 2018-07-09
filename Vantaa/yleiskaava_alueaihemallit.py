# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 10:55:04 2018

REQUIREMENTS
============
Python 3.5 or newer
GeoPandas
Pandas
Gensim

DESCRIPTION
===========
This script reads in lemmatized comments from Vantaa Masterplan 2020 map questionnaire
per area and does topic modelling for each of them and saves the results. Topic
modelling is done using the Gensim LDA pipeline. This script outputs pickled
dataframes and CSV tables, but also saves the resulting dictionaries, corpora
and models.

USAGE
=====
1. Make sure all filepaths fit your current system.
2. Make sure the texts in all files are are in a field called 'Fin_lemmas' or change that part of
the script to fit whatever your field name is. The field name for text must be
exactly the same across all GeoPackage files.
3. Run the script
4. Go grab coffee
5. Return to your topic modelling results

QUOTE OF THE SCRIPT
===================
"Look on my works, ye mighty, and despair!"
    - Percy Bysshe Shelley

@author: Tuomas Väisänen, Vantaan kaupunki 2018
'"""

import geopandas as gpd
import pandas as pd
from collections import defaultdict
import warnings
warnings.filterwarnings(action='ignore', category=UserWarning, module='gensim')
import logging
logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
from gensim import corpora, models

# Initializing the filepaths
av = r'C:\GIS\Vantaa\points\Lems\Lems_Aviapolis.gpkg'
hk = r'C:\GIS\Vantaa\points\Lems\Lems_Hakunila.gpkg'
kv = r'C:\GIS\Vantaa\points\Lems\Lems_Kivisto.gpkg'
kk = r'C:\GIS\Vantaa\points\Lems\Lems_Koivukyla.gpkg'
ko = r'C:\GIS\Vantaa\points\Lems\Lems_Korso.gpkg'
my = r'C:\GIS\Vantaa\points\Lems\Lems_Myyrmaki.gpkg'
tk = r'C:\GIS\Vantaa\points\Lems\Lems_Tikkurila.gpkg'

# Reading the files in
aviapolis = gpd.read_file(av)
hakunila = gpd.read_file(hk)
kivisto = gpd.read_file(kv)
koivukyla = gpd.read_file(kk)
korso = gpd.read_file(ko)
myyrmaki = gpd.read_file(my)
tikkurila = gpd.read_file(tk)

# Reading in finnish stopwords, listing and setting them
swfp = r'C:\GIS\Vantaa\stopwords\stopwords_suomi.txt'
sWords = pd.read_csv(swfp, names=('a'), header=None)
stopwords = sWords['a'].values.tolist()
stopwords.extend(['alue','voida'])
stoplist = set(stopwords)

# Creating a dataframe list and a name list for topic model creation
alueet = [aviapolis, hakunila, kivisto, koivukyla, korso, myyrmaki, tikkurila]
nimet = ['aviapolis', 'hakunila', 'kivisto', 'koivukyla', 'korso', 'myyrmaki', 'tikkurila']

# Initializing LDA topic model from gensim
lda = models.ldamodel.LdaModel

# For loop that tokenizes, models topics and saves results
for i in range(len(alueet)):
    '''
    First, the texts are tokenized
    '''
    docs = alueet[i]['Fin_lemmas'].tolist()
    texts = [[word for word in document.lower().split() if word not in stoplist] for document in docs]
    frequency = defaultdict(int)
    for text in texts:
        for token in text:
            frequency[token] += 1
    texts = [[token for token in text if frequency[token] > 1] for text in texts]
    tokentexts = pd.Series(texts)
    alueet[i]['tokens'] = tokentexts.values
    '''
    Second, dictionaries and corpora are formed with gensim and saved
    '''
    dictionary = corpora.Dictionary(alueet[i]['tokens'])
    dictionary.save(r'C:\GIS\Vantaa\Aihemallit\dictionaries\{}_dictionary.dict'.format(nimet[i]))
    alueet[i]['bow'] = alueet[i]['tokens'].map(dictionary.doc2bow)
    corpus = alueet[i]['bow'].tolist()
    corpora.MmCorpus.serialize(r'C:\GIS\Vantaa\Aihemallit\corpora\{}_corpus.mm'.format(nimet[i]), corpus)
    '''
    Third, areal topic models are created and saved, this part might take a 30 minutes 
    '''
    alue_lda = lda(corpus, num_topics=10, id2word=dictionary, iterations=6500, passes=300)
    alueCm = models.coherencemodel.CoherenceModel(texts=texts, model=alue_lda, corpus=corpus, dictionary=dictionary, coherence='c_v', topn=10)
    topicCoh = alueCm.get_coherence_per_topic()
    aluetopics = alue_lda.show_topics(num_topics=10, num_words=10, formatted=True)
    alue_lda.save(r'C:\GIS\Vantaa\Aihemallit\models\{}_10topics_6500iters_300pass_38coh.model'.format(nimet[i]))
    '''
    Fourth, areal topic model results are appended into corresponding dataframes
    '''
    alueet[i]['topic_results'] = [list(alue_lda[x]) for x in alueet[i]['bow']]
    # Extracting the top 3 most prevalent topics
    top3_res = []
    for row, tops in alueet[i]['topic_results'].iteritems():
        top3 = sorted(tops,key=lambda x: x[1], reverse=True)[:3] # Sorting Top 3 topics from best to worst
        top3_res.append(top3)
    alueet[i]['top3_topics'] = top3_res
    # Extracting the best topic and probability
    top_res = []
    top_prob = []
    for row, tops in alueet[i]['topic_results'].iteritems():
        top = sorted(tops,key=lambda x: x[1], reverse=True)[0]
        top_res.append(top[0])
        top_prob.append(top[1])
    alueet[i]['best_topic'] = top_res
    alueet[i]['best_prob'] = top_prob
    # Saving the dataframe as pickle
    picklepath = r'C:\GIS\Vantaa\Aihemallit\pickles\{}_vanyk_topics.pkl'.format(nimet[i])
    alueet[i].to_pickle(picklepath)
    # Saving the point results as CSV table to be joined by ID into GeoPackage
    selected = alueet[i][['id','Fin_lemmas','best_topic','best_prob']]
    topfp = r'C:\GIS\Vantaa\Aihemallit\topics\{}_points_topics.txt'.format(nimet[i])
    selected.to_csv(topfp, sep='\t', encoding='utf-8')
    '''
    Fifth, most contributing words-per-topic and corresponding coherences are
    extracted and then saved into corresponding CSV table
    '''
    # Extracting the topic words
    results = []
    for t in aluetopics:
        words = t[1]
        words = words.split(' + ')
        words = [w.split('*')[1].strip('"') for w in words]
        results.append(words)
    # Adding the corresponding coherences for each topic
    for pos in range(len(results)):
        results[pos].extend([topicCoh[pos]])
    # Saving the 10 topics with words and coherences into a csv file
    topicDf = pd.DataFrame.from_records(results)
    topicsave = r'C:\GIS\Vantaa\Aihemallit\topics\10_{}_topics_with_coherences.txt'.format(nimet[i])
    topicDf.to_csv(topicsave, encoding='utf-8', sep='\t')
    