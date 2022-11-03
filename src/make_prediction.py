'''
(2) Predictions

This script will apply your model on new data and return the m6A predictions.

Input:
1. direct RNA-Seq data, processed by m6Anet
note: the pre-trained model from the first script needs to be present

Output:
1. predicted m6A sites
'''

import sys
import pandas as pd
import numpy as np
import pickle
import json

from sklearn import preprocessing
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler

from sklearn.ensemble import RandomForestClassifier

try:
    model = pickle.load(open('../data/model.sav', 'rb'))
except:
    print('Error: please run train_model.py before running this py file')

try:
    path = sys.argv[1]
    with open(path, "r") as read_file:
        data = [json.loads(l) for l in read_file]
except:
    print('Error: please input a valid direct RNA-Seq data path e.g. ../data/test_data.json')
    sys.exit(1)


output_path = sys.argv[2]


## read_json.ipynb ##

# return the key and value of the dictionary separately
# input: dictionary
# output: key, value
def split(x):
    key = list(x.keys())[0]
    value = x.get(key)
    return key, value

# expand the nested dictionary
transcript_list, position_list, sevenmers_list, reads_list = [], [], [], []
for line in data:
    t, dic = split(line)
    p, dic = split(dic)
    s, reads = split(dic)
    
    for read in reads:
        transcript_list.append(t)
        position_list.append(p)
        sevenmers_list.append(s)
        reads_list.append(read)

df = pd.DataFrame.from_dict({'transcript_id': transcript_list, 
                             'transcript_position': position_list, 
                             'sevenmers': sevenmers_list, 
                             'reads': reads_list})

df1 = pd.DataFrame(df['reads'].to_list())
df3 = pd.concat([df.reset_index(), df1], axis=1)
df3 = df3.rename(columns = {0:'dwelling_time_1', 1:'sd_current_1', 2:'mean_current_1',
                            3:'dwelling_time_2', 4:'sd_current_2', 5:'mean_current_2',
                            6:'dwelling_time_3', 7:'sd_current_3', 8:'mean_current_3'})

final_df = df3.drop(columns=['reads', 'index'])


## feature_engineering.ipynb ##

# input_path = '../data/data.csv'
# output_path = '../data/grouped_data.csv'

df = final_df
# get the difference of each features 
df['diff_dwelling_time_1'] = df['dwelling_time_2'] - df['dwelling_time_1']
df['diff_dwelling_time_2'] = df['dwelling_time_3'] - df['dwelling_time_2']
df['diff_sd_current_1'] = df['sd_current_2'] - df['sd_current_1']
df['diff_sd_current_2'] = df['sd_current_3'] - df['sd_current_2']
df['diff_mean_current_1'] = df['mean_current_2'] - df['mean_current_1']
df['diff_mean_current_2'] = df['mean_current_3'] - df['mean_current_2']
# aggregate and find min, max, median, std for each features

lst = ['transcript_id', 'transcript_position', 'sevenmers']

grouped_df = df.groupby(by = lst).agg(
    {'dwelling_time_1': [min, max, 'median', 'std'],
     'sd_current_1': [min, max, 'median', 'std'],
     'mean_current_1': [min, max, 'median', 'std'],
     'dwelling_time_2': [min, max, 'median', 'std'],
     'sd_current_2': [min, max, 'median', 'std'],
     'mean_current_2': [min, max, 'median', 'std'],
     'dwelling_time_3': [min, max, 'median', 'std'],
     'sd_current_3': [min, max, 'median', 'std'],
     'mean_current_3': [min, max, 'median', 'std'],
     'diff_dwelling_time_1': [min, max, 'median', 'std'],
     'diff_dwelling_time_2': [min, max, 'median', 'std'],
     'diff_sd_current_1': [min, max, 'median', 'std'],
     'diff_sd_current_2': [min, max, 'median', 'std'],
     'diff_mean_current_1': [min, max, 'median', 'std'],
     'diff_mean_current_2': [min, max, 'median', 'std'],
    }).reset_index()

# rename the columns
grouped_df2 = grouped_df
grouped_df2.columns = ["_".join(x) for x in np.ravel(grouped_df.columns)]
grouped_df2 = grouped_df.rename(columns = {'transcript_id_': 'transcript_id', 
                                           'transcript_position_': 'transcript_position',
                                           'sevenmers_': 'sevenmers'})

# find the relative position of each read in each transcript
grouped_df2['transcript_position'] = grouped_df2['transcript_position'].astype(int)
grouped_df2['relative_position'] = grouped_df2.groupby('transcript_id')['transcript_position'].transform(lambda x: (x - x.min())/(x.max()-x.min()))
# note: have NAs because there's transcripts with only one position
# fill the NAs with 0
grouped_df2['relative_position'] = grouped_df2['relative_position'].fillna(0)
# split the sevenmers into seven columns
order_df = pd.DataFrame(grouped_df2['sevenmers'].str.split('').to_list())[[x for x in range(1, 8)]].rename(
    columns = {1: 'order_1', 2:'order_2', 
               3: 'order_3', 4:'order_4', 
               5: 'order_5', 6:'order_6',
               7: 'order_7'})
grouped_df3 = pd.concat([grouped_df2, order_df], axis = 1)

# order_4 is always A, order_5 is always C --> drop
grouped_df3 = grouped_df3.drop(columns=['order_4', 'order_5'])

# find the number of occurrence of a letter in a word
# input: str word, str letter
# output: int
def find(word, letter):
    res = 0
    for i in word:
        if i==letter:
            res += 1
    return res
# count the A,C,G,T in the sevenmers
grouped_df3['count_A'] = grouped_df3['sevenmers'].map(lambda x: find(x, 'A'))
grouped_df3['count_C'] = grouped_df3['sevenmers'].map(lambda x: find(x, 'C'))
grouped_df3['count_G'] = grouped_df3['sevenmers'].map(lambda x: find(x, 'G'))
grouped_df3['count_T'] = grouped_df3['sevenmers'].map(lambda x: find(x, 'T'))


## resample_normalise.ipynb ##

df = grouped_df3
features_nominal = ['order_1', 'order_2', 'order_3', 'order_6', 'order_7']
df[features_nominal] = df[features_nominal].astype('category')

# normalise the numerical columns
# input: df
# output: normalised df
def normalise(df):
    numerical_columns = df.select_dtypes(include=['int64', 'float64']).columns
    string_columns = df.select_dtypes(include=['object', 'category']).columns
    
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(df[numerical_columns])
    
    df_normalised = pd.DataFrame(x_scaled)
    df_normalised.columns = numerical_columns
    
    final_df = pd.concat([df[string_columns].reset_index(), df_normalised], axis=1)
    final_df = final_df.drop(columns = ['index'])
    return final_df

df_norm = normalise(df)
df_norm['transcript_position'] = df['transcript_position']

## feature_selection.ipynb ##
# one hot encode the categories
features_nominal = ['order_1', 'order_2', 'order_3', 'order_6', 'order_7']
df_norm = pd.get_dummies(df_norm, columns = features_nominal)

selected = pickle.load(open('../data/features.info', 'rb'))
df_norm = df_norm[selected]

# make predictions

y_pred = model.predict_proba(df_norm)[:, 1]

def join_pred(df, y_pred):
    new_df = pd.concat([df[['transcript_id', 'transcript_position']], pd.DataFrame(y_pred)], axis=1)
    new_df = new_df.rename(columns = {0: 'score'})
    return new_df

y_df = join_pred(df, y_pred)

# Export to csv
y_df.to_csv(output_path, index=False)
print()
print(f'Excecution successful! The prediction can now be found in {output_path}')
print()

print(y_df.head())