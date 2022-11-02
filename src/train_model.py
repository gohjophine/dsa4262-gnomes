'''
(1) Model training

This script will use training data to learn a model.

Input:
1. direct RNA-Seq data, processed by m6Anet
2. m6A labels

Output:
1. the trained model
2. features used by the model
'''

import warnings
warnings.filterwarnings("ignore")

import sys
import pandas as pd
import numpy as np
import pickle
import json

from sklearn.model_selection import train_test_split, GroupShuffleSplit 
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler

from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier

from xgboost import XGBClassifier

from factor_analyzer import FactorAnalyzer

import itertools
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score, confusion_matrix, mean_squared_error, auc, roc_curve

try:
    path1 = sys.argv[1]
    path2 = sys.argv[2]
    with open(path1, "r") as read_file:
        data = [json.loads(l) for l in read_file]
    df2 = pd.read_csv(path2)
except:
    print('Error: please input a valid direct RNA-Seq data and m6A labels paths e.g. ../data/data.json ../data/data.info')
    sys.exit(1)

output_path = '../data/trial_model.sav'

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

# Join with Labels (data.csv)
data = df2
grouped_df3 = grouped_df3.merge(data, how='left', on=['transcript_id', 'transcript_position'])
gene_id_col = grouped_df3.pop('gene_id')
grouped_df3.insert(0, 'gene_id', gene_id_col)


## resample_normalise (dataset 0 sub).ipynb ##
# split dataset into X y train test, based on gene_id
# input: df, split_size
# output: train df, test df
def split_tt(df, split_size=0.2):
    splitter = GroupShuffleSplit(test_size=split_size, n_splits=1, random_state=42)
    split = splitter.split(df, groups=df['gene_id'])
    train_inds, test_inds = next(split)
    train = df.iloc[train_inds]
    test = df.iloc[test_inds]
    
    y_train = train['label']
    X_train = train.drop(['label', 'sevenmers'], axis = 1)
    y_test = test['label']
    X_test = test.drop(['label', 'sevenmers'], axis = 1)
    
    return X_train, y_train, X_test, y_test

# oversample and undersample such that ratio of minority to majority samples becomes 3:4
# input: df, df (X_train, y_train)
# output: df, df (resampled version)
def resample(X_train, y_train):
    # define oversampling strategy so that ratio of minority samples to majority samples is 1:2
    oversample = RandomOverSampler(sampling_strategy=0.5, random_state=42)
    X_train_over, y_train_over = oversample.fit_resample(X_train, y_train)
    
    # define undersampling strategy so that the ratio of minority to majority samples becomes 3:4
    under = RandomUnderSampler(sampling_strategy=0.75, random_state=42)
    X_train_under, y_train_under = under.fit_resample(X_train_over, y_train_over)
    return X_train_under, y_train_under

df = grouped_df3
features_nominal = ['order_1', 'order_2', 'order_3', 'order_6', 'order_7']
df[features_nominal] = df[features_nominal].astype('category')

# full processed data set
y_train = df['label']
X_train = df.drop(['label', 'sevenmers'], axis = 1)
X_train, y_train = resample(X_train, y_train)
X_train = X_train.drop(columns=['gene_id', 'transcript_id'])

# split processed data set into train set and test set
# keeping gene_id for cross validation in later section
X_train_id, y_train_id, X_test_id, y_test_id = split_tt(df)
X_train_id, y_train_id = resample(X_train_id, y_train_id)


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

# on full data
X_train_norm = normalise(X_train)
X_train_norm['transcript_position'] = X_train['transcript_position']
df_norm_d0 = normalise(df) # normalise for dataset 0 

# on train-test set
X_train_id = normalise(X_train_id)
X_test_id = normalise(X_test_id)


## feature_selection.ipynb ##

# one hot encode the categories
features_nominal = ['order_1', 'order_2', 'order_3', 'order_6', 'order_7']
df_norm_d0 = pd.get_dummies(df_norm_d0, columns = features_nominal)
X_train = pd.get_dummies(X_train_norm, columns = features_nominal)
X_train_id = pd.get_dummies(X_train_id, columns=features_nominal)
X_test_id = pd.get_dummies(X_test_id, columns=features_nominal)


impt_feat = []

## XGB feature importance
xgb = XGBClassifier(eval_metric = ['logloss'], use_label_encoder=False, random_state=42)
xgb.fit(X_train, y_train.values.ravel())

feats = {} # a dict to hold feature_name: feature_importance
for feature, importance in zip(X_train.columns, xgb.feature_importances_):
    feats[feature] = importance # add the name/value pair 

importances = pd.DataFrame(feats.items(), columns=['Feature', 'Importance'])
importances = importances.sort_values(by = ['Importance'], ascending = False)
impt_feat.extend(importances.Feature.iloc[0:30].tolist())

## PCA dimensionality reduction
# Remove categorical features
df_pca = df.drop(columns = ['label', 'sevenmers', 'gene_id', 'transcript_id', 'order_1', 'order_2', 'order_3', 'order_6', 'order_7'])
fa = FactorAnalyzer(n_factors = 10, method = 'principal', rotation='varimax')
fa.fit(df_pca)

eigenvalues, _ = fa.get_eigenvalues()
variances = fa.get_factor_variance()
def evaluate_pcs(num_of_pcs,data):
    def encode_vals(x):
        if x <= -0.7 or x >= 0.7:
            return x
        else:
            return("")    
    # REMARK: we use 'principal' method and 'varimax' rotation in the FactorAnalyzer function.
    f = FactorAnalyzer(n_factors=num_of_pcs, method = 'principal',rotation='varimax')
    f.fit(data)
    loadings = pd.DataFrame(f.loadings_).set_index(data.columns)
    loadings = loadings.applymap(encode_vals)
    loadingcols= list(loadings.columns)
    newcols = {}
    for i in loadingcols:
        newcols[i] = "PC" + str(i+1)
    loadings.rename(columns = newcols,inplace=True)
    return loadings

# The following function generates the rotation matrix. Recall that we use
# this matrix to determine if the PCs generated are easily understandable and appropriate.
# The argument "num_of_pcs" specifies, the number of PCs we wish to generate.

PCA_df = evaluate_pcs(9,df_pca)
value = []
which = lambda lst:list(np.where(lst)[0])
for i in PCA_df.columns:
    value.extend(which(PCA_df[i] == ''))
    
value = pd.DataFrame(value).rename(columns = {0: 'rowno'})
val_counts = pd.DataFrame(value.value_counts()).reset_index()
drop_cols = val_counts[val_counts[0] == 9]['rowno']
keep_df = df.drop(columns = list(PCA_df.index[drop_cols]))
impt_feat.extend(keep_df.columns)

## RFE Recursive Feature Elimination
logreg = LogisticRegression(max_iter=1000, random_state=42)
rfe3 = RFE(logreg, n_features_to_select=30)
rfe3 = rfe3.fit(X_train, y_train.values.ravel())

cols_keep = X_train.columns.values[rfe3.support_] # cols remaining
impt_feat.extend(cols_keep)

## Feature Importance using Random Forest
forest1 = RandomForestClassifier(random_state = 42, n_jobs= -1)
forest1.fit(X_train,y_train.values.ravel())
feats = {}
for feature, importance in zip(X_train.columns, forest1.feature_importances_):
    feats[feature] = importance

importances = pd.DataFrame(feats.items(), columns=['Feature', 'Importance'])
importances = importances.sort_values(by = ['Importance'], ascending = False)
impt_feat.extend(importances.Feature.iloc[0: 30].tolist())

# select feature counts >= 2
features = pd.DataFrame(impt_feat).rename(columns = {0: 'feat'})
feature_count = pd.DataFrame(features.value_counts()).reset_index()
impt_feat = feature_count[feature_count[0] >= 2]['feat']

# Remove Collinear Features
# https://stackoverflow.com/questions/29294983/how-to-calculate-correlation-between-all-columns-and-remove-highly-correlated-on
def remove_collinear_features(df_model, target_var, threshold, verbose, final_features):
    '''
    Objective:
        Remove collinear features in a dataframe with a correlation coefficient
        greater than the threshold and which have the least correlation with the target (dependent) variable. Removing collinear features can help a model 
        to generalize and improves the interpretability of the model.

    Inputs: 
        df_model: features dataframe
        target_var: target (dependent) variable
        threshold: features with correlations greater than this value are removed
        verbose: set to "True" for the log printing

    Output: 
        dataframe that contains only the non-highly-collinear features
    '''

    # Calculate the correlation matrix
    corr_matrix = df_model.drop(target_var, 1).corr()
    iters = range(len(corr_matrix.columns) - 1)
    drop_cols = []
    dropped_feature = ""

    # Iterate through the correlation matrix and compare correlations
    for i in iters:
        for j in range(i+1): 
            item = corr_matrix.iloc[j:(j+1), (i+1):(i+2)]
            col = item.columns
            row = item.index
            val = abs(item.values)

            # If correlation exceeds the threshold
            if val >= threshold:
                # Print the correlated features and the correlation value
                if verbose:
                    print(col.values[0], "|", row.values[0], "|", round(val[0][0], 2))
                col_value_corr = df_model[col.values[0]].corr(df_model[target_var])
                row_value_corr = df_model[row.values[0]].corr(df_model[target_var])
                if verbose:
                    print("{}: {}".format(col.values[0], np.round(col_value_corr, 3)))
                    print("{}: {}".format(row.values[0], np.round(row_value_corr, 3)))
                if col_value_corr < row_value_corr:
                    drop_cols.append(col.values[0])
                    dropped_feature = "dropped: " + col.values[0]
                else:
                    drop_cols.append(row.values[0])
                    dropped_feature = "dropped: " + row.values[0]
                if verbose:
                    print(dropped_feature)
                    print("-----------------------------------------------------------------------------")

    # Drop one of each pair of correlated columns
    drops = set(drop_cols)
    df_model = df_model.drop(columns=drops)
    final_features = final_features.extend(df_model.columns.tolist())

X_train = X_train[impt_feat]
df = pd.concat([X_train, y_train], axis=1)
final_features = []
remove_collinear_features(df, 'label', 0.9, False, final_features)
final_features.remove('label')

# print("final features: \n",final_features)
X_train = X_train[final_features]
X_train_id2 = X_train_id[final_features]
X_test_id2 = X_test_id[final_features]
df_norm_d0_2 = df_norm_d0[final_features]

## rf_model.ipynb ##
# baseline model with default parameters
forest = RandomForestClassifier(random_state = 42, n_jobs= -1)
forest.fit(X_train_id2, y_train_id)
# test metrics
rf_y_pred = forest.predict(X_test_id2)
y_predict_prob = forest.predict_proba(X_test_id2)[:, 1]

ra = metrics.roc_auc_score(y_test_id, y_predict_prob)
pa = metrics.average_precision_score(y_test_id, y_predict_prob)


# use a manually-defined cross-validation method to tune hyperparameters: max_depth and n_estimators
# input: number of folds, list of features to use, training set df, training labels df,
#           minimum max_depth, maximum max_depth, step size of max_depth,
#           minimum n_estimators, maximum n_estimators, step size of n_estimators
# output: df containing different combinations of parameters and result metrics
def manual_fold(k, features_lst, x_df, y_df, min_depth, max_depth, step_depth, min_trees, max_trees, step_trees):
    # use GroupShuffleSplit to group by gene_id and split into k folds
    gss = GroupShuffleSplit(n_splits=k, random_state = 42)
    folds = gss.split(X=x_df, groups=x_df['gene_id'])
    
    # get row indice of folds
    train_lst = []
    test_lst = []
    for train, test in folds:
        train_lst.append(train)
        test_lst.append(test)
    
    # select relevant features for model
    x_df = x_df[features_lst]
    
    # grid for hyperparameter tuning
    max_d = np.arange(min_depth, max_depth+1, step_depth).tolist()
    n_trees = np.arange(min_trees, max_trees+1, step_trees).tolist()
    a = [max_d, n_trees]
    space = list(itertools.product(*a))
    
    # hyperparameter tuning using cross validation
    # store parameters and result metrics into lists
    d = []
    t =[]
    accuracy = []
    precision = []
    recall = []
    roc_auc = []
    pr_auc = []
    
    for param in space:
        # fit model
        depth = param[0]
        trees = param[1]
        forest_cv = RandomForestClassifier(random_state = 42, n_jobs= -1, max_depth = depth, n_estimators = trees)
        
        acc = []
        pre = []
        rec = []
        rocauc = []
        prauc = []
        
        for i in range(k):

            cv_x_train = x_df.iloc[train_lst[i]]
            cv_x_test = x_df.iloc[test_lst[i]]

            cv_y_train = y_df.iloc[train_lst[i]]
            cv_y_test = y_df.iloc[test_lst[i]]

            forest_cv.fit(cv_x_train, cv_y_train.values.ravel())

            # test metrics
            rf_y_pred = forest_cv.predict(cv_x_test)

            a = metrics.accuracy_score(cv_y_test, rf_y_pred)
            p = metrics.precision_score(cv_y_test,rf_y_pred)
            r = metrics.recall_score(cv_y_test, rf_y_pred)

            y_predict_prob = forest_cv.predict_proba(cv_x_test)[:, 1]
            ra = metrics.roc_auc_score(cv_y_test, y_predict_prob)
            pa = metrics.average_precision_score(cv_y_test, y_predict_prob)

            acc.append(a)
            pre.append(p)
            rec.append(r)
            rocauc.append(ra)
            prauc.append(pa)
        
        d.append(depth)
        t.append(trees)
        accuracy.append(sum(acc)/k)
        precision.append(sum(pre)/k)
        recall.append(sum(rec)/k)
        roc_auc.append(sum(rocauc)/k)
        pr_auc.append(sum(prauc)/k)
        
    # make dataframe containing all combinations of parameters and their results
    dictionary = {"depth":d, "trees":t, "accuracy":accuracy, "precision":precision, "recall":recall, "roc-auc":roc_auc, "pr-auc":pr_auc}
    result = pd.DataFrame(dictionary)
    return result

# the search-grid used here has been already been refined 
res_df = manual_fold(5, final_features, X_train_id, y_train_id, 16, 18, 1, 120, 150, 10)

# obtain the best 5 models scored on the average between roc-auc and pr-auc
res_df['score'] = (res_df['roc-auc'] + res_df['pr-auc'])/2
best_five = res_df.sort_values('score', ascending = False).head().reset_index()

# fit models using these best parameters and select the one that performs the best on y-train
best_ind = 0
best_d = 0
best_n = 0
best_score = 0
# fit models using these best parameters and select the one that performs the best on y-train
for i in range(5):
    d = best_five['depth'].loc[i]
    n = best_five['trees'].loc[i]

    tuned_forest = RandomForestClassifier(random_state = 42, n_jobs= -1, max_depth = d, n_estimators = n)
    tuned_forest.fit(X_train_id2, y_train_id.values.ravel())
    
    y_predict_prob = tuned_forest.predict_proba(X_test_id2)[:, 1]
    ra = metrics.roc_auc_score(y_test_id, y_predict_prob)
    pa = metrics.average_precision_score(y_test_id, y_predict_prob)

    if (ra+pa)/2 > best_score:
        best_score = (ra+pa)/2
        best_ind = i
        best_d = d
        best_n = n

# fit model with best average score using best_d and best_n on full training data
forest = RandomForestClassifier(random_state = 42, n_jobs= -1, max_depth = best_d, n_estimators = best_n)
forest.fit(X_train, y_train.values.ravel())

# get predicted values on train set
y_predict_prob = forest.predict_proba(df_norm_d0_2)[:, 1]
ra = metrics.roc_auc_score(df_norm_d0['label'], y_predict_prob)
pa = metrics.average_precision_score(df_norm_d0['label'], y_predict_prob)

pickle.dump(forest, open(output_path, 'wb'))
pickle.dump(final_features, open('../data/trial_features.info', 'wb'))

print()
print(f'Execution successful! The model can now be found in {output_path}')
print()
