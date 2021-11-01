# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 10:36:53 2021

@author: TK

         >>>> Norlys ET -- Quant exercise <<<<

One of your colleagues have discovered a trading strategy that looks to be quite profitable. The strategy works by opening a position on the DAH (Day Ahead) market, and then trade this out ID (intraday). It’s however starting to decline in performance, so you’re tasked with coming up with a new strategy. 

The new strategy shall be based on the forecasted wind and solar production as well as the forecasted consumption. The original strategy also contains derived features (an example could be difference in solar level from rolling two-week average). To compare the two strategies directly, you’ve been given a file (data_for_exercise.csv), that contains the following columns. 

data_input.columns = [ts, wind, solar, cons, spot, market] 

> ts: Timestamp [-] 
> wind, solar, cons: Forecasted wind, solar and consumption [MW]. 
> spot: The price we’re able to open the positions at [EUR/MWh]. 
> market: The price we’re able to close the positions at [EUR/MWh]

NOTE: Python and module versions as follows...
> Python: 3.8.8 
> scipy: 1.7.1
> numpy: 1.21.2
> matplotlib: 3.4.3
> pandas: 1.3.3
> sklearn: 0.24.2
"""
#%% load modules etc
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#from sci-kit-learn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import mean_squared_error
# ML models
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVC
from sklearn.svm import SVR

#%% Load dataset
filename = "data_new.csv"
df = pd.read_csv(filename, parse_dates=['ts'])

# set time column to index
df = df.set_index('ts')

# summary check
print(df.describe())


#%% resample dfs to coarser time scales - may be useful for plotting/comparison
daily_df = df.resample('D').mean()
weekly_df = df.resample('W').mean()
fortnight_df = df.resample('2W').mean()

# Derived features: difference from rolling n-day aves
def rolling_diff(dataframe, name , n):
    roll_ave = dataframe[name].rolling(n, center=True).mean()
    diff = dataframe[name] - roll_ave
    return diff
    
# e.g. rolling 14 day aves for energy forecasts, as per the strategy doc
solar_14d = df['solar'].rolling(14, center=True).mean()
wind_14d = df['wind'].rolling(14, center=True).mean()
cons_14d = df['cons'].rolling(14, center=True).mean()


#%% new columns of interest?
df0 = df.copy() # make a copy of original df
df['dwind_14d'] = rolling_diff(df,'wind',14)
df['dsolar_14d'] = rolling_diff(df,'solar',14)
df['dcons_14d'] = rolling_diff(df,'cons',14)
# df['total'] = df['solar'] + df['wind'] # total energy
# df['dE'] = df['cons'] - df['total'] # difference in energy
df['dP'] = df['market'] - df['spot'] # difference in price
# df['Sp+'] = np.where(df['spot'] > 0, 1, 0)
# df['Ma+'] = np.where(df['market'] > 0, 1, 0)
df['dP+'] = np.where(df['market'] > df['spot'], 1., 0.)
#NOTE: if market > spot, we make profit. So dP+ = 1 is the goal.

# summary
print(df.describe())


#%% 2. Explore data: summary plots
# box and whisker plots (from original dataframe df0)
df0.plot(kind='box', subplots=True, sharex=False, sharey=False)
plt.show()

# scatter plot matrix with histograms
pd.plotting.scatter_matrix(df)
plt.show()


def correlated_mat(dataframe, plot=False):
    corr = dataframe.corr()
    if plot:
        # sns.set(rc={'figure.figsize': (10, 10)})
        sns.heatmap(corr, cmap=plt.cm.RdBu,
                    vmin=-1.0, vmax=1.0, annot=True, linewidths=.7)
        plt.xticks(rotation=60, size=10)
        plt.yticks(rotation=60, size=10)
        plt.title('Correlation Matrix', size=20)
        plt.show()

correlated_mat(df, plot=True)

#%% seaborn plots

# choose variables of interest
sns.jointplot(x='dsolar_14d', y='market', data=df)

#recall dP+ = 1 if spot > market price, 0 if spot < market
sns.FacetGrid(df, hue='dP+') \
   .map(plt.scatter, 'solar', 'wind') \
   .add_legend()
plt.show()

sns.pairplot(df, hue='dP+')
plt.show()

sns.FacetGrid(df, hue='dP+') \
   .map(sns.kdeplot, 'wind') \
   .add_legend()
plt.show()

#%% time series
#energy forecasts
for name in ['cons', 'wind', 'solar']:
    df[name].plot(ylabel='MWh',alpha = 0.4)
    daily_df[name].plot(legend=True)
    
#energy forecasts
name1 = ['cons', 'wind', 'solar']
df[name].plot(ylabel='MWh',alpha = 0.4)
df[name].plot(legend=True)

#prices
for name in ['market', 'spot']:
    df[name].plot(ylabel='EUR/MWh',alpha = 0.4)
    daily_df[name].plot(legend=True)

#compare e.g. wind and market price
plt.figure()   
ax = df['market'].plot(color='r',ylabel='EUR/MWh',alpha = 0.4)
daily_df['market'].plot(color='r',legend=True)
df['wind'].plot(color='b',secondary_y=True,ylabel='MWh',alpha = 0.4)
daily_df['wind'].plot(color='b',secondary_y=True,legend=True)
ax.set_ylabel('EUR/MWh')
ax.right_ax.set_ylabel('MWh')
plt.show()

# Autocorrelation (temporal)
pd.plotting.autocorrelation_plot(df['wind'])
pd.plotting.autocorrelation_plot(df['dP'])

#%% 3. MODEL BUILDING - predicting dP+
#choose predictors and predictand
# predictors = ['solar', 'wind', 'cons']
predictors = ['solar', 'wind', 'cons','dsolar_14d', 'dwind_14d', 'dcons_14d']
predictand = ['dP+'] # Binary Y/N: market price > spot price.

df = df.dropna() # drop nans if they are present (e.g., due to rolling aves)

X = df[predictors].values # predictors 
y = df[predictand].values # predictand y = f(X)
y = np.ravel(y)
X_train, X_validation, Y_train, Y_validation = train_test_split(X, y, test_size=0.2, random_state=0)

#%% sklearn ML models -- i.e., the f in y = f(X). 
models = []
models.append(('LogR', LogisticRegression(solver='liblinear', multi_class='ovr')))
# models.append(('LinR', LinearRegression()))
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('KNN', KNeighborsClassifier(n_neighbors=10)))
models.append(('CART', DecisionTreeClassifier()))
models.append(('NB', GaussianNB()))
# models.append(('SVM', SVC(gamma='auto')))

#%% evaluate each model in turn
results = []
names = []
for name, model in models:
 	kfold = StratifiedKFold(n_splits=10, random_state=1, shuffle=True)
 	cv_results = cross_val_score(model, X_train, Y_train, cv=kfold, scoring='accuracy')
 	results.append(cv_results)
 	names.append(name)
 	print('%s: %f (%f)' % (name, cv_results.mean(), cv_results.std()))
     
#%% Compare Algorithms
plt.boxplot(results, labels=names)
plt.title('Algorithm Comparison')
plt.show()

#%% Make predictions on validation dataset
# model = SVC(gamma='auto')
# model = LogisticRegression()
model = models[2]
print('Predictions using %s model...' % names[2])
model[1].fit(X_train, Y_train)
predictions = model[1].predict(X_validation)

#%% Evaluate predictions
print(accuracy_score(Y_validation, predictions))
print(confusion_matrix(Y_validation, predictions))
print(classification_report(Y_validation, predictions))    
#%% 
'''
#%% CONTINUOUS PREDICTIONS..?

#%% 2. MODEL BUILDING - predicting dP
#choose predictors and predictand
predictors = ['solar', 'wind', 'cons']
predictand = ['market'] 


X = df[predictors].values # predictors 
y = df[predictand].values # predictand y = f(X)
y = np.ravel(y)

X_train, X_validation, Y_train, Y_validation = train_test_split(X, y, test_size=0.2, random_state=7)

#%% sklearn ML models -- i.e., the f in y = f(X). 
models = []

models.append(('LinR', LinearRegression()))
models.append(('Las', Lasso()))
models.append(('Ridge', Ridge(alpha=.5)))
models.append(('KNR', KNeighborsRegressor()))
models.append(('DTR', DecisionTreeRegressor()))
models.append(('RFR', RandomForestRegressor()))
models.append(('SVR', SVR()))

#%% evaluate each model in turn
results = []
names = []
for name, model in models:
 	kfold = KFold(n_splits=10)
 	cv_results = cross_val_score(model, X_train, Y_train, cv = kfold)
 	results.append(cv_results)
 	names.append(name)
 	print('%s: %f (%f)' % (name, cv_results.mean(), cv_results.std()))
     
#%% Compare Algorithms
plt.boxplot(results, labels=names)
plt.title('Algorithm Comparison')
plt.show()

#%% Make predictions on validation dataset
# model = SVC(gamma='auto')
# model = LogisticRegression()
model = models[3]
print('Predictions using %s model...' % model[0])
model[1].fit(X_train, Y_train)
predictions = model[1].predict(X_validation)

print(model[1].score(X_validation, Y_validation), mean_squared_error(Y_validation, predictions))
'''
