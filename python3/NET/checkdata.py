# -*- coding: utf-8 -*-
"""

@author: TK

<checkdata.py>: inital data upload-- clean and tidy up etc for further analysis. Output: new .csv to be uploaded elsewhere.

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

#%% Load dataset
filename = "data_for_exercise.csv"
df = pd.read_csv(filename)
print(df.info())
print(df.index) #check index
print(df.columns) #check default column names
print(df.shape) # shape
print(df.head(2)) # first rows

# after inspection, import as ; sep and comma decimal?
df = pd.read_csv(filename, sep=';', decimal=',')
print(df.index)
print(df.columns) #check default column names
print(df.shape) # shape
print(df.dtypes) #data type?
print(df.head(2)) # first rows

# NOTE: some numbers separated by commas, some by points. Need to make data type consistent (numeric).

#%% explore & tidy up data:
    
print(df.isnull().sum()) #missing data?
if sum(df.isnull().sum()) == 0:
    print('No missing data: proceed!')
else:
    sys.exit('Some missing data: explore further before proceeding...')
    
# column names 
cols = ['ts', 'solar', 'wind', 'cons', 'market', 'spot']      
df.columns = cols #rename columns

# set time column to datetime and then to index
df['ts'] = pd.to_datetime(df['ts'], utc=True)
print(df['ts'].min(), df['ts'].max())
df = df.set_index('ts')

# make sure var values are numeric
df[['solar','wind','cons']] = df[['solar','wind','cons']].apply(pd.to_numeric)
print(df.dtypes)

#Check head again:
print(df.head(2))

# summary
print(df.describe())
#Looks reasonable...

#%% save new data file
fn_new = 'data_new.csv'
df.to_csv(fn_new)
print('Data tidied and saved as...', fn_new)
