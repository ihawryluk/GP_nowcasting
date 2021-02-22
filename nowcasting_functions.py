# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 21:02:47 2021

@author: iwona
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pystan
import datetime
import arviz
import sklearn.metrics
import pickle

def test(test):
    print(test + test)
    
def fit_model(raw_data, precompiled=True, modelname='model_4comp_SE', date_nowcast=None,
              maxD=10, iters=1000, warmup=400, chains=4, adapt_delta=0.9,
              max_treedepth=12, seed=9876,
              pickle_run=False, save=True, savepath=''):
    
    # Prepare the data for nowcast
    data = raw_data.copy()
    data['Release_index'] = data.Release.astype('category').cat.codes
    data = data[data['Date'] >= '2020-06-30'] # cut off early days as they are less relevant
    data = data[data['Release'] >= '2020-06-30'] # cut off early days as they are less relevant
    
    if date_nowcast is None:
        date_nowcast = data['Date'].values[-1]
    
    delays_data, triangle = nowcasting_prep(raw_data, date_nowcast, maxD)

    # Get the pre-compiled model
    if precompiled:
        model = get_model(modelname)
    else:
        #model = pystan.StanModel(file=modelname)
        model = modelname
        
    # Run the model
    fit = fit_single_model(model, triangle, iters =iters,
              warmup=warmup, chains=chains, adapt_delta=adapt_delta,
              max_treedepth=max_treedepth, seed=seed)
    
    # Save all the samples
    if pickle_run:
        with open(savepath + model + '_' + date_nowcast + ".pkl", "wb") as f:
            pickle.dump({'fit' : fit}, f, protocol=-1)
            
    # Get and saved processed results
    if save:
        savefile = savepath + modelname + '_' + date_nowcast + '.csv'
        
    results = get_results_df(fit, savefile=savepath)
    
    return fit, results

##############################################################################
# DATA PREPARATION FUNCTIONS
##############################################################################

def sum_two_values(val1, val2):
    """Sums two values, treats NaN as 0 if one value is NaN, returns NaN if both values are NaN.
    Used for summing vectors of deaths."""
    if np.isnan(val1) & np.isnan(val2):
        return np.nan
    if np.isnan(val2):
        return val1
    return val1 + val2

def sum_two_vectors(vec1, vec2):
    """Used for getting a total number of deaths per day reported."""
    vec1 = vec1.values
    vec2 = vec2.values
    # if both values are NaN, sum should give NaN
    # if only one is NaN, this NaN should be treated as 0
    for i in range(len(vec1)):
        vec1[i] = sum_two_values(vec1[i], vec2[i])
    return vec1

#def get_weeks(date_now, date_last, all_weeks):
#    """Get number of weeks elapsed from date_last to date_now"""
#    date_now = datetime.datetime.strptime(date_now,'%Y-%m-%d')
#    weeks_from_end = int((date_last - date_now).days / 7)
#    return all_weeks - weeks_from_end
    
def get_weeks_vector(n):
    week_vec = []
    week=0
    counter=0
    for i in range(n):
        week_vec.append(week)
        counter += 1
        if counter == 7:
            counter=0
            week = week+1
    return week_vec
    
def create_reporting_triangle(raw_data):
    """Creates a reporting triangle from SIVEP data binned by date and release"""
    # create and empty Dataframe to store the reporting triangle
    index = raw_data.Date.unique()  # all Dates in the dataset
    columns = np.sort(raw_data.Release_index.unique()) # maximum number of weeks (releases)
    triangle = pd.DataFrame(index = index, columns = columns)
    
    # fill in the delays info
    for i in triangle.index:  # iterate over all dates
        current_date = raw_data[raw_data['Date'] == i]  # get reported deaths for this date for every release
        all_deaths = current_date.Deaths.values
        # take the first number of deaths appearing -- this is 0 delay
        triangle.loc[i,0] = all_deaths[0]
        next_releases = current_date.Release_index.values
        first_release = next_releases[0]
        for j in range(1,len(all_deaths)):
            triangle.loc[i,next_releases[j]-first_release] = all_deaths[j] - all_deaths[j-1]
#     triangle[triangle < 0] = 0
    return triangle

def sum_up_delays_above_max(triangle, maxD=10):
    """Sums up all deaths reported with delay larger than max_delay in the reporting triangle"""
    triangle_max_D = triangle.copy()
    max_delays_reported = int(triangle_max_D.columns[-1])
    for c in range(maxD+1, max_delays_reported+1):
        triangle_max_D.loc[:,maxD] = sum_two_vectors(triangle_max_D.loc[:,maxD],
                                                        triangle_max_D.loc[:,c])
        triangle_max_D.drop(columns = [c], inplace = True)
    return triangle_max_D

def mask_bottom_triangle(unmasked_triangle, val=0):
    """Adds a mask to obtain a correct reporting triangle;
    issue arising due to the fact that releases are not in equal weekly intervals"""
    nrow = unmasked_triangle.shape[0]
    ncol = unmasked_triangle.shape[1]
    masked_triangle = unmasked_triangle.copy()
    
    for c in range(0, ncol-1):
        masked_triangle.iloc[nrow-c-1,c+1] = val
    return masked_triangle

def test_reporting_triangle(triangle):
    """Tests whether the reporting triangle has a correct structure.
    Returns True if the structure is correct and False if it is not."""
    test_triangle = triangle.copy()  
    test_triangle.replace(0, np.nan, inplace=True) # might cause issues if 0s are outside of the missing information
    nans_per_row = test_triangle.isnull().sum(axis=1)

    for i in range(1,len(nans_per_row.values)):
        if nans_per_row[i] - nans_per_row[i-1] > 1:
            return False
        if nans_per_row[i] > 0:
            if nans_per_row[i] - nans_per_row[i-1] != 1:
                return False
    return True    

def mask_below_diagonal(triangle):
    """Makes all values below nan value to be nan"""
    nrow = triangle.shape[0]
    ncol = triangle.shape[1]
    masked_triangle = triangle.copy()
    
    for c in range(1, ncol):
        nan_idx = np.argwhere(np.isnan(masked_triangle.iloc[:,c].values))
        if nan_idx:
            for r in range(nan_idx[0][0], nrow):
                masked_triangle.iloc[r,c] = np.nan
    return masked_triangle

    
def nowcasting_prep(data, now_date, maxD=10):
    # Prepare the data for nowcasting
    data_filtered = data[data['Release'] <= now_date].copy() # remove the data above nowcasting date

    # 1. Change the Release index to start from 0
    releases = np.sort(data_filtered.Release_index.unique())
    for i in range(len(releases)):
        data_filtered.loc[data_filtered["Release_index"] == releases[i], "Release_index"] = i

    # 2. Create a reporting triangle for the validation data
    delays_data = create_reporting_triangle(data_filtered)

    # 3. Sum up the deaths above max_D (maximum delay to be used)
    delays_data = sum_up_delays_above_max(delays_data, maxD)

    # 4. Bin the data into weeks
    first_date = datetime.datetime.strptime(delays_data.index[0],'%Y-%m-%d')
    last_date = datetime.datetime.strptime(delays_data.index[-1],'%Y-%m-%d')
    delays_data['week'] = get_weeks_vector(len(delays_data))

    #all_weeks = int((last_date - first_date).days / 7)  #+1?

    #for i in delays_data.index:
    #   delays_data.loc[i,'week'] = get_weeks(i, last_date, all_weeks)

    delays_data_weekly = delays_data.groupby(['week']).sum(min_count=1)

    # 5. Mask some values
    # this is done to keep a correct structure of the reporting triangle
    # this issue is cause because the columns are 'releases' and rows are 'weeks', and the release tend to be weekly
    # with some shorter/longer periods in the beginning of August
    delays_data_weekly = mask_bottom_triangle(delays_data_weekly)

    # 5.5 change the begative values into 0s
    delays_data_weekly[delays_data_weekly < 0] = 0


    # 6. Add a total number of deaths per date reported
    delays_data_weekly["all_deaths"] = delays_data_weekly.sum(axis=1)
    
    triangle = delays_data_weekly.iloc[:,:-1].copy()
    triangle = mask_bottom_triangle(triangle, np.nan)
    triangle = mask_below_diagonal(triangle)
    triangle.replace(np.nan, 10000000, inplace=True)
    triangle = triangle.astype('int32')

    return delays_data_weekly, triangle


##############################################################################
# MODEL FITTING FUNCTIONS
##############################################################################
# def get_model(modelname):
#     """gets a precompiled stan model"""
#     with open("paper_models/backtesting/stan_models/" + modelname + ".pkl", "rb") as f:
#         data_dict = pickle.load(f)
#     model = data_dict['model']
#     return model

def get_model(modelname):
    """gets a precompiled stan model"""
    with open("stan_models/" + modelname + ".pkl", "rb") as f:
        data_dict = pickle.load(f)
    model = data_dict['model']
    return model


def fit_single_model(model, triangle, iters =1000,
              warmup=500, chains=4, adapt_delta=0.9,
              max_treedepth=12, seed=12345):
    """fits model to the reporting triangle"""
    
    T =  triangle.shape[0]
    D = triangle.shape[1]
    x = list(range(0,T))
    time = list(range(0,T))
    delay = list(range(0,D))
    
    data_stan = {'T': T, 'D': D, 'n': triangle, 'x': x, 'time' : time, 'delay' : delay}
    
    fit = model.sampling(data=data_stan,iter=iters,
                                   warmup=warmup, 
                                   chains=chains,
                                   thin=1, 
                                   n_jobs=-1, 
                                   control=dict(adapt_delta = adapt_delta, max_treedepth=max_treedepth), 
                                   seed = seed)
    return fit


def fit_to_all_dates(model, model_str, triangles, dates_to_test):
    """fits the model to each of the dates in the dates_to_test list"""
    print('Fitting started for', model_str)
    fit_all = {}
    for d in dates_to_test:
        n = triangles[d]
        try:
            fit = fit_single_model(model, n)
            fit_all.update({d: fit})
            print('Fit finished for date', d)
        except:
            print('Fitting failed for', model_str, 'date', str(d))   
    return fit_all

##############################################################################
# SAVING MODELS OUTPUTS
##############################################################################

def get_results_df(fit, savefile = ''):
    fit_df_GP = fit.extract()['sum_n_predict']
    predicted_GP_mean = np.mean(fit_df_GP, axis = 0)
    predicted_GP_median = np.median(fit_df_GP,axis = 0)
    predicted_q025_GP = np.quantile(fit_df_GP,0.025,axis=0)
    predicted_q975_GP = np.quantile(fit_df_GP,0.975,axis=0)
    predicted_q25_GP = np.quantile(fit_df_GP,0.25,axis=0)
    predicted_q75_GP = np.quantile(fit_df_GP,0.75,axis=0)
    
    nowcasted_data = pd.DataFrame(columns = ['week','mean', 'median', 'q025', 'q975', 'q25', 'q75'])
    nowcasted_data['week'] = range(len(predicted_GP_mean))
    nowcasted_data['mean'] = predicted_GP_mean.astype(int)
    nowcasted_data['median'] = predicted_GP_median.astype(int)
    nowcasted_data['q025'] = predicted_q025_GP.flatten().astype(int)
    nowcasted_data['q975'] = predicted_q975_GP.flatten().astype(int)
    nowcasted_data['q25'] = predicted_q25_GP.flatten().astype(int)
    nowcasted_data['q75'] = predicted_q75_GP.flatten().astype(int)

    if savefile:
        nowcasted_data.to_csv(savefile + '.csv', index=False)
    return nowcasted_data

def save_all_results(fits, savefile_start):
    for k in fits.keys():
        tmp = get_results_df(fits[k], savefile = savefile_start + '_' + k)
        
def get_all_results(fits):
    results = {}
    for k in fits.keys():
        tmp = get_results_df(fits[k], savefile = '') # do not save to csv
        results.update({k: tmp})
    return results

##############################################################################
# CALCULATING ERRORS
##############################################################################

def nowcasting_error(deaths_true, deaths_nowcasted, weeks_test, method):
    """method: 'MSE', 'MAE', 'RMSE'
    weeks_test: how many last weeks use for testing"""
    n = len(deaths_nowcasted)
    deaths_true = deaths_true[:n] # doesn't matter what is later
    deaths_true = deaths_true[-weeks_test:] # we only test last few weeks
    deaths_nowcasted = deaths_nowcasted[-weeks_test:]
    
    if method == 'MSE':
        return sklearn.metrics.mean_squared_error(y_true = deaths_true, y_pred = deaths_nowcasted, squared = True)
    elif method == 'RMSE':
        return sklearn.metrics.mean_squared_error(y_true = deaths_true, y_pred = deaths_nowcasted, squared = False)
    elif method == 'MAE':
        return sklearn.metrics.mean_absolute_error(y_true = deaths_true, y_pred = deaths_nowcasted)
    else:
        print('Wrong method')
        return -1000
    
def nowcasting_error_weighted(deaths_true, deaths_nowcasted, deaths_raw,
                              weeks_test, method):
    """method: 'MSE', 'MAE', 'RMSE'
    weeks_test: how many last weeks use for testing"""
    n = len(deaths_nowcasted)
    deaths_true = deaths_true[:n] # doesn't matter what is later
    deaths_true = deaths_true[-weeks_test:] # we only test last few weeks
    deaths_nowcasted = deaths_nowcasted[-weeks_test:]
    deaths_raw = deaths_raw[-weeks_test:]
    
    # calculate the weights
    weights = abs(deaths_true - deaths_raw)

    if method == 'MSE':
        return sklearn.metrics.mean_squared_error(y_true = deaths_true, 
                                                  y_pred = deaths_nowcasted, 
                                                  sample_weight = weights,
                                                  squared = True)
    elif method == 'RMSE':
        return sklearn.metrics.mean_squared_error(y_true = deaths_true, 
                                                  y_pred = deaths_nowcasted,
                                                  sample_weight = weights,
                                                  squared = False)
    elif method == 'MAE':
        return sklearn.metrics.mean_absolute_error(y_true = deaths_true,
                                                   y_pred = deaths_nowcasted,
                                                   sample_weight = weights)

    else:
        print('Wrong method')
        return -1000