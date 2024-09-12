%% VAR-IMPLIED FORECASTS FOR HISTORICAL EVOLUTION
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/tomyc/Dropbox (MIT)/mp_modelcnfctls/code/github_public/varplus';
vintage = '';
task = '/var_inputs';

addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage task '/_data/main'])

save_results = 1;

cd([path vintage task]);

%% DATA

% import data

data_table = readtable('_data_cmw.csv',detectImportOptions('_data_cmw.csv'));
data = table2array(data_table);
date = data(:,1);

series = data_table.Properties.VariableNames;
series = series(2:end);

% raw macro outcomes

gdp       = data(:,2);
unemp     = data(:,3);
ffr       = 4 * data(:,4);
infl      = 400 * data(:,5);
inv       = data(:,6);
cons      = data(:,7);
lab       = data(:,8);
lab_share = data(:,9);
lab_prod  = data(:,10);
tfp       = data(:,11);

% macro outcome transformations

gdp       = 100 * stat_transform(gdp,1);
inv       = 100 * stat_transform(inv,1);
cons      = 100 * stat_transform(cons,1);
lab       = 100 * stat_transform(lab,1);
lab_share = 100 * stat_transform(lab_share,1);
lab_prod  = 100 * stat_transform(lab_prod,1);
tfp       = 100 * stat_transform(tfp,1);

% dates

startdate = 1960;
enddate   = 2019.75;

startdate_indx = find(date == startdate);
enddate_indx = find(date == enddate);

% collect VAR inputs

vardata = [unemp gdp inv cons lab tfp lab_prod lab_share infl ffr]; 
series_names = {'Unemployment', 'Output', 'Investment', 'Consumption','Hours',...
    'TFP', 'Labor Productivity', 'Labor Share', 'Inflation', 'FFR'};

% positions of the variables of interest

pi_pos = 9;
y_pos  = 2;
i_pos  = 10;

vardata = vardata(startdate_indx:enddate_indx,:);
date    = date(startdate_indx:enddate_indx,:);

% de-trend

vardata_orig = vardata;

det_X = [ones(length(vardata),1),(1:1:length(vardata))'];
[det_coeff,vardata] = ls_detrend(vardata,2);

%% SETTINGS

% VAR specification

n_lags     = 2;                    % number of lags
constant   = 0;                    % constant?
IRF_hor    = 250;
n_draws    = 1000;
n_y        = size(vardata,2);

% historical episode of interest

dates_fcst = [2008.75:0.25:2012];
n_dates    = length(dates_fcst);
vars_fcst  = [pi_pos y_pos i_pos];
n_vars     = length(vars_fcst);

n_dates_fcst = length(dates_fcst);
n_vars_fcst  = length(vars_fcst);

fcst_length = IRF_hor; % this is for computing appropiate dimensions in the code, won't be shown.
fcst_lag    = 10; %how many lags of the data you want to plot.

series_names = series_names(vars_fcst);

%% VAR ESTIMATION

%----------------------------------------------------------------
% Estimate Reduced-Form VAR
%----------------------------------------------------------------

T = size(vardata,1) - n_lags;
[B_draws,Sigma_draws,B_OLS,Sigma_OLS] = bvar_fn(vardata,n_lags,constant,n_draws);

%----------------------------------------------------------------
% OLS Point Estimate Forecasts
%----------------------------------------------------------------

var_forecasts_OLS = NaN(fcst_length,n_vars,n_dates);

for i_date = 1:n_dates_fcst

% get date at which we compute the forecast
fcst_date = find(date == dates_fcst(i_date));

% given that date, use the VAR to compute forecasts
temp = forecast_fn(vardata,n_lags,constant,B_OLS,fcst_date,fcst_length,2); 

% Store the VAR implied forecasts at that date
var_forecasts_OLS(:,:,i_date) = temp(:,vars_fcst);

end

% this is for plotting purposes only.

var_history = vardata(:,vars_fcst);

pi_history = var_history(:,1);
y_history  = var_history(:,2);
i_history  = var_history(:,3);

%----------------------------------------------------------------
% Posterior Forecast Distribution
%----------------------------------------------------------------

% not used in the paper.

var_forecasts_draws = NaN(fcst_length,n_vars,n_draws,n_dates);

for i_date = 1:n_dates_fcst

var_forecasts_all = NaN(fcst_length,n_y,n_draws);

fcst_date = find(date == dates_fcst(i_date));

for i_draw = 1:n_draws

    B = B_draws(:,:,i_draw);

    var_forecasts_all(:,:,i_draw) = forecast_fn(vardata,n_lags,constant,B,fcst_date,fcst_length,2);

end

var_forecasts_all = var_forecasts_all(:,vars_fcst,:);
var_forecasts_draws(:,:,:,i_date) = var_forecasts_all;

end

var_forecasts_lb  = squeeze(quantile(var_forecasts_draws,0.16,3));
var_forecasts_med = squeeze(quantile(var_forecasts_draws,0.5,3));
var_forecasts_ub  = squeeze(quantile(var_forecasts_draws,0.84,3));

%% SAVE RESULTS

if save_results == 1

    cd([path vintage task '/_results']);
    
    save fcst_evol_results var_forecasts_OLS var_history det_coeff det_X dates_fcst fcst_lag date series_names
    
    cd([path vintage task]);

end