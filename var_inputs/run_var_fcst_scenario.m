%% VAR-IMPLIED FORECASTS FOR SCENARIO ANALYSIS
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
enddate   = 2021.25;

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

fcst_date  = find(date == enddate); % date at which you are making the forecast, here  I start at 2021Q2
fcst_hor   = 20; % how many quarters out you want to forecast for plotting
vars_fcst  = [pi_pos y_pos i_pos]; % ordering of the variables to forecast
n_vars     = length(vars_fcst);

fcst_length = IRF_hor; % this is for computing appropiate dimensions in the code, won't be shown
fcst_lag    = 10; % how many lags of the data you want to plot

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

var_forecasts_OLS = forecast_fn(vardata,n_lags,constant,B_OLS,fcst_date,fcst_length,1);
var_forecasts_OLS = var_forecasts_OLS(:,vars_fcst);

% this will be used for plotting only.

var_history = vardata(:,vars_fcst);

pi_history = var_history(:,1);
y_history  = var_history(:,2);
i_history  = var_history(:,3);

%----------------------------------------------------------------
% Posterior Forecast Distribution
%----------------------------------------------------------------

% not used in the paper

var_forecasts_draws = NaN(fcst_length,n_y,n_draws);

for i_draw = 1:n_draws

    B = B_draws(:,:,i_draw);
    var_forecasts_draws(:,:,i_draw) = forecast_fn(vardata,n_lags,constant,B,fcst_date,fcst_length,1);

end

var_forecasts_draws = var_forecasts_draws(:,vars_fcst,:);

var_forecasts_lb  = squeeze(quantile(var_forecasts_draws,0.16,3));
var_forecasts_med = squeeze(quantile(var_forecasts_draws,0.5,3));
var_forecasts_ub  = squeeze(quantile(var_forecasts_draws,0.84,3));

%% COLLECT RESULTS

%----------------------------------------------------------------
% Add History to Forecasts
%----------------------------------------------------------------

% extend det_X

det_trend = [(1:1:length(vardata))'; repmat(length(vardata),fcst_hor,1)];
det_X_ext = [ones(length(vardata)+fcst_hor,1),det_trend];

% inflation

pi_forecasts_lb  = [pi_history(fcst_date-fcst_lag:fcst_date);var_forecasts_lb(1:fcst_hor,1,1)];
pi_forecasts_med = [pi_history(fcst_date-fcst_lag:fcst_date);var_forecasts_med(1:fcst_hor,1,1)];
pi_forecasts_ub  = [pi_history(fcst_date-fcst_lag:fcst_date);var_forecasts_ub(1:fcst_hor,1,1)];

pi_history = pi_history + det_X * det_coeff(:,pi_pos);
pi_history = [pi_history;NaN(fcst_hor,1)];

pi_forecasts_lb  = pi_forecasts_lb + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);
pi_forecasts_med = pi_forecasts_med + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);
pi_forecasts_ub  = pi_forecasts_ub + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);

% output

y_forecasts_lb  = [y_history(fcst_date-fcst_lag:fcst_date);var_forecasts_lb(1:fcst_hor,2,1)];
y_forecasts_med = [y_history(fcst_date-fcst_lag:fcst_date);var_forecasts_med(1:fcst_hor,2,1)];
y_forecasts_ub  = [y_history(fcst_date-fcst_lag:fcst_date);var_forecasts_ub(1:fcst_hor,2,1)];

y_history = y_history + det_X * det_coeff(:,y_pos);
y_history = [y_history;NaN(fcst_hor,1)];

y_forecasts_lb  = y_forecasts_lb + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);
y_forecasts_med = y_forecasts_med + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);
y_forecasts_ub  = y_forecasts_ub + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);

% interest rates

i_forecasts_lb  = [i_history(fcst_date-fcst_lag:fcst_date);var_forecasts_lb(1:fcst_hor,3,1)];
i_forecasts_med = [i_history(fcst_date-fcst_lag:fcst_date);var_forecasts_med(1:fcst_hor,3,1)];
i_forecasts_ub  = [i_history(fcst_date-fcst_lag:fcst_date);var_forecasts_ub(1:fcst_hor,3,1)];

i_history = i_history + det_X * det_coeff(:,i_pos);
i_history = [i_history;NaN(fcst_hor,1)];

i_forecasts_lb  = i_forecasts_lb + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);
i_forecasts_med = i_forecasts_med + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);
i_forecasts_ub  = i_forecasts_ub + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);

% extend date

date_ext = [date(end):0.25:date(end)+0.25*fcst_hor]';
date_ext = [date(1:end-1);date_ext];

%----------------------------------------------------------------
% Save Results
%----------------------------------------------------------------

if save_results == 1

    cd([path vintage task '/_results']);
    
    save fcst_scenario_results var_forecasts_OLS var_history det_coeff det_X det_X_ext ... 
        fcst_lag fcst_hor fcst_date date date_ext series_names pi_pos y_pos i_pos
    
    cd([path vintage task]);

end