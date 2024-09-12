%% VAR-IMPLIED FORECASTS FOR COMPARISON AGAINST SPF
% Tomas Caravello, Alisdair McKay, Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

task = '/var_inputs';

%% DATA

% import data

data_table = readtable([path vintage task  '/_data/spf/_data_spf_comparison.csv']);
data = table2array(data_table);
date = data(:,1);

series = data_table.Properties.VariableNames;
series = series(2:end);

% raw macro outcomes

gdp       = data(:,2);
unemp     = data(:,3);
tbill     = data(:,4);
infl      = 400 * data(:,5);
inv       = data(:,6);
cons      = data(:,7);
lab       = data(:,8);
lab_share = data(:,9);
lab_prod  = data(:,10);
tfp       = data(:,11);

data_orig = [infl gdp tbill];

% macro outcome transformations

[gdp      , gdp_trend]       = regcyc(gdp);
[inv      , inv_trend]       = regcyc(inv);
[cons     , cons_trend]      = regcyc(cons);
[lab      , lab_trend]       = regcyc(lab);
[lab_share, lab_share_trend] = regcyc(lab_share);
[lab_prod , lab_prod_trend]  = regcyc(lab_prod);
[tfp      , tfp_trend]       = regcyc(tfp);

% dates

startdate = 1960;
enddate   = 2019.75;

startdate_indx = find(date == startdate);
enddate_indx = find(date == enddate);

% collect VAR inputs

vardata = [unemp gdp inv cons lab tfp lab_prod lab_share infl tbill]; 
series_names = {'Unemployment', 'Output', 'Investment', 'Consumption','Hours',...
    'TFP', 'Labor Productivity', 'Labor Share', 'Inflation', 'TBill'};

trends = zeros(size(vardata));
trends(:,2) = gdp_trend;
trends(:,3) = inv_trend;
trends(:,4) = cons_trend;
trends(:,5) = lab_trend;
trends(:,6) = tfp_trend;
trends(:,7) = lab_prod_trend;
trends(:,8) = lab_share_trend;

vardata = vardata(startdate_indx:enddate_indx,:);
date    = date(startdate_indx:enddate_indx,:);
trends  = trends(startdate_indx:enddate_indx,:);
T       = size(vardata,1);
data_orig = data_orig(startdate_indx:enddate_indx,:);

% de-trend

vardata_orig = vardata;

det_X = [ones(length(vardata),1),(1:1:length(vardata))'];
[det_coeff,vardata] = ls_detrend(vardata,2);

trends = trends + [ones(T,1),(1:1:T)'] * det_coeff;

%% SETTINGS

% VAR specification

n_lags     = 2;                    % number of lags
constant   = 0;                    % constant?
IRF_hor    = 20;
n_draws    = 20;
n_y        = size(vardata,2);

% historical episode of interest

dates_fcst = [1968.75:0.25:2019.75];
n_dates    = length(dates_fcst);
vars_fcst  = [9 2 10]; % inflation, GDP, interest rate
n_vars     = length(vars_fcst);

n_dates_fcst = length(dates_fcst);
n_vars_fcst  = length(vars_fcst);

fcst_length = IRF_hor;
fcst_lag    = 10;

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

    fcst_date = find(date == dates_fcst(i_date));
    
    temp = forecast_fn(vardata,n_lags,constant,B_OLS,fcst_date,fcst_length,2);
    var_forecasts_OLS(:,:,i_date) = temp(:,vars_fcst);

end

%----------------------------------------------------------------
% Further Transformations
%----------------------------------------------------------------

% add back trends

for i_date = 1:n_dates_fcst
    for j = 1:n_vars
        for h = 1:fcst_length
            fcst_date = find(date == dates_fcst(i_date));
            fdate = fcst_date+h-1;
            if fdate > size(trends,1)
                tr = nan;
            else
                tr = trends(fdate,vars_fcst(j));
            end
            var_forecasts_OLS(h,j,i_date) = var_forecasts_OLS(h,j,i_date) + tr;
        end
    end
end

% express real gdp as cumulative log change from current period

j = 2;
assert(strcmp(series_names{j},'Output'))
for i_date = 1:n_dates_fcst
    for h = 2:fcst_length
        var_forecasts_OLS(h,j,i_date) = 100*(var_forecasts_OLS(h,j,i_date) - var_forecasts_OLS(1,j,i_date));
    end
end

% convert to 2-d array shape T x 12
% [inflation 1-4 quarters ahead,   gdp 1-4 qtrs,   tbill  1-4 qtrs]
vararr = zeros(n_dates_fcst,4*n_vars_fcst);
for i_date = 1:n_dates_fcst
    for j = 1:n_vars
        for h = 2:5
            vararr(i_date,h-1 + (j-1)*4) = var_forecasts_OLS(h,j,i_date);
        end
    end
end

%% LOAD SPF

spfpgdp = readtable([path vintage task '/_data/spf/' 'spf_meanGrowth'],sheet="PGDP");
spfrgdp = readtable([path vintage task '/_data/spf/' 'spf_meanLevel'],sheet="RGDP");
spftbil = readtable([path vintage task  '/_data/spf/' 'spf_meanLevel'],sheet="TBILL");

assert(all(spfpgdp.YEAR == spftbil.YEAR))
assert(all(spfpgdp.QUARTER == spftbil.QUARTER))
assert(all(spfpgdp.YEAR == spfrgdp.YEAR))
assert(all(spfpgdp.QUARTER == spfrgdp.QUARTER))

spf = [spfpgdp removevars(spfrgdp,{'YEAR','QUARTER'}) removevars(spftbil,{'YEAR','QUARTER'})];
spf.date = spf.YEAR + 0.25*(spf.QUARTER-1);

% express real gdp as cumulative growth from current quarter
spf.RGDP3 = 100*(log(spf.RGDP3) - log(spf.RGDP2));
spf.RGDP4 = 100*(log(spf.RGDP4) - log(spf.RGDP2));
spf.RGDP5 = 100*(log(spf.RGDP5) - log(spf.RGDP2));
spf.RGDP6 = 100*(log(spf.RGDP6) - log(spf.RGDP2));

% convert to array
spfarr = [table2array(spf(ismember (spf.date , dates_fcst),{'dpgdp3','dpgdp4','dpgdp5','dpgdp6',...
                    'RGDP3','RGDP4','RGDP5','RGDP6',...
                    'TBILL3','TBILL4','TBILL5','TBILL6'}))];

%% CREATE ARRAY OF ACTUALS

actualarr = zeros(n_dates_fcst,12);
for i_date = 1:n_dates_fcst  
    for h = 2:5
        fcst_date = find(date == dates_fcst(i_date));
        fdate = fcst_date+h-1;
        if fdate > size(data_orig,1)
            ac = NaN(1,3);
        else
            ac = data_orig(fdate,:);
            ac(2) = 100*(ac(2) - data_orig(fcst_date,2)); % gdp is cumulative log change
        end
        for j = 1:n_vars
            actualarr(i_date,h-1+(j-1)*4) = ac(j);
        end
    end
end


%% COMPUTE MSEs

assert(all(size(actualarr) == size(vararr)))
assert(all(size(actualarr) == size(spfarr)))

start_comp = 1981.5;
end_comp   = 2007.5;
dates_comp = start_comp:0.25:end_comp;

vararr = vararr(ismember(dates_fcst, dates_comp),:);
spfarr = spfarr(ismember(dates_fcst, dates_comp),:);
actualarr = actualarr(ismember(dates_fcst, dates_comp),:);

varmse = mean((actualarr - vararr).^2 , 1, "omitnan");
spfmse = mean((actualarr - spfarr).^2 , 1, "omitnan");

rownames ={};
for s = {'infl ','GDP  ','Tbill'}
    for h = 1:4
        rownames = [rownames,[s{1} ' ' num2str(h) ' qtr ahead']];
    end
end

mse = array2table([varmse' spfmse'],...
    'VariableNames',{'VAR','SPF'},'RowNames',rownames);

%% FINAL TABLE

cd([path vintage task])

f = fopen(fullfile('_results', strcat('spf_compare', '.tex')), 'w'); % open file for writing

fprintf(f,"\\begin{table}\n");
fprintf(f,"\\centering\n");
fprintf(f,"\\begin{tabular}{lccccc}\n");
fprintf(f,"\\toprule \\toprule");
fprintf(f,"\\textbf{Variable} & \\multicolumn{2}{c}{\\textbf{1-quarter ahead} } && \\multicolumn{2}{c}{ \\textbf{4-quarter ahead}} \\\\ \n");
fprintf(f,"\\cline{2-3}  \\cline{5-6} \n");
fprintf(f,"& \\textbf{VAR} & \\textbf{SPF}  && \\textbf{VAR} & \\textbf{SPF} \\\\ \n");
fprintf(f,strcat("\\midrule GDP & ", num2str(varmse(5),3),"  & ",num2str(spfmse(5),3)," && ",num2str(varmse(8),3),"  & ",num2str(spfmse(8),3)," \\\\ \n"));
fprintf(f,strcat("Inflation & ", num2str(varmse(1),3),"  & ",num2str(spfmse(1),3)," && ",num2str(varmse(4),3),"  & ",num2str(spfmse(4),3)," \\\\ \n"));
fprintf(f,strcat("T-Bill rate & ", num2str(varmse(9),3),"  & ",num2str(spfmse(9),3)," && ",num2str(varmse(12),3),"  & ",num2str(spfmse(12),3)," \\\\ \\bottomrule \n"));
fprintf(f,"\\end{tabular} \n");
fprintf(f,"\\caption{Mean squared error of 1- and 4-quarter ahead forecasts.} \n");
fprintf(f,"\\label{tab:spf_compare} \n");
fprintf(f,"\\end{table} \n");

fclose(f);