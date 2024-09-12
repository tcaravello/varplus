%% DATA

% import data

data_table = readtable([path vintage task  '/_data/main/_data_cmw.csv'],detectImportOptions([path vintage task  '/_data/main/_data_cmw.csv']));
data = table2array(data_table);
date = data(:,1);

series = data_table.Properties.VariableNames;
series = series(2:end);

% raw macro outcomes

gdp       = data(:,2);
unemp     = data(:,3);
ffr       = data(:,4);
infl      = 400 * data(:,5);
inv       = data(:,6);
cons      = data(:,7);
lab       = data(:,8);
lab_share = data(:,9);
lab_prod  = data(:,10);
tfp       = data(:,11);

% macro outcome transformations
[gdp      , gdp_trend]       = regcyc(gdp);
[inv      , inv_trend]       = regcyc(inv);
[cons     , cons_trend]      = regcyc(cons);
[lab      , lab_trend]       = regcyc(lab);
[lab_share, lab_share_trend] = regcyc(lab_share);
[lab_prod , lab_prod_trend]  = regcyc(lab_prod);
[tfp      , tfp_trend]       = regcyc(tfp);

%% COLLECT VAR INPUTS

vardata = [unemp gdp inv cons lab tfp lab_prod lab_share infl ffr]; 
n_y     = size(vardata,2);

% trim the data
startdate_indx = find(date == startdate-2);
enddate_indx = find(date == enddate);

vardata = vardata(startdate_indx:enddate_indx,:);
date = date(startdate_indx:enddate_indx,:);

series_names = {'Unemployment', 'Output', 'Investment', 'Consumption','Hours',...
    'TFP', 'Labor Productivity', 'Labor Share', 'Inflation', 'ffr'};

% de-trend
[det_coeff,vardata] = ls_detrend(vardata,2);

%% DATA FORMATTING

n_var = n_y;                   % number of endogenous variables
m     = n_y*n_lags + constant; % number of exogenous variables

yt = vardata;
T  = size(yt,1);
xt = zeros(T,n_y*n_lags+constant);
for i=1:n_lags
    xt(:,n_y*(i-1)+1:n_y*i) = [NaN(i+HORIZON-1,n_y);vardata(1:end-i-(HORIZON-1),:)];
end
if constant==1
    xt(:,n_y*n_lags+1)=ones(T,1);
elseif constant == 2
    xt(:,n_y*n_lags+1)=ones(T,1);
    xt(:,n_y*n_lags+2)=[1:1:T]';
elseif constant == 3
    xt(:,n_y*n_lags+1)=ones(T,1);
    xt(:,n_y*n_lags+2)=[1:1:T]';
    xt(:,n_y*n_lags+3)=([1:1:T].^2)';
end

startdate_indx = find(date == startdate);
enddate_indx = find(date == enddate);

Y = yt(startdate_indx:enddate_indx,:); % T by nvar matrix of observations
X = xt(startdate_indx:enddate_indx,:); % T by (nvar*nlag+1) matrix of regressors

%% LOAD SW FACTORS

outcomes_of_interest = [2 9 10];

swf = readmatrix([path vintage task  '/_data/stockwatson2016/SW_factors.csv']);
sw_dates = swf(:,1);
sw_date_start = find(sw_dates == startdate -0.25*HORIZON ); % align timing with lag=HORIZON in xt 
sw_date_end =   find(sw_dates == enddate   -0.25*HORIZON );
sw_fac = swf(sw_date_start:sw_date_end   ,  2:end); 

%% COMPUTE R^2

R2 = zeros(n_y,2);
for j = 1:n_y
    rss1 = sum( (Y(:,j) - X * (X \ Y(:,j))).^2  );
    rss2 = sum( (Y(:,j) - [X sw_fac] * ([X sw_fac] \ Y(:,j))).^2  );
    tss = sum((Y(:,j)).^2);
    R2(j,:) = 1-[rss1 rss2]/tss;
end

R2Table(:,:,HORIZON) =  R2(outcomes_of_interest,:);