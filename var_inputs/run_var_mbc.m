%% MBC SHOCK IRFs
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

startdate = find(date == 1960);
enddate = find(date == 2019.75);

% collect VAR inputs

vardata = [unemp gdp inv cons lab tfp lab_prod lab_share infl ffr]; 
series_names = {'Unemployment', 'Output', 'Investment', 'Consumption','Hours',...
    'TFP', 'Labor Productivity', 'Labor Share', 'Inflation', 'FFR'};

vardata = vardata(startdate:enddate,:);

% de-trend

vardata = detrend(vardata);

%% SETTINGS

% VAR specification

n_lags     = 4;                    % number of lags
constant   = 0;                    % constant?
IRF_hor    = 250;
n_draws    = 1000;
n_y        = size(vardata,2);

% MBC shock identification

grdd = 1024;
wmin = 2*pi/32;
wmax = 2*pi/6;

idx_shock = 1;

%% VAR ESTIMATION

%----------------------------------------------------------------
% Estimate Reduced-Form VAR
%----------------------------------------------------------------

T = size(vardata,1) - n_lags;
[B_draws,Sigma_draws,B_OLS,Sigma_OLS] = bvar_fn(vardata,n_lags,constant,n_draws);

%----------------------------------------------------------------
% MBC Shock IRFs
%----------------------------------------------------------------

% auxiliary objects

MYtmp = [eye(n_y) zeros(n_y,n_y*(n_lags-1))];
MY    = MYtmp(1:n_y, :);

% some placeholders

Qmat    = NaN(n_draws,n_y);          % impulse vector
IRFsr   = NaN(n_draws,n_y*IRF_hor);  % IRFs

% draw from posterior

for i_draw = 1:n_draws

Sigma_u   = Sigma_draws(:,:,i_draw);
B         = B_draws(:,:,i_draw);

bench_rot = chol(Sigma_u,'lower');
tmp = B';
DYN = tmp(:, 1:(n_y*n_lags));
MX  = [DYN;eye(n_y*(n_lags - 1)) zeros(n_y*(n_lags - 1), n_y)];
ME  = [bench_rot;zeros(n_y*(n_lags - 1), n_y)];

Itmp1 = NaN(IRF_hor*n_y, n_y); % IRF
Itmp2 = NaN(n_y, n_y);  % IRF at test horizon
Itmp3 = NaN(IRF_hor * n_y, n_y); % Cumulative squared IRF
for kk=1:n_y
    [tmpYY, tmpXX]  = comp_irf(ME, MX, MY, IRF_hor, kk);
    ttmp = tmpYY(1:IRF_hor,:);
    Itmp1(:,kk) = ttmp(:);
    Itmp2(:,kk) = tmpYY(1,:)';
    vtmp = tmpYY(1:IRF_hor, :);
    Itmp3(:,kk) = vtmp(:);
end

[eigvec, eigval] = funcfdq(MX,MY,bench_rot,wmin,wmax,idx_shock,1,grdd); % extract shock
tst = Itmp2(idx_shock,:)*eigvec;
sgn = sign(tst(1));
Qmat(i_draw,:) =  sgn* eigvec;
IRFsr(i_draw,:) = sgn*(Itmp1 * eigvec);

end

% collect IRFs

airf = reshape(mean(IRFsr,1),IRF_hor,n_y);
mirf = reshape(quantile(IRFsr,0.5), IRF_hor, n_y);
lirf = reshape(quantile(IRFsr,0.16), IRF_hor, n_y);
sirf = reshape(quantile(IRFsr,0.84), IRF_hor, n_y);

%% SAVE RESULTS

IS_MBC  = mirf;

if save_results == 1

    cd([path vintage task '/_results']);
    
    save mbc_results IS_MBC series_names
    
    cd([path vintage task]);
    
end