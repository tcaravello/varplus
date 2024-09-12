%% MONETARY POLICY SHOCK IRFs: ARUOBA-DRECHSEL SHOCK
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
infl      = 100 * 4 * data(:,5);
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
unemp     = stat_transform(unemp,1);

% monetary shocks

ff4     = data(:,12);
onrun10 = data(:,13);
bm_iv   = ff4 + onrun10;

ir_1 = data(:,14);
ir_2 = data(:,15);
ir_3 = data(:,16);

rr_1 = data(:,17);
rr_2 = data(:,18);
rr_3 = data(:,19);
rr_4 = data(:,20);

gk_var_1 = data(:,21);
gk_var_2 = data(:,22);

mp1_tc      = data(:,23);
ff4_tc      = data(:,24);
ed2_tc      = data(:,25);
ed3_tc      = data(:,26);
ed4_tc      = data(:,27);
ff1_vr      = data(:,28);
ff4_vr      = data(:,29);
ed2_vr      = data(:,30);
ff1_gkgreen = data(:,31);
ff4_gkgreen = data(:,32);
ed2_gkgreen = data(:,33);

mar_raw    = data(:,34);
mar_svariv = data(:,35);
ad_shock   = data(:,36);
lpcom = data(:,37);

% dates

startdate = find(date == 1969);
enddate = find(date == 2006.75);

% collect VAR inputs

vardata = [ad_shock gdp infl ffr];
vardata = vardata(startdate:enddate,:);
disp('Note: I am setting missing values of the IV to 0.')
vardata(isnan(vardata)) = 0;

% series names

series_names = {'AD Shock','Output','Inflation','Interest Rate'};

%% SETTINGS

% VAR specification

n_lags     = 2;                    % number of lags
constant   = 2;                    % constant?
IRF_hor    = 100;
n_draws    = 1000;
n_y        = size(vardata,2);

% identification

IS.shock_pos = [1];
n_shocks     = length(IS.shock_pos);

% outcomes of interest

IS.y_pos = [2 3 4];

% placeholders

IS.IRF     = NaN(IRF_hor,n_y,n_shocks,n_draws);
IS.IRF_med = NaN(IRF_hor,n_y,n_shocks);
IS.IRF_lb  = NaN(IRF_hor,n_y,n_shocks);
IS.IRF_ub  = NaN(IRF_hor,n_y,n_shocks);
IS.IRF_OLS = NaN(IRF_hor,n_y,n_shocks);
IS.IRF_variance = NaN(IRF_hor,n_y,n_shocks);

%% VAR ESTIMATION

%----------------------------------------------------------------
% Estimate Reduced-Form VAR
%----------------------------------------------------------------

T = size(vardata,1) - n_lags;
[B_draws,Sigma_draws,B_OLS,Sigma_OLS] = bvar_fn(vardata,n_lags,constant,n_draws);

%----------------------------------------------------------------
% OLS IRFs
%----------------------------------------------------------------

% extract VAR inputs
    
Sigma_u   = Sigma_OLS;
B         = B_OLS;

% benchmark rotation

bench_rot = chol(Sigma_u,'lower');

% Wold IRFs

IRF_Wold = zeros(n_y,n_y,IRF_hor); % row is variable, column is shock
IRF_Wold(:,:,1) = eye(n_y);

for l = 1:IRF_hor
    
    if l < IRF_hor
        for j=1:min(l,n_lags)
            IRF_Wold(:,:,l+1) = IRF_Wold(:,:,l+1) + B(1+(j-1)*n_y:j*n_y,:)'*IRF_Wold(:,:,l-j+1);
        end
    end
    
end

W = bench_rot;

% get IRFs

IRF_OLS = NaN(n_y,n_y,IRF_hor);
for i_hor = 1:IRF_hor
    IRF_OLS(:,:,i_hor) = IRF_Wold(:,:,i_hor) * W;
end

% collect results

for i_shock = 1:length(IS.shock_pos)
    IS.IRF_OLS(:,:,i_shock) = squeeze(IRF_OLS(:,IS.shock_pos(i_shock),:))';
end

%----------------------------------------------------------------
% Confidence Set
%----------------------------------------------------------------

for i_draw = 1:n_draws
    
% extract VAR inputs
    
Sigma_u   = Sigma_draws(:,:,i_draw);
B         = B_draws(:,:,i_draw);

% benchmark rotation

bench_rot = chol(Sigma_u,'lower');

% Wold IRFs

IRF_Wold = zeros(n_y,n_y,IRF_hor); % row is variable, column is shock
IRF_Wold(:,:,1) = eye(n_y);

for l = 1:IRF_hor
    
    if l < IRF_hor
        for j=1:min(l,n_lags)
            IRF_Wold(:,:,l+1) = IRF_Wold(:,:,l+1) + B(1+(j-1)*n_y:j*n_y,:)'*IRF_Wold(:,:,l-j+1);
        end
    end
    
end

W = bench_rot;

% get IRFs

IRF_draw = NaN(n_y,n_y,IRF_hor);
for i_hor = 1:IRF_hor
    IRF_draw(:,:,i_hor) = IRF_Wold(:,:,i_hor) * W;
end

% collect results

for i_shock = 1:length(IS.shock_pos)
    IS.IRF(:,:,i_shock,i_draw) = squeeze(IRF_draw(:,IS.shock_pos(i_shock),:))';
end

end

for ii=1:IRF_hor
    for jj=1:n_y
        for kk = 1:n_shocks
            IS.IRF_med(ii,jj,kk) = quantile(IS.IRF(ii,jj,kk,:),0.5);
            IS.IRF_lb(ii,jj,kk) = quantile(IS.IRF(ii,jj,kk,:),0.16);
            IS.IRF_ub(ii,jj,kk) = quantile(IS.IRF(ii,jj,kk,:),0.84);
            IS.IRF_variance(ii,jj,kk) = var(IS.IRF(ii,jj,kk,:));
        end
    end
end

Y_m_emp         = IS.IRF_OLS(:,2,1);
Y_m_var_emp     = IS.IRF_variance(:,2,1);
Y_m_lb_emp      = IS.IRF_lb(:,2,1);
Y_m_ub_emp      = IS.IRF_ub(:,2,1);
Y_m_draws_emp   = IS.IRF(:,2,1,:);

Pi_m_emp        = IS.IRF_OLS(:,3,1);
Pi_m_var_emp    = IS.IRF_variance(:,3,1);
Pi_m_lb_emp     = IS.IRF_lb(:,3,1);
Pi_m_ub_emp     = IS.IRF_ub(:,3,1);
Pi_m_draws_emp  = IS.IRF(:,3,1,:);

R_n_m_emp       = IS.IRF_OLS(:,4,1);
R_n_m_var_emp   = IS.IRF_variance(:,4,1);
R_n_m_lb_emp    = IS.IRF_lb(:,4,1);
R_n_m_ub_emp    = IS.IRF_ub(:,4,1);
R_n_m_draws_emp = IS.IRF(:,4,1,:);

% save results

if save_results == 1

    cd([path vintage task '/_results'])
    
    save IRFs_ad_results Y_m_emp Y_m_var_emp Y_m_lb_emp Y_m_ub_emp Y_m_draws_emp...
            Pi_m_emp Pi_m_var_emp Pi_m_lb_emp Pi_m_ub_emp Pi_m_draws_emp...
            R_n_m_emp R_n_m_var_emp R_n_m_lb_emp R_n_m_ub_emp R_n_m_draws_emp
    
    cd([path vintage task]);

end