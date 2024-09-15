%% FIND POSTERIOR MODE AND GET IRFS FOR ANY MODEL
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/tomyc/Dropbox (MIT)/mp_modelcnfctls/code/github_public/varplus';
vintage = '';
task = '/model_estim';

addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage task '/_results'])
addpath([path vintage '/var_inputs/_results'])

cd([path vintage task]);

%% SET-UP

%----------------------------------------------------------------
% Model
%----------------------------------------------------------------

model      = '/rank'; %type '/rank' for RANK and'/hank' for HANK
behavioral = 0; % 1 if you want to estimate the behavioral parameters, 0 if you want to do ratex

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

global T T_use n_shock

save_results = 0; % 1 if you want to save the results

load_hank = 1; % only applies if model = hank. 1: load Jacobians from previous run
               % 0: compute steady state and jacobians from scratch, takes a few seconds

cov_mat = 0; % 1 if you want to use diagonal, 0 if you want to use non-diagonal
             % in the paper we use non-diagonal because it better represents the informativeness of the data to do model selection

n_shock = 25;       % number of news shocks to use
T_use = 25;         % horizon to match IRFs
T = T_use + 275;    % truncation horizon for model solution
T_save = 200;       % horizon used to save IRF matrices

%% IMPORTS AND PREPARATIONS

%----------------------------------------------------------------
% Empirical Targets
%----------------------------------------------------------------

get_empirical_targets;

%----------------------------------------------------------------
% Calibration
%----------------------------------------------------------------

get_calibration_general;

%----------------------------------------------------------------
% Jacobians
%----------------------------------------------------------------

get_baseline_jacobians;

%----------------------------------------------------------------
% Priors
%----------------------------------------------------------------

get_priors_general;

%% OBTAIN POSTERIOR MODE

%----------------------------------------------------------------
% Find Mode
%----------------------------------------------------------------

if behavioral == 1

if strcmp(model, '/rank')
     x0 = [0.71    0.86    0.99999    0.946    0.99999    5.3    0.47    m_d_mean     m_f_mean];
elseif strcmp(model, '/hank')
     x0 = [0.9598    0.85    0.99999    0.95     0.99999    5.3219    0.47    m_d_mean    m_f_mean];
end

lower_bound = [0,0,0,0,0,0,0,0,0]; 
upper_bound = [1,1,1,1,1,1,1,1,1]-(1e-13);
upper_bound(6) = inf;

else
    if strcmp(model, '/rank')
        x0 = [0.73    0.85   0.99999    0.74   0.99999    5.3    0.47];
    elseif strcmp(model, '/hank')
        x0 = [0.9654    0.9530    0.999999    0.8763    0.999999    5.7686    0.4618];
    end
lower_bound = [0,0,0,0,0,0,0]; 
upper_bound = [1,1,1,1,1,1,1]-(1e-13);
upper_bound(6) = inf;

end

% use low levels of tolerance for illustration purpose only.

options_min = optimset('Display','iter','MaxFunEvals',1000,'tolFun',1e-4, 'TolX', 1e-4);

% get posterior mode
[param_sol,posterior_mode,exit_flag,output] = fminsearchcon(@(x) solve_model(x,target,target_Sigma_inv,model),x0, lower_bound,upper_bound, [],[],[],options_min);

%----------------------------------------------------------------
% Collect Parameters
%----------------------------------------------------------------

kappa_p = (1-beta*param_sol(2))*(1-param_sol(2))/param_sol(2);
kappa_w = (1-beta*param_sol(4))*(1-param_sol(4))/param_sol(4);

if length(param_sol)>7
    beta_p = beta*param_sol(2)*param_sol(9)*(1+kappa_p/(1-beta*param_sol(2)*param_sol(9)));
    beta_w = beta*param_sol(4)*param_sol(9)*(1+kappa_w/(1-beta*param_sol(4)*param_sol(9)));
else
    beta_p = beta;
    beta_w = beta;
end

%----------------------------------------------------------------
% Final IRFs
%----------------------------------------------------------------

[error_fit, Pi_m_model, Y_m_model, R_n_m_model, m_fit_model,...
    Pi_m_fit, Y_m_fit, R_n_m_fit, C_m_model, Inv_m_model, C_m_fit, Inv_m_fit] = solve_model(param_sol,target,target_Sigma_inv,model);
 
m_fit_model_2 = [m_fit_model(1:end-1);0]; 
Pi_m_fit_2 = Pi_m_model(1:T_use,1:T_use) * m_fit_model_2;
Y_m_fit_2 = Y_m_model(1:T_use,1:T_use) * m_fit_model_2;
R_n_m_fit_2 = R_n_m_model(1:T_use,1:T_use) * m_fit_model_2;

%% SAVE POSTERIOR MODE

cd([path vintage task '/_results']);

if save_results ==1

if behavioral == 1

    if strcmp(model, '/rank')

        Pi_m_rank = Pi_m_model;
        Y_m_rank = Y_m_model;
        R_n_m_rank = R_n_m_model; 
        m_fit_rank = m_fit_model;
        
        Pi_m_base_rank = Pi_m_rank(1:T_save,1:T_save); 
        Y_m_base_rank = Y_m_rank(1:T_save,1:T_save); 
        I_m_base_rank = R_n_m_rank(1:T_save,1:T_save);
        
        save params_rank_mode_behav T T_use target_Sigma_inv param_sol posterior_mode Pi_m_rank Y_m_rank R_n_m_rank m_fit_rank
        save rank_base_mode_behav Pi_m_base_rank Y_m_base_rank I_m_base_rank

    elseif strcmp(model, '/hank')

        Pi_m_hank = Pi_m_model;
        Y_m_hank = Y_m_model;
        R_n_m_hank = R_n_m_model; 
        m_fit_hank = m_fit_model;
        
        Pi_m_base_hank = Pi_m_hank(1:T_save,1:T_save); 
        Y_m_base_hank = Y_m_hank(1:T_save,1:T_save); 
        I_m_base_hank = R_n_m_hank(1:T_save,1:T_save);
        
        save params_hank_mode_behav T T_use target_Sigma_inv param_sol posterior_mode Pi_m_hank Y_m_hank R_n_m_hank m_fit_hank
        save hank_base_mode_behav Pi_m_base_hank Y_m_base_hank I_m_base_hank

    end

    else

    if strcmp(model, '/rank')

        Pi_m_rank = Pi_m_model;
        Y_m_rank = Y_m_model;
        R_n_m_rank = R_n_m_model; 
        m_fit_rank = m_fit_model;
        
        Pi_m_base_rank = Pi_m_rank(1:T_save,1:T_save); 
        Y_m_base_rank = Y_m_rank(1:T_save,1:T_save); 
        I_m_base_rank = R_n_m_rank(1:T_save,1:T_save);
        
        save params_rank_mode T T_use target_Sigma_inv param_sol posterior_mode Pi_m_rank Y_m_rank R_n_m_rank m_fit_rank
        save rank_base_mode Pi_m_base_rank Y_m_base_rank I_m_base_rank

    elseif strcmp(model, '/hank')

        Pi_m_hank = Pi_m_model;
        Y_m_hank = Y_m_model;
        R_n_m_hank = R_n_m_model; 
        m_fit_hank = m_fit_model;
        
        Pi_m_base_hank = Pi_m_hank(1:T_save,1:T_save); 
        Y_m_base_hank = Y_m_hank(1:T_save,1:T_save); 
        I_m_base_hank = R_n_m_hank(1:T_save,1:T_save);
        
        save params_hank_mode T T_use target_Sigma_inv param_sol posterior_mode Pi_m_hank Y_m_hank R_n_m_hank m_fit_hank
        save hank_base_mode Pi_m_base_hank Y_m_base_hank I_m_base_hank

    end

end

end

cd([path vintage task]);