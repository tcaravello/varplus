%% SAMPLE FROM POSTERIOR FOR ANY MODEL
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:declareGlobalBeforeUse')

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

model      = '/hank'; % '/rank' for RANK, '/hank' for HANK
behavioral = 0; % 1 if you want to estimate the behavioral parameters, 0 if you wanna do ratex. 

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

save_results = 0; % 1 if you want to save results

load_hank = 0; % 0 if you want to compute HANK jacobians from scratch, 1 if you want to load them. 

load_init_state = 0; % 1 if you want to load the initial state (parameters and variance matrix)

global T T_use n_shock

n_shock = 25;       % Number of news shocks to use
T_use = 25;         % Horizon to match IRFs
T = T_use + 275;    % Truncation horizon for model solution
T_save = 200;        % Horizon used to save IRF matrices.

cov_mat = 0; %1 if you want to use diagonal, 0 if you want to use non-diagonal.

% set a low number of draws for exposition purposes. For this to work in
% practise you need a larger number of N_adapt and N_burn. 

N_adapt = 1000; % number of draws to estimate draw covariance matrix. If you wanna turn off set N_adapt = 0
N_burn = N_adapt + 500; % number of draws to burn
N_keep = 100;  % number of matrix draws to keep.
keep_every = 5; %this means that we store 1 out of keep_every draws
N_draws = N_burn + keep_every*N_keep; % total number of draws.
n_iter_print = 100; % print time every n_iter_print iters. Set equal to N_draws to turn off.

var_scale = 0.25;
var_scale_intial_draws = 0.6;

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

%% SAMPLE FROM POSTERIOR

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

settings_MH = struct();
settings_MH.N_adapt = N_adapt;
settings_MH.N_burn = N_burn;
settings_MH.N_keep = N_keep;
settings_MH.keep_every = keep_every;
settings_MH.N_draws = N_draws;
settings_MH.n_iter_print = n_iter_print;
settings_MH.T_use = T_use;
settings_MH.T_save = T_save;
settings_MH.start_chol = load_init_state;
 
if load_init_state == 1

    % load states from a previous run in order to avoid the adaptive and burning steps
    if strcmp(model, '/rank')
        load rank_MH_init
    elseif strcmp(model, '/hank')
        load hank_MH_init
    end
    % use user-provided intial values.
    init_val = init_val_use;
    var_mat_draws_init = chol_use;
else

if behavioral == 1

    % initial values are set at the posterior mode, but other would also
    % work. Initial var-cov matrix is diagonal with the prior sds as
    % suggested Fernandez-Villaverde (2016), sec 12.2.2

    if strcmp(model, '/rank')
        init_val = [0.7527    0.9108   0.99999    0.9054   0.99999    5.5372    0.4881    m_d_mean    m_f_mean];
        var_mat_draws_init = diag([h_sd, calvo_p_sd , 0, calvo_w_sd , 0, kappa_sd, psi_uti_sd, 0, m_f_sd ])*var_scale_intial_draws/var_scale;
    elseif strcmp(model, '/hank')
        init_val = [0.9341    0.9153    0.9999    0.9077    0.9999    5.1929    0.4652    m_d_mean    m_f_mean];
        var_mat_draws_init = diag([0.1, calvo_p_sd , 0, calvo_w_sd , 0, kappa_sd, psi_uti_sd, 0, m_f_sd ])*var_scale_intial_draws/var_scale;
    end

else

    if strcmp(model, '/rank')
        init_val = [0.7527    0.9108   0.99999    0.9054   0.99999    5.5372    0.4881];
        var_mat_draws_init = diag([h_sd, calvo_p_sd , 0, calvo_w_sd , 0, kappa_sd, psi_uti_sd])*var_scale_intial_draws/var_scale;
    elseif strcmp(model, '/hank')
        init_val = [0.9341    0.94   0.99999    0.88    0.99999    5.1929    0.4652];
        var_mat_draws_init = diag([0.1, calvo_p_sd , 0, calvo_w_sd, 0, kappa_sd, psi_uti_sd])*var_scale_intial_draws/var_scale;
    end

end

end

settings_MH.init_val = init_val;
settings_MH.var_scale = var_scale;
settings_MH.var_mat_draws_init = var_mat_draws_init;

%----------------------------------------------------------------
% MH Chain
%----------------------------------------------------------------

[Pi_m_collector, Y_m_collector,R_n_m_collector, m_fit_collector,...
   param_collector, log_posterior_collector,acceptance_collector, ...
   Pi_m_fit_collector, Y_m_fit_collector,R_n_m_fit_collector,chol_use] = run_MH_full_chain(settings_MH,init_val,var_scale,var_mat_draws_init,target,target_Sigma_inv,model);

settings_MH.chol_use = chol_use;
acc_rate = mean(acceptance_collector);

% report mean acceptance:

fprintf('Mean acceptance is %1.4f \n',mean(acceptance_collector));
fprintf('Pre-adj is %1.4f \n',mean(acceptance_collector(1:N_adapt)));
fprintf('Post-adj is %1.4f \n',mean(acceptance_collector(N_adapt+1:end)));

%% SAVE RESULTS

if save_results == 1

cd([path vintage task '/_results'])  

if behavioral == 1
     if strcmp(model, '/rank')
    
    save rank_draws_main_behav Pi_m_collector Y_m_collector R_n_m_collector m_fit_collector
    save rank_draws_other_behav param_collector log_posterior_collector m_fit_collector Pi_m_fit_collector Y_m_fit_collector R_n_m_fit_collector...
         acceptance_collector settings_MH
    
    elseif strcmp(model, '/hank')
    save hank_draws_main_behav Pi_m_collector Y_m_collector R_n_m_collector m_fit_collector
    save hank_draws_other_behav param_collector log_posterior_collector m_fit_collector Pi_m_fit_collector Y_m_fit_collector R_n_m_fit_collector...
         acceptance_collector settings_MH
    
     end

else

    if strcmp(model, '/rank')
    
    save rank_draws_main Pi_m_collector Y_m_collector R_n_m_collector m_fit_collector
    save rank_draws_other param_collector log_posterior_collector m_fit_collector Pi_m_fit_collector Y_m_fit_collector R_n_m_fit_collector...
         acceptance_collector settings_MH
    
    elseif strcmp(model, '/hank')
    save hank_draws_main Pi_m_collector Y_m_collector R_n_m_collector m_fit_collector
    save hank_draws_other param_collector log_posterior_collector m_fit_collector Pi_m_fit_collector Y_m_fit_collector R_n_m_fit_collector...
         acceptance_collector settings_MH
    
    end

end

cd([path vintage task])

end