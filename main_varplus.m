%% REPLICATION FILES: "EVALUATING POLICY COUNTERFACTUALS: A VAR-PLUS APPROACH"
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:declareGlobalBeforeUse')

path = '/Users/tomyc/Dropbox (MIT)/mp_modelcnfctls/code/github_public/varplus';
%path = '/Users/christianwolf/Dropbox/Research/mp_modelcnfctls/code/replication';
% path = '/Users/iragm01/Library/CloudStorage/Dropbox/mp_modelcnfctls/code/replication';

vintage = '';

addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage '/var_inputs'])
addpath([path vintage '/var_inputs/_results'])
addpath([path vintage '/model_estim'])
addpath([path vintage '/applications/second_moments'])
addpath([path vintage '/applications/hist_scenario'])
addpath([path vintage '/applications/hist_evol'])
addpath([path vintage '/invertibility/sw_sequence'])

%% THEORETICAL INVERTIBILITY ANALYSIS (FIGURE 1)

% counterfactuals

get_2mom_ninvwold_all;

% forecasting

get_2mom_fcstvars;

%% APPLICATIONS

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

% sample

indic_early  = 0; % early sample?

% target shocks and extrapolation

indic_1shock = 0; % use only empirical target?

indic_RE    = 0; % only RE models?
indic_behav = 0; % only behavioral models?
indic_joint = 1; % all models?

% target counterfactual rule

cnfctl_0y       = 0; % output gap targeting
cnfctl_0pi      = 0; % inflation targeting
cnfctl_0ib      = 0; % nominal rate peg
cnfctl_tylr     = 0; % Taylor rule
cnfctl_ngdp     = 0; % NGDP targeting
cnfctl_ibtarget = 0; % rate target
cnfctl_optpol   = 1; % dual mandate

%----------------------------------------------------------------
% Figures 2-4
%----------------------------------------------------------------

plot_model_irfs;

%----------------------------------------------------------------
% Figures 5-6
%----------------------------------------------------------------

get_cnfctl_stats;
get_cnfctl_mbc;

%----------------------------------------------------------------
% Figure 7
%----------------------------------------------------------------

get_historical_evol;

%----------------------------------------------------------------
% Figure 8
%----------------------------------------------------------------

% RE

indic_RE    = 1; % only RE models?
indic_behav = 0; % only behavioral models?
indic_joint = 0; % all models?

get_historical_scenario;

% behavioral

indic_RE    = 0; % only RE models?
indic_behav = 1; % only behavioral models?
indic_joint = 0; % all models?

get_historical_scenario;

% joint

indic_RE    = 0; % only RE models?
indic_behav = 0; % only behavioral models?
indic_joint = 1; % all models?

get_historical_scenario;

%----------------------------------------------------------------
% Figure 9
%----------------------------------------------------------------

decompose_realrates_brank;

%----------------------------------------------------------------
% Figures D.1 - D.4
%----------------------------------------------------------------

% second moments, early sample

indic_1shock = 0; % use only empirical target?
indic_early  = 1; % early sample?

get_cnfctl_stats;

% second moments, 1 shock

indic_1shock = 1; % use only empirical target?
indic_early  = 0; % early sample?

get_cnfctl_stats;

% historical evolution, 1 shock

indic_1shock = 1; % use only empirical target?
indic_early  = 0; % early sample?

get_historical_evol;

% historical scenario, 1 shock

indic_1shock = 1; % use only empirical target?
indic_early  = 0; % early sample?

get_historical_scenario;

%----------------------------------------------------------------
% Table 4.1
%----------------------------------------------------------------

get_posterior_probs;

%----------------------------------------------------------------
% Table C.1
%----------------------------------------------------------------

run_var_spf_fcst_compare;

%----------------------------------------------------------------
% Table C.2
%----------------------------------------------------------------

run_var_swfactors;

%----------------------------------------------------------------
% Table C.4
%----------------------------------------------------------------

get_param_post;