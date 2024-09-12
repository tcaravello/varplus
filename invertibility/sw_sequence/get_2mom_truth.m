%% GET TRUE COUNTERFACTUAL SECOND MOMENTS
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/christianwolf/Dropbox/Research/mp_modelcnfctls/code/replication';
vintage = '/24_09_03';
task = '/invertibility';

addpath([path vintage task '/_auxiliary_functions'])
addpath([path vintage task '/sw_sequence/_subroutines'])
addpath([path vintage task '/sw_sequence/_get_mp_irfs'])

cd([path vintage task '/sw_sequence']);

%% SETTINGS

%----------------------------------------------------------------
% IS Characterization
%----------------------------------------------------------------

settings.VMA_hor          = 500; % maximal horizon for VMA representation
settings.VAR_poplaglength = 250; % lag length in population VAR
settings.IRF_hor          = 30;  % horizon for computed IRFs
settings.IS_draws         = 200; % number of draws from the identified set
settings.weight_hor       = settings.VMA_hor; % horizon for weight polynomial P(L)

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [205/255 205/255 205/255];
settings.colors.orange = [204/255 102/255 0/255];

settings.colors.list = [settings.colors.black;settings.colors.grey;settings.colors.orange];

%% VAR REPRESENTATION

%----------------------------------------------------------------
% Solve SW Model
%----------------------------------------------------------------

% model run

dynare SW_Model noclearall
clean_folder_SW

disp('I have solved and simulated the model.')

disp('Collecting model properties...')

% get law of motion for all model variables

SW_model.decision = decision(2:end,:);

% specify observables

SW_model.obs = [5 4 19]; % (r,y,pi)

% size indicators

SW_model.n_y   = size(SW_model.obs,2);
SW_model.n_eps = M_.exo_nbr;
SW_model.n_s   = M_.nspred;

% ABCD representations

SW_model.ABCD = ABCD_fun_SW(SW_model);

% delete superfluous variables

clean_workspace_SW

%----------------------------------------------------------------
% Get Population IRFs + FVDs + Shock Sequences
%----------------------------------------------------------------

[SW_model.IRF,SW_model.M,SW_model.tot_weights] = pop_analysis(SW_model,settings);

disp('...done!')

%----------------------------------------------------------------
% Get VAR Representation
%----------------------------------------------------------------

disp('Getting the VAR representation...')

VAR = popVAR(SW_model,settings);

disp('...done!')

%% PREDICT SECOND MOMENTS

%----------------------------------------------------------------
% Inputs
%----------------------------------------------------------------

% monetary policy causal effects

load sw_mp_irfs

% counterfactual rule

crpi = 1.5;
crr  = 0.8;
cry  = 0.5;
crdy = 0;

A_pi = zeros(settings.VMA_hor,settings.VMA_hor);
for t = 1:settings.VMA_hor
    A_pi(t,t) = -(1-crr) * crpi;
end

A_y = zeros(settings.VMA_hor,settings.VMA_hor);
for t = 1:settings.VMA_hor
    A_y(t,t) = -(1-crr) * (cry + crdy);
    if t > 1
        A_y(t,t-1) = (1-crr) * crdy;
    end
end

A_i = zeros(settings.VMA_hor,settings.VMA_hor);
for t = 1:settings.VMA_hor
    A_i(t,t) = 1;
    if t > 1
        A_i(t,t-1) = -crr;
    end
end

%----------------------------------------------------------------
% Counterfactual IRFs
%----------------------------------------------------------------

SW_model.IRF_cnfctl = NaN(settings.VMA_hor,3,SW_model.n_eps);

for i_eps = 1:SW_model.n_eps
    m_star = zeros(settings.VMA_hor,1);
    if i_eps == 1
        m_star(1) = -0.2290;
        for t = 2:settings.VMA_hor
            m_star(t) = 0.1999 * m_star(t-1);
        end
    end
    m_seq_aux = -(A_pi * Pi_m + A_y * Y_m + A_i * I_m)^(-1) * (A_pi * SW_model.IRF(:,3,i_eps) ...
        + A_y * SW_model.IRF(:,2,i_eps) + A_i * SW_model.IRF(:,1,i_eps) + m_star);
    SW_model.IRF_cnfctl(:,3,i_eps) = SW_model.IRF(:,3,i_eps) + Pi_m * m_seq_aux;
    SW_model.IRF_cnfctl(:,2,i_eps) = SW_model.IRF(:,2,i_eps) + Y_m * m_seq_aux;
    SW_model.IRF_cnfctl(:,1,i_eps) = SW_model.IRF(:,1,i_eps) + I_m * m_seq_aux;
end

%----------------------------------------------------------------
% Predicted Second Moments
%----------------------------------------------------------------

% variances

var_base   = zeros(SW_model.n_y,SW_model.n_y);
var_cnfctl = zeros(SW_model.n_y,SW_model.n_y);
for t = 1:settings.VMA_hor
    var_base = var_base + squeeze(SW_model.IRF(t,:,:)) * squeeze(SW_model.IRF(t,:,:))';
    var_cnfctl = var_cnfctl + squeeze(SW_model.IRF_cnfctl(t,:,:)) * squeeze(SW_model.IRF_cnfctl(t,:,:))';
end

% spectral density

omega_lb = (2*pi)/32;
omega_ub = (2*pi)/6;

n_omega = 1e2;
omega_grid = linspace(omega_lb,omega_ub,n_omega);

s_y_base_fn   = @(omega) 1/(2*pi) * sd_fun(omega,SW_model.IRF);
s_y_cnfctl_fn = @(omega) 1/(2*pi) * sd_fun(omega,SW_model.IRF_cnfctl);

s_y_base   = NaN(n_omega,3);
s_y_cnfctl = NaN(n_omega,3);
for i_omega = 1:n_omega
    omega_val = omega_grid(i_omega);
    
    s_y_base_val = s_y_base_fn(omega_val);
    s_y_base(i_omega,:) = diag(s_y_base_val);
    
    s_y_cnfctl_val = s_y_cnfctl_fn(omega_val);
    s_y_cnfctl(i_omega,:) = diag(s_y_cnfctl_val);
end

clear s_y_val

save s_y_true s_y_cnfctl s_y_base