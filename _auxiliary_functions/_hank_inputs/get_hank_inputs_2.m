%% GET HANK INPUTS
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

warning('off','MATLAB:dispatcher:nameConflict')
path_2 = [path vintage '/_auxiliary_functions/_hank_inputs']; 
addpath(genpath([path_2 '/_aux']))
addpath([path_2 '/_income_process'])
addpath([path_2 '/_steady_state'])
addpath([path_2 '/_jacobians'])

cd(path_2);

%% COMPUTE STEADY STATE

disp('I am solving for the steady state of the model.')

get_ss

disp('Done!')

%% COMPUTE JACOBIANS

disp('I am computing the PE Jacobians.')

get_jacobians

disp('Done!')

%% LOAD INPUTS

load inputs_hank.mat
