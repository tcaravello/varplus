%% GET SW MODEL MONETARY POLICY CAUSAL EFFECTS
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

addpath([path vintage task '/sw_sequence/_get_mp_irfs/_auxiliary_functions'])

cd([path vintage task '/sw_sequence/_get_mp_irfs']);

%% MODEL PARAMETERS & SETTINGS

% fixed parameters

global ctou clandaw cg curvp curvw

ctou    = 0.025;
clandaw	= 1.5;
cg		= 0.18;
curvp	= 10;
curvw	= 10;

% estimated parameters

csadjcost  = 6.3144;
csigma	   = 1.2679;
chabb      = 0.8056;
cprobw     = 0.7668;
csigl      = 2.5201;
cprobp     = 0.5304;
cindw      = 0.5345;
cindp      = 0.1779;
czcap      = 0.3597;
cfc	       = 1.6670;
crpi       = 3;
crr        = 0;
cry        = 0;
crdy       = 0;
constepinf = 0.6365;
constebeta = 0.1126;
constelab  = 1.3263;
ctrend     = 0.5113;
cgy	       = 0.5881;
calfa	   = 0.2024;

param = [csadjcost,csigma,chabb,cprobw,csigl,cprobp,cindw,...
    cindp,czcap,cfc,crpi,crr,cry,crdy,...
    constepinf,constebeta,constelab,ctrend,cgy,calfa];

n_param = length(param);

% solution settings

T = 500;

%% BASE MODEL SOLUTION

[Pi_m,Y_m,I_m] = sw_sol_fn(param,T);

%% SAVE RESULTS

save sw_mp_irfs Pi_m Y_m I_m