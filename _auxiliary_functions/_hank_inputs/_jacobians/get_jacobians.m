%% GET HANK INPUTS: JACOBIANS
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

cd([path_2 '/_steady_state'])

%% PREPARATIONS

%----------------------------------------------------------------
% Global Variables
%----------------------------------------------------------------

% aggregate parameters

global beta_hat gamma probdeath wealth_0_pos ...
     l_tax TransY_ratio BY_ratio
 
% household parameters

global a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap

% steady-state quantities

global C_SS Y_SS Trans_SS D_SS W_SS L_SS A_SS Pi_SS R_n_SS R_b_SS eis ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

global grid_a spliorder states states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

load param_agg
load param_households
load SS
load aux_use

%----------------------------------------------------------------
% Re-definitions
%----------------------------------------------------------------

r_b_SS = R_b_SS - 1;

%----------------------------------------------------------------
% Time Horizon
%----------------------------------------------------------------

global T

%% GET ALL JACOBIANS

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

global step

step = 10^(-4);

%----------------------------------------------------------------
% Prep Run
%----------------------------------------------------------------

global c_opt_vec_SS Qt_big_SS

w_seq     = zeros(T,1);
l_seq     = zeros(T,1);
tau_l_seq = zeros(T,1);
d_H_seq   = zeros(T,1);
r_n_seq   = zeros(T,1);
pi_seq    = zeros(T,1);
tau_seq   = zeros(T,1);

solve_hh_problem

c_opt_vec_SS = c_opt_t(:,:,end);

Qt_SS     = Qt;
Qt_big_SS = NaN(n_y*n_a,n_y*n_a);

for i_yT = 1:n_yT
    Qt_big_SS(:,1+(i_yT-1)*(n_a*n_yP):i_yT*(n_a*n_yP)) = repmat(Qt_SS,n_yT,1) * yT_dist(i_yT);
end
%----------------------------------------------------------------
% Labor income
%----------------------------------------------------------------

% Labor income here is in log deviations already. i.e the jacobian ouputed
% is d log(C)/d log(Y). 

vars = zeros(7,1);
vars(1) = step;

[Y_c_w,D_w] = YD_fn(vars);
[C_w,A_H_w] = jac_fn(Y_c_w,D_w);

% get various other jacobians just by rescaling

C_l = C_w;
A_H_l = A_H_w;
C_tau_l = -l_tax/(1-l_tax) * C_l;
A_H_tau_l =  -l_tax/(1-l_tax) * A_H_l;
C_d = D_SS/(W_SS*L_SS*(1-l_tax)) * C_l;
A_H_d = D_SS/(W_SS*L_SS*(1-l_tax)) * A_H_l;

%----------------------------------------------------------------
% Inflation
%----------------------------------------------------------------

vars = zeros(7,1);
vars(6) = step;

[Y_c_pi,D_pi] = YD_fn(vars);

[C_pi,A_H_pi] = jac_fn(Y_c_pi,D_pi);

%----------------------------------------------------------------
% Realized real interest rate
%----------------------------------------------------------------

C_r = - C_pi;
A_H_r = - A_H_pi;

%----------------------------------------------------------------
% Transfer
%----------------------------------------------------------------

vars = zeros(5,1);
vars(7) = step;

[Y_c_tau,D_tau] = YD_fn(vars);

[C_tau,A_H_tau] = jac_fn(Y_c_tau,D_tau);

%% Get marginal utilities (use aggregate consumption)

Lambda_w     = -(1/eis)*C_w;
Lambda_l     = -(1/eis)*C_l;
Lambda_r     = -(1/eis)*C_r;
Lambda_d     = -(1/eis)*C_d;
Lambda_tau   = -(1/eis)*C_tau;
Lambda_tau_l = -(1/eis)*C_tau_l;

beta_hank = beta;

%% SAVE RESULTS

cd(path_2)

save inputs_hank Lambda_w Lambda_l Lambda_r Lambda_d Lambda_tau Lambda_tau_l ...
    C_w C_l C_r C_d C_tau C_tau_l A_H_w A_H_l A_H_r A_H_d A_H_tau A_H_tau_l ...
    beta_hat r_b_SS l_tax C_SS W_SS L_SS D_SS A_SS Trans_SS