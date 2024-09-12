function [M_c,M_b] = jac_fn(Y_c,D);

%% PREPARATIONS: GLOBAL VARIABLES

% aggregate parameters

global beta beta_hat gamma probdeath wealth_0_pos ...
     kappa ...
     y_tax TransY_ratio BY_ratio rho_b lambda_pi lambda_y ...
     rho_s rho_d sigma_s sigma_d
 
% household parameters

global a_lb n_y n_yP n_yT grid_y grid_yP grid_yT y_dist yP_dist yT_dist Pi_y Pi_yP Pi_yT annuity_gap

% steady-state quantities

global C_SS Y_SS Trans_SS Pi_SS R_n_SS R_b_SS A_SS ...
    r_b_grid r_b_SS mutilde_SS c_opt_SS ap_opt_SS lambda_SS lambda_vec_SS lambdafull_SS

% other quantities

global grid_a spliorder states states_yP states_a Phi_yP Emat_yP fspaceerga fspace ...
    n_a n_s a_min a_max

% settings

global T step

% auxiliary variables

global c_opt_vec_SS Qt_big_SS

%% MAIN COMPUTATIONS

%----------------------------------------------------------------
% Construct P
%----------------------------------------------------------------

P_c = NaN(n_a*n_y,T-1);

for i_yT = 1:n_yT
    P_c(1+(i_yT-1)*(n_a*n_yP):i_yT*(n_a*n_yP),1) = c_opt_vec_SS(:,i_yT);
end

for t = 2:T-1
    P_c(:,t) = Qt_big_SS * P_c(:,t-1);
end

P_b = NaN(n_a*n_y,T-1);
P_b(:,1) = repmat(grid_a',n_y,1);
for t = 2:T
    P_b(:,t) = Qt_big_SS * P_b(:,t-1);
end

%----------------------------------------------------------------
% Construct F
%----------------------------------------------------------------

F_c   = [Y_c;P_c'*D];
F_b   = P_b'*D;

%----------------------------------------------------------------
% Construct M
%----------------------------------------------------------------

M_c = F_c;
for t = 2:T
    M_c(2:end,t) = M_c(2:end,t) + M_c(1:end-1,t-1);
end
M_c = M_c/C_SS;

M_b = F_b;
for t = 2:T
    M_b(2:end,t) = M_b(2:end,t) + M_b(1:end-1,t-1);
end
M_b = M_b/A_SS;