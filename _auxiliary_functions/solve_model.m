function [error_fit, varargout] = solve_model(param,target,target_Sigma_inv,model)

% given parameter values, this function solves the model and computes the
% quasi-likelihood.

if strcmp(model, '/rank')
    global prior_h
else
    global prior_theta
end

global prior_kappa prior_calvo_p prior_zeta_p prior_calvo_w prior_zeta_w prior_psi_uti prior_m_d prior_m_f

global T_use n_shock

%----------------------------------------------------------------
% Get Model-implied IRFs
%----------------------------------------------------------------

if nargout > 8
    [Pi_m_model, Y_m_model, R_n_m_model, C_m_model, Inv_m_model] = compute_irfs_model(param,model);
else
    [Pi_m_model, Y_m_model, R_n_m_model] = compute_irfs_model(param,model);
end

%----------------------------------------------------------------
% Get Loss
%----------------------------------------------------------------

% this part uses the news shocks to get the Fitted IRF.

% first, get the 1:H columns of the jacobians, and stack them.
A_model = [Pi_m_model(1:T_use,1:n_shock); Y_m_model(1:T_use,1:n_shock); R_n_m_model(1:T_use,1:n_shock)];

% initialize a vector of zeros for the full tilde{v}^{*}
m_fit = zeros(T_use,1);

% fill in the first H components with the "best fit" shock.
% Notice that this is essentially a GLS formula:
% where A_model = X,  target_Sigma_inv = Sigma^{-1} and target = Y.
m_fit(1:n_shock) = ( A_model' * target_Sigma_inv * A_model)\ ... 
     (A_model' * target_Sigma_inv * target);

% get the best fit values as Jacobian(1:H,1:H) * tilde{v}^{*}
Pi_m_fit = Pi_m_model(1:T_use,1:T_use) * m_fit;
Y_m_fit  = Y_m_model(1:T_use,1:T_use) * m_fit;
R_n_m_fit  = R_n_m_model(1:T_use,1:T_use) * m_fit;

if length(param)>7 % if param vector is long, then we are using behavioral models.
    behav_prior =  log(prior_m_d(param(8)))+log(prior_m_f(param(9)));
else % otherwise are not so ignore this term.
    behav_prior = 0;
end

% use loop because RANK has habit parameter "h", whereas HANK has info
% stickiness parameter "theta"
if strcmp(model, '/rank')
    error_fit = 0.5*(([Pi_m_fit ; Y_m_fit; R_n_m_fit] - target)' *  target_Sigma_inv * ([Pi_m_fit ; Y_m_fit; R_n_m_fit] - target)) ...
   - log(prior_h(param(1))) - log(prior_calvo_p(param(2))) - log(prior_zeta_p(param(3))) ...
   - log(prior_calvo_w(param(4))) - log(prior_zeta_w(param(5))) - log(prior_kappa(param(6))) ...
   - log(prior_psi_uti(param(7)))- behav_prior;
else
   error_fit = 0.5*(([Pi_m_fit ; Y_m_fit; R_n_m_fit] - target)' *  target_Sigma_inv * ([Pi_m_fit ; Y_m_fit; R_n_m_fit] - target)) ...
   - log(prior_theta(param(1))) - log(prior_calvo_p(param(2))) - log(prior_zeta_p(param(3))) ...
   - log(prior_calvo_w(param(4))) - log(prior_zeta_w(param(5))) - log(prior_kappa(param(6))) ...
   - log(prior_psi_uti(param(7))) - behav_prior;
end

% depending on how many variables you put when calling the function, you
% can obtain more outputs.

if nargout==8
    varargout{1} = Pi_m_model;
    varargout{2} = Y_m_model;
    varargout{3} = R_n_m_model;
    varargout{4} = m_fit;
    varargout{5} = Pi_m_fit;
    varargout{6} = Y_m_fit;
    varargout{7} = R_n_m_fit;
elseif nargout== 10
    varargout{1} = Pi_m_model;
    varargout{2} = Y_m_model;
    varargout{3} = R_n_m_model;
    varargout{4} = m_fit;
    varargout{5} = Pi_m_fit;
    varargout{6} = Y_m_fit;
    varargout{7} = R_n_m_fit;
    varargout{8} = C_m_model;
    varargout{9} = Inv_m_model;
elseif nargout == 12
    varargout{1} = Pi_m_model;
    varargout{2} = Y_m_model;
    varargout{3} = R_n_m_model;
    varargout{4} = m_fit;
    varargout{5} = Pi_m_fit;
    varargout{6} = Y_m_fit;
    varargout{7} = R_n_m_fit;
    varargout{8} = C_m_model;
    varargout{9} = Inv_m_model;
    C_m_fit = C_m_model(1:T_use,1:T_use) * m_fit;
    Inv_m_fit =  Inv_m_model(1:T_use,1:T_use) * m_fit;
    varargout{10} = C_m_fit;
    varargout{11} = Inv_m_fit;
end    

end

%----------------------------------------------------------------
% Auxiliary Functions
%----------------------------------------------------------------

function  varargout = compute_irfs_model(param,model) 

global rho_tr  phi_pi  phi_y  phi_dy

global T L_mat D_mat Id_mat

J_s = get_transformed_J_s(param,model); % get transformed params.

% model solution

step = 1;
% each column of this matrix is a different shock (i.e pi_shocks 1 to T, w
% shocks 1 to T and Y shocks 1 to t). Each row is the response of a target
% equation.
A_GE = NaN(3*T,3*T);

for i_deriv = 1:3
    guess_seq_deriv = zeros(3*T,T);
    m_seq_deriv     = zeros(T,T);
    guess_seq_deriv(1+(i_deriv-1)*T:i_deriv*T,:) = eye(T)*step;
    excess_demand_up = excess_demand_fn_given_Js(guess_seq_deriv,m_seq_deriv,J_s,param);
    A_GE(:,1+(i_deriv-1)*T:i_deriv*T) = (excess_demand_up)/step;
end

% monetary shock wedge matrix

guess_seq_deriv = zeros(3*T,T);
m_seq_deriv     = eye(T)*step;
excess_demand_up = excess_demand_fn_given_Js(guess_seq_deriv,m_seq_deriv,J_s,param);
A_m = (excess_demand_up)/step;

% model solution

sol_m = -A_GE\A_m;

Pi_m_model = sol_m(1:T,:);
W_m_model  = sol_m(T+1:2*T,:);
Y_m_model  = sol_m(2*T+1:3*T,:);

%----------------------------------------------------------------
% Final IRF Matrices
%----------------------------------------------------------------

% nominal rates
R_n_m_model = (Id_mat - rho_tr * L_mat)\(phi_pi * Pi_m_model + (phi_y * Id_mat  +  phi_dy * D_mat) * Y_m_model + Id_mat)*(1-rho_tr);

if nargout == 3
    varargout{1} = Pi_m_model;
    varargout{2} = Y_m_model;
    varargout{3} = R_n_m_model;
elseif nargout == 5 % also compute consumption and investment

    global Y_SS C_SS I_SS beta F_mat alpha delta

    calvo_p = param(2); % PC slope
    zeta_p = param(3); % PC intertia
    calvo_w = param(4); % Wage PC slope
    zeta_w = param(5); % Wage PC intertia
    kappa = param(6); %investment adjustment cost.
    psi_uti = param(7); % capacity utilization
    zeta = psi_uti/(1-psi_uti); %capacity utilization parameter.

    if length(param)>7
        m_d = param(8);
        m_f = param(9);
    else
        m_d = 1;
        m_f = 1;
    end

    kappa_p = (1-beta*calvo_p)*(1-calvo_p)/calvo_p;
    kappa_w = (1-beta*calvo_w)*(1-calvo_w)/calvo_w;

    beta_p = beta*calvo_p*m_f*(1+kappa_p/(1-beta*calvo_p*m_f));
    beta_w = beta*calvo_w*m_f*(1+kappa_w/(1-beta*calvo_w*m_f));

    % marginal costs

    p_I_seq = 1/kappa_p * (Id_mat - zeta_p * L_mat - beta_p * (F_mat - zeta_p * Id_mat)) * Pi_m_model;

    % hours worked

    l_seq = p_I_seq + Y_m_model - W_m_model;
    
    % capacity utilization
    
    u_seq = 1/zeta * (p_I_seq + (1-alpha)/alpha * (l_seq -  Y_m_model));
    
    % capital and investment
    
    k_seq_lag = 1/alpha * (W_m_model- p_I_seq - alpha * u_seq + alpha * l_seq);
    k_seq = [k_seq_lag(2:T,:); zeros(1,T)] ;
    Inv_m_model = 1/delta * (k_seq - (1-delta) * k_seq_lag);
    C_m_model   = (Y_SS*Y_m_model - I_SS * Inv_m_model)/C_SS;

    varargout{1} = Pi_m_model;
    varargout{2} = Y_m_model;
    varargout{3} = R_n_m_model;
    varargout{4} = C_m_model;
    varargout{5} = Inv_m_model;
end

end

function excess_demand = excess_demand_fn_given_Js(guess_seq,shock_seq,J_s,param)

global beta varphi  alpha ...
        delta delta_1 delta_2 eta rho_tr phi_pi phi_y phi_dy tau_b_fix tau_b_dist ...
        tau_l_rate r_SS L_SS I_SS Y_SS W_SS D_SS B_SS A_SS Tau_SS;

global T L_mat F_mat D_mat Id_mat

%----------------------------------------------------------------
% Collect Inputs
%----------------------------------------------------------------

calvo_p = param(2); % PC slope
zeta_p = param(3); % PC intertia
calvo_w = param(4); % Wage PC slope
zeta_w = param(5); % Wage PC intertia
kappa = param(6); %investment adjustment cost.
psi_uti = param(7); % capacity utilization
zeta = psi_uti/(1-psi_uti); %capacity utilization parameter.

if length(param)>7
    m_d = param(8);
    m_f = param(9);
else
    m_d = 1;
    m_f = 1;
end

kappa_p = (1-beta*calvo_p)*(1-calvo_p)/calvo_p;
kappa_w = (1-beta*calvo_w)*(1-calvo_w)/calvo_w;

beta_p = beta*calvo_p*m_f*(1+kappa_p/(1-beta*calvo_p*m_f));
beta_w = beta*calvo_w*m_f*(1+kappa_w/(1-beta*calvo_w*m_f)); %we use here the same behavioral parameter as the households.

Lambda_w     = J_s.Lambda_w;
Lambda_l     = J_s.Lambda_l;
Lambda_r     = J_s.Lambda_r ;
Lambda_d     = J_s.Lambda_d ;
Lambda_tau   = J_s.Lambda_tau ;
Lambda_tau_l = J_s.Lambda_tau_l ;

C_w     = J_s.C_w;
C_l     = J_s.C_l;
C_r     = J_s.C_r ;
C_d     = J_s.C_d ;
C_tau   = J_s.C_tau ;
C_tau_l = J_s.C_tau_l;

A_H_w     = J_s.A_H_w;
A_H_l     = J_s.A_H_l;
A_H_r     = J_s.A_H_r ;
A_H_d     = J_s.A_H_d ;
A_H_tau   = J_s.A_H_tau ;
A_H_tau_l = J_s.A_H_tau_l ;

n_cols  = size(guess_seq,2);
pi_seq  = guess_seq(1:T,:);
w_seq   = guess_seq(T+1:2*T,:);
y_seq   = guess_seq(2*T+1:3*T,:);

m_seq   = shock_seq;

%----------------------------------------------------------------
% Get Outcomes
%----------------------------------------------------------------

% marginal costs

p_I_seq = 1/kappa_p * (Id_mat - zeta_p * L_mat - beta_p * (F_mat - zeta_p * Id_mat)) * pi_seq;

% nominal rates

r_n_seq = (Id_mat - rho_tr * L_mat)\(phi_pi * pi_seq + (phi_y * Id_mat  +  phi_dy * D_mat) * y_seq + m_seq)*(1-rho_tr);

% hours worked

l_seq = p_I_seq + y_seq - w_seq;

% capacity utilization

u_seq = 1/zeta * (p_I_seq + (1-alpha)/alpha * (l_seq - y_seq));

% capital and investment

k_seq_lag = 1/alpha * (w_seq - p_I_seq - alpha * u_seq + alpha * l_seq);
k_seq = [k_seq_lag(2:T,:); zeros(1,n_cols)] ;
i_seq = 1/delta * (k_seq - (1-delta) * k_seq_lag);

% Behavioral frictions: just added m_d. To turn off set m_d = 1

% capital price 
q_seq = kappa * (1+1/(1+r_SS)) * i_seq - kappa * 1/(1+r_SS) *m_d*[i_seq(2:end,:);zeros(1,n_cols)] - kappa * [zeros(1,n_cols);i_seq(1:end-1,:)];

% added the 1/m_d to account for the behavioral friction in investment..
k_seq_supply = -1/((1-alpha) * (1 - 1/(1+r_SS) * (1-delta))) * (q_seq/m_d + r_n_seq - [pi_seq(2:end,:);zeros(1,n_cols)] - 1/(1+r_SS) * (1-delta) * [q_seq(2:end,:);zeros(1,n_cols)] ...
   - (1-1/(1+r_SS) * (1-delta)) * ([p_I_seq(2:end,:);zeros(1,n_cols)] - (1-alpha) * [u_seq(2:end,:);zeros(1,n_cols)] + (1-alpha) * [l_seq(2:end,:);zeros(1,n_cols)]));

i_seq_supply = 1/delta * (k_seq_supply - (1-delta) * [zeros(1,n_cols);k_seq_supply(1:end-1,:)]);

% dividends
d_seq = (Y_SS * y_seq - W_SS * L_SS * (w_seq + l_seq) - I_SS * i_seq)/ D_SS;

% financial intermediary

a_I_seq = (Id_mat - (1-delta_2)*(1+r_SS) * L_mat)\((1-delta_1)*D_SS * d_seq);

d_H_seq = delta_1 * d_seq + (delta_2/D_SS) * (1 + r_SS) * [zeros(1,n_cols);a_I_seq(1:end-1,:)];

% asset market

p_bond_seq = - (Id_mat-(1-eta)/(1+r_SS) * F_mat)\r_n_seq;

r_seq      = - pi_seq - [zeros(1,n_cols);p_bond_seq(1:end-1,:)] + (1-eta)/(1+r_SS) * p_bond_seq;

% fiscal policy

b_seq = (Id_mat - (1+r_SS- tau_l_rate * tau_b_dist - tau_b_fix) * L_mat)\((1 + r_SS) * r_seq - (tau_l_rate * W_SS * L_SS / B_SS) * (w_seq + l_seq));

tau_l_seq = (B_SS)/(W_SS * L_SS) * tau_b_dist * [zeros(1,n_cols); b_seq(1:T-1,:)];

tau_seq = (B_SS /Tau_SS) * tau_b_fix * [zeros(1,n_cols); b_seq(1:T-1,:)];

% consumption & savings

lambda_seq = Lambda_w * w_seq + Lambda_l * l_seq + Lambda_r * r_seq + Lambda_d * d_H_seq + Lambda_tau * tau_seq + Lambda_tau_l * tau_l_seq;
c_seq      = C_w * w_seq + C_l * l_seq + C_r * r_seq + C_d * d_H_seq + C_tau * tau_seq + C_tau_l * tau_l_seq;
a_H_seq    = A_H_w * w_seq + A_H_l * l_seq + A_H_r * r_seq + A_H_d * d_H_seq + A_H_tau * tau_seq + A_H_tau_l * tau_l_seq;

% labor supply

pi_w_seq = w_seq - [zeros(1,n_cols);w_seq(1:T-1,:)] + pi_seq;

chi_seq = 1/kappa_w * ((pi_w_seq - zeta_w * [zeros(1,n_cols);pi_seq(1:T-1,:)]) - beta_w * ([pi_w_seq(2:T,:);zeros(1,n_cols)] - zeta_w * pi_seq));

l_seq_supply = varphi * (chi_seq + (w_seq + lambda_seq) - tau_l_rate/(1-tau_l_rate) * tau_l_seq);

%----------------------------------------------------------------
% Check Accuracy
%----------------------------------------------------------------

%asset market clearing

excess_demand_1 = A_SS * a_H_seq + a_I_seq - B_SS * b_seq;

% labor market-clearing

excess_demand_2 = l_seq - l_seq_supply;

% capital market-clearing

excess_demand_3 = [i_seq(1:T-1,:) - i_seq_supply(1:T-1,:);k_seq_lag(1,:)];

% collect all wedges

excess_demand = [excess_demand_1;excess_demand_2;excess_demand_3];

end

function J_s = get_transformed_J_s(param,model)

J_s = struct();

if strcmp(model, '/rank')

global beta eis Id_mat F_mat L_mat T C_SS A_SS W_SS L_SS tau_l_rate Tau_SS D_SS r_SS

h = param(1);
if length(param)>7
    m_d = param(8);
else
    m_d = 1;
end

gamma = 1/eis;

% first, build FI J_Zs:

h_vec = (1-h.^((1:T)'))/(1-h);
h_vec_a = 1-beta*h - (1-beta)*h.^((2:T+1)');
beta_vec = beta.^((0:T-1));

Forward_sum_mat = ((Id_mat-beta * F_mat)\Id_mat - Id_mat);

%----------------------------------------------------------------
% FI After tax income
%----------------------------------------------------------------

C_z_FI = ((1-beta)*(1-beta*h)/C_SS)*h_vec*beta_vec;
A_H_z_FI = (beta/(A_SS*(1-h)))*h_vec_a*beta_vec -Forward_sum_mat/A_SS;
Lambda_z_FI = -(gamma/((1-beta*h)*(1-h)))*(Id_mat - h * L_mat ) * C_z_FI + (gamma/((1-beta*h)*(1-h)))*beta*h*(F_mat *C_z_FI - h*C_z_FI);

% now build H_rs

FR_mat = ((1-beta*h)*(1-h)/(gamma*beta*h))*((Id_mat- (beta*h*F_mat))\Id_mat - Id_mat);
FS_mat_aux = zeros(T,T);

for j=1:T
    FS_mat_aux(j+1:end,j) = h_vec(1:end-j);
end

FS_mat =  FS_mat_aux *FR_mat;
psi_T = -((1-beta)*(1-beta*h)*(1-h))/((1-beta^(T))*(1-beta*h)-h*(1-beta)*(1-(beta*h)^(T)));
c_end = (((1-psi_T*((1-h^(T))*beta^(T))/((1-h)*(1-beta))))^(-1))*(psi_T * ((1-h^(T))/(1-h)) * beta_vec * FS_mat + FS_mat(T,:));

H_c_r = psi_T * h_vec * (beta_vec*FS_mat + c_end *(beta^(T))/(1-beta)) + FS_mat;
beta_vec_2 =  beta.^(((T-1):-1:0)');
H_a_r = (C_SS/A_SS)*(Forward_sum_mat*H_c_r + (beta/(1-beta))*beta_vec_2*c_end); 

%----------------------------------------------------------------
% FI Interest rates
%----------------------------------------------------------------

C_r_FI = H_c_r + C_z_FI * (1+r_SS)*A_SS;
A_H_r_FI = H_a_r + A_H_z_FI * (1+r_SS)*A_SS;
Lambda_r_FI = -(gamma/((1-beta*h)*(1-h)))*(Id_mat - h * L_mat ) * C_r_FI  + (gamma/((1-beta*h)*(1-h)))*beta*h*(F_mat *C_r_FI - h*C_r_FI);

if m_d > 0.99999 % for numerical stability
% use full information
    C_z_use = C_z_FI; 
    A_H_z_use = A_H_z_FI;
    Lambda_z_use = Lambda_z_FI; 
    C_r_use =  C_r_FI;
    A_H_r_use = A_H_r_FI ;
    Lambda_r_use = Lambda_r_FI;
else
    C_z_use = add_cognitive_disc(C_z_FI,m_d); 
    A_H_z_use = add_cognitive_disc(A_H_z_FI,m_d);
    Lambda_z_use = add_cognitive_disc(Lambda_z_FI,m_d); 
    C_r_use =  add_cognitive_disc(C_r_FI,m_d);
    A_H_r_use = add_cognitive_disc(A_H_r_FI,m_d) ;
    Lambda_r_use = add_cognitive_disc(Lambda_r_FI,m_d);
end

%----------------------------------------------------------------
% Interest rates
%----------------------------------------------------------------

J_s.C_r      = C_r_use;
J_s.A_H_r    = A_H_r_use;
J_s.Lambda_r = Lambda_r_use;

%----------------------------------------------------------------
% Wages
%----------------------------------------------------------------

J_s.Lambda_w = W_SS * L_SS * (1-tau_l_rate) * Lambda_z_use; 
J_s.C_w = W_SS * L_SS * (1-tau_l_rate) * C_z_use ;
J_s.A_H_w = W_SS * L_SS * (1-tau_l_rate)* A_H_z_use ;

%----------------------------------------------------------------
% Labor income
%----------------------------------------------------------------

J_s.Lambda_l = J_s.Lambda_w;
J_s.C_l      = J_s.C_w;
J_s.A_H_l    = J_s.A_H_w;

%----------------------------------------------------------------
% Labor taxes
%----------------------------------------------------------------

J_s.Lambda_tau_l =  -W_SS * L_SS * tau_l_rate * Lambda_z_use; 
J_s.C_tau_l = -W_SS * L_SS * tau_l_rate * C_z_use; 
J_s.A_H_tau_l = -W_SS * L_SS * tau_l_rate * A_H_z_use; 

%----------------------------------------------------------------
% Dividends
%----------------------------------------------------------------

J_s.Lambda_d =  D_SS * Lambda_z_use; 
J_s.C_d = D_SS * C_z_use; 
J_s.A_H_d = D_SS*  A_H_z_use; 

%----------------------------------------------------------------
% Lump sum taxes
%----------------------------------------------------------------

J_s.Lambda_tau =  -Tau_SS * Lambda_z_use; 
J_s.C_tau = -Tau_SS * C_z_use; 
J_s.A_H_tau = -Tau_SS * A_H_z_use; 

else

global Lambda_w  Lambda_l  Lambda_r  Lambda_d  Lambda_tau  Lambda_tau_l ...
    C_w  C_l  C_r  C_d  C_tau  C_tau_l  A_H_w  A_H_l  A_H_r  A_H_d  A_H_tau  A_H_tau_l tau_l_rate eis

theta = param(1);
if length(param)>7
    m_d = param(8);
else
    m_d = 1;
end

if m_d > 0.99999 % use this for numerical stability.

    theta = param(1);
    
    J_s.C_w = add_sticky(C_w, theta);
    J_s.C_l = J_s.C_w;
    J_s.C_r = add_sticky(C_r, theta);
    J_s.C_d = add_sticky(C_d, theta);
    J_s.C_tau = add_sticky(C_tau, theta);
    J_s.C_tau_l = -(tau_l_rate/(1-tau_l_rate))*J_s.C_w;
    
    J_s.A_H_w = add_sticky(A_H_w, theta);
    J_s.A_H_l = J_s.A_H_w;
    J_s.A_H_r = add_sticky(A_H_r, theta);
    J_s.A_H_d = add_sticky(A_H_d, theta);
    J_s.A_H_tau = add_sticky(A_H_tau, theta);
    J_s.A_H_tau_l = -(tau_l_rate/(1-tau_l_rate))*J_s.A_H_w;
    
    J_s.Lambda_w = -(1/eis)*J_s.C_w;
    J_s.Lambda_l = -(1/eis)*J_s.C_l;
    J_s.Lambda_r = -(1/eis)*J_s.C_r;
    J_s.Lambda_d = -(1/eis)*J_s.C_d;
    J_s.Lambda_tau = -(1/eis)*J_s.C_tau;
    J_s.Lambda_tau_l = -(1/eis)*J_s.C_tau_l;

else % if we want to add cognitive discounting

    J_s.C_w = add_sticky(add_cognitive_disc(C_w,m_d), theta);
    J_s.C_l = J_s.C_w;
    J_s.C_r = add_sticky(add_cognitive_disc(C_r,m_d), theta);
    J_s.C_d = add_sticky(add_cognitive_disc(C_d,m_d), theta);
    J_s.C_tau = add_sticky(add_cognitive_disc(C_tau,m_d), theta);
    J_s.C_tau_l = -(tau_l_rate/(1-tau_l_rate))*J_s.C_w;
    
    J_s.A_H_w = add_sticky(add_cognitive_disc(A_H_w,m_d), theta);
    J_s.A_H_l = J_s.A_H_w;
    J_s.A_H_r = add_sticky(add_cognitive_disc(A_H_r,m_d), theta);
    J_s.A_H_d = add_sticky(add_cognitive_disc(A_H_d, m_d) , theta);
    J_s.A_H_tau = add_sticky(add_cognitive_disc(A_H_tau, m_d) , theta);
    J_s.A_H_tau_l = -(tau_l_rate/(1-tau_l_rate))*J_s.A_H_w;
    
    J_s.Lambda_w = -(1/eis)*J_s.C_w;
    J_s.Lambda_l = -(1/eis)*J_s.C_l;
    J_s.Lambda_r = -(1/eis)*J_s.C_r;
    J_s.Lambda_d = -(1/eis)*J_s.C_d;
    J_s.Lambda_tau = -(1/eis)*J_s.C_tau;
    J_s.Lambda_tau_l = -(1/eis)*J_s.C_tau_l;

end

end

end

function J_sticky = add_sticky(J,theta)

global T

J_sticky = zeros(T, T);
J_sticky(:,1) = J (:,1);
J_sticky(1,2:end) = (1-theta) * J (1,2:end);

% build the rest of J_sticky recurisevely

for t = 2:T
    for s = 2:T
        J_sticky(t,s) = theta * J_sticky(t-1,s-1) + (1-theta) * J (t,s); 
    end
end

end

function J_cd = add_cognitive_disc(J,m_d)

global T

J_cd = zeros(T, T);
J_cd(:,1) = J(:,1);
for j=2:T
    J_cd (:,j) = m_d^(j-1) * (J(:,j)- [0;J(1:end-1,j-1)]) + [0; J_cd(1:end-1,j-1)];
end

end