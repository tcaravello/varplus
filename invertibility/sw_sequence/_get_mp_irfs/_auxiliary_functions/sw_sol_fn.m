function [Pi_m,Y_m,I_m] = sw_sol_fn(param,T);

global ctou clandaw cg curvp curvw

%----------------------------------------------------------------
% Settings
%----------------------------------------------------------------

step = 1;

%----------------------------------------------------------------
% GE Solution
%----------------------------------------------------------------

% updating matrix

A_GE = NaN(3*T,3*T);

for i_deriv = 1:3*T
    guess_seq_deriv = zeros(3*T,1);
    m_seq_deriv     = zeros(T,1);
    guess_seq_deriv(i_deriv,1) = guess_seq_deriv(i_deriv,1) + step;
    excess_demand_up = excess_demand_fn(guess_seq_deriv,m_seq_deriv,param,T);
    A_GE(:,i_deriv) = (excess_demand_up)/step;
end

% monetary shock wedge matrix

A_m = NaN(3*T,T);

for i_deriv = 1:T
    guess_seq_deriv = zeros(3*T,1);
    m_seq_deriv     = zeros(T,1);
    m_seq_deriv(i_deriv,1) = m_seq_deriv(i_deriv,1) + step;
    excess_demand_up = excess_demand_fn(guess_seq_deriv,m_seq_deriv,param,T);
    A_m(:,i_deriv) = (excess_demand_up)/step;
end

% model solution

sol_m = -A_GE^(-1) * A_m;

rk_m = sol_m(1:T,:);
kp_m = sol_m(T+1:2*T,:);
pi_m = sol_m(2*T+1:3*T,:);

%----------------------------------------------------------------
% Final IRF Matrices
%----------------------------------------------------------------

% inflation

Pi_m = pi_m;

% nominal rate

A_r = NaN(T,3*T);

for i_deriv = 1:3*T
    guess_seq_deriv = zeros(3*T,1);
    m_seq_deriv     = zeros(T,1);
    guess_seq_deriv(i_deriv,1) = guess_seq_deriv(i_deriv,1) + step;
    r_up = r_fn(guess_seq_deriv,m_seq_deriv,param,T);
    A_r(:,i_deriv) = (r_up)/step;
end

I_m = A_r(:,1:T) * rk_m + A_r(:,T+1:2*T) * kp_m + A_r(:,2*T+1:3*T) * pi_m;

% output

A_y = NaN(T,3*T);

for i_deriv = 1:3*T
    guess_seq_deriv = zeros(3*T,1);
    m_seq_deriv     = zeros(T,1);
    guess_seq_deriv(i_deriv,1) = guess_seq_deriv(i_deriv,1) + step;
    y_up = y_fn(guess_seq_deriv,m_seq_deriv,param,T);
    A_y(:,i_deriv) = (y_up)/step;
end

Y_m = A_y(:,1:T) * rk_m + A_y(:,T+1:2*T) * kp_m + A_y(:,2*T+1:3*T) * pi_m;