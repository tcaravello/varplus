% stabilize output

if cnfctl_0y == 1

A_pi = zeros(T,T);
A_y  = eye(T);
A_i  = zeros(T,T);

wedge = zeros(T,1);

end

% stabilize inflation

if cnfctl_0pi == 1

A_pi = eye(T);
A_y  = zeros(T,T);
A_i  = zeros(T,T);

wedge = zeros(T,1);

end

% stabilize interest rates

if cnfctl_0ib == 1

A_pi = zeros(T,T);
A_y  = zeros(T,T);
A_i = eye(T);

wedge = zeros(T,1);

end

% Taylor rule

if cnfctl_tylr == 1

rho_ib = 0.8;
phi_pi = 1.5;
phi_y  = 0.5;

A_pi = (1-rho_ib) * (-phi_pi * eye(T));
A_y  = (1-rho_ib) * (-phi_y * eye(T));
A_i = eye(T);
for t = 2:T
    A_i(t,t-1) = -rho_ib;
end

wedge = zeros(T,1);

end

% Nominal GDP targeting

if cnfctl_ngdp == 1

warning('Lags are probably not implemented correctly.')

A_pi = eye(T);
A_y  = eye(T);
for i = 2:T
    A_y(i,i-1) = -1;
end
A_i = zeros(T,T);

wedge = zeros(T,1);

end

% interest rate target

if cnfctl_ibtarget == 1

A_pi = zeros(T,T);
A_y  = zeros(T,T);
A_i  = eye(T);

wedge = zeros(T,1);

end

% optimal dual mandate

if cnfctl_optpol == 1

% Loss function = \sum_{t=0}^\infty \beta^t [ \lambda_\pi \pi_t^2 + \lambda_y y_t^2  ...
                        % + \lambda_i i_t^2  + \lambda_{di} (i_t - i_{t-1})^2 ]
lambda_pi   = 1;
lambda_y    = 1;
lambda_i    = 0;  
lambda_di   = 1;

beta = 1; % discount factor

W = zeros(T,T);

for t = 1:T 
    W(t,t) = beta^(t-1);
end

W_pi = lambda_pi * W;
W_y  = lambda_y * W;
W_i  = lambda_i * W;

W_i  = W_i + lambda_di * W * (1+beta);
for t = 2:T
    W_i(t-1,t) = -lambda_di * beta^(t-1);
    W_i(t,t-1) = -lambda_di * beta^(t-1);
end
   
end