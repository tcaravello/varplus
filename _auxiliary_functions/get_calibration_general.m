global beta eis  varphi  alpha ...
   delta  delta_1  delta_2  eta  rho_tr  phi_pi  phi_y  phi_dy  tau_b_fix  tau_b_dist ...
tau_l_rate  r_SS  L_SS  C_SS  I_SS  Y_SS  W_SS  D_SS  B_SS  A_SS  Tau_SS  G_SS T F_mat L_mat D_mat Id_mat;

beta    = 0.99;
eis     = 0.5;
varphi  = 0.5;
alpha	= 0.36;
delta	= 0.025;
delta_1	= 0.2;
delta_2	= 0.05;
eta	    = 0.2;
rho_tr	= 0.85;
phi_pi	= 2;
phi_y	= 0.25;
phi_dy	= 0.3;

tau_b_fix = 0.0; % adjustment via lump sum taxes
tau_b_dist = 0.15; % adjustment via distortionary taxes
tau_l_rate = 0.3;
Tau_Y = -0.05; % with a minus, these are transfers
B_Y_SS = 1.04;

r_SS	= 1/beta - 1;
L_SS	= 1;
K_SS	= (alpha/(r_SS + delta))^(1/(1-alpha)) * L_SS;
Y_SS	= K_SS^alpha * L_SS^(1-alpha);
I_SS	= delta * K_SS;
W_SS	= (1-alpha) * K_SS^alpha * L_SS^(-alpha);
D_SS	= Y_SS - W_SS * L_SS - I_SS;
B_SS    = B_Y_SS * Y_SS;
A_SS    = B_SS; 
Tau_SS  = Tau_Y * Y_SS;
G_SS    = tau_l_rate*(1-alpha)*Y_SS + Tau_SS - r_SS * A_SS;
C_SS	= Y_SS - I_SS - G_SS;

% create operator matrices

% forward
F_mat = zeros(T,T);
F_mat(1:end-1,2:end) = eye(T-1);

% lag
L_mat = zeros(T,T);
L_mat(2:end,1:end-1) = eye(T-1);

% difference operator
D_mat = eye(T) - L_mat;

Id_mat = eye(T);
