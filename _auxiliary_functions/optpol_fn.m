function [pi_z_optpol,y_z_optpol,i_z_optpol,nu_z_optpol] = optpol_fn(W_pi,W_y,W_i,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

% get optimal shock sequence

nu_z_optpol = - (Pi_m' * W_pi * Pi_m + Y_m' * W_y * Y_m + I_m' * W_i * I_m)^(-1) ...
    * (Pi_m' * W_pi * pi_z + Y_m' * W_y * y_z + I_m' * W_i * i_z);

% get outcomes

pi_z_optpol = pi_z + Pi_m * nu_z_optpol;
y_z_optpol  = y_z + Y_m * nu_z_optpol;
i_z_optpol  = i_z + I_m * nu_z_optpol;

end