function [pi_z_cnfctl,y_z_cnfctl,i_z_cnfctl,nu_z_cnfctl] = cnfctl_fn(A_pi,A_y,A_i,wedge,...
    Pi_m,Y_m,I_m,pi_z,y_z,i_z);

% get optimal shock sequence

nu_z_cnfctl = - ((A_pi * Pi_m + A_y * Y_m + A_i * I_m)' * (A_pi * Pi_m + A_y * Y_m + A_i * I_m))^(-1) ...
                * ((A_pi * Pi_m + A_y * Y_m + A_i * I_m)' * (A_pi * pi_z + A_y * y_z + A_i * i_z - wedge));

% get outcomes

pi_z_cnfctl = pi_z + Pi_m * nu_z_cnfctl;
y_z_cnfctl  = y_z + Y_m * nu_z_cnfctl;
i_z_cnfctl  = i_z + I_m * nu_z_cnfctl;