function [pi_cnfctl,y_cnfctl,i_cnfctl,m_aux] = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,...
    Pi_m,Y_m,I_m,pi_x,y_x,i_x,hist_indic)

% get optimal shock sequence

m_aux = - ((A_pi * Pi_m + A_y * Y_m + A_i * I_m)' * (A_pi * Pi_m + A_y * Y_m + A_i * I_m))^(-1) ...
                * ((A_pi * Pi_m + A_y * Y_m + A_i * I_m)' * (A_pi * pi_x + A_y * y_x + A_i * i_x - wedge));

% get outcomes

pi_cnfctl = pi_x + Pi_m * m_aux;
y_cnfctl  = y_x + Y_m * m_aux;
i_cnfctl  = i_x + I_m * m_aux;

% adjust for historical evolution

if hist_indic == 2

    pi_cnfctl = pi_cnfctl(1);
    y_cnfctl  = y_cnfctl(1);
    i_cnfctl  = i_cnfctl(1);

end