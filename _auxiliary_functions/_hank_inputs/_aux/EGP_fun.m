function [c_opt,ap_opt,mutilde_upd] = EGP_fun(mutilde_opt,r_b_grid,W,L,D,Trans,...
    beta,gamma,l_tax,states,grid_a,Emat_yP,grid_yT,yT_dist);

% preparations

n_a  = size(grid_a,2);
n_yP = size(states,1)/n_a;
n_yT = length(grid_yT);

a_lb = min(grid_a);

% consumption and savings on endogenous grid
    
c_endo = repmat((beta .* mutilde_opt).^(-1/gamma),1,n_yT);
a_endo = 1./(1+r_b_grid) .* (states(:,1) + c_endo ...
    - (states(:,2) * grid_yT') .* ((1-l_tax).*W.*L + D) - Trans);

% interpolate

for i_yT = 1:n_yT
    a_endo_temp = reshape(a_endo(:,i_yT),n_a,n_yP);
    V_temp = interpOneD_vec(repmat(grid_a',1,n_yP),a_endo_temp);
    ap_temp(:,i_yT) = V_temp * repmat(grid_a',n_yP,1);
    
end

% deal with constraint

constr = (ap_temp < a_lb);
ap_opt = constr .* a_lb + (1-constr) .* ap_temp;

% get final consumption

c_opt = (1 + r_b_grid) .* states(:,1) + (states(:,2) * grid_yT') .* ((1-l_tax).*W.*L + D) ...
    + Trans - ap_opt;

% get updated mu_tilde

mutilde_upd = Emat_yP * (((1 + r_b_grid) .* c_opt.^(-gamma)) * yT_dist');

end