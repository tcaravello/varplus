% tuning parameters if covariance matrix is non-diagonal

beta_cov_mat = 1;
theta_1_cov_mat = 1;
theta_2_cov_mat = 1;
bandwidth_cov_mat = 8;

% import targets

global Y_m_emp Y_m_var_emp Pi_m_emp Pi_m_var_emp R_n_m_emp R_n_m_var_emp Y_m_lb_emp Y_m_ub_emp Pi_m_lb_emp Pi_m_ub_emp R_n_m_lb_emp R_n_m_ub_emp

load IRFs_ad_results

target = [Pi_m_emp(1:T_use)/4; Y_m_emp(1:T_use); R_n_m_emp(1:T_use)/4];

% diagonal

Sigma_var_hat_diag = diag([Pi_m_var_emp(1:T_use)/16;Y_m_var_emp(1:T_use);R_n_m_var_emp(1:T_use)/16]);
target_Sigma_inv_diag = inv(diag([Pi_m_var_emp(1:T_use)/16;Y_m_var_emp(1:T_use);R_n_m_var_emp(1:T_use)/16]));

% non-Diagonal

Psi_i_mat = [Pi_m_draws_emp(1:T_use,:)/4;Y_m_draws_emp(1:T_use,:);R_n_m_draws_emp(1:T_use,:)/4];
Psi_i_mat_mean = mean(Psi_i_mat,2);
M_draws = size(Psi_i_mat,2);

Psi_i_var = zeros(3*T_use,3*T_use);
for n_draw = 1:M_draws
    Psi_i_var = Psi_i_var + (1/M_draws)*((Psi_i_mat(:,n_draw) -Psi_i_mat_mean)  * (Psi_i_mat(:,n_draw) - Psi_i_mat_mean)') ; 
end

% apply finite sample correction as in Christiano et al (2010)

% build triangular kernel

kernel_aux_1 = 1:-1/bandwidth_cov_mat:0;
kernel_aux_2 = [flip(kernel_aux_1(2:end)),kernel_aux_1];
mat_kernel_1 = zeros(T_use,3*T_use);
for t_mat = 1:T_use
mat_kernel_1(t_mat,t_mat:t_mat+2*bandwidth_cov_mat) = kernel_aux_2;
end
mat_kernel_same_var = mat_kernel_1(1:T_use,bandwidth_cov_mat+1:bandwidth_cov_mat+T_use).^theta_1_cov_mat;
mat_kernel_other_var = beta_cov_mat*mat_kernel_1(1:T_use,bandwidth_cov_mat+1:bandwidth_cov_mat+T_use).^theta_2_cov_mat;

adj_mat_var = [mat_kernel_same_var, mat_kernel_other_var, mat_kernel_other_var;
    mat_kernel_other_var, mat_kernel_same_var, mat_kernel_other_var;
    mat_kernel_other_var, mat_kernel_other_var, mat_kernel_same_var];

Sigma_var_hat = Psi_i_var.*adj_mat_var;
target_Sigma_inv_non_diag = inv(Sigma_var_hat);

% pick which one to use

if cov_mat == 1

    target_Sigma_inv = target_Sigma_inv_diag;

else

    target_Sigma_inv = target_Sigma_inv_non_diag;

end