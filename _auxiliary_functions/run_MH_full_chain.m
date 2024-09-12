function [Pi_m_collector, Y_m_collector,R_n_m_collector, m_fit_collector,...
   param_collector, log_posterior_collector,acceptance_collector, ...
   Pi_m_fit_collector, Y_m_fit_collector,R_n_m_fit_collector,chol_use] = run_MH_full_chain(settings_MH,init_val,var_scale,var_mat_draws_init,target,target_Sigma_inv,model)

% unpack settings

N_adapt      = settings_MH.N_adapt;
N_burn       = settings_MH.N_burn;
N_keep       = settings_MH.N_keep;
keep_every   = settings_MH.keep_every;
N_draws      = settings_MH.N_draws;
n_iter_print = settings_MH.n_iter_print;
T_use        = settings_MH.T_use;
T_save       = settings_MH.T_save;
start_chol   = settings_MH.start_chol; % if start_chol = 0 use diagonal prior. Otherwise, takes var_mat_draws_init as the cholesky matrix to use draws.

% placeholders

param_collector = zeros(length(init_val),N_draws);
Pi_m_collector  = zeros(T_save,T_save,N_keep);
Y_m_collector   = zeros(T_save,T_save,N_keep);
R_n_m_collector = zeros(T_save,T_save,N_keep);

m_fit_collector  = zeros(T_use,N_keep);

Pi_m_fit_collector  = zeros(T_use,N_keep);
Y_m_fit_collector   = zeros(T_use,N_keep);
R_n_m_fit_collector = zeros(T_use,N_keep);

log_posterior_collector = zeros(1,N_draws);
acceptance_collector    = zeros(N_draws,1);

param_old_orig = init_val(:);

% define a few auxiliary functions to transform parameters, in order to ensure bounds

[posterior_f_old, Pi_m_old, Y_m_old, R_n_m_old,...
    m_fit_old, Pi_m_fit_old,Y_m_fit_old, R_n_m_fit_old] = solve_model(param_old_orig,target,target_Sigma_inv,model);

draws_scale_mat = var_mat_draws_init;

var_step = var_scale * draws_scale_mat;

if length(param_old_orig) == 7
    lower_bound = [0,0,0,0,0,0,0]'; 
    upper_bound = [1,1,1,1,1,100,1]';
elseif  length(param_old_orig) == 9
    lower_bound = [0,0,0,0,0,0,0,0,0]'; 
    upper_bound = [1,1,1,1,1,100,1,1,1]';
end

% predraw normals and uniforms

normal_draws = randn(length(param_old_orig),N_draws);
uniform_draws = rand(N_draws,1);

% main loop

tic
for n = 1:N_draws
   % updating step 
    param_new_orig = param_old_orig + var_step*normal_draws(:,n);
    % reject if infeasible:
    reject = sum(param_new_orig<lower_bound)+sum(param_new_orig>upper_bound);
    if reject == 0
    % evaluate new posterior
    [posterior_f_new, Pi_m_new, Y_m_new, R_n_m_new,...
        m_fit_new, Pi_m_fit_new,Y_m_fit_new, R_n_m_fit_new] = solve_model(param_new_orig,target,target_Sigma_inv,model); 

    % acceptance step
    if uniform_draws(n)< exp(-posterior_f_new+posterior_f_old)
        % accept
        acceptance_collector(n) = 1;
        param_old_orig  = param_new_orig;
        posterior_f_old = posterior_f_new;
        Pi_m_old = Pi_m_new;
        Y_m_old = Y_m_new;
        R_n_m_old = R_n_m_new;
        
        m_fit_old = m_fit_new;
        Pi_m_fit_old = Pi_m_fit_new;
        Y_m_fit_old = Y_m_fit_new;
        R_n_m_fit_old = R_n_m_fit_new;
    else
        acceptance_collector(n) = 0;
    end
    else
        acceptance_collector(n) = 0;
    end
    % store draws
        param_collector(:,n) = param_old_orig;
        log_posterior_collector(n) = posterior_f_old;
    if n>N_burn && floor((n-N_burn)/keep_every) == (n-N_burn)/keep_every %if above the burning limit, store results
        Pi_m_collector(:,:,(n-N_burn)/keep_every) = Pi_m_old(1:T_save, 1:T_save);
        Y_m_collector(:,:,(n-N_burn)/keep_every) = Y_m_old(1:T_save, 1:T_save);
        R_n_m_collector(:,:,(n-N_burn)/keep_every) = R_n_m_old(1:T_save, 1:T_save);
        m_fit_collector(:,(n-N_burn)/keep_every) = m_fit_old;

        Pi_m_fit_collector(:,(n-N_burn)/keep_every) = Pi_m_fit_old;
        Y_m_fit_collector(:,(n-N_burn)/keep_every) = Y_m_fit_old;
        R_n_m_fit_collector(:,(n-N_burn)/keep_every) = R_n_m_fit_old;
    end
    if floor(n/n_iter_print) == n/n_iter_print
        fprintf('Time elapsed is %2.4f. Starting iter %6.0f \n',toc, n)
        fprintf('Acceptance rate for last iters is %1.4f\n',mean(acceptance_collector(n+1-n_iter_print:n)))
        tic
    end
    % adapt covariance matrix
    if n == N_adapt
        index_change = find(diag(var_mat_draws_init) ~=0);
        cov_aux = cov(param_collector(index_change,1:N_adapt)');
        chol_aux = chol(cov_aux);
        chol_use = zeros(size(var_mat_draws_init));
        for i = 1:length(index_change)
           for j = 1:length(index_change)
                chol_use(index_change(i),index_change(j)) = chol_aux(i,j);
            end
        end
        var_step = var_scale * chol_use;
    end

end

if exist('chol_use','var') == 0
    chol_use = draws_scale_mat;
end

end
