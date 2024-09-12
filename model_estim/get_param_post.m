%% ANALYZE SAMPLING RESULTS
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

task = '/model_estim';

addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage '/var_inputs/_results'])
addpath([path vintage '/suff_stats/ratex']);
addpath([path vintage '/suff_stats/behavioral']);
addpath([path vintage '/suff_stats/mix']);
addpath([path vintage task '/_results'])

cd([path vintage task]);

%% PRODUCE ESTIMATION TABLE

get_priors_general;

n_models = 4;

f = fopen(fullfile('_results', strcat('param_prior_post', '.tex')), 'w'); % open file for writing

fprintf(f,'\\begin{table} \n');
fprintf(f,'\\centering \n');
fprintf(f,'\\begin{tabular}{|cc|ccc|c|ccccc|} \\hline \n');
fprintf(f,'& & \\multicolumn{3}{|c|}{Prior} & & \\multicolumn{5}{|c|}{Posterior} \\\\ \n');
fprintf(f,'Model & Parameter & Dist. & Mean & St. Dev & & Mode & Mean & Median &  5 percent & 95 percent \\\\ \\hline \n');

% load each model, compute the statistics, in order to get the full table

for i_model = 1:n_models

    % RANK

    if i_model == 1
    load rank_draws_main
    load rank_draws_other
    load params_rank_mode
    model = '\rank';
    model_name = 'RANK - RE';

    % HANK

    elseif i_model == 2
    load hank_draws_main
    load hank_draws_other
    load params_hank_mode
    model = '\hank';
    model_name = 'HANK - RE';

    % B-RANK

    elseif i_model == 3
    load rank_draws_main_behav
    load rank_draws_other_behav
    load params_rank_mode_behav
    model = '\rank';
    model_name = 'RANK - CD';

    % B-HANK

    elseif i_model == 4
    load hank_draws_main_behav
    load hank_draws_other_behav
    load params_hank_mode_behav
    model = '\hank';
    model_name = 'HANK - CD';

    end

    % unpack settings

    N_adapt      = settings_MH.N_adapt;
    N_burn       = settings_MH.N_burn;
    N_keep       = settings_MH.N_keep;
    keep_every   = settings_MH.keep_every;
    N_draws      = settings_MH.N_draws;
    n_iter_print = settings_MH.n_iter_print;
    T_use        = settings_MH.T_use;
    T_save       = settings_MH.T_save;
    init_val     = settings_MH.init_val;
    var_scale    =settings_MH.var_scale;
    var_mat_draws_init = settings_MH.var_mat_draws_init;
    
    param_sd = std(param_collector(1:6,N_burn+1:end),0,2);
    param_uti = param_collector(7,N_burn+1:end)./(1-param_collector(7,N_burn+1:end));
    param_uti_sd = std(param_uti);
    param_mean = mean(param_collector(1:7,N_burn+1:end),2);

    % compute parameter moments
    
    start_2 = N_burn+1;
    end_2 = N_draws;
    param_quants = quantile(param_collector',[0.5 0.05 0.95]);
    param_means = mean(param_collector');

    % create table for estimated parameters
    
    % habit/info
    if strcmp(model,'\rank')
        fprintf(f,[model_name ' & $h$ & Beta & $%1.2f$ & $%1.2f$ & & $ %1.4f$ & $ %1.4f$ & $ %1.4f$& $ %1.4f$ & $ %1.4f$ \\\\ \n '],...
            [h_mean, h_sd, param_sol(1), param_means(1), param_quants(:,1)']);
    elseif strcmp(model,'\hank')
        fprintf(f,[model_name ' & $\\theta$ & Beta & $%1.2f$ & $%1.2f$ & & $ %1.4f$ & $ %1.4f$ & $ %1.4f$& $ %1.4f$& $ %1.4f$ \\\\ \n '],...
            [theta_mean, theta_sd, param_sol(1), param_means(1), param_quants(:,1)']);
    end

    % calvo_p
    fprintf(f,'& $\\theta_p$ & Beta & $%1.2f$ & $%1.2f$ & & $ %1.4f$ & $ %1.4f$ & $ %1.4f$& $ %1.4f$ & $ %1.4f$ \\\\ \n ',...
        [calvo_p_mean, calvo_p_sd, param_sol(2), param_means(2), param_quants(:,2)']);

    % calvo_w
    fprintf(f,'& $\\theta_w$ & Beta & $%1.2f$ & $%1.2f$ & & $ %1.4f$ & $ %1.4f$ & $ %1.4f$& $ %1.4f$ & $ %1.4f$ \\\\ \n ',...
        [calvo_w_mean, calvo_w_sd, param_sol(4), param_means(4), param_quants(:,4)']);

    % investment adjustment cost
    fprintf(f,'& $\\kappa$ & Normal & $%1.2f$ & $%1.2f$ & & $ %1.4f$ & $ %1.4f$ & $ %1.4f$& $ %1.4f$ & $ %1.4f$ \\\\ \n ',...
        [kappa_mean, kappa_sd, param_sol(6), param_means(6), param_quants(:,6)']);

    % capacity utilization
    fprintf(f,'& $\\psi$ & Beta & $%1.2f$ & $%1.2f$ & & $ %1.4f$ & $ %1.4f$ & $ %1.4f$& $ %1.4f$ & $ %1.4f$ \\\\ \\hline \n ',...
        [psi_uti_mean, psi_uti_sd, param_sol(7), param_means(7), param_quants(:,7)']);

end
fprintf(f,'\\end{tabular} \n');
fprintf(f,'\\end{table} \n');

fclose(f);