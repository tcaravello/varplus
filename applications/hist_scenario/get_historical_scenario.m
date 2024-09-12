%% HISTORICAL SCENARIO UNDER COUNTERFACTUAL POLICY RULE
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

experiment = '/applications/hist_scenario';

save_fig = 1;

addpath([path vintage '/suff_stats/ratex']);
addpath([path vintage '/suff_stats/behavioral']);
addpath([path vintage '/suff_stats/mix']);
addpath([path vintage '/_auxiliary_functions'])
addpath([path vintage '/var_inputs/_results']);

cd([path vintage experiment]);

%% IMPORTS & SETTINGS

%----------------------------------------------------------------
% Policy Shock Sufficient Statistics
%----------------------------------------------------------------

% import

if indic_1shock == 0

    import_suffstats

elseif indic_1shock == 1

    import_suffstats_1shock

end

% sizes

T       = size(Pi_m_draws,1);
n_draws = size(Pi_m_draws,3);

%----------------------------------------------------------------
% Historical Evolution Inputs
%----------------------------------------------------------------

load fcst_scenario_results

pi_base_x = var_forecasts_OLS(1:T,1);
y_base_x  = var_forecasts_OLS(1:T,2);
i_base_x  = var_forecasts_OLS(1:T,3);
r_base_x  = i_base_x - [pi_base_x(2:T);0];

pi_history = var_history(:,1);
y_history  = var_history(:,2);
i_history  = var_history(:,3);
r_history  = i_history - [pi_history(2:end);0];

clear var_forecasts_OLS var_history

%----------------------------------------------------------------
% Specify Counterfactual Rule
%----------------------------------------------------------------

set_cnfctl_rule;

%----------------------------------------------------------------
% Shock Space
%----------------------------------------------------------------

if indic_1shock == 0
    shock_max = T; % set = T for all shocks
elseif indic_1shock == 1
    shock_max = 1;
end

if shock_max < T
    disp('Note: I am not using all shocks.')
end

%% CONSTRUCT COUNTERFACTUAL HISTORICAL SCENARIO

%----------------------------------------------------------------
% Settings & Placeholders
%----------------------------------------------------------------

pi_cnfctl    = zeros(fcst_hor,n_draws);
y_cnfctl     = zeros(fcst_hor,n_draws);
i_cnfctl     = zeros(fcst_hor,n_draws);
r_cnfctl     = zeros(fcst_hor,n_draws);
m_cnfctl     = zeros(fcst_hor,n_draws);

n_models = 2;

if indic_behav == 1
    model_names = {'B-RANK','B-HANK'};
else
    model_names = {'RANK','HANK'};
end

pi_cnfctl_models = zeros(fcst_hor,n_models,n_draws);
y_cnfctl_models  = zeros(fcst_hor,n_models,n_draws);
i_cnfctl_models  = zeros(fcst_hor,n_models,n_draws);
r_cnfctl_models  = zeros(fcst_hor,n_models,n_draws);
m_cnfctl_models  = zeros(fcst_hor,n_models,n_draws);

%----------------------------------------------------------------
% Posterior Draws
%----------------------------------------------------------------

for i_draw = 1:n_draws
    
% causal effects
    
Pi_m = Pi_m_draws(:,1:shock_max,i_draw);
Y_m  = Y_m_draws(:,1:shock_max,i_draw);
I_m  = I_m_draws(:,1:shock_max,i_draw);

% some preparations

if cnfctl_optpol == 1
    A_pi  = Pi_m' * W_pi;
    A_y   = Y_m' * W_y;
    A_i   = I_m' * W_i;
    wedge = lambda_di * I_m' * [i_history(end); zeros(T-1,1)];
elseif cnfctl_tylr == 1
    wedge(1) = rho_ib * i_history(end);
end

% "wedge" term accounts for initial conditions if the rule includes lags.
% for example, the Taylor rule includes lagged interest rates.

% construct the historical counterfactual

pi_x = pi_base_x;
y_x  = y_base_x;
i_x  = i_base_x;

[pi_cnfctl_tmp, y_cnfctl_tmp, i_cnfctl_tmp, m_cnfctl_tmp] ...
        = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m,Y_m,I_m,pi_x,y_x,i_x,1);

pi_cnfctl(:,i_draw) = pi_cnfctl_tmp(1:fcst_hor);
y_cnfctl(:,i_draw)  = y_cnfctl_tmp(1:fcst_hor);
i_cnfctl(:,i_draw)  = i_cnfctl_tmp(1:fcst_hor);
r_cnfctl(:,i_draw)  = i_cnfctl_tmp(1:fcst_hor) - pi_cnfctl_tmp(2:fcst_hor+1);

if indic_1shock == 1
    m_cnfctl(1:shock_max,i_draw) = m_cnfctl_tmp;
else
    m_cnfctl(:,i_draw)  = m_cnfctl_tmp(1:fcst_hor);
end

clear pi_cnfctl_tmp y_cnfctl_tmp i_cnfctl_tmp m_cnfctl_tmp

end

% percentiles

pi_cnfctl_lb  = quantile(pi_cnfctl,0.16,2);
pi_cnfctl_med = quantile(pi_cnfctl,0.5,2);
pi_cnfctl_ub  = quantile(pi_cnfctl,0.84,2);

y_cnfctl_lb  = quantile(y_cnfctl,0.16,2);
y_cnfctl_med = quantile(y_cnfctl,0.5,2);
y_cnfctl_ub  = quantile(y_cnfctl,0.84,2);

i_cnfctl_lb  = quantile(i_cnfctl,0.16,2);
i_cnfctl_med = quantile(i_cnfctl,0.5,2);
i_cnfctl_ub  = quantile(i_cnfctl,0.84,2);

r_cnfctl_lb  = quantile(r_cnfctl,0.16,2);
r_cnfctl_med = quantile(r_cnfctl,0.5,2);
r_cnfctl_ub  = quantile(r_cnfctl,0.84,2);

%----------------------------------------------------------------
% Individual Models
%----------------------------------------------------------------

for i_model = 1:n_models
    
for i_draw = 1:n_draws
    
% causal effects
    
if i_model == 1

    Pi_m = Pi_m_rank_draws(:,1:shock_max,i_draw);
    Y_m  = Y_m_rank_draws(:,1:shock_max,i_draw);
    I_m  = I_m_rank_draws(:,1:shock_max,i_draw);

elseif i_model == 2

    Pi_m = Pi_m_hank_draws(:,1:shock_max,i_draw);
    Y_m  = Y_m_hank_draws(:,1:shock_max,i_draw);
    I_m  = I_m_hank_draws(:,1:shock_max,i_draw);

end

% some preparations

if cnfctl_optpol == 1
    A_pi  = Pi_m' * W_pi;
    A_y   = Y_m' * W_y;
    A_i   = I_m' * W_i;
    wedge = lambda_di * I_m' * [i_history(end); zeros(T-1,1)];
elseif cnfctl_tylr == 1
    wedge(1) = rho_ib * i_history(end);
end

% construct the historical counterfactual

pi_x = pi_base_x;
y_x  = y_base_x;
i_x  = i_base_x;

[pi_cnfctl_tmp, y_cnfctl_tmp, i_cnfctl_tmp, m_cnfctl_tmp] ...
        = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m,Y_m,I_m,pi_x,y_x,i_x,1);

pi_cnfctl_models(:,i_model,i_draw) = pi_cnfctl_tmp(1:fcst_hor);
y_cnfctl_models(:,i_model,i_draw)  = y_cnfctl_tmp(1:fcst_hor);
i_cnfctl_models(:,i_model,i_draw)  = i_cnfctl_tmp(1:fcst_hor);
r_cnfctl_models(:,i_model,i_draw)  = i_cnfctl_tmp(1:fcst_hor) - pi_cnfctl_tmp(2:fcst_hor+1);

if indic_1shock == 1
    m_cnfctl_models(1:shock_max,i_model,i_draw) = m_cnfctl_tmp;
else
    m_cnfctl_models(:,i_model,i_draw)  = m_cnfctl_tmp(1:fcst_hor);
end

clear pi_cnfctl_tmp y_cnfctl_tmp i_cnfctl_tmp m_cnfctl_tmp

end

end

pi_cnfctl_models = quantile(pi_cnfctl_models,0.5,3);
y_cnfctl_models  = quantile(y_cnfctl_models,0.5,3);
i_cnfctl_models  = quantile(i_cnfctl_models,0.5,3);
r_cnfctl_models  = quantile(r_cnfctl_models,0.5,3);
m_cnfctl_models  = quantile(m_cnfctl_models,0.5,3);

%% PLOT RESULTS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% color settings

settings.colors.black  = [0 0 0];
settings.colors.grey   = [150/255 150/255 150/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.list = [settings.colors.black;settings.colors.grey;settings.colors.orange];
settings.colors.models = [196/255 174/255 120/255; ... % beige
                            204/255 0/255 0/255; ... % red
                            102/255 178/255 255/255;...% blue
                            32/255, 119/255, 34/255 ]; %orange

% figure spacing

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

%----------------------------------------------------------------
% Add History to Forecasts
%----------------------------------------------------------------

% inflation

pi_cnfctl_lb     = [pi_history(fcst_date-fcst_lag:fcst_date);pi_cnfctl_lb];
pi_cnfctl_med    = [pi_history(fcst_date-fcst_lag:fcst_date);pi_cnfctl_med];
pi_cnfctl_ub     = [pi_history(fcst_date-fcst_lag:fcst_date);pi_cnfctl_ub];
pi_cnfctl_models = [repmat(pi_history(fcst_date-fcst_lag:fcst_date),1,n_models);pi_cnfctl_models];
pi_base_var      = [pi_history(fcst_date-fcst_lag:fcst_date);pi_base_x(1:fcst_hor)];

pi_cnfctl_lb     = pi_cnfctl_lb + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);
pi_cnfctl_med    = pi_cnfctl_med + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);
pi_cnfctl_ub     = pi_cnfctl_ub + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);
pi_cnfctl_models = pi_cnfctl_models + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);
pi_base_var      = pi_base_var + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,pi_pos);
pi_hist_plot     = [pi_base_var(1:fcst_lag+1);NaN(fcst_hor,1)];

% output

y_cnfctl_lb     = [y_history(fcst_date-fcst_lag:fcst_date);y_cnfctl_lb];
y_cnfctl_med    = [y_history(fcst_date-fcst_lag:fcst_date);y_cnfctl_med];
y_cnfctl_ub     = [y_history(fcst_date-fcst_lag:fcst_date);y_cnfctl_ub];
y_cnfctl_models = [repmat(y_history(fcst_date-fcst_lag:fcst_date),1,n_models);y_cnfctl_models];
y_base_var      = [y_history(fcst_date-fcst_lag:fcst_date);y_base_x(1:fcst_hor)];

y_cnfctl_lb     = y_cnfctl_lb + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);
y_cnfctl_med    = y_cnfctl_med + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);
y_cnfctl_ub     = y_cnfctl_ub + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);
y_cnfctl_models = y_cnfctl_models + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);
y_base_var      = y_base_var + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,y_pos);
y_hist_plot     = [y_base_var(1:fcst_lag+1);NaN(fcst_hor,1)];

% nominal interest rates

i_cnfctl_lb     = [i_history(fcst_date-fcst_lag:fcst_date);i_cnfctl_lb];
i_cnfctl_med    = [i_history(fcst_date-fcst_lag:fcst_date);i_cnfctl_med];
i_cnfctl_ub     = [i_history(fcst_date-fcst_lag:fcst_date);i_cnfctl_ub];
i_cnfctl_models = [repmat(i_history(fcst_date-fcst_lag:fcst_date),1,n_models);i_cnfctl_models];
i_base_var      = [i_history(fcst_date-fcst_lag:fcst_date);i_base_x(1:fcst_hor)];

i_cnfctl_lb     = i_cnfctl_lb + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);
i_cnfctl_med    = i_cnfctl_med + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);
i_cnfctl_ub     = i_cnfctl_ub + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);
i_cnfctl_models = i_cnfctl_models + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);
i_base_var      = i_base_var + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos);
i_hist_plot     = [i_base_var(1:fcst_lag+1);NaN(fcst_hor,1)];

% real interest rates

r_cnfctl_lb     = [r_history(fcst_date-fcst_lag:fcst_date);r_cnfctl_lb];
r_cnfctl_med    = [r_history(fcst_date-fcst_lag:fcst_date);r_cnfctl_med];
r_cnfctl_ub     = [r_history(fcst_date-fcst_lag:fcst_date);r_cnfctl_ub];
r_cnfctl_models = [repmat(r_history(fcst_date-fcst_lag:fcst_date),1,n_models);r_cnfctl_models];
r_base_var      = [r_history(fcst_date-fcst_lag:fcst_date);r_base_x(1:fcst_hor)];

r_cnfctl_lb     = r_cnfctl_lb + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos) ...
                    - [det_X_ext(fcst_date-fcst_lag+1:fcst_date+fcst_hor,:); det_X_ext(end,:)] * det_coeff(:,pi_pos);
r_cnfctl_med    = r_cnfctl_med + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos) ...
                    - [det_X_ext(fcst_date-fcst_lag+1:fcst_date+fcst_hor,:); det_X_ext(end,:)] * det_coeff(:,pi_pos);
r_cnfctl_ub     = r_cnfctl_ub + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos) ...
                    - [det_X_ext(fcst_date-fcst_lag+1:fcst_date+fcst_hor,:); det_X_ext(end,:)] * det_coeff(:,pi_pos);
r_cnfctl_models = r_cnfctl_models + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos) ...
                    - [det_X_ext(fcst_date-fcst_lag+1:fcst_date+fcst_hor,:); det_X_ext(end,:)] * det_coeff(:,pi_pos);
r_base_var      = r_base_var + det_X_ext(fcst_date-fcst_lag:fcst_date+fcst_hor,:) * det_coeff(:,i_pos) ...
                    - [det_X_ext(fcst_date-fcst_lag+1:fcst_date+fcst_hor,:); det_X_ext(end,:)] * det_coeff(:,pi_pos);
r_hist_plot     = [r_base_var(1:fcst_lag+1);NaN(fcst_hor,1)];

%----------------------------------------------------------------
% Historical Scenario Plot
%----------------------------------------------------------------

if indic_RE == 1
    cd([path vintage experiment '/_results/RE']);
elseif indic_behav == 1
    cd([path vintage experiment '/_results/behav']);
elseif indic_joint == 1
    cd([path vintage experiment '/_results/joint']);
end

figure

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
plot(0,0,'-','Color',settings.colors.black,'LineWidth',4)
hold on
plot(0,0,'-','Color',settings.colors.grey,'LineWidth',4)
hold on
plot(0,0,'-','Color',settings.colors.blue,'LineWidth',4)
hold on
jbfill(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),(pi_cnfctl_lb)',(pi_cnfctl_ub)',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),pi_cnfctl_med,'-','Color',settings.colors.blue,'LineWidth',4)
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),pi_base_var,'--','Color',settings.colors.grey,'LineWidth',4)
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),pi_hist_plot,'-','Color',settings.colors.black,'LineWidth',4)
hold on
xlim([date_ext(fcst_date-fcst_lag),date_ext(fcst_date+fcst_hor)])
set(gcf,'color','w')
title(series_names(1),'interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
% ylabel('\% Deviation','interpreter','latex','FontSize',20)
legend({'History','Forecast','Counterfct''l'},'Location','Southeast','fontsize',18,'interpreter','latex')
grid on
hold off

subplot(1,3,2)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),(y_cnfctl_lb)',(y_cnfctl_ub)',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),y_cnfctl_med,'-','Color',settings.colors.blue,'LineWidth',4)
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),y_base_var,'--','Color',settings.colors.grey,'LineWidth',4)
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),y_hist_plot,'-','Color',settings.colors.black,'LineWidth',4)
hold on
xlim([date_ext(fcst_date-fcst_lag),date_ext(fcst_date+fcst_hor)])
set(gcf,'color','w')
title(series_names(2),'interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
% ylabel('\% Deviation','interpreter','latex','FontSize',20)
grid on
hold off

subplot(1,3,3)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(3);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
jbfill(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),(i_cnfctl_lb)',(i_cnfctl_ub)',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),i_cnfctl_med,'-','Color',settings.colors.blue,'LineWidth',4)
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),i_base_var,'--','Color',settings.colors.grey,'LineWidth',4)
hold on
plot(date_ext(fcst_date-fcst_lag):0.25:date_ext(fcst_date+fcst_hor),i_hist_plot,'-','Color',settings.colors.black,'LineWidth',4)
hold on
xlim([date_ext(fcst_date-fcst_lag),date_ext(fcst_date+fcst_hor)])
set(gcf,'color','w')
title(series_names(3),'interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
% ylabel('\% Deviation','interpreter','latex','FontSize',20)
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [0.5*pos(1) 0.5*pos(2) 2.25*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

if save_fig == 1
if indic_1shock == 0
    if cnfctl_optpol == 1
        print('scenario_2021_optpol','-dpng');
    elseif cnfctl_tylr == 1
        print('scenario_2021_tylr','-dpng');
    end
elseif indic_1shock == 1
    if cnfctl_optpol == 1
        print('scenario_2021_optpol_1shock','-dpng');
    end
end
end
cd([path vintage experiment]);