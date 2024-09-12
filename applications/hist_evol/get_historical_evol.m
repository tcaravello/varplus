%% HISTORICAL EVOLUTION UNDER COUNTERFACTUAL POLICY RULE
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

experiment = '/applications/hist_evol';

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

load fcst_evol_results

fcst_date = find(date == dates_fcst(1));

pi_base_x = squeeze(var_forecasts_OLS(1:T,1,:));
y_base_x  = squeeze(var_forecasts_OLS(1:T,2,:));
i_base_x  = squeeze(var_forecasts_OLS(1:T,3,:));

pi_history = var_history(:,1);
y_history  = var_history(:,2);
i_history  = var_history(:,3);

clear var_forecasts_OLS var_history

%----------------------------------------------------------------
% Specify Counterfactual Rule
%----------------------------------------------------------------

set_cnfctl_rule

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

%% CONSTRUCT COUNTERFACTUAL HISTORICAL EVOLUTION

%----------------------------------------------------------------
% Settings & Placeholders
%----------------------------------------------------------------

fcst_hor = size(dates_fcst,2);
n_models = 0;

if indic_behav == 1
    model_names = {'B-RANK','B-HANK'};
else
    model_names = {'RANK','HANK'};
end

pi_cnfctl    = zeros(fcst_hor,n_draws);
y_cnfctl     = zeros(fcst_hor,n_draws);
i_cnfctl     = zeros(fcst_hor,n_draws);

pi_x_cnfctl_OLS  = zeros(fcst_hor,fcst_hor);
y_x_cnfctl_OLS   = zeros(fcst_hor,fcst_hor);
i_x_cnfctl_OLS   = zeros(fcst_hor,fcst_hor);

pi_cnfctl_models    = zeros(fcst_hor,n_models,n_draws);
y_cnfctl_models     = zeros(fcst_hor,n_models,n_draws);
i_cnfctl_models     = zeros(fcst_hor,n_models,n_draws);

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
    A_pi = Pi_m' * W_pi;
    A_y  = Y_m' * W_y;
    A_i  = I_m' * W_i;
end
    
% historical counterfactual: initial date (use the forecasted path and the lagged interest rate)

t_simul = 1;

pi_x  = pi_base_x(:,t_simul);
y_x   = y_base_x(:,t_simul);
i_x   = i_base_x(:,t_simul);

if cnfctl_optpol == 1
    wedge = lambda_di * I_m' * [i_history(fcst_date-1); zeros(T-1,1)];
elseif cnfctl_tylr == 1
    wedge = zeros(shock_max,1);
    wedge(1) = rho_ib * i_history(fcst_date-1);
else
    wedge = zeros(shock_max,1);
end

[pi_path, y_path, i_path] = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m,Y_m,I_m,pi_x,y_x,i_x,1);

pi_cnfctl(t_simul,i_draw) = pi_path(1);
y_cnfctl(t_simul,i_draw)  = y_path(1);
i_cnfctl(t_simul,i_draw)  = i_path(1);

% historical counterfactual: subsequent dates

for t_simul = 2:fcst_hor

    % get the forecast revisions
    pi_x  = pi_base_x(:,t_simul) - [pi_base_x(2:end,t_simul-1);0];
    y_x   = y_base_x(:,t_simul)  - [y_base_x(2:end,t_simul-1);0];
    i_x   = i_base_x(:,t_simul)  - [i_base_x(2:end,t_simul-1);0];
    
    wedge = zeros(shock_max,1);

    % construct counterfactual response to forecast revision
    [d_pi, d_y, d_i] = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m,Y_m,I_m,pi_x,y_x,i_x,1);

    pi_path = [pi_path(2:end);0] + d_pi;
    y_path  = [y_path(2:end);0] + d_y;
    i_path  = [i_path(2:end);0] + d_i;

    pi_cnfctl(t_simul,i_draw) = pi_path(1);
    y_cnfctl(t_simul,i_draw)  = y_path(1);
    i_cnfctl(t_simul,i_draw)  = i_path(1);

end

end

clear Pi_m Y_m I_m

pi_cnfctl_lb  = quantile(pi_cnfctl,0.16,2);
pi_cnfctl_med = quantile(pi_cnfctl,0.5,2);
pi_cnfctl_ub  = quantile(pi_cnfctl,0.84,2);

y_cnfctl_lb  = quantile(y_cnfctl,0.16,2);
y_cnfctl_med = quantile(y_cnfctl,0.5,2);
y_cnfctl_ub  = quantile(y_cnfctl,0.84,2);

i_cnfctl_lb  = quantile(i_cnfctl,0.16,2);
i_cnfctl_med = quantile(i_cnfctl,0.5,2);
i_cnfctl_ub  = quantile(i_cnfctl,0.84,2);
%----------------------------------------------------------------
% OLS Point Estimates
%----------------------------------------------------------------

% some preparations

% keep only the first "shock_max" shocks.

Pi_m_base = Pi_m_base(:,1:shock_max); 
Y_m_base = Y_m_base(:,1:shock_max);
I_m_base = I_m_base(:,1:shock_max,:);

if cnfctl_optpol == 1
    A_pi = Pi_m_base' * W_pi;
    A_y  = Y_m_base' * W_y;
    A_i  = I_m_base' * W_i;
end
    
m_cnfctl_aux = zeros(T,T);

% historical counterfactual: initial date (use the forecasted path and the lagged interest rate)

t_simul = 1;

pi_x  = pi_base_x(:,t_simul);
y_x   = y_base_x(:,t_simul);
i_x   = i_base_x(:,t_simul);

if cnfctl_optpol == 1
    wedge = lambda_di * I_m_base' * [i_history(fcst_date-1); zeros(T-1,1)];
elseif cnfctl_tylr == 1
    wedge = zeros(shock_max,1);
    wedge(1) = rho_ib * i_history(fcst_date-1);
else
    wedge = zeros(shock_max,1);
end

[pi_path, y_path, i_path] = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m_base,Y_m_base,I_m_base,pi_x,y_x,i_x,1);

pi_x_cnfctl_OLS(:,t_simul) = pi_path(1:fcst_hor);
y_x_cnfctl_OLS(:,t_simul)  = y_path(1:fcst_hor);
i_x_cnfctl_OLS(:,t_simul)  = i_path(1:fcst_hor);

% historical counterfactual: subsequent dates

for t_simul = 2:fcst_hor

    % get the forecast revisions
    pi_x  = pi_base_x(:,t_simul) - [pi_base_x(2:end,t_simul-1);0];
    y_x   = y_base_x(:,t_simul)  - [y_base_x(2:end,t_simul-1);0];
    i_x   = i_base_x(:,t_simul)  - [i_base_x(2:end,t_simul-1);0];
    wedge = zeros(shock_max,1);

   % construct counterfactual response to forecast revision
    [d_pi, d_y, d_i] = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m_base,Y_m_base,I_m_base,pi_x,y_x,i_x,1);

    pi_path = [pi_path(2:end);0] + d_pi;
    y_path = [y_path(2:end);0] + d_y;
    i_path = [i_path(2:end);0] + d_i;

    pi_x_cnfctl_OLS(:,t_simul) = pi_path(1:fcst_hor);
    y_x_cnfctl_OLS(:,t_simul)  = y_path(1:fcst_hor);
    i_x_cnfctl_OLS(:,t_simul)  = i_path(1:fcst_hor);

end

%----------------------------------------------------------------
% Individual Models
%----------------------------------------------------------------

if n_models>0

for i_model = 1:n_models
    
for i_draw = 1:n_draws
    
% get policy shock IRFs

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
    A_pi = Pi_m' * W_pi;
    A_y  = Y_m' * W_y;
    A_i  = I_m' * W_i;
end
    
% historical counterfactual: initial date (use the forecasted path and the lagged interest rate)

t_simul = 1;

pi_x  = pi_base_x(:,t_simul);
y_x   = y_base_x(:,t_simul);
i_x   = i_base_x(:,t_simul);

if cnfctl_optpol == 1
    wedge = lambda_di * I_m' * [i_history(fcst_date-1); zeros(T-1,1)];
elseif cnfctl_tylr == 1
    wedge = zeros(shock_max,1);
    wedge(1) = rho_ib * i_history(fcst_date-1);
else
    wedge = zeros(shock_max,1);
end

[pi_path, y_path, i_path] = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m,Y_m,I_m,pi_x,y_x,i_x,1);

pi_cnfctl_models(t_simul,i_model,i_draw) = pi_path(1);
y_cnfctl_models(t_simul,i_model,i_draw)  = y_path(1);
i_cnfctl_models(t_simul,i_model,i_draw)  = i_path(1);

% historical counterfactual: subsequent dates

for t_simul = 2:fcst_hor

    % get the forecast revisions
    pi_x  = pi_base_x(:,t_simul) - [pi_base_x(2:end,t_simul-1);0];
    y_x   = y_base_x(:,t_simul)  - [y_base_x(2:end,t_simul-1);0];
    i_x   = i_base_x(:,t_simul)  - [i_base_x(2:end,t_simul-1);0];
    
    wedge = zeros(shock_max,1);

    % construct counterfactual response to forecast revision
    [d_pi, d_y, d_i] = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m,Y_m,I_m,pi_x,y_x,i_x,1);

    pi_path = [pi_path(2:end);0] + d_pi;
    y_path  = [y_path(2:end);0] + d_y;
    i_path  = [i_path(2:end);0] + d_i;

    pi_cnfctl_models(t_simul,i_model,i_draw) = pi_path(1);
    y_cnfctl_models(t_simul,i_model,i_draw)  = y_path(1);
    i_cnfctl_models(t_simul,i_model,i_draw)  = i_path(1);

end

end

end

clear Pi_m Y_m I_m

pi_cnfctl_models = quantile(pi_cnfctl_models,0.5,3);
y_cnfctl_models  = quantile(y_cnfctl_models,0.5,3);
i_cnfctl_models  = quantile(i_cnfctl_models,0.5,3);
end

%% PLOT RESULTS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% color settings

settings.colors.black  = [0 0 0];
settings.colors.grey   = [150/255 150/255 150/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.green = [37/255 152/255 14/255];
settings.colors.navyblue = [0/255 0/255 50/255];
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

% show interest rate forecasts at each date?

show_fcst = 0;

%----------------------------------------------------------------
% Add History to Forecasts
%----------------------------------------------------------------

% inflation

pi_cnfctl_lb  = [pi_history(fcst_date-fcst_lag:fcst_date-1);pi_cnfctl_lb];
pi_cnfctl_med = [pi_history(fcst_date-fcst_lag:fcst_date-1);pi_cnfctl_med];
pi_cnfctl_ub  = [pi_history(fcst_date-fcst_lag:fcst_date-1);pi_cnfctl_ub];
if n_models>0
    pi_cnfctl_models = [repmat(pi_history(fcst_date-fcst_lag:fcst_date-1),1,n_models);pi_cnfctl_models];
end
pi_cnfctl_lb  = pi_cnfctl_lb + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,9);
pi_cnfctl_med = pi_cnfctl_med + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,9);
pi_cnfctl_ub  = pi_cnfctl_ub + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,9);
if n_models>0
    pi_cnfctl_models = pi_cnfctl_models + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,9);
end

pi_history = pi_history + det_X * det_coeff(:,9);

% output

y_cnfctl_lb  = [y_history(fcst_date-fcst_lag:fcst_date-1);y_cnfctl_lb];
y_cnfctl_med = [y_history(fcst_date-fcst_lag:fcst_date-1);y_cnfctl_med];
y_cnfctl_ub  = [y_history(fcst_date-fcst_lag:fcst_date-1);y_cnfctl_ub];
if n_models>0
    y_cnfctl_models = [repmat(y_history(fcst_date-fcst_lag:fcst_date-1),1,n_models);y_cnfctl_models];
end
y_cnfctl_lb  = y_cnfctl_lb + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,2);
y_cnfctl_med = y_cnfctl_med + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,2);
y_cnfctl_ub  = y_cnfctl_ub + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,2);
if n_models>0
    y_cnfctl_models = y_cnfctl_models + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,2);
end

y_history = y_history + det_X * det_coeff(:,2);

% interest rates

i_cnfctl_lb  = [i_history(fcst_date-fcst_lag:fcst_date-1);i_cnfctl_lb];
i_cnfctl_med = [i_history(fcst_date-fcst_lag:fcst_date-1);i_cnfctl_med];
i_cnfctl_ub  = [i_history(fcst_date-fcst_lag:fcst_date-1);i_cnfctl_ub];
if n_models>0
    i_cnfctl_models = [repmat(i_history(fcst_date-fcst_lag:fcst_date-1),1,n_models);i_cnfctl_models];
end
i_cnfctl_lb  = i_cnfctl_lb  + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,10);
i_cnfctl_med = i_cnfctl_med + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,10);
i_cnfctl_ub  = i_cnfctl_ub  + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,10);
if n_models>0
    i_cnfctl_models = i_cnfctl_models + det_X(fcst_date-fcst_lag:fcst_date+fcst_hor-1,:) * det_coeff(:,10);
end

i_history = i_history + det_X * det_coeff(:,10);

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
plot(0,0,'-','Color',settings.colors.blue,'LineWidth',4)
hold on
jbfill(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),(pi_cnfctl_lb)',(pi_cnfctl_ub)',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),pi_cnfctl_med,'-','Color',settings.colors.blue,'LineWidth',4)
hold on
plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),pi_history(fcst_date-fcst_lag:fcst_date+fcst_hor-1),'-','Color',settings.colors.black,'LineWidth',4)
hold on
if n_models>0
for i_model = 1:n_models
    plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),pi_cnfctl_models(:,i_model),'-','Color',settings.colors.models(i_model,:),'LineWidth',4)
    hold on
end
end
xlim([date(fcst_date-fcst_lag),date(fcst_date+fcst_hor-1)])
set(gcf,'color','w')
title(series_names(1),'interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
% ylabel('\% Deviation','interpreter','latex','FontSize',20)
legend('History','Counterfct''l','Location','Southwest','fontsize',18,'interpreter','latex')
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
jbfill(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),(y_cnfctl_lb)',(y_cnfctl_ub)',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),y_cnfctl_med,'-','Color',settings.colors.blue,'LineWidth',4)
hold on
plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),y_history(fcst_date-fcst_lag:fcst_date+fcst_hor-1),'-','Color',settings.colors.black,'LineWidth',4)
hold on
if n_models>0
for i_model = 1:n_models
    plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),y_cnfctl_models(:,i_model),'-','Color',settings.colors.models(i_model,:),'LineWidth',4)
    hold on
end
end
xlim([date(fcst_date-fcst_lag),date(fcst_date+fcst_hor-1)])
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
jbfill(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),(i_cnfctl_lb)',(i_cnfctl_ub)',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),i_cnfctl_med,'-','Color',settings.colors.blue,'LineWidth',4)
% plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),i_cnfctl_OLS,'-','Color',settings.colors.black,'LineWidth',4)
hold on
plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),i_history(fcst_date-fcst_lag:fcst_date+fcst_hor-1),'-','Color',settings.colors.black,'LineWidth',4)
hold on
if n_models>0
for i_model = 1:n_models
    plot(date(fcst_date-fcst_lag):0.25:date(fcst_date+fcst_hor-1),i_cnfctl_models(:,i_model),'-','Color',settings.colors.models(i_model,:),'LineWidth',4)
    hold on
end
end
if show_fcst == 1
for i_fcst = 1:fcst_hor
    plot(date(fcst_date+i_fcst-1):0.25:date(fcst_date+fcst_hor-1),i_fcst_cnfctl_OLS(1:end-i_fcst+1,i_fcst),'--','Color',settings.colors.black,'LineWidth',1)
    hold on
end
end
hold on
xlim([date(fcst_date-fcst_lag),date(fcst_date+fcst_hor-1)])
set(gcf,'color','w')
title(series_names(3),'interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
% ylabel('\% Deviation','interpreter','latex','FontSize',20)
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [0.5*pos(1) pos(2) 2.25*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

if save_fig == 1
if indic_1shock == 0
    if cnfctl_optpol == 1
        print('evol_gr_optpol','-dpng');
    elseif cnfctl_tylr == 1
        print('evol_gr_tylr','-dpng');
    end
elseif indic_1shock == 1
    if cnfctl_optpol == 1
        print('evol_gr_optpol_1shock','-dpng');
    end
end
end

cd([path vintage experiment]);