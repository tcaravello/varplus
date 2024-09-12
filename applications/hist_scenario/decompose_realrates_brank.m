%% DECOMPOSE REAL RATE CAUSAL EFFECTS
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

% load posterior mode

load rank_base_mode

Pi_m_rank = 4 * Pi_m_base_rank;
Y_m_rank  = Y_m_base_rank;
I_m_rank  = 4 * I_m_base_rank;

load rank_base_mode_behav

Pi_m_brank = 4 * Pi_m_base_rank;
Y_m_brank  = Y_m_base_rank;
I_m_brank  = 4 * I_m_base_rank;

clear Pi_m_base_rank Y_m_base_rank I_m_base_rank

% time horizon

T = size(Pi_m_rank,1);

% inflation horizon response

pi_hor_indx = 4;

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

set_cnfctl_rule

%% CONSTRUCT COUNTERFACTUAL HISTORICAL SCENARIO

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

pi_cnfctl   = NaN(T,2);
y_cnfctl    = NaN(T,2);
i_cnfctl    = NaN(T,2);
r_cnfctl    = NaN(T,2);
pi_r_cnfctl = NaN(T,2);

%----------------------------------------------------------------
% Counterfactuals
%----------------------------------------------------------------

for i_model = 1:2
    
% causal effect matrices

if i_model == 1

    Pi_m = Pi_m_rank;
    Y_m  = Y_m_rank;
    I_m  = I_m_rank;

elseif i_model == 2

    Pi_m = Pi_m_brank;
    Y_m  = Y_m_brank;
    I_m  = I_m_brank;

end

R_m  = I_m - [Pi_m(2:end,:); [0 Pi_m(end,1:end-1)]];
Pi_r = Pi_m / R_m;

if i_model == 1
    Pi_r_rank = Pi_r;
elseif i_model == 2
    Pi_r_brank = Pi_r;
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

[pi_cnfctl_tmp, y_cnfctl_tmp, i_cnfctl_tmp] ...
        = cnfctl_pred_fn(A_pi,A_y,A_i,wedge,Pi_m,Y_m,I_m,pi_x,y_x,i_x,1);

pi_cnfctl(:,i_model)  = pi_cnfctl_tmp;
y_cnfctl(:,i_model)   = y_cnfctl_tmp;
i_cnfctl(:,i_model)   = i_cnfctl_tmp;
r_cnfctl(:,i_model)   = (i_cnfctl_tmp - i_x) - [pi_cnfctl_tmp(2:end) - pi_x(2:end);0];
pi_r_cnfctl(:,i_model)   = (Pi_r(pi_hor_indx,:) .* r_cnfctl(:,i_model)')';

end

%% PLOT RESULTS

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% color settings

settings.colors.black  = [0 0 0];
settings.colors.grey   = [150/255 150/255 150/255];
settings.colors.orange = [204/255 102/255 0/255];
settings.colors.green  = [37/255 152/255 14/255];
settings.colors.beige  = [196/255 174/255 120/255];
settings.colors.blue   = [102/255 178/255 255/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.purple = [160/255 32/255 240/255];
settings.colors.list = [settings.colors.black;settings.colors.grey;settings.colors.orange];
settings.colors.models = [196/255 174/255 120/255; ... % beige
                            204/255 0/255 0/255; ... % red
                            102/255 178/255 255/255]; % blue

% figure spacing

plotwidth = 0.265;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

%----------------------------------------------------------------
% Interest Rate Decomposition Plot
%----------------------------------------------------------------

cd([path vintage experiment '/_results/RE']);

max_hor = 100;
xdates  = date_ext(fcst_date)+0.25*(0:max_hor-1);

x_min = min(xdates);
x_max = 2037;

figure

subplot(1,3,1)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
plot(xdates,100 * r_cnfctl(1:max_hor,1),'-','Color',settings.colors.beige,'LineWidth',4)
hold on
plot(xdates,100 * r_cnfctl(1:max_hor,2),'-','Color',settings.colors.blue,'LineWidth',4)
hold on
set(gcf,'color','w')
title('Change in Real Rate $r$','interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
xlim([x_min x_max])
ylabel('Basis Points','interpreter','latex','FontSize',20)
ylim([-100 125])
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
plot(xdates,Pi_r_rank(pi_hor_indx,1:max_hor),'-','Color',settings.colors.beige,'LineWidth',4)
hold on
plot(xdates,Pi_r_brank(pi_hor_indx,1:max_hor),'-','Color',settings.colors.blue,'LineWidth',4)
hold on
set(gcf,'color','w')
title('Response of $\pi_{2022Q2}$ to $r$','interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
xlim([x_min x_max])
ylabel('\% Deviation','interpreter','latex','FontSize',20)
ylim([-0.5 0.25])
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
jbfill(xdates,(zeros(max_hor,1))',100 * pi_r_cnfctl(1:max_hor,1)',...
    settings.colors.beige,settings.colors.beige,0,0.5);
hold on
jbfill(xdates,(zeros(max_hor,1))',100 * pi_r_cnfctl(1:max_hor,2)',...
    settings.colors.blue,settings.colors.blue,0,0.5);
hold on
set(gcf,'color','w')
title('Contribution of $r$ to $\pi_{2022Q2}$','interpreter','latex','fontsize',24)
xlabel('Date','interpreter','latex','FontSize',20)
xlim([x_min x_max])
ylabel('\% Deviation','interpreter','latex','FontSize',20)
ylim([-3 3])
legend({'RANK','B-RANK'},'Location','Northeast','fontsize',18,'interpreter','latex')
grid on
hold off

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
    print('decompose_pir','-dpng');
end

cd([path vintage experiment]);