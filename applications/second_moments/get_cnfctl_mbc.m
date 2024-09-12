%% COUNTERFACTUAL POLICY IRFs AFTER MBC SHOCK
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

experiment = '/applications/second_moments';

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
% MBC IRFs
%----------------------------------------------------------------

load mbc_results

mbc_base  = IS_MBC(1:T,[9 2 10]); % variables: (pi, y, i)

series_names = series_names([9 2 10]);

n_y      = size(mbc_base,2);

clear IS_MBC

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

%% CONSTRUCT MBC SHOCK COUNTERFACTUAL

%----------------------------------------------------------------
% Counterfactual Posterior Draws
%----------------------------------------------------------------

mbc_cnfctl = NaN(T,n_y,n_draws);

for i_draw = 1:n_draws
    
% get policy shock IRFs

Pi_m = Pi_m_draws(:,1:shock_max,i_draw);
Y_m  = Y_m_draws(:,1:shock_max,i_draw);
I_m  = I_m_draws(:,1:shock_max,i_draw);

% get the corresponding base sequences

pi_mbc = mbc_base(:,1);
y_mbc  = mbc_base(:,2);
i_mbc  = mbc_base(:,3);

% find best fit to counterfactual rule

if cnfctl_optpol == 0

[pi_mbc_cnfctl,y_mbc_cnfctl,i_mbc_cnfctl] = cnfctl_fn(A_pi,A_y,A_i,wedge,...
    Pi_m,Y_m,I_m,pi_mbc,y_mbc,i_mbc);

elseif cnfctl_optpol == 1

[pi_mbc_cnfctl,y_mbc_cnfctl,i_mbc_cnfctl] = optpol_fn(W_pi,W_y,W_i,...
    Pi_m,Y_m,I_m,pi_mbc,y_mbc,i_mbc);

end

mbc_cnfctl(:,1,i_draw) = pi_mbc_cnfctl;
mbc_cnfctl(:,2,i_draw) = y_mbc_cnfctl;
mbc_cnfctl(:,3,i_draw) = i_mbc_cnfctl;

end

%----------------------------------------------------------------
% Construct Percentiles
%----------------------------------------------------------------

mbc_cnfctl_lb  = quantile(mbc_cnfctl,0.16,3);
mbc_cnfctl_med = quantile(mbc_cnfctl,0.5,3);
mbc_cnfctl_ub  = quantile(mbc_cnfctl,0.84,3);

%----------------------------------------------------------------
% Individual Models
%----------------------------------------------------------------

n_models          = 2;
mbc_cnfctl_models = NaN(T,n_y,n_models);

for i_model = 1:n_models
    
pi_mbc_cfnctl_tmp = NaN(T,n_draws);
y_mbc_cfnctl_tmp  = NaN(T,n_draws);
i_mbc_cfnctl_tmp  = NaN(T,n_draws);
    
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

% get the corresponding base sequences

pi_mbc = mbc_base(:,1);
y_mbc  = mbc_base(:,2);
i_mbc  = mbc_base(:,3);

% find best fit to counterfactual rule

if cnfctl_optpol == 0

[pi_mbc_cnfctl,y_mbc_cnfctl,i_mbc_cnfctl] = cnfctl_fn(A_pi,A_y,A_i,wedge,...
    Pi_m,Y_m,I_m,pi_mbc,y_mbc,i_mbc);

elseif cnfctl_optpol == 1

[pi_mbc_cnfctl,y_mbc_cnfctl,i_mbc_cnfctl] = optpol_fn(W_pi,W_y,W_i,...
    Pi_m,Y_m,I_m,pi_mbc,y_mbc,i_mbc);

end

pi_mbc_cfnctl_tmp(:,i_draw) = pi_mbc_cnfctl;
y_mbc_cfnctl_tmp(:,i_draw)  = y_mbc_cnfctl;
i_mbc_cfnctl_tmp(:,i_draw)  = i_mbc_cnfctl;

end

mbc_cnfctl_models(:,1,i_model) = mean(pi_mbc_cfnctl_tmp,2);
mbc_cnfctl_models(:,2,i_model) = mean(y_mbc_cfnctl_tmp,2);
mbc_cnfctl_models(:,3,i_model) = mean(i_mbc_cfnctl_tmp,2);

end

clear pi_mbc_cfnctl_tmp y_mbc_cfnctl_tmp i_mbc_cfnctl_tmp

%% PLOT COUNTERFACTUAL IRFs

%----------------------------------------------------------------
% General Settings
%----------------------------------------------------------------

% plot horizon

IRF_hor_plot = 30;

% color settings

settings.colors.black  = [0 0 0];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.models = [196/255 174/255 120/255; ... % beige
                            204/255 0/255 0/255; ... % red
                            102/255 178/255 255/255]; % blue

% plot size

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

%----------------------------------------------------------------
% MBC Shock
%----------------------------------------------------------------

if indic_RE == 1
    cd([path vintage experiment '/_results/RE']);
elseif indic_behav == 1
    cd([path vintage experiment '/_results/behav']);
elseif indic_joint == 1
    cd([path vintage experiment '/_results/joint']);
end

figure

for i_y = 1:n_y

subplot(1,3,i_y)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(i_y);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
plot(0,0,':','Color',settings.colors.black,'LineWidth',4)
hold on
plot(0,0,'Color',settings.colors.blue,'LineWidth',4)
hold on
jbfill(0:1:IRF_hor_plot,(mbc_cnfctl_lb(1:IRF_hor_plot+1,i_y))',(mbc_cnfctl_ub(1:IRF_hor_plot+1,i_y))',...
    settings.colors.lblue,settings.colors.lblue,0,1);
hold on
plot(0:1:IRF_hor_plot,mbc_base(1:IRF_hor_plot+1,i_y),':','Color',settings.colors.black,'LineWidth',4)
hold on
plot(0:1:IRF_hor_plot,mbc_cnfctl_med(1:IRF_hor_plot+1,i_y),'Color',settings.colors.blue,'LineWidth',4)
hold on
% xlim([1 IRF_hor_plot])
% ylim([-2 2])
% yticks([-2 -1 0 1 2])
set(gcf,'color','w')
title(series_names(i_y),'interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
if i_y == 1
    ylabel('\% Deviation','interpreter','latex','FontSize',20)
end
if i_y == 2
    legend({'Data','Counterfct''l'},'Location','Southeast','fontsize',18,'interpreter','latex','NumColumns',2)
end
grid on
hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

if save_fig == 1
if indic_1shock == 0
    if cnfctl_optpol == 1
        print('mbc_optpol','-dpng');
    end
elseif indic_1shock == 1
    if cnfctl_optpol == 1
        print('mbc_optpol_1shock','-dpng');
    end
end
end

cd([path vintage experiment]);