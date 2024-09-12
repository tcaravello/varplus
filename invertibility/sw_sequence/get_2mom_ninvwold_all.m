%% PREDICT SECOND MOMENTS IN SMETS-WOUTERS
% Tomas Caravello, Alisdair McKay, and Christian Wolf
% this version: 09/03/2024

%% HOUSEKEEPING

task = '/invertibility';

save_fig = 1;

addpath([path vintage task '/_auxiliary_functions'])
addpath([path vintage task '/sw_sequence/_subroutines'])
addpath([path vintage task '/sw_sequence/_get_mp_irfs'])

cd([path vintage task '/sw_sequence']);

%% SETTINGS

%----------------------------------------------------------------
% IS Characterization
%----------------------------------------------------------------

settings.VMA_hor          = 500; % maximal horizon for VMA representation
settings.VAR_poplaglength = 250; % lag length in population VAR
settings.IRF_hor          = 30;  % horizon for computed IRFs
settings.IS_draws         = 200; % number of draws from the identified set
settings.weight_hor       = settings.VMA_hor; % horizon for weight polynomial P(L)

%----------------------------------------------------------------
% Colors
%----------------------------------------------------------------

settings.colors.black  = [0 0 0];
settings.colors.grey   = [205/255 205/255 205/255];
settings.colors.orange = [204/255 102/255 0/255];

settings.colors.list = [settings.colors.black;settings.colors.grey;settings.colors.orange];

%----------------------------------------------------------------
% List of Observables
%----------------------------------------------------------------

n_obs    = 4;
obs_list = cell(n_obs,1);

obs_list{1} = [5 4 19]; % (r,y,pi)
obs_list{2} = [5 4 19 21]; % (r,y,pi,lab)
obs_list{3} = [5 4 19 17 18]; % (r,y,pi,c,i)
obs_list{4} = [5 4 19 17 18 20 21]; % (r,y,pi,c,i,w,lab)

%----------------------------------------------------------------
% Placeholders
%----------------------------------------------------------------

omega_lb = (2*pi)/32;
omega_ub = (2*pi)/6;

n_omega = 1e2;
omega_grid = linspace(omega_lb,omega_ub,n_omega);

var_pred_all = NaN(3,3,n_obs);
s_y_pred_all = NaN(n_omega,3,n_obs);

%% PREDICT SECOND MOMENTS

for i_obs = 1:n_obs

%----------------------------------------------------------------
% Solve SW Model
%----------------------------------------------------------------

% model run

pause(3) % use this pause to avoid dynare error when overwriting files

dynare SW_Model noclearall
clean_folder_SW

% get law of motion for all model variables

SW_model.decision = decision(2:end,:);

% specify observables

SW_model.obs = obs_list{i_obs};

% size indicators

SW_model.n_y   = size(SW_model.obs,2);
SW_model.n_eps = M_.exo_nbr;
SW_model.n_s   = M_.nspred;

% ABCD representations

SW_model.ABCD = ABCD_fun_SW(SW_model);

% delete superfluous variables

clean_workspace_SW

%----------------------------------------------------------------
% Get Population IRFs + FVDs + Shock Sequences
%----------------------------------------------------------------

[SW_model.IRF,SW_model.M,SW_model.tot_weights] = pop_analysis(SW_model,settings);

%----------------------------------------------------------------
% Get VAR Representation
%----------------------------------------------------------------

VAR = popVAR(SW_model,settings);

%----------------------------------------------------------------
% Prepare Inputs
%----------------------------------------------------------------

% monetary policy causal effects

load sw_mp_irfs

% counterfactual rule

crpi = 1.5;
crr  = 0.8;
cry  = 0.5;
crdy = 0;

A_pi = zeros(settings.VMA_hor,settings.VMA_hor);
for t = 1:settings.VMA_hor
    A_pi(t,t) = -(1-crr) * crpi;
end

A_y = zeros(settings.VMA_hor,settings.VMA_hor);
for t = 1:settings.VMA_hor
    A_y(t,t) = -(1-crr) * (cry + crdy);
    if t > 1
        A_y(t,t-1) = (1-crr) * crdy;
    end
end

A_i = zeros(settings.VMA_hor,settings.VMA_hor);
for t = 1:settings.VMA_hor
    A_i(t,t) = 1;
    if t > 1
        A_i(t,t-1) = -crr;
    end
end

% adjust Wold IRFs

VAR.Wold_rot = chol(VAR.Sigma_u,'lower');
for t = 1:settings.VMA_hor
    VAR.IRF_Wold(:,:,t) = VAR.IRF_Wold(:,:,t) * VAR.Wold_rot;
end
VAR.IRF_Wold = permute(VAR.IRF_Wold,[3 1 2]);

%----------------------------------------------------------------
% Get Counterfactual IRFs
%----------------------------------------------------------------

% rotate the Wold innovations

VAR.IRF_Wold_cnfctl = NaN(settings.VMA_hor,3,SW_model.n_y);

for i_y = 1:SW_model.n_y
    m_seq_aux = -(A_pi * Pi_m + A_y * Y_m + A_i * I_m)^(-1) * (A_pi * VAR.IRF_Wold(:,3,i_y) ...
        + A_y * VAR.IRF_Wold(:,2,i_y) + A_i * VAR.IRF_Wold(:,1,i_y));
    VAR.IRF_Wold_cnfctl(:,3,i_y) = VAR.IRF_Wold(:,3,i_y) + Pi_m * m_seq_aux;
    VAR.IRF_Wold_cnfctl(:,2,i_y) = VAR.IRF_Wold(:,2,i_y) + Y_m * m_seq_aux;
    VAR.IRF_Wold_cnfctl(:,1,i_y) = VAR.IRF_Wold(:,1,i_y) + I_m * m_seq_aux;
end

% add a policy shock

m_star = zeros(settings.VMA_hor,1);
m_star(1) = -0.2290;
for t = 2:settings.VMA_hor
    m_star(t) = 0.1999 * m_star(t-1);
end
m_seq_aux = -(A_pi * Pi_m + A_y * Y_m + A_i * I_m)^(-1) * m_star;
VAR.IRF_Wold_cnfctl(:,3,SW_model.n_y+1) = Pi_m * m_seq_aux;
VAR.IRF_Wold_cnfctl(:,2,SW_model.n_y+1) = Y_m * m_seq_aux;
VAR.IRF_Wold_cnfctl(:,1,SW_model.n_y+1) = I_m * m_seq_aux;

%----------------------------------------------------------------
% Predicted Second Moments
%----------------------------------------------------------------

% variances

var_pred = zeros(3,3);
for t = 1:settings.VMA_hor
    var_pred = var_pred + squeeze(VAR.IRF_Wold_cnfctl(t,:,:)) * squeeze(VAR.IRF_Wold_cnfctl(t,:,:))';
end

var_pred_all(:,:,i_obs) = var_pred;

% spectral density

s_y_fn = @(omega) 1/(2*pi) * sd_fun(omega,VAR.IRF_Wold_cnfctl);

s_y_pred = NaN(n_omega,3);
for i_omega = 1:n_omega
    omega_val = omega_grid(i_omega);
    s_y_val = s_y_fn(omega_val);
    s_y_pred(i_omega,:) = diag(s_y_val);
end

s_y_pred_all(:,:,i_obs) = s_y_pred;

%----------------------------------------------------------------
% Clean-Up
%----------------------------------------------------------------

clear A_i A_pi A_y crdy crpi crr cry I_m m_seq_aux m_star Pi_m s_y_fn s_y_pred s_y_val t var_pred Y_m

end

clear i_obs i_omega i_y omega_val

%% PLOT RESULTS

%----------------------------------------------------------------
% Import Truth
%----------------------------------------------------------------

load s_y_true

%----------------------------------------------------------------
% Re-Scale
%----------------------------------------------------------------

% defining scaling coefficients

scale = NaN(3,1);

for i_y = 1:3
    scale(i_y) = max(s_y_base(1,i_y),s_y_cnfctl(1,i_y));
end

% do the re-scaling

for i_y = 1:3
    s_y_base(:,i_y)       = s_y_base(:,i_y) ./ scale(i_y);
    s_y_cnfctl(:,i_y)     = s_y_cnfctl(:,i_y) ./ scale(i_y);
    s_y_pred_all(:,i_y,:) = s_y_pred_all(:,i_y,:) ./ scale(i_y);
end

%----------------------------------------------------------------
% Figure Settings
%----------------------------------------------------------------

% color settings

settings.colors.black  = [0 0 0];
settings.colors.dgrey  = [130/255 130/255 130/255];
settings.colors.pink  = [255/255 192/255 203/255];
settings.colors.blue   = [116/255 158/255 178/255];
settings.colors.lblue  = 0.25 * settings.colors.blue + 0.75 * [1 1 1];
settings.colors.models = [30/255 129/255 176/255; ...
                            226/255 135/255 67/255; ...
                            135/255 62/255 35/255; ...
                            118/255 181/255 197/255];

% line specification

settings.line_specs = {'-.', '--', ':', '-.'};

% plot size

plotwidth = 0.27;
gapsize = 0.05;
gapsize_edges = (1-3*plotwidth-2*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth, gapsize_edges + 2 * gapsize + 2 * plotwidth];

% names

series_names = {'Interest Rate', 'Output', 'Inflation'};

%----------------------------------------------------------------
% Figure
%----------------------------------------------------------------

cd([path vintage task '/sw_sequence/_results']);

figure()

for i_y = 1:3

subplot(1,3,i_y)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(i_y);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
plot(omega_grid,s_y_base(:,i_y),'-','Color',settings.colors.pink,'LineWidth',5)
hold on
plot(omega_grid,s_y_cnfctl(:,i_y),'-','Color',settings.colors.black,'LineWidth',5)
hold on
for i_obs = 1:n_obs
    plot(omega_grid,s_y_pred_all(:,i_y,i_obs),settings.line_specs{i_obs},'Color',settings.colors.models(i_obs,:),'LineWidth',5)
    hold on
end
xlim([min(omega_grid) max(omega_grid)])
ylim([0 1])
yticks([0:0.2:1])
set(gcf,'color','w')
title(series_names(i_y),'interpreter','latex','fontsize',24)
xlabel('Frequency','interpreter','latex','FontSize',20)
if i_y == 3
    legend({'Old Rule','Cnfctl','Y,$\Pi$,R','+ L','+ I \& C','All'},'Location','Northeast','fontsize',18,'interpreter','latex','NumColumns',2,'Orientation','horizontal')
end
grid on
hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1
    print('cnfctl_2mom_pred','-dpng');
end

cd([path vintage task '/sw_sequence']);