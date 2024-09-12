%% FORECAST UNCERTAINTY IN SMETS-WOUTERS
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
settings.IRF_hor          = 30; % horizon for computed IRFs
settings.IS_draws         = 200; % number of draws from the identified set
settings.weight_hor       = settings.VMA_hor; % horizon for weight polynomial P(L)
settings.T_fcst           = 200; % maximal horizon for forecast uncertainty

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

fcst_var = zeros(settings.T_fcst+1,3,n_obs);

%% COMPUTE FORECAST UNCERTAINTY

for i_obs = 1:n_obs

%----------------------------------------------------------------
% Solve SW Model
%----------------------------------------------------------------

% model run

dynare SW_Model noclearall nolog
pause(5)
clean_folder_SW

disp('I have solved and simulated the model.')

disp('Collecting model properties...')

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

disp('...done!')

%----------------------------------------------------------------
% Get VAR Representation
%----------------------------------------------------------------

disp('Getting the VAR representation...')

VAR = popVAR(SW_model,settings);

disp('...done!')

%----------------------------------------------------------------
% Get Uncertainty
%----------------------------------------------------------------

disp('Getting the forecast uncertainty...')

VAR_coeff = VAR.VAR_coeff;
Sigma_u   = VAR.Sigma_u;
n_y       = SW_model.n_y;
p         = VAR.laglength;

A_KF = zeros(n_y*p,n_y*p);
for i = 1:p
    A_KF(1:n_y,1+(i-1)*n_y:i*n_y) = VAR_coeff(1+(i-1)*n_y:i*n_y,:)';
end
A_KF(n_y+1:n_y*p,1:n_y*(p-1)) = eye(n_y*(p-1));

B_KF = zeros(n_y*p,n_y);
B_KF(1:n_y,:) = Sigma_u^(0.5);

T_KF = settings.T_fcst;

var_list = zeros(n_y*p,n_y*p,T_KF);
var_list(:,:,1) = B_KF * B_KF';
for t = 2:T_KF
    var_list(:,:,t) = A_KF * var_list(:,:,t-1) * A_KF' + B_KF * B_KF';
end

var_list = var_list(1:3,1:3,:);

for i_y = 1:3
    fcst_var(2:end,i_y,i_obs) = squeeze(var_list(i_y,i_y,:));
end

clear VAR_coeff Sigma_u n_y p VMA_hor A_KF B_KF T_KF

disp('...done!')

end

%% PLOT RESULTS

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

% horizon

T_plot = 100;

%----------------------------------------------------------------
% Re-Scale
%----------------------------------------------------------------

fcst_var_scale = NaN(settings.T_fcst+1,3,n_obs);

for i_y = 1:3
    for i_obs = 1:n_obs
        fcst_var_scale(:,i_y,i_obs) = fcst_var(:,i_y,i_obs) ./ fcst_var(:,i_y,4);
    end
end

fcst_var_scale(1,:,:) = 1;

%----------------------------------------------------------------
% Figure
%----------------------------------------------------------------

cd([path vintage task '/sw_sequence/_results']);

figure

for i_y = 1:3

subplot(1,3,i_y)
pos = get(gca, 'Position');
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex')
pos(1) = left_pos(i_y);
pos(3) = plotwidth;
set(gca,'Position', pos)
hold on
for i_obs = 1:n_obs
    plot(0,0,settings.line_specs{i_obs},'Color',settings.colors.models(i_obs,:),'LineWidth',5)
    hold on
end
hold on
for i_obs = 1:n_obs
    plot(1:1:T_plot,fcst_var_scale(2:T_plot+1,i_y,i_obs),settings.line_specs{i_obs},'Color',settings.colors.models(i_obs,:),'LineWidth',5)
    hold on
end
xlim([1 T_plot])
ylim([1 1.2])
yticks([1:0.05:1.2])
set(gcf,'color','w')
title(series_names(i_y),'interpreter','latex','fontsize',24)
xlabel('Horizon','interpreter','latex','FontSize',20)
% if i_y == 1
%     ylabel('Rel. Forecast Variance','interpreter','latex','FontSize',20)
% end
if i_y == 3
    legend({'Y,$\Pi$,R','+ L','+ I \& C','All'},'Location','Northeast','fontsize',18,'interpreter','latex','NumColumns',2,'Orientation','horizontal')
end
grid on
hold off

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 2.25*pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
if save_fig == 1;
    print('cnfctl_2mom_fcst_vars','-dpng');
end
cd([path vintage task '/sw_sequence']);