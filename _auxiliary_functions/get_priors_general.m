% habit parameter in RANK
global h_a h_b prior_h

h_mean = 0.7;
h_sd = 0.1;

% we always need var< mu(1-mu)
h_a = h_mean*(h_mean*(1-h_mean)/(h_sd^2)-1);
h_b = (1-h_mean)*(h_mean*(1-h_mean)/(h_sd^2)-1);

prior_h = @(x) betapdf(x,h_a,h_b);

global theta_a theta_b prior_theta

% sticky info parameter in HANK
theta_mean = 0.7;
theta_sd = 0.2;

% we always need var< mu(1-mu)
theta_a = theta_mean*(theta_mean*(1-theta_mean)/(theta_sd^2)-1);
theta_b = (1-theta_mean)*(theta_mean*(1-theta_mean)/(theta_sd^2)-1);

prior_theta = @(x) betapdf(x,theta_a,theta_b);

global kappa_mean kappa_sd prior_kappa calvo_p_a  calvo_p_b  prior_calvo_p...
        zeta_p_a zeta_p_b prior_zeta_p calvo_w_a  calvo_w_b  prior_calvo_w zeta_w_a zeta_w_b prior_zeta_w ...
        psi_uti_a psi_uti_b prior_psi_uti m_d_a m_d_b prior_m_d m_f_a m_f_b prior_m_f

% Calvo price parameter

calvo_p_mean = 2/3;
calvo_p_sd = 0.2;

calvo_p_a = calvo_p_mean*(calvo_p_mean*(1-calvo_p_mean)/(calvo_p_sd^2)-1);
calvo_p_b = (1-calvo_p_mean)*(calvo_p_mean*(1-calvo_p_mean)/(calvo_p_sd^2)-1);

prior_calvo_p = @(x) betapdf(x,calvo_p_a,calvo_p_b);

% price inertia (numerical way to essentially fix it at 1)
zeta_p_mean = 0.999;
zeta_p_sd = 0.001;

% we always need var< mu(1-mu)
zeta_p_a = zeta_p_mean*(zeta_p_mean*(1-zeta_p_mean)/(zeta_p_sd^2)-1);
zeta_p_b = (1-zeta_p_mean)*(zeta_p_mean*(1-zeta_p_mean)/(zeta_p_sd^2)-1);

prior_zeta_p = @(x) betapdf(x,zeta_p_a,zeta_p_b);

% Calvo wage parameter
calvo_w_mean = 2/3;
calvo_w_sd = 0.2;

calvo_w_a = calvo_w_mean*(calvo_w_mean*(1-calvo_w_mean)/(calvo_w_sd^2)-1);
calvo_w_b = (1-calvo_w_mean)*(calvo_w_mean*(1-calvo_w_mean)/(calvo_w_sd^2)-1);

prior_calvo_w = @(x) betapdf(x,calvo_w_a,calvo_w_b);

% wage inertia (numerical way to essentially fix it at 1)
zeta_w_mean = 0.999;
zeta_w_sd = 0.001;

% we always need var< mu(1-mu)
zeta_w_a = zeta_w_mean*(zeta_w_mean*(1-zeta_w_mean)/(zeta_w_sd^2)-1);
zeta_w_b = (1-zeta_w_mean)*(zeta_w_mean*(1-zeta_w_mean)/(zeta_w_sd^2)-1);

prior_zeta_w = @(x) betapdf(x,zeta_w_a,zeta_w_b);

% investment adjustment cost.
kappa_mean = 5;
kappa_sd = 1.5;
prior_kappa = @(x) normpdf(x,kappa_mean,kappa_sd);

% capacity utilization (this is paramterized by psi_uti in SW, with zeta = psi/(1-psi))
psi_uti_mean = 0.5;
psi_uti_sd = 0.15;

% we always need var< mu(1-mu)
psi_uti_a = psi_uti_mean*(psi_uti_mean*(1-psi_uti_mean)/(psi_uti_sd^2)-1);
psi_uti_b = (1-psi_uti_mean)*(psi_uti_mean*(1-psi_uti_mean)/(psi_uti_sd^2)-1);
prior_psi_uti = @(x) betapdf(x,psi_uti_a,psi_uti_b);

% demand side (consumption, investment) behavioral parameter

m_d_mean = 0.999999; %set it to 1. 
m_d_sd = 0.00001;

m_d_a = m_d_mean*(m_d_mean*(1-m_d_mean)/(m_d_sd^2)-1);
m_d_b = (1-m_d_mean)*(m_d_mean*(1-m_d_mean)/(m_d_sd^2)-1);

prior_m_d = @(x) betapdf(x,m_d_a,m_d_b);

% pricing side (NKPC, wage PC) behavioral parameter

m_f_mean = 0.65;
m_f_sd = 0.00001;

m_f_a = m_f_mean*(m_f_mean*(1-m_f_mean)/(m_f_sd^2)-1);
m_f_b = (1-m_f_mean)*(m_f_mean*(1-m_f_mean)/(m_f_sd^2)-1);
prior_m_f = @(x) betapdf(x,m_f_a,m_f_b);