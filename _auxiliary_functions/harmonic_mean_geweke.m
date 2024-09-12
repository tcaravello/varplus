function [marginal_lik, quantiles_use] = harmonic_mean_geweke(param_collector, log_posterior_collector, settings_HM, model)

global prior_zeta_p prior_zeta_w

if size(param_collector,1) == 7


index_use = settings_HM.index_use;
tau = settings_HM.tau ;
target_Sigma_inv = settings_HM.target_Sigma_inv;
N_IRFS = settings_HM.N_IRFS;

param_p = param_collector(3,1);
param_w = param_collector(5,1);

log_posterior_use = log_posterior_collector + (N_IRFS/2)*log(2*pi) - 0.5*log(det(target_Sigma_inv)) + log(prior_zeta_p(param_p)) + log(prior_zeta_w(param_w));

n_param = length(index_use);
N_use = length(log_posterior_collector);
param_use = param_collector(index_use,:);
chol_cov = chol(cov(param_use'))\eye(n_param);
f_exp = chol_cov*(param_use - mean(param_use,2));
quantiles_use = quantile(f_exp',[0.025 0.16 0.84 0.975]); % I compute this as a diagnostic for non-normality.

f_exp_sq = f_exp.^2;
threshold = chi2inv(tau,n_param);
indicator_use = (sum(f_exp_sq,1)<=threshold);

%build f vector

f_vec = (1/tau)*((2*pi)^(-n_param/2))*(det(chol_cov))*exp(-0.5*sum(f_exp_sq,1)).*indicator_use;

%posterior_use_inv = exp(-log_posterior_use);
%inside_sum = f_vec.*exp(log_posterior_use);

% compute harmonic mean. Notice that log_posterior_use already has a minus
% in front.

marginal_lik = 1/mean(f_vec.*exp(log_posterior_use));

elseif size(param_collector,1) == 9

index_use = settings_HM.index_use;
tau = settings_HM.tau ;
target_Sigma_inv = settings_HM.target_Sigma_inv;
N_IRFS = settings_HM.N_IRFS;

param_p = param_collector(3,1);
param_w = param_collector(5,1);

if length(index_use) == 5

param_m_d = param_collector(8,1);
param_m_f = param_collector(9,1);

% create the beta functions used for the prior of the behavioral parameters

%Demand side (consumption, investment) behavioral parameter

m_d_mean = param_m_d; %this function encodes sort the uncertainty we have.. 
m_d_sd = 0.00001;

%m_d_mean = 0.999;
%m_d_sd = 0.001;

m_d_a = m_d_mean*(m_d_mean*(1-m_d_mean)/(m_d_sd^2)-1);
m_d_b = (1-m_d_mean)*(m_d_mean*(1-m_d_mean)/(m_d_sd^2)-1);

prior_m_d = @(x) betapdf(x,m_d_a,m_d_b);

%Pricing side (NKPC, wage PC) behavioral parameter

m_f_mean = param_m_f;
m_f_sd = 0.00001;


m_f_a = m_f_mean*(m_f_mean*(1-m_f_mean)/(m_f_sd^2)-1);
m_f_b = (1-m_f_mean)*(m_f_mean*(1-m_f_mean)/(m_f_sd^2)-1);

prior_m_f = @(x) betapdf(x,m_f_a,m_f_b);



log_posterior_use = log_posterior_collector + (N_IRFS/2)*log(2*pi) - 0.5*log(det(target_Sigma_inv)) ...
    + log(prior_zeta_p(param_p)) + log(prior_zeta_w(param_w)) + log(prior_m_d(param_m_d)) + log(prior_m_f(param_m_f));

   elseif length(index_use) == 7

log_posterior_use = log_posterior_collector + (N_IRFS/2)*log(2*pi) - 0.5*log(det(target_Sigma_inv)) ...
    + log(prior_zeta_p(param_p)) + log(prior_zeta_w(param_w));
    end

n_param = length(index_use);
N_use = length(log_posterior_collector);
param_use = param_collector(index_use,:);
chol_cov = chol(cov(param_use'))\eye(n_param);
f_exp = chol_cov*(param_use - mean(param_use,2));
quantiles_use = quantile(f_exp',[0.025 0.16 0.84 0.975]); % I compute this as a diagnostic for non-normality.

f_exp_sq = f_exp.^2;
threshold = chi2inv(tau,n_param);
indicator_use = (sum(f_exp_sq,1)<=threshold);

%build f vector

f_vec = (1/tau)*((2*pi)^(-n_param/2))*(det(chol_cov))*exp(-0.5*sum(f_exp_sq,1)).*indicator_use;

posterior_use_inv = exp(-log_posterior_use);
inside_sum = f_vec.*exp(log_posterior_use);

% compute harmonic mean. Notice that log_posterior_use already has a minus
% in front.

marginal_lik = 1/mean(f_vec.*exp(log_posterior_use));


end
end