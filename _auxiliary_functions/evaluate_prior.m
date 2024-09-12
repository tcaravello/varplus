function out = evaluate_prior(param,model)

if strcmp(model, 'rank')
    global prior_h
else
    global prior_theta
end

global prior_kappa prior_calvo_p prior_zeta_p prior_calvo_w prior_zeta_w prior_psi_uti prior_m_d prior_m_f
global prior_gamma_nkpc prior_omega_is prior_kappa_nkpc prior_eis prior_h_teq prior_i_gamma prior_frisch

if strcmp(model, 'teq')
if length(param)>4
    behav_prior =  log(prior_m_f(param(5)));
else
    behav_prior = 0;
end

elseif strcmp(model, 'teq_simple')
    behav_prior = 0;
else

if length(param)>7
    behav_prior =  log(prior_m_d(param(8)))+log(prior_m_f(param(9)));
else
    behav_prior = 0;
end
end

if strcmp(model, 'rank')
    error_fit =...
   - log(prior_h(param(1))) - log(prior_calvo_p(param(2))) - log(prior_zeta_p(param(3))) ...
   - log(prior_calvo_w(param(4))) - log(prior_zeta_w(param(5))) - log(prior_kappa(param(6))) ...
   - log(prior_psi_uti(param(7)))-behav_prior;
elseif strcmp(model, 'teq_simple')
error_fit =...
    - log(prior_gamma_nkpc(param(1))) - log(prior_omega_is(param(2))) - log(prior_kappa_nkpc(param(3))) ...
   - log(prior_eis(param(4)));
elseif strcmp(model, 'teq')
     error_fit =  ...
   - log(prior_h_teq(param(1))) - log(prior_calvo_p(param(2)))-log(prior_i_gamma(param(3)))...
   -log(prior_frisch(param(4)))-behav_prior;
else
   error_fit =...
   - log(prior_theta(param(1))) - log(prior_calvo_p(param(2))) - log(prior_zeta_p(param(3))) ...
   - log(prior_calvo_w(param(4))) - log(prior_zeta_w(param(5))) - log(prior_kappa(param(6))) ...
   - log(prior_psi_uti(param(7))) - behav_prior;
end

out = exp(-error_fit);

end