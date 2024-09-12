%% GET OUTCOMES

% get marginal costs

mc_seq = ((cfc-1)*curvp+1)/((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp) * ...
    ((1+cbetabar*cgamma*cindp) * pi_seq - cbetabar*cgamma*[pi_seq(2:T);0] - cindp * [0;pi_seq(1:T-1)]);

% get wages

w_seq = 1/(1-calfa) * (mc_seq - calfa * rk_seq);

% get capacity utilization

zcap_seq = (1/(czcap/(1-czcap)))* rk_seq;

% get effective capital

k_seq = [0;kp_seq(1:T-1)] + zcap_seq;

% get labor

lab_seq = rk_seq - w_seq + k_seq;

% get investment

inve_seq = 1/cikbar * (kp_seq - (1-cikbar) * [0;kp_seq(1:T-1)]);

% get the price of capital

pk_seq = (cgamma^2*csadjcost) * ((1+cbetabar*cgamma) * inve_seq - [0;inve_seq(1:T-1)] - cbetabar*cgamma*[inve_seq(2:T);0]);

% get nominal rates

r_seq = -pk_seq + [pi_seq(2:T);0] + (crk/(crk+(1-ctou)))*[rk_seq(2:T);0] + ((1-ctou)/(crk+(1-ctou)))*[pk_seq(2:T);0];

% get output

y_seq = cfc*(calfa*k_seq+(1-calfa)*lab_seq);

% get consumption

c_seq = 1/ccy * (y_seq - ciy*inve_seq - crkky * zcap_seq);

% get monetary policy wedge

ms_seq = zeros(T,1);
ms_seq(1) = 0.2290 * m_seq(1);
for t = 2:T
    ms_seq(t) = 0.1999 * ms_seq(t-1) + 0.2290 * m_seq(t);
end