function y_seq = y_fn(guess_seq,shock_seq,param,T);

global ctou clandaw cg curvp curvw

%% COLLECT INPUTS

% parameters

csadjcost  = param(1);
csigma     = param(2);
chabb      = param(3);
cprobw     = param(4);
csigl      = param(5);
cprobp     = param(6);
cindw      = param(7);
cindp      = param(8);
czcap      = param(9);
cfc        = param(10);
crpi       = param(11);
crr        = param(12);
cry        = param(13);
crdy       = param(14);
constepinf = param(15);
constebeta = param(16);
constelab  = param(17);
ctrend     = param(18);
cgy        = param(19);
calfa      = param(20);

% several further auxiliary parameters

cpie	= 1+constepinf/100;
cgamma	= 1+ctrend/100 ;
cbeta	= 1/(1+constebeta/100);

clandap	 = cfc;
cbetabar = cbeta*cgamma^(-csigma);
cr		 = cpie/(cbeta*cgamma^(-csigma));
crk		 = (cbeta^(-1))*(cgamma^csigma) - (1-ctou);
cw 		 = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
cikbar	 = (1-(1-ctou)/cgamma);
cik		 = (1-(1-ctou)/cgamma)*cgamma;
clk		 = ((1-calfa)/calfa)*(crk/cw);
cky		 = cfc*(clk)^(calfa-1);
ciy		 = cik*cky;
ccy		 = 1-cg-cik*cky;
crkky	 = crk*cky;
cwhlc	 = (1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
cwly	 = 1-crk*cky;
conster	 = (cr-1)*100;

% guessed sequences

rk_seq  = guess_seq(1:T,1);
kp_seq  = guess_seq(T+1:2*T,1);
pi_seq  = guess_seq(2*T+1:3*T,1);

m_seq   = shock_seq;

%% GET OUTCOMES

get_aggregates

% y_seq = mc_seq;