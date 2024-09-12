function [B_draws,Sigma_draws,B_OLS,Sigma_OLS] = bvar_fn(data,n_lags,constant,n_draws);

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_var = size(data,2);                % number of endogenous variables
m     = n_var*n_lags + constant;     % number of exogenous variables

yt = data(n_lags+1:end,:);
T  = size(yt,1);
xt = zeros(T,n_var*n_lags+constant);
for i=1:n_lags
    xt(:,n_var*(i-1)+1:n_var*i) = data((n_lags-(i-1)):end-i,:) ;
end
if constant==1
    xt(:,n_var*n_lags+1)=ones(T,1);
elseif constant == 2
    xt(:,n_var*n_lags+1)=ones(T,1);
    xt(:,n_var*n_lags+2)=[1:1:T]';
elseif constant == 3
    xt(:,n_var*n_lags+1)=ones(T,1);
    xt(:,n_var*n_lags+2)=[1:1:T]';
    xt(:,n_var*n_lags+3)=([1:1:T].^2)';
end

% write data in Zellner (1971, pp 224-227) notation
Y = yt; % T by nvar matrix of observations
X = xt; % T by (nvar*nlag+1) matrix of regressors

%----------------------------------------------------------------
% Prior for Reduced Form Parameters
%----------------------------------------------------------------

nnuBar              = 0;
OomegaBarInverse    = zeros(m);
PpsiBar             = zeros(m,n_var);
PphiBar             = zeros(n_var);

%----------------------------------------------------------------
% Posterior for Reduced-Form Parameters
%----------------------------------------------------------------

nnuTilde            = T +nnuBar;
OomegaTilde         = (X'*X  + OomegaBarInverse)\eye(m);
OomegaTildeInverse  =  X'*X  + OomegaBarInverse;
PpsiTilde           = OomegaTilde*(X'*Y + OomegaBarInverse*PpsiBar);
PphiTilde           = Y'*Y + PphiBar + PpsiBar'*OomegaBarInverse*PpsiBar - PpsiTilde'*OomegaTildeInverse*PpsiTilde;
PphiTilde           = (PphiTilde'+PphiTilde)*0.5;

%----------------------------------------------------------------
% Draw from Posterior
%----------------------------------------------------------------

% definitios used to store orthogonal-reduced-form draws, volume elements, and unnormalized weights
B_draws         = NaN(m,n_var,n_draws); % reduced-form lag parameters
Sigma_draws     = NaN(n_var,n_var,n_draws); % reduced-form covariance matrices

% definition to facilitate the draws from B|Sigma
cholOomegaTilde = chol(OomegaTilde)'; % this matrix is used to draw B|Sigma below

% draws from posterior
for i_draw = 1:n_draws
    
    Sigmadraw     = iwishrnd(PphiTilde,nnuTilde);
    cholSigmadraw = chol(Sigmadraw)';
    Bdraw         = kron(cholSigmadraw,cholOomegaTilde)*randn(m*n_var,1) + reshape(PpsiTilde,n_var*m,1);
    Bdraw         = reshape(Bdraw,n_var*n_lags+constant,n_var);
    
    % store reduced-form draws
    B_draws(:,:,i_draw)     = Bdraw;
    Sigma_draws(:,:,i_draw) = Sigmadraw;
    
end

%----------------------------------------------------------------
% Collect Results
%----------------------------------------------------------------

B_OLS = PpsiTilde;
Sigma_OLS = PphiTilde/T;