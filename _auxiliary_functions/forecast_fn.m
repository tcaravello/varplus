function var_forecasts = forecast_fn(data,n_lags,constant,B,fcst_date,fcst_length,hist_indic);

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_var = size(data,2); 

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
u = Y-X*B;

% augment X to contain all data, write in companion form.

%xTT = [Y(end,:),Y(end-1,:)];

xTT = zeros(1,n_var*n_lags);
for ll = 1:n_lags
    xTT(1,(ll-1)*n_var+1:ll*n_var) = Y(end+1-ll,:); 
end

if constant==1
    xTT(:,n_var*n_lags+1)=1;
elseif constant == 2
    xTT(:,n_var*n_lags+1)=1;
    xTT(:,n_var*n_lags+2)=T+1;
elseif constant == 3
    xt(:,n_var*n_lags+1)=1;
    xt(:,n_var*n_lags+2)=T+1;
    xt(:,n_var*n_lags+3)=(T+1)^2;
end

X = [X;xTT];

%----------------------------------------------------------------
% Construct Forecasts
%----------------------------------------------------------------

fcst_date = fcst_date - n_lags + 1;

% Determininstic terms.

X_t = NaN(fcst_date+fcst_length,size(X,2));
if constant==1
    X_t(:,n_var*n_lags+1)=ones(fcst_date+fcst_length,1);
elseif constant == 2
    X_t(:,n_var*n_lags+1)=ones(fcst_date+fcst_length,1);
    X_t(:,n_var*n_lags+2)=[1:1:fcst_date+fcst_length]';
elseif constant == 3
    X_t(:,n_var*n_lags+1)=ones(fcst_date+fcst_length,1);
    X_t(:,n_var*n_lags+2)=[1:1:fcst_date+fcst_length]';
    X_t(:,n_var*n_lags+3)=([1:1:fcst_date+fcst_length].^2)';
end
X_t(1:fcst_date,:) = X(1:fcst_date,:);

% construct forecast recursively: Carry forecasts for y_{t},y_{t-1}, etc,
for i = 1:fcst_length
    X_t(fcst_date+i,1:n_var) = X_t(fcst_date+i-1,:) * B; %set $E_{t}[y_{t+h}] = B E_{t}[y_{t+h-1}] using the VAR equation.
    for j = 2:n_lags % and then for the lags, just use identities 
        X_t(fcst_date+i,n_var*(j-1)+1:n_var*j) = X_t(fcst_date+i-1,1+n_var*(j-2):n_var*(j-1));
    end
end

% whether you want to include the current period or not.
if hist_indic == 1
    var_forecasts = X_t(fcst_date+1:fcst_date+fcst_length,1:n_var);
elseif hist_indic == 2
    var_forecasts = X_t(fcst_date:fcst_date+fcst_length-1,1:n_var);
end