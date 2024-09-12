function [cond_var,cond_var_1] = Kalman_filter(A,B,C,y)

% number of observations and variables
T = size(y,1);
n_s = size(A,1);

% initial conditions
st = zeros(n_s,1);
Pt = zeros(n_s,n_s);

% filtering
for t=1:T

    yt=y(t,:)';

    % prediction equations    
    st_1 = A * st;    
    Pt_1 = A * Pt * A' + B * B';

    % innovations 
    yt_1 = C*st_1;
    vt = yt-yt_1;

    % updating equations       
    Kt = Pt_1 * C' * (C * Pt_1 * C')^(-1);
    st = st_1 + Kt * vt;   
    Pt = Pt_1 - Kt * C * Pt_1;
    
end

cond_var   = Pt;
cond_var_1 = Pt_1;

end