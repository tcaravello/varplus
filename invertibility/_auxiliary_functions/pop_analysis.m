function [IRF,M,tot_weights] = pop_analysis(model,settings);

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_s   = model.n_s;
n_eps = model.n_eps;
n_y   = model.n_y;
A     = model.ABCD.A;
B     = model.ABCD.B;
C     = model.ABCD.C;
D     = model.ABCD.D;

VMA_hor = settings.VMA_hor;

%----------------------------------------------------------------
% IRFs
%----------------------------------------------------------------

IRF = NaN(VMA_hor,n_y,n_eps);

for shock = 1:n_eps
    for i = 1:VMA_hor
        if i == 1
            resp = D';
            IRF(1,:,shock) = resp(shock,:);
        elseif i == 2
            resp = (C* B)';
            IRF(2,:,shock) = resp(shock,:);
        else
            resp = (C * A^(i-2) * B)';
            IRF(i,:,shock) = resp(shock,:);
        end
    end
end

%----------------------------------------------------------------
% Link Reduced-Form Errors and Structural Shocks
%----------------------------------------------------------------

% derive Sigma and K

Sigma_states = eye(n_s);
K_states = zeros(n_s,n_y);

dist = 1;
tol = 10^(-10);
relax = 0.9;

while dist >= tol
    Sigma_upd = (A - K_states*C)*Sigma_states*(A-K_states*C)' + B*B' + K_states*(D*D')*K_states' - B * D'*K_states' - K_states * D * B';
    K_upd = (A * Sigma_upd * C' + B * D') * (C * Sigma_upd * C' + D * D')^(-1);
    
    dist1 = max(max(abs(Sigma_upd - Sigma_states)));
    dist2 = max(max(abs(K_upd - K_states)));
    
    Sigma_states = relax * Sigma_states + (1-relax) * Sigma_upd;
    K_states = relax * K_states + (1-relax) * K_upd;
    
    dist = max(dist1, dist2);
end

% get the M matrices

aux = [A, zeros(n_s,n_s); K_states * C, A - K_states * C];
aux_lag = NaN(n_y,n_eps,VMA_hor);

for i = 1:VMA_hor
    aux_lag(:,:,i) = [C, -C] * aux^(i-1) * [B; K_states*D];
end

M = NaN(n_y,n_eps,VMA_hor+1);

M(:,:,1) = D;

for i = 2:VMA_hor+1
    M(:,:,i) = aux_lag(:,:,i-1); % u = M(L) * epsilon
end

%----------------------------------------------------------------
% Total Available Weights
%----------------------------------------------------------------

tot_weights = NaN(n_eps,VMA_hor+1);

tot_weights(:,1) = diag(D' * (C * Sigma_states * C' + D * D')^(-1) * D);
for i = 2:VMA_hor+1
    tot_weights(:,i) = diag(([C, -C] * aux^(i-2) * [B; K_states*D])' ...
        * (C * Sigma_states * C' + D * D')^(-1) * ([C, -C] * aux^(i-2) * [B; K_states*D]));
end

end