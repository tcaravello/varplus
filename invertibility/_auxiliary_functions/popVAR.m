function VAR_pop = popVAR(model,settings)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

n_s = model.n_s;
n_y = model.n_y;
A   = model.ABCD.A;
B   = model.ABCD.B;
C   = model.ABCD.C;
D   = model.ABCD.D;

VAR_pop.laglength = settings.VAR_poplaglength;
VAR_laglength     = VAR_pop.laglength;
VMA_hor           = settings.VMA_hor;

%----------------------------------------------------------------
% VAR in Macro Aggregates
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

% get VAR

VAR_pop.VAR_coeff = NaN(VAR_laglength*n_y+1,n_y);

for i = 1:VAR_laglength
    VAR_pop.VAR_coeff(1+(i-1)*n_y:i*n_y,:) = (C * (A - K_states * C)^(i-1) * K_states)';
end

VAR_pop.Sigma_u = C * Sigma_states * C' + D * D';

%----------------------------------------------------------------
% Wold IRFs
%----------------------------------------------------------------

VAR_pop.IRF_Wold = zeros(n_y,n_y,VMA_hor);
VAR_pop.IRF_Wold(:,:,1) = eye(n_y);

for l = 1:VMA_hor
    
    if l<VMA_hor
        for j=1:min(l,VAR_laglength)
            VAR_pop.IRF_Wold(:,:,l+1) = VAR_pop.IRF_Wold(:,:,l+1) + VAR_pop.VAR_coeff(1+(j-1)*n_y:j*n_y,:)'*VAR_pop.IRF_Wold(:,:,l-j+1);
        end
    end
    
end

end