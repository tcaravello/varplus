% Computes IRFs
function[YY, XX] = comp_irf(ME, MX, MY, nrep, kk)
[nx, ne] = size(ME);
eyes = eye(ne);
XX = zeros(nx, nrep);
XX(:,1) = ME*eyes(:,kk);
for tt=2:nrep
    XX(:,tt) = MX * XX(:, tt - 1);
end

YY = (MY*XX)';
XX = XX';

end


