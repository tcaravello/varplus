%% Computes  q, the eigenvector associated to the largest eigenvalue in the matrix Theta from equation 1 in Angeletos. Translated from his coe
function[eigvec, eigval] =  funcfdq(MX, MY,S,wmin,wmax,idx,weight,grid)
    % Computes the impulse vector for a specific set of variables
    freqs   = linspace(0.0,2*pi,grid);
    F       = zeros(size(freqs));
    F( ( freqs >= wmin )& (freqs <= wmax)) = 1;
    F( ( freqs >= (2*pi - wmax))& (freqs <= (2*pi - wmin))) = 1;
    ne      = size(S,2);
    nx      = size(MX,1);
    ME      = [S(:,1:ne);zeros(nx-ne,ne)];
    zi      = exp(-(1i)*freqs);
    r2pi    = 1/(2*pi);
    nidx    = length(idx);
    VD      = zeros(ne,ne);
    for ii=1:nidx
        sp      = complex(zeros(grid,1));
        sp2     = complex(zeros(grid,ne*ne));

        for gp=1:grid
            if F(gp)==1
                fom     = MY(idx(ii),:)*((eye(nx)-MX*zi(gp))\ME);
                tmp     = r2pi*(fom*fom');
                tmp     = F(gp)*tmp;
                sp(gp)  = tmp;
                tmp     = r2pi*(fom'*fom);
                tmp     = F(gp)*tmp;
                sp2(gp,:) = tmp(:)';
            end
        end
        sp(isnan(sp))   = 0.0;
        sp2(isnan(sp2)) = 0.0;
        VTtmp           = 2*pi*real(ifft(sp));
        VDtmp           = 2*pi*real(ifft(sp2));
        VD              = VD+weight(ii)*reshape(VDtmp(1,:)/VTtmp(1),ne,ne);
    end
    [P,D]   = eig(VD);
    [val, ii]     = max(max(abs(D)));
    if size(find(D == val), 1) > 1
        print("Multiple identical (maximum) eigenvalues found")
        eigvec = P;
        eigval = D;
    else
        eigvec =  real(P(:,ii));
        eigval = val;
    end
end