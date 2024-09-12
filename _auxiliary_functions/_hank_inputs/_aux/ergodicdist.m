function q=ergodicdist(Q,alg)
% ErgodicDist - Ergodic distribution of discrete Markov Chains
%
% q=ergodicdist(Q)
%
% Input: Q -   right stochastic matrix (COLUMNS sum to one): the prob.
%              of reaching state j from state i is P(i,j).
%        alg - choice of method: 1: iterative, 2: direct. (optional;
%              default = 1)
%
% Output: q - (nx1) vector containing the ergodic distribution
%
% Author: Marco Maffezzoli. Ver. 1.0.1, 11/2012.
%

if nargin==1
    alg=1;
end
h=size(Q,1);
switch alg
    case 1
        q=sparse(1,h);
        q(1,1)=1;
        %q = q+1/h ; 
        dif=1; 
        while dif>1e-8
            z=q*Q;
            dif=norm(z-q);
            q=z;
        end
        q=q';
    case 2
        [q,~] = eigs(Q',1);
        q = q'./sum(q);
        dif = 1;
        while dif>1e-8
            z=q*Q;
            dif=norm(z-q);
            q=z;
        end
        q=q';
    otherwise
        error('Please set alg = 1 or 2!')
end

end