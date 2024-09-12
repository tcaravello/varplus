function [Beta,Res] = ls_detrend(Y,const)

T = size(Y,1);
if const == 1
    X = ones(T,1);
elseif const == 2
    X = [ones(T,1),(1:1:T)'];
end

Beta = X \ Y;
Res = Y - X * Beta;

end