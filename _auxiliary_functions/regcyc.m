function [yreg,ypred] = regcyc(y)
% Matlab code to calculate cyclical component based on 2-year-ahead forecast error from linear regression as recommended in
%      James D. Hamilton, "Why You Should Never Use the Hodrick-Prescott Filter"
%      Review of Economics and Statistics, forthcoming
% input:  y = (T x 1) vector of data, tth element is observation for date t
% output  yreg = (T x 1) vector, tth element is cyclical component for date t

T = size(y,1);
yreg = NaN(T,1);
ypred = NaN(T,1);
h = 8;    % default for quarterly data and 2-year horizon
p =  4;   % default for quarterly data (number of lags in regression)

% construct X matrix of lags
X = ones(T-p-h+1,1);
for j = 1:p
   X = [X y(p-j+1:T-h-j+1,1)];
end

% do OLS regression
b = inv(X'*X)*X'*y(p+h:T,1);
yreg(p+h:T,1) = y(p+h:T,1) - X*b;
ypred(p+h:T,1) = X*b;
