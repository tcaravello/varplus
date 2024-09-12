function y = winsorize(x,varargin)
% Winsorising extreme values. 
% 
% Syntax: 	Y=WINSORISING(X,W)
%           X - data matrix or vector
%           For vectors, WINSORISING(X) is the winsorized X array. 
%           For matrices, WINSORISING(X) is a matrix containing the 
%           winsorized element from each column. 
% 
%           W - Amount of winsoritazion (90 by default). If you set W=90 this 
%           means that the remaing 10% (0-5th percentile and 95-100th
%           percentile) will be substituted.
% 
%      Example: 
% 
% x=[92 19 101 58 103 91 26 78 10 13 0 101 86 85 15 89 89 25 2 41];
% 
%           Calling on Matlab the function: winsorizing(x)
% 
%           Answer is:
% 
% y=92 19 101 58 102 91 26 78 10 13 1 101 86 85 15 89 89 25 2 41
% 
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
% 
% To cite this file, this would be an appropriate format:
% Cardillo G. (2011). WINSORISING: WINSORISING Data
% http://www.mathworks.com/matlabcentral/fileexchange/32327

%Input error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan'}));
addOptional(p,'w',90, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',100}));
parse(p,x,varargin{:});
x=p.Results.x; w=p.Results.w;
clear p

rw=(100-w)/2;
b=[rw 100-rw];
lb=prctile(x,min(b)); ub=prctile(x,max(b)); %set lower and upper bound

if isvector(x)
    y=x; y(y<lb)=lb; y(y>ub)=ub; %winsorising using logical indexing
else
    y=zeros(size(x));
    for I=1:size(x,2)
        c=x(:,I);
        c(c<lb(I))=lb(I); c(c>ub(I))=ub(I);
        y(:,I)=c;
    end
end