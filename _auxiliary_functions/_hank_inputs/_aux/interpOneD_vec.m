function V = interpOneD_vec(x_fine,x_knot)
% Compute Linear interpolation over 1d grid. It extrapolates outside of
% knot points
%
% by SeHyoun Ahn, March 2017
%
% PARAMETERS:
%    x_fine = x grid to interpolate to
%    x_knot = x knot points to interpolate from
%
% OUTPUT:
%    V = matrix giving interpolation values
%
% EXAMPLE:
%    x_fine = linspace(-1,1,300)';
%    x_knot = linspace(0,2,20)';
%    V = interpOneD(x_fine,x_knot);
%    plot(V*(x_knot.^5+x_knot*0.5));
%

n_x_fine = size(x_fine,1);
n_x_knot = size(x_knot,1);
n_epsi = size(x_fine,2);
loc = sum((repmat(x_fine,1,1,n_x_knot) - permute(repmat(x_knot,1,1,n_x_fine),[3 2 1]))>=0,3);
% loc = min(loc,length(x_knot)-1);
loc = min(loc,size(x_knot,1)-1);
loc = max(loc,1);

x_fine_vec = x_fine(:);
x_knot_vec = x_knot(:);
offset = repmat([0:1:n_epsi-1] * n_x_knot,n_x_fine,1);
loc_temp = loc + offset;
loc_vec = loc_temp(:);
t = (x_fine_vec-x_knot_vec(loc_vec))./(x_knot_vec(loc_vec+1)-x_knot(loc_vec));
ind_x = 1:length(x_fine_vec);

i = repmat(ind_x',2,1);
j = [loc_vec;loc_vec+1];
v = [(1-t);t];

V = sparse(i,j,v,length(x_fine_vec),length(x_knot_vec));