%clear all;
function [H] = make_jinc(P, scale)

x = [-floor(P/2) : 1 : floor(P/2)];
y = x;

[xx,yy]=meshgrid(x,y);

r = sqrt(xx.^2+yy.^2);

z = ones(size(r));
k = find(r);
z(k) = besselj(1, 2*pi*r(k)*scale) ./ ( 2*pi*r(k)*scale );
z(ceil(P/2), ceil(P/2)) = 0.5;

H = z / sum(sum(z));
