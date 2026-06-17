function plot_pmesh_multi(N, meshVel, meshWall, u)
[z,~]=zwgl(N-1);
NelxV = length(meshVel)-1;
NelyV = length(meshVel)-1;
nxV = NelxV*(N-2) + 1;
nyV = NelyV*(N-2) + 1;
U = zeros(nyV,nxV);
X = zeros(nyV,nxV);
Y = zeros(nyV,nxV);
n = N-1;
n2 = n^2;
figure; hold on;
for ey = 1:NelyV
for ex = 1:NelxV
    e = (ey-1)*NelxV + ex;
    idx = (e-1)*n2 + (1:n2);
    uel = reshape(u(idx),n,n);
    x0 = meshVel(ex);
    x1 = meshVel(ex+1);
    y0 = meshVel(ey);
    y1 = meshVel(ey+1);
    xloc = (x1-x0)/2*z + (x0+x1)/2;
    yloc = (y1-y0)/2*z + (y0+y1)/2;
    [YY,XX] = meshgrid(yloc,xloc);
    surf(XX,YY,uel);
end
end

% shading interp;
view(2);
colorbar;
