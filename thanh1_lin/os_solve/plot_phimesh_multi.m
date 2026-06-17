function plot_phimesh_multi(N, meshVel, meshWall, u)
meshPhi  = [-meshWall(end:-1:1), meshVel(2:end-1), meshWall];
[z,~]=zwgll(N);
Nelx = length(meshPhi)-1;
Nely = length(meshPhi)-1;
nx = Nelx*N + 1;
ny = Nely*N + 1;
U = zeros(ny,nx);
X = zeros(ny,nx);
Y = zeros(ny,nx);
n = N+1;
n2 = n^2;
figure; hold on;
for ey = 1:Nely
for ex = 1:Nelx
    e = (ey-1)*Nelx + ex;
    idx = (e-1)*n2 + (1:n2);
    uel = reshape(u(idx),n,n);
    x0 = meshPhi(ex);
    x1 = meshPhi(ex+1);
    y0 = meshPhi(ey);
    y1 = meshPhi(ey+1);
    xloc = (x1-x0)/2*z + (x0+x1)/2;
    yloc = (y1-y0)/2*z + (y0+y1)/2;
    [YY,XX] = meshgrid(yloc,xloc);
    surf(XX,YY,uel);
end
end

% shading interp;
view(2);
colorbar;
