% Flow in z-dir
% Magnetic field in x-dir

format long;

N = 10;

% Re = 3e4;
% Ha = 100;
% c = 1.0;
% W = 0.1;


Re = 1e4;
Ha = 15;
c = 1.0;
W = 0.1;

alpha = 1.0; % Wave number in z


rho = 1;
mu = 1/Re;
sigma = Ha^2/Re;
sigma_w = c*sigma/W;
By = 1;

meshVel  = [-1.0, -0.5, 0.0, 0.5, 1.0];
meshWall = [0.0, 1.0]*W + 1.0;



[U_base,Phi_base,fz] = solve_steady_wall_multi(N, meshVel, meshWall, mu, sigma, sigma_w, W, By);
[uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin_wall_multi(N, meshVel, meshWall, rho, mu, sigma,sigma_w, W, U_base, Phi_base, By, alpha);

plot_vmesh_multi(N, meshVel, meshWall, U_base)
plot_vmesh_multi(N, meshVel, meshWall, real(uvec))

rho
mu
sigma
sigma_w
fz
gamma_r = real(gamma)
gamma_i = imag(gamma)

% save("ha10n10.mat")
output_nek


