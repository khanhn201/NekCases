A = csvread('one_period.csv');

t = A(:,1);
u = A(:,2)/100.0; % cm to m


D = 8e-3; % D=8mm
Qvol = 0.25/1000/60; % L/min -> m^3/s
Area = pi*D^2/4;

[t, idx] = sort(t);
u = u(idx);
[t, ia] = unique(t);
u = u(ia);
% u = u(1:end/3);
% t = t(1:end/3);
t = -t/min(t) + 1;
plot(t,u);
figure;

N = 1024;
tu = linspace(min(t), max(t), N);
uu = interp1(t, u, tu, 'spline');

% Subtract mean flow rate
T = tu(end) - tu(1);
umean = trapz(tu, uu)/T;
uu = uu/umean*Qvol/Area;

% uu = uu + Qvol/Area;

% ----------------------------------------
% FFT
% ----------------------------------------
U = fft(uu)/N;
nmodes = 20;
U20 = zeros(size(U));
U20(1:nmodes) = U(1:nmodes);
U20(end-nmodes+2:end) = U(end-nmodes+2:end);

u_rec = real(ifft(U20*N));

% ----------------------------------------
% Plot original and reconstructed waveform
% ----------------------------------------
% figure;
% plot(tu, uu, 'b', 'linewidth', 2);
% hold on;
% plot(tu, u_rec, 'r--', 'linewidth', 2);
%
% xlabel('Time');
% ylabel('Velocity');
%
% legend('Original', '20-mode reconstruction');
%
% title('FFT Reconstruction using First 20 Modes');
%
% grid on;


% Export
U = fft(uu)/N;
nmodes = 20;

a = zeros(nmodes,1);
b = zeros(nmodes,1);

a(1) = real(U(1));   % mean term

for k = 1:nmodes
    a(k) = 2*real(U(k+1));
    b(k) = -2*imag(U(k+1));
end

fid = fopen('cl.fft','w');

jmax = 1;
imax = 1;
start = 0;
endt = 1;
vis = 0.0035;
rho = 1060;
radread = 1.0;
freqbpm = 60;
womp = 1;

itype = 1;
qmean = Qvol*1000/60;

fprintf(fid,'%d %d %g %g %g %g %g %g %g\n', ...
        jmax, imax, start, endt, vis, rho, radread, freqbpm, womp);

fprintf(fid,'%d %d %g\n', itype, nmodes, qmean);

for k = 1:nmodes
    fprintf(fid,'%e %e\n', a(k), b(k));
end

fclose(fid);



% Reconstruction
u_rec = real(U(1))*ones(size(tu));
u_rec = 0*ones(size(tu));

T = tu(end) - tu(1);
t0 = tu(1);

for k = 1:nmodes
    theta = 2*pi*k*tu;

    u_rec = u_rec ...
          + a(k)*cos(theta) ...
          + b(k)*sin(theta);
end

figure;
plot(tu, uu, 'b', 'linewidth', 2);
hold on;
plot(tu, u_rec, 'r--', 'linewidth', 2);

xlabel('Time');
ylabel('Velocity');

legend('Original', '20-mode reconstruction');

title('Fourier Series Reconstruction');

grid on;

