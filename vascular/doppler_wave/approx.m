% time vector (one cycle)
t = linspace(0, 1, 1000);

% systolic forward flow
As = 1.0;
ts = 0.20;
ss = 0.10;
sys = As * exp(-((t - ts).^2) / (2*ss^2));

% early diastolic reversal (negative)
Ar = 0.45;
tr = 0.45;
sr = 0.12;
rev = -Ar * exp(-((t - tr).^2) / (2*sr^2));

% late diastolic forward flow
Ad = 0.15;
td = 0.75;
sd = 0.15;
dia = Ad * exp(-((t - td).^2) / (2*sd^2));

% total waveform
f = sys + rev + dia;

% plot
figure;
plot(t, f, 'k', 'linewidth', 2);
grid on;
xlabel('Time (normalized cardiac cycle)');
ylabel('Flow');
title('Approximate Doppler-like Flow Waveform');

% reference line
hold on;
plot(t, 0*t, '--k');
hold off;
