%% Matched Filter DEMO.


%% Sinusoidal vs. Chirp Signal
clc;
clear;

t_stop = 2e-2; Ntime = 1e3;
t_interval = 5e-2; Ninterval = Ntime * t_stop / t_interval;

dt = t_interval/Ntime;
t_p = 0:dt:t_stop-dt;
t = 0:dt:t_interval-dt;
f0 = 0.1e3;
f1 = 0.6e3;

Nd = 300;
s_chirp = zeros([1, Ntime]); s_sin = zeros([1, Ntime]);
template_chirp = chirp(t_p, f0, t_stop, f1, "linear", 0, "complex");
template_sin = cos(2*pi*(f0+f1)/2*t_p) + 1j * sin(2*pi*(f0+f1)/2*t_p);
s_chirp(Nd:Nd+Ninterval-1) = template_chirp;
s_sin(Nd:Nd+Ninterval-1) = template_sin;

figure;
subplot(2, 1, 1);
plot(t, real(s_chirp), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(s_chirp), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title("original chirp signal");
subplot(2, 1, 2);
plot(t, real(s_sin), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(s_sin), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title("original sinusoidal signal (In-phase, Quadrature)");

% noise addition
Pnoise = 3e-1;
n_chirp = s_chirp + sqrt(Pnoise) .* 1/sqrt(2) * (randn([1, Ntime]) + 1j * randn([1, Ntime]));
n_sin = s_sin + sqrt(Pnoise) .* 1/sqrt(2) * (randn([1, Ntime]) + 1j * randn([1, Ntime]));

figure;
subplot(2, 1, 1);
plot(t, real(n_chirp), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(n_chirp), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title(sprintf("chirp signal with AWGN noise, SNR = %.2f dB", 10 * log10(1/Pnoise)));
subplot(2, 1, 2);
plot(t, real(n_sin), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(n_sin), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title(sprintf("sinusoidal signal with AWGN noise, SNR = %.2f dB", 10 * log10(1/Pnoise)));

% matched filtering

template = conj(fliplr(template_chirp));
m_chirp = conv(n_chirp, template, "valid");
template = conj(fliplr(template_sin));
m_sin = conv(n_sin, template, "valid");

figure;
subplot(2, 1, 1);
plot(abs(m_chirp), "LineWidth", 1.5, "LineStyle", "-", "Color", "red");
grid on;
title(sprintf("matched filtered signal from chirp signal, SNR = %.2f dB", 10 * log10(1/Pnoise)));
subplot(2, 1, 2);
plot(abs(m_sin), "LineWidth", 1.5, "LineStyle", "-", "Color", "red");
grid on;
title(sprintf("matched filtered signal from sinusoidal signal, SNR = %.2f dB", 10 * log10(1/Pnoise)));



%% Sinusoidal vs. Chirp Signal
clc;
clear;

t_stop = 2e-2; Ntime = 1e3;
t_interval = 5e-2; Ninterval = Ntime * t_stop / t_interval;

dt = t_interval/Ntime;
t_p = 0:dt:t_stop-dt;
t = 0:dt:t_interval-dt;
f0 = 0.1e3;
f1 = 0.6e3;

Nd = 300;
Nd_o = 380;
s_chirp = zeros([1, Ntime]); s_sin = zeros([1, Ntime]);
template_chirp = chirp(t_p, f0, t_stop, f1, "linear", 0, "complex");
template_sin = cos(2*pi*(f0+f1)/2*t_p) + 1j * sin(2*pi*(f0+f1)/2*t_p);
s_chirp(Nd:Nd+Ninterval-1) = template_chirp;
s_chirp(Nd_o:Nd_o+Ninterval-1) = s_chirp(Nd_o:Nd_o+Ninterval-1) + 0.8 * template_chirp;
s_sin(Nd:Nd+Ninterval-1) = template_sin;
s_sin(Nd_o:Nd_o+Ninterval-1) = s_sin(Nd_o:Nd_o+Ninterval-1) + 0.8 * template_sin;

figure;
subplot(2, 1, 1);
plot(t, real(s_chirp), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(s_chirp), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title("original chirp signal");
subplot(2, 1, 2);
plot(t, real(s_sin), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(s_sin), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title("original sinusoidal signal (In-phase, Quadrature)");

% noise addition
Pnoise = 5e-1;
n_chirp = s_chirp + sqrt(Pnoise) .* 1/sqrt(2) * (randn([1, Ntime]) + 1j * randn([1, Ntime]));
n_sin = s_sin + sqrt(Pnoise) .* 1/sqrt(2) * (randn([1, Ntime]) + 1j * randn([1, Ntime]));

figure;
subplot(2, 1, 1);
plot(t, real(n_chirp), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(n_chirp), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title(sprintf("chirp signal with AWGN noise, SNR = %.2f dB", 10 * log10(1/Pnoise)));
subplot(2, 1, 2);
plot(t, real(n_sin), "LineWidth", 1.5, "LineStyle", "-");
hold on;
plot(t, imag(n_sin), "LineWidth", 1.5, "LineStyle", ":", "Color", "black");
grid on;
title(sprintf("sinusoidal signal with AWGN noise, SNR = %.2f dB", 10 * log10(1/Pnoise)));

% matched filtering

template = conj(fliplr(template_chirp));
m_chirp = conv(n_chirp, template, "valid");
template = conj(fliplr(template_sin));
m_sin = conv(n_sin, template, "valid");

figure;
subplot(2, 1, 1);
plot(abs(m_chirp), "LineWidth", 1.5, "LineStyle", "-", "Color", "red");
grid on;
title(sprintf("matched filtered signal from chirp signal, SNR = %.2f dB", 10 * log10(1/Pnoise)));
subplot(2, 1, 2);
plot(abs(m_sin), "LineWidth", 1.5, "LineStyle", "-", "Color", "red");
grid on;
title(sprintf("matched filtered signal from sinusoidal signal, SNR = %.2f dB", 10 * log10(1/Pnoise)));