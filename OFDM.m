%% OFDM and FFTs
clear; 
clc;
% LTE spec (Only for referencing practica values)
% Subcarrier Spacing: 15kHz
% # of Subcarriers: 12
% sample/symbol: 2048

Ns = 64; % # of OFDM symbols within a signal frame
Nc = 24; % # of subcarriers
Df = 15e3; % subcarrier spacing
Ts = 1/Df; % OFDM symbol duration
B = Nc * Df; % Total BW
Dtau = 1/B;

Nsample = 2048;

% generate OFDM signal
% 1. Comm. data generation
data = ones(Ns, Nc); 

% s(t) = \sum_{m=1}^{Ns} \sum_{n=1}^{Nc} d_{mn} e^{j 2\pi f_n t} rect(~~)
% 1. signal symbol sampled s((m-1)Ts + u) (0 \leq u \leq Ts)
sample_rate = Ts / Nsample;
carrier_fs = Df * (0:1:(Nc-1));

ofdm_signal = zeros(Nsample * Ns, 1);

for m=1:Ns
    for n=1:Nc 
    end
end

%% Delay-Doppler profile
tau = 20e-6; % 20us
doppler_f = 550; % 550Hz
v_delay = exp(-1j*2*pi*Df*tau*(0:(Nc-1)));
v_f = exp(-1j*2*pi*doppler_f*Ts*(0:(Ns-1)));

%
h = 100 * (randn+1j*randn)/sqrt(2); % channel coeff mtx
Y_tilde = h * (v_f' * v_delay);

%
Y_dd = ifft( fft(Y_tilde, [], 1 ), [], 2 ) / sqrt(Ns*Nc);
Y_dd = fftshift(Y_dd,1);

%
fd_axis  = (-Ns/2:Ns/2-1)/(Ns*Ts);   % Hz
tau_axis = (0:Nc-1)*Dtau; % s

%
imagesc(tau_axis*1e6, fd_axis, abs(Y_dd)); axis xy;
xlabel('\tau (µs)'); ylabel('f_D (Hz)');
title('Delay–Doppler profile'); colorbar;

[~,idx] = max(abs(Y_dd(:)));
[k,l] = ind2sub(size(Y_dd), idx);
f_d_est = fd_axis(k);        % ≈ 180 Hz
tau_est = tau_axis(l);  


%% 
taus = [24e-6; 50e-6; 53e-6]; % 20us
doppler_fs = [550; 170; -1200]; % 550Hz
v_delays = exp(-1j*2*pi*Df*taus*(0:(Nc-1)));
v_fs = exp(-1j*2*pi*doppler_fs*Ts*(0:(Ns-1)));

%
h = 10 * (randn+1j*randn)/sqrt(2); % channel coeff mtx
Y_tilde = h * (v_fs' * v_delays);

%
Y_dd = ifft( fft(Y_tilde, [], 1 ), [], 2 ) / sqrt(Ns*Nc);
Y_dd = fftshift(Y_dd,1);

%
fd_axis  = (-Ns/2:Ns/2-1)/(Ns*Ts);   % Hz
tau_axis = (0:Nc-1)*Dtau; % s

%
imagesc(tau_axis*1e6, fd_axis, abs(Y_dd)); axis xy;
xlabel('Delay $\tau$ ($\mu$s)', 'Interpreter', 'latex', 'FontSize', 14); 
ylabel('Doppler $f_D$ (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
title('Delay-Doppler profile', 'Interpreter', 'latex', 'FontSize', 15); 
colorbar;

ax = gca;
ax.FontSize = 13;

[~,idx] = max(abs(Y_dd(:)));
[k,l] = ind2sub(size(Y_dd), idx);
f_d_est = fd_axis(k);        % ≈ 180 Hz
tau_est = tau_axis(l);  


%% Helper Functions