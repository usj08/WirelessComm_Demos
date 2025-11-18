%% OFDM and FFTs
clear; 
clc;
% LTE spec (Only for referencing practica values)
% Subcarrier Spacing: 15kHz
% # of Subcarriers: 12
% sample/symbol: 2048

Ns = 32; % # of OFDM symbols within a signal frame
Nc = 12 * 2^2; % # of subcarriers
Df = 15e3; % subcarrier spacing
Ts = 1/Df; % OFDM symbol duration
B = Nc * Df; % Total BW
Dtau = 1/B;

Nsample = 2048;
noise_pow = 1;

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
doppler_f = 1800; % 1.8kHz
v_delay = exp(-1j*2*pi*Df*tau*(0:(Nc-1)));
v_f = exp(-1j*2*pi*doppler_f*Ts*(0:(Ns-1)));

%
h = 100 * (randn+1j*randn)/sqrt(2); % channel coeff mtx
Y_tilde = h * (v_f' * v_delay);
Y_tilde = Y_tilde + noise_pow * (randn+1j*randn)/sqrt(2);

%
Y_dd = ifft( fft(Y_tilde, [], 1 ), [], 2 ) / sqrt(Ns*Nc);
Y_dd = fftshift(Y_dd,1);

%
fd_axis  = (-Ns/2:Ns/2-1)/(Ns*Ts);   % Hz
tau_axis = (0:Nc-1)*Dtau; % s

% 2D profile
Z = 20*log10(abs(Y_dd)/max(abs(Y_dd(:))));
Z = max(Z, -40);   % 선택: 다이내믹 레인지 제한

[Tau, Fd] = meshgrid(tau_axis, fd_axis);  % X,Y 그리드

figure;
surf(Tau*1e6, Fd, zeros(size(Z)), Z, ...   % Z높이는 0, 색은 Z
     'EdgeColor','none');
view(2);               % 위에서 내려보기 = 2D
shading flat;          % 격자선 숨기기
colormap turbo;        % 색상맵(취향대로)
colorbar;
xlabel('\tau (\mus)');
ylabel('f_D (Hz)');
title('Delay–Doppler (2D surf color only)');
set(gca,'YDir','normal');  % imagesc처럼 아래→위 증가 방향

[~,idx] = max(abs(Y_dd(:)));
[k,l] = ind2sub(size(Y_dd), idx);
f_d_est = fd_axis(k);        % ≈ 180 Hz
tau_est = tau_axis(l);  


% 3D surf
Zlin = abs(Y_dd);                        
ZdB  = 20*log10(Zlin / max(Zlin(:)));     
ZdB  = max(ZdB, -40);                    

% surf grid
[Tau, Fd] = meshgrid(tau_axis, fd_axis);   

figure;
surf(Tau*1e6, Fd, ZdB, ZdB, 'EdgeColor','none');     

shading interp;           
colormap turbo;            
colorbar;                   
view(35, 40);               
xlabel('$\tau (\mu s)$', 'FontSize', 15, 'Interpreter', 'latex');
ylabel('$f_D$ (Hz)', 'FontSize', 15, 'Interpreter', 'latex');
zlabel('Amplitude (dB)', 'FontSize', 15, 'Interpreter', 'latex');
title('Delay–Doppler profile', 'FontSize', 15, 'Interpreter', 'latex');

% 추정된 피크 위치 표시(옵션)
hold on;
plot3(tau_est*1e6, f_d_est, interp2(Tau*1e6, Fd, ZdB, tau_est*1e6, f_d_est), ...
      'ko', 'MarkerFaceColor','w', 'MarkerSize',6);

%% 
taus = [24e-6; 50e-6; 38e-6]; % 20us
doppler_fs = [2200; 170; -1200]; % 550Hz
v_delays = exp(-1j*2*pi*Df*taus*(0:(Nc-1)));
v_fs = exp(-1j*2*pi*doppler_fs*Ts*(0:(Ns-1)));

%
h = 100 * (randn+1j*randn)/sqrt(2); % channel coeff mtx
Y_tilde = h * (v_fs' * v_delays);

%
Y_dd = ifft( fft(Y_tilde, [], 1 ), [], 2 ) / sqrt(Ns*Nc);
Y_dd = fftshift(Y_dd,1);

%
fd_axis  = (-Ns/2:Ns/2-1)/(Ns*Ts);   % Hz
tau_axis = (0:Nc-1)*Dtau; % s

%
colormap turbo;  
imagesc(tau_axis*1e6, fd_axis, 20*log10(abs(Y_dd))); axis xy;
xlabel('Delay $\tau$ ($\mu$s)', 'Interpreter', 'latex', 'FontSize', 14); 
ylabel('Doppler $f_D$ (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
title('Delay-Doppler profile', 'Interpreter', 'latex', 'FontSize', 15); 
colorbar;

ax = gca;
ax.FontSize = 13;

% 3D surf
Zlin = abs(Y_dd);                    
ZdB  = 20*log10(Zlin);     
ZdB  = max(ZdB, -40);                     

% surf grid
[Tau, Fd] = meshgrid(tau_axis, fd_axis);   

figure;
surf(Tau*1e6, Fd, ZdB, ZdB, 'EdgeColor','none');     

shading interp;           
colormap turbo;            
colorbar;                   
view(35, 40);               
xlabel('$\tau (\mu s)$', 'FontSize', 15, 'Interpreter', 'latex');
ylabel('$f_D$ (Hz)', 'FontSize', 15, 'Interpreter', 'latex');
zlabel('Amplitude (dB)', 'FontSize', 15, 'Interpreter', 'latex');
title('Delay–Doppler profile (surf)', 'FontSize', 15, 'Interpreter', 'latex');

%% Helper Functions