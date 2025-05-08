%% MIMO Beamforming Implementation

%% Settings & Parameters
clc;
clear all;

% Rx Side
Lr = 4; % normalized Rx antenna array length
nr = 1; % # of Rx antennas
Dr = Lr/nr; % normalized antenna spacing by wavelength

% Tx Side


%% 1. Plot of the Signature Alignment Magnitude
phi = 0:0.001:2*pi; % angular axis
phi_init = pi/3;
omega = cos(phi) - cos(phi_init); % directional cosine ref. to phi_0

M = 1 / nr * exp(1j*pi*Dr.*omega*(nr-1)).*(sin(pi*Lr*omega))./(sin(pi*Lr*omega/nr));
figure;
polarplot(phi, abs(M), "Color", 'b');
title(['Lr= ', num2str(Lr), ', nr= ', num2str(nr)]);


%% 2. Decomposition into Orthogonal Spatial Signatures
intarr = (0:nr-1)';
er = 1/sqrt(nr) * exp(-1j*2*pi*Dr.*intarr.*omega);

figure;
for i=0:nr-1
    subplot(fix(nr/4), min([nr, 4]), i+1);
    omega = cos(phi) - i/Lr;
    M = 1 / nr * exp(1j*pi*Dr.*omega*(nr-1)).*(sin(pi*Lr*omega))./(sin(pi*Lr*omega/nr));
    polarplot(phi, abs(M), "Color", 'b');
    rlim([0 1]);
    title(['Basis ', num2str(i)]);
end


%% 3. Angular Windows
modulus = 1/Dr;
dir_cos_arr = 0:(1/Lr):(nr/Lr);

for i=1:nr 
    if (dir_cos_arr(i) == 1)

    end
end


