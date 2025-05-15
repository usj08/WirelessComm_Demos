%% MIMO Beamforming Pattern

%% Settings & Parameters
clc;
clear all;

% Rx Side
% Lr: normalized Rx antenna array length
% nr: # of Rx antennas
% Dr: normalized antenna spacing by wavelength

% Tx Side


% 1. Plot of the Signature Alignment Magnitude
phi = 0:0.001:2*pi; % angular axis
phi_init = pi/3; % arriving direction of signal (main lobe)
omega = cos(phi) - cos(phi_init); % directional cosine ref. to phi_0

nr_list = [2 4 8 16 32 64];
Lr_list = [2 4 8 16];
N_nr = length(nr_list);
N_Lr = length(Lr_list); 

for i_Lr = 1:N_Lr
    figure;
    for i_nr = 1:N_nr
        nr = nr_list(i_nr);
        Lr = Lr_list(i_Lr);
        Dr = Lr / nr;
        subplot(2, N_nr/2, i_nr);
        M = 1 / nr * exp(1j*pi*Dr.*omega*(nr-1)).*(sin(pi*Lr*omega))./(sin(pi*Lr*omega/nr));
        polarplot(phi, abs(M), "Color", 'b');
        title(['$n_r=$ ', num2str(nr), ...
            '$, \Delta_r=$ ', num2str(rats(Lr/nr))], 'Interpreter', 'latex', 'FontSize', 13);
    end
    sgtitle(['Spatial Signature of $L_r=$ ', num2str(Lr), ', Main lobe at $\phi = \pi/$', num2str(pi/phi_init)], 'Interpreter', 'latex', 'FontSize', 15);
end 

%% 2. Decomposition into Orthogonal Spatial Signatures
Lr = 8;
nr = 16;
Dr = Lr / nr;
phi = 0:0.001:2*pi; % angular axis
phi_init = pi/3; % arriving direction of signal (main lobe)
omega = cos(phi) - cos(phi_init); % directional cosine ref. to phi_0

figure;
M = 1 / nr * exp(1j*pi*Dr.*omega*(nr-1)).*(sin(pi*Lr*omega))./(sin(pi*Lr*omega/nr));
polarplot(phi, abs(M), "Color", 'r');

intarr = (0:nr-1)';
er = 1/sqrt(nr) * exp(-1j*2*pi*Dr.*intarr.*omega);
figure;
for i=0:nr-1
    subplot(max(1, fix(nr/4)), max(1, min([nr, 4])), i+1);
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


