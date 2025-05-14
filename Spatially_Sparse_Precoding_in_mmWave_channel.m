%% Orthogonal Matching Pursuit Implementation

%% Settings & Parameters
clc;
clear;

% Tx/Rx side
NtRF = 4;
NrRF = 4; % Ns < NtRF << Nt. Ns < NtRF << Nt
Ns = 2; 
Nt = 64;
Nr = 16; 
Ncl = 8;
Nray = 10;

% If Ncl*Nray << min(Nr, Nt), becomes very accurate (Why?)
k = 2*pi/0.01; % 2pi/lambda. typically lambda = 3*10^8 / N*10^9 ~~ 0.3/N
d = 0.05; % inter-element spacing. typically lambda/2

% A. Rx SIDE 
% 1. Array response vectors (UPA) and Channel Formulation

% Channel Instantiation (p.1501)
% ::::::1. Cluster Gain (alpha)::::::
% Eq.(4) is baseline but refer to the explanation below Eq. (4);
% We follow to set gamma based on this explanation.
gamma = Nt*Nr;
clusterGain = (randn([Ncl, Nray]) + 1j * randn([Ncl, Nray])) * sqrt(gamma / (2*Ncl*Nray));

% ::::::2. AoA, AoD sampling for channel::::::
% for simplicity, alpha is set to meet equal average cluster power
% sigma_i^2 = gamma / NclNray. Summing all power gains equal to NtNr.
% angles ~ Laplacian, mu(phi_cl,theta_cl) and std(phi_spread,theta_spread).
% Laplacian: f(x|mu, b) = 1/2b exp(-|x-mu|/b) with mean=mu, std=\sqrt{2}b
phi_cl_rx = 2*pi*rand([Ncl, 1]); theta_cl_rx = 2*pi*rand([Ncl, 1]);
phi_cl_tx = 2*pi*rand([Ncl, 1]); theta_cl_tx = 2*pi*rand([Ncl, 1]);
phi_rx = zeros([Ncl, Nray]); theta_rx = zeros([Ncl, Nray]); 
phi_tx = zeros([Ncl, Nray]); theta_tx = zeros([Ncl, Nray]); 

% We assume spread is same both in rx and tx for simplicity (but revisable)
phi_spread = 7.5/180 * pi; theta_spread = 7.5/180 * pi;

for i = 1:Ncl
    phi_rx(i, :) = laplacian_sample(phi_cl_rx(i), phi_spread/sqrt(2), Nray);
    theta_rx(i, :) = laplacian_sample(theta_cl_rx(i), theta_spread/sqrt(2), Nray);

    phi_tx(i, :) = laplacian_sample(phi_cl_tx(i), phi_spread/sqrt(2), Nray);
    theta_tx(i, :) = laplacian_sample(theta_cl_tx(i), theta_spread/sqrt(2), Nray);
end 

% ::::::3. Element Gain(Gamma_r, Gamma_t)::::::
% Elements gain for Rx/Tx are 0 or 1 according to whether angles belong to 
% the ranged sector. A sector angle range is simplified to be
% omni-directional (So set all Gains to 1; see p.1509 in paper)
RxElemGain = ones(Ncl, Nray); TxElemGain = ones(Ncl, Nray); 

% ((UNCOMMENT THESE LINES TO DEFINE SECTOR RANGE FOR ELEMENT GAIN.))
% phi_min = -pi; phi_max = pi; theta_min = -pi; theta_max = pi;
% 
% % Wrap angles to bound in [-pi, pi] using modular for convenience in 
% % computing Element Gains.
% phi_rx  = mod(phi_rx + pi, 2*pi) - pi; 
% theta_rx  = mod(theta_rx + pi, 2*pi) - pi;
% phi_tx  = mod(phi_tx + pi, 2*pi) - pi;
% theta_tx  = mod(theta_tx + pi, 2*pi) - pi;
% for i = 1:Ncl
%     for l = 1:Nray
%         if (phi_rx(i, l) > phi_min) && (phi_rx(i, l) < phi_max)
%             RxElemGain(i, l) = 1;
%         end
%         if (phi_tx(i, l) > phi_min) && (phi_tx(i, l) < phi_max)
%             TxElemGain(i, l) = 1;
%         end
%     end
% end

% ::::::4. Generate Channel according to Clustered Channel Model::::::
% Be ware of normalizations.
ch = zeros([Nr, Nt]);
At = zeros([Nt, Ncl*Nray]);
Ar = zeros([Nr, Ncl*Nray]);

Wt = sqrt(Nt); Ht = sqrt(Nt); % W * H = Nt
Wr = sqrt(Nr); Hr = sqrt(Nr);

for i = 1:Ncl
    for l = 1:Nray
        ar = UPA(Wr, Hr, k, d, phi_rx(i, l), theta_rx(i, l));
        at = UPA(Wt, Ht, k, d, phi_tx(i, l), theta_tx(i, l));
        ch = ch + clusterGain(i, l) * RxElemGain(i, l) * TxElemGain(i, l) * ar * at';
        At(:,(Nray)*(i-1)+l) = at; 
        Ar(:,(Nray)*(i-1)+l) = ar;
    end
end

disp('::::::::: SANITY CHECK :::::::::');
disp(['E[|H|^2] = ', num2str(sum(abs(ch).^2, "all"))]);
disp(['Nt*Nr    = ', num2str(Nt*Nr)]);


% 2. Orthogonal Matching Pursuit
[U, S, V] = svd(ch);
Fopt = V(:, 1:Ns);
FRF = zeros(Nt, NtRF);
FBB = zeros(NtRF, Ns);
Fres = Fopt;

for i=1:NtRF
    proj = At' * Fres;
    [m, k] = max(diag(proj * proj'));
    FRF(:, i) = At(:, k);
    FBB = pinv(FRF(:, 1:i)) * Fopt;
    Fres = (Fopt - FRF(:, 1:i) * FBB) / norm(Fopt - FRF(:, 1:i) * FBB, "fro");
end
FBB = sqrt(Ns) * FBB / norm(FRF * FBB, "fro");
F_design = FRF * FBB;

% figure;
% surf(abs(F_design));
% title("Hybrid Precoder from Orthogonal Matching Pursuit");
% 
% figure;
% surf(abs(Fopt));
% title("Hybrid Precoder from Optimal Unconstrained solving");

% B. Rx Side & Signal Recovery
% 1. 
Ndata = 20;
s = (randi(2, [Ns, Ndata])-1) / sqrt(Ns/2); % Note that E[ss']=1/Ns. assumed binary bit stream 1 or 0.

noise_var = 1; rho = 1; % default noise_var: 1. noise_var = 10 means SNR = -10dB
noise = sqrt(noise_var / 2) * (randn([Nr, Ndata]) + 1j * randn([Nr, Ndata]));
y_tx = sqrt(rho) * ch * F_design * s + noise;

% MMSE Optimal Combiner W
Wmmse = (1/sqrt(rho) * inv(F_design' * (ch' * ch) * F_design + noise_var * Ns / rho) * (F_design' * ch'))';
Wres = Wmmse;
Eyy = rho/Ns * ch * (F_design * F_design') * ch' + noise_var * eye(Nr);
WRF = zeros(Nr, 1);
WBB = zeros(NrRF, Ns);

for i=1:NrRF 
    proj = Ar' * Eyy * Wres;
    [m, k] = max(diag(proj * proj'));
    if (i==1)
        WRF = Ar(:, k);
    else
        WRF = [WRF, Ar(:, k)];
    end
    WBB = inv(WRF' * Eyy * WRF) * (WRF' * Eyy * Wmmse);
    Wres = (Wmmse - WRF * WBB) / norm(Wmmse - WRF * WBB, "fro");
end

W_design = WRF * WBB;


% sparse signal recovery
s_recover = abs(W_design' * y_tx);

Nscale = 50; % only for plotting
s_plot = zeros(Ns, Ndata*Nscale);
s_recover_plot = zeros(Ns, Ndata*Nscale);


% Sanity check 
for i = 1:Ns
    s_plot(i,:) = repelem(s(i,:), Nscale) * sqrt(Ns/2); % energy normalized
    s_recover_plot(i,:) = repelem(s_recover(i,:), Nscale) * sqrt(Ns/2);
end

figure;
for i = 1:Ns
    subplot(Ns, 1, i);
    plot(1:Ndata*Nscale, s_plot(i,:), "Color", "green");
    title(["Datastream ", num2str(i)]);
end
% 
% figure;
% for i = 1:Ns
%     subplot(Ns, 1, i);
%     plot(1:Ndata*Nscale, s_recover_plot(i,:), "Color", "blue");
%     title(["Recovered Signal (Before Threshold Detector)", num2str(i)]);
% end

figure;
for i = 1:Ns
    subplot(Ns, 1, i);
    plot(1:Ndata*Nscale, s_recover_plot(i,:)>0.5, "Color", "red");
    title(["Recovered Signal (After Threshold Detector)", num2str(i)]);
end

% 2. Spectral Efficiency
SNR = 10 .^ (-4:0.5:0); % rho / (Ns * noise_var) % dB
SNRdB = 10 * log(SNR) * 1/log(10);
Rres = zeros(length(SNR), 1);
Ropt = zeros(length(SNR), 1);

i = 1;
for snr=SNR
    R = log(det(eye(Ns) + snr * inv(W_design' * W_design) * (W_design' * ch * F_design * F_design' * ch' * W_design))) * 1 / log(2);
    Ro = log(det(eye(Ns) + snr * inv(Wmmse' * Wmmse) *  (Wmmse' * ch * Fopt * Fopt' * ch' * Wmmse))) * 1 / log(2);
    Rres(i) = real(R);
    Ropt(i) = real(Ro);
    i = i + 1;
end

figure;
plot(SNRdB, Rres, "Color", "Blue", "LineWidth", 1.2, "Marker", "o");

hold on;
plot(SNRdB, Ropt, "Color", "Red", "LineWidth", 1.2, "Marker", ">");
xlabel("SNR (dB)");
ylabel("Spectral Efficiency (bits/s/Hz)");
grid on;
legend("Sparse Precoding & Combining", "Optimal Unconstrained Precoding");



%% Functions.
function y = UPA(W, H, k, d, phi, theta)
    y = zeros(W, H);
    for w=1:W
        for h=1:H
            y(w, h) = exp(1j*k*d*((w-1)*sin(phi)*sin(theta)+(h-1)*cos(theta)));
        end
    end
    y = y * 1/sqrt(W*H);
    y = reshape(y, [W*H, 1]);
end

function samples = laplacian_sample(mu, b, N)
    % mu: location parameter
    % b : scale parameter
    % N : number of samples
    u = rand(N, 1) - 0.5;  % Uniform(-0.5, 0.5)
    samples = mu - b * sign(u) .* log(1 - 2 * abs(u));
end


