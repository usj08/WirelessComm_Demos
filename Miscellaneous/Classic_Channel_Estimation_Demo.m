%% Classic Channel Estimation Demo
% Author: Seongwook Jung
clc;
clear;

%%
% Params
Nr = 64;
Nt = 16; 
Pt = 1;

% Possibly Variable Params for Simulation
Np = 12;
SNR_dB_list = -15:5:15; 
MSEs_1 = zeros([length(SNR_dB_list), 1]);
MSEs_2 = zeros([length(SNR_dB_list), 1]);
i = 1; 
num_iter = 100;
Delta = 0.5;

for SNR_dB = SNR_dB_list
    SNR = 10^(SNR_dB/10);
    N0 = Pt / SNR;
    for iter = 1:1:num_iter
        % Channel Generation
        % h_{ij} follows (normalized) Rayleigh fading CN(0, 1)
        H = 1/sqrt(2) * (randn(Nr, Nt)+1j*randn(Nr, Nt));

        % Pilot (Randomly drawn from QPSK symbol)
        QPSK_symbols = sqrt(Pt/2) * [1+1j -1+1j -1-1j 1-1j];
        pilots = QPSK_symbols(randi(4, Nt, Np));
        W = (randn(Nr, Np)+1j*randn(Nr, Np)) * sqrt(N0/2);
        Y = H * pilots + W;
        
        % vectorize for LS (i.e. Y=HP+W is transformed into y=Ah+w)
        y = reshape(Y, [Nr*Np, 1]);
        h = reshape(H, [Nr*Nt, 1]);
        w = reshape(W, [Nr*Np, 1]);
        A = kron(transpose(pilots), eye(Nr));
        
        [h_vec, H_v, A_T, A_R] = virtualChannelModel(H, Nt, Nr, Delta);

        % change into virtual channel model ver;
        A_v = A * kron(conj(A_T), A_R); % 
        y_v = A_v * h_vec + w;

        % 1. LS (Least-Square) ONLY FOR N_p > N_t.
        %h_est1 = (A'*A)\(A'*y);
        %MSEs_1(i) = MSEs_1(i) + mean(abs(h-h_est1).^2) / mean(abs(h).^2);

        % 2. LMMSE with Known Channel Statiscal Model 
        % h_hat = C_h A' * (A C_h A' + N0 I)^{-1} y
        h_est2 = A' * ((A * A' + N0 * eye(Nr*Np))\y);
        MSEs_2(i) = MSEs_2(i) + mean(abs(h-h_est2).^2) / mean(abs(h).^2);
    end

    MSEs_1(i) = MSEs_1(i) / num_iter;
    MSEs_2(i) = MSEs_2(i) / num_iter;
    i = i+1;
end

% MSE vs. SNR_dB
figure;
plot(SNR_dB_list, 10 * log10(MSEs_1), LineWidth=1.0, Marker='>', Color="b", DisplayName="Least Square");
hold on;
plot(SNR_dB_list, 10 * log10(MSEs_2), LineWidth=1.0, Marker='diamond', Color="g", DisplayName="LMMSE");
xlabel("SNR [dB]");
ylabel("NMSE [dB]");
legend();
title("NMSE of Channel Estimation Schemes");
grid on;
fa