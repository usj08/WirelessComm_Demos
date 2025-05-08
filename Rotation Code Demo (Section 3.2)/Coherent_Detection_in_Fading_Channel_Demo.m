%% Coherent Detection in fading channel Demo
% Author: Seongwook Jung

% Description:
% channel h_l (2,) ~ CN(0, 1) (Rayleigh fading)
% w_l ~ CN(0, N0) (noise)
% Here, N0 = a^2 / SNR
% y_l = h_l x_l + w_l (l=1, 2)
clc; 
clear;

%%
% simulation settings
SNR_dB = -10:5:20;
SNRs = 10 .^ (SNR_dB / 10);
num_iter = 100000;
iterations = 1:1:num_iter;
error_counts = zeros(length(SNR_dB), 1);
error_probs_theory = 1/2 * (1 - sqrt(SNRs ./ (1 + SNRs)));

% TODO: update probability matrix
% TODO: make ML detector
% TODO: plot graphs

SNR_i = 0; 
L = 1;
a = 1;

for SNR = SNRs
    SNR_i = SNR_i + 1;
    
    % codewords of M codes of length L (L, M)
    % each symbol is +- a
    codewords = [-a a];

    for iter = iterations
        % channel settings (h, w)
        h = 1/sqrt(2) * (randn(L,1) + 1i*randn(L,1));
        N0 = a^2 / SNR;
        w = sqrt(N0/2) * (randn(L,1) + 1i*randn(L,1));
        
        % sent signal and received signal
        answer_index = randi([1, 2]);
        x = codewords(1, answer_index);
        y = h .* x + w;

        % detector function
        SS = real(h' / abs(h) * y); % r = |h|x+z
        U = abs(h) * codewords;
        if (SS < 0) 
            index = 1;
        else 
            index = 2;
        end 

        
        if answer_index ~= index
            error_counts(SNR_i) = error_counts(SNR_i) + 1;
        end 
    end
end

% plot the results
error_probs = error_counts ./ num_iter;

plot(SNR_dB, error_probs, '-ro', 'LineWidth', 1.2, 'MarkerSize', 6, 'DisplayName', 'error probability');
hold on;
plot(SNR_dB, error_probs_theory, '--bs', 'LineWidth', 1.2, 'MarkerSize', 6, 'DisplayName', 'error probability (theory)');
legend show;

xlabel('SNR (dB)');
ylabel('error probability p_e');
title('Coherent Detection in Fading Channel (iterations = 10000)');

grid on;
hold off;