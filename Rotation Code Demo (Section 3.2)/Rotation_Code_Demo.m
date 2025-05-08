%% Rotation Code Demo
% Author: Seongwook Jung

% Brief Description of Channel (See Section 3.1.2 in Tse book for details):
% channel h_l ~ CN(0, 1) (Normalized Rayleigh fading)
% w_l ~ CN(0, N0) (noise)
% Here, N0 = a^2 / SNR
% y_l = h_l x_l + w_l (l=1, 2)
clc; 
clear;

%%
% simulation settings
SNR_dB = 0:5:25;
SNRs = 10 .^ (SNR_dB / 10);
thetas = [0, pi/8, 1/2*atan(2), pi/4];
theta_star = atan(1/2); % this is for tightened bound
num_iter = 10^5;
iterations = 1:1:num_iter;
error_counts = zeros(length(SNR_dB), length(thetas));
error_probs_theory = zeros(length(SNR_dB), length(thetas));
error_probs_theory_tighter = zeros(length(SNR_dB), length(thetas));

SNR_i = 0; 
a = 10;

for SNR = SNRs
    SNR_i = SNR_i + 1;
    theta_i = 0;
    N0 = a^2 / SNR;

    for theta = thetas
        theta_i = theta_i + 1;
        % Rotation Matrix R
        R = [cos(theta) -sin(theta);
             sin(theta) cos(theta)];
        % codewords of M codes of length L (L, M)
        % each symbol is +- a
        
        codewords = R * [a a; -a a; -a -a; a -a]';
        diff = 1/a .* (codewords(:, 1) - codewords(:, 2:end)); % (L, M-1)
        coeffs = 1 + SNR/4 * abs(diff.^2);
        coeffs_tighter = coeffs + 1; % this is because Re[h]^2 ~ exp(2)
        prod_dist = abs(diff(1, :) .* diff(2, :)).^2;
        min_prod_dist = min(prod_dist);
     
        % choose the option among 1, 2 for theory value
        % 1. TIGHTER BOUND MODE
        theta_th = atan(1/2);
        
        TERM_B = 1/(coeffs_tighter(1, 1) * (coeffs_tighter(1, 1) + coeffs_tighter(2, 1)));
        if (theta < theta_th)
            TERM_ADD = 1/(coeffs_tighter(1, 3) * (coeffs_tighter(1, 3) + coeffs_tighter(2, 3)));
        else 
            gamma = (2*tan(theta)-1)/(tan(theta)*(tan(theta)+2));
            TERM_ADD_1 = 1/(coeffs_tighter(1, 3) * (coeffs_tighter(1, 3) + coeffs_tighter(2, 3)));
            TERM_ADD_2 = -1/(coeffs_tighter(1, 3) * (coeffs_tighter(1, 3)/gamma + coeffs_tighter(2, 3)));
            TERM_ADD_3 = 1/(coeffs_tighter(1, 2) * (coeffs_tighter(1, 2)/gamma + coeffs_tighter(2, 2)));
            TERM_ADD = TERM_ADD_1 + TERM_ADD_2 + TERM_ADD_3;
        end
        %error_probs_theory_tighter(SNR_i, theta_i) = TERM_B + TERM_ADD;
        % 

        % 2. UNION BOUND MODE
        error_probs_theory(SNR_i, theta_i) = 1/prod(coeffs(:, 1)) + 1/prod(coeffs(:, 2)) + 1/prod(coeffs(:, 3));

        for iter = iterations
            % channel settings (h, w)
            L = 2;
            h = 1/sqrt(2) * (randn(L,1) + 1i*randn(L,1));
            
            w = sqrt(N0/2) * (randn(L,1) + 1i*randn(L,1));
            
            % sent signal and received signal
            answer_index = randi([1, 4]);
            x = codewords(:, answer_index);
            y = h .* x + w;

            % detector function
            U = h .* codewords;
            dists = vecnorm(y - U);
            [~, index] = min(dists);

            % coherent combining with h'/|h|^2 causes BIG error when |h|->0

            if answer_index ~= index
                error_counts(SNR_i, theta_i) = error_counts(SNR_i, theta_i) + 1;
            end
     
        end
    end
end

% plot the results
error_probs = error_counts ./ num_iter;
colors = [1 0 0; 0 1 0; 0 0 1; 0 0.5 0.5; 0.5 0.5 0];
% subplot(1, 2, 1);
for theta_i = 1:length(thetas)
    plot(SNR_dB, error_probs(:, theta_i), '-o', 'LineWidth', 1.2, 'MarkerSize', 6, ...
        'DisplayName', ['error probability ', num2str(thetas(theta_i) / pi), '\pi'], 'Color', colors(theta_i, :));
    hold on;
    plot(SNR_dB, error_probs_theory(:, theta_i), '-->', 'LineWidth', 1.2, 'MarkerSize', 6, ...
        'DisplayName', ['theoretical ', num2str(thetas(theta_i) / pi), '\pi'],'Color', colors(theta_i, :));
end 
set(gca, 'YScale', 'log');
ylim([10^(-4),10]);
legend show;
xlabel('SNR (dB)');
ylabel('error probability p_e (log)');
title('Coherent Detection in Fading Channel (iterations = 10^5)');
grid on;
hold off;

% subplot(1, 2, 2);
% for theta_i = 1:length(thetas)
%     plot(SNR_dB, error_probs(:, theta_i), '-o', 'LineWidth', 1.2, 'MarkerSize', 6, ...
%         'DisplayName', ['error probability ', num2str(thetas(theta_i) / pi), '\pi'], 'Color', colors(theta_i, :));
%     hold on;
%     plot(SNR_dB, error_probs_theory_tighter(:, theta_i), '-->', 'LineWidth', 1.2, 'MarkerSize', 6, ...
%         'DisplayName', ['theoretical ', num2str(thetas(theta_i) / pi), '\pi'],'Color', colors(theta_i, :));
% end 
% set(gca, 'YScale', 'log');
% ylim([10^(-4),10]);
% legend show;
% xlabel('SNR (dB)');
% ylabel('error probability p_e (log)');
% title('Coherent Detection in Fading Channel (tighter bounds, iterations = 10^5)');

grid on;
hold off;