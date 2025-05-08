%% 4. Opportunistic Beamforming + PF; Effect on Slow Fading
clear; clc;

% Parameters
K_list   = [1 2 4 8 16 32];
Kmax     = max(K_list);
nt       = 2;                          % Tx antennas (Smaller -> closer performance to coherent one)
gamma    = 10^(0/10);                  % 0 dB SNR
slot_len = 1.67e-3;
t_c      = 1.7;                        % PF adaptation time constant
Nslots   = 2e4;                        % Number of time slots (channel is coherent on Nslots)
alpha_pf = slot_len / t_c;
num_iter = 20;

% Preallocate accumulators
SE_obf = zeros(size(K_list));
SE_cbf = 0;

for iter = 1:num_iter
    fprintf('Iter %d/%d\n', iter, num_iter);
    
    % (1) Draw one static block‐fading channel per iteration
    H_full = (randn(nt,Kmax) + 1j*randn(nt,Kmax)) / sqrt(2);  % nt×Kmax

    % (3) Coherent BF (oracle best‐user using full CSI)
    sumRate_cbf = 0;
    for n = 1:Nslots
        rates = zeros(Kmax,1);
        for k = 1:Kmax
            h_k = H_full(:,k);
            q_k = h_k / norm(h_k);             % matched filter
            gain = abs(h_k' * q_k)^2;          % = ||h_k||^2
            rates(k) = log2(1 + gamma * gain);
        end
        sumRate_cbf = sumRate_cbf + rates(1);
    end
    SE_cbf = SE_cbf + sumRate_cbf / Nslots;

    % (2) Opportunistic BF (PF) on the same static channel
    for idx = 1:length(K_list)
        K = K_list(idx);
        T = ones(K,1);
        sumRate_obf = 0;
        H_k = H_full(:,1:K);               % nt×K static sub‐channel
        
        sel_slots = 0;

        for n = 1:Nslots
            q = randn(nt, 1) + 1j * randn(nt, 1);
            q = q / norm(q);
            gain   = abs(H_k' * q).^2;     % K×1
            R_all  = log2(1 + gamma * gain);

            [~, k_star] = max(R_all ./ T);
            %sumRate_obf = sumRate_obf + R_all(k_star);
            if (k_star == 1)
                sel_slots = sel_slots + 1;
                sumRate_obf = sumRate_obf + R_all(k_star);
            end

            % PF update
            T = (1 - alpha_pf) * T;
            T(k_star) = T(k_star) + alpha_pf * R_all(k_star);
        end

        SE_obf(idx) = SE_obf(idx) + sumRate_obf / sel_slots;
    end
end

% Final average over iterations
SE_cbf = (SE_cbf / num_iter) * ones(size(K_list));
SE_obf = SE_obf / num_iter;

% (4) Plot as before
figure;
plot(K_list, SE_cbf, '--ks', 'LineWidth',2, 'DisplayName','Coherent BF');
hold on;
plot(K_list, SE_obf, '-o',  'LineWidth',2, 'DisplayName','Opportunistic BF (PF)');
xlabel('Number of users \(K\)', 'Interpreter','latex');
ylabel('Average throughput (bps/Hz)', 'Interpreter','latex');
title('\textbf{Spectral Efficiency vs. Number of Users}','Interpreter','latex');
legend('Location','southeast','FontSize',12);
grid on;
set(gca,'FontSize',13,'XTick',K_list);

% (5) plot channels
Nslots_small = 50;
Htest = H_full(:, 1:2);
gains = zeros(Nslots, 2);
for n = 1:Nslots
    q = randn(nt, 1) + 1j * randn(nt, 1);
    q = q / norm(q);
    gain = abs(Htest' * q).^2;
    gains(n, :) = gain;
end

% Plot Channel
figure;
subplot(2,1,1);
plot(1:Nslots_small, sqrt(gains(1:Nslots_small, 1)), '-o', 'LineWidth',1.5,'Color', 'r', 'DisplayName','User 1 (After Opp. BF)');
hold on;
plot(1:Nslots_small, repmat(abs(Htest(1, 1)), 1, Nslots_small), '--', 'Color', 'r', 'LineWidth',1,'DisplayName','User 1 (Before Opp. BF)');
xlabel('Time slots', 'Interpreter','latex');
ylabel('Channel Strength', 'Interpreter','latex');
title('\textbf{Channel Strength after Opportunistic Beamforming}','Interpreter','latex', 'FontSize',14);
legend('Location','southeast','FontSize',12);
grid on;

subplot(2,1,2);
plot(1:Nslots_small, sqrt(gains(1:Nslots_small, 2)), '->', 'LineWidth',1.5, 'Color', 'b', 'DisplayName','User 2 (After Opp. BF)');
hold on;
plot(1:Nslots_small, repmat(abs(Htest(1, 2)), 1, Nslots_small), '--', 'Color', 'b','LineWidth',1,'DisplayName','User 2 (Before Opp. BF)');
xlabel('Time slots', 'Interpreter','latex');
ylabel('Channel Strength', 'Interpreter','latex');
title('\textbf{Channel Strength after Opportunistic Beamforming}','Interpreter','latex', 'FontSize',14);
legend('Location','southeast','FontSize',12);
grid on;
