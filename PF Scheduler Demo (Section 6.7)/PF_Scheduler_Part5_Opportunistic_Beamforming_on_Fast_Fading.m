%% 5. Effect on Fast Fading
clear; clc;

% Parameters
K_list    = [1 2 4 8 16 32 64];   % user counts
num_K     = length(K_list);
num_iter  = 10;                   % Monte Carlo iterations

t_c_slots = 100;                  % PF time constant ≃100 slots
alpha_pf  = 1 / t_c_slots;    
Nslots    = 5e3;                  % number of fast‐fading slots

gamma     = 10^(0/10);            % 0 dB SNR
kappa     = 10;                   % Ricean K‐factor

% Throughput accumulator
TP = zeros(5, num_K);

for iter = 1:num_iter
  for idx = 1:num_K
    K = K_list(idx);

    % Scenario 1: Fast Rayleigh, N_t=1, PF
    T1 = ones(K,1); sumR1 = 0;
    for n = 1:Nslots
      H = (randn(1,K) + 1j*randn(1,K)) / sqrt(2);
      q = 1;                                 % scalar beam
      g = abs(H' * q).^2;                   
      R = log2(1 + gamma * g);
      [~, k] = max(R ./ T1);
      sumR1 = sumR1 + R(k);
      T1    = (1 - alpha_pf)*T1; T1(k) = T1(k) + alpha_pf*R(k);
    end
    TP(1,idx) = TP(1,idx) + sumR1 / Nslots;

    % Scenario 2: Fast Rayleigh, N_t=2, PF with dumb‐antenna BF
    T2 = ones(K,1); sumR2 = 0;
    for n = 1:Nslots
      H = (randn(2,K) + 1j*randn(2,K)) / sqrt(2);
      alpha_vec = rand(2,1);                   % U[0,1] power weights
      theta_vec = 2*pi * rand(2,1);            % U[0,2π] phases
      q = sqrt(alpha_vec) .* exp(1j*theta_vec);
      q = q / norm(q);                 % normalize beam

      g = abs(H' * q).^2;
      R = log2(1 + gamma * g);
      [~, k] = max(R ./ T2);
      sumR2 = sumR2 + R(k);
      T2 = (1 - alpha_pf)*T2; T2(k) = T2(k) + alpha_pf*R(k);
    end
    TP(2,idx) = TP(2,idx) + sumR2 / Nslots;

    % Scenario 3: Fast Rician, N_t=1, PF with LOS changing each slot
    T3 = ones(K,1); sumR3 = 0;
    for n = 1:Nslots
      los1 = exp(1j*2*pi*rand(1,K)).';         % draw LOS phases per slot
      H = sqrt(kappa/(kappa+1))*los1.' ...
        + (randn(1,K) + 1j*randn(1,K)) / sqrt(2*(kappa+1));
      q = 1;                                % scalar beam
      g = abs(H' * q).^2;
      R = log2(1 + gamma * g);
      [~, k] = max(R ./ T3);
      sumR3 = sumR3 + R(k);
      T3    = (1 - alpha_pf)*T3; T3(k) = T3(k) + alpha_pf*R(k);
    end
    TP(3,idx) = TP(3,idx) + sumR3 / Nslots;

    % Scenario 4: Fast Rician, N_t=2, PF with dumb‐antenna BF and LOS each slot
    T4 = ones(K,1); sumR4 = 0;
    for n = 1:Nslots
      LOS = exp(1j*2*pi*rand(2,K));
      H   = sqrt(kappa/(kappa+1))*LOS ...
          + (randn(2,K) + 1j*randn(2,K)) / sqrt(2*(kappa+1));
      alpha_vec = rand(2,1);
      theta_vec = 2*pi * rand(2,1);
      q = sqrt(alpha_vec) .* exp(1j*theta_vec);
      q = q / norm(q);

      g = abs(H' * q).^2;
      R = log2(1 + gamma * g);
      [~, k] = max(R ./ T4);
      sumR4 = sumR4 + R(k);
      T4 = (1 - alpha_pf)*T4; T4(k) = T4(k) + alpha_pf*R(k);
    end
    TP(4,idx) = TP(4,idx) + sumR4 / Nslots;

    % Scenario 5 (EXTRA): Fast Rician, N_t=4, PF with dumb‐antenna BF and LOS each slot
    T5 = ones(K,1); sumR5 = 0;
    for n = 1:Nslots
      LOS = exp(1j*2*pi*rand(4,K));
      H   = sqrt(kappa/(kappa+1))*LOS ...
          + (randn(4,K) + 1j*randn(4,K)) / sqrt(2*(kappa+1));
      alpha_vec = rand(4,1);
      theta_vec = 2*pi * rand(4,1);
      q = sqrt(alpha_vec) .* exp(1j*theta_vec);
      q = q / norm(q);

      g = abs(H' * q).^2;
      R = log2(1 + gamma * g);
      [~, k] = max(R ./ T5);
      sumR5 = sumR5 + R(k);
      T5 = (1 - alpha_pf)*T5; T5(k) = T5(k) + alpha_pf*R(k);
    end
    TP(5,idx) = TP(5,idx) + sumR5 / Nslots;
  end
end

% Average over iterations
TP = TP / num_iter;

% Plot results
figure; hold on;
styles = {'-s','-o','-d','-^', '->'};
labels = {
  'Rayleigh, Nt=1, PF', ...
  'Rayleigh, Nt=2, PF+Opp. BF', ...
  'Ricean,   Nt=1, PF', ...
  'Ricean,   Nt=2, PF+Opp. BF', ...
  'Ricean,   Nt=4, PF+Opp. BF'
};

for s = 1:5
  plot(K_list, TP(s,:), styles{s}, 'LineWidth',2, 'DisplayName',labels{s});
end
xlabel('Number of users $K$',       'Interpreter','latex');
ylabel('Average throughput (bps/Hz)','Interpreter','latex');
title('\textbf{Fast Fading Total Throughput with Opportunistic BF}', ...
      'Interpreter','latex');
legend('Location','southeast');
grid on;
set(gca,'FontSize',13,'XTick',K_list);