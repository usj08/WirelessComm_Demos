%% Proportional Fair Scheduler DEMO.
% This Demo contains all parts that are described in Tse book Section 6.7
% (Multiuser Diversity: System Aspect);
% Author: Seongwook Jung

% Part 1. Comparison of Various Schedulers in Asymmetric Environment
% Part 2. Evolution of average throughput per user based on requested rates
% Part 3. Various scenarios: fixed, low-mobility, high-mobility
% Part 4. Opportunistic Beamforming + Proportional Fair Scheduling @ Slow Fading

clc;
clear;

%% 1. Comparison of Various Schedulers in Asymmetric Environment
% Detail: K=4, each with different SNR ([0, 5, 10, 20] dB)
clc;
clear;

% Simulation parameters
slot_len = 1.67e-3;    
Nslots   = 10000;       
K        = 4;           
t_c      = 1000;      

% Asymmetric channel environment scenario
meanSNR_dB = [0; 5; 10; 20];
meanSNR    = 10.^(meanSNR_dB/10);

% Channel model: Rayleigh fading
R = zeros(K, Nslots);
for n = 1:Nslots
    h2 = exprnd(1, [K,1]);              % |h|^2 ~ Exp(1)
    R(:,n) = log2(1 + meanSNR .* h2);   % {R_k}
end

% 1) Round-Robin
serviceSlots_RR = cell(K,1);
for n = 1:Nslots
    k = mod(n-1, K) + 1;
    serviceSlots_RR{k}(end+1) = n;
end
latRR = zeros(K,1);
for k = 1:K
    latRR(k) = mean(diff(serviceSlots_RR{k})) * slot_len;
end

% 2) Max-Rate
serviceSlots_MR = cell(K,1);
for n = 1:Nslots
    [~, k] = max(R(:,n));
    serviceSlots_MR{k}(end+1) = n;
end
latMR = zeros(K,1);
for k = 1:K
    latMR(k) = mean(diff(serviceSlots_MR{k})) * slot_len;
end

% 3) Proportional-Fair
T = ones(K,1) * mean(R(:,1));  % Same initialization T_k[1] across users
alpha = 1/t_c;
serviceSlots_PF = cell(K,1);
for n = 1:Nslots
    PFmetric = R(:,n) ./ T;
    [~, k] = max(PFmetric);
    serviceSlots_PF{k}(end+1) = n;
    % EMA
    T = (1-alpha)*T;
    T(k) = T(k) + alpha * R(k,n);
end
latPF = zeros(K,1);
for k = 1:K
    latPF(k) = mean(diff(serviceSlots_PF{k})) * slot_len;
end

% Plots
users = 1:K;
figure('Name','All Schedulers','NumberTitle','off'); hold on;
bar(users - 0.25, latRR, 0.25, 'FaceColor',[.7 .7 .7], 'DisplayName','RR');
bar(users       , latMR, 0.25, 'FaceColor',[.3 .3 .8], 'DisplayName','Max-Rate');
bar(users + 0.25, latPF, 0.25, 'FaceColor',[.2 .8 .2], 'DisplayName','PF');
xticks(users);
xlabel('User Index','FontSize',14);
ylabel('Inter-service Latency (s)','FontSize',14);
title('Per‐User Latency: RR vs Max-Rate vs PF','FontSize',16);
legend('Location','northwest','FontSize',12);
set(gca,'FontSize',12);
grid on;

figure('Name','RR vs PF','NumberTitle','off'); hold on;
bar(users - 0.125, latRR, 0.25, 'FaceColor',[.7 .7 .7], 'DisplayName','RR');
bar(users + 0.125, latPF, 0.25, 'FaceColor',[.2 .8 .2], 'DisplayName','PF');
xticks(users);
xlabel('User Index','FontSize',14);
ylabel('Inter-service Latency (s)','FontSize',14);
title('Per‐User Latency: RR vs PF Only','FontSize',16);
legend('Location','northwest','FontSize',12);
set(gca,'FontSize',12);
grid on;

% 1) Average Throughput.
thrptRR = zeros(K,1);
thrptMR = zeros(K,1);
thrptPF = zeros(K,1);
for k = 1:K
    thrptRR(k) = sum( R(k, serviceSlots_RR{k}) ) * slot_len / (Nslots*slot_len);
    thrptMR(k) = sum( R(k, serviceSlots_MR{k}) ) * slot_len / (Nslots*slot_len);
    thrptPF(k) = sum( R(k, serviceSlots_PF{k}) ) * slot_len / (Nslots*slot_len);
end

% 2) Total Throughput
totalRR = sum(thrptRR);
totalMR = sum(thrptMR);
totalPF = sum(thrptPF);

% 3) Jain's Fairness Index
slotCount_RR = cellfun(@numel, serviceSlots_RR);   % K×1
slotCount_MR = cellfun(@numel, serviceSlots_MR);
slotCount_PF = cellfun(@numel, serviceSlots_PF);

% Jain’s fairness: (sum(n_k)^2) / (K * sum(n_k.^2))
fairRR = (sum(slotCount_RR)^2) / (K * sum(slotCount_RR.^2));
fairMR = (sum(slotCount_MR)^2) / (K * sum(slotCount_MR.^2));
fairPF = (sum(slotCount_PF)^2) / (K * sum(slotCount_PF.^2));

% 4) Trade-off Scatter Plot
% Trade-off Scatter Plot with custom markers/colors
figure; hold on;
scatter(fairRR, totalRR, 100, 'o', 'MarkerEdgeColor','r',  'MarkerFaceColor','r', 'DisplayName','RR');
scatter(fairMR, totalMR, 100, 's', 'MarkerEdgeColor','b',  'MarkerFaceColor','b', 'DisplayName','Max-Rate');
scatter(fairPF, totalPF, 100, '^', 'MarkerEdgeColor','g',  'MarkerFaceColor','g', 'DisplayName','PF');

xlabel('Jain''s Fairness Index','FontSize',14);
ylabel('Total Throughput (normalized)','FontSize',14);
title('Throughput–Fairness Trade-off','FontSize',16);
legend('Location','best','FontSize',12);
grid on;
set(gca,'FontSize',12);

% 1) Product of throughputs 
P_RR = prod(thrptRR);    % \prod_k T_k for Round-Robin
P_MR = prod(thrptMR);    % for Max-Rate
P_PF = prod(thrptPF);    % for Proportional-Fair

% 2) figure
figure;
bar([P_RR, P_MR, P_PF], 0.6);
set(gca,'XTickLabel',{'RR','Max-Rate','PF'}, 'FontSize',12);
set(gca,'YScale','log')
ylabel('Product of T_k','FontSize',14);
title('Utility across schedulers','FontSize',16);
grid on;



%% 2. Evolution of average throughput per user based on requested rates
slot_len   = 1.67e-3;           % [s]
Nslots     = 10000;             % # of slots
K          = 4;
t_c        = 1000;              % PF EMA Window Size
meanSNR_dB = [0;5;10;20];
meanSNR    = 10.^(meanSNR_dB/10);

% 1) Generate per-slot rates R(k,n)
R = zeros(K,Nslots);
for n = 1:Nslots
    h2 = exprnd(1,[K,1]);               % |h|^2 ~ Exp(1)
    R(:,n) = log2(1 + meanSNR .* h2);
end

% 2) Run PF scheduler, record selection and T history
alpha = 1/t_c;
T     = mean(R(:,1))*ones(K,1);        % Init T_k
T_hist = zeros(K,Nslots);              % {T_k} history per slot
sel    = zeros(1,Nslots);              % selected user idx

for n = 1:Nslots
    metric   = R(:,n) ./ T;
    [~, k_idx]  = max(metric);             % k_idx = argmax_k R_k/T_k
    sel(n)   = k_idx;
    % update EMA
    T        = (1-alpha)*T;
    T(k_idx)    = T(k_idx) + alpha*R(k_idx,n);
    T_hist(:,n) = T;
end

plot_range = 3000:3000+255;   
K = 4;
colors = lines(K);

figure('Name','Time‐Aligned Comparison','NumberTitle','off','Position',[100 100 1400 350]);

% 1) Requested rates R_k[n]
subplot(3,1,1);
hold on;
for k = 1:K
    plot(plot_range, R(k,plot_range), '-','Color',colors(k,:),'LineWidth',0.8);
end
xlabel('Slot index','FontSize',12);
ylabel('Requested rate R_k[n]','FontSize',12);
title('Per-slot Rates','FontSize',14);
legend('User1','User2','User3','User4','Location','northeast','FontSize',9);
xlim([plot_range(1), plot_range(end)]);
grid on;

% 2) Scheduled user per slot (scatter)
subplot(3,1,2);
hold on;
for k = 1:K
    idx = plot_range(sel(plot_range)==k);
    scatter(idx, k*ones(size(idx)), 15, ...
            'MarkerEdgeColor',colors(k,:),...
            'MarkerFaceColor',colors(k,:));
end
xlabel('Slot index','FontSize',12);
ylabel('Scheduled user','FontSize',12);
title('User Selection','FontSize',14);
yticks(1:K); ylim([0.5, K+0.5]);
xlim([plot_range(1), plot_range(end)]);
grid on;

% 3) Evolution of average throughput T_k[n]
subplot(3,1,3);
hold on;
for k = 1:K
    plot(plot_range, T_hist(k,plot_range), '-','Color',colors(k,:),'LineWidth',1);
end
xlabel('Slot index','FontSize',12);
ylabel('EMA throughput T_k[n]','FontSize',12);
title('Avg Throughput Evolution','FontSize',14);
legend('User1','User2','User3','User4','Location','east','FontSize',9);
xlim([plot_range(1), plot_range(end)]);
grid on;


%% 3. Various scenarios: fixed, low-mobility, high-mobility
clear; clc;

% Params
num_iter    = 4;                
t_c         = 1.7;               % PF time constant [s]
slot_len    = 1.67e-3;           % slot duration [s]

SNR_dB      = 0;
gamma       = 10^(SNR_dB/10);
B           = 1.25e6;            % 1.25 MHz
M           = 8;                % Sum‐of‐Sinusoids path count

Nslots      = 2e5;
K_list      = [1 2:2:16];        % sweep user count
Kmax        = max(K_list);

% Doppler & Constants
fc        = 2e9; 
c0        = 3e8;
v_low     = 3/3.6;            
v_high    = 30/3.6;          
v_ext     = 70/3.6;
fD_low    = v_low  * fc/c0;    
fD_high   = v_high * fc/c0;   
fD_fixed  = 2;                
Kappa     = 5;                

t = (0:Nslots-1) * slot_len;   % time vector

sumR_fixed    = zeros(size(K_list));
sumR_low      = zeros(size(K_list));
sumR_high     = zeros(size(K_list));
sumR_ext      = zeros(size(K_list));

for iter = 1:num_iter
    tic;
    % 1) Create Channel
    h_fixed_full = zeros(Kmax, Nslots);
    h_low_full   = zeros(Kmax, Nslots);
    h_high_full  = zeros(Kmax, Nslots);
    h_ext_full = zeros(Kmax, Nslots);
    
    disp(['Iteration = ' num2str(iter) '/' num2str(num_iter)]);

    fprintf('Creating Channels... ');
    for u = 1:Kmax
        % Doppler effect (Clarke's Model) implemented with Jakes (sum-of-sinusoids)
        % 1. Fixed Rician 
        phi_LOS = 2*pi*rand;
        alpha_f = 2*pi*rand(M,1);
        phi_f   = 2*pi*rand(M,1);
        arg_f   = 2*pi*(fD_fixed * cos(alpha_f)) * t + phi_f;
        h_diff  = sum(exp(1j*arg_f),1)/sqrt(M);
        h_fixed_full(u,:) = sqrt(Kappa/(Kappa+1))*exp(1j*phi_LOS) + sqrt(1/(Kappa+1))*h_diff;
        h_fixed_full(u,:) = h_fixed_full(u, :) / sqrt(mean(abs(h_fixed_full(u, :)).^2));

        % 2. Low‐mobility Rayleigh (3km/h)
        alpha_l = 2*pi*rand(M,1);
        phi_l   = 2*pi*rand(M,1);
        arg_l   = 2*pi*(fD_low * cos(alpha_l)) * t + phi_l;
        h_low_full(u,:) = sum(exp(1j*arg_l),1)/sqrt(M);
        h_low_full(u,:) = h_low_full(u, :) / sqrt(mean(abs(h_low_full(u, :)).^2));

        % 3. High‐mobility Rayleigh (30km/h)
        alpha_h = 2*pi*rand(M,1);
        phi_h   = 2*pi*rand(M,1);
        arg_h   = 2*pi*(fD_high * cos(alpha_h)) * t + phi_h;
        h_high_full(u,:) = sum(exp(1j*arg_h),1)/sqrt(M);
        h_high_full(u,:) = h_high_full(u, :) / sqrt(mean(abs(h_high_full(u, :)).^2));

        % 4. High‐mobility Rayleigh (70km/h)
        fD_ext = v_ext * fc / c0;
        alpha_ext = 2*pi*rand(M,1);
        phi_ext   = 2*pi*rand(M,1);
        arg_ext   = 2*pi*(fD_ext * cos(alpha_ext)) * t + phi_ext;
        h_temp_ext = sum(exp(1j*arg_ext),1)/sqrt(M);
        h_ext_full(u,:) = h_temp_ext / sqrt(mean(abs(h_temp_ext).^2));
    end
    disp('Done.');

    if (iter == 1)
        % plot channels for visualization (sanity check, just for K=1)
        startSlot = 1;
        endSlot   = 512;               
        slots     = startSlot:endSlot;  
        t_win     = (slots-1)*slot_len;
        figure;
        plot(t_win, abs(h_fixed_full(1,slots).^2), 'b', 'LineWidth', 1); hold on;
        plot(t_win, abs(h_low_full(1,slots).^2),   'r', 'LineWidth', 1);
        plot(t_win, abs(h_high_full(1,slots).^2),  'g', 'LineWidth', 1);
        xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'Latex');
        ylabel('$|h|^2$', 'FontSize', 14, 'Interpreter', 'Latex');
        title('\textbf{Channel Samples}', 'FontSize', 14, 'Interpreter', 'Latex');
        legend('Fixed Rician','3km/h Rayleigh','30km/h Rayleigh');
        grid on;
    end

    % 2) PF simulation (iterate across user numbers)
    fprintf('Simulating Proportional-Fair... ');
    Rtot_fixed  = zeros(size(K_list));
    Rtot_low    = zeros(size(K_list));
    Rtot_high   = zeros(size(K_list));
    Rtot_ext = zeros(size(K_list));
    for iK = 1:length(K_list)
        K = K_list(iK);

        T_f = log2(1+gamma)*ones(K,1);
        T_l = T_f;
        T_h = T_f;
        T_e = T_f;
        
        sumBits_f = 0;
        sumBits_l = 0;
        sumBits_h = 0;
        sumBits_e = 0;
        
        h_fixed = h_fixed_full(1:K,:);
        h_low   = h_low_full(1:K,:);
        h_high  = h_high_full(1:K,:);
        h_ext = h_ext_full(1:K,:);

        for n = 1:Nslots
            % IS-856 has feedback delay of 3.33ms (2 time slots)
            nn = max(1, n-2);
            hf = h_fixed(:,n);
            hl = h_low(:,n);
            hh = h_high(:,n);
            he = h_ext(:,n);

            hf_predict = h_fixed(:,nn);
            hl_predict = h_low(:,nn);
            hh_predict = h_high(:,nn);
            he_predict = h_ext(:,nn);

            Rf_actual = log2(1 + gamma * abs(hf).^2);
            Rl_actual = log2(1 + gamma * abs(hl).^2);
            Rh_actual = log2(1 + gamma * abs(hh).^2);
            Re_actual = log2(1 + gamma * abs(he).^2);

            Rf_predict = log2(1 + gamma * abs(hf_predict).^2);
            Rl_predict = log2(1 + gamma * abs(hl_predict).^2);
            Rh_predict = log2(1 + gamma * abs(hh_predict).^2);
            Re_predict = log2(1 + gamma * abs(he_predict).^2);

            [~, kf] = max(Rf_predict./T_f);
            [~, kl] = max(Rl_predict./T_l);
            [~, kh] = max(Rh_predict./T_h);
            [~, ke] = max(Re_predict ./ T_e);

            sumBits_f = sumBits_f + Rf_actual(kf)*slot_len*B;
            sumBits_l = sumBits_l + Rl_actual(kl)*slot_len*B;
            sumBits_h = sumBits_h + Rh_actual(kh)*slot_len*B;
            sumBits_e = sumBits_e + Re_actual(ke)*slot_len*B;

            alpha_pf = slot_len/t_c;

            T_f = (1-alpha_pf)*T_f;
            T_l = (1-alpha_pf)*T_l;
            T_h = (1-alpha_pf)*T_h;
            T_e = (1-alpha_pf)*T_e; 

            T_f(kf) = T_f(kf) + alpha_pf*Rf_predict(kf);
            T_l(kl) = T_l(kl) + alpha_pf*Rl_predict(kl);
            T_h(kh) = T_h(kh) + alpha_pf*Rh_predict(kh);
            T_e(ke) = T_e(ke) + alpha_pf*Re_predict(ke);
        end

        totalTime       = Nslots * slot_len;
        Rtot_fixed(iK)  = sumBits_f/1e3/totalTime;
        Rtot_low(iK)    = sumBits_l/1e3/totalTime;
        Rtot_high(iK)   = sumBits_h/1e3/totalTime;
        Rtot_ext(iK)    = sumBits_e / 1e3 / totalTime;
    end
    disp('Done.');
    
    % Total thruput
    sumR_fixed = sumR_fixed + Rtot_fixed;
    sumR_low   = sumR_low   + Rtot_low;
    sumR_high  = sumR_high  + Rtot_high;
    sumR_ext  = sumR_ext  + Rtot_ext;

    toc;
end

% 3) Mean across iterations
avgR_fixed = sumR_fixed / num_iter;
avgR_low   = sumR_low   / num_iter;
avgR_high  = sumR_high  / num_iter;
avgR_ext  = sumR_ext / num_iter;

% 4) Plot
figure;
plot(K_list, avgR_fixed, '-s','LineWidth',2,'DisplayName','Fixed (Rician, \kappa=5)'); hold on;
plot(K_list, avgR_low,   '-o','LineWidth',2,'DisplayName','Low-mobility (Rayleigh, 3 km/h)'); 
plot(K_list, avgR_high,  '-d','LineWidth',2,'DisplayName','High-mobility (Rayleigh, 30 km/h)');
plot(K_list, avgR_ext, '-^','LineWidth',2,'DisplayName','High-mobility (Rayleigh, 70 km/h)');
xlabel('Number of users \(K\)','FontSize',16,'Interpreter','latex');
ylabel('Total throughput (kbps)','FontSize',16,'Interpreter','latex');
title('\textbf{PF Scheduler: Total Throughput vs } $K$','FontSize',18,'Interpreter','latex');
legend('Location','southeast','FontSize',14);
set(gca,'FontSize',14,'XLim',[1 16],'XTick',K_list);
grid on; hold off;


%% 4. Opportunistic Beamforming + PF; Effect on Slow Fading
clear; clc;

% Parameters
K_list   = [1 2 4 8 16 32];
Kmax     = max(K_list);
nt       = 2;                          % Tx antennas (Smaller -> closer performance to coherent one)
gamma    = 10^(0/10);                  % 0 dB SNR
slot_len = 1.67e-3;
t_c      = 1.7;                        % PF adaptation time constant
Nslots   = 1e4;                        % Number of time slots (channel is coherent on Nslots)
alpha_pf = slot_len / t_c;
num_iter = 10;

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
plot(1:Nslots_small, repmat(sqrt(sum(abs(Htest(:, 1)).^2)), 1, Nslots_small), '--', 'Color', 'r', 'LineWidth',1,'DisplayName','User 1 (Before Opp. BF)');
xlabel('Time slots', 'Interpreter','latex');
ylabel('Channel Strength', 'Interpreter','latex');
title('\textbf{Channel Strength after Opportunistic Beamforming}','Interpreter','latex', 'FontSize',14);
legend('Location','southeast','FontSize',12);
grid on;

subplot(2,1,2);
plot(1:Nslots_small, sqrt(gains(1:Nslots_small, 2)), '->', 'LineWidth',1.5, 'Color', 'b', 'DisplayName','User 2 (After Opp. BF)');
hold on;
plot(1:Nslots_small, repmat(sqrt(sum(abs(Htest(:, 2)).^2)), 1, Nslots_small), '--', 'Color', 'b','LineWidth',1,'DisplayName','User 2 (Before Opp. BF)');
xlabel('Time slots', 'Interpreter','latex');
ylabel('Channel Strength', 'Interpreter','latex');
title('\textbf{Channel Strength after Opportunistic Beamforming}','Interpreter','latex', 'FontSize',14);
legend('Location','southeast','FontSize',12);
grid on;


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


%% MISCELLANEOUS (Playground, test)
clc; clear;
Kmax = 1;
Nslots = 1e5;

num_iter    = 20;                
t_c         = 1.6;               % PF time constant [s]
slot_len    = 1.67e-3;           % slot duration [s]
update_slots= round(t_c/slot_len);

SNR_dB      = 0;
gamma       = 10^(SNR_dB/10);
B           = 1.25e6;            % 1.25 MHz
M           = 80;                % Sum‐of‐Sinusoids path count

% Doppler & Constants
fc        = 2e9; 
c0        = 3e8;
v_low     = 2.5/3.6;            
v_high    = 30/3.6;           
fD_low    = v_low  * fc/c0;    
fD_high   = v_high * fc/c0;   
fD_fixed  = 2;                
Kappa     = 5;                

t = (0:Nslots-1) * slot_len;

h_fixed_test = zeros(Kmax, Nslots); 
h_low_test = zeros(Kmax, Nslots);
h_high_test = zeros(Kmax, Nslots);
for u = 1:Kmax
    % Doppler effect (Clarke's Model) implemented with Jakes (sum-of-sinusoids)
    % 1. Fixed Rician 
    phi_LOS = 2*pi*rand;
    alpha_f = 2*pi*rand(M,1);
    phi_f   = 2*pi*rand(M,1);
    arg_f   = 2*pi*(fD_fixed * cos(alpha_f)) * t + phi_f;
    h_diff  = sum(exp(1j*arg_f),1)/sqrt(M);
    h_fixed_test(u,:) = sqrt(Kappa/(Kappa+1))*exp(1j*phi_LOS) + sqrt(1/(Kappa+1))*h_diff;

    % 2. Low‐mobility Rayleigh (3km/h)
    alpha_l = 2*pi*rand(M,1);
    phi_l   = 2*pi*rand(M,1);
    arg_l   = 2*pi*(fD_low * cos(alpha_l)) * t + phi_l;
    h_low_test(u,:) = sum(exp(1j*arg_l),1)/sqrt(M);

    % 3. High‐mobility Rayleigh (30km/h)
    alpha_h = 2*pi*rand(M,1);
    phi_h   = 2*pi*rand(M,1);
    arg_h   = 2*pi*(fD_high * cos(alpha_h)) * t + phi_h;
    h_high_test(u,:) = sum(exp(1j*arg_h),1)/sqrt(M);
end

disp(['Average power of fixed channel: ' num2str(mean(abs(h_fixed_test).^2))]);
disp(['Average power of low-mobility channel: ' num2str(mean(abs(h_low_test).^2))]);
disp(['Average power of high-mobility channel: ' num2str(mean(abs(h_high_test).^2))]);


% plot channels for visualization (sanity check, just for K=1)
startSlot = 1;
endSlot   = 512;               
slots     = startSlot:endSlot;  
t_win     = (slots-1)*slot_len;
figure;
plot(t_win, abs(h_fixed_test(1,slots)), 'b', 'LineWidth', 1); hold on;
plot(t_win, abs(h_low_test(1,slots)),   'r', 'LineWidth', 1);
plot(t_win, abs(h_high_test(1,slots)),  'g', 'LineWidth', 1);
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'Latex');
ylabel('$|h|$', 'FontSize', 14, 'Interpreter', 'Latex');
title('\textbf{Channel Samples}', 'FontSize', 14, 'Interpreter', 'Latex');
legend('Fixed Rician','2.5 km/h Rayleigh','30 km/h Rayleigh', 'FontSize', 12);
grid on;