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