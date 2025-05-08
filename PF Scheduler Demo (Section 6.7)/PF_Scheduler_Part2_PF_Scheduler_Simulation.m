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