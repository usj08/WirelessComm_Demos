%% Capacity of MIMO fading Channel (Ch 8.1 ~ Ch 8.2)
clear;
clc;

%%
NtNr_list = {[1 1], [2 4], [4 4], [1 8], [8 8]};
SNR_dB_list = -10:5:40;
SNR_list = 10.^(SNR_dB_list/10);

% system model: Rayleigh of all entries being i.i.d CN(0, 1)

C = zeros([length(NtNr_list), length(SNR_dB_list)]);
num_iter = 2048;

for i_NtNr = 1:length(NtNr_list)
    NtNr = NtNr_list{i_NtNr};
    Nt = NtNr(1, 1);
    Nr = NtNr(1, 2);
    for i_SNR = 1:length(SNR_list)
        SNR = SNR_list(i_SNR);
        for iter = 1:num_iter
            H = 1/sqrt(2) * (randn([Nt, Nr]) + 1j * randn([Nt, Nr]));
            [U, S, V] = svd(H, 'econ', 'vector');
            C(i_NtNr, i_SNR) = C(i_NtNr, i_SNR) + sum(log2(1+SNR/Nt * S.^2));
        end
        C(i_NtNr, i_SNR) = C(i_NtNr, i_SNR) / num_iter;
    end
end

% Figure
figure;

% you may add to canvas more (for more plots)
color_canvas = {'k' 'r' 'r', 'b' 'b', }; 
line_canvas = {'--', '-.', '-', '-.', '-', };
marker_canvas = {'None', '>', 'o', '>', 'o'};

for i_NtNr = 1:length(NtNr_list)
    NtNr = NtNr_list{i_NtNr}; Nt = NtNr(1, 1); Nr = NtNr(1, 2);
    plot(SNR_dB_list, C(i_NtNr, :), 'LineWidth', 1.2, 'Color', color_canvas{i_NtNr}, ...
        'LineStyle', line_canvas{i_NtNr}, 'Marker', marker_canvas{i_NtNr}, 'DisplayName', ...
        ['$n_t = $', num2str(Nt), ', $n_r = $', num2str(Nr)]);
    hold on;
end
xlabel('SNR (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Capacity (bits/s/Hz)', 'FontSize', 12, 'Interpreter', 'latex');
title('{\bf Capacity of MIMO Fast fading channel}', 'FontSize', 14, 'Interpreter', 'latex');
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on;
ax.FontName = 'Arial';
legend('FontSize', 13, 'Interpreter','latex');

%% Low SNR Regime
NtNr_list = {[1 1], [1 4], [4 4], [1 8], [8 8]}; % [1 1] is baseline
SNR_dB_list = -30:5:10;
SNR_list = 10.^(SNR_dB_list/10);

% system model: Rayleigh of all entries being i.i.d CN(0, 1)

C = zeros([length(NtNr_list), length(SNR_dB_list)]);
num_iter = 4096;

% baseline: NtNr = [1 1]
i_NtNr = 1; NtNr = NtNr_list{i_NtNr}; Nt = NtNr(1, 1); Nr = NtNr(1, 2);
for i_SNR = 1:length(SNR_list)
    SNR = SNR_list(i_SNR);
    for iter = 1:num_iter
        H = 1/sqrt(2) * (randn([Nt, Nr]) + 1j * randn([Nt, Nr]));
        [U, S, V] = svd(H, 'econ', 'vector');
        C(i_NtNr, i_SNR) = C(i_NtNr, i_SNR) + sum(log2(1+SNR/Nt * S.^2));
    end
    C(i_NtNr, i_SNR) = C(i_NtNr, i_SNR) / num_iter;
end

for i_NtNr = 2:length(NtNr_list)
    NtNr = NtNr_list{i_NtNr};
    Nt = NtNr(1, 1); Nr = NtNr(1, 2);
    for i_SNR = 1:length(SNR_list)
        SNR = SNR_list(i_SNR);
        for iter = 1:num_iter
            H = 1/sqrt(2) * (randn([Nt, Nr]) + 1j * randn([Nt, Nr]));
            [U, S, V] = svd(H, 'econ', 'vector');
            C(i_NtNr, i_SNR) = C(i_NtNr, i_SNR) + sum(log2(1+SNR/Nt * S.^2));
        end
        C(i_NtNr, i_SNR) = C(i_NtNr, i_SNR) / num_iter / C(1, i_SNR);
    end
end

% Figure
figure;

% you may add to canvas more (for more plots)
color_canvas = {'k' 'r' 'r', 'b' 'b', }; 
line_canvas = {'--', '-.', '-', '-.', '-', };
marker_canvas = {'None', '>', 'o', '>', 'o'};

for i_NtNr = 2:length(NtNr_list)
    NtNr = NtNr_list{i_NtNr}; Nt = NtNr(1, 1); Nr = NtNr(1, 2);
    plot(SNR_dB_list, C(i_NtNr, :), 'LineWidth', 1.2, 'Color', color_canvas{i_NtNr}, ...
        'LineStyle', line_canvas{i_NtNr}, 'Marker', marker_canvas{i_NtNr}, 'DisplayName', ...
        ['$n_t = $', num2str(Nt), ', $n_r = $', num2str(Nr)]);
    hold on;
end
xlabel('SNR (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$$\frac{C_{n_t n_r}}{C_{11}}$$', 'Rotation', 0, 'FontSize', 14, 'Interpreter', 'latex');
title('{\bf Capacity of MIMO Fast fading channel (Low SNR regime)}', 'FontSize', 14, 'Interpreter', 'latex');
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on;
ax.FontName = 'Arial';
legend('FontSize', 13, 'Interpreter','latex');

%% Quarter Circle Law
% We focus on n * n channel
n_list = [32 64 128 1024];
nn = length(n_list);
edges = 0:0.1:2;

for i_n = 1:nn
    subplot(2, nn/2, i_n);
    n = n_list(i_n);
    H = 1/sqrt(2 * n) * (randn([n, n]) + 1j * randn([n, n]));
    [U, S, V] = svd(H, 'econ', 'vector');
    hist = histogram(S, 'BinEdges', edges, 'FaceColor', [0 0.8 0.8]);
    set(gca,'FontSize',14);
    grid on;
    title(['$n=$ ', num2str(n)], 'FontSize', 16, 'Interpreter','latex');
end

sgtitle('Empirical Distribution of Singular Value of $$\mathbf{H}/\sqrt{n}$$', ...
    'FontSize', 16, 'Interpreter', 'latex');

%% Large Antenna Array Regime
n_list = [2 4];
SNR_dB_list = -10:5:30;
SNR_list = 10.^(SNR_dB_list/10);

% system model: Rayleigh of all entries being i.i.d CN(0, 1)

C = zeros([length(n_list), length(SNR_dB_list)]);
num_iter = 2048;

for i_n = 1:length(n_list)
    n = n_list(i_n);
    for i_SNR = 1:length(SNR_list)
        SNR = SNR_list(i_SNR);
        for iter = 1:num_iter
            H = 1/sqrt(2) * (randn([n, n]) + 1j * randn([n, n]));
            [U, S, V] = svd(H, 'econ', 'vector');
            C(i_n, i_SNR) = C(i_n, i_SNR) + sum(log2(1+SNR/n * S.^2));
        end
        C(i_n, i_SNR) = C(i_n, i_SNR) / num_iter;
    end
end

F_temp = (sqrt(4*SNR_list + 1)-1).^2;
c_star = 2*log2(1+SNR_list-1/4.*(F_temp)) - log2(exp(1)) ./ (4*SNR_list) .* F_temp;

% Figure
figure;

% you may add to canvas more (for more plots)
color_canvas = {'r' 'b' 'm', }; 
line_canvas = {'-.', '-.', '-', };
marker_canvas = {'>', 'o', '+', };

% c^*(SNR) (baseline approx)
plot(SNR_dB_list, c_star, 'LineWidth', 1.2, 'Color', 'black', ...
    'LineStyle', '--', 'Marker', 'None', 'DisplayName', ...
    ['$c^*(\mathrm{SNR})$']);
hold on;

% exact rate
for i_n = 1:length(n_list)
    n = n_list(i_n);
    plot(SNR_dB_list, C(i_n, :)/n, 'LineWidth', 1.2, 'Color', color_canvas{i_n}, ...
        'LineStyle', line_canvas{i_n}, 'Marker', marker_canvas{i_n}, 'DisplayName', ...
        sprintf('$\\frac{1}{%d} C_{%d%d} $', n, n, n));
    hold on;
end
xlabel('SNR (dB)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Rate (bits/s/Hz)', 'FontSize', 12, 'Interpreter', 'latex');
title('{\bf Capacity approximation: $C_{nn} \approx n ~ c^*(\mathrm{SNR})$}', 'FontSize', 14, 'Interpreter', 'latex');
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on;
ax.FontName = 'Arial';
legend('FontSize', 13, 'Interpreter','latex');


%% MIMO, SIMO, MISO w.r.t n
n_list = 1:1:15;
SNR_dB = 0; % fix SNR state
SNR = 10.^(SNR_dB/10);

% system model: Rayleigh of all entries being i.i.d CN(0, 1)

C = zeros([3, length(n_list)]); % 3 set of capacity; MIMO, SIMO, MISO
num_iter = 2048;

% I. MIMO
for i_n = 1:length(n_list)
    n = n_list(i_n);
    for iter = 1:num_iter
        H = 1/sqrt(2) * (randn([n, n]) + 1j * randn([n, n]));
        [U, S, V] = svd(H, 'econ', 'vector');
        C(1, i_n) = C(1, i_n) + sum(log2(1+SNR/n * S.^2));
    end
    C(1, i_n) = C(1, i_n) / num_iter;
end

% II. SIMO
% y = hx + w, x = scalar
for i_n = 1:length(n_list)
    n = n_list(i_n);
    for iter = 1:num_iter
        h = 1/sqrt(2) * (randn([n, 1]) + 1j * randn([n, 1]));
        [U, S, V] = svd(h, 'econ', 'vector');
        C(2, i_n) = C(2, i_n) + sum(log2(1+SNR * S.^2));
    end
    C(2, i_n) = C(2, i_n) / num_iter;
end


% III. MISO
% y = h*x + w, y = scalar, h : (1, N) vector
for i_n = 1:length(n_list)
    n = n_list(i_n);
    for iter = 1:num_iter
        h = 1/sqrt(2) * (randn([1, n]) + 1j * randn([1, n]));
        [U, S, V] = svd(h, 'econ', 'vector');
        C(3, i_n) = C(3, i_n) + sum(log2(1+SNR/n * S.^2));
    end
    C(3, i_n) = C(3, i_n) / num_iter;
end

% Figure
figure;

% you may add to canvas more (for more plots)
color_canvas = {'r' 'b' 'm', }; 
line_canvas = {'-.', '-.', '-', };
marker_canvas = {'>', 'o', '+', };
struct_canvas = {'MIMO', 'SIMO', 'MISO'};

% exact rate
for i = 1:3
    plot(n_list, C(i, :), 'LineWidth', 1.2, 'Color', color_canvas{i}, ...
        'LineStyle', line_canvas{i}, 'Marker', marker_canvas{i}, ...
        'DisplayName', struct_canvas{i});
    hold on;
end

xlabel('$n$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Rate (bits/s/Hz)', 'FontSize', 12, 'Interpreter', 'latex');
title('{\bf Capacity Comparison Across MIMO, SIMO, MISO}', 'FontSize', 14, 'Interpreter', 'latex');
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
grid on;
ax.FontName = 'Arial';
legend('FontSize', 13, 'Interpreter','latex');
