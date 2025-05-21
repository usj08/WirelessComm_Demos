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