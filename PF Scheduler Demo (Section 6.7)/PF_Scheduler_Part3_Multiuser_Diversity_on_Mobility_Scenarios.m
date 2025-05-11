%% 3. Various scenarios: fixed, low-mobility, high-mobility
clear; clc;

% Params
num_iter    = 4;                
t_c         = 1.7;               % PF time constant [s]
slot_len    = 1.67e-3;           % slot duration [s]

SNR_dB      = 0;
gamma       = 10^(SNR_dB/10);
B           = 1.25e6;            % 1.25 MHz
M           = 64;                % Sum‐of‐Sinusoids path count

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