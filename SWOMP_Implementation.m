%% SW-OMP Implementation
% OBJECTIVE:
% At frequency-selective channel, we exploit SW-OMP to estimate the
% channel.


%% Settings & Parameters
clc;
clear all;

% Tx/Rx side
Lt = 1; Lr = 4; % Ns < NtRF << Nt. Ns < NrRF << Nr. (Default: 1, 4)
Ns = 1; % # of data streams, Here, Ns=NtRF(Lt). (Default: 1)
Nr = 32; Nt = 32; % Default: 32, 32

Nq = 2; % # of Quantization bits that phase shifter will exploit. (Def: 2)

K = 16; % # of OFDM subcarriers (Default: 16)
Gr = 64; Gt = 64; % Grid Size (Dictionary Size, Default: 64, 64)

%% 1. Array response vectors (ULA) and Channel Formulation

% Iteration (NSME Averaged For Various Channels)
iters = 50;
SNRdB = -15:5:10;
NMSEs = zeros(length(SNRdB), 1);
NMSEs_SS = zeros(length(SNRdB), 1);
for iter=1:iters
    % Channel Formulation (Refer to Section II. and IV.)
    L = 4; % # of independent paths
    Nc = 4; % # of delay taps of the channel
    Ts = 1/1760*1e-6; % sampling rate [s] from IEEE. 802.11ad
    
    % Bandlimiting filter p_rc; RC(beta, Ts, t) (Raised-Cosine Filter with roll off=0.8)
    
    % AoA, AoD ~ Unif(0, pi) - OFF-GRID
    % AoA = pi*rand([L, 1]); AoD = pi*rand([L, 1]);

    % AoA, AoD in angles on dict - ON-GRID
    dr = 2*pi/Gr; dt = 2*pi/Gt;
    r_grid = 0:dr:2*pi-dr;
    t_grid = 0:dt:2*pi-dt;
    AoA = r_grid(randi(Gr, [L, 1]));
    AoD = t_grid(randi(Gt, [L, 1]));

    % Delay ~ Unif(0, (Nc-1)Ts)
    tau = (Nc-1)*Ts*rand([L, 1]);
    % (Average) Path Loss
    rho_L = 1e4;
    % Channel Gain s.t. alpha_l ~ CN(0, 1) Per Path
    gain = sqrt(1/2) * (randn([L, 1]) + 1j * randn([L, 1]));
   
    
    % Channel mtx. Per Delay Tap (H0, H1, ... , H_Nc-1)
    Hd = zeros([Nr, Nt, Nc]); % Hd: (Nr, Nt, Nc)
    cf = sqrt(Nt*Nr/L/rho_L);
    for d = 0:Nc-1
        for l = 1:L
            Hd(:, :, d+1) = Hd(:, :, d+1) + cf*gain(l)*ULA(Nr, AoA(l))*ULA(Nt, AoD(l))';
        end
    end
    %RC(beta, Ts, d*Ts-tau(l)). randn()+1 is just a random factor for
    %different d (d=0, 1, ... , Nc-1)
    
    % Store Frequency Version (OFDM With K-subcarriers)
    Hf = zeros(Nr, Nt, K);
    for k = 1:K
        for d = 0:Nc-1
            Hf(:, :, k) = Hf(:, :, k) + 1/sqrt(Nc) * Hd(:, :, d+1) * exp(-1j*2*pi*k*d/K);
        end
    end
    
    % Sanity Check (Check E[|H[k]|^2] ~= NtNr/rho_L)
    [norm(Hf(:, :, :), "Fro")^2/K, Nt*Nr/rho_L]
    
    % Quantized Angle
    dq = 2*pi/(2^Nq);
    AQ = 0:dq:2*pi-dq;
    
    % Dictionary mtx; A_R and A_T
    dr = 2*pi/Gr; dt = 2*pi/Gt;
    nr = 0:Nt-1; nt = 0:Nr-1; 
    r_grid = 0:dr:2*pi-dr;
    t_grid = 0:dt:2*pi-dt;
    Ar_dict = 1/sqrt(Nr)*exp(1j*pi*nr'*cos(r_grid));
    At_dict = 1/sqrt(Nt)*exp(1j*pi*nt'*cos(t_grid)); 
    Psi = kron(conj(At_dict), Ar_dict);
    
    % 2. Simultaneously Weighted-Orthogonal Matching Pursuit (SW-OMP)
    % 2-1. Training Phase
    M = 80; % # of training frames
    
    % Training precoder / combiner for m-th frame
    F_tr = zeros(Nt, Lt, M); W_tr = zeros(Nr, Lr, M);
    
    % Measurement Mtx. independent of subcarriers
    Phi = zeros(M*Lr, Nt*Nr);
    yf = zeros(M*Lr, K); % [y[1], y[2], ... , y[k]]
    
    % t^(m)[k] = 1 for simplicity, first.
    
    % FOR SNRs
    SNRdB = -15:5:10;
    power = 1;
    N = power / rho_L ./ 10 .^ (SNRdB / 10);
    i=1;
    
    for Pnoise = N
        % compute per frame for SW-OMP (BUILDUPs)
        for m=1:M
            Ftemp = 1/sqrt(Nt) * reshape(exp(-1j*AQ(randi([1, 2^Nq], Nt, Lt))), [Nt, Lt]); % F_tr(m)
            Wtemp = 1/sqrt(Nr) * reshape(exp(-1j*AQ(randi([1, 2^Nq], Nr, Lr))), [Nr, Lr]); % W_tr(m)
            q = (2*randi([0, 1], [Lt, 1])-1) * sqrt(power) / sqrt(Ns); % q(m) . t = 1 for simplicity.
            F_tr(:, :, m) = Ftemp;
            W_tr(:, :, m) = Wtemp;
            Phi((m-1)*Lr+1:m*Lr, :) = kron(transpose(q) * transpose(Ftemp), Wtemp');
            for k=1:K
                noise = sqrt(Pnoise) * (1/sqrt(2) * randn([Nr, 1]) + 1j * 1/sqrt(2) * randn([Nr, 1]));
                yf((m-1)*Lr+1:m*Lr, k) = Wtemp'*(Hf(:, :, k) * Ftemp * q) + Wtemp'*noise; % noise ~ CN(0, Pnoise)
            end
        end
        
        % Noise Covariance
        Cw = zeros([M*Lr, M*Lr]); % Noise variance = sigma^2 * blkdiag(W1'W1, W2'W2, ... , W'MWM)
        for m=1:M
            Cw((m-1)*Lr+1:m*Lr, (m-1)*Lr+1:m*Lr) = W_tr(:, :, m)'*W_tr(:, :, m);
        end
        %Cw = Cw * Pnoise; (No need to do this since I didn't multiply Pnoise in the for loop)
        Dw = chol(Cw); % Cholesky Decomposition Cw = Dw' * Dw (Returns Upper-Triangular mtx.)
        Gamma = Phi * Psi; % Equivalent Measurement(Sensing) mtx.
        
        % 2-2. Channel Estimation Phase; Exploiting SW-OMP
        % Complexity: Approximately, O(K*M*Lr*(GrGt))
        cf = zeros(Gt*Gr, K);
        x_supp = [];
        supp_idx = [];
        whiten = inv(Dw)';
        GammaW = whiten*Gamma;
    
        % Initialize residual vectors
        ywf = whiten * yf;
        rf = ywf;
        ths = Pnoise; % epsilon that will stop loop
        MSE = ths + 1;
        k = 1; % counter
        
        while (MSE > ths) % (MSE > ths) is one option, another is (k < L+1)
            % Distributed Correlation
            cf = GammaW' * rf; % cf: (GtGr, K); c[k] : (GtGr, 1)
            [m, p] = max(sum(abs(cf), 2)); % p is argmax
            
            % Update the Current Guess of the common support
            supp_idx = union(supp_idx, p); 
        
            % Weighted-LS
            x_supp = pinv(GammaW(:, supp_idx)) * ywf;
            
            % Update Residual
            rf = ywf - GammaW(:, supp_idx)*x_supp; % Gammaw(:, supp_idx) * x_recon = yw[k] estimator
        
            % MSE
            MSE = 1/(K*M*Lr) * trace(rf' * rf);

            % counter 
            k = k + 1;
        end
        
        % Some Results with SNR
        % Sparse Channel Recovery
        length(supp_idx)  % ONLY FOR WHEN Condition MSE>ths is used
        Hhat = reshape(Psi(:, supp_idx) * x_supp, [Nr, Nt, K]);
        
        % Perf Metrics; NMSE (Normalized Mean Squared Error)
        num = 0; denom = 0;
        for k=1:K
            num = num + norm(Hhat(:, :, k) - Hf(:, :, k), "Fro")^2;
            denom = denom + norm(Hf(:, :, k), "Fro")^2;
        end
        NMSE = num / denom;
        NMSEs(i) = NMSEs(i) + 10 * log(NMSE) / log(10);

        % 2-3. Channel Estimating Phase; SS-SW-OMP + thresholding
        % Let's choose Kp << K for major carriers.
        Kp = 2; % default: 2~4 when K = 16.
        beta = 0.025*Pnoise;
        K_ind = [];
        yf_copy = yf;
        for ki=1:Kp
            [m, p] = max(sum(abs(yf_copy).^2, 1));
            K_ind = union(K_ind, p);
            yf_copy(:, p) = zeros(M*Lr, 1);
        end
    
        cf = zeros(Gt*Gr, Kp);
        x_supp = [];
        supp_idx = [];
        whiten = inv(Dw)';
        GammaW = whiten*Gamma;
    
        % Initialize residual vectors
        ywf = whiten * yf;
        rf = ywf;
        ths = Pnoise; % epsilon that will stop loop
        MSE = ths + 1;
        k2 = 1;
        while (MSE > ths)
            % Distributed Correlation
            cf = GammaW' * rf(:, K_ind); % cf: (GtGr, Kp); c[k] : (GtGr, 1)
            [m, p] = max(sum(abs(cf), 2)); % p is argmax
            
            % Update the Current Guess of the common support
            supp_idx = union(supp_idx, p); 
        
            % Weighted-LS
            x_supp = pinv(GammaW(:, supp_idx)) * ywf;
            
            % Update Residual
            rf = ywf - GammaW(:, supp_idx)*x_supp; % Gammaw(:, supp_idx) * x_recon = yw[k] estimator
        
            % MSE
            MSE = 1/(K*M*Lr) * trace(rf' * rf);
            
            % counter
            k2 = k2 + 1;
        end
    
        % Pruning Low-Power Paths
        Pstar = max(sum(abs(x_supp).^2, 2));
        p_av = sum(abs(x_supp).^2, 2)/K;
        prune_idx = [];
        for i_=1:length(supp_idx)
            if (p_av(i_) > beta * Pstar)
                prune_idx = union(prune_idx, i_);
            end
        end
        
        x_supp = x_supp(prune_idx, :);
        % Sparse Channel Recovery
        Hhat = reshape(Psi(:, supp_idx(prune_idx)) * x_supp, [Nr, Nt, K]);
        
        % Perf Metrics; NMSE (Normalized Mean Squared Error)
        num = 0; denom = 0;
        for k=1:K
            num = num + norm(Hhat(:, :, k) - Hf(:, :, k), "Fro")^2;
            denom = denom + norm(Hf(:, :, k), "Fro")^2;
        end

        NMSE = num / denom;
        NMSEs_SS(i) = NMSEs_SS(i) + 10 * log(NMSE) / log(10);
        i = i + 1;
    end
end

% Averaging over iters * channels
NMSEs = NMSEs / iters;
NMSEs_SS = NMSEs_SS / iters;

figure;
plot(SNRdB, NMSEs, "Marker", "o", "Color", "r", "LineWidth", 1.2);
hold on;
plot(SNRdB, NMSEs_SS, "Marker", "x", "Color", "b", "LineWidth", 1.2);
legend('SW-OMP', 'SS-SW-OMP/Th');
grid on;
xlabel('SNR (dB)');
ylabel('NMSE (dB)');


%% Functions.
function y = ULA(N, theta)
    n = 0:N-1;
    y = exp(1j.*n'*pi*cos(theta)) / sqrt(N);
end

function y = RC(beta, Ts, t)
% Raised-Cosine Filter in time domain. beta = roll-off factor
    if (abs(t) == Ts / (2 * beta))
        y = pi / (4 * Ts) * sinc(1/(2*beta));
    else
        y = 1/Ts * sinc(t/Ts) * cos(pi*beta*t/Ts) / (1 - (2*beta*t/Ts)^2);
    end
end



