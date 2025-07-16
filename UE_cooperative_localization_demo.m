%% UE-Cooperative Sensing (Single Target Localization) DEMO
% Note that AoA can be of +- 180 deg diff because of 1-to-1 mapping of sin
% function. This incurs no problem in triangulation process.
clc;
clear;

% Params
N = 2; % mono-static sensing
Ncl = 1; Nray = 1;
lambda = 0.01; % f_c = 3GHz
k = 2*pi/lambda; % 2pi/lambda. typically lambda = 3*10^8 / N*10^9 ~~ 0.3/N
d = lambda/2; % inter-element spacing. typically lambda/2
num_users = [2 4 8 16 32];
Niter = 128;
Nsample = 2*N; % for MUSIC
xb = 50; yb = 50; % for grid world
theta_spread = 0.1/180 * pi; % 0.1 deg std error
SNR_dB = 20;
SNR = 10^(SNR_dB/10);
AoA_min = -90 * pi/180;
AoA_max = 90 * pi/180;
mean_pos_err = zeros(length(num_users), 1);

for ki = 1:length(num_users)
    K = num_users(ki);
    disp(['K: ', num2str(K)]);
    for ni = 1:Niter
        if (mod(ni, 100) == 0)
            disp(['iteration: ', num2str(ni), '/', num2str(Niter)]);
        end
        %% I. Grid World Generation
        % [-xb, +xb] * [-yb, +yb] rectangular map
        
        % I-1. position and UE target
        pos_target = [2*xb*(rand(1)-0.5) 2*yb*(rand(1)-0.5)];

        % we only assume UE can scan 0~180deg for simplicity
        x_range = xb + pos_target(1);
        pos_UE = [x_range*rand(K, 1)-xb 2*yb*(rand(K, 1)-0.5)];

        dx = pos_target(1) - pos_UE(:,1); dy = pos_target(2) - pos_UE(:,2);

        % I-2. ground-truth AoA
        AoA = atan2(dy, dx);

        %% II. Channel Generation "per UE"
        ch = zeros([N, N, K]);

        for kk = 1:K
            % II-1. Cluster Gain (alpha)
            gamma = N^2;
            clusterGain = (randn([Ncl, Nray]) + 1j * randn([Ncl, Nray])) * sqrt(gamma / (2*Ncl*Nray));
            
            % II-2. AoA and AoD
            % angles ~ Laplacian, mu(phi_cl,theta_cl) and std(phi_spread,theta_spread).
            % Laplacian: f(x|mu, b) = 1/2b exp(-|x-mu|/b) with mean=mu, std=\sqrt{2}b
            theta_cl_tx = AoA(kk);
            theta_tx = zeros([Ncl, Nray]); 
            
            for i = 1:Ncl
                theta_tx(i, :) = laplacian_sample(theta_cl_tx(i), theta_spread/sqrt(2), Nray);
            end 
            
            % II-3. Element Gain(Gamma_r, Gamma_t)
            % Elements gain for Rx/Tx are 0 or 1 according to whether angles belong to 
            % the ranged sector. A sector angle range is simplified to be
            % omni-directional (So set all Gains to 1)
            RxElemGain = ones(Ncl, Nray); TxElemGain = ones(Ncl, Nray); 
            % II-4. Generate one-user channel
            for i = 1:Ncl
                for l = 1:Nray
                    ar = ULA(N, k, d, theta_tx(i, l));
                    ch(:,:,kk) = ch(:,:,kk) + clusterGain(i, l) * RxElemGain(i, l) * TxElemGain(i, l) * (ar * ar');
                end
            end
        end
        
        
        %% III. MUSIC (from sensing results {y(k)}, X = AF + W)
        
        AoA_est = zeros(K, 1);
        Ntheta = 256; % spectrum resolution
        dtheta = (AoA_max - AoA_min) / Ntheta;
        theta_cont = linspace(AoA_min, AoA_max - dtheta, Ntheta);
        array_cont = zeros(N, Ntheta);
        D = 1; % we don't need to estimate D
        
        for i = 1:Ntheta
            array_cont(:, i) = ULA(N, k, d, theta_cont(i));
        end 
        
        for kk = 1:K
            % III-1. Cov matrix generation (one may use spatial smoothing)
            X = zeros(N, Nsample);
            F = (randn(N, Nsample) + 1j * randn(N, Nsample)) * sqrt(SNR/2);
            W = (randn(N, Nsample) + 1j * randn(N, Nsample)) / sqrt(2);
            X = ch(:,:,kk) * F + W;
            S = (X * X') / Nsample;

            % III-2. Eigendecomposition
            S0 = eye(N);
            [eigvecs, eigvals] = eig(S, S0, 'chol');
            [eigvals, idx] = sort(diag(eigvals), 'descend');
            eigvecs = eigvecs(:, idx);
            eigval_min = min(eigvals);
            
            % III-3. Obtain MUSIC Spectrum and estimate AoAs
            noise_eigvecs = eigvecs(:, D+1:end);
            music_spectrum = zeros(Ntheta, 1);
            
            for n = 1:Ntheta
                vv = noise_eigvecs' * array_cont(:, n);
                music_spectrum(n) = 1 / (vv' * vv);
            end

            plot(music_spectrum);

            [~,locs,~,~] = findpeaks(music_spectrum, ...
                'MinPeakProminence',0.7, ...
                'MinPeakDistance',100, ...
                'NPeaks',1);
            if (length(locs) < 1) 
                continue
            end

            AoA_est(kk) = theta_cont(locs);
        end

        % IV. Triangulation
        [p_hat, ~, ~] = LS_bearing(pos_UE, AoA_est*180/pi); 
        
        mean_pos_err(ki) = mean_pos_err(ki) + norm(pos_target - p_hat');

    end
end

mean_pos_err = mean_pos_err / Niter;


% Plot

figure;
plot(num_users, mean_pos_err, 'LineWidth', 1.5, 'Color', 'b', 'LineStyle', '-.');
grid on;
xlabel('number of UE $K$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('mean position error (m)', 'Interpreter', 'latex', 'FontSize', 14);
title('mean position error vs. number of UEs, SNR = 20dB', 'Interpreter', 'latex', 'FontSize', 14);


%% Functions.

function [p_hat, resid, condA] = LS_bearing(S, ang_loc_deg, heading_deg)
% TRIANGULATE A POINT FROM MULTIPLE BEARING MEASUREMENTS (no wrapping).
%
% Inputs
%   S           : Kx2 sensor coordinates [x_k, y_k].
%   ang_loc_deg : Kx1 measured local angles (deg) relative to each sensor's broadside.
%                 e.g., MUSIC AoA in [-90, +90], but *no wrap assumed/required*.
%   heading_deg : Kx1 sensor broadside heading (deg) w.r.t. global +x axis.
%                 If omitted/empty -> zeros(K,1).
%
% Outputs
%   p_hat  : [x; y] LS estimate (same units as S).
%   resid  : RMS perpendicular residual to the K bearing lines (meters).
%   condA  : condition number of A (geometry diagnostic).
%
% Notes
%   - No angle wrapping/mod performed.
%   - If geometry is ill-conditioned (condA >> 1), estimate may be unstable.
%   - Use pinv for robustness.

    if nargin < 3 || isempty(heading_deg), heading_deg = zeros(size(ang_loc_deg)); end

    ang_loc_deg = ang_loc_deg(:);
    heading_deg = heading_deg(:);

    K = size(S,1);
    if numel(ang_loc_deg) ~= K || numel(heading_deg) ~= K
        error('Length mismatch: sensors=%d, ang=%d, heading=%d', K, numel(ang_loc_deg), numel(heading_deg));
    end

    % Global absolute bearing of each line-of-bearing
    theta_glob_deg = heading_deg + ang_loc_deg;  % no wrap needed

    % Line coefficients: a*x + b*y = c
    a = -sind(theta_glob_deg);
    b =  cosd(theta_glob_deg);
    A = [a b];
    c = a.*S(:,1) + b.*S(:,2);

    % Diagnostics
    condA = cond(A);

    % Solve (robust to near-singularity)
    p_hat = pinv(A) * c;   % 2x1

    % Perpendicular residuals
    r = A*p_hat - c;          % signed
    nrm = hypot(a,b);         % should be 1 but compute for safety
    d = r ./ nrm;             % distances
    resid = sqrt(mean(d.^2)); % RMS

end


function y = UPA(W, H, k, d, phi, theta)
    y = zeros(W, H);
    for w=1:W
        for h=1:H
            y(w, h) = exp(1j*k*d*((w-1)*sin(phi)*sin(theta)+(h-1)*cos(theta)));
        end
    end
    y = y * 1/sqrt(W*H);
    y = reshape(y, [W*H, 1]);
end

function y = ULA(N, k, d, theta)
    y = zeros(N, 1);
    for n=1:N
        y(n) = exp(1j*k*d*((n-1)*sin(theta)));
    end
    y = y * 1/sqrt(N);
end

function samples = laplacian_sample(mu, b, N)
    % mu: location parameter
    % b : scale parameter
    % N : number of samples
    u = rand(N, 1) - 0.5;  % Uniform(-0.5, 0.5)
    samples = mu - b * sign(u) .* log(1 - 2 * abs(u));
end
