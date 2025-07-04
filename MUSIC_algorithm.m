%% MUSIC Reproduction
% ref: Schmidt, Ralph. "Multiple emitter location and signal parameter estimation."
% IEEE transactions on antennas and propagation 34.3 (1986): 276-280.
% Here, only considered white noise, but colored noise is easily adaptable;
% Just do whitening process and its over.

% S0 is pre-calculated or pre-known, but this can be done in practice by
% pilot-estimation stage. i.e. with known F.

% I. SYSTEM MODEL
clc;
clear;

lambda = 0.01;
k = 2*pi/lambda; % 2pi/lambda. typically lambda = 3*10^8 / N*10^9 ~~ 0.3/N
d = lambda/2; % inter-element spacing. typically lambda/2

Nray = 4; 
Nr = 32;
Nsample = 2048;

% steering vector continuum a(theta) where \theta \in [AoA_min, AoA_max]
% should be one-to-one mapping to the range of \theta to discern AoAs.
AoA_min = 0 * pi/180;
AoA_max = 180 * pi/180;

AoAs =  (AoA_max-AoA_min) * rand(Nray, 1) + AoA_min; % [AoA_min, AoA_max]
A = zeros(Nr, Nray);
for i = 1:Nray
    A(:, i) = ULA(Nr, k, d, AoAs(i));
end 

% My own signal. Complex number is also possible. 
% F should be size of (Nrays, Nsamples)
powers = randi(4, [Nray, 1]);
F = diag(sqrt(powers)) * randn(Nray, Nsample); 

powers

sigma = 1e-2;
W = sigma * 1/sqrt(2) * (randn(Nr, Nsample) + 1j * randn(Nr, Nsample));
X = A * F + W; % (Nr, Nsamples)


sprintf('rank(A) = %d \nrank(X) = %d', rank(A), rank(X))



% II. Apply MUSIC Algorithm

% 1. Construct S with samples.
S = zeros(Nr, Nr);
S0 = eye(Nr); % S0 should be made prior to analyze eigenstructure.

for n = 1:Nsample 
    S = S + X(:, n) * X(:, n)';
end
S = S ./ Nsample;


% 2. Calculate eigenstructure of S in metric of S_0.
% Solve Sv = lambda S0 v by pre- and post-multiplying S_0^{-1/2}.
% This is actually done by 'eig(A, B)' function, so no need to implement.
[eigvecs, eigvals] = eig(S, S0, 'chol');

% sort eigvecs and eigvals
[eigvals, idx] = sort(diag(eigvals), 'descend');
eigvecs = eigvecs(:, idx);
eigval_min = min(eigvals);

Ntheta = 256;
dtheta = (AoA_max - AoA_min) / Ntheta;
theta_cont = linspace(AoA_min, AoA_max - dtheta, Ntheta);
theta_display = theta_cont * 180 / pi;
array_cont = zeros(Nr, Ntheta);

for i = 1:Ntheta
    array_cont(:, i) = ULA(Nr, k, d, theta_cont(i));
end 

% 2-1. Estimate D = M - N
% a problem to think: how to distinguish N? anyways...
% Exploit the statistics of W. (We already assumed we know S0)

% [idx, C] = kmeans(log(eigvals), 2);
% [~, I] = min(C);
% sum(idx == I)
% D = Nr - sum(idx == I, 'all');

% chi-square LR test
eigvals'
D = 1;
eps = 1e-15;
reject_level = 0.05;
anomaly_flag = 0;

while ~anomaly_flag
    % LR ~ \chi^2 (Nsample(M-D))
    df = Nsample*(Nr-D);
    LR = Nsample / sigma^2 * sum(eigvals(D+1:end));

    pvalue = chi2cdf(LR, df, 'upper');

    if (pvalue > reject_level)
        anomaly_flag = 1;
        break;
    end 
    D = D + 1;
end


sprintf('estimated D = %d', D)

noise_eigvecs = eigvecs(:, D+1:end);

% 3. Evaluate P(theta)
spectrum = zeros(Ntheta, 1);

for n = 1:Ntheta
    vv = noise_eigvecs' * array_cont(:, n);
    spectrum(n) = 1 / (vv' * vv);
end

AoAs_deg = AoAs * 180 / pi;

figure;
plot(theta_display, 10*log10(spectrum), 'LineWidth', 1.25);
grid on;
title('Spectrum against AoA', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Angle of Arrival (AoA) $\theta^\circ$', 'Interpreter', 'latex', ...
    'FontSize', 13);
ylabel('$P_{MU}(\theta)$ (dB)', 'Interpreter', 'latex', ...
    'FontSize', 13);
xlim([AoA_min * 180 / pi, AoA_max * 180 / pi]);
ax = gca;
ax.FontSize = 12;

hold on;
for degree = AoAs_deg'
    xline(degree, '--r', sprintf('%.2f $^{\\circ}$', degree), 'Interpreter', 'latex', ...
        'Fontsize', 13);
    hold on;
end 

% Estimating other params?

P_hat = real(pinv(A) * (S - eigval_min * S0) * pinv(A)');
P_hat % Compare with F.


%% A. Rank Deficiency without Spatial Smoothing

% only limits to single observation
X = A * F(:, 1) + W(:, 1);

% II. Apply MUSIC Algorithm

% 1. Construct S with samples.
S = zeros(Nr, Nr);
S0 = eye(Nr); % S0 should be made prior to analyze eigenstructure.

for n = 1:Nsample 
    S = S + X * X';
end
S = S ./ Nsample;

% 2. Calculate eigenstructure of S in metric of S_0.
% Solve Sv = lambda S0 v by pre- and post-multiplying S_0^{-1/2}.
% This is actually done by 'eig(A, B)' function, so no need to implement.
[eigvecs, eigvals] = eig(S, S0, 'chol');

% sort eigvecs and eigvals
[eigvals, idx] = sort(diag(eigvals), 'descend');
eigvecs = eigvecs(:, idx);
eigval_min = min(eigvals);

Ntheta = 256;
dtheta = (AoA_max - AoA_min) / Ntheta;
theta_cont = linspace(AoA_min, AoA_max - dtheta, Ntheta);
theta_display = theta_cont * 180 / pi;
array_cont = zeros(Nr, Ntheta);

for i = 1:Ntheta
    array_cont(:, i) = ULA(Nr, k, d, theta_cont(i));
end 

% 2-1. Estimate D = M - N
% a problem to think: how to distinguish N? anyways...
% Exploit the statistics of W. (We already assumed we know S0)

% [idx, C] = kmeans(log(eigvals), 2);
% [~, I] = min(C);
% sum(idx == I)
% D = Nr - sum(idx == I, 'all');

% chi-square LR test
eigvals'
D = 1;
eps = 1e-15;
reject_level = 0.05;
anomaly_flag = 0;

while ~anomaly_flag
    % LR ~ \chi^2 (Nsample(M-D))
    df = Nsample*(Nr-D);
    LR = Nsample / sigma^2 * sum(eigvals(D+1:end));

    pvalue = chi2cdf(LR, df, 'upper');

    if (pvalue > reject_level)
        anomaly_flag = 1;
        break;
    end 
    D = D + 1;
end


sprintf('estimated D = %d', D)

noise_eigvecs = eigvecs(:, D+1:end);

% 3. Evaluate P(theta)
spectrum_def = zeros(Ntheta, 1);

for n = 1:Ntheta
    vv = noise_eigvecs' * array_cont(:, n);
    spectrum_def(n) = 1 / (vv' * vv);
end

AoAs_deg = AoAs * 180 / pi;


figure;
plot(theta_display, 10*log10(spectrum_def), 'LineWidth', 1.25);
grid on;
title('Spectrum against AoA (Spatial Smoothing)', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Angle of Arrival (AoA) $\theta^\circ$', 'Interpreter', 'latex', ...
    'FontSize', 13);
ylabel('$P_{MU}(\theta)$ (dB)', 'Interpreter', 'latex', ...
    'FontSize', 13);
xlim([AoA_min * 180 / pi, AoA_max * 180 / pi]);
ax = gca;
ax.FontSize = 12;

hold on;
for degree = AoAs_deg'
    xline(degree, '--r', sprintf('%.2f $^{\\circ}$', degree), 'Interpreter', 'latex', ...
        'Fontsize', 13);
    hold on;
end 

%% B. Rank Augmentation via Spatial Smoothing

% II. Apply MUSIC Algorithm together with SS

% 1. Construct S with samples.
L = 8;
S = zeros(Nr-L+1, Nr-L+1);
S0 = eye(Nr-L+1); % S0 should be made prior to analyze eigenstructure.

for l = 1:L 
    S = S + X(l:Nr-L+l) * X(l:Nr-L+l)';
end

S = S ./ L;

% 2. Calculate eigenstructure of S in metric of S_0.
% Solve Sv = lambda S0 v by pre- and post-multiplying S_0^{-1/2}.
% This is actually done by 'eig(A, B)' function, so no need to implement.
[eigvecs, eigvals] = eig(S, S0, 'chol');

% sort eigvecs and eigvals
[eigvals, idx] = sort(diag(eigvals), 'descend');
eigvecs = eigvecs(:, idx);
eigval_min = min(eigvals);

Ntheta = 256;
dtheta = (AoA_max - AoA_min) / Ntheta;
theta_cont = linspace(AoA_min, AoA_max - dtheta, Ntheta);
theta_display = theta_cont * 180 / pi;
array_cont = zeros(Nr-L+1, Ntheta);

for i = 1:Ntheta
    array_cont(:, i) = ULA(Nr-L+1, k, d, theta_cont(i));
end 

% 2-1. Estimate D = M - N
% a problem to think: how to distinguish N? anyways...
% Exploit the statistics of W. (We already assumed we know S0)

% [idx, C] = kmeans(log(eigvals), 2);
% [~, I] = min(C);
% sum(idx == I)
% D = Nr - sum(idx == I, 'all');

% chi-square LR test
eigvals'
D = 1;
eps = 1e-15;
reject_level = 0.05;
anomaly_flag = 0;

while ~anomaly_flag
    % LR ~ \chi^2 (Nsample(M-D))
    df = Nsample*(Nr-D);
    LR = Nsample / sigma^2 * sum(eigvals(D+1:end));

    pvalue = chi2cdf(LR, df, 'upper');

    if (pvalue > reject_level)
        anomaly_flag = 1;
        break;
    end 
    D = D + 1;
end


sprintf('estimated D = %d', D)

noise_eigvecs = eigvecs(:, D+1:end);

% 3. Evaluate P(theta)
spectrum_ss = zeros(Ntheta, 1);

for n = 1:Ntheta
    vv = noise_eigvecs' * array_cont(:, n);
    spectrum_ss(n) = 1 / (vv' * vv);
end

AoAs_deg = AoAs * 180 / pi;


figure;
plot(theta_display, 10*log10(spectrum), 'LineWidth', 1.25, 'Color', [0 0 0], ...
    'DisplayName', 'MUSIC with plenty data');
grid on;
title(sprintf('Spectrum against AoA (Spatial Smoothing), $M = %d$, $D = %d$', Nr, D), 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Angle of Arrival (AoA) $\theta^\circ$', 'Interpreter', 'latex', ...
    'FontSize', 13);
ylabel('$P_{MU}(\theta)$ (dB)', 'Interpreter', 'latex', ...
    'FontSize', 13);
xlim([AoA_min * 180 / pi, AoA_max * 180 / pi]);
ax = gca;
ax.FontSize = 12;

hold on;
plot(theta_display, 10*log10(spectrum_def), 'LineWidth', 1.25, 'Color', 'r', ...
    'DisplayName', 'MUSIC(1 sample)');

hold on;
plot(theta_display, 10*log10(spectrum_ss), 'LineWidth', 1.25, 'Color', 'b', ...
    'DisplayName', 'MUSIC(1 sample), Spatial Smoothing');

legend;
lgd = legend; 
lgd.FontName = 'Times New Roman';   
lgd.FontSize = 13;   
lgd.AutoUpdate = 'off';

hold on;
for degree = AoAs_deg'
    xline(degree, '--r', sprintf('%.2f $^{\\circ}$', degree), 'Interpreter', 'latex', ...
        'Fontsize', 13);
    hold on;
end 

%% Helper Function.
% array steering vector a(theta).
function y = ULA(N, k, d, theta)
    y = zeros(N, 1);
    for n=1:N
        y(n) = 1/sqrt(N) * exp(1j*k*d*((n-1)*cos(theta)));
    end
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