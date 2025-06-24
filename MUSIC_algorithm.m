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

Nray = 6; 
Nr = 16;
Nsample = 256;

AoAs = pi * rand(Nray, 1) - pi/2; % range: -pi/2, pi/2
A = zeros(Nr, Nray);
for i = 1:Nray
    A(:, i) = ULA(Nr, k, d, AoAs(i));
end 
rank(A) 

% My own signal. Complex number is also possible. 
% F should be size of (Nrays, Nsamples)
F = diag([1 10 5 1 1 1]) * randn(Nray, Nsample); 

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
theta_interval = pi/Ntheta;
theta_cont = linspace(-pi/2, pi/2, Ntheta);
theta_display = theta_cont * 180 / pi;
array_cont = zeros(Nr, Ntheta);

for i = 1:Ntheta
    array_cont(:, i) = ULA(Nr, k, d, theta_cont(i));
end 

% 2-1. Estimate D = M - N
% a problem to think: how to distinguish N? anyways...
% This time, heuristic.
D = sum(eigvals > eigval_min * 20);
sprintf('estimated D = %d', D)

Nnoise = Nr - D;
noise_eigvecs = eigvecs(:, end-(Nnoise-1):end);

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
xlim([-90, 90]);
ax = gca;
ax.FontSize = 12;

hold on;
for angle = AoAs_deg'
    xline(angle, '--r', sprintf('%.2f $^{\\circ}$', angle), 'Interpreter', 'latex', ...
        'Fontsize', 13);
    hold on;
end 

% Estimating other params?

P_hat = real(pinv(A) * (S - eigval_min * S0) * pinv(A)');
P_hat % Compare with F.




%% Helper Function.
% array steering vector a(theta).
function y = ULA(N, k, d, theta)
    y = zeros(N, 1);
    for n=1:N
        y(n) = 1/sqrt(N) * exp(1j*k*d*((n-1)*sin(theta)));
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