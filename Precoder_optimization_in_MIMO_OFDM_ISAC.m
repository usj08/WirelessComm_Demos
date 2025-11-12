clc; clear;

% Params
K = 16;
M = 64;
Ngrid = 180;
P_tot = 20;
P_per = P_tot / M;
pow_noise = 1e-2;

Niter = 64;
alphas = [0 0.03 0.1 0.3 0.5 1 2];
Na = numel(alphas);

Nch = 100;

% Early-stopping params
tol_abs  = 1e-6;   % absolute tolerance on original objective J
min_iter = 5;      % do at least this many outer iterations before stopping

% Steering matrix A (ULA) on full grid
d = 0.5;
k = 2*pi;
theta_full = linspace(-pi/2, pi/2, Ngrid);
A = zeros(M, Ngrid);
for g = 1:Ngrid
    A(:,g) = ULA(M, k, d, theta_full(g));
end

% Storage over channels
obj_orig_all  = zeros(Niter, Na, Nch);
obj_comm_all  = zeros(Niter, Na, Nch);
obj_sense_all = zeros(Niter, Na, Nch);

% Line-search params
beta = 0.5;
mu0  = 1e-1;

for ch = 1:Nch
    fprintf('Channel %d / %d\n', ch, Nch);

    % Communication channel: i.i.d. Rayleigh
    Hc = (randn(K,M)+1j*randn(K,M))/sqrt(2);

    % Inits
    W0 = (randn(M,K)+1j*randn(M,K))/sqrt(2);
    rad = sqrt(P_per);
    for m = 1:M
        nm = norm(W0(m,:),2);
        if nm > 0
            W0(m,:) = W0(m,:) * (rad / nm);
        else
            W0(m,1) = rad;
        end
    end
    y0 = randn(K,1)/sqrt(K);
    gamma0 = randn(K,1)/sqrt(K);

    for ia = 1:Na
        alpha = alphas(ia);

        W = W0; y = y0; gamma = gamma0;

        [J, comm_obj, sense_pen] = compute_obj(Hc, W, A, alpha, pow_noise);
        obj_orig_all(1, ia, ch)  = J;
        obj_comm_all(1, ia, ch)  = comm_obj;
        obj_sense_all(1, ia, ch) = sense_pen;

        stopped = false;

        for n = 2:Niter
            % FP-style gamma, y updates
            hw_diag = diag(Hc * W);
            b = real(conj(y) .* hw_diag);
            gamma = 0.25 * ( b + sqrt(b.^2 + 4) ).^2 - 1;

            HW = Hc * W;
            denom = sum(abs(HW).^2, 2) + pow_noise;
            y = 2 * sqrt(1 + gamma) .* hw_diag ./ denom;

            % Sensing penalty terms (full grid)
            U = A' * W;
            Z = U * U';
            Off = Z - diag(diag(Z));

            % Comm gradient
            Q = Hc' * diag(abs(y).^2) * Hc;
            G = Hc' * diag(conj(y) .* sqrt(1 + gamma));

            % Total gradient
            grad = -2 * Q * W + 2 * G - 4 * alpha * ( A * (Off * (A' * W)) );

            % Projected gradient ascent with monotone acceptance on original J
            mu = mu0;
            J_old = J;
            while true
                W_try = W + mu * grad;

                % Per-antenna equality projection: ||row_m||_2 = sqrt(P_per)
                for m = 1:M
                    nm = norm(W_try(m,:), 2);
                    if nm > 0
                        W_try(m,:) = W_try(m,:) * (rad / nm);
                    else
                        W_try(m,1) = rad;
                    end
                end

                [J_try, comm_obj, sense_pen] = compute_obj(Hc, W_try, A, alpha, pow_noise);

                if J_try >= J_old
                    W = W_try; J = J_try;
                    break;
                else
                    mu = beta * mu;
                    if mu < 1e-9
                        break;
                    end
                end
            end

            obj_orig_all(n, ia, ch)  = J;
            obj_comm_all(n, ia, ch)  = comm_obj;
            obj_sense_all(n, ia, ch) = sense_pen;

            % Early stopping on original objective
            if n >= min_iter
                dJ = abs(obj_orig_all(n, ia, ch) - obj_orig_all(n-1, ia, ch));
                if dJ <= tol_abs
                    % fill the rest with the last value for smooth plots
                    if n < Niter
                        obj_orig_all(n+1:end, ia, ch)  = J;
                        obj_comm_all(n+1:end, ia, ch)  = comm_obj;
                        obj_sense_all(n+1:end, ia, ch) = sense_pen;
                    end
                    stopped = true;
                    break;
                end
            end
        end

        if stopped
            continue;
        end
    end
end

% Channel-average for plotting
obj_orig  = mean(obj_orig_all, 3);
obj_comm  = mean(obj_comm_all, 3);
obj_sense = mean(obj_sense_all, 3);

% Plot objectives over iterations
cmap = lines(Na);
markers = {'o','s','^','d','v','>','<','p','h'};
linestyles = {'-','--',':','-.'};

figure; hold on;
for ia = 1:Na
    mk = markers{mod(ia-1, numel(markers))+1};
    ls = linestyles{mod(ia-1, numel(linestyles))+1};
    plot(1:Niter, obj_orig(:,ia), 'Color', cmap(ia,:), ...
         'LineStyle', ls, 'Marker', mk, 'LineWidth', 1.6, 'MarkerSize', 4, ...
         'DisplayName', sprintf('\\alpha = {%g}', alphas(ia)));
end
grid on; xlabel('Iteration');
ylabel('Objective Value', 'Interpreter', 'latex');
title(sprintf('Original Objective Value $\\Sigma \\log_2(1+\\mathrm{SINR}_k) - \\alpha\\,\\Phi$ vs. iteration (avg over %d channels)', Nch), 'Interpreter', 'latex');
legend('Location','best'); hold off;

figure; hold on;
for ia = 1:Na
    mk = markers{mod(ia-1, numel(markers))+1};
    ls = linestyles{mod(ia-1, numel(linestyles))+1};
    plot(1:Niter, obj_comm(:,ia), 'Color', cmap(ia,:), ...
         'LineStyle', ls, 'Marker', mk, 'LineWidth', 1.6, 'MarkerSize', 4, ...
         'DisplayName', sprintf('\\alpha = {%g}', alphas(ia)));
end
grid on; xlabel('Iteration');
ylabel('Comm objective Value', 'Interpreter', 'latex');
title( sprintf('Comm Objective Value $\\sum_k \\log_2(1+\\mathrm{SINR}_k)$ vs. iteration (avg over %d channels)', Nch), 'Interpreter','latex' );
legend('Location','best'); hold off;

figure; hold on;
for ia = 1:Na
    mk = markers{mod(ia-1, numel(markers))+1};
    ls = linestyles{mod(ia-1, numel(linestyles))+1};
    plot(1:Niter, obj_sense(:,ia), 'Color', cmap(ia,:), ...
         'LineStyle', ls, 'Marker', mk, 'LineWidth', 1.6, 'MarkerSize', 4, ...
         'DisplayName', sprintf('\\alpha = {%g}', alphas(ia)));
end
grid on; xlabel('Iteration');
ylabel('Sense Penalty Value', 'Interpreter', 'latex');
title(sprintf('Sense Penalty Value $\\Phi = ||offdiag(A^H W W^H A)||_F^2$ vs. iteration (avg over %d channels)', Nch), 'Interpreter', 'latex');
legend('Location','best'); hold off;

% ---------------- helper functions ----------------
function [J, comm_obj, sense_pen] = compute_obj(Hc, W, A, alpha, pow_noise)
    HW     = Hc * W;
    sig2   = abs(diag(HW)).^2;
    interf = sum(abs(HW).^2, 2) - sig2;
    SINR   = sig2 ./ max(interf + pow_noise, eps);
    comm_obj = sum(log2(1 + SINR));
    U = A' * W;
    Z = U * U';
    Off = Z - diag(diag(Z));
    sense_pen = norm(Off,'fro')^2;
    J = comm_obj - alpha * sense_pen;
end

function y = ULA(N, k, d, theta)
    y = zeros(N, 1);
    for n=1:N
        y(n) = (1/sqrt(N)) * exp(-1j*k*d*((n-1)*sin(theta)));
    end
end
