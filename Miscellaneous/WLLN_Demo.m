% WLLN Demo
% For i.i.d {X_i}s, to what value (X_1 X_2 ... X_n)^(1/n) converge?

X = [1, 2, 3];
pX = [1/2 1/4 1/4];
EX = sum(X .* pX); % the product will not converge to this
L = 2^sum(log2(X) .* pX); % exact convergence; 2^E[log X] Due to WLLN

N_exp_list = 1:1:6;
N_list = 10.^N_exp_list;
result_list = ones([1, length(N_list)]);
std_list = ones([1, length(N_list)]);


num_iters = 100;
num_acc = 10;
i = 1;

for N = N_list
    sampled_result_list = ones([1, num_iters]);
    for iter = 1:1:num_iters
        samples = randsample(X, N, true, pX);
        
        % due to exploding value, we accumulate step by step (w/ acc_duration)
        num_subsamples = N / num_acc; 
        for s = 1:1:num_subsamples
            subsamples = samples((s-1)*num_acc+1:s*num_acc);
            sampled_result_list(iter) = sampled_result_list(iter) * power(prod(subsamples), 1/N);
        end
    end
    result_list(i) = mean(sampled_result_list);
    std_list(i) = std(sampled_result_list);
    i = i+1;
end

figure;
errorbar(N_exp_list, result_list, std_list, Marker='>', Color='b', LineWidth=1, DisplayName='n-th root of X_1 X_2 ... X_n');
hold on;
plot(N_exp_list, repelem(EX, length(N_list)), Marker='o', Color='r', DisplayName='E[X]');
plot(N_exp_list, repelem(L, length(N_list)), Marker='x', Color='g', DisplayName='L');

xlabel("log # of samples log_{10} N");
ylabel("n-th root of X_1 X_2 ... X_n");
legend;
grid on;