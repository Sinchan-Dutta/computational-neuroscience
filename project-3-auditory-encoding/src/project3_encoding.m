clc;
close all;
clear;

load('data_cn_project_iii_a26.mat');

%% Visualising the stimulus
figure(1);
plot(Stimulus);
title('Stimulus');

%% Q1 
% Using Matlab Built-in Function
[R, lags] = xcorr(Stimulus, 50, 'unbiased');
figure(2);
plot(lags, R, 'LineWidth', 1.5);
title('Autocorrelation function of s(t)')

% From Scratch
autocorr = zeros(1, 2*50+1);
for i = 0:50
    mult = sum(Stimulus(i+1:end).*Stimulus(1:end-i));
    autocorr(51+i) = mult/20000;
end

for i = 1:50
    autocorr(i) = autocorr(end-i+1);
end

figure(3);
plot(-50:50, autocorr, 'LineWidth', 1.5);
title('Autocorrelation function of s(t)')

%% Q2
% PSTH Calculation
edges = 0:1e-3:20;
PSTH_count = zeros(4, length(edges)-1);
for i = 1:4
    for j = 1:50
        [N, bin] = histcounts(All_Spike_Times{i, j}, edges);
        PSTH_count(i, :) = PSTH_count(i, :) + N;
    end
end

PSTH_rate = PSTH_count / (50 * 1e-3);
PSTH_predict = PSTH_rate(:, 15001:end);

timeBins = edges(1:end-1) + diff(edges)/2;

figure(4);
for k = 1:4
    subplot(2, 2, k);
    plot(timeBins, PSTH_rate(k, :));
    xlabel('Time (s)');
    ylabel('Spiking Rate (Spikes/sec)');
    title(['Neuron No. ', num2str(k)]);
end
sgtitle('PSTH Rate');

%% Q3
% Poisson or Non-Poisson
time_scales = [10, 20, 50, 100, 200, 500];

for i = 1:4
    figure(i+4);
    for j = 1:length(time_scales)
        edges = 0:time_scales(j)*1e-3:20;
        num_bins = length(edges)-1;
        spikeCounts = zeros(50, num_bins);

        for k = 1:50
            [N, ~] = histcounts(All_Spike_Times{i, k}, edges);
            spikeCounts(k, :) = N;
        end

        mean_vals = mean(spikeCounts, 1);
        var_vals = var(spikeCounts, 0, 1);

        subplot(2, 3, j);
        scatter(mean_vals, var_vals, 'filled');
        maximum_val = max([mean_vals, var_vals]);
        hold on;
        plot([0, maximum_val], [0, maximum_val], '--', 'Color', 'red', 'LineWidth', 1.5);
        xlabel('Mean Spike Count');
        ylabel('Variance');
        title([num2str(time_scales(j)), ' ms bins']);
    end
    sgtitle(['Neuron no. ', num2str(i)]);
end

% ISI and CV plot
edges = 0:1e-3:1.5;
ISI_s = cell(4, 1);
for i = 1:4
    for j = 1:50
        ISI_s{i} = [ISI_s{i}, diff(All_Spike_Times{i, j}(:)')];
    end
end

ISI = zeros(4, length(edges)-1);
for i = 1:4
    [N, ~] = histcounts(ISI_s{i}, edges);
    ISI(i, :) = N;
end

means = zeros(1, 4);
variances = zeros(1, 4);
for i = 1:4
    means(1, i) = mean(ISI_s{i});
    variances(1, i) = var(ISI_s{i});
end
CV = sqrt(variances) ./ means;

timeBins = edges(1:end-1) + diff(edges)/2;
figure(9);
for k = 1:4
    subplot(2, 2, k);
    bar(timeBins, ISI(k, :), 'EdgeColor', 'none');
    xlabel('Spike Interval (s)');
    ylabel('Number of Spikes');
    title(['Neuron No. ', num2str(k)]);
    legend(['Mean = ', num2str(means(k)), '  Variance = ', num2str(variances(k)), '  CV = ', num2str(CV(k))]);
end
sgtitle("ISI");

figure(10);
bar(1:4, CV, 0.2, 'EdgeColor', 'b', 'FaceColor', 'b', 'DisplayName', 'CV');
yline(1, '--', 'Color', 'r', displayName = 'y = 1');
xlabel('Neuron No. ');
ylabel('CV');
legend();
title('CV Plot');

%% Q4 (Spike Triggered Average)
STA = zeros(4, 50);

for k = 1:4
    spike_indices = cell(1);
    for run = 1:50
        spike_indices{1} = [spike_indices{1}, ceil(1000 * All_Spike_Times{k, run}(:)')];
    end
    spike_indices{1} = spike_indices{1}(spike_indices{1} >= 50 & spike_indices{1} <= 15000);
    no_of_spikes = length(spike_indices{1});
    for position = spike_indices{1}
        STA(k, :) = STA(k, :) + Stimulus(position - 49 : position);
    end
    STA(k, :) = flip(STA(k, :)) / no_of_spikes;
end

figure(11);
for k = 1:4
    subplot(2, 2, k);
    plot(0:49, STA(k, :), 'LineWidth', 1.5);
    title(['Neuron No. ', num2str(k)]);
end
sgtitle('Raw STA');

% Correction
R_tau = autocorr(51:100);
Cxx = toeplitz(R_tau);

h_corrected = (Cxx \ STA')';

figure(12);
for k = 1:4
    subplot(2, 2, k);
    plot(0:49, h_corrected(k, :), 'LineWidth', 1.5);
    title(['Neuron No. ', num2str(k)]);
end
sgtitle('Corrected filter');

%% Q5
x = Stimulus(1:15000);
estimated_y = zeros(4, 14951);
for k = 1:4
    y_hat = conv(x, h_corrected(k, :));
    estimated_y(k, :) = y_hat(50:15000);
end

figure(13);
for k = 1:4
    y_hat = estimated_y(k, :);
    lambda = PSTH_rate(k, 50:15000);
    subplot(2, 2, k);
    scatter(y_hat, lambda, "filled");
    title(['Neuron No. ', num2str(k)]);
    xlabel('Estimated y(t)');
    ylabel('Measured \lambda(t)');
end
sgtitle('Measured \lambda(t) vs Estimated y(t)');

% Binning
sigmoid = @(p, x) p(1) ./ (1 + exp(-p(2)*(x - p(3)))) + p(4);
opts = optimoptions('lsqcurvefit', 'Display', 'off');
% sigmoid = @(p, x) p(1) + p(2) * x + p(3) * x.^2 + p(4) * x.^3;
curves = cell(4, 1);
figure(14);
num_bins = 30;
for k = 1:4
    y_hat = estimated_y(k, :);
    lambda = PSTH_rate(k, 50:15000);
    [bin_assigned, edges] = discretize(y_hat, num_bins);
    mean_y = zeros(1, num_bins);
    mean_lambda = zeros(1, num_bins);

    for b = 1:num_bins
        idx = (bin_assigned == b);
        if any(idx)
            mean_y(b) = mean(y_hat(idx));
            mean_lambda(b) = mean(lambda(idx));
        end
    end
    initial_guess = [1 1 0 0];
    curves{k} = lsqcurvefit(sigmoid, initial_guess, mean_y, mean_lambda, [], [], opts);

    subplot(2, 2, k);
    scatter(mean_y, mean_lambda, "filled");
    hold on;
    plot(mean_y, sigmoid(curves{k}, mean_y), '--');
    title(['Neuron No. ', num2str(k)]);
    xlabel('Mean y');
    ylabel('Mean \lambda');
end
sgtitle('Mean \lambda vs Mean y');

%% Q6
x_ = Stimulus(15001:end);
estimated_y = zeros(4, 4951);
actual_y = zeros(4, 4951);
figure(15);
for k = 1:4
    y_hat = conv(x_, h_corrected(k, :));
    estimated_y(k, :) = sigmoid(curves{k}, y_hat(50:5000));
    actual_y(k, :) = PSTH_rate(k, 15050:end);
    R = corrcoef(estimated_y(k, :), actual_y(k, :));
    r = R(1, 2);
    
    subplot(2, 2, k);
    scatter(actual_y(k, :), estimated_y(k, :), 'filled');
    xlabel('Actual Measured Value');
    ylabel('Predicted Value');
    text(0.65, 0.9, ['r^2 = ', num2str(r^2)], 'Units', 'normalized');
    title(['Neuron', num2str(k)]);
end
sgtitle('Predicted vs Actual');

%% Pruning
performance = zeros(4, 50);
h_final = zeros(4, 50);
for k = 1:4
    [sorted_h, idx] = sort(abs(h_corrected(k, :)), 'ascend');
    h = h_corrected(k, :);
    for i = 0:49
        if i >= 1
            h(idx(1:i)) = 0;
        end
        y_hat = conv(x, h);
        estimated_y_ = y_hat(50:15000);
        lambda = PSTH_rate(k, 50:15000);

        initial_guess = [max(lambda) 1 mean(estimated_y_) min(lambda)];
        curves{k} = lsqcurvefit(sigmoid, initial_guess, estimated_y_, lambda, [], [], opts);

        y_hat = conv(x_, h);
        estimated_y(k, :) = sigmoid(curves{k}, y_hat(50:5000));
        actual_y(k, :) = PSTH_rate(k, 15050:end);
        R = corrcoef(estimated_y(k, :), actual_y(k, :));
        r = R(1, 2);
        performance(k, 50-i) = r^2;
    end
    % From observation I decided to keep 44, 35, 24, 43 parameters respectively for h(t)
    keepParameters = [44, 35, 24, 43];
    h_final(k, :) = h_corrected(k, :);
    h_final(k, idx(1:50 - keepParameters(k))) = 0;
end

figure(16);
for k = 1:4
    subplot(2, 2, k);
    plot(1:50, performance(k, :));
    xlabel('No. of parameters');
    ylabel('r^2');
    title(['Neuron', num2str(k)]);
end
sgtitle('No. of parameters vs performance');

figure(17);
for k = 1:4
    subplot(2, 2, k);
    plot(1:50, h_final(k, :));
    title(['Neuron', num2str(k)]);
end
sgtitle('Filters with pruning');

figure(18);
for k = 1:4
    Nfft = 512;
    Fs = 1000;
    h_fft = fft(h_final(k, :), Nfft);
    H_mag = abs(h_fft(1:Nfft/2+1));
    f_s = (0:Nfft/2) * (Fs / Nfft);
    subplot(2, 2, k);
    plot(f_s, H_mag, 'LineWidth', 1.5);
    xlabel('Freq (Hz)');
    ylabel('Magnitude');
    title(['Neuron', num2str(k)]);
end
sgtitle('Filter Response');

%% Maximally Informative Dimension (Track A)
len_v = 100;
mid_bins = 25;
num_epochs = 20;
n_restarts = 5;
rng(0);

train_stimulus = Stimulus(1:15000);
num_of_train_samples = length(train_stimulus) - len_v + 1;
X = zeros(num_of_train_samples, len_v);
for i = 1:num_of_train_samples
    X(i, :) = flip(train_stimulus(i:i + len_v - 1));
end

test_stimulus = Stimulus(15001:end);
num_of_test_samples = length(test_stimulus) - len_v + 1;
X_test = zeros(num_of_test_samples, len_v);
for i = 1:num_of_test_samples
    X_test(i, :) = flip(test_stimulus(i:i + len_v - 1));
end

spike_indices = cell(4, 1);
for k = 1:4
    for run = 1:50
        spike_times = ceil(1000 * All_Spike_Times{k, run}(:)');
        spike_times = spike_times(spike_times >= len_v & spike_times <= 15000);
        spike_indices{k} = [spike_indices{k}, spike_times - (len_v - 1)];
        % Row i of X holds the stimulus ending at ms (i + len_v - 1), so a spike in ms n belongs to row n - len_v + 1.
    end
end

STA100 = zeros(4, len_v);
for k = 1:4
    STA100(k, :) = mean(X(spike_indices{k}, :), 1);
end

v_MID = zeros(4, len_v);
MI_history = zeros(4, num_epochs);
MI_sta = zeros(1, 4);
restart_MI = zeros(4, n_restarts);
restart_cos = zeros(4, n_restarts);
for k = 1:4
    MI_sta(k) = mid_grad(STA100(k, :) / norm(STA100(k, :)), X, spike_indices{k}, mid_bins);
    [v_MID(k, :), MI_history(k, :)] = mid_run(STA100(k, :), X, spike_indices{k}, mid_bins, num_epochs);
    for r = 1:n_restarts
        [vr, hist_r] = mid_run(randn(1, len_v), X, spike_indices{k}, mid_bins, num_epochs);
        restart_MI(k, r) = hist_r(end);
        restart_cos(k, r) = abs(vr * v_MID(k, :)');
    end
end

figure(19);
for k = 1:4
    subplot(2, 2, k);
    plot(1:num_epochs, MI_history(k, :));
    xlabel('epoch');
    ylabel('MI');
    title(['Neuron', num2str(k)]);
end
sgtitle('Epoch vs MI');

figure(20);
for k = 1:4
    subplot(2, 2, k);
    sta_k = STA100(k, :) / norm(STA100(k, :));
    plot(0:len_v-1, sta_k, 'LineWidth', 1.5);
    hold on;
    plot(0:len_v-1, v_MID(k, :) * sign(v_MID(k, :) * sta_k'), '--', 'LineWidth', 1.5);
    xlabel('Lag (ms)');
    title(['Neuron', num2str(k)]);
    legend('STA', 'MID');
end
sgtitle('STA vs MID vector (100 ms)');

figure(21);
for k = 1:4
    subplot(2, 2, k);
    bar([MI_sta(k), MI_history(k, end), restart_MI(k, :)]);
    set(gca, 'XTickLabel', [{'STA', 'MID'}, arrayfun(@(r) ['rand', num2str(r)], 1:n_restarts, 'UniformOutput', false)]);
    ylabel('MI (bits/spike)');
    title(['Neuron', num2str(k)]);
end
sgtitle('MI at STA start, MID optimum and random restarts');

% Output nonlinearity and prediction (1 vector)
r2_sta = zeros(1, 4);
r2_mid = zeros(1, 4);
pred_mid = cell(4, 1);
actual_mid = cell(4, 1);
for k = 1:4
    lambda_tr = PSTH_rate(k, len_v:15000)';
    actual_mid{k} = PSTH_rate(k, 15000+len_v:end)';

    v = v_MID(k, :)';
    p = fit_nonlinearity(X * v, lambda_tr, num_bins);
    pred_mid{k} = sig_fun(p, X_test * v);
    r2_mid(k) = r_squared(pred_mid{k}, actual_mid{k});

    sta_k = STA100(k, :)' / norm(STA100(k, :));
    p = fit_nonlinearity(X * sta_k, lambda_tr, num_bins);
    r2_sta(k) = r_squared(sig_fun(p, X_test * sta_k), actual_mid{k});
end

figure(22);
for k = 1:4
    subplot(2, 2, k);
    scatter(actual_mid{k}, pred_mid{k}, 'filled');
    xlabel('Actual Measured Value');
    ylabel('Predicted Value');
    text(0.65, 0.9, ['r^2 = ', num2str(r2_mid(k))], 'Units', 'normalized');
    title(['Neuron', num2str(k)]);
end
sgtitle('MID prediction (1 vector)');

% Pruning of the MID vector
perf_mid = zeros(4, len_v);
v_MID_final = v_MID;
for k = 1:4
    lambda_tr = PSTH_rate(k, len_v:15000)';
    [~, idx] = sort(abs(v_MID(k, :)), 'ascend');
    for i = 0:len_v-1
        h = v_MID(k, :)';
        if i >= 1
            h(idx(1:i)) = 0;
        end
        p = fit_nonlinearity(X * h, lambda_tr, num_bins);
        perf_mid(k, len_v-i) = r_squared(sig_fun(p, X_test * h), actual_mid{k});
    end
    full_r2 = perf_mid(k, len_v);
    n_keep = len_v;
    if full_r2 > 0.01
        for np = len_v-1:-1:1
            if perf_mid(k, np) >= 0.98 * full_r2
                n_keep = np;
            else
                break;
            end
        end
    end
    v_MID_final(k, idx(1:len_v-n_keep)) = 0;
end

figure(23);
for k = 1:4
    subplot(2, 2, k);
    plot(1:len_v, perf_mid(k, :));
    xlabel('No. of parameters');
    ylabel('r^2');
    title(['Neuron', num2str(k)]);
end
sgtitle('MID: No. of parameters vs performance');

% Two vectors
mid_bins_2d = 14;
C0 = cov(X);
r2_two = zeros(1, 4);
MI_two = zeros(1, 4);
MI_two_shuffled = zeros(1, 4);
V_two = zeros(len_v, 2, 4);
for k = 1:4
    w = accumarray(spike_indices{k}(:), 1, [num_of_train_samples 1]);
    lambda_tr = PSTH_rate(k, len_v:15000)';
    v1 = v_MID(k, :)';

    mu = (X' * w) / sum(w);
    Xc = X - mu';
    dC = (Xc' * (Xc .* w)) / sum(w) - C0;
    [Ev, Dv] = eig(dC);
    ev = diag(Dv);
    if abs(ev(end)) > abs(ev(1))
        v2 = Ev(:, end);
    else
        v2 = Ev(:, 1);
    end
    v2 = v2 - v1 * (v1' * v2);
    v2 = v2 / norm(v2);

    [V, MI_two(k)] = mid2_optimize([v1, v2], X, w, mid_bins_2d);
    V_two(:, :, k) = V;
    MI_two_shuffled(k) = mid2_info(V, X, circshift(w, 7000), mid_bins_2d);

    [~, f2] = fit_nonlinearity2(X * V, lambda_tr, 12);
    r2_two(k) = r_squared(f2(X_test * V), actual_mid{k});
end

figure(24);
for k = 1:4
    subplot(2, 2, k);
    plot(0:len_v-1, V_two(:, 1, k), 'LineWidth', 1.5);
    hold on;
    plot(0:len_v-1, V_two(:, 2, k), 'LineWidth', 1.5);
    xlabel('Lag (ms)');
    title(['Neuron', num2str(k), '  r^2 = ', num2str(r2_two(k))]);
    legend('v_1', 'v_2');
end
sgtitle('Two MID vectors');

fprintf('Neuron | r2 STA(100ms) | r2 MID 1 vec | r2 MID 2 vec | MI 1D | MI 2D | MI 2D (shuffled control)\n');
for k = 1:4
    fprintf('%6d | %13.3f | %12.3f | %12.3f | %5.4f | %5.3f | %.3f\n', k, r2_sta(k), r2_mid(k), r2_two(k), MI_history(k, end), MI_two(k), MI_two_shuffled(k));
end

%% Local functions
function [MI, grad] = mid_grad(v, X, spike_rows, nb)
    x_all = X * v';
    sd = std(x_all);
    edges = linspace(-4*sd, 4*sd, nb + 1);
    centers = edges(1:end-1) + diff(edges)/2;
    bin_all = min(max(floor((x_all + 4*sd) / (8*sd/nb)) + 1, 1), nb);
    bin_spike = bin_all(spike_rows);

    P_all = accumarray(bin_all(:), 1, [nb 1])' / numel(bin_all);
    P_spike = accumarray(bin_spike(:), 1, [nb 1])' / numel(bin_spike);
    ok = P_all > 0 & P_spike > 0;
    MI = sum(P_spike(ok) .* log2(P_spike(ok) ./ P_all(ok)));

    f_x = (P_spike + 1e-8) ./ (P_all + 1e-8);
    f_x = conv(f_x, [1 2 3 2 1] / 9, 'same');
    df_x = gradient(f_x, centers);

    w_spike = df_x(bin_spike) ./ f_x(bin_spike);
    w_all = df_x(bin_all);
    term1 = (w_spike(:)' * X(spike_rows, :)) / numel(bin_spike);
    term2 = (w_all(:)' * X) / numel(bin_all);
    grad = term1 - term2;
    grad = grad - v * (grad * v');
end

function [v, MI_history] = mid_run(v0, X, spike_rows, nb, num_epochs)
    v = v0 / norm(v0);
    [MI, grad] = mid_grad(v, X, spike_rows, nb);
    lr = 1;
    MI_history = zeros(1, num_epochs);
    for epoch = 1:num_epochs
        v_new = v + lr * grad;
        v_new = v_new / norm(v_new);
        [MI_new, grad_new] = mid_grad(v_new, X, spike_rows, nb);
        if MI_new > MI
            v = v_new;
            MI = MI_new;
            grad = grad_new;
            lr = lr * 1.2;
        else
            lr = lr / 2;
        end
        MI_history(epoch) = MI;
    end
end

function [I, G] = mid2_info(V, X, w, nb)
    P = X * V;
    sd = std(P);
    bw = 7 * sd / nb;
    b = min(max(floor((P + 3.5*sd) ./ bw) + 1, 1), nb);
    k = (b(:, 2) - 1) * nb + b(:, 1);
    K = nb * nb;
    p = accumarray(k, 1, [K 1]) / size(P, 1);
    ps = accumarray(k, w, [K 1]) / sum(w);
    ok = p > 0 & ps > 0;
    I = sum(ps(ok) .* log2(ps(ok) ./ p(ok)));
    if nargout < 2
        return;
    end
    rho = zeros(K, 1);
    rho(p > 0) = ps(p > 0) ./ p(p > 0);
    rho = reshape(rho, nb, nb);
    ker = exp(-(-2:2)'.^2 / (2 * 0.8^2));
    ker = ker / sum(ker);
    rho = conv2(conv2(rho, ker, 'same'), ker', 'same');
    [d_ax2, d_ax1] = gradient(rho, bw(2), bw(1));
    S = sparse(1:numel(k), k, 1, numel(k), K);
    Sw = spdiags(w, 0, numel(w), numel(w)) * S;
    Bm = (X' * S) ./ max(sum(S, 1), 1);
    Am = (X' * Sw) ./ max(w' * S, eps);
    G = zeros(size(X, 2), 2);
    G(:, 1) = (Am - Bm) * (p .* d_ax1(:));
    G(:, 2) = (Am - Bm) * (p .* d_ax2(:));
    G = G - V .* sum(G .* V, 1);
end

function [V, I] = mid2_optimize(V0, X, w, nb)
    V = V0 ./ sqrt(sum(V0.^2, 1));
    [I, G] = mid2_info(V, X, w, nb);
    step = 1;
    for it = 1:200
        improved = false;
        for tries = 1:12
            Vn = V + step * G;
            Vn = Vn ./ sqrt(sum(Vn.^2, 1));
            In = mid2_info(Vn, X, w, nb);
            if In > I
                improved = true;
                break;
            end
            step = step / 2;
        end
        if ~improved
            break;
        end
        V = Vn;
        [I, G] = mid2_info(V, X, w, nb);
        step = step * 1.5;
    end
end

function out = sig_fun(p, x)
    out = p(1) ./ (1 + exp(-p(2) * (x - p(3)))) + p(4);
end

function p = fit_nonlinearity(y_hat, lambda, nb)
    bin_assigned = discretize(y_hat, nb);
    mean_y = nan(1, nb);
    mean_lambda = nan(1, nb);
    for b = 1:nb
        idx = (bin_assigned == b);
        if any(idx)
            mean_y(b) = mean(y_hat(idx));
            mean_lambda(b) = mean(lambda(idx));
        end
    end
    valid = ~isnan(mean_y);
    initial_guess = [max(lambda), 1/std(y_hat), median(y_hat), min(lambda)];
    opts = optimoptions('lsqcurvefit', 'Display', 'off');
    p = lsqcurvefit(@sig_fun, initial_guess, mean_y(valid), mean_lambda(valid), [], [], opts);
end

function [p, f] = fit_nonlinearity2(P, lambda, nb)
    b1 = discretize(P(:, 1), nb);
    b2 = discretize(P(:, 2), nb);
    k = (b2 - 1) * nb + b1;
    ks = unique(k);
    n_bins = numel(ks);
    mean_P = zeros(n_bins, 2);
    mean_lambda = zeros(n_bins, 1);
    cnt = zeros(n_bins, 1);
    for q = 1:n_bins
        idx = (k == ks(q));
        mean_P(q, :) = mean(P(idx, :), 1);
        mean_lambda(q) = mean(lambda(idx));
        cnt(q) = sum(idx);
    end
    prod_sig = @(p, XY) p(1) ./ ((1 + exp(-(XY(:, 1) - p(2)) / p(3))) .* (1 + exp(-(XY(:, 2) - p(4)) / p(5)))) + p(6);
    model = @(p, XY) sqrt(cnt) .* prod_sig(p, XY);
    initial_guess = [max(lambda), 0, std(P(:, 1))/2, 0, std(P(:, 2))/2, min(lambda)];
    opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 20000);
    p = lsqcurvefit(model, initial_guess, mean_P, sqrt(cnt) .* mean_lambda, [], [], opts);
    f = @(XY) prod_sig(p, XY);
end

function r2 = r_squared(a, b)
    R = corrcoef(a(:), b(:));
    r2 = R(1, 2)^2;
    if isnan(r2)
        r2 = 0;
    end
end