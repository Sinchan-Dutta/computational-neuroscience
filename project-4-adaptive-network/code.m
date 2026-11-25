clc;
close all;
clear;

%% Q1: inhomogeneous Poisson, PSTH, RMSE
dt = 1e-3;
lambda_t = 40 + 10 * rand(1, 1000);

repetitions = 640;
spike_timings = cell(repetitions, 1);
spike_trains = zeros(repetitions, 1000);

for i = 1:640
    [spike_trains(i, :), spike_timings{i}] = generate_poisson(lambda_t, dt);
end

PSTHs = zeros(7, 1000);
RMSE_errors = zeros(1, 7);
for k = 0:6
    m = 2^k * 10;
    PSTHs(k+1, :) = sum(spike_trains(1:m, :), 1) / (m * dt);
    RMSE_errors(k+1) = sqrt(sum((PSTHs(k+1, :) - lambda_t).^2) / 1000);
end

figure(1);
plot(10*2.^(0:6), RMSE_errors, 'LineWidth', 1.5);
xlabel('Number of repetitions');
ylabel('RMSE');
title('PSTH estimation error versus repetitions');
grid on;

figure(2);
plot((0:999)*dt, lambda_t, 'k', 'LineWidth', 1); hold on;
plot((0:999)*dt, PSTHs(7, :), 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('sp/s');
legend('\lambda(t)', 'PSTH (320 repetitions)');
title('PSTH versus \lambda(t)');

%% Q2: ODDBALL, verNOp, tau_ir = 5000 ms
seq = repmat('S', 1, 15);
seq(8) = 'D';
edges = (0:10:4500)*dt;
PSTH_sp = zeros(1, length(edges)-1);
PSTH_l4 = zeros(1, length(edges)-1);
plastic = 0;
w0 = [0.02; 0.02; 0.11];
[rate_S, rate_D] = generate_rates(seq, 50, 250);
for k = 1:50
    [train_s, ~] = generate_poisson(rate_S, dt);
    [train_d, ~] = generate_poisson(rate_D, dt);
    [~, spk_sp, x_s, x_d] = sp_neuron(train_s, train_d, 5000, dt);
    [~, spk_l4, ~] = l4_neuron(train_s, train_d, spk_sp, x_s, x_d, w0, plastic, 5000, dt);
    spike_times_sp = (find(spk_sp) - 1) * dt;
    spike_times_l4 = (find(spk_l4) - 1) * dt;
    [sp, ~] = histcounts(spike_times_sp, edges);
    [l4, ~] = histcounts(spike_times_l4, edges);
    PSTH_sp = PSTH_sp + sp;
    PSTH_l4 = PSTH_l4 + l4;
end

PSTH_sp = PSTH_sp / (50 * 10 * dt);
PSTH_l4 = PSTH_l4 / (50 * 10 * dt);

ctr = (edges(1:end-1) + edges(2:end))/2;
figure(3);
subplot(2,1,1); bar(ctr, PSTH_sp, 1); title('SP PSTH'); ylabel('sp/s');
subplot(2,1,2); bar(ctr, PSTH_l4, 1); title('L4 PSTH'); ylabel('sp/s'); xlabel('Time (s)');
for ax = 1:2
    subplot(2,1,ax); hold on;
    xline((0:14)*0.3, ':k');
    xline(7*0.3, 'r', 'LineWidth', 1.5);
end
sgtitle('\tau_{ir} = 5000 ms');

%% Q3: tau_ir = 1000, 3000, 10000 ms
tau_res = [1000, 3000, 10000];

for l = 1:3
    PSTH_sp = zeros(1, length(edges)-1);
    PSTH_l4 = zeros(1, length(edges)-1);
    for k = 1:50
        [train_s, ~] = generate_poisson(rate_S, dt);
        [train_d, ~] = generate_poisson(rate_D, dt);
        [~, spk_sp, x_s, x_d] = sp_neuron(train_s, train_d, tau_res(l), dt);
        [~, spk_l4, ~] = l4_neuron(train_s, train_d, spk_sp, x_s, x_d, w0, plastic, tau_res(l), dt);
        spike_times_sp = (find(spk_sp) - 1) * dt;
        spike_times_l4 = (find(spk_l4) - 1) * dt;
        [sp, ~] = histcounts(spike_times_sp, edges);
        [l4, ~] = histcounts(spike_times_l4, edges);
        PSTH_sp = PSTH_sp + sp;
        PSTH_l4 = PSTH_l4 + l4;
    end

    PSTH_sp = PSTH_sp / (50 * 10 * dt);
    PSTH_l4 = PSTH_l4 / (50 * 10 * dt);

    ctr = (edges(1:end-1) + edges(2:end))/2;
    figure(3+l);
    subplot(2,1,1); bar(ctr, PSTH_sp, 1); title('SP PSTH'); ylabel('sp/s');
    subplot(2,1,2); bar(ctr, PSTH_l4, 1); title('L4 PSTH'); ylabel('sp/s'); xlabel('Time (s)');
    for ax = 1:2
        subplot(2,1,ax); hold on;
        xline((0:14)*0.3, ':k');
        xline(7*0.3, 'r', 'LineWidth', 1.5);
    end
    sgtitle(sprintf('\\tau_{ir} = %d ms', tau_res(l)));
end

%% Q4: 60 min, p(D) = 0.1, with plasticity (chunked)
dt = 1e-3;
n_stim = 12000;
chunk_stim = 100;
n_chunks = n_stim/chunk_stim;
steps_chunk = chunk_stim*300;
ds = 10;
tau_ir = 5000;
pD = 0.1;
w0 = [0.02; 0.02; 0.11];
seq = repmat('S', 1, n_stim);
seq(rand(1, n_stim) < pD) = 'D';
n_steps = n_stim*300;
w_hist = zeros(3, n_steps/ds);
spk_sp_all = false(1, n_steps);
spk_l4_all = false(1, n_steps);
train_s_all = false(1, n_steps);
train_d_all = false(1, n_steps);
st_sp = [];
st_l4 = [];
for c = 1:n_chunks
    idx = (c-1)*chunk_stim + 1 : c*chunk_stim;
    t_idx = (c-1)*steps_chunk + (1:steps_chunk);
    [rate_S, rate_D] = generate_rates(seq(idx), 50, 250);
    train_s = generate_poisson(rate_S, dt);
    train_d = generate_poisson(rate_D, dt);
    [~, spk_sp, x_s, x_d, st_sp] = sp_neuron(train_s, train_d, tau_ir, dt, st_sp);
    [~, spk_l4, w, st_l4] = l4_neuron(train_s, train_d, spk_sp, x_s, x_d, w0, 1, tau_ir, dt, st_l4);
    w_hist(:, (c-1)*steps_chunk/ds + (1:steps_chunk/ds)) = w(:, ds:ds:end);
    spk_sp_all(t_idx) = spk_sp;
    spk_l4_all(t_idx) = spk_l4;
    train_s_all(t_idx) = train_s;
    train_d_all(t_idx) = train_d;
    save('q4_state.mat', 'st_sp', 'st_l4', 'c', 'seq');
end
save('q4_results.mat', 'seq', 'w_hist', 'spk_sp_all', 'spk_l4_all', 'train_s_all', 'train_d_all', 'ds', 'dt', '-v7.3');

t_min = (1:size(w_hist, 2))*ds*dt/60;
figure(7);
hw = plot(t_min, w_hist', 'LineWidth', 1.2);
legend(hw, {'S \rightarrow L4','D \rightarrow L4','SP \rightarrow L4'}, 'Location','east');
xlabel('Time (min)');
ylabel('Synaptic weight');

%% Q4: verNOp ODDBALL responses at time points detected from the weight curves
S4 = load('q4_results.mat', 'w_hist', 'ds', 'dt');
tt = (1:size(S4.w_hist, 2))*S4.ds*S4.dt/60;
w_s = movmean(S4.w_hist, 1000, 2);
first_t = @(mask) tt(find(mask, 1));
t_rise     = first_t(max(w_s(1:2, :), [], 1) > 0.05);
t_sp_fall  = first_t(w_s(3, :) < 0.10 & tt > 1);
t_D_sat    = first_t(w_s(2, :) > 0.38);
t_S_sat    = first_t(w_s(1, :) > 0.38);
t_sp_floor = first_t(w_s(3, :) < 0.001);
t_pts = unique(round([0.5, t_rise, t_sp_fall, t_D_sat, t_S_sat, t_sp_floor, t_sp_floor + 5], 2));
if numel(t_pts) < 5
    warning('Fewer than 5 time points detected; set t_pts by hand.');
end
figure(7); hold on;
h = xline(t_pts, ':k', compose('%.2f min', t_pts));
va = {'top','middle','bottom'};
for k = 1:numel(h)
    h(k).LabelVerticalAlignment = va{mod(k,3)+1};
end
legend(hw, {'S \rightarrow L4','D \rightarrow L4','SP \rightarrow L4'}, 'Location','east');
tau_ir = 5000;
plastic = 0;
edges = (0:10:4500)*dt;
ctr = (edges(1:end-1) + edges(2:end))/2;
seq = repmat('S', 1, 15);
seq(8) = 'D';
[rate_S, rate_D] = generate_rates(seq, 50, 250);

resp_sp = zeros(numel(t_pts), 15);
resp_l4 = resp_sp;
for p = 1:numel(t_pts)
    i_w = min(max(round(t_pts(p)*60/(S4.ds*S4.dt)), 1), size(S4.w_hist, 2));
    w_p = S4.w_hist(:, i_w);
    PSTH_sp = zeros(1, length(edges)-1);
    PSTH_l4 = zeros(1, length(edges)-1);
    for k = 1:50
        [train_s, ~] = generate_poisson(rate_S, dt);
        [train_d, ~] = generate_poisson(rate_D, dt);
        [~, spk_sp, x_s, x_d] = sp_neuron(train_s, train_d, tau_ir, dt);
        [~, spk_l4, ~] = l4_neuron(train_s, train_d, spk_sp, x_s, x_d, w_p, plastic, tau_ir, dt);
        [sp, ~] = histcounts((find(spk_sp) - 1) * dt, edges);
        [l4, ~] = histcounts((find(spk_l4) - 1) * dt, edges);
        PSTH_sp = PSTH_sp + sp;
        PSTH_l4 = PSTH_l4 + l4;
    end
    PSTH_sp = PSTH_sp / (50 * 10 * dt);
    PSTH_l4 = PSTH_l4 / (50 * 10 * dt);

    on_t = (0:14)*0.3;
    for j = 1:15
        idx = ctr >= on_t(j) & ctr < on_t(j) + 0.1;
        resp_sp(p, j) = sum(PSTH_sp(idx))*0.01;
        resp_l4(p, j) = sum(PSTH_l4(idx))*0.01;
    end

    figure(10+p);
    subplot(2,1,1); bar(ctr, PSTH_sp, 1); title('SP PSTH'); ylabel('sp/s');
    subplot(2,1,2); bar(ctr, PSTH_l4, 1); title('L4 PSTH'); ylabel('sp/s'); xlabel('Time (s)');
    for ax = 1:2
        subplot(2,1,ax); hold on;
        xline((0:14)*0.3, ':k');
        xline(7*0.3, 'r', 'LineWidth', 1.5);
    end
    sgtitle(sprintf('t = %.2f min, w = [%.3f %.3f %.3f]', t_pts(p), w_p(1), w_p(2), w_p(3)));
end

% Summary: spikes per trial (0-100 ms after onset) vs stimulus number
figure(30);
for a = 1:2
    subplot(1,2,a); hold on;
    if a == 1, R = resp_sp; ttl = 'SP'; else, R = resp_l4; ttl = 'L4'; end
    hr = plot(1:15, R', '-o'); xline(8, 'r');
    xlabel('Stimulus number'); ylabel('Spikes per trial (0-100 ms)'); title(ttl);
    legend(hr, compose('%.2f min', t_pts), 'Location', 'best');
end

%% Q5: 60 min, p(D) = 0.5, with plasticity (chunked)
dt = 1e-3;
n_stim = 12000;
chunk_stim = 100;
n_chunks = n_stim/chunk_stim;
steps_chunk = chunk_stim*300;
ds = 10;
tau_ir = 5000;
pD = 0.5;
w0 = [0.02; 0.02; 0.11];
seq = repmat('S', 1, n_stim);
seq(rand(1, n_stim) < pD) = 'D';
n_steps = n_stim*300;
w_hist = zeros(3, n_steps/ds);
spk_sp_all = false(1, n_steps);
spk_l4_all = false(1, n_steps);
train_s_all = false(1, n_steps);
train_d_all = false(1, n_steps);
st_sp = [];
st_l4 = [];
for c = 1:n_chunks
    idx = (c-1)*chunk_stim + 1 : c*chunk_stim;
    t_idx = (c-1)*steps_chunk + (1:steps_chunk);
    [rate_S, rate_D] = generate_rates(seq(idx), 50, 250);
    train_s = generate_poisson(rate_S, dt);
    train_d = generate_poisson(rate_D, dt);
    [~, spk_sp, x_s, x_d, st_sp] = sp_neuron(train_s, train_d, tau_ir, dt, st_sp);
    [~, spk_l4, w, st_l4] = l4_neuron(train_s, train_d, spk_sp, x_s, x_d, w0, 1, tau_ir, dt, st_l4);
    w_hist(:, (c-1)*steps_chunk/ds + (1:steps_chunk/ds)) = w(:, ds:ds:end);
    spk_sp_all(t_idx) = spk_sp;
    spk_l4_all(t_idx) = spk_l4;
    train_s_all(t_idx) = train_s;
    train_d_all(t_idx) = train_d;
    save('q5_state.mat', 'st_sp', 'st_l4', 'c', 'seq');
end
save('q5_results.mat', 'seq', 'w_hist', 'spk_sp_all', 'spk_l4_all', 'train_s_all', 'train_d_all', 'ds', 'dt', '-v7.3');

t_min = (1:size(w_hist, 2))*ds*dt/60;
figure(8);
plot(t_min, w_hist', 'LineWidth', 1.2);
legend('S\rightarrow L4', 'D\rightarrow L4', 'SP\rightarrow L4');
xlabel('Time (min)');
ylabel('Synaptic weight');

%% Functions
function [spike_train, spike_timing] = generate_poisson(rate, dt)
    p_t = 1 - exp(-rate * dt);
    spike_train = rand(size(p_t)) < p_t;
    spike_timing = (find(spike_train) - 1) * dt;
end

function [rate_S, rate_D] = generate_rates(seq, duration, gap)
    rate_S = zeros(1, length(seq) * (duration + gap));
    rate_D = zeros(1, length(seq) * (duration + gap));
    for k = 1:length(seq)
        st = (k-1)*(duration + gap) + 1;
        if seq(k) == 'S'
            rate_S(st : st + duration - 1) = 10;
            rate_D(st : st + duration - 1) = 2.5;
        else
            rate_S(st : st + duration - 1) = 2.5;
            rate_D(st : st + duration - 1) = 10;
        end
        rate_S(st + duration : st + duration + gap - 1) = 0.5;
        rate_D(st + duration : st + duration + gap - 1) = 0.5;
    end
end

function [V_out, spk, x_s, x_d, st] = sp_neuron(train_s, train_d, tau_ir, dt, st)
    N = length(train_s);
    tau_syn = 10e-3;
    tau_ref = 2e-3;
    ws = 0.2;
    wd = 0.2;
    if nargin < 5 || isempty(st)
        st.V = 0;
        st.g_s = 0;
        st.g_d = 0;
        st.x_s = [1; 0; 0];
        st.x_d = [1; 0; 0];
        st.refr = zeros(1, 21);
    end
    V = zeros(1, N + 1);
    g_s = zeros(1, N + 1);
    g_d = zeros(1, N + 1);
    x_s = zeros(3, N + 1);
    x_d = zeros(3, N + 1);
    V(1) = st.V;
    g_s(1) = st.g_s;
    g_d(1) = st.g_d;
    x_s(:, 1) = st.x_s;
    x_d(:, 1) = st.x_d;
    refr = zeros(1, N + 21);
    refr(1:21) = st.refr;
    spk = zeros(1, N);
    refract_amount = 5*exp(-(0:20)*dt/tau_ref);
    for k = 2:N + 1
        M_s = train_s(k-1);
        M_d = train_d(k-1);
        x_s(:, k) = rk4(x_s(:, k-1), M_s, 0.9*dt, tau_ir*dt, 10*dt, dt);
        x_d(:, k) = rk4(x_d(:, k-1), M_d, 0.9*dt, tau_ir*dt, 10*dt, dt);
        g_s(k) = g_s(k-1)*exp(-dt/tau_syn) + M_s;
        g_d(k) = g_d(k-1)*exp(-dt/tau_syn) + M_d;
        V(k) = 0.9*V(k-1) + g_s(k-1)*ws*x_s(2, k-1) + g_d(k-1)*wd*x_d(2, k-1);
        if V(k) - refr(k) >= 0.05
            spk(k-1) = 1;
            refr(k:k+20) = refr(k:k+20) + refract_amount;
        end
    end
    st.V = V(N+1);
    st.g_s = g_s(N+1);
    st.g_d = g_d(N+1);
    st.x_s = x_s(:, N+1);
    st.x_d = x_d(:, N+1);
    st.refr = refr(N+1:N+21);
    V_out = V(2:end) - refr(2:N+1);
end

function [V_l4, spk_l4, w, st] = l4_neuron(train_s, train_d, spk_sp, x_s, x_d, w0, plastic, tau_ir, dt, st)
    N = length(train_s);
    tau_syn = 10e-3;
    tau_ref = 2e-3;
    a_ltp = 0.015;
    tau_ltp = 13e-3;
    a_ltd = 0.021;
    tau_ltd = 20e-3;
    w_min = [0.0001; 0.0001; 0.0001];
    w_max = [0.4; 0.4; 0.11];
    if nargin < 10 || isempty(st)
        st.V = 0;
        st.g_s = 0;
        st.g_d = 0;
        st.g_sp = 0;
        st.x_sp = [1; 0; 0];
        st.refr = zeros(1, 21);
        st.w = w0(:);
        st.pre_tr = zeros(3, 1);
        st.post_tr = 0;
    end
    V = zeros(1, N + 1);
    g_s = zeros(1, N + 1);
    g_d = zeros(1, N + 1);
    g_sp = zeros(1, N + 1);
    x_sp = zeros(3, N + 1);
    w = zeros(3, N + 1);
    V(1) = st.V;
    g_s(1) = st.g_s;
    g_d(1) = st.g_d;
    g_sp(1) = st.g_sp;
    x_sp(:, 1) = st.x_sp;
    w(:, 1) = st.w;
    pre_tr = st.pre_tr;
    post_tr = st.post_tr;
    refr = zeros(1, N + 21);
    refr(1:21) = st.refr;
    spk_l4 = zeros(1, N);
    refract_amount = 5*exp(-(0:20)*dt/tau_ref);
    pre = [train_s; train_d; spk_sp];
    for k = 2:N + 1
        M_s = train_s(k-1);
        M_d = train_d(k-1);
        M_sp = spk_sp(k-1);
        x_sp(:, k) = rk4(x_sp(:, k-1), M_sp, 0.9*dt, tau_ir*dt, 27*dt, dt);
        g_s(k) = g_s(k-1)*exp(-dt/tau_syn) + M_s;
        g_d(k) = g_d(k-1)*exp(-dt/tau_syn) + M_d;
        g_sp(k) = g_sp(k-1)*exp(-dt/tau_syn) + M_sp;
        V(k) = 0.9*V(k-1) + g_s(k-1)*w(1, k-1)*x_s(2, k-1) ...
                          + g_d(k-1)*w(2, k-1)*x_d(2, k-1) ...
                          + g_sp(k-1)*w(3, k-1)*x_sp(2, k-1);
        if V(k) - refr(k) >= 0.05
            spk_l4(k-1) = 1;
            refr(k:k+20) = refr(k:k+20) + refract_amount;
        end
        w(:, k) = w(:, k-1);
        if plastic
            post = spk_l4(k-1);
            pre_tr = pre_tr*exp(-dt/tau_ltp);
            post_tr = post_tr*exp(-dt/tau_ltd);
            dw = a_ltp*post*pre_tr - a_ltd*pre(:, k-1)*post_tr;
            w(:, k) = min(max(w(:, k-1).*(1 + dw), w_min), w_max);
            pre_tr = pre_tr + pre(:, k-1);
            post_tr = post_tr + post;
        end
    end
    st.V = V(N+1);
    st.g_s = g_s(N+1);
    st.g_d = g_d(N+1);
    st.g_sp = g_sp(N+1);
    st.x_sp = x_sp(:, N+1);
    st.refr = refr(N+1:N+21);
    st.w = w(:, N+1);
    st.pre_tr = pre_tr;
    st.post_tr = post_tr;
    V_l4 = V(2:end) - refr(2:N+1);
    w = w(:, 2:end);
end

function next_step = rk4(x, M, tau_re, tau_ir, tau_ei, dt)
    D = [-M/tau_re, 0, 1/tau_ir;
           M/tau_re, -1/tau_ei, 0;
           0, 1/tau_ei, -1/tau_ir];
    k1 = D*x;
    k2 = D * (x + 0.5 * dt * k1);
    k3 = D * (x + 0.5 * dt * k2);
    k4 = D * (x + dt * k3);
    next_step = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end