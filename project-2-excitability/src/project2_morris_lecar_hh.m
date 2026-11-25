%% ========================================================================
%  Computational Neuroscience (EC60007) — Project II, Part I
%  Morris-Lecar Model: Phase-plane analysis, equilibria, action
%  potentials, threshold, bistability, firing rate, and manifolds.
%  =======================================================================
clc;
clear;
close all;
warning('off', 'all');

%% ===================== Q2: Model setup & nullclines =====================
% Parameters as in Rinzel & Ermentrout (R&E), Chapter 7
gCa = 4.4;    % Max Ca2+ conductance (mS/cm^2)
gK  = 8.0;    % Max K+ conductance   (mS/cm^2)
gL  = 2.0;    % Leak conductance     (mS/cm^2)
VCa = 120;    % Ca2+ reversal potential (mV)
VK  = -84;    % K+ reversal potential   (mV)
VL  = -60;    % Leak reversal potential (mV)
C   = 20;     % Membrane capacitance (uF/cm^2)
phi = 0.02;   % Recovery time-scale factor
V1  = -1.2;
V2  = 18;
V3  = 2;
V4  = 30;
I   = 0;      % Applied current (uA/cm^2)

% Steady-state / time-constant functions
m_inf  = @(V) 0.5*(1 + tanh((V-V1)/V2));
w_inf  = @(V) 0.5*(1 + tanh((V-V3)/V4));
tau_w  = @(V) 1./cosh((V-V3)/(2*V4));

% Governing ODEs. State vector X = [w; V]
dw_dt = @(w, V, phi) phi*(w_inf(V)-w)./tau_w(V);
dV_dt = @(w, V, I)   (-gCa*m_inf(V).*(V-VCa) - gK*w.*(V-VK) - gL*(V-VL) + I)/C;
dX_dt = @(X, phi, I) [dw_dt(X(1),X(2),phi); dV_dt(X(1),X(2),I)];

options_fsolve = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

% --- Nullclines ---
V_range = [-80, 60];
v_grid  = linspace(V_range(1), V_range(2), 250);
w_nullcline = w_inf(v_grid);
V_nullcline = (gCa*m_inf(v_grid).*(v_grid-VCa) + gL*(v_grid-VL) - I) ./ (-gK*(v_grid-VK));

figure(1);
plot(v_grid, V_nullcline, 'r', v_grid, w_nullcline, 'b', 'LineWidth', 1.5);
hold on;

%% ===================== Q2: Equilibrium point (two methods) ==============
% Method 1: solve nullcline intersection (1-D root find on V)
w_minus_V = @(v) w_inf(v) - (gCa*m_inf(v).*(v-VCa) + gL*(v-VL)) ./ (-gK*(v-VK));
eq_V = fzero(w_minus_V, 10);
eq_w = w_inf(eq_V);
plot(eq_V, eq_w, 'd', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 5);

% Method 2: solve dX/dt = 0 directly for both states simultaneously
X0 = [0.5, -40];
F = @(X) dX_dt(X, phi, I);
eq_X = fsolve(F, X0, options_fsolve);
plot(eq_X(2), eq_X(1), 's', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

xlabel('V'); ylabel('w');
title('Nullclines of the Morris-Lecar Model');
legend('V-nullcline', 'w-nullcline', sprintf('Eq Point (%.5f, %.5f)', eq_V, eq_w));

% --- Vector field (quiver) ---
w_vec = linspace(-0.1, 1.1, 25);
V_vec = linspace(-80, 60, 30);
[w_grid, V_grid] = meshgrid(w_vec, V_vec);

dw = dw_dt(w_grid, V_grid, phi) * 100;   % scale w so axes are comparable (V ~ O(100), w ~ O(1))
dV = dV_dt(w_grid, V_grid, I);

mag = sqrt(dw.^2 + dV.^2);
mag(mag == 0) = 1;                       % avoid divide-by-zero
dw_n = dw ./ mag * 0.01;
dV_n = dV ./ mag;

L = 1.5;  % arrow length scale
q = quiver(V_grid, w_grid, L*dV_n, L*dw_n, ...
    'Color', [0.2 0.55 0.2], 'LineWidth', 1.4, 'MaxHeadSize', 0.0004, 'AutoScale', 'off');
q.Annotation.LegendInformation.IconDisplayStyle = 'off';  % keep quiver out of legend
xlim([-80 60]); ylim([-0.1 1]);

%% ===================== Q3: Jacobian & eigenvalues at equilibrium ========
fprintf('Equilibrium point is w = %f, V = %f\n', eq_w, eq_V);

h1 = 1e-6; h2 = 1e-6;   % finite-difference steps for V and w respectively
k = dw_dt(eq_w, eq_V, phi);
l = dV_dt(eq_w, eq_V, I);
J = [(dw_dt(eq_w+h2, eq_V, phi)-k)/h2, (dw_dt(eq_w, eq_V+h1, phi)-k)/h1;
     (dV_dt(eq_w+h2, eq_V, I)-l)/h2,   (dV_dt(eq_w, eq_V+h1, I)-l)/h1];

[~, eigvals] = eig(J);
fprintf('Eigenvalues at equilibrium: %f+i%f, %f+i%f\n', ...
    real(eigvals(1,1)), imag(eigvals(1,1)), real(eigvals(2,2)), imag(eigvals(2,2)));
%% ===================== Q5: Action potential (phi = 0.02) ================
% Minimum current for an AP: phi=0.04 -> ~70.857 uA, phi=0.02 -> ~61.43 uA, phi=0.01 -> ~54.791 uA
phi = 0.02;
X0 = eq_X;
tspan = [0 500];

figure(2);
plot(v_grid, V_nullcline, 'r', 'DisplayName', 'V nullcline'); hold on;
plot(v_grid, w_nullcline, 'b', 'DisplayName', 'w nullcline'); hold on;

for I = 61.43:0.001:61.431
    [t, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan, X0);

    figure(2);
    plot(X(:,2), X(:,1), 'LineWidth', 2, 'DisplayName', sprintf('I = %.3f, phi = %.2f', I, phi));

    figure(3);
    plot(t, X(:,2), 'LineWidth', 2, 'DisplayName', sprintf('I = %.3f, phi = %.2f', I, phi));
    hold on;
end

figure(2);
plot(eq_X(2), eq_X(1), 'd', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', ...
    'MarkerSize', 5, 'DisplayName', sprintf('Starting point (%.3f, %.3f)', eq_X(2), eq_X(1)));
ylabel('w'); xlabel('Membrane Potential (V)');
title(sprintf('Membrane Potential vs w (phi = %.2f)', phi));
grid on; legend();
xlim([-80 60]); ylim([-0.1 0.5]);

figure(3);
xlabel('time'); ylabel('Membrane Potential (V)');
title(sprintf('Time vs Membrane Potential (phi = %.2f)', phi));
grid on; legend();

%% ===================== Q5: Comparison across phi values =================
X0 = eq_X;
tspan = [0 500];

figure(4);
plot(v_grid, V_nullcline, 'r', 'DisplayName', 'V nullcline'); hold on;
plot(v_grid, w_nullcline, 'b', 'DisplayName', 'w nullcline'); hold on;

for I = 80
    for phi = [0.04, 0.02, 0.01]
        [t, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan, X0);

        figure(4);
        plot(X(:,2), X(:,1), 'LineWidth', 2, 'DisplayName', sprintf('I = %.1f, phi = %.2f', I, phi));
        hold on;

        figure(5);
        plot(t, X(:,2), 'LineWidth', 2, 'DisplayName', sprintf('I = %.1f, phi = %.2f', I, phi));
        hold on;
    end

    % Equilibrium point for the last phi used above (phi = 0.01), for reference
    X1 = [0.1, -30];
    F = @(X) dX_dt(X, phi, I);
    eq_X1 = fsolve(F, X1, options_fsolve);
end

figure(4);
plot(eq_X(2), eq_X(1), 'd', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', ...
    'MarkerSize', 5, 'DisplayName', sprintf('Starting point (%.3f, %.3f)', eq_X(2), eq_X(1)));
plot(eq_X1(2), eq_X1(1), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerSize', 5, 'DisplayName', sprintf('New Eq point (%.3f, %.3f)', eq_X1(2), eq_X1(1)));
ylabel('w'); xlabel('Membrane Potential (V)');
title('Membrane Potential vs w'); grid on; legend();
xlim([-80 60]); ylim([-0.1 0.5]);

figure(5);
xlabel('time'); ylabel('Membrane Potential (V)');
title('Time vs Membrane Potential'); grid on; legend();

%% ===================== Q6: Depolarizing current pulses (Iext = 0) =======
phi = 0.02;
I = 0;
tspan = [0 500];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15);

% Initial V-depolarizations to test (chosen to bracket the threshold precisely)
dV = [45, 45.883521, 45.883522, 45.883523, 45.88352385, 45.883524, 46, 50];

figure(6);
plot(v_grid, V_nullcline, 'r', 'DisplayName', 'V nullcline', 'LineWidth', 1.5); hold on;
plot(v_grid, w_nullcline, 'b', 'DisplayName', 'w nullcline', 'LineWidth', 1.5); hold on;

for k = 1:length(dV)
    X0 = eq_X + [0, dV(k)];
    [t, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan, X0, options);

    figure(6);
    h_line = plot(X(:,2), X(:,1), 'LineWidth', 2, 'DisplayName', sprintf('V_i = %.8f', X0(2)));
    col = h_line.Color;
    plot(X0(2), X0(1), 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
        'DisplayName', sprintf('Starting point (%.3f, %.3f)', X0(2), X0(1)));

    figure(7);
    plot(t, X(:,2), 'LineWidth', 2, 'DisplayName', sprintf('V_i = %.3f, phi = %.2f', X0(2), phi));
    hold on;
end

figure(6);
plot(eq_X(2), eq_X(1), 'd', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', ...
    'MarkerSize', 5, 'DisplayName', sprintf('Eq point (%.3f, %.3f)', eq_X(2), eq_X(1)));
ylabel('w'); xlabel('Membrane Potential (V)');
title(sprintf('Membrane Potential vs w (phi = %.2f)', phi));
grid on; legend();
xlim([-80 60]); ylim([-0.1 0.5]);

figure(7);
xlabel('time'); ylabel('Membrane Potential (V)');
title(sprintf('Time vs Membrane Potential (phi = %.2f)', phi));
grid on; legend();

%% ===================== Q6: Pinpointing the threshold ====================
dV = [45.883522:1e-8:45.88352385, ...
      45.883523856:1e-10:45.883523865, ...
      45.883523865:1e-8:45.883525];

n_pts  = length(dV);
values = zeros(1, n_pts);   % peak V reached
times  = zeros(1, n_pts);   % time at which peak V occurs

for k = 1:n_pts
    X0 = eq_X + [0, dV(k)];
    [t, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan, X0, options);
    [max_V, id] = max(X(:,2));
    times(k)  = t(id);
    values(k) = max_V;
end

figure(8);
plot(dV + eq_X(2), values);
ylabel('V_{max}'); xlabel('V_{initial}');
title('max V vs V_{initial}');

figure(9);
plot(dV + eq_X(2), times);
ylabel('time taken to reach maximum voltage'); xlabel('V_{initial}');
title('time vs V_{initial}');

%% ===================== Q7: Response to Iext = 86 uA/cm^2 ================
I = 86;
phi = 0.02;
X0 = eq_X;   % IC 1: equilibrium point for Iext = 0

% IC 2: equilibrium point for Iext = 86
X2 = [0.1, -30];
F = @(X) dX_dt(X, phi, I);
eq_X2 = fsolve(F, X2, options_fsolve);

% Nullclines recomputed for Iext = 86
v_grid = linspace(-80, 60, 250);
w_nullcline_86 = w_inf(v_grid);
V_nullcline_86 = (gCa*m_inf(v_grid).*(v_grid-VCa) + gL*(v_grid-VL) - I) ./ (-gK*(v_grid-VK));

figure(10);
plot(v_grid, V_nullcline_86, 'r', 'DisplayName', 'V nullcline', 'LineWidth', 1.5); hold on;
plot(v_grid, w_nullcline_86, 'b', 'DisplayName', 'w nullcline', 'LineWidth', 1.5); hold on;

[~, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan, X0, options);
h1 = plot(X(:,2), X(:,1), 'LineWidth', 2, 'DisplayName', sprintf('I_{ext} = %.1f uA', I));

[~, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan, eq_X2, options);
plot(X(:,2), X(:,1), 'LineWidth', 2, 'DisplayName', sprintf('Initial point (%.3f, %.3f)', eq_X2(2), eq_X2(1)));

% IC 3: off-equilibrium starting point
X0_3 = [0.17; -27.9];
[~, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan, X0_3, options);
h2 = plot(X(:,2), X(:,1), 'LineWidth', 2, 'DisplayName', sprintf('Initial point (%.3f, %.3f)', X0_3(2), X0_3(1)));

plot(eq_X(2), eq_X(1), 'o', 'MarkerFaceColor', h1.Color, 'MarkerEdgeColor', h1.Color, ...
    'MarkerSize', 5, 'DisplayName', sprintf('Initial point (%.3f, %.3f)', eq_X(2), eq_X(1)));
plot(eq_X2(2), eq_X2(1), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerSize', 5, 'DisplayName', sprintf('Starting point (%.3f, %.3f)', eq_X2(2), eq_X2(1)));
plot(X0_3(2), X0_3(1), 'o', 'MarkerFaceColor', h2.Color, 'MarkerEdgeColor', h2.Color, ...
    'MarkerSize', 5, 'DisplayName', sprintf('Starting point (%.3f, %.3f)', X0_3(2), X0_3(1)));

grid on; ylabel('w'); xlabel('Membrane Potential (V)');
title(sprintf('Response to I_{ext} = %.1f uA', I));
legend();
xlim([-80 60]); ylim([-0.1 0.5]);

% --- Characterize the Iext = 86 equilibrium point ---
h1 = 1e-5; h2 = 1e-7;
eq_w = eq_X2(1);
eq_V = eq_X2(2);
k = dw_dt(eq_w, eq_V, phi);
l = dV_dt(eq_w, eq_V, I);
J = [(dw_dt(eq_w+h2, eq_V, phi)-k)/h2, (dw_dt(eq_w, eq_V+h1, phi)-k)/h1;
     (dV_dt(eq_w+h2, eq_V, I)-l)/h2,   (dV_dt(eq_w, eq_V+h1, I)-l)/h1];

[~, eigvals] = eig(J);
fprintf('Eigenvalues at equilibrium point (%.3f, %.3f): %f  %f\n', eq_X2, diag(eigvals));

%% ===================== Q8: Unstable Periodic Orbit (backward time) ======
figure(11);
plot(v_grid, V_nullcline_86, 'r', 'DisplayName', 'V nullcline', 'LineWidth', 1.5); hold on;
plot(v_grid, w_nullcline_86, 'b', 'DisplayName', 'w nullcline', 'LineWidth', 1.5); hold on;

X0_UPO = [0.17; -27.9];

% Integrating -dX/dt forward = integrating dX/dt backward in time
[~, X_UPO] = ode45(@(t, X) -dX_dt(X, phi, I), [0, 1000], X0_UPO, options);
plot(X_UPO(:,2), X_UPO(:,1), 'LineWidth', 2, 'DisplayName', 'UPO Trajectory');

plot(eq_X2(2), eq_X2(1), 'd', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
    'MarkerSize', 5, 'DisplayName', sprintf('Eq point (%.3f, %.3f)', eq_X2(2), eq_X2(1)));
plot(X0_UPO(2), X0_UPO(1), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', ...
    'MarkerSize', 5, 'DisplayName', sprintf('Starting point (%.3f, %.3f), backward trajectory', X0_UPO(2), X0_UPO(1)));

grid on; ylabel('w'); xlabel('Membrane Potential (V)');
title(sprintf('Response to I_{ext} = %.1f uA', I));
legend();
xlim([-80 60]); ylim([-0.1 0.5]);

%% ===================== Q8: UPO acts as a true threshold =================
UPO = X_UPO(16000:end, :);

figure(12);
plot(v_grid, V_nullcline_86, 'r', 'DisplayName', 'V nullcline', 'LineWidth', 1.5); hold on;
plot(v_grid, w_nullcline_86, 'b', 'DisplayName', 'w nullcline', 'LineWidth', 1.5); hold on;
plot(UPO(:,2), UPO(:,1), 'k', 'DisplayName', 'UPO', 'LineWidth', 1.5); hold on;

% Perturb a point on the UPO by an infinitesimal epsilon, one step out, one step in
epsilon = [0; 1e-9];
X3  = UPO(3700, :)';
X31 = X3 + epsilon;   % just outside the UPO
X32 = X3 - epsilon;   % just inside the UPO
tspan_upo = [0 800];

[~, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan_upo, X31, options);
h1 = plot(X(:,2), X(:,1), 'LineWidth', 1.5, 'DisplayName', 'Point Outside UPO orbit');
plot(X31(2), X31(1), 'd', 'MarkerFaceColor', h1.Color, 'MarkerEdgeColor', h1.Color, ...
    'MarkerSize', 5, 'DisplayName', sprintf('Starting point (%.6f, %.6f)', X31(2), X31(1)));

[~, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan_upo, X32, options);
h2 = plot(X(:,2), X(:,1), 'LineWidth', 1.5, 'DisplayName', 'Point Inside UPO orbit');
plot(X32(2), X32(1), 'o', 'MarkerFaceColor', h2.Color, 'MarkerEdgeColor', h2.Color, ...
    'MarkerSize', 3, 'DisplayName', sprintf('Starting point (%.6f, %.6f)', X32(2), X32(1)));

grid on; ylabel('w'); xlabel('Membrane Potential (V)');
title(sprintf('Response to I_{ext} = %.1f uA', I));
legend();
xlim([-80 60]); ylim([-0.1 0.5]);

%% ===================== Q9: Eigenvalues at Iext = 80, 86, 90 =============
phi = 0.02;
Is = [80, 86, 90];
X0 = [0.5, -35];

for m = 1:length(Is)
    I = Is(m);
    F = @(X) dX_dt(X, phi, I);
    eq_X_m = fsolve(F, X0, options_fsolve);

    h1 = 1e-5; h2 = 1e-7;
    eq_w = eq_X_m(1);
    eq_V = eq_X_m(2);
    k = dw_dt(eq_w, eq_V, phi);
    l = dV_dt(eq_w, eq_V, I);
    J = [(dw_dt(eq_w+h2, eq_V, phi)-k)/h2, (dw_dt(eq_w, eq_V+h1, phi)-k)/h1;
         (dV_dt(eq_w+h2, eq_V, I)-l)/h2,   (dV_dt(eq_w, eq_V+h1, I)-l)/h1];

    [~, eigvals] = eig(J);
    fprintf('Eigenvalues at equilibrium for I = %.0f uA/cm^2: %f + %fi, %f + %fi\n', ...
        I, real(eigvals(1,1)), imag(eigvals(1,1)), real(eigvals(2,2)), imag(eigvals(2,2)));
end

%% ===================== Q9: Firing rate vs applied current ===============
Is = 80:100;
X0 = [0.5, -35];
firing_rates = zeros(size(Is));

for m = 1:length(Is)
    I = Is(m);
    F = @(X) dX_dt(X, phi, I);
    eq_X_m = fsolve(F, X0, options_fsolve);

    % Simulation length / steady-state cutoff grows near the bifurcation
    % point (~I=90-92) where transients take longer to settle
    if I < 90
        tspan_fr = [0 500];  R = 250;
    elseif I == 90
        tspan_fr = [0 15000]; R = 10000;
    elseif I == 91
        tspan_fr = [0 7000];  R = 5000;
    elseif I == 92
        tspan_fr = [0 5000];  R = 3000;
    else
        tspan_fr = [0 4000];  R = 2500;
    end

    [t, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan_fr, eq_X_m, options);

    steady_response = t > R;
    steady_V = X(steady_response, 2);
    steady_t = t(steady_response);

    [peaks, locs] = findpeaks(steady_V, steady_t, 'MinPeakHeight', 35);
    if length(peaks) >= 2
        T = mean(diff(locs));   % mean inter-spike interval (ms)
        f = 1000/T;             % convert to Hz
    else
        f = 0;                  % no sustained firing
    end
    firing_rates(m) = f;
end

figure(13);
plot(Is, firing_rates);
grid on;
xlabel('I_{ext} uA/cm^2'); ylabel('Firing Rate (Hz)');
title('Firing Rate for different I_{ext}');

%% ===================== Q10: New parameter set, 3 equilibria & manifolds =
gCa = 4; gK = 8.0; gL = 2.0;
VCa = 120; VK = -84; VL = -60;
C = 20; phi = 0.0667;
V1 = -1.2; V2 = 18; V3 = 12; V4 = 17.4;
I = 30;

m_inf = @(V) 0.5*(1 + tanh((V-V1)/V2));
w_inf = @(V) 0.5*(1 + tanh((V-V3)/V4));
tau_w = @(V) 1./cosh((V-V3)/(2*V4));
dw_dt = @(w, V, phi) phi*(w_inf(V)-w)./tau_w(V);
dV_dt = @(w, V, I)   (-gCa*m_inf(V).*(V-VCa) - gK*w.*(V-VK) - gL*(V-VL) + I)/C;
dX_dt = @(X, phi, I) [dw_dt(X(1),X(2),phi); dV_dt(X(1),X(2),I)];

% --- Nullclines ---
v_grid = linspace(-80, 60, 400);
w_nullcline = w_inf(v_grid);
V_nullcline = (gCa*m_inf(v_grid).*(v_grid-VCa) + gL*(v_grid-VL) - I) ./ (-gK*(v_grid-VK));

figure(14); clf;
plot(v_grid, V_nullcline, 'r', 'LineWidth', 1.5); hold on;
plot(v_grid, w_nullcline, 'b', 'LineWidth', 1.5);

% --- Find the three equilibrium points ---
options_fsolve = optimoptions('fsolve', 'Display', 'off');
X0_s = [0.01, -45;   % left
        0.15, -20;   % middle
        0.50,  10];  % right

eq_Xs = zeros(3, 2);
for i = 1:3
    eq_Xs(i, :) = fsolve(@(X) [dw_dt(X(1),X(2),phi); dV_dt(X(1),X(2),I)], X0_s(i,:), options_fsolve);
end
[~, idx] = sort(eq_Xs(:,2));
eq_Xs = eq_Xs(idx, :);   % sorted by V: left -> middle -> right

% --- Characterize each equilibrium (stability via eigenvalues) ---
saddle_idx = 0;
fprintf('\n--- Equilibrium Point Stability Analysis ---\n');
for m = 1:3
    eq_w = eq_Xs(m,1);
    eq_V = eq_Xs(m,2);

    h1 = 1e-5; h2 = 1e-7;
    k = dw_dt(eq_w, eq_V, phi);
    l = dV_dt(eq_w, eq_V, I);
    J = [(dw_dt(eq_w+h2, eq_V, phi)-k)/h2, (dw_dt(eq_w, eq_V+h1, phi)-k)/h1;
         (dV_dt(eq_w+h2, eq_V, I)-l)/h2,   (dV_dt(eq_w, eq_V+h1, I)-l)/h1];

    [eig_vecs, eig_vals] = eig(J);
    l1 = eig_vals(1,1); l2 = eig_vals(2,2);

    fprintf('Eq %d (V=%.2f, w=%.4f): lambda1 = %.4f + %.4fi, lambda2 = %.4f + %.4fi\n', ...
        m, eq_V, eq_w, real(l1), imag(l1), real(l2), imag(l2));

    if real(l1) * real(l2) < 0
        saddle_idx  = m;
        saddle_vecs = eig_vecs;
        saddle_vals = [l1; l2];
    end
end

plot(eq_Xs(1,2), eq_Xs(1,1), 'og', 'MarkerFaceColor', 'g', 'MarkerSize', 5);  % left  (stable node/focus)
plot(eq_Xs(2,2), eq_Xs(2,1), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5);  % middle (saddle)
plot(eq_Xs(3,2), eq_Xs(3,1), 'om', 'MarkerFaceColor', 'm', 'MarkerSize', 5);  % right (unstable)

tspan = [0 300];

% --- Stable/unstable manifolds of the saddle point ---
if saddle_idx > 0
    eps_val = 1e-5;   % small perturbation off the saddle along each eigenvector
    saddle_pt = eq_Xs(saddle_idx, :)';

    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'Events', @stopOutside);

    for v_i = 1:2
        v_dir = saddle_vecs(:, v_i);
        val   = real(saddle_vals(v_i));

        if val > 0
            % Unstable manifold: integrate forward in time
            [~, sol1] = ode45(@(t, X) dX_dt(X, phi, I), tspan, saddle_pt + eps_val*v_dir, opts);
            [~, sol2] = ode45(@(t, X) dX_dt(X, phi, I), tspan, saddle_pt - eps_val*v_dir, opts);
            plot(sol1(:,2), sol1(:,1), 'm--', 'LineWidth', 1.5);
            plot(sol2(:,2), sol2(:,1), 'm--', 'LineWidth', 1.5);
        else
            % Stable manifold: integrate backward in time
            dX_backward = @(t, X) -dX_dt(X, phi, I);
            [~, sol1] = ode45(dX_backward, tspan, saddle_pt + eps_val*v_dir, opts);
            [~, sol2] = ode45(dX_backward, tspan, saddle_pt - eps_val*v_dir, opts);
            plot(sol1(:,2), sol1(:,1), 'k-', 'LineWidth', 1.5);
            plot(sol2(:,2), sol2(:,1), 'k-', 'LineWidth', 1.5);
        end
    end
end

% --- Vector field (quiver) ---
w_vec = linspace(-0.1, 1.1, 25);
V_vec = linspace(-80, 60, 30);
[w_grid, V_grid] = meshgrid(w_vec, V_vec);

dw = dw_dt(w_grid, V_grid, phi) * 100;
dV = dV_dt(w_grid, V_grid, I);
mag = sqrt(dw.^2 + dV.^2);
mag(mag == 0) = 1;
dw_n = dw ./ mag * 0.01;
dV_n = dV ./ mag;

L = 1.5;
q = quiver(V_grid, w_grid, L*dV_n, L*dw_n, ...
    'Color', [0.2 0.55 0.2], 'LineWidth', 1.4, 'MaxHeadSize', 0.0004, 'AutoScale', 'off');
q.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlim([-80 60]); ylim([-0.05 0.8]);
xlabel('Membrane Potential V (mV)'); ylabel('Recovery Variable w');
title('Morris-Lecar Phase Plane with Manifolds (I_{ext} = 30)');
legend('V-nullcline', 'w-nullcline', 'Eq 1 (Stable Node/Focus)', 'Eq 2 (Saddle Point)', ...
    'Eq 3 (Unstable)', 'Stable Manifold', 'Stable Manifold', 'Unstable Manifold', 'Unstable Manifold');
grid on;

%% ===================== Q11: Equilibria vs Iext (30-50) ==================

Is = [30:0.5:39.5, 39.5:0.1:39.9, 39.92:0.02:39.98, 40:50];
options_fsolve = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Display', 'off');
resid_tol = 1e-8;   % max allowed ||F(X)|| to accept a solution as a true equilibrium
dedup_tol = 1e-3;   % treat two roots closer than this (in [w,V]) as the same point

X0_guesses = [0.05, -40;   % left / resting branch
              0.05, -20;   % middle / saddle branch
              0.30,   0];  % right / upper branch

figure(15); clf;

for I = Is
    F = @(X) dX_dt(X, phi, I);
    found_pts = [];   % accepted [w, V] equilibria for this I

    for i = 1:3
        [X, fval, exitflag] = fsolve(F, X0_guesses(i,:), options_fsolve);

        if exitflag > 0 && norm(fval) < resid_tol
            is_dup = false;
            for p = 1:size(found_pts,1)
                if norm(found_pts(p,:) - X) < dedup_tol
                    is_dup = true;
                    break;
                end
            end
            if ~is_dup
                found_pts = [found_pts; X]; %#ok<AGROW>
            end
        end
    end

    % --- Classify stability of each accepted equilibrium (as in Q9/Q10) ---
    for p = 1:size(found_pts,1)
        eq_w = found_pts(p,1);
        eq_V = found_pts(p,2);

        h1 = 1e-5; h2 = 1e-7;
        k = dw_dt(eq_w, eq_V, phi);
        l = dV_dt(eq_w, eq_V, I);
        J = [(dw_dt(eq_w+h2, eq_V, phi)-k)/h2, (dw_dt(eq_w, eq_V+h1, phi)-k)/h1;
             (dV_dt(eq_w+h2, eq_V, I)-l)/h2,   (dV_dt(eq_w, eq_V+h1, I)-l)/h1];
        eigvals = eig(J);
        re = real(eigvals);

        if all(re < 0)
            color = 'g';   % stable node/focus
        elseif all(re > 0)
            color = 'm';   % unstable node/focus
        else
            color = 'k';   % saddle
        end

        subplot(1, 2, 1); hold on;
        plot(I, eq_V, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'MarkerSize', 4);

        subplot(1, 2, 2); hold on;
        plot(I, eq_w, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'MarkerSize', 4);
    end
end

subplot(1, 2, 1);
xlabel('I_{ext} (\muA/cm^2)'); ylabel('Equilibrium V (mV)');
title('Equilibrium V vs I_{ext} (green=stable, black=saddle, magenta=unstable)');
grid on;

subplot(1, 2, 2);
xlabel('I_{ext} (\muA/cm^2)'); ylabel('Equilibrium w');
title('Equilibrium w vs I_{ext} (green=stable, black=saddle, magenta=unstable)');
grid on;

%% ===================== Q11: Firing rate vs Iext (30-50) =================

Is_fr = [30:39, 39.1:0.25:39.85, 39.9, 39.95, 39.97:0.01:39.99, 40:50];
firing_rates_11 = zeros(size(Is_fr));
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

X0_track = X0_guesses(2,:);   % branch-tracking guess, updated by continuation below

for idx = 1:length(Is_fr)
    I = Is_fr(idx);
    F = @(X) dX_dt(X, phi, I);

    [X_eq, fval, exitflag] = fsolve(F, X0_track, options_fsolve);
    if exitflag > 0 && norm(fval) < resid_tol
        X0 = X_eq;
        X0_track = X_eq;   % carry forward for next I
    else
        X0 = X0_track;      % equilibrium gone past the bifurcation -- start nearby
    end

    tspan_fr = [0 3000];
    R = 1000;   % discard initial transient before this time when counting spikes

    [t, X] = ode45(@(t, X) dX_dt(X, phi, I), tspan_fr, X0, options);

    steady_response = t > R;
    steady_V = X(steady_response, 2);
    steady_t = t(steady_response);

    [peaks, locs] = findpeaks(steady_V, steady_t, 'MinPeakHeight', 25);
    if length(peaks) >= 2
        T = mean(diff(locs));
        firing_rates_11(idx) = 1000/T;
    else
        firing_rates_11(idx) = 0;
    end
end

figure(16);
plot(Is_fr, firing_rates_11, '-o');
xlabel('I_{ext} (\muA/cm^2)'); ylabel('Firing Rate (Hz)');
title('Firing rate vs I_{ext}');
grid on;

%% ========================================================================
%  Computational Neuroscience (EC60007) — Project II, Part II
%  Hodgkin-Huxley Model: action potentials, threshold, equilibria,
%  reduced models, anode-break excitation.
%  =======================================================================
% Units: milliVolt, milliSecond, microAmpere/cm^2, milliSiemens/cm^2,
% microFarad/cm^2 (a consistent unit set, as established for the MLE).

%% ===================== Q12: Rate functions & model setup ================
% alpha_n and alpha_m are defined at the end, as local functions

beta_n = @(V) 0.125*exp(-(V+60)/80);
beta_m = @(V) 4*exp(-(V+60)/18);

alpha_h = @(V) 0.07*exp(-(V+60)/20);
beta_h  = @(V) 1./(exp(-(V+30)/10) + 1);

C = 1.0; GL_bar = 0.3; GK_bar = 36; GNa_bar = 120;
EK = -72; ENa = 55; v_rest = -60;

m_rest = alpha_m(v_rest)/(alpha_m(v_rest)+beta_m(v_rest));
n_rest = alpha_n(v_rest)/(alpha_n(v_rest)+beta_n(v_rest));
h_rest = alpha_h(v_rest)/(alpha_h(v_rest)+beta_h(v_rest));

%% ===================== Q13: Fit EL so resting potential is -60 mV =======
GK_rest  = GK_bar*n_rest^4;
GNa_rest = GNa_bar*m_rest^3*h_rest;
EL = v_rest + (GK_rest*(v_rest-EK) + GNa_rest*(v_rest-ENa))/GL_bar;

% Governing ODEs. State vector Y = [V; m; n; h]
dm_dt = @(m, V) alpha_m(V)*(1-m) - beta_m(V)*m;
dn_dt = @(n, V) alpha_n(V)*(1-n) - beta_n(V)*n;
dh_dt = @(h, V) alpha_h(V)*(1-h) - beta_h(V)*h;
dV_dt = @(m, n, h, V, I) (-GK_bar*n^4*(V-EK) - GNa_bar*m^3*h*(V-ENa) - GL_bar*(V-EL) + I)/C;

dY_dt = @(X, I) [ dV_dt(X(2), X(3), X(4), X(1), I);
                  dm_dt(X(2), X(1));
                  dn_dt(X(3), X(1));
                  dh_dt(X(4), X(1)) ];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15);   % stiff-safe tolerances for ode15s

% --- Action potential with Iext = 10 uA/cm^2 ---
I1 = 10;
tspan = [0 50];
Y0 = [v_rest; m_rest; n_rest; h_rest];
[t1, Y1] = ode15s(@(t, Y) dY_dt(Y, I1), tspan, Y0, options);

figure(17);
plot(t1, Y1(:,1), 'LineWidth', 1.5);
xlabel('time'); ylabel('Membrane Potential (V)');
title('Action Potentials using Hodgkin Huxley Model');
legend(sprintf('I = %.1fuA/cm^2', I1));
grid on;

% --- m, n, h over one full AP period ---
figure(18);
idx = (t1 > 20 & t1 < 35.7);   % one full AP, -60mV -> spike -> back to -60mV (found by inspection)
plot(t1(idx), Y1(idx,2), 'LineWidth', 1.5, 'DisplayName', 'm'); hold on;
plot(t1(idx), Y1(idx,3), 'LineWidth', 1.5, 'DisplayName', 'n');
plot(t1(idx), Y1(idx,4), 'LineWidth', 1.5, 'DisplayName', 'h');
xlabel('time'); ylabel('m / n / h');
legend(); grid on;
title('How m, n, h vary over a full action potential');

%% ===================== Q14: Stability at rest & threshold ===============
I = 0;
tspan = [0 100];
Y0 = [v_rest; m_rest; n_rest; h_rest];

% Jacobian via forward finite differences on the 4-D system
eps_val = 1e-7;
N = length(Y0);
J = zeros(N, N);
f0 = dY_dt(Y0, I);

for j = 1:N
    Y_perturb = Y0;
    Y_perturb(j) = Y_perturb(j) + eps_val;
    f_perturb = dY_dt(Y_perturb, I);
    J(:, j) = (f_perturb - f0) / eps_val;
end

eigen_values = eig(J);
if all(real(eigen_values) < 0)
    fprintf('The system is Stable at rest with I_ext = 0\n');
end

% --- Threshold search: sweep initial V, find max V reached ---
del_V = [linspace(-54, -53.5, 20), linspace(-53.5, -53.4925, 40), ...
         linspace(-53.4925, -53.4924, 100), linspace(-53.4924, -53.47, 40), ...
         linspace(-53.47, -53, 20)];
V_max = zeros(size(del_V));

figure(19);
for i = 1:length(del_V)
    Y0_ = [del_V(i); m_rest; n_rest; h_rest];
    [~, Y] = ode15s(@(t, Y) dY_dt(Y, 0), tspan, Y0_, options);
    V_max(i) = max(Y(:,1));
end

plot(del_V, V_max);
ylabel('V_{max}'); xlabel('V_{initial}');
title('V_{max} vs V_{initial} to find threshold');
% Threshold identified at approximately V_initial = -53.49245 mV

%% ===================== Q15: Equilibria vs Iext (8-12) ===================
Is = sort(unique([linspace(8, 12, 50), linspace(9.62, 9.88, 30)]));   % refined near the transition
options_fsolve = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Display', 'off');
resid_tol = 1e-8;
dedup_tol = 1e-3;

Y0_guess = [v_rest, m_rest, n_rest, h_rest];

figure(20); clf; hold on;
N = length(Y0_guess);
J = zeros(N, N);

for I = Is
    F = @(Y) dY_dt(Y, I);
    Y0 = fsolve(F, Y0_guess, options_fsolve);

    f0 = dY_dt(Y0, I);
    for j = 1:N
        Y_perturb = Y0;
        Y_perturb(j) = Y_perturb(j) + eps_val;
        f_perturb = dY_dt(Y_perturb, I);
        J(:, j) = (f_perturb - f0) / eps_val;
    end

    eigvals = eig(J);
    if all(real(eigvals) < 0)
        color = 'g';   % stable
    else
        color = 'm';   % unstable
    end

    plot(I, Y0(1), 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'MarkerSize', 4);
end

xlabel('I_{ext} (\muA/cm^2)'); ylabel('Equilibrium V (mV)');
title('Equilibrium points vs I_{ext} (green=stable, magenta=unstable)');
grid on;

%% ===================== Q16: Reduced (V, n) model ========================
% m activates much faster than V, n, h, so it is replaced by its
% instantaneous steady state m_inf(V) rather than integrated explicitly.
% h is fixed at its resting value here (this reduction keeps n as the
% only slow recovery variable, with h held constant rather than slaved
% to n).
m_reduced = @(V) alpha_m(V)/(alpha_m(V)+beta_m(V));
dV_dt_reduced1 = @(n, V, I) (-GK_bar*n^4*(V-EK) - GNa_bar*m_reduced(V)^3*h_rest*(V-ENa) - GL_bar*(V-EL) + I)/C;

dX_dt1 = @(X, I) [dV_dt_reduced1(X(2), X(1), I); dn_dt(X(2), X(1))];

I = 0;
X0 = [v_rest; n_rest] + [4.5; 0];   % depolarizing displacement from rest
tspan = [0 50];
[t_reduced, X_reduced] = ode15s(@(t, X) dX_dt1(X, I), tspan, X0, options);
V_reduced = X_reduced(:,1);

figure(21);
plot(t_reduced, V_reduced, 'LineWidth', 1.5);
xlabel('time'); ylabel('Membrane Potential (V)');
title('Action Potentials using Hodgkin Huxley Reduced Model (h = h_{rest})');
legend(sprintf('I = %.1f uA/cm^2, V_{initial} = %.2f', I, X0(1)));
grid on;

%% ===================== Q17: Anode-break excitation ======================
% Hyperpolarizing pulse (I = -3 uA/cm^2, 20 ms), then released to I = 0
I = -3;
tspan2 = [0 20];
Y0 = [v_rest; m_rest; n_rest; h_rest];
[t2, Y2] = ode15s(@(t, Y) dY_dt(Y, I), tspan2, Y0, options);

I = 0;
tspan3 = [20 50];
Y2_end = Y2(end, :);
[t3, Y3] = ode15s(@(t, Y) dY_dt(Y, I), tspan3, Y2_end, options);

t_total = [t2; t3(2:end)];
Y = [Y2; Y3(2:end, :)];

figure(22);
plot(t_total, Y(:,1), 'r', 'LineWidth', 1.5); hold on;
xlabel('Time (ms)'); ylabel('Membrane Potential (V)');
title('Generation of Action Potential due to Anode Break');
grid on;

%% ===================== Q18 (preliminary): (V, m) reduced systems ========
% Two reduced (V, m) systems, with n and h fixed either at their resting
% values (Case 1) or at their values at the end of the anodal pulse
% (Case 2) -- initial exploration ahead of the full Q18 phase-plane below.
I = 0;
dV_dt_reduced1 = @(m, V, I) (-GK_bar*n_rest^4*(V-EK) - GNa_bar*m^3*h_rest*(V-ENa) - GL_bar*(V-EL) + I)/C;

h_ = Y2_end(3);
n_ = Y2_end(4);
dV_dt_reduced2 = @(m, V, I) (-GK_bar*n_^4*(V-EK) - GNa_bar*m^3*h_*(V-ENa) - GL_bar*(V-EL) + I)/C;

dX_dt1 = @(X, I) [dV_dt_reduced1(X(2), X(1), I); dm_dt(X(2), X(1))];
dX_dt2 = @(X, I) [dV_dt_reduced2(X(2), X(1), I); dm_dt(X(2), X(1))];

X0_guess = [v_rest, m_rest];

F = @(X) dX_dt1(X, I);
X10 = fsolve(F, X0_guess, options_fsolve);
f10 = dX_dt1(X10, I);
J1 = zeros(2, 2);
for j = 1:2
    X_perturb = X10;
    X_perturb(j) = X_perturb(j) + eps_val;
    J1(:, j) = (dX_dt1(X_perturb, I) - f10) / eps_val;
end

F = @(X) dX_dt2(X, I);
X20 = fsolve(F, X0_guess, options_fsolve);
f20 = dX_dt2(X20, I);
J2 = zeros(2, 2);
for j = 1:2
    X_perturb = X20;
    X_perturb(j) = X_perturb(j) + eps_val;
    J2(:, j) = (dX_dt2(X_perturb, I) - f20) / eps_val;
end

tspan = [0 50];
[~, X1] = ode15s(@(t, X) dX_dt1(X, I), tspan, X10, options);
[~, X2] = ode15s(@(t, X) dX_dt2(X, I), tspan, X20, options);

figure(23);
subplot(1, 2, 1);
plot(X1(:,1), X1(:,2), 'LineWidth', 1.5, 'DisplayName', sprintf('I = %.1f uA/cm^2', I)); hold on;
plot(X10(2), X10(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5, 'DisplayName', 'Stable Point');
ylabel('m'); xlabel('Membrane Potential (V)');
title('Reduced Model 1'); legend(); grid on;

subplot(1, 2, 2);
plot(X2(:,1), X2(:,2), 'LineWidth', 1.5, 'DisplayName', sprintf('I = %.1f uA/cm^2', I)); hold on;
plot(X20(1), X20(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5, 'DisplayName', 'Stable Point');
ylabel('m'); xlabel('Membrane Potential (V)');
title('Reduced Model 2'); legend(); grid on;

%% ===================== Q18: (V, m) phase-plane, anode-break threshold ====
I_post = 0;

% Case 1: n, h fixed at resting values (-60 mV)
n_case1 = n_rest;
h_case1 = h_rest;

% Case 2: n, h fixed at their values at the end of the hyperpolarizing pulse
n_case2 = Y2_end(3);   % n de-activated
h_case2 = Y2_end(4);   % h de-inactivated

% Reduced 2-D system. State vector X = [V; m]
dVdt_1 = @(V, m) (-GK_bar*(n_case1^4)*(V-EK) - GNa_bar*(m^3)*h_case1*(V-ENa) - GL_bar*(V-EL) + I_post) / C;
dVdt_2 = @(V, m) (-GK_bar*(n_case2^4)*(V-EK) - GNa_bar*(m^3)*h_case2*(V-ENa) - GL_bar*(V-EL) + I_post) / C;
dmdt_red = @(V, m) alpha_m(V)*(1-m) - beta_m(V)*m;

dXdt_1 = @(X) [dVdt_1(X(1), X(2)); dmdt_red(X(1), X(2))];
dXdt_2 = @(X) [dVdt_2(X(1), X(2)); dmdt_red(X(1), X(2))];

options_fsolve = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'Display', 'off');
eps_val = 1e-7;

% --- Case 1 equilibria ---
X1_eq = fsolve(dXdt_1, [v_rest; m_rest], options_fsolve);
J1 = zeros(2, 2); f10 = dXdt_1(X1_eq);
for j = 1:2
    X_p = X1_eq; X_p(j) = X_p(j) + eps_val;
    J1(:, j) = (dXdt_1(X_p) - f10) / eps_val;
end
eig1 = eig(J1);

X1_eq2 = fsolve(dXdt_1, [40, 1], options_fsolve);
J12 = zeros(2, 2); f10 = dXdt_1(X1_eq2);
for j = 1:2
    X_p = X1_eq2; X_p(j) = X_p(j) + eps_val;
    J12(:, j) = (dXdt_1(X_p) - f10) / eps_val;
end
eig12 = eig(J12);

% --- Case 2 equilibria ---
X2_eq = fsolve(dXdt_2, [40; 1], options_fsolve);
J2 = zeros(2, 2); f20 = dXdt_2(X2_eq);
for j = 1:2
    X_p = X2_eq; X_p(j) = X_p(j) + eps_val;
    J2(:, j) = (dXdt_2(X_p) - f20) / eps_val;
end
eig2 = eig(J2);

X2_eq2 = fsolve(dXdt_2, [40, 1], options_fsolve);
J22 = zeros(2, 2); f20 = dXdt_2(X2_eq2);
for j = 1:2
    X_p = X2_eq2; X_p(j) = X_p(j) + eps_val;
    J22(:, j) = (dXdt_2(X_p) - f20) / eps_val;
end
eig22 = eig(J22);

fprintf('=== Q18: Equilibrium Point Characterization ===\n');
fprintf('Case 1 Eq1: V = %.2f mV, m = %.4f | Eigenvalues: %.2f, %.2f\n', X1_eq(1), X1_eq(2), eig1(1), eig1(2));
fprintf('Case 1 Eq2: V = %.2f mV, m = %.4f | Eigenvalues: %.2f, %.2f\n', X1_eq2(1), X1_eq2(2), eig12(1), eig12(2));
fprintf('Case 2 Eq1: V = %.2f mV, m = %.4f | Eigenvalues: %.2f, %.2f\n\n', X2_eq(1), X2_eq(2), eig2(1), eig2(2));
fprintf('Case 2 Eq2: V = %.2f mV, m = %.4f | Eigenvalues: %.2f, %.2f\n\n', X2_eq2(1), X2_eq2(2), eig22(1), eig22(2));

% --- Nullclines for both cases ---
V_grid = linspace(-80, 60, 300);
m_null = alpha_m(V_grid) ./ (alpha_m(V_grid) + beta_m(V_grid));   % m-nullcline (shared shape)

m_null_V1 = zeros(size(V_grid));
m_null_V2 = zeros(size(V_grid));
for k = 1:length(V_grid)
    Vk = V_grid(k);
    term1_1 = (GK_bar*(n_case1^4)*(Vk-EK) + GL_bar*(Vk-EL) - I_post) / (-GNa_bar*h_case1*(Vk-ENa));
    term1_2 = (GK_bar*(n_case2^4)*(Vk-EK) + GL_bar*(Vk-EL) - I_post) / (-GNa_bar*h_case2*(Vk-ENa));
    m_null_V1(k) = nthroot(term1_1, 3);
    m_null_V2(k) = nthroot(term1_2, 3);
end

tspan = [0 30];
[t1, X1] = ode15s(@(t, X) dXdt_1(X), tspan, X1_eq, options);
[t2, X2] = ode15s(@(t, X) dXdt_2(X), tspan, [Y2_end(1), Y2_end(2)], options);

% --- Phase planes ---
figure(24); clf;

subplot(1, 2, 1); hold on;
plot(V_grid, m_null, 'k--', 'LineWidth', 1.5, 'DisplayName', 'm-nullcline');
plot(V_grid, m_null_V1, 'r-', 'LineWidth', 1.5, 'DisplayName', 'V-nullcline');
plot(X1(:,1), X1(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'evolution');
plot(X1_eq(1), X1_eq(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5, 'DisplayName', 'Stable Rest Eq');
plot(X1_eq2(1), X1_eq2(2), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 5, 'DisplayName', 'Stable Eq 2');
xlabel('Membrane Potential V (mV)'); ylabel('m');
title('Q18 Case 1: Fixed at Rest (-60 mV)');
ylim([0 1.1]); xlim([-80 60]);
legend('Location', 'northwest'); grid on;

subplot(1, 2, 2); hold on;
plot(V_grid, m_null, 'k--', 'LineWidth', 1.5, 'DisplayName', 'm-nullcline');
plot(V_grid, m_null_V2, 'r-', 'LineWidth', 1.5, 'DisplayName', 'V-nullcline');
plot(X2(:,1), X2(:,2), 'b--', 'LineWidth', 1.5, 'DisplayName', 'evolution');
plot(X2_eq2(1), X2_eq2(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5, 'DisplayName', 'Stable Eq');
plot(Y2_end(1), Y2_end(2), 'm*', 'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'State at Anode Release');
xlabel('Membrane Potential V (mV)'); ylabel('m');
title('Q18 Case 2: Fixed at End of Anodal Stimulus');
ylim([0 1.1]); xlim([-80 60]);
legend('Location', 'northwest'); grid on;

% --- Time traces for both cases ---
figure(25); clf;

subplot(1, 2, 1); hold on;
plot(t1, X1(:,1), 'b-', 'LineWidth', 1.5);
ylabel('Membrane Potential V (mV)'); xlabel('time');
title('Q18 Case 1: Fixed at Rest (-60 mV)');
grid on;

subplot(1, 2, 2); hold on;
plot(t2, X2(:,1), 'b-', 'LineWidth', 1.5);
ylabel('Membrane Potential V (mV)'); xlabel('time');
title('Q18 Case 2: Fixed at End of Anodal Stimulus');
grid on;

%% ===================== Local functions (Part I & II) ====================
% Local functions defined at the end
function [value, isterminal, direction] = stopOutside(~, X)
    % Event function: stops integration once the trajectory leaves the
    % plotted phase-space window (used for the saddle manifolds above).
    in_bounds = (X(2) >= -90 && X(2) <= 70 && X(1) >= -0.1 && X(1) <= 1.1);
    value = in_bounds - 0.5;
    isterminal = 1;
    direction = 0;
end

function a = alpha_n(V)
    a = zeros(size(V));
    mask = abs(V + 50) < 1e-7;
    a(mask)  = 0.1;
    a(~mask) = -0.01*(V(~mask)+50)./(exp(-(V(~mask)+50)/10)-1);
end

function b = alpha_m(V)
    b = zeros(size(V));
    mask = abs(V + 35) < 1e-7;
    b(mask)  = 1;
    b(~mask) = -0.1*(V(~mask)+35)./(exp(-(V(~mask)+35)/10)-1);
end