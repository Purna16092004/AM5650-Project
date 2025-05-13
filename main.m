%Purna Ananthkrishnan - 20th April 2025
%Indian Institute of Technology Madras
%Coupling of structure and wake oscillators in vortex-induced vibrations by M.L. Facchinetti, E. de Langre, and F. Biolley.

clc; clear; close all;
format long g;


%% Parameters from the paper (Section 4.2)
St = 0.2;
epsilon = 0.3;
A = 12;
gamma = 0.8;
xi = 3.1e-3;
M = 2e-4;
mu = 0.05/M;

%% Figure 5

% Simulation setup
Ur_vec_up = linspace(3.5, 6.5, 150);
Ur_vec_down = fliplr(Ur_vec_up);
tspan = [0 5000];
dt = 0.1;
t = tspan(1):dt:tspan(2);

omega1_up = zeros(size(Ur_vec_up));
omega1_down = zeros(size(Ur_vec_down));
omega2_up = zeros(size(Ur_vec_up));
omega2_down = zeros(size(Ur_vec_down));
omega3_up = zeros(size(Ur_vec_up));
omega3_down = zeros(size(Ur_vec_down));

function [freq, amp] = dominant_frequency(t, y)
    Fs = 1 / mean(diff(t));
    y = y - mean(y);
    L = length(y);

    N = 2^nextpow2(4*L);  % Zero padding
    Y = fft(y, N);
    P2 = abs(Y / N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    f = Fs * (0:(N/2)) / N;

    [amp, idx] = max(P1);
    if idx > 1 && idx < length(P1)
        alpha = P1(idx - 1);
        beta = P1(idx);
        gamma = P1(idx + 1);
        p = (gamma - alpha) / (2 * (2 * beta - alpha - gamma));
        freq = f(idx) + p * (f(2) - f(1));
    else
        freq = f(idx);
    end
end

function [t, Y] = ode4(odefun, tspan, y0)
    t = tspan;
    Y = zeros(length(t), length(y0));
    Y(1, :) = y0;

    for i = 1:length(t) - 1
        h = t(i+1) - t(i);
        k1 = odefun(t(i), Y(i,:)');
        k2 = odefun(t(i) + h/2, Y(i,:)' + h*k1/2);
        k3 = odefun(t(i) + h/2, Y(i,:)' + h*k2/2);
        k4 = odefun(t(i) + h, Y(i,:)' + h*k3);
        Y(i+1, :) = Y(i,:) + h*((k1 + 2*k2 + 2*k3 + k4)/6)';
    end
end


% figure 5a

for k = 1:2
    if k == 1
        Ur_vec = Ur_vec_up;
    else
        Ur_vec = Ur_vec_down;
    end

    omega = zeros(size(Ur_vec));

    for i = 1:length(Ur_vec)
        Ur = Ur_vec(i);
        delta = 1 / (St * Ur);
        d_coeff = 2 * xi * delta + gamma / mu;

        if i == 1
            y0 = [0.1; 0; 2; 0];
        else
            y0 = Y(end, :);
        end

        [~, Y] = ode4(@(t,y) viv_odea(t, y, delta, d_coeff, epsilon, A, M), t, y0);
        Y1 = Y;

        % Analyze last 40% of signal
        idx_start = floor(0.6 * length(t));
        idx = idx_start:length(t);
        y_signal = Y(idx, 1);

        if all(abs(diff(y_signal)) < 1e-6) || any(isnan(y_signal))
            omega(i) = NaN;
            continue;
        end

        try
            [freq, ~] = dominant_frequency(t(idx), y_signal);
            omega(i) = 2 * pi * freq;  % just ω, in rad/s
            fprintf('Ur = %.2f → omega = %.4f\n', Ur, omega(i));
        catch
            warning('FFT failed at Ur = %.2f', Ur);
            omega(i) = NaN;
        end
    end

    if k == 1
        omega1_up = omega;
    else
        omega1_down = omega;
    end
end


% Plotting
figure;
hold on;
plot(Ur_vec_up, omega1_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, omega1_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('Angular Frequency, \omega (rad/s)');
xlim([3.5 6.5]);
ylim([0.75 1.25]);
legend('Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('\omega as a Function of U_r — Displacement Coupling');
drawnow;

% figure 5b
for k = 1:2
    if k == 1
        Ur_vec = Ur_vec_up;
    else
        Ur_vec = Ur_vec_down;
    end

    omega = zeros(size(Ur_vec));

    for i = 1:length(Ur_vec)
        Ur = Ur_vec(i);
        delta = 1 / (St * Ur);
        d_coeff = 2 * xi * delta + gamma / mu;

        if i == 1
            y0 = [0.1; 0; 2; 0];
        else
            y0 = Y(end, :);
        end

        [~, Y] = ode4(@(t,y) viv_odeb(t, y, delta, d_coeff, epsilon, A, M), t, y0);
        Y2 = Y;

        % Analyze last 40% of signal
        idx_start = floor(0.6 * length(t));
        idx = idx_start:length(t);
        y_signal = Y(idx, 1);

        if all(abs(diff(y_signal)) < 1e-6) || any(isnan(y_signal))
            omega(i) = NaN;
            continue;
        end

        try
            [freq, ~] = dominant_frequency(t(idx), y_signal);
            omega(i) = 2 * pi * freq;  % just ω, in rad/s
            fprintf('Ur = %.2f → omega = %.4f\n', Ur, omega(i));
        catch
            warning('FFT failed at Ur = %.2f', Ur);
            omega(i) = NaN;
        end
    end

    if k == 1
        omega2_up = omega;
    else
        omega2_down = omega;
    end
end

% Plotting
figure;
hold on;
plot(Ur_vec_up, omega2_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, omega2_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('Angular Frequency, \omega (rad/s)');
xlim([3.5 6.5]);
ylim([0.75 1.25]);
legend('Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('\omega as a Function of U_r — Velocity Coupling');
drawnow;

% figure 5c
for k = 1:2
    if k == 1
        Ur_vec = Ur_vec_up;
    else
        Ur_vec = Ur_vec_down;
    end

    omega = zeros(size(Ur_vec));

    for i = 1:length(Ur_vec)
        Ur = Ur_vec(i);
        delta = 1 / (St * Ur);
        d_coeff = 2 * xi * delta + gamma / mu;

        if i == 1
            y0 = [0.1; 0; 2; 0];
        else
            y0 = Y(end, :);
        end

        [~, Y] = ode4(@(t,y) viv_odec(t, y, delta, d_coeff, epsilon, A, M), t, y0);
        Y3 = Y;

        % Analyze last 40% of signal
        idx_start = floor(0.6 * length(t));
        idx = idx_start:length(t);
        y_signal = Y(idx, 1);

        if all(abs(diff(y_signal)) < 1e-6) || any(isnan(y_signal))
            omega(i) = NaN;
            continue;
        end

        try
            [freq, ~] = dominant_frequency(t(idx), y_signal);
            omega(i) = 2 * pi * freq;  % just ω, in rad/s
            fprintf('Ur = %.2f → omega = %.4f\n', Ur, omega(i));
        catch
            warning('FFT failed at Ur = %.2f', Ur);
            omega(i) = NaN;
        end
    end

    if k == 1
        omega3_up = omega;
    else
        omega3_down = omega;
    end
end

% Plotting
figure;
hold on;
plot(Ur_vec_up, omega3_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, omega3_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('Angular Frequency, \omega (rad/s)');
xlim([3.5 6.5]);
ylim([0.8 1.2]);
legend('Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('\omega as a Function of U_r — Acceleration Coupling');
drawnow;

%% Figures 6 & 7

Ur_vec_up = linspace(3.5, 6.5, 150);
Ur_vec_down = fliplr(Ur_vec_up);

q01_up = zeros(length(Ur_vec_up));
q01_down = zeros(length(Ur_vec_up));

q02_up = zeros(length(Ur_vec_up));
q02_down = zeros(length(Ur_vec_up));

q03_up = zeros(length(Ur_vec_up));
q03_down = zeros(length(Ur_vec_up));

y01_up = zeros(length(Ur_vec_up));
y01_down = zeros(length(Ur_vec_up));

y02_up = zeros(length(Ur_vec_up));
y02_down = zeros(length(Ur_vec_up));

y03_up = zeros(length(Ur_vec_up));
y03_down = zeros(length(Ur_vec_up));

for i = 1:length(Ur_vec_up)

    Ur = Ur_vec_up(i);
    delta = 1/(St*Ur);

    q01_up(i) = wake_amplitude(omega1_up(i),delta,1);
    q02_up(i) = wake_amplitude(omega2_up(i),delta,2);
    q03_up(i) = wake_amplitude(omega3_up(i),delta,3);

    q01_down(i) = wake_amplitude(omega1_down(i),delta,1);
    q02_down(i) = wake_amplitude(omega2_down(i),delta,2);
    q03_down(i) = wake_amplitude(omega3_down(i),delta,3);

    y01_up(i) = structure_amplitude(omega1_up(i),delta,q01_up(i));
    y02_up(i) = structure_amplitude(omega2_up(i),delta,q02_up(i));
    y03_up(i) = structure_amplitude(omega3_up(i),delta,q03_up(i));

    y01_down(i) = structure_amplitude(omega1_down(i),delta,q01_down(i));
    y02_down(i) = structure_amplitude(omega2_down(i),delta,q02_down(i));
    y03_down(i) = structure_amplitude(omega3_down(i),delta,q03_down(i));

end

q01_up = smoothdata(q01_up, 'movmean', 1);
q01_down = smoothdata(q01_down, 'movmedian', 10);
q01_down = flip(q01_down);

q02_up = smoothdata(q02_up, 'gaussian', 5);
q02_down = filloutliers(q02_down, 'pchip');
q02_down = flip(q02_down);

q03_up = smoothdata(q03_up, 'lowess', 9);
q03_down = smoothdata(q03_down, 'movmean', 9);
q03_down = flip(q03_down);

y01_up = smoothdata(y01_up, 'movmean', 5);
y01_down = smoothdata(y01_down, 'movmedian', 15);
y01_down = flip(y01_down);

y02_up = smoothdata(y02_up, 'movmean', 5);
y02_down = smoothdata(y02_down, 'movmedian', 10);


y03_up = smoothdata(y03_up, 'movmean', 5);
y03_down = smoothdata(y03_down, 'movmean', 10);


%Plotting
figure;
hold on;
plot(Ur_vec_up, q01_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, q01_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('q0');
xlim([3.5 6.5]);
ylim([1 3]);
h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);
legend([h1, h2], 'Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('Variation of q_0 with U_ — Displacement Coupling');
drawnow;

figure;
hold on;
plot(Ur_vec_up, q02_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, q02_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('q0');
xlim([3.5 6.5]);
ylim([1 3]);
h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);
legend([h1, h2], 'Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('Variation of q_0 with U_ — Velocity Coupling');
drawnow;

figure;
hold on;
plot(Ur_vec_up, q03_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, q03_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('q_0');
xlim([3.5 6.5]);
ylim([1 3]);
h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);
legend([h1, h2], 'Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('Variation of q_0 with U_r — Acceleration Coupling');
drawnow;

y01_up(y01_up < 1e-6) = NaN;
y01_down(y01_down < 1e-6) = NaN;
figure;
hold on;
plot(Ur_vec_up, y01_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, y01_down, 'k--', 'LineWidth', 1.5);
h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);
legend([h1, h2], 'Increasing U_r', 'Decreasing U_r');
xlabel('Reduced Velocity, U_r');
ylabel('y_0');
xlim([3.5 6.5]);
ylim([0 0.1]);
grid on;
set(gca, 'FontSize', 12);
title('Variation of y_0 with U_r — Displacement Coupling');
drawnow;

y02_up(y02_up < 1e-6) = NaN;
y02_down(y02_down < 1e-6) = NaN;
figure;
hold on;
plot(Ur_vec_up, y02_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, y02_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('y_0');
xlim([3.5 6.5]);
ylim([0 0.1]);
h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);
legend([h1, h2], 'Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('Variation of y_0 with U_r — Velocity Coupling');
drawnow;
y03_up(y03_up < 1e-6) = NaN;
y03_down(y03_down < 1e-6) = NaN;
figure;
hold on;
plot(Ur_vec_up, y03_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, y03_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('y_0');
xlim([3.5 6.5]);
ylim([0 0.1]);
h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);
legend([h1, h2], 'Increasing U_r', 'Decreasing U_r');
grid on;
set(gca, 'FontSize', 12);
title('Variation of y_0 with U_r — Acceleration Coupling');
drawnow;


%% Novel Result 

% Reduced velocity sweep
Ur_vec_up = linspace(3.5, 6.5, 100);
Ur_vec_down = fliplr(Ur_vec_up);

% Time
tspan = [0 4000];
dt = 0.1;
t = tspan(1):dt:tspan(2);
y0 = [0.1; 0; 0; 0];  % Initial condition

omega_up = zeros(size(Ur_vec_up));
omega_down = zeros(size(Ur_vec_down));

% Sweep UR up and down
for k = 1:2
    if k == 1
        Ur_vec = Ur_vec_up;
    else
        Ur_vec = Ur_vec_down;
    end
    omega = zeros(size(Ur_vec));

    for i = 1:length(Ur_vec)
        Ur = Ur_vec(i);
        delta = 1 / (St * Ur);

        if i == 1
            y0 = [0.1; 0; 0; 0];
        else
            y0 = Y(end,:);
        end

        % Run simulation
        [~, Y] = ode4(@(t,y) viv_nonlinear(t, y, xi, gamma, mu, delta, M, epsilon), t, y0);

        % Frequency analysis
        idx = round(0.6 * length(t)):length(t);
        y_signal = Y(idx,1);
        if all(abs(diff(y_signal)) < 1e-6) || any(isnan(y_signal))
            omega(i) = NaN;
            continue;
        end
        [freq, ~] = dominant_frequency(t(idx), y_signal);
        omega(i) = 2 * pi * freq;
    end

    if k == 1
        omega_up = omega;
    else
        omega_down = omega;
    end
end

% Plot
figure;
hold on;
plot(Ur_vec_up, omega_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, omega_down, 'k--', 'LineWidth', 1.5);
xlabel('Reduced Velocity, U_r');
ylabel('Angular Frequency \omega (rad/s)');
xlim([3.5 6.5]);
ylim([0.95 1.05]);
title('\omega as a Function of U_r — Nonlinear Coupling');
grid on;
set(gca, 'FontSize', 12);

q0_up = zeros(length(Ur_vec_up));
q0_down = zeros(length(Ur_vec_up));

y0_up = zeros(length(Ur_vec_up));
y0_down = zeros(length(Ur_vec_up));


for i = 1:length(Ur_vec_up)

    Ur = Ur_vec_up(i);
    delta = 1/(St*Ur);

    q0_up(i) = wake_amplitude(omega_up(i),delta,1);
    q0_down(i) = wake_amplitude(omega_down(i),delta,1);
    y0_up(i) = structure_amplitude(omega_up(i),delta,q0_up(i));
    y0_down(i) = structure_amplitude(omega_down(i),delta,q0_down(i));

end

q0_up = smoothdata(q0_up, 'movmean', 1);
q0_down = smoothdata(q0_down, 'movmedian', 10);
q0_down = flip(q0_down);
y0_up = smoothdata(y0_up, 'movmean', 5);
y0_down = smoothdata(y0_down, 'movmedian', 15);
y0_down = flip(y0_down);


%Plotting
y0_up(y0_up < 1e-12) = NaN;
y0_down(y0_down < 1e-12) = NaN;
figure;
hold on;
plot(Ur_vec_up, y0_up, 'k-', 'LineWidth', 1.5);
plot(Ur_vec_down, y0_down, 'k--', 'LineWidth', 1.5);
h1 = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.5);
legend([h1, h2], 'Increasing U_r', 'Decreasing U_r');
xlabel('Reduced Velocity, U_r');
ylabel('y_0');
xlim([3.5 6.5]);
ylim([0 0.1]);
grid on;
set(gca, 'FontSize', 12);
title('Variation of y_0 with U_r — Nonlinear Coupling');
drawnow;


