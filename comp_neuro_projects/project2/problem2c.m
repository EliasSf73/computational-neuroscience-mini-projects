% HW2Prob2c_CoupledLIF_Noise.m
% -------------------------------------------------------------------------
% Simulates two LIF neurons with independent Gaussian noise added to
% the input current. Finds the noise standard deviation (sigma) needed
% to achieve approximately 5 Hz standard deviation in instantaneous
% firing rate. Plots the resulting noisy spike trains.
% Assumes NO synaptic coupling (g_peak = 0).
% -------------------------------------------------------------------------

clear; clc; close all;

fprintf('--- LIF Neurons with Independent Noise ---\n');

%% === Parameters & Setup ===

% --- LIF Neuron Parameters (consistent) ---
Cm = 1.0;         % uF/cm^2
gL = 0.1;         % mS/cm^2
tau_m = Cm / gL;     % ms
V_rest = -65.0;   % mV  <- Using original V_rest for consistency w/ report
V_th = -55.0;     % mV
V_reset = -75.0;  % mV  <- Using corrected reset
t_ref = 2.0;      % ms

% --- Simulation Parameters ---
dt = 0.1;         % Time step (ms)
T_sim_final = 1000; % Final simulation duration for plot (ms)
T_sim_search = 5000;% Longer duration for statistics during search (ms)
time_vector_final = 0:dt:T_sim_final;
n_steps_final = length(time_vector_final);

% --- Noise Parameters ---
target_std_dev_rate = 5.0; % Target standard deviation of instantaneous rate (Hz)
Ie_base = 1.0430; % Constant baseline current from 2a (~25 Hz average)
sigma_found = NaN; % Value of sigma to be found

fprintf('Base LIF Params: tau_m=%.1f, V_rest=%.1f, V_th=%.1f, V_reset=%.1f, t_ref=%.1f\n', ...
    tau_m, V_rest, V_th, V_reset, t_ref);
fprintf('Using baseline Ie = %.4f uA/cm^2\n', Ie_base);
fprintf('Target std dev of inst. rate: %.1f Hz\n', target_std_dev_rate);
fprintf('--------------------------------------------------\n');


%% === Function to Simulate ONE Noisy LIF & Calculate Rate SD ===
% Simulates for T_duration, returns std dev of instantaneous rate
function rate_sd = calculate_noisy_rate_sd(Ie, sigma, T_duration, dt, Cm, gL, V_rest, V_th, V_reset, t_ref, tau_m)
    numSteps = round(T_duration / dt);
    V = V_rest;
    t_last_spike = -inf;
    spike_times = [];

    for i = 1:numSteps
        t_current = i * dt;
        time_since_last = t_current - t_last_spike;

        if time_since_last < t_ref
            V = V_reset;
        else
            % Add independent Gaussian noise to current
            I_noise = sigma * randn(); % Generate noise for this step
            I_total = Ie + I_noise;

            % Correct LIF update equation
            dV = ((V_rest - V) + (I_total / gL)) / tau_m;
            V = V + dt * dV;

            if V >= V_th
                spike_times(end+1) = t_current;
                V = V_reset;
                t_last_spike = t_current;
            end
        end
    end

    % Calculate rate standard deviation from ISIs
    if length(spike_times) > 1
        ISIs = diff(spike_times);
        if ~isempty(ISIs) & ISIs > 0 % Avoid division by zero if ISIs are bad
             rates_inst = 1000 ./ ISIs; % Convert ISIs (ms) to rates (Hz)
             rate_sd = std(rates_inst);
        else
             rate_sd = 0; % Or NaN, if no valid ISIs
        end
    else
        rate_sd = 0; % Not enough spikes to calculate SD
    end
end


%% === Find Sigma (Noise Strength) ===
% Using simple iterative search for clarity (Bisection could also be used)
fprintf('Searching for sigma to achieve ~%.1f Hz rate SD...\n', target_std_dev_rate);

sigma_low = 0; rate_sd_low = 0;
sigma_high = 0.5; % Initial guess for upper bound (adjust if needed)
rate_sd_high = calculate_noisy_rate_sd(Ie_base, sigma_high, T_sim_search, dt, Cm, gL, V_rest, V_th, V_reset, t_ref, tau_m);
fprintf('  Initial Check: sigma=%.3f -> Rate SD=%.2f Hz\n', sigma_high, rate_sd_high);

% Adjust initial sigma_high if necessary
iter_adjust = 0;
max_iter_adjust = 10;
while rate_sd_high < target_std_dev_rate && iter_adjust < max_iter_adjust
    fprintf('  Rate SD too low. Increasing sigma_high...\n');
    sigma_high = sigma_high * 2; % Double the guess
    rate_sd_high = calculate_noisy_rate_sd(Ie_base, sigma_high, T_sim_search, dt, Cm, gL, V_rest, V_th, V_reset, t_ref, tau_m);
    fprintf('  Check: sigma=%.3f -> Rate SD=%.2f Hz\n', sigma_high, rate_sd_high);
    iter_adjust = iter_adjust + 1;
end
if iter_adjust == max_iter_adjust
     error('Could not find sigma_high that gives Rate SD >= target.');
end

% Simple iterative refinement (like manual bisection steps)
tolerance_sigma = 0.005; % How close sigma needs to be
max_iter_search = 20;
iter = 0;
while (sigma_high - sigma_low) > tolerance_sigma && iter < max_iter_search
    sigma_mid = (sigma_low + sigma_high) / 2;
    rate_sd_mid = calculate_noisy_rate_sd(Ie_base, sigma_mid, T_sim_search, dt, Cm, gL, V_rest, V_th, V_reset, t_ref, tau_m);
    fprintf('  Iter %d: Test sigma=%.4f -> Rate SD=%.2f Hz\n', iter+1, sigma_mid, rate_sd_mid);

    if rate_sd_mid < target_std_dev_rate
        sigma_low = sigma_mid; % Need more noise
        rate_sd_low = rate_sd_mid;
    else
        sigma_high = sigma_mid; % Have enough or too much noise
        rate_sd_high = rate_sd_mid;
    end
    iter = iter + 1;
end

% Choose the sigma that gave the rate SD closest to target
if abs(rate_sd_low - target_std_dev_rate) < abs(rate_sd_high - target_std_dev_rate)
    sigma_found = sigma_low;
    final_rate_sd = rate_sd_low;
else
    sigma_found = sigma_high;
    final_rate_sd = rate_sd_high;
end


fprintf('--------------------------------------------------\n');
fprintf('Search complete (%d iterations).\n', iter);
fprintf('Using final noise sigma = %.4f uA/cm^2\n', sigma_found);
fprintf('This yields ~%.2f Hz standard deviation in instantaneous rate.\n', final_rate_sd);
fprintf('--------------------------------------------------\n');


%% === Simulate TWO Independent Noisy LIF Neurons (1 second) ===
% Uses Ie_base and the found sigma_found

V = [V_rest; V_rest];           % Voltages for N1, N2
t_last_spike = [-inf; -inf];    % Last spike times for N1, N2
spike_times_noisy = {[], []};   % Cell array to store spike times

fprintf('Running final 1s simulation for two noisy neurons...\n');
tic;
for i = 1:n_steps_final
    t = time_vector_final(i);

    % Update each neuron independently
    for neuron_idx = 1:2
        time_since_last = t - t_last_spike(neuron_idx);

        if time_since_last < t_ref
            V(neuron_idx) = V_reset; % Clamp during refractory
        else
            % Add independent Gaussian noise
            I_noise = sigma_found * randn();
            I_total = Ie_base + I_noise;

            % Integrate voltage
            dV = ((V_rest - V(neuron_idx)) + (I_total / gL)) / tau_m;
            V(neuron_idx) = V(neuron_idx) + dt * dV;

            % Check for spike
            if V(neuron_idx) >= V_th
                spike_times_noisy{neuron_idx}(end+1) = t; % Record spike time
                V(neuron_idx) = V_reset;                  % Reset voltage
                t_last_spike(neuron_idx) = t;             % Update last spike time
            end
        end
    end % End neuron loop
end % End time loop
toc;
fprintf('Simulation complete.\n');

% Calculate final average rates
rate1 = length(spike_times_noisy{1}) / (T_sim_final/1000);
rate2 = length(spike_times_noisy{2}) / (T_sim_final/1000);
fprintf('Avg Firing Rates: N1=%.1f Hz, N2=%.1f Hz\n', rate1, rate2);
fprintf('--------------------------------------------------\n');

%% === Plot Spike Trains (Raster Plot) ===
figure;
hold on;
colors = {'b', 'r'};
for neuron_idx = 1:2
    if ~isempty(spike_times_noisy{neuron_idx})
        plot(spike_times_noisy{neuron_idx}, neuron_idx * ones(size(spike_times_noisy{neuron_idx})), ...
             '|', 'Color', colors{neuron_idx}, 'MarkerSize', 10, 'LineWidth', 1.5);
    end
end
hold off;
ylim([0.5, 2.5]);
set(gca, 'YTick', [1, 2], 'YTickLabel', {'Neuron 1', 'Neuron 2'});
xlabel('Time (ms)'); ylabel('Neuron');
title(sprintf('Noisy LIF Spike Trains (Ie=%.2f, \\sigma=%.4f)', Ie_base, sigma_found));
xlim([0, T_sim_final]); grid on;

fprintf('Raster plot generated for noisy neurons.\n');