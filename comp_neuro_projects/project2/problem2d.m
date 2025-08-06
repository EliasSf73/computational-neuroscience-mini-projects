% HW2Prob2d_CoupledLIF_WithNoise.m
% -------------------------------------------------------------------------
% Simulates two mutually coupled LIF neurons WITH independent Gaussian noise.
% Repeats the scenarios from part (b) (Excitatory/Simultaneous,
% Inhibitory/Alternating) using the same parameters, but adds noise
% found in part (c) to observe its effect on firing patterns.
% -------------------------------------------------------------------------

clear; clc; close all;

fprintf('--- Coupled LIF Neuron Simulation with Noise ---\n');

%% === Parameters & Setup ===

% --- LIF Neuron Parameters (consistent) ---
Cm = 1.0;         % uF/cm^2
gL = 0.1;         % mS/cm^2
tau_m = Cm / gL;     % ms
V_rest = -65.0;   % mV
V_th = -55.0;     % mV
V_reset = -75.0;  % mV
t_ref = 2.0;      % ms

% --- Simulation Parameters ---
dt = 0.1;         % Time step (ms)
T_sim = 1000;     % Total simulation time (ms) = 1 second
time_vector = 0:dt:T_sim;
n_steps = length(time_vector);

% --- Noise Parameter (from part c) ---
sigma_found = 0.9102; % Noise standard deviation (uA/cm^2)

fprintf('Base LIF Params: tau_m=%.1f, V_rest=%.1f, V_th=%.1f, V_reset=%.1f, t_ref=%.1f\n', ...
    tau_m, V_rest, V_th, V_reset, t_ref);
fprintf('Simulation: T=%d ms, dt=%.2f ms\n', T_sim, dt);
fprintf('Noise level: sigma = %.4f uA/cm^2\n', sigma_found);
fprintf('--------------------------------------------------\n');


%% === Scenario (i) B: Excitatory Synapses + Noise - Target: Simultaneous Firing? ===
fprintf('\n--- Running Scenario: Excitatory + Noise - Simultaneous? ---\n');

% -- Parameters for this scenario (same synaptic/Ie as noise-free case) --
E_syn_exc = 0;      % Excitatory reversal potential
g_peak_exc_sync = 0.15; % Synaptic strength (from part b)
tau_syn_exc = 5.0;   % Synaptic decay time (from part b)
Ie_exc_sync = 1.05;  % Base Input current (from part b)

fprintf('Params: E_syn=%.1f mV, g_peak=%.3f mS/cm^2, tau_syn=%.1f ms, Ie_base=%.3f uA/cm^2\n', ...
        E_syn_exc, g_peak_exc_sync, tau_syn_exc, Ie_exc_sync);

% -- Initialization --
V = [V_rest; V_rest];
t_last_spike = [-inf; -inf];
g_syn = [0; 0];
spike_times_exc_noisy = {[], []};

% -- Simulation Loop with Noise --
for i = 1:n_steps
    t = time_vector(i);

    % Decay synaptic conductances
    g_syn = g_syn * (1 - dt/tau_syn_exc);
    g_syn = max(0, g_syn);

    % Calculate total input current including independent noise for each neuron
    noise_input = sigma_found * randn(2, 1); % Generate 2 independent noise values
    I_total = Ie_exc_sync + noise_input;     % Add noise to base current

    % Update each neuron
    for neuron_idx = 1:2
        % Synaptic current for this neuron
        I_syn = g_syn(neuron_idx) * (V(neuron_idx) - E_syn_exc);

        % Refractory check
        time_since_last = t - t_last_spike(neuron_idx);
        is_refractory = (time_since_last < t_ref);

        if is_refractory
            V(neuron_idx) = V_reset;
        else
            % Integrate voltage (LIF + Synapse + Noise)
            current_drive = I_total(neuron_idx); % Use noisy current for this neuron
            dV = ( -gL*(V(neuron_idx) - V_rest) - I_syn + current_drive ) / Cm;
            V(neuron_idx) = V(neuron_idx) + dt * dV;

            % Check for spike
            if V(neuron_idx) >= V_th
                V(neuron_idx) = V_reset;
                t_last_spike(neuron_idx) = t;
                spike_times_exc_noisy{neuron_idx}(end+1) = t;
                other_neuron = 3 - neuron_idx;
                g_syn(other_neuron) = g_syn(other_neuron) + g_peak_exc_sync;
            end
        end
    end
end

% -- Plotting for this scenario --
figure;
hold on; colors = {'b', 'r'};
plot(spike_times_exc_noisy{1}, 1 * ones(size(spike_times_exc_noisy{1})), '|', 'Color', colors{1}, 'MarkerSize', 10, 'LineWidth', 1.5);
plot(spike_times_exc_noisy{2}, 2 * ones(size(spike_times_exc_noisy{2})), '|', 'Color', colors{2}, 'MarkerSize', 10, 'LineWidth', 1.5);
hold off; ylim([0.5, 2.5]); set(gca, 'YTick', [1, 2], 'YTickLabel', {'Neuron 1', 'Neuron 2'});
xlabel('Time (ms)'); ylabel('Neuron'); xlim([0, T_sim]); grid on;
title(sprintf('Excitatory Coupling + Noise (\\sigma=%.3f) - Simultaneous Firing?', sigma_found));
fprintf(' -> Scenario complete. Check plot for effect of noise on synchrony.\n');
fprintf('--------------------------------------------------\n');


%% === Scenario (ii) A: Inhibitory Synapses + Noise - Target: Alternating Firing? ===
fprintf('\n--- Running Scenario: Inhibitory + Noise - Alternating? ---\n');

% -- Parameters for this scenario (same synaptic/Ie as noise-free case) --
E_syn_inh = -80.0;    % Inhibitory reversal potential
g_peak_inh_alt = 0.05; % Synaptic strength (from part b - successful alternation)
tau_syn_inh = 5.0;    % Synaptic decay time (from part b - successful alternation)
Ie_inh_alt = 2.0;     % Base Input current (from part b - successful alternation)

fprintf('Params: E_syn=%.1f mV, g_peak=%.3f mS/cm^2, tau_syn=%.1f ms, Ie_base=%.3f uA/cm^2\n', ...
        E_syn_inh, g_peak_inh_alt, tau_syn_inh, Ie_inh_alt);

% -- Initialization --
V = [V_rest; V_rest + 0.1]; % Slight asymmetry
t_last_spike = [-inf; -inf];
g_syn = [0; 0];
spike_times_inh_noisy = {[], []};

% -- Simulation Loop with Noise --
for i = 1:n_steps
    t = time_vector(i);

    % Decay synaptic conductances
    g_syn = g_syn * (1 - dt/tau_syn_inh);
    g_syn = max(0, g_syn);

    % Calculate total input current including independent noise
    noise_input = sigma_found * randn(2, 1);
    I_total = Ie_inh_alt + noise_input;

    % Update each neuron
    for neuron_idx = 1:2
        I_syn = g_syn(neuron_idx) * (V(neuron_idx) - E_syn_inh);
        time_since_last = t - t_last_spike(neuron_idx);
        is_refractory = (time_since_last < t_ref);

        if is_refractory
            V(neuron_idx) = V_reset;
        else
            current_drive = I_total(neuron_idx);
            dV = ( -gL*(V(neuron_idx) - V_rest) - I_syn + current_drive ) / Cm;
            V(neuron_idx) = V(neuron_idx) + dt * dV;

            if V(neuron_idx) >= V_th
                V(neuron_idx) = V_reset;
                t_last_spike(neuron_idx) = t;
                spike_times_inh_noisy{neuron_idx}(end+1) = t;
                other_neuron = 3 - neuron_idx;
                g_syn(other_neuron) = g_syn(other_neuron) + g_peak_inh_alt;
            end
        end
    end
end

% -- Plotting for this scenario --
figure;
hold on;
plot(spike_times_inh_noisy{1}, 1 * ones(size(spike_times_inh_noisy{1})), '|', 'Color', colors{1}, 'MarkerSize', 10, 'LineWidth', 1.5);
plot(spike_times_inh_noisy{2}, 2 * ones(size(spike_times_inh_noisy{2})), '|', 'Color', colors{2}, 'MarkerSize', 10, 'LineWidth', 1.5);
hold off; ylim([0.5, 2.5]); set(gca, 'YTick', [1, 2], 'YTickLabel', {'Neuron 1', 'Neuron 2'});
xlabel('Time (ms)'); ylabel('Neuron'); xlim([0, T_sim]); grid on;
title(sprintf('Inhibitory Coupling + Noise (\\sigma=%.3f) - Alternating Firing?', sigma_found));
fprintf(' -> Scenario complete. Check plot for effect of noise on alternation.\n');
fprintf('--------------------------------------------------\n');


