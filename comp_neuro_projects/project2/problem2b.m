% HW2Prob2b_CoupledLIF.m
% -------------------------------------------------------------------------
% Simulates two mutually coupled Leaky Integrate-and-Fire (LIF) neurons.
% Explores excitatory and inhibitory coupling to achieve simultaneous
% and alternating firing patterns, respectively.
% -------------------------------------------------------------------------

clear; clc; close all;

fprintf('--- Coupled LIF Neuron Simulation ---\n');

%% === Parameters & Setup ===

% --- LIF Neuron Parameters (consistent units: ms, mV, nA/cm^2) ---
Cm = 1.0;         % Membrane Capacitance (uF/cm^2 == mS*ms/cm^2)
gL = 0.1;         % Leak conductance (mS/cm^2)
tau_m = Cm / gL;     % Membrane time constant (ms)
V_rest = -65.0;   % Resting potential (mV)
V_th = -55.0;     % Threshold potential (mV)
V_reset = -75.0;  % Reset potential (mV) - Corrected from previous examples
t_ref = 2.0;      % Refractory period (ms)

% --- Simulation Parameters ---
dt = 0.1;         % Time step (ms)
T_sim = 1000;     % Total simulation time (ms) = 1 second
time_vector = 0:dt:T_sim;
n_steps = length(time_vector);

% --- Synapse Base Parameters (they will be adjusted per scenario) ---
tau_syn = 5.0;    % Synaptic decay time constant (ms) 
g_peak = 0.1;     % Peak synaptic conductance increase per spike (mS/cm^2) - our initial guess
E_syn = 0;        % Synaptic reversal potential (mV) - will be set per scenario

% --- Input Current ---
% we will start near the value found in 2a, tho we may need adjustment per scenario
Ie_base = 1.05;   % Base constant input current (nA/cm^2 or uA/cm^2 if Cm=1)

fprintf('Base LIF Params: tau_m=%.1f, V_rest=%.1f, V_th=%.1f, V_reset=%.1f, t_ref=%.1f\n', ...
    tau_m, V_rest, V_th, V_reset, t_ref);
fprintf('Simulation: T=%d ms, dt=%.2f ms\n', T_sim, dt);
fprintf('--------------------------------------------------\n');


%% === Scenario (i) B: Excitatory Synapses - Target: Simultaneous Firing ===
fprintf('\n--- Running Scenario: Excitatory - Simultaneous ---\n');

% -- Parameters for this scenario --
E_syn_exc = 0;      % Excitatory reversal potential
g_peak_exc_sync = 0.15; % Synaptic strength (tuned for synchrony, might need adjustment)
tau_syn_exc = 5.0;   % Synaptic decay time
Ie_exc_sync = 1.05;  % Input current ()

fprintf('Params: E_syn=%.1f mV, g_peak=%.3f mS/cm^2, tau_syn=%.1f ms, Ie=%.3f uA/cm^2\n', ...
        E_syn_exc, g_peak_exc_sync, tau_syn_exc, Ie_exc_sync);

% -- Initialization --
V = [V_rest; V_rest];           % Voltages for N1, N2
t_last_spike = [-inf; -inf];    % Last spike times for N1, N2
g_syn = [0; 0];                 % Synaptic conductance received by N1, N2
spike_times_exc_sync = {[], []};% Cell array to store spike times

% -- Simulation Loop --
for i = 1:n_steps
    t = time_vector(i);
    Ie_now = Ie_exc_sync; % Constant current for this scenario

    % Decay synaptic conductances (linear approximation)
    g_syn = g_syn * (1 - dt/tau_syn_exc);
    g_syn = max(0, g_syn); % Ensure conductance doesn't go negative

    % Update each neuron
    for neuron_idx = 1:2
        % Calculate synaptic current for this neuron
        I_syn = g_syn(neuron_idx) * (V(neuron_idx) - E_syn_exc);

        % Check refractory period
        time_since_last = t - t_last_spike(neuron_idx);
        is_refractory = (time_since_last < t_ref);

        if is_refractory
            V(neuron_idx) = V_reset; % Clamp during refractory
        else
            % Integrate voltage (LIF + Synapse) - Correct leak term sign
            dV = ( -gL*(V(neuron_idx) - V_rest) - I_syn + Ie_now ) / Cm;
            V(neuron_idx) = V(neuron_idx) + dt * dV;

            % Check for spike
            if V(neuron_idx) >= V_th
                V(neuron_idx) = V_reset;              % Reset voltage
                t_last_spike(neuron_idx) = t;       % Update last spike time
                spike_times_exc_sync{neuron_idx}(end+1) = t; % Record spike time

                % Update *other* neuron's synaptic conductance
                other_neuron = 3 - neuron_idx; % Quick way to get 2 if j=1, 1 if j=2
                g_syn(other_neuron) = g_syn(other_neuron) + g_peak_exc_sync;
            end % End spike check
        end % End refractory check
    end % End neuron loop
end % End time loop

% -- Plotting for this scenario --
figure;
hold on;
colors = {'b', 'r'};
plot(spike_times_exc_sync{1}, 1 * ones(size(spike_times_exc_sync{1})), '|', 'Color', colors{1}, 'MarkerSize', 10, 'LineWidth', 1.5);
plot(spike_times_exc_sync{2}, 2 * ones(size(spike_times_exc_sync{2})), '|', 'Color', colors{2}, 'MarkerSize', 10, 'LineWidth', 1.5);
hold off;
ylim([0.5, 2.5]); set(gca, 'YTick', [1, 2], 'YTickLabel', {'Neuron 1', 'Neuron 2'});
xlabel('Time (ms)'); ylabel('Neuron'); xlim([0, T_sim]); grid on;
title(sprintf('Excitatory Coupling (g_{peak}=%.2f) - Simultaneous Firing?', g_peak_exc_sync));
fprintf(' -> Scenario complete. Check plot for synchrony.\n');
fprintf('--------------------------------------------------\n');


%% === Scenario (ii) A: Inhibitory Synapses - Target: Alternating Firing ===
fprintf('\n--- Running Scenario: Inhibitory - Alternating ---\n');

% -- Parameters for this scenario --
E_syn_inh = -80.0;    % Inhibitory reversal potential (below V_reset)
g_peak_inh_alt = 0.05; % Synaptic strength (tuned for alternation)
tau_syn_inh = 5.0;   % Slower decay can help stabilize alternation
Ie_inh_alt = 2;     % Input current (needs to be strong enough to overcome inhibition)

fprintf('Params: E_syn=%.1f mV, g_peak=%.3f mS/cm^2, tau_syn=%.1f ms, Ie=%.3f uA/cm^2\n', ...
        E_syn_inh, g_peak_inh_alt, tau_syn_inh, Ie_inh_alt);

% -- Initialization --
% Introducing slight asymmetry to encourage alternation
V = [V_rest; V_rest + 0.1];     % Start V2 slightly higher
t_last_spike = [-inf; -inf];
g_syn = [0; 0];
spike_times_inh_alt = {[], []};

% -- Simulation Loop (Identical logic, different parameters) --
for i = 1:n_steps
    t = time_vector(i);
    Ie_now = Ie_inh_alt; % Use current for this scenario

    % Decay synaptic conductances
    g_syn = g_syn * (1 - dt/tau_syn_inh);
    g_syn = max(0, g_syn);

    % Update each neuron
    for neuron_idx = 1:2
        I_syn = g_syn(neuron_idx) * (V(neuron_idx) - E_syn_inh); % Use inhibitory E_syn
        time_since_last = t - t_last_spike(neuron_idx);
        is_refractory = (time_since_last < t_ref);

        if is_refractory
            V(neuron_idx) = V_reset;
        else
            dV = ( -gL*(V(neuron_idx) - V_rest) - I_syn + Ie_now ) / Cm;
            V(neuron_idx) = V(neuron_idx) + dt * dV;
            if V(neuron_idx) >= V_th
                V(neuron_idx) = V_reset;
                t_last_spike(neuron_idx) = t;
                spike_times_inh_alt{neuron_idx}(end+1) = t;
                other_neuron = 3 - neuron_idx;
                g_syn(other_neuron) = g_syn(other_neuron) + g_peak_inh_alt; % Use inhibitory g_peak
            end
        end
    end
end

% -- Plotting for this scenario --
figure;
hold on;
plot(spike_times_inh_alt{1}, 1 * ones(size(spike_times_inh_alt{1})), '|', 'Color', colors{1}, 'MarkerSize', 10, 'LineWidth', 1.5);
plot(spike_times_inh_alt{2}, 2 * ones(size(spike_times_inh_alt{2})), '|', 'Color', colors{2}, 'MarkerSize', 10, 'LineWidth', 1.5);
hold off;
ylim([0.5, 2.5]); set(gca, 'YTick', [1, 2], 'YTickLabel', {'Neuron 1', 'Neuron 2'});
xlabel('Time (ms)'); ylabel('Neuron'); xlim([0, T_sim]); grid on;
title(sprintf('Inhibitory Coupling (g_{peak}=%.2f) - Alternating Firing?', g_peak_inh_alt));
fprintf(' -> Scenario complete. Check plot for alternation.\n');
fprintf('--------------------------------------------------\n');


