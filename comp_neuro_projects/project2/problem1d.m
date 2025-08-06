% Homewrok2: Problem 1d_LIF_ResponseFunction.m
% -------------------------------------------------------------------------
% Purpose:
% Calculates and plots the response function (spike count vs. injected current)
% for a simple Leaky Integrate-and-Fire (LIF) neuron model.
% Compares its response profile to the Hodgkin-Huxley model.
% -------------------------------------------------------------------------

clear; clc; close all; % Start fresh

fprintf('--- LIF Neuron Response Function Calculation ---\n');

%% === LIF Model Parameters ===
% Use parameters comparable to HH where applicable, but typical for LIF
Cm = 1.0;          % Membrane Capacitance (uF/cm^2) - Same as HH
gL = 0.3;          % Leak Conductance (mS/cm^2) - Same as HH leak
V_rest = -65.0;    % Resting Potential (mV) - More typical resting value
V_thresh = -55.0;  % Spike Threshold (mV)
V_reset = -70.0;   % Reset Potential after spike (mV)
t_ref = 1.0;       % Absolute Refractory Period (ms)

% Calculate membrane time constant (optional, but informative)
tau_m = Cm / gL; % ms
fprintf('LIF Parameters: Cm=%.1f, gL=%.1f, V_rest=%.1f, V_thresh=%.1f, V_reset=%.1f, t_ref=%.1f, tau_m=%.2f\n', ...
    Cm, gL, V_rest, V_thresh, V_reset, t_ref, tau_m);

%% === simulation setup ===
dt = 0.01;                     % Simulation time step (ms)
t_start = 0;                   % Simulation start time (ms)
t_end = 500;                   % Simulation end time (ms)
time_vector = t_start:dt:t_end;% Vector of time points
stim_onset_time = 100;         % Stimulus onset time (ms)
stim_duration_sec = (t_end - stim_onset_time) / 1000; % Duration for rate calc (s)

%% === Range of Input Currents ===
% range might need to be adjusted here 
Ic_min = 0;
Ic_max = 70;       % 
num_Ic_steps = 21; % 
Ic_values = linspace(Ic_min, Ic_max, num_Ic_steps); % Vector of currents

% an array to store results
spike_counts_per_Ic_lif = zeros(size(Ic_values));

fprintf('Testing %d current values from %.2f to %.2f uA/cm^2...\n', num_Ic_steps, Ic_min, Ic_max);
progress_timer = tic;

%% === Loop Through Each Current Value ===
for k = 1:length(Ic_values)
    current_Ic = Ic_values(k);

    % --- Reset Initial Neuron State ---
    membrane_potential = V_rest; % Start at rest
    time_since_last_spike = inf; % Not in refractory period initially

    % --- Reset Spike Counting Variables ---
    current_spike_count = 0;

    % --- Run the LIF Simulation for this Ic value ---
    for current_time = time_vector

        % Check if currently in refractory period
        is_refractory = (time_since_last_spike < t_ref);

        % Update time since last spike
        time_since_last_spike = time_since_last_spike + dt;

        % If refractory, hold voltage at reset and do nothing else
        if is_refractory
            membrane_potential = V_reset; % Clamp voltage during refractory
            continue; % Skip integration and threshold check
        end

        % --- Integration Step (if not refractory) ---
        % 1. Determine External Input Current
        if current_time >= stim_onset_time
            external_current = current_Ic;
        else
            external_current = 0;
        end

        % 2. Update Membrane Potential using Forward Euler
        % dV/dt = (-gL*(V - V_rest) + I_ext) / Cm
        V = membrane_potential; % Use temp variable for clarity
        dV = ((-gL * (V - V_rest) + external_current) / Cm);
        membrane_potential = V + dV * dt;

        % --- Check for Spike ---
        if membrane_potential >= V_thresh
            % Spike occurred!
            membrane_potential = V_reset; % Reset voltage
            time_since_last_spike = 0;    % Reset refractory timer

            % Count spike ONLY if it occurs during the stimulus period
            if current_time >= stim_onset_time
                current_spike_count = current_spike_count + 1;
            end
        end % End spike check

    end % End simulation loop for one Ic value

    % --- Store the result for this Ic ---
    spike_counts_per_Ic_lif(k) = current_spike_count;
    % Print progress
    % fprintf('  Ic = %.4f uA/cm^2 --> Spikes = %d\n', current_Ic, current_spike_count);

end % end the loop over Ic_values

fprintf('Finished calculations in %.2f seconds.\n', toc(progress_timer));
fprintf('--------------------------------------------------\n');

%% === Plot the LIF Response Function ===
figure;
plot(Ic_values, spike_counts_per_Ic_lif, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
xlabel('Injected Current I_c (\muA/cm^2)');
ylabel(sprintf('Spike Count (%d-%d ms)', stim_onset_time, t_end));
title('Leaky Integrate-and-Fire (LIF) Neuron Response Function');
grid on;
xlim([Ic_min - 0.2, Ic_max + 0.2]); % 