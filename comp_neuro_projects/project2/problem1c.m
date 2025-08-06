% Homework 2: Problem1c_ResponseFunction.m
% -------------------------------------------------------------------------
% Purpose:
% Calculates and plots the response function (spike count vs. injected current)
% for the Hodgkin-Huxley model. Investigates behavior at high currents.
% Uses simple threshold crossing for spike counting and the linear
% approximation for Exponential Euler on gating variables.
% -------------------------------------------------------------------------

clear; clc; close all; % Start fresh

fprintf('--- Hodgkin-Huxley Response Function Calculation ---\n');

%% === Model Parameters (Consistent) ===
g_max = [36, 120, 0.3];        % [gK_max, gNa_max, gLeak] (mS/cm^2)
E_rev = [-12, 115, 10.613];    % [EK, ENa, ELeak] (mV)
Cm = 1.0;                      % Membrane Capacitance (uF/cm^2)
gK_max = g_max(1); gNa_max = g_max(2); gLeak = g_max(3);
EK = E_rev(1); ENa = E_rev(2); ELeak = E_rev(3);

%% === Simulation Setup (Consistent) ===
dt = 0.01;                     % Simulation time step (ms)
t_start = 0;                   % Simulation start time (ms)
t_end = 500;                   % Simulation end time (ms)
time_vector = t_start:dt:t_end;% Vector of time points
stim_onset_time = 100;         % Stimulus onset time (ms)
stim_duration_sec = (t_end - stim_onset_time) / 1000; % Duration for rate calc (s)
spike_detection_threshold = 0; % Voltage threshold for spike counting (mV)

%% === Range of Input Currents ===
% defining the range and steps for the injected current Ic
Ic_min = 0;
Ic_max = 120;        % range can be updated
num_Ic_steps = 21;  % Number of steps (e.g., steps of 2 uA/cm^2)
Ic_values = linspace(Ic_min, Ic_max, num_Ic_steps); % Vector of currents

% an array to store results (spike counts for each Ic)
spike_counts_per_Ic = zeros(size(Ic_values));

fprintf('Testing %d current values from %.2f to %.2f uA/cm^2...\n', num_Ic_steps, Ic_min, Ic_max);
progress_timer = tic; % Start a timer for the loop

%% === Loop Through Each Current Value ===
for k = 1:length(Ic_values)
    current_Ic = Ic_values(k);

    % --- Reset Initial Neuron State for each simulation ---
    membrane_potential = -10.0;      % Initial V (mV)
    gating_vars = [0, 0, 1];         % Initial [n, m, h]

    % --- reset spike counting variables for this Ic ---
    current_spike_count = 0;
    voltage_previous_step = membrane_potential;
    below_threshold = true; % State for threshold crossing detection

    % --- Run the HH Simulation for this Ic value ---
    for current_time = time_vector

        % 1. Determine External Input Current
        if current_time >= stim_onset_time
            external_current = current_Ic; % Use the current Ic from the outer loop
        else
            external_current = 0;
        end

        % 2. Calculate Rate Constants (Alpha & Beta) - with singularity handling
        V = membrane_potential;
        alpha_rates = zeros(1,3); beta_rates = zeros(1,3);
        if abs(V - 10) < 1e-7, alpha_rates(1) = 0.1; else, alpha_rates(1) = (10-V) / (100 * (exp((10-V)/10) - 1)); end
        beta_rates(1) = 0.125 * exp(-V/80);
        if abs(V - 25) < 1e-7, alpha_rates(2) = 1.0; else, alpha_rates(2) = (25-V) / (10 * (exp((25-V)/10) - 1)); end
        beta_rates(2) = 4 * exp(-V/18);
        alpha_rates(3) = 0.07 * exp(-V/20);
        beta_rates(3) = 1 / (exp((30-V)/10) + 1);

        % 3. Update Gating Variables (Using Linear Approx of Exp Euler)
        tau_gates = 1.0 ./ (alpha_rates + beta_rates + eps);
        x_inf_gates = alpha_rates .* tau_gates;
        gating_vars = gating_vars .* (1 - dt ./ tau_gates) + x_inf_gates .* (dt ./ tau_gates);
        gating_vars = max(0, min(1, gating_vars)); % Clamp [0, 1]
        n = gating_vars(1); m = gating_vars(2); h = gating_vars(3);

        % 4. Calculate Ionic Conductances
        gK = gK_max * n^4; gNa = gNa_max * m^3 * h; gL = gLeak;
        conductances = [gK, gNa, gL];

        % 5. Calculate Ionic Currents
        ionic_currents = conductances .* (V - E_rev); % [IK, INa, IL]

        % 6. Update Membrane Potential
        dV = (external_current - sum(ionic_currents)) / Cm;
        membrane_potential = V + dt * dV;

        % --- Spike Counting (Simple Threshold Crossing) ---
        % Count only during the stimulus period [stim_onset_time, t_end]
        if current_time >= stim_onset_time
            V_now = membrane_potential;
            % Detect upward crossing if previously below threshold
            if below_threshold && voltage_previous_step < spike_detection_threshold && V_now >= spike_detection_threshold
                current_spike_count = current_spike_count + 1;
                below_threshold = false; % Now above threshold, wait to go below
            end
            % Reset flag if voltage drops below threshold again
            if ~below_threshold && V_now < spike_detection_threshold
                below_threshold = true; % Ready for next spike crossing
            end
            voltage_previous_step = V_now; % Update previous V
        else
             % Still need to update previous V even before stimulus starts
             voltage_previous_step = membrane_potential;
        end
        % --- End Spike Counting ---

    end % End simulation loop for one Ic value

    % --- Store the result for this Ic ---
    spike_counts_per_Ic(k) = current_spike_count;
   
    % fprintf('  Ic = %.4f uA/cm^2 --> Spikes = %d\n', current_Ic, current_spike_count);

end % End loop over Ic_values

fprintf('Finished calculations in %.2f seconds.\n', toc(progress_timer));
fprintf('--------------------------------------------------\n');

%% === Plot the Response Function ===
figure; % New figure for the response curve
plot(Ic_values, spike_counts_per_Ic, 'bs-', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
xlabel('Injected Current I_c (\muA/cm^2)');
ylabel(sprintf('Spike Count (%d-%d ms)', stim_onset_time, t_end));
title('Hodgkin-Huxley Neuron Response Function (Spike Count)');
grid on;
xlim([Ic_min - 1, Ic_max + 1]); % 

% ------- Ends Here -----











% --------------------------------------------------\n');

%% ==: Simulate and Plot Example of High Current Behavior ===
Ic_high_test = 120.0; % Choose a significantly higher current likely to cause block
fprintf('Running simulation for a high current example: Ic = %.1f uA/cm^2...\n', Ic_high_test);

% --- Reset and Run Simulation ---
membrane_potential = -10.0; gating_vars = [0, 0, 1];
V_trace_high_Ic = zeros(size(time_vector)); % Store trace
time_index = 1;
% Re-run the simulation loop logic just for this high current
for current_time = time_vector
    if current_time >= stim_onset_time, external_current = Ic_high_test; else, external_current = 0; end
    V = membrane_potential; alpha_rates = zeros(1,3); beta_rates = zeros(1,3);
    if abs(V - 10) < 1e-7, alpha_rates(1) = 0.1; else, alpha_rates(1) = (10-V) / (100 * (exp((10-V)/10) - 1)); end; beta_rates(1) = 0.125 * exp(-V/80);
    if abs(V - 25) < 1e-7, alpha_rates(2) = 1.0; else, alpha_rates(2) = (25-V) / (10 * (exp((25-V)/10) - 1)); end; beta_rates(2) = 4 * exp(-V/18);
    alpha_rates(3) = 0.07 * exp(-V/20); beta_rates(3) = 1 / (exp((30-V)/10) + 1);
    tau_gates = 1.0 ./ (alpha_rates + beta_rates + eps); x_inf_gates = alpha_rates .* tau_gates;
    gating_vars = gating_vars .* (1 - dt ./ tau_gates) + x_inf_gates .* (dt ./ tau_gates); gating_vars = max(0, min(1, gating_vars));
    n = gating_vars(1); m = gating_vars(2); h = gating_vars(3);
    gK = gK_max * n^4; gNa = gNa_max * m^3 * h; gL = gLeak; conductances = [gK, gNa, gL];
    ionic_currents = conductances .* (V - E_rev);
    dV = (external_current - sum(ionic_currents)) / Cm; membrane_potential = V + dt * dV;
    V_trace_high_Ic(time_index) = membrane_potential; % Store V
    time_index = time_index + 1;
end

% --- Plot the High Current Trace ---
figure; % New figure for the high current example
plot(time_vector, V_trace_high_Ic, 'r-', 'LineWidth', 1.5);
hold on;
plot([stim_onset_time stim_onset_time], ylim, 'k--', 'LineWidth', 1); % Stimulus onset line
plot(xlim, [spike_detection_threshold spike_detection_threshold], 'g:', 'LineWidth', 1); % Show threshold
hold off;
xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');
title(sprintf('High Current Example: V(t) with Ic = %.1f µA/cm^2', Ic_high_test));
legend('V(t)', 'Ic ON', 'Spike Threshold', 'Location', 'best');
grid on;
ylim([-80 130]); % Keep consistent y-axis

fprintf('Generated plot showing behavior at high injected current.\n');