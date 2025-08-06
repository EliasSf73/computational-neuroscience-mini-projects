% Homework2: Problem1_b.m
% -------------------------------------------------------------------------
% Purpose:
% This script simulates the Hodgkin-Huxley model with a specific input
% current (Ic) and then implements methods to detect the precise timing (τ)
% of each action potential (spike). It compares two methods: peak detection
% and simple threshold crossing. Finally, it visualizes the detected spike
% times (using peak detection) on the voltage trace plot.
%
% Based on HW2Prob1_a_Verification.m
% -------------------------------------------------------------------------

clear; clc; close all; % start fresh

fprintf('--- Hodgkin-Huxley Spike Timing Detection Script ---\n');

%% === Model Parameters (Consistent with previous step) ===
g_max = [36, 120, 0.3];        % [gK_max, gNa_max, gLeak] (mS/cm^2)
E_rev = [-12, 115, 10.613];    % [EK, ENa, ELeak] (mV)
Cm = 1.0;                      % Membrane Capacitance (uF/cm^2)
gK_max = g_max(1); gNa_max = g_max(2); gLeak = g_max(3);
EK = E_rev(1); ENa = E_rev(2); ELeak = E_rev(3);

%% === Simulation Setup ===
dt = 0.01;                     % Simulation time step (ms)
t_start = 0;                   % Simulation start time (ms)
t_end = 500;                   % Simulation end time (ms)
time_vector = t_start:dt:t_end;% Vector of time points
stim_onset_time = 100;         % Stimulus onset time (ms)
stim_duration_sec = (t_end - stim_onset_time) / 1000; % Duration for rate calc (s)

%% === Input Current (Use the value verified previously) ===
Ic_value = 6.2148;             % Input current amplitude (uA/cm^2) - Verified to give ~25Hz

fprintf('Using Input Current: Ic = %.4f uA/cm^2\n', Ic_value);
fprintf('Simulation Time Span: [%.1f, %.1f] ms\n', t_start, t_end);
fprintf('Stimulus Applied From: %.1f ms onwards\n', stim_onset_time);
fprintf('--------------------------------------------------\n');

%% === Initial Neuron State ===
membrane_potential = -10.0;      % Initial V (mV)
gating_vars = [0, 0, 1];         % Initial [n, m, h]

%% === Simulation Data Storage ===
voltage_trace = zeros(size(time_vector)); % Array to store V(t)
time_index = 1;

%% === Run the Simulation ===
fprintf('Simulation running...\n');
tic; % Start timer

% Loop through each time step
for current_time = time_vector

    % --- 1. Determine External Input Current ---
    if current_time >= stim_onset_time
        external_current = Ic_value;
    else
        external_current = 0;
    end

    % --- 2. Calculate gating variable rate constants (alpha & beta) ---
    V = membrane_potential;
    alpha_rates = zeros(1,3); beta_rates = zeros(1,3);

    % Rate constants for n (K activation) - Handle V=10 singularity
    if abs(V - 10) < 1e-7, alpha_rates(1) = 0.1; else, alpha_rates(1) = (10-V) / (100 * (exp((10-V)/10) - 1)); end
    beta_rates(1) = 0.125 * exp(-V/80);

    % Rate constants for m (Na activation) - Handle V=25 singularity
    if abs(V - 25) < 1e-7, alpha_rates(2) = 1.0; else, alpha_rates(2) = (25-V) / (10 * (exp((25-V)/10) - 1)); end
    beta_rates(2) = 4 * exp(-V/18);

    % Rate constants for h (Na inactivation)
    alpha_rates(3) = 0.07 * exp(-V/20);
    beta_rates(3) = 1 / (exp((30-V)/10) + 1);

    % --- 3. Update Gating Variables (Exponential Euler) ---
    tau_gates = 1.0 ./ (alpha_rates + beta_rates + eps);
    x_inf_gates = alpha_rates .* tau_gates;
    gating_vars = gating_vars .* (1 - dt ./ tau_gates) + x_inf_gates .* (dt ./ tau_gates);
    gating_vars = max(0, min(1, gating_vars)); % Clamp [0, 1]
    n = gating_vars(1); m = gating_vars(2); h = gating_vars(3);

    % --- 4. Calculate Ionic Conductances ---
    gK = gK_max * n^4; gNa = gNa_max * m^3 * h; gL = gLeak;
    conductances = [gK, gNa, gL];

    % --- 5. Calculate Ionic Currents ---
    ionic_currents = conductances .* (V - E_rev); % [IK, INa, IL]

    % --- 6. Update Membrane Potential (Forward Euler) ---
    dV = (external_current - sum(ionic_currents)) / Cm;
    membrane_potential = V + dt * dV;

    % --- 7. Record Data ---
    voltage_trace(time_index) = membrane_potential;
    time_index = time_index + 1;

end % End simulation loop
simulation_duration = toc;
fprintf('Simulation finished in %.2f seconds.\n', simulation_duration);
fprintf('--------------------------------------------------\n');

%% === Spike Timing Detection ===
% Now that we have the full voltage_trace, we can analyze it to find spike times.

fprintf('Detecting spike times...\n');

% --- Method 1: Peak Detection (Our primary method) ---
% Definition: Spike time τ is the time when the voltage reaches its maximum
% value during an action potential.
% Approach:
% 1. Find upward crossings of a voltage threshold.
% 2. After a crossing, find the maximum voltage value before the voltage
%    starts significantly decreasing (e.g., drops below threshold or peak is passed).
% 3. Record the time corresponding to that maximum voltage.
% 4. Implement a refractory mechanism to avoid detecting the same spike multiple times.

spike_times_peak = [];       % Array to store spike times found by peak method
spike_indices_peak = [];     % Array to store corresponding indices
peak_voltages = [];          % Store the peak voltage values
spike_detection_threshold = 0; % Use 0 mV threshold
in_spike_peak = false;       % State flag: Are we currently above threshold looking for a peak?
current_peak_V = -Inf;       % Voltage of the peak found so far within the current spike
current_peak_idx = -1;       % Index of the peak found so far

% Find the index corresponding to the stimulus onset time
stim_onset_index = find(time_vector >= stim_onset_time, 1);

% Iterate through the voltage trace *after* stimulus onset
for i = stim_onset_index:length(voltage_trace)
    V_now = voltage_trace(i);
    V_prev = voltage_trace(i-1); % Need voltage from the previous step

    % Detect upward threshold crossing *if not already looking for a peak*
    if ~in_spike_peak && V_prev < spike_detection_threshold && V_now >= spike_detection_threshold
        in_spike_peak = true;       % Start looking for the peak of this spike
        current_peak_V = V_now;     % Initialize peak search
        current_peak_idx = i;
    end

    % If we are looking for a peak
    if in_spike_peak
        % If current voltage is higher than the peak found so far, update the peak
        if V_now > current_peak_V
            current_peak_V = V_now;
            current_peak_idx = i;
        % If voltage starts dropping significantly after finding a potential peak
      
        % A simple check is if V starts decreasing after having crossed the threshold.
        elseif V_now < V_prev && current_peak_idx ~= -1 % Check V starts decreasing and we found a potential peak idx
            % We passed the peak. Record the time of the *maximum* voltage found.
            spike_indices_peak(end+1) = current_peak_idx;
            spike_times_peak(end+1) = time_vector(current_peak_idx);
            peak_voltages(end+1) = current_peak_V;

            % Reset for the next spike - must go below threshold before finding another
            in_spike_peak = false;
            current_peak_V = -Inf;
            current_peak_idx = -1;
        end
    end

     % A condition to reset if V drops below threshold even if peak wasn't formally declared
     % this handles cases where spike might not fully repolarize quickly
     if in_spike_peak && V_now < spike_detection_threshold
           % If we were tracking a peak, record it now as voltage dropped below threshold
           if current_peak_idx ~= -1
               spike_indices_peak(end+1) = current_peak_idx;
               spike_times_peak(end+1) = time_vector(current_peak_idx);
               peak_voltages(end+1) = current_peak_V;
           end
           % Reset state
           in_spike_peak = false;
           current_peak_V = -Inf;
           current_peak_idx = -1;
     end

end

fprintf(' -> Method 1 (Peak Detection): Found %d spikes.\n', length(spike_times_peak));
% Display the first few detected peak times
if ~isempty(spike_times_peak)
    fprintf('    Spike times (τ_peak): %.2f', spike_times_peak(1));
    for k=2:min(5, length(spike_times_peak)) % Print up to first 5
        fprintf(', %.2f', spike_times_peak(k));
    end
    if length(spike_times_peak) > 5, fprintf('...'); end
    fprintf(' ms\n');
end


% --- Method 2: Simple Threshold Crossing (Alternative for comparison) ---
% Definition: Spike time τ is the time when the voltage first crosses a
% specific threshold value in the upward direction.
% Approach:
% 1. Iterate through the voltage trace.
% 2. Detect when V(i-1) < threshold and V(i) >= threshold.
% 3. Record the time t(i).
% 4. Implement a simple refractory state: don't detect another crossing
%    until the voltage has dropped below the threshold again.

spike_times_threshold = [];  % Array to store spike times
below_threshold = true;      % State flag: Is the voltage currently below threshold?

% Iterate through the voltage trace *after* stimulus onset
for i = stim_onset_index:length(voltage_trace)
    V_now = voltage_trace(i);
    V_prev = voltage_trace(i-1);

    % Detect upward crossing only if we were previously below threshold
    if below_threshold && V_prev < spike_detection_threshold && V_now >= spike_detection_threshold
        spike_times_threshold(end+1) = time_vector(i); % Record time of crossing
        below_threshold = false; % Now we are above threshold, ignore further crossings
    end

    % Reset the flag when voltage drops below threshold again
    if ~below_threshold && V_now < spike_detection_threshold
        below_threshold = true; % Ready to detect the next spike crossing
    end
end
fprintf(' -> Method 2 (Threshold Crossing): Found %d spikes.\n', length(spike_times_threshold));
% Display the first few detected threshold times
if ~isempty(spike_times_threshold)
    fprintf('    Spike times (τ_thresh): %.2f', spike_times_threshold(1));
    for k=2:min(5, length(spike_times_threshold)) % Print up to first 5
        fprintf(', %.2f', spike_times_threshold(k));
    end
    if length(spike_times_threshold) > 5, fprintf('...'); end
    fprintf(' ms\n');
end


%% === Discussion: Comparison of Methods ===
fprintf('--------------------------------------------------\n');
fprintf('Discussion: Spike Timing Detection Methods\n');
fprintf('1. Peak Detection:\n');
fprintf('   - Pros: Defines spike time based on a clear physiological event (max depolarization).\n');
fprintf('           Less sensitive to the exact value of the detection threshold.\n');
fprintf('           Often considered a more accurate representation of the spike event time.\n');
fprintf('   - Cons: More computationally complex (requires post-simulation analysis and search).\n');
fprintf('           Requires careful logic to identify the peak correctly and avoid double counting.\n');
fprintf('2. Simple Threshold Crossing:\n');
fprintf('   - Pros: Very simple and computationally fast (can be done during simulation).\n');
fprintf('   - Cons: Highly sensitive to the chosen threshold value (different threshold -> different τ).\n');
fprintf('           The time doesn''t correspond to a unique physiological point, just the start of the upswing.\n');
fprintf('           More susceptible to noise if threshold is set too low.\n');
fprintf('Conclusion: Peak detection is generally preferred for precise spike timing analysis,\n');
fprintf('            while threshold crossing is simpler for basic spike counting or triggering events.\n');
fprintf('            Note that τ_peak will always occur slightly *after* τ_threshold for the same spike.\n');
fprintf('--------------------------------------------------\n');

%% === Plotting Results with Spike Timings (Using Peak Detection) ===

figure; % Create a new figure window
plot(time_vector, voltage_trace, 'b-', 'LineWidth', 1.5, 'DisplayName', 'V(t)'); % Plot voltage trace
hold on; % Keep the plot active

% Check if any spikes were detected by the peak method
if ~isempty(spike_times_peak)
    % Plot markers at the *peak* of each detected spike
    plot(spike_times_peak, peak_voltages, 'ro', ...
         'MarkerFaceColor', 'r', ...
         'MarkerSize', 7, ...
         'DisplayName', 'Spike Peak Time (τ_{peak})');
else
    fprintf('No spike peaks detected for plotting.\n');
end

% Add visual aids: stimulus onset line and frequency text
ylim_current = ylim; % Get y-axis limits after plotting V(t)
plot([stim_onset_time stim_onset_time], ylim_current, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off'); % Dashed line for Ic ON
text(stim_onset_time + 5, ylim_current(2)*0.9, 'I_c ON', 'Color', 'k', 'FontSize', 10);

% Calculate frequency based on peak detection count for annotation
num_spikes_peak = length(spike_times_peak);
if stim_duration_sec > 0
    frequency_peak = num_spikes_peak / stim_duration_sec;
else
    frequency_peak = 0;
end
text_content = sprintf('Freq (%.0f-%.0fms, Peaks): %.2f Hz', stim_onset_time, t_end, frequency_peak);
text(time_vector(end)*0.55, ylim_current(2)*0.8, text_content, 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

hold off; % Release the plot hold

% Add final plot details
xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');
title_str = sprintf('HH Voltage Trace with Detected Spike Peak Times (Ic = %.4f µA/cm^2)', Ic_value);
title(title_str);
grid on;
legend('show', 'Location', 'SouthEast'); % Show legend
ylim([-80 130]); % Set consistent y-axis limits

fprintf('Plot generated showing V(t) and detected spike peak times.\n');