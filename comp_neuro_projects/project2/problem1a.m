% Homework2:  problem1_a.m
% Hodgkin-Huxley simulation with constant current injection for 25 Hz firing

% Parameters
g = [36, 120, 0.3];            % Max conductances (mS/cm^2): K, Na, Leak (Based on usage in code)
E = [-12, 115, 10.613];        % Reversal potentials (mV): K, Na, Leak
Cm = 1;                        % Membrane Capacitance (uF/cm^2) - Making explicit
dt = 0.01;                     % Time step (ms)
tspan = 0:dt:500;              % Simulation time (ms)
stim_onset = 100;              % Time when stimulus starts (ms)
stim_duration_sec = (tspan(end) - stim_onset) / 1000; % Duration for rate calc (s)

% searching for Ic that causes ~25 Hz firing in the [100, 500]ms window
Ic_low = 2;
Ic_high = 18;
target_spikes = 25 * stim_duration_sec; % CORRECTED version: Target spikes for 25 Hz in 400ms
fprintf('Target number of spikes (%.0fms window): %.1f\n', stim_duration_sec*1000, target_spikes);
tolerance = 0.01; % Tolerance value for Ic search

% --- Bisection Search ---
num_search_iter = 0;
max_search_iter = 100;
while (Ic_high - Ic_low) > tolerance && num_search_iter < max_search_iter
    Ic_mid = (Ic_low + Ic_high) / 2;
    % Passing stim_onset to simulate function for correct counting
    spikes = simulate(Ic_mid, g, E, Cm, dt, tspan, stim_onset);
    fprintf(' Search Iter %d: Ic = %.4f uA/cm^2 -> Spikes = %d\n', num_search_iter+1, Ic_mid, spikes);
    if spikes < target_spikes
        Ic_low = Ic_mid;
    else
        Ic_high = Ic_mid;
    end
    num_search_iter = num_search_iter + 1;
end

Ic_opt = (Ic_low + Ic_high) / 2;
fprintf('--------------------------------------------------\n');
fprintf('Search complete.\n');
fprintf('Optimal Ic found for ~%.1f Hz firing: %.4f µA/cm^2\n', 25, Ic_opt);
fprintf('--------------------------------------------------\n');


% --- Running final simulation to record V(t) and verify spike count ---
V = -10; % the original initial V from  lecture material
x = [0, 0, 1]; % the original initial x = [n, m, h] from lecture material
V_trace = zeros(size(tspan));
i = 1;

% Variables for final verification
final_spike_count = 0;
prev_V_final = V; % Initial previous V for spike counting
spike_threshold = 0; % Threshold for spike detection

fprintf('Running final simulation with Ic = %.4f uA/cm^2...\n', Ic_opt);
for t = tspan
    % Setting external current
    if t >= stim_onset
        I_ext = Ic_opt;
    else
        I_ext = 0;
    end

    % --- Alpha and Beta functions with division-by-zero handling ---
    Alpha = zeros(1,3); % Preallocating Alpha vector
    Beta = zeros(1,3);  % Preallocating Beta vector

    % Alpha/Beta for n (x(1)) - K channel activation
    V_eff_n = V - 10;
    if abs(V_eff_n) < 1e-7 % Handle V near 10 for alpha_n
         Alpha(1) = 0.1; % L'Hopital's rule result ( just to be careful)
    else
         Alpha(1) = (10-V)/(100*(exp((10-V)/10)-1));
    end
    Beta(1) = 0.125*exp(-V/80);

    % Alpha/Beta for m (x(2)) - Na channel activation
    V_eff_m = V - 25;
    if abs(V_eff_m) < 1e-7 % Handle V near 25 for alpha_m
         Alpha(2) = 1.0; % L'Hopital's rule result
    else
         Alpha(2) = (25-V)/(10*(exp((25-V)/10)-1));
    end
    Beta(2) = 4*exp(-V/18);

    % Alpha/Beta for h (x(3)) - Na channel inactivation
    Alpha(3) = 0.07*exp(-V/20);
    Beta(3) = 1/(exp((30-V)/10)+1);
    % --- End Alpha/Beta calculation ---

    % Time constants and steady states
    tau = 1./(Alpha + Beta);
    x_0 = Alpha .* tau;

    % Exponential Euler integration for gating variables
    x = (1 - dt ./ tau) .* x + (dt ./ tau) .* x_0;
    % Clamp gating variables to physically realistic range [0,1]
    x = max(0, min(1, x));

    % Conductances ( g(1)=gK, g(2)=gNa, g(3)=gL)
    gnmh(1) = g(1) * x(1)^4;       % K conductance (using n=x(1))
    gnmh(2) = g(2) * x(2)^3 * x(3); % Na conductance (using m=x(2), h=x(3))
    gnmh(3) = g(3);                % Leak conductance

    % Ionic currents
    I = gnmh .* (V - E); % I = [IK, INa, IL]

    % Membrane potential update (using explicit Cm)
    V = V + (dt/Cm) * (I_ext - sum(I));

    % Record voltage trace
    V_trace(i) = V;
    i = i + 1; % Increment index

    % --- Spike Counting for Final Verification ---
    % detecting spikes, during the stimulus period,  using threshold crossing
    if t >= stim_onset && prev_V_final < spike_threshold && V >= spike_threshold
        final_spike_count = final_spike_count + 1;
    end
    prev_V_final = V; % Update previous V for next iteration's check
    % --- End Spike Counting ---

end % End of simulation loop

% calculating final frequency
if stim_duration_sec > 0
    final_frequency = final_spike_count / stim_duration_sec;
else
    final_frequency = 0;
end

fprintf('Final simulation complete.\n');
fprintf('Spikes detected in [%.0f, %.0f] ms window: %d\n', stim_onset, tspan(end), final_spike_count);
fprintf('Resulting frequency: %.2f Hz\n', final_frequency);
fprintf('--------------------------------------------------\n');


% --- Plotting ---
figure;
plot(tspan, V_trace, 'b');
hold on;
% stimulus onset
plot([stim_onset stim_onset], ylim, 'k--', 'LineWidth', 1);
text(stim_onset + 5, max(ylim)*0.9, 'I_c ON', 'Color', 'k');
% frequency text
text(tspan(end)*0.6, max(ylim)*0.8, sprintf('Freq (%.0f-%.0fms): %.2f Hz', stim_onset, tspan(end), final_frequency), 'FontSize', 10, 'BackgroundColor', 'w');
hold off;
xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');
title(sprintf('Membrane Potential with Ic = %.4f µA/cm^2', Ic_opt));
grid on;
ylim([-80 130]); % I might adjust ylim based on observation of voltage range

% --- Local function required by the bisection search  above---
% arguments Cm and stim_onset passed now
function spike_count = simulate(Ic, g, E, Cm, dt, tspan, stim_onset)
    V = -10;                   % Initial membrane potential (mV)
    x = [0, 0, 1];             % Gating variables: n, m, h
    spike_count = 0;
    prev_V = V;
    spike_threshold = 0; % Threshold for spike detection

    for t = tspan
        % Set external current
        if t >= stim_onset
            I_ext = Ic;
        else
            I_ext = 0;
        end

        % --- Alpha and Beta functions with division-by-zero handling ---
        Alpha = zeros(1,3);
        Beta = zeros(1,3);

        % Alpha/Beta for n (x(1))
        V_eff_n = V - 10;
        if abs(V_eff_n) < 1e-7
             Alpha(1) = 0.1;
        else
             Alpha(1) = (10-V)/(100*(exp((10-V)/10)-1));
        end
        Beta(1) = 0.125*exp(-V/80);

        % Alpha/Beta for m (x(2))
        V_eff_m = V - 25;
        if abs(V_eff_m) < 1e-7
             Alpha(2) = 1.0;
        else
             Alpha(2) = (25-V)/(10*(exp((25-V)/10)-1));
        end
        Beta(2) = 4*exp(-V/18);

        % Alpha/Beta for h (x(3))
        Alpha(3) = 0.07*exp(-V/20);
        Beta(3) = 1/(exp((30-V)/10)+1);
        % --- End Alpha/Beta calculation ---


        % Time constants and steady states
        % integrating small epsilon to prevent division by zero if Alpha+Beta is zero
        tau = 1./(Alpha + Beta + eps);
        x_0 = Alpha .* tau;

        % Exponential Euler integration for gating variables
        x = (1 - dt ./ tau) .* x + (dt ./ tau) .* x_0;
        x = max(0, min(1, x)); % Clamp

        % Conductances
        gnmh(1) = g(1) * x(1)^4;       % K
        gnmh(2) = g(2) * x(2)^3 * x(3); % Na
        gnmh(3) = g(3);                % Leak

        % Ionic currents
        I_ion = gnmh .* (V - E); % [IK, INa, IL]

        % Membrane potential update
        V = V + (dt/Cm) * (I_ext - sum(I_ion));

        % --- Spike Detection ---
        % Counting spikes only during the stimulus period using threshold crossing
        if t >= stim_onset && prev_V < spike_threshold && V >= spike_threshold
            spike_count = spike_count + 1;
        end
        prev_V = V; % Update previous V for next iteration
    end % End loop over time
end % End simulate function