%% LIF neuron with noise: multiple trials
clear; clc;
rng(35);

% parameters
dt = 0.1;         % time step in ms
T = 1000;         % simulation duration in ms (1 sec)
t_array = 0:dt:T;
nSteps = length(t_array);

Erest = -65;      % resting membrane potential in mV
tau = 20;         % membrane time constant in ms
Eth = -55;        % spike threshold in mV
R = 10;           % membrane resistance in Ohm
I_inj = 0;        % zero constant injected current

% Set of sigma_noise values to test
sigma_vals = [0, 5, 8,8.5, 9, 10];
nSigmas = length(sigma_vals);

% Number of trials for each sigma_noise
nTrials = 100;

% matrices for firing rates and spike times
firing_rates = zeros(nSigmas, nTrials);  % rows for each sigma_noise, columns for each trial
all_spike_times = cell(nSigmas, nTrials);  % store spike times for each trial

% Loop over each sigma_noise value and each trial
for j = 1:nSigmas
    sigma_noise = sigma_vals(j);
    for trial = 1:nTrials
 
        v = Erest * ones(1, nSteps);  % initialize membrane potential
        spike_times = [];
        for i = 1:nSteps-1
            I_noise = sigma_noise * randn;  % noise current sample from N(0, sigma_noise^2)
            v(i+1) = v(i) + (dt/tau) * (Erest - v(i) + R * I_noise);
            if v(i+1) >= Eth
                spike_times = [spike_times, t_array(i)];  % recording spike time
                v(i+1) = Erest;  % reset membrane potential
            end
        end
        all_spike_times{j, trial} = spike_times;
        num_spikes = length(spike_times);
        firing_rates(j, trial) = num_spikes / (T/1000);  % convert ms to seconds for Hz
    end
    % average firing rate for current sigma_noise value
    avg_rate = mean(firing_rates(j, :));
    fprintf('For sigma_noise = %g, Average firing rate = %.2f Hz\n', sigma_vals(j), avg_rate);
end

%% Plot 1: average firing rate vs. sigma_noise
figure;
% Compute mean and standard deviation across trials for each sigma_noise value
mean_rates = mean(firing_rates, 2);
std_rates = std(firing_rates, 0, 2);
errorbar(sigma_vals, mean_rates, std_rates, 'o-', 'LineWidth', 2);
xlabel('\sigma_{noise} (mA)');
ylabel('Average Firing Rate (Hz)');
title('Average Firing Rate vs. Noise Standard Deviation');
grid on;






























% figure;
% hold on;
% For each trial, plot vertical lines at each spike time, with the trial number on the y-axis.
% for trial = 1:nTrials
%     spike_times = all_spike_times{sel_index, trial};
%     for k = 1:length(spike_times)
%         % Draw a vertical line at the spike time, spanning trial-0.4 to trial+0.4
%         line([spike_times(k) spike_times(k)], [trial-0.8, trial+0.8], 'Color', 'k', 'LineWidth', 1);
%     end
% end
% hold off;
% xlabel('Time (ms)');
% ylabel('Trial');
% title(['Raster Plot for \sigma_{noise} = ' num2str(selected_sigma) ' mA']);
% xlim([0 T]);
% ylim([0 nTrials+1]);
