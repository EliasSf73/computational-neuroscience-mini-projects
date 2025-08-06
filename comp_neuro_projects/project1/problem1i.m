%%  comparing neuron behavior at I_weak = 1.5 mA ( value between 0 and 2 mA)
%  simulation parameters
dt = 0.1;          % time step in ms
T = 1000;          % simulation duration in ms (1 second)
t_array = 0:dt:T;
nSteps = length(t_array);
Erest = -65;       % resting potential (mV)
Eth = -55;         % threshold (mV)
tau = 20;          % membrane time constant (ms)
R = 10;            % membrane resistance (Ohm)
I_weak = 1.2;      % chosen input current (mA)

nTrials = 5;       % number of repeated trials

%% case1: Simulation witout background noise
sigma_noisesnf = 0;  % noise-free scenario

all_spike_times_nf = cell(nTrials,1);
firing_rates_nf = zeros(1, nTrials);

for trial = 1:nTrials
    rng(trial); % fixed seed for reproducibility for each trial
    v = Erest * ones(1, nSteps);
    spike_times = [];
    for i = 1:nSteps-1
        I_noise = sigma_noisesnf * randn;  % always 0 here sice sigma_noise=0
        v(i+1) = v(i) + (dt/tau)*(Erest - v(i) + R*(I_weak + I_noise));
        if v(i+1) >= Eth
            spike_times = [spike_times, t_array(i)];
            v(i+1) = Erest;  % reset after spike
        end
    end
    all_spike_times_nf{trial} = spike_times;
    firing_rates_nf(trial) = length(spike_times) / (T/1000);
end

%% case2: Sismulation with background noise
sigma_noise_n = 8.5;  % background noise value found from part (d)

all_spike_times_n = cell(nTrials,1);
firing_rates_n = zeros(1, nTrials);

for trial = 1:nTrials
    rng(trial+100);  % using a different seed offset for noise condition
    v = Erest * ones(1, nSteps);
    spike_times = [];
    for i = 1:nSteps-1
        I_noise = sigma_noise_n * randn;
        v(i+1) = v(i) + (dt/tau)*(Erest - v(i) + R*(I_weak + I_noise));
        if v(i+1) >= Eth
            spike_times = [spike_times, t_array(i)];
            v(i+1) = Erest;
        end
    end
    all_spike_times_n{trial} = spike_times;
    firing_rates_n(trial) = length(spike_times) / (T/1000);
end

% calculating mean firing rates 
mean_rate_nf = mean(firing_rates_nf);
mean_rate_n = mean(firing_rates_n);
fprintf('Noise-free condition: Mean firing rate = %.2f Hz\n', mean_rate_nf);
fprintf('With noise: Mean firing rate = %.2f Hz\n', mean_rate_n);

%% plot raster plots for comparison
figure;

% case1: raster plot for noise-free condition
subplot(1,2,1);
hold on;
for trial = 1:nTrials
    spike_times = all_spike_times_nf{trial};
    for s = 1:length(spike_times)
        line([spike_times(s) spike_times(s)], [trial-0.4, trial+0.4], 'Color', 'b', 'LineWidth', 1.5);
    end
end
hold off;
xlabel('Time (ms)');
ylabel('Trial');
title('Noise-Free (I_{weak} = 1.5 mA)');
xlim([0, T]);
ylim([0, nTrials+1]);

% case2: raster plot for noise condition
subplot(1,2,2);
hold on;
for trial = 1:nTrials
    spike_times = all_spike_times_n{trial};
    for s = 1:length(spike_times)
        line([spike_times(s) spike_times(s)], [trial-0.4, trial+0.4], 'Color', 'r', 'LineWidth', 1.5);
    end
end
hold off;
xlabel('Time (ms)');
ylabel('Trial');
title('With Noise (I_{weak} = 1.2 mA, \sigma_{noise}=8.5)');
xlim([0, T]);
ylim([0, nTrials+1]);

sgtitle('Raster Plots Comparing Neuron Behavior at I_{weak} = 1.2 mA');

