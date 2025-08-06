% define parameters
% rng(9);  % fixed seed for reproducibility

dt = 0.1;             % time step (ms)
T = 1000;             % total simulation time (ms) = 1 second
t_array = 0:dt:T;
nSteps = length(t_array);

Erest = -65;          % resting potential (mV)
Eth = -55;            % threshold (mV)
tau = 20;             % membrane time constant (ms)
R = 10;               % membrane resistance (Ohm)
I_inj = 0;            % no constant injected current
sigma_noise = 8.5;  % from part (d)   
nTrials=50;




% run multiple trials and then calculate mean firing rate for each I_c
function mean_rate=measure_rate(I_c,sigma_noise,Erest,Eth,tau,R,dt,T,nTrials);
% array for each trial's firing rate
rates= zeros(1,nTrials);
for trial=1:nTrials
    rng(trial); 
     v = Erest * ones(1, round(T/dt)+1);
        spike_count = 0;
        
        % Euler's integration
        for i = 1:length(v)-1
            I_noise = sigma_noise * randn; % random noise
            % updating LIF voltage
            v(i+1) = v(i) + (dt/tau)*( (Erest - v(i)) + R*(I_c + I_noise) );
            
            % threshold crossing check
            if v(i+1) >= Eth
                spike_count = spike_count + 1;
                v(i+1) = Erest; % reset
            end
        end
        
        % Ccalculate the firing rate for this trial (spikes/sec)
        rates(trial) = spike_count / (T/1000); 
    end
    
    % return the mean across all trials
    mean_rate = mean(rates);3
end



% 10 values for I_c from 0 to Imax 
Imax = 2; 
Ic_vals = linspace(0, Imax, 10);

%  arrays for average firing rates for each condition
mean_rates_no_noise = zeros(1, length(Ic_vals));
mean_rates_noise = zeros(1, length(Ic_vals));

%% calculate 'response function' for each I_c Value
% case1: without noise (sigma_noise = 0) and case2: with noise (sigma_noise = 8.5)
for idx = 1:length(Ic_vals)
    current = Ic_vals(idx);
    % measure average firing rate of case1: without noise
    mean_rates_no_noise(idx) = measure_rate(current, 0, Erest, Eth, tau, R, dt, T, nTrials);
    % Measure average firing rate of case2: with noise (sigma_noise = 8.5)
    mean_rates_noise(idx) = measure_rate(current, 8.5, Erest, Eth, tau, R, dt, T, nTrials);
end

%% Plot the "response functions"
figure;
plot(Ic_vals, mean_rates_no_noise, 'bo-', 'LineWidth', 2);
hold on;
plot(Ic_vals, mean_rates_noise, 'ro-', 'LineWidth', 2);
xlabel('Injected Current I_c (mA)');
ylabel('Average Firing Rate (Hz)');
title('Neuron Response Function (f-I Curve) with and without Background Noise');
legend('No Noise (\sigma_{noise} = 0)', 'With Noise (\sigma_{noise} = 8.5)');
grid on;
hold off;