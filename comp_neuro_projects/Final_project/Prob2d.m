%% Problem 2d: Generate 5 Poisson Spike Trains (30 Hz) in 10 ms Bins

% We model spike generation as a (discrete‐time) Poisson process
% with mean rate λ = 30 Hz, dt = 10 ms ⇒ expected p = λ·dt = 0.3 probability of a spike per bin.

clear; close all; clc;

%% 1) Parameters
dt  = 10;           % ms per bin
T   = 100;          % total window (ms)
n   = T/dt;         % number of bins (10)
p   = 0.3;          % probability of ≥1 spike in a bin
M   = 5;            % number of example trains to generate

%% 2) Generate binary spike‐train matrix (M × n)
spikes = rand(M, n) < p;    % each entry is 1 with prob p, else 0

%% 3) Display the binary patterns
disp('Five example Poisson spike patterns (rows):');
disp(spikes);

%% 4) Raster plot of the example trains
figure('Name','Example Poisson Spike Trains (30 Hz)','NumberTitle','off');
hold on;
for i = 1:M
    % Convert bin indices to time: we mark the start of each bin
    spike_bins = find(spikes(i,:));
    spike_times = (spike_bins - 1) * dt;  
    % Plot as vertical tick marks
    plot(spike_times, i * ones(size(spike_times)), 'k|', 'MarkerSize',12);
end
hold off;
xlabel('Time (ms)');
ylabel('Trial index');
title('Five Example Poisson Spike Trains at 30 Hz');
ylim([0.5, M + 0.5]);
yticks(1:M);
xlim([0, T]);
grid on;
