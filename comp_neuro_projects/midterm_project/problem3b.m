%%  poisson spike-train under time-varying rate r(t)
clear;clc;
% fundamental parameters
N=5; % 5 trials
dt=1e-3;
t= 0:dt:1;
% sinusoidal parameters from problem 3a
c1= 20;
c2= 4*pi;
c3=pi/2;
c4=30;
% compute r(t) --> utilizing vectorized operation
r= c1*sin(c2*t-c3)+c4;

% 3) Generate spikes ---------------------------------------------
% Preallocate: rows=time bins, cols=trials
spikes= false(numel(t),N);
rng(0,'twister');             % for reproducibility
for ti=1:numel(t)
    p_current= r(ti)*dt;
    spikes(ti,:)= rand(1,N)<p_current;
end

%%  -------------------plot r(t) and raster ----------------------------------

figure('Position',[200 200 600 400]);
subplot(2,1,1);
plot(t, r, 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('r(t) [Hz]');
title('Time‑varying firing rate from Part 3a');
grid on;

subplot(2,1,2);
hold on;


for trial = 1:N
    % extract spike times for this trial
    ts = t(spikes(:,trial));
    % plot each spike as a vertical tick
    plot(ts, trial*ones(size(ts)), 'k|', 'MarkerSize', 8);
end
ylim([0.5 N+0.5]);
xlabel('Time (s)'); ylabel('Trial');
title(sprintf('Raster plot: Poisson spikes (N = %d trials)',N));
set(gca,'YDir','reverse');  % so Trial 1 is on top
box on;
