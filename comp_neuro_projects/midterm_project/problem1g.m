%% --------------------Problem 1g: autocorrelation of spikes with varying firing rate----------- 
r0      = 50;            % mean rate (Hz)
A       = 50;            % modulation amplitude (Hz)
Tosc    = 0.1;           % oscillation period (s)
Tmax    = 60;            % simulate 60 seconds
dt      = 1e-4;          % 0.1 ms bins
tvec    = 0:dt:Tmax;     % time axis

% instantaneous rate vector
r_t = r0 + A*cos(2*pi*tvec/Tosc);

% thinning method:
% (a) generate candidates at max rate:
rmax = r0 + A;           % 100 Hz
pmax = rmax*dt;
cand = rand(size(tvec)) < pmax;  % candidate spikes

% (b) thin by keeping with prob r(t)/rmax
thin_p = r_t / rmax;
keep   = rand(size(tvec)) < thin_p;
is_spike = cand & keep;

% extract spike times
t_spikes = tvec(is_spike);


% choose lag window ±80 ms
maxLag = 0.16;                  
edges  = -maxLag:dt:maxLag;      
lags   = [];

for i = 1:length(t_spikes)
  d    = t_spikes - t_spikes(i);
  sel  = (d>=-maxLag)&(d<=+maxLag)&(d~=0);
  lags = [lags; d(sel).'];
end

counts  = histcounts(lags, edges);
centers = edges(1:end-1) + dt/2;

figure('Color','w');
bar(1e3*centers, counts, 'k','EdgeColor','none');
xlabel('Lag (ms)'), ylabel('Count')
title('Autocorrelogram: cosine‐modulated Poisson, r(t)=50+50 cos(2\pi t/0.1)')
xlim(1e3*[-maxLag, +maxLag])
