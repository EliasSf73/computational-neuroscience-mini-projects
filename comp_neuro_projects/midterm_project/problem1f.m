% Parameters
r      = 50;            % target rate (Hz)
T      = 100;           % total duration (s)
maxLag = 0.08;          % +/- 80 ms window
dtLag  = 1e-3;          % bin width for the correlogram
%% --------------------Problem 1f: autocorrelation of spikes with constant firing rate-----------

% 1. Draw ISI from exponential of 1/r rate
spikes=[];
t=0;
while t<T
    isi=-log(rand)/r;  % exactly equivalent to exprnd(1/r)
    t=t+isi;
    if t<T
        spikes(end+1)=t;
    end
end


% 2) Compute all pairwise non‑zero lags within ±maxLag
lags=[]; % array to collect inter-spike differences
for i=1:numel(spikes)
    d= spikes-spikes(i); % Subtract the reference time spikes(i) from every spike time in spikes
    sel= (d>= -maxLag)& (d<= maxLag) & (d~=0);% sel true if d in [-maxlag,maxlag] but !=0
     lags= [lags; d(sel).'];   
end

% 3) Histogram those lags
edges   = -maxLag : dtLag : +maxLag;
counts  = histcounts(lags, edges);
centers = edges(1:end-1) + dtLag/2;

% 4) Plot
figure('Color','w');
bar(centers*1e3, counts, 'k');
xlabel('Lag (ms)');
ylabel('Count');
title('Autocorrelogram (ISI‑based Poisson, r=50Hz)');
xlim(1e3*maxLag*[-1 1]);
grid on;

