%% Problem 1b: Find Bursting Parameters @ 50 Hz + Extension
% (1) Grid‐search for (a,b,c,d) producing bursting at ≈50 Hz, I=10 mA
% (2) Confirm with V(t) & ISI histogram
% (3) extension: F–B curve (bursts/sec vs. I)

clear; close all;

%% 1) Simulation settings
dt = 1;            % ms
T  = 500;          % total time (ms)
t  = 0:dt:T;       
N  = numel(t);
I0 = 10;           % drive for search (mA)

%% 2) Parameter grid (bursting regime around Izhikevich’s values)
a_list = [0.02, 0.03, 0.04];
b_list = [0.15, 0.20, 0.25];
c_list = [-60, -55, -50];
d_list = [2, 4, 6, 8];

best_cv   = -inf;
best_err  = inf;
best_par  = [];

% Loop through every candidate
for a = a_list
  for b = b_list
    for c = c_list
      for d = d_list
        %--- simulate one run ---
        V = -65*ones(1,N);
        u = b*V;
        spikes = [];
        for k = 2:N
          dV = 0.04*V(k-1)^2 + 5*V(k-1) + 140 - u(k-1) + I0;
          du = a*(b*V(k-1) - u(k-1));
          V(k) = V(k-1) + dt*dV;
          u(k) = u(k-1) + dt*du;
          if V(k) >= 30
            spikes(end+1) = t(k);    %#ok<AGROW>
            V(k) = c;                
            u(k) = u(k) + d;         
          end
        end
        
        %--- metrics: mean rate & ISI‐CV ---
        f_avg = numel(spikes)/(T/1000);       % Hz
        ISIs  = diff(spikes);
        if isempty(ISIs), continue; end
        cv_isi = std(ISIs)/mean(ISIs);        % bursting high CV
        
        %--- pick best: close to 50 Hz & max CV ---
        err = abs(f_avg - 50);
        if (err < 5) && (cv_isi > best_cv)   % allow ±5 Hz window
          best_cv  = cv_isi;
          best_err = err;
          best_par = [a, b, c, d, f_avg, cv_isi];
        end
      end
    end
  end
end

% Unpack best
a = best_par(1);  b = best_par(2);  c = best_par(3);  d = best_par(4);
f_avg = best_par(5);  cv_isi = best_par(6);
fprintf('Found bursting params: a=%.2f, b=%.2f, c=%d, d=%d → f≈%.1f Hz, CV≈%.2f\n', ...
        a,b,c,d,f_avg,cv_isi);

%% 3) Confirm bursting with best set
V = -65*ones(1,N);  
u = b*V;            
spikes = [];
for k = 2:N
  dV = 0.04*V(k-1)^2 + 5*V(k-1) + 140 - u(k-1) + I0;
  du = a*(b*V(k-1) - u(k-1));
  V(k) = V(k-1) + dt*dV;
  u(k) = u(k-1) + dt*du;
  if V(k) >= 30
    spikes(end+1) = t(k);
    V(k) = c;
    u(k) = u(k) + d;
  end
end
ISIs = diff(spikes);

figure('Name','Bursting (Confirmed)','NumberTitle','off');
subplot(2,1,1);
plot(t, V, 'b', 'LineWidth',1);
xlabel('Time (ms)'); ylabel('V (mV)');
title(sprintf('Intrinsic Bursting (f≈%.1f Hz, CV≈%.2f)', f_avg, cv_isi));

subplot(2,1,2);
histogram(ISIs,20);
xlabel('ISI (ms)'); ylabel('Count');
title('ISI Histogram (Bursting)');

%% 4) extension: F–B curve (burst rate vs. I)
I_list = 8:0.5:12;  
burst_thr = 15;      % ms to separate bursts
burst_rate = zeros(size(I_list));

for idx = 1:numel(I_list)
  I = I_list(idx);
  
  % simulate with best (a,b,c,d)
  V = -65*ones(1,N);
  u = b*V;
  sp = [];
  for k = 2:N
    dV = 0.04*V(k-1)^2 + 5*V(k-1) + 140 - u(k-1) + I;
    du = a*(b*V(k-1) - u(k-1));
    V(k) = V(k-1) + dt*dV;
    u(k) = u(k-1) + dt*du;
    if V(k) >= 30
      sp(end+1) = t(k);
      V(k) = c;
      u(k) = u(k) + d;
    end
  end
  
  % count bursts
  ISIsp = diff(sp);
  burst_starts = [1, find(ISIsp > burst_thr)+1];
  burst_rate(idx) = numel(burst_starts)/(T/1000);
end

figure('Name','Burst Rate vs. I','NumberTitle','off');
plot(I_list, burst_rate, '-o', 'LineWidth',1.2);
xlabel('Input current I (mA)');
ylabel('Burst rate (bursts/s)');
title('F–B Curve: Bursting Frequency vs. I');
grid on;
