%% Problem 1c: Synaptic Discrimination of Bursting vs. Tonic Input
% Goal: find a synaptic weight w such that
%  - bursting input (from 1b) triggers a postsynaptic spike
%  - fast-spiking input (from 1a) does NOT trigger a spike
% We model: 
%   * EPSC:    dI_syn/dt = -I_syn/tau_syn  + w·sum(delta(t - t_pre))
%   * LIF:     dV/dt     = -(V-V_rest)/tau_m + I_syn/C
% Euler, dt=1 ms, T=500 ms

clear; close all;

%% 1) Common sim parameters
dt      = 1;           % ms
T       = 500;         % ms
t       = 0:dt:T;      
N       = numel(t);
I_drive = 10;          % mA (for presyn model)

%% 2) Generate presynaptic spike trains
% 2a) Fast-spiking train (a=0.05,b=0.15,c=-75,d=2)
a1=0.05; b1=0.15; c1=-75; d1=2;
V1 = -65*ones(1,N); u1 = b1*V1; spikes_fs = [];
for i=2:N
    dV = 0.04*V1(i-1)^2 + 5*V1(i-1) + 140 - u1(i-1) + I_drive;
    du = a1*(b1*V1(i-1) - u1(i-1));
    V1(i)=V1(i-1)+dt*dV;  u1(i)=u1(i-1)+dt*du;
    if V1(i)>=30
      spikes_fs(end+1)=t(i);
      V1(i)=c1;  u1(i)=u1(i)+d1;
    end
end

% 2b) Bursting train  (a=0.02,b=0.15,c=-50,d=2)
a2=0.02; b2=0.15; c2=-50; d2=2;
V2 = -65*ones(1,N); u2 = b2*V2; spikes_bu = [];
for i=2:N
    dV = 0.04*V2(i-1)^2 + 5*V2(i-1) + 140 - u2(i-1) + I_drive;
    du = a2*(b2*V2(i-1) - u2(i-1));
    V2(i)=V2(i-1)+dt*dV;  u2(i)=u2(i-1)+dt*du;
    if V2(i)>=30
      spikes_bu(end+1)=t(i);
      V2(i)=c2;  u2(i)=u2(i)+d2;
    end
end

%% 3) Postsynaptic LIF + synapse params
V_rest   = -65;        % mV
V_thresh = -40;        % mV
V_reset  = V_rest;     % mV
tau_m    = 20;         % ms
C         = 1;         % arb. units
tau_syn  = 10;         % ms
w_list   = linspace(0.1,5,100);  % search range for w

best_w = NaN;
% search synaptic weight w
for w = w_list
  % simulate FS input
  [spk_fs_post, ~] = run_LIF_with_syn(w, tau_syn, spikes_fs);
  % simulate bursting input
  [spk_bu_post, ~] = run_LIF_with_syn(w, tau_syn, spikes_bu);
  
  if isempty(spk_fs_post) && ~isempty(spk_bu_post)
    best_w = w;
    break
  end
end
fprintf('Chosen synaptic weight: w = %.2f\n', best_w);

%% 4) Confirm & plot both cases
[spk_fs_post, V_fs_post] = run_LIF_with_syn(best_w, tau_syn, spikes_fs);
[spk_bu_post, V_bu_post] = run_LIF_with_syn(best_w, tau_syn, spikes_bu);

figure('Name','Pre vs Post: Fast-Spiking','NumberTitle','off');
subplot(2,1,1);
stem(spikes_fs, ones(size(spikes_fs))*V_thresh+5, 'k','Marker','none');
hold on; plot(t, V_fs_post, 'b','LineWidth',1);
hold off;
xlabel('Time (ms)'); ylabel('mV');
title('Fast-Spiking Input → No Postsynaptic Spike');

subplot(2,1,2);
stem(spikes_bu, ones(size(spikes_bu))*V_thresh+5, 'k','Marker','none');
hold on; plot(t, V_bu_post, 'r','LineWidth',1);
hold off;
xlabel('Time (ms)'); ylabel('mV');
title('Bursting Input → Postsynaptic Spike');

%% 5) Summary stats: postsynaptic firing rates
fr_fs_post = numel(spk_fs_post) / (T/1000);    % Hz, should be 0
fr_bu_post = numel(spk_bu_post) / (T/1000);    % Hz, >0

fprintf('Postsynaptic FR (fast):    %.1f Hz\n', fr_fs_post);
fprintf('Postsynaptic FR (bursting): %.1f Hz\n', fr_bu_post);

%% 6) Bar plot for the report
figure('Name','Postsynaptic Firing Rate Comparison','NumberTitle','off');
bar([fr_fs_post, fr_bu_post], 0.5);
set(gca, 'XTickLabel', {'Fast-Spiking','Bursting'}, 'FontSize',12);
ylabel('Firing rate (Hz)', 'FontSize',12);
title('Postsynaptic Firing Rate', 'FontSize',14);
ylim([0, max([fr_bu_post,1])*1.2]);   % give some headroom



%% --- helper function ---
function [spike_times_post, V] = run_LIF_with_syn(w, tau_syn, pre_spikes)
  % Simulate LIF neuron with exponential synapse (dt=1ms, T=500ms)
  dt = 1; T = 500; t = 0:dt:T; N = numel(t);
  V = -65 * ones(1,N);
  I_syn = zeros(1,N);
  spike_times_post = [];
  next_pre = 1;
  
  for i = 2:N
    % update I_syn
    I_syn(i) = I_syn(i-1) - (dt/tau_syn)*I_syn(i-1);
    if next_pre <= numel(pre_spikes) && t(i)==pre_spikes(next_pre)
      I_syn(i) = I_syn(i) + w;
      next_pre = next_pre + 1;
    end
    
    % update V
    dV = -(V(i-1)-(-65))/20 + I_syn(i)/1;
    V(i) = V(i-1) + dt*dV;
    
    % check spike
    if V(i) >= -40
      spike_times_post(end+1) = t(i); %#ok<AGROW>
      V(i) = -65;  
      
      % reset
    end
  end
end
