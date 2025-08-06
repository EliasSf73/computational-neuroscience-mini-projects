%%  poisson spike-train under time-varying rate r(t)
clear;clc;
% fundamental parameters
N=5; % 5 trials
dt=1e-3;
t= 0:dt:1;
% sinusoidal parameters f
c1= 20;
c2= 4*pi;
c3=pi/2;
c4=30;

nbins= round(1/dt); % number of bins for each trial
% derived parameters
rng(41);

for i in t
    r=c1*sin(c2*i-c3)+c4;  % r(t)
    p= r*dt;
    spk= rand(nbins,N)<p; % spike matrix
% bins on rows and trials on columns then compare each value to p
% spk will have 0s(no spike) and 1s(spike) as entries

% firing rate in in each trial (column sum)
num_spk= sum(spk,1);
rate=num_spk/T; % spike rate in Hz
