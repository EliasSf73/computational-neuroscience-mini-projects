% %  problem 1e. Poisson generator via exponential ISIs
% compare speed and statistics  with problem 1a of bin(Bernoulli) method
clear; clc;
% experiment parameters
N=10;
r=5;
T=0.1;
dt= 1e-3;
n_bins= round(T/dt);
rng(41);

%% ---------- (A) Δt–Bernoulli method (copy of part‑a) -------------------
tic
spk_bin=rand(n_bins,N)< r*dt;
rateA= sum(spk_bin,1)/T;
timeA=toc;

%% ---------- (B) ISI method --------------------------------------------
tic
spk_isi=false(n_bins,N);
for tr=1:N
    t=0;
    ix=1;
    while t<T
        isi= -log(rand)/r;    % exactly equivalent to exprnd(1/r)
        t=t+isi;
        if t>T;break;end
        ix=floor(t/dt)+1;
        spk_isi(ix,tr)=1;
    end
end

rateB= sum(spk_isi,1)/T;
timeB=toc;

%% ---------- statistics & timing ---------------------------------------
fprintf('---  Part e results  --------------------------------\n')
fprintf('Method A (Δt Bernoulli):  μ = %4.2f Hz,  σ = %4.2f Hz,  time = %.2e s\n',...
        mean(rateA),std(rateA),timeA);
fprintf('Method B (ISI draw)    :  μ = %4.2f Hz,  σ = %4.2f Hz,  time = %.2e s\n',...
        mean(rateB),std(rateB),timeB);


%% optional visual check: raster plots ----------------------------------
figure('Color','w');
subplot(2,1,1)
imagesc(spk_bin'); colormap(gray); axis xy
title('Δt–Bernoulli (part a)');
ylabel('trial'), xticklabels([]);

subplot(2,1,2)
imagesc(spk_isi'); colormap(gray); axis xy
title('ISI‑based (part e)');
xlabel('time bin (1 ms)'), ylabel('trial');