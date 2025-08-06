
%% ------Problem 2d: ---------Recovering RF from L--------

% GETTING RF FROM PROBLEM 2A
%%-------Problem 2a: -------Prepare a simplified model RF of asimple cell in V
clear; clf; rng(0,'twister');
% basic parameters
sig_x=5; % Gaussian width along carrier axis
sig_y=5; % Gaussian width along orthogonal axis
k=0.5; % spatial frequency (cycles/unti)
phi= 0; % phase offset
theta= pi/6; %title angle from vertical

% the grid
x= -25:25; y=-25:25;
[xx,yy]=meshgrid(x,y);

% coordinate rotation: make x' and y' by rotating the coordinates with ø
xp= xx*cos(theta)+ yy*sin(theta);
yp= -xx*sin(theta)+ yy*cos(theta);

% Gabor RF
envelope= exp(-(xp.^2)/(2*sig_x^2)-(yp.^2)/(2*sig_y^2)); % localize the filter in space
norm_factor= 1/(2*pi*sig_y*sig_x); % normalization factor
RF= norm_factor*envelope.*cos(k*xp-phi);

%% GETTING STIM FROM PROBLEM 2B

Nwhite   = 200;
Nblack   = 200;
numTrials = 5000;
STIM = zeros(numel(y),numel(x),numTrials);

for j = 1:numTrials
    frame = zeros(size(xx));               % background = 0
    inds  = randperm(numel(frame),Nwhite+Nblack);
    widx  = inds(1:Nwhite);
    bidx  = inds(Nwhite+1:end);
    frame(widx) = +1;                      % white dots
    frame(bidx) = -1;                      % black dots
    STIM(:,:,j) = frame;
end




%%  From Problem 2c: COMBINE RF AND STIM TO GET LINEAR RESPONSE
L = zeros(numTrials,1);
for j = 1:numTrials
    Sj    = STIM(:,:,j);
    L(j)  = sum( RF(:) .* Sj(:) );        % dot-product
end

%% USING THE RESULTS FOR PROBLEM 2D
% extract top 20% 
pct=prctile(L,90); % 80th percentile of L
sel=find(L>=pct); % indices of trials in the top 20%
% average those top 20%
RF_reconst= mean(STIM(:,:,sel),3); % averaging over the 3rd dimension

% displaying the reconstructed RF


figure('Color','w');
imagesc(x, y, RF_reconst);        % show matrix as image
axis image xy;                 % keep aspect ratio and y‑axis ascending
colormap(jet(256));    % flipped jet gives blue→white→red
c = colorbar;
c.Label.String = 'Reconstructed amplitude';
xlabel('x (deg)');
ylabel('y (deg)');
title(sprintf('Reconstructed RF (top %d%% of STIMs)', numel(sel)/numTrials*100));
