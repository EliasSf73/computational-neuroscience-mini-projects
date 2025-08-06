
%% ------Problem 2c: calculating Linear Response--------

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

Nwhite   = 30;
Nblack   = 30;
numTrials = 200;
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




%%  COMBINE RF AND STIM TO GET LINEAR RESPONSE
L = zeros(numTrials,1);
for j = 1:numTrials
    Sj    = STIM(:,:,j);
    L(j)  = sum( RF(:) .* Sj(:) );        % dot-product
end

figure('Name','L histogram','Color','w');
histogram(L,20,'Normalization','probability');
xlabel('Linear response L_j');
ylabel('P(L)');
title('Histogram of linear responses (200 trials)');
grid on;
