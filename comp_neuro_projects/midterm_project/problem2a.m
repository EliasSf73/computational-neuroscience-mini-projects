%%-------Problem 2a: -------Prepare a simplified model RF of asimple cell in V

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

% plot

figure('Color','w');
imagesc(x,y,RF);
axis image xy;           % equal scaling, origin bottom‑left
colormap(jet); 
colorbar;
title('Gabor RF: \sigma_x=5,\sigma_y=5, k=0.4, \theta=30°');
xlabel('x (deg)'), ylabel('y (deg)');