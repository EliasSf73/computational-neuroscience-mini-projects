%% Problem 3d: Error Surface E(a,b) Heatmap + Extensions
% Sweep (a,b) over a grid, compute the AND‐misclassification cost E(a,b),
% plot the heatmap, compute the fraction of perfect‐AND lines, and show
% perceptron‐learning dynamics in (a,b)-space.

clear; close all; clc;

%% 1) Define AND inputs and target outputs
inputs  = [0 0;  0 1;  1 0;  1 1];  % [r1,r2] rows
desired = all(inputs,2);           % AND truth: [0;0;0;1]

%% 2) Grid sweep for E(a,b)
a_min = -3; a_max =  3;   % slope range
b_min = -1; b_max =  3;   % intercept range
Ngrid = 200;              % resolution

a_vals = linspace(a_min,a_max,Ngrid);
b_vals = linspace(b_min,b_max,Ngrid);
E      = zeros(Ngrid,Ngrid);

for ia = 1:Ngrid
  a = a_vals(ia);
  for ib = 1:Ngrid
    b = b_vals(ib);
    preds    = inputs(:,2) > (a*inputs(:,1) + b);
    E(ib,ia)= sum(abs(preds - desired));
  end
end

%% 3) Plot heatmap of misclassification cost
figure('Name','3d: Error Heatmap','NumberTitle','off');
imagesc(a_vals, b_vals, E);
set(gca,'YDir','normal');
colormap(flipud(hot)); 
colorbar('Ticks',0:4,'TickLabels',0:4);
xlabel('slope a'); ylabel('intercept b');
title('Misclassification Error E(a,b) for AND Gate');
hold on;
% overlay zero‐error contour
contour(a_vals,b_vals,E,[0 0],'c-','LineWidth',2);
hold off;

%% 4) Extension 1: Fraction of (a,b)-space yielding perfect AND
frac_AND = sum(E(:)==0)/numel(E);
fprintf('Fraction of grid with E=0: %.4f (%.2f%%)\n', frac_AND, frac_AND*100);

% Annotate on the heatmap
hold on;
text(a_min+0.2, b_max-0.2, ...
     sprintf('Zero‐error region = %.2f%%', frac_AND*100), ...
     'Color','w','FontSize',12,'FontWeight','bold');
hold off;

%% 5) Extension 2: Perceptron‐Learning Dynamics in (a,b)-space
% We include bias Θ and learn on the four AND examples.

eta   = 0.1;        % learning rate
w     = [0;0];      % initial weights [w1;w2]
Theta = 0;          % initial threshold
traj  = [];         % to store (a,b) trajectory

for epoch = 1:50
  for i = 1:4
    x = inputs(i,:)';     % column vector [r1;r2]
    y = desired(i);       % target 0 or 1
    pred = (w'*x > Theta);
    delta = eta*(y - pred);
    w     = w + delta*x;       % update weights
    Theta = Theta - delta;     % update threshold
  end
  % record equivalent (a,b)
  a_learn = -w(1)/w(2);
  b_learn =  Theta/w(2);
  traj(end+1,:) = [a_learn, b_learn];
end

% Overlay learning trajectory on the heatmap
figure('Name','3d: Learning Trajectory','NumberTitle','off');
imagesc(a_vals, b_vals, E);
set(gca,'YDir','normal');
colormap(flipud(hot));
hold on;
plot(traj(:,1), traj(:,2), 'w.-','LineWidth',1.5,'MarkerSize',10);
plot(traj(1,1),traj(1,2),'go','MarkerSize',8,'MarkerFaceColor','g');    % start
plot(traj(end,1),traj(end,2),'ks','MarkerSize',8,'MarkerFaceColor','k');% end
colorbar('Ticks',0:4,'TickLabels',0:4);
xlabel('slope a'); ylabel('intercept b');
title('Perceptron Learning Trajectory in (a,b)-space');
legend('E(a,b)','Learning path','Start','Converged','Location','best');
hold off;
