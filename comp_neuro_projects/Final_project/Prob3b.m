%% Problem 3b: Equivalence Between Line y = a x + b and (w1,w2,Θ)
% We show that choosing a line r2 = a*r1 + b is equivalent to selecting
% weights w1, w2 and threshold Θ for the perceptron, and then decide r_out.

clear; close all; clc;

%% 1) Choose line parameters (slope a, intercept b)
% These define the decision boundary r2 = a*r1 + b
a = -1;    
b = 1.5;   

% Convert back to perceptron parameters (up to scale)
w2    = 1;        % we fix w2=1 for simplicity
w1    = -a;       % slope relation: a = -w1/w2 ⇒ w1 = -a*w2 = -a
Theta = b;        % intercept relation: b = Θ/w2 ⇒ Θ = b*w2 = b

fprintf('Derived perceptron params: w1=%.1f, w2=%.1f, Θ=%.1f\n\n', ...
         w1, w2, Theta);

%% 2) Test on the four AND inputs
inputs   = [0 0; 0 1; 1 0; 1 1];   % [r1 r2] pairs
N = size(inputs,1);

% (a) perceptron rule: output = (w1*r1 + w2*r2) > Θ
net_sum   = inputs(:,1)*w1 + inputs(:,2)*w2;
out_perc  = net_sum > Theta;

% (b) line rule: output = (r2 > a*r1 + b)
boundary  = a * inputs(:,1) + b;
out_line  = inputs(:,2)  > boundary;

% AND truth
out_and   = all(inputs,2);

%% 3) Display table of results
disp(' r1  r2 |  sum>Θ  |  r2>a*r1+b  |  AND_expected');
disp('--------------------------------------------');
for i = 1:N
    fprintf('  %d    %d    |   %d     |     %d       |      %d\n', ...
            inputs(i,1), inputs(i,2), ...
            out_perc(i), out_line(i), out_and(i));
end

%% 4) Plot inputs and decision boundary
figure('Name','3b: Decision Boundary','NumberTitle','off');
hold on;
% plot AND inputs: black filled = 1,0; open circles = 0
idx_pos = find(out_and);
idx_neg = find(~out_and);
plot(inputs(idx_neg,1), inputs(idx_neg,2), 'wo','MarkerEdgeColor','k','MarkerSize',10);
plot(inputs(idx_pos,1), inputs(idx_pos,2), 'ks','MarkerFaceColor','k','MarkerSize',10);

% plot boundary line over [0,1] range
xg = linspace(-0.2,1.2,100);
yg = a*xg + b;
plot(xg, yg, 'r--','LineWidth',2);

xlabel('r_1'); ylabel('r_2');
title('Perceptron AND: Decision boundary r_2 = a,r_1 + b');
legend('Class=0','Class=1','Decision Boundary','Location','Best');
xlim([-0.1 1.1]); ylim([-0.1 1.6]);
grid on;
hold off;
