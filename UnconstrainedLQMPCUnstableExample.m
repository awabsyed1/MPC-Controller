% NOTE that this code has been created to support the revision slides and to aid in your understanding
% It is NOT a one size fits all script and as a result will not work on many problems/can give incorrect answers
% It is simply a tool to help you understand the concepts better
% It is highly recommended that you create your own scripts from scratch to ensure you understand the concepts properly.
%% Problem Setup
clear
clc

% System is reachable
A = [2.5 0;0 0.9];
B = [1;2];
n = size(A,2); % number of states
m = size(B,2); % number of inputs

Q = eye(2); % (Q^1/2,A) is observable
R = 1;
P = Q; % we don't yet have a proper way to choose this
N = 2;

[F,G] = predict_mats(A,B,N);
[H,L,M] = cost_mats(F,G,Q,R,P);

x0 = [2;2]; % Intial Condition

S = -H\L;
uOpt1 = S*x0;
uOpt2 = quadprog(H,L*x0);
[uOpt1 uOpt2] % display them side by side. SAME RESULT

%% To find Unconstrained LQ-MPC controller
clc
KN = S(1:m,:); % First m rows of S

max(abs(eig(A+B*KN))) % check system stability (should be less than 1 for a stable system)
%% Simulate system
clc
% close all

xS = []; % Array to store states
uS = []; % Array to store inputs

x0 = [10;10]; % Intial Condition
nk = 15; % Number of simulation steps
x = x0;

for k = 1:nk
    % Get control
    u = KN*x;
    
    x = A*x + B*u; % Update state
    
    % Store results in last column of respective arrays
    xS(:,k) = x;
    uS(:,k) = u;
end

% Cleanup arrays
xS = [x0 xS]; % Add initial condition to array
uS = [uS uS(:,end)]; % Replicate final column to show results properly when plotting with stairs

% Plot results
figure
plot(0:nk,xS);
grid on

figure
stairs(0:nk,uS);
grid on









