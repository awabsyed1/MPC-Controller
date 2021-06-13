% NOTE that this code has been created to support the revision slides and to aid in your understanding
% It is NOT a one size fits all script and as a result will not work on many problems/can give incorrect answers
% It is simply a tool to help you understand the concepts better
% It is highly recommended that you create your own scripts from scratch to ensure you understand the concepts properly.
%% Problem Setup
clear
clc

% System is reachable
A = [1.5 0;0 0.9];
B = [0.5; 0.6];
n = size(A,2); % 2
m = size(B,2); % 1

% Constraints
xMax = 10; xMin = -10;
uMax = 5; uMin = -5;
Px = [eye(n);-eye(n)];
qx = [xMax xMax -xMin -xMin]';
Pu = [eye(m);-eye(m)];
qu = [uMax -uMin]';

Q = 1e3*eye(2);  % Observable because rank(obsv(A,sqrtm(Q))) = n = 2
R = 1;
N = 5;

% Same as normal constraints (Even if P satisfies Lyapunov there are no guarantees)
% P = 0*Q; PxN = Px; qxN = qx; 

% Terminal Equality constraints
% P = 0*Q; PxN = [eye(n);-eye(n)]; qxN = zeros(2*n,1); 

% Deadbeat Mode-2
K = -acker(A,B,[0 0]); 
P = dlyap((A+B*K)',Q + K'*R*K);
Phi = A+B*K; % Closed loop transition matrix for Mode-2
PxN = kron(eye(n),[Px;Pu*K])*[Phi^0;Phi^1];
qxN = kron(ones(n,1),[qx;qu]);

[F,G] = predict_mats(A,B,N);
[H,L,M] = cost_mats(F,G,Q,R,P);
[Pc, qc, Sc] = constraint_mats(F,G,Pu,qu,Px,qx,PxN,qxN);
%% Simulate
clc
close all

xS = []; % Array to store states
uS = []; % Array to store inputs

x0 = [-1;8]; % Initial Condition
nk = 30; % Number of simulation steps
x = x0;

for k = 1:nk
    % Get control
    [uOpt,fVal,flag] = quadprog(H,L*x,Pc,qc+Sc*x);
    
    if flag < 1 % infeasible
        break; % exit loop
    end
    
    u = uOpt(1:m,:); % Control is the first m rows
    x = A*x + B*u; % Update state
    
    % Store results in last column of respective arrays
    xS(:,k) = x;
    uS(:,k) = u;
end

% Cleanup arrays
xS = [x0 xS]; % Add initial condition to array
uS = [uS uS(:,end)]; % Replicate final column to show results properly when plotting with stairs

% Plot results
plotRng = 0:(size(xS,2)-1);
figure
plot(plotRng,xS);
grid on

figure
stairs(plotRng,uS);
grid on

%% Estimate Feasibility Region (NOTE THIS IF FOR THEOREMS IN L15 and 16) BECAUSE OF RECURSIVE FEASIBILITY GUARANTEES
clc
close all

for x1 = -10:0.5:10
    for x2 = -10:0.5:10
        x = [x1; x2];
        
        % solve the MPC problem
        [Uopt, fval, flag] = quadprog(H,L*x,Pc,qc+Sc*x);
        
        if flag >= 1 % feasible
            plot(x1,x2,'b.','MarkerSize',15);
            hold on
        end
 
    end
end

axis([-10 10 -10 10]);
grid on
