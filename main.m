%%
clear; clc;
% System Matrices 
A = [7/3 -7/4 0 ; 1 0 0 ;0 0 3/4];
B = [2 ; 0 ; 0]; 
c = 33/100; 
C = [-c 1 1];
D = 0;
N = 4;
Ts = 1;
% Identifying eigenvalues of Q matrices 
syms s 
equ = (0.5776-s)*((1-s)^2-1)+2*0.76*(0.76-0.76*(1-s));
solve(equ);
% Dimensions 
n = size(A,2); % number of states
m = size(B,2); % number of inputs
%  Tuning Parameters 
Q = C'*C; % *Q_design;
% Q = eye(3); 
P = 0*eye(3);  %Q;
R = 10; 

check_ABQR(A,B,Q,R); % Checking the stability of the system
%-----------------------------Unconstrained MPC-------------------------%
% % cost matrices
[F,G] = predict_mats(A,B,N); 
[H,L,M] = cost_mats(F,G,Q,R,P);

% optimal policy
S = -H\L;

% optimal feedback law
KN = S(1,:);

% form closed-loop system
Phi = A+B*KN;

% stability check of the closed loop 
rho = max(abs(eig(Phi)));

if rho >= 1
    display('System with terminal P is not stable')
else
    display('System with terminal P is stable')
end
%%
% Plot to examine the correlation between Horizon Length and tuning
%   parameter R, 
figure (1)
for N=1:25
    for R=1:10000
        [F,G] = predict_mats(A,B,N);
        [H,L,M] = cost_mats(F,G,Q,R,P);
        S = -H\L;
        KN = S(1,:);
        rho = max(abs(eig(A+B*KN)));    % check if the eigenvalue is stable
        if rho <1
            semilogy(N,R,'.r');
            xlim([0 25]);
            hold on
        end
    end
end
xlabel('Horizon Length') 
ylabel('Rho')
title('Horizon Length vs Rho') 
% print('-clipboard','-dmeta')
%%
%-----------------------Unconstrained MPC Response----------------------%
% System Response to set initial condition 
x0 = [10;10;10]; % Intial Condition [10;10;10] or [1;1;1]
nk = 15; % Number of simulation steps
x = x0;
uOpt1 = S*x0;

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
figure (2) % State Respone 
plot(0:nk,xS);
legend('x_1','x_2','x_3')
title('Change in State per time step (k)')
xlabel('Simulation Steps') 
ylabel ('States')
grid on
print('-clipboard','-dmeta')

figure (3) % Input Response 
stairs(0:nk,uS);
title('Change in input')
xlabel('Simulation Steps') 
grid on
print('-clipboard','-dmeta')
%%
% ---------------------------Constrained MPC-----------------------------%
uMax = 1; uMin = -1; % Input Constraint
Px = [eye(n);-eye(n)];
% qx = [xMax xMax -xMin -xMin]';
Pu = [eye(m);-eye(m)];
qu = [uMax -uMin]';

PxN = Px;% kron(eye(n),[Px;Pu*K])*[Phi^0;Phi^1];
% qxN = qx % kron(ones(n,1),[qx;qu]);

% -------Constrained Terminal Equality MPC controller --------------%
% Terminal Equality constraints
% P = 0*Q; PxN = [eye(n);-eye(n)]; qxN = zeros(2*n,1); 

% Terminal equality 
% [F,G] = predict_mats(A,B,N);
% [H,L,M] = cost_mats(F,G,Q,R,P);
% [Pc, qc, Sc] = constraint_mats(F,G,Pu,qu,[],[],PxN,qxN);

%-------------Constrained Terminal Inequality MPC Contoller-------------%
xS = []; % Array to store states
uS = []; % Array to store inputs

x0 = [1;1;1]; % Initial Condition
nk = 30; % Number of simulation steps
x = x0;

% Deadbeat Mode-2
K = -dlqr(A,B,Q,R) % [0 0 0]);  % Deadbeat 
P = dlyap((A+B*K)',Q + K'*R*K);
Phi = A+B*K; % Closed loop transition matrix for Mode-2
rho = max(abs(eig(Phi)));
if rho >= 1
    display('System with terminal P is not stable')
else
    display('System with terminal P is stable')
end
 px = [eye(n); eye(n)];
%  px = [C;-C]; 
 xmin = [-10;-10;-10];
xmax = [10;10;10]
qx = [xmax ; -xmin];   
Maux = [];  %Extended Mode-1 constraints based on K
for i =0:n-1    %N
    Maux = [Maux;(A+B*K)^(i)];
end 
Mm = kron(eye(n),[px; Pu*K]);      
pxf = Mm*Maux; 
qxf = kron(ones(3,1),[qx; qu]);

PxN = kron(eye(n),[Px;Pu*K])*[Phi^0;Phi^1;Phi^2];
qxN = kron(ones(n,1),[qx;qu]);

[F,G] = predict_mats(A,B,N);
[H,L,M] = cost_mats(F,G,Q,R,P);
[Pc, qc, Sc] = constraint_mats(F,G,Pu,qu,px,qx,PxN,qxN);

% PxN = kron(eye(n),[Px;Pu*K])*[Phi^0;Phi^1];
% qxN = kron(ones(n,1),[qx;qu]);
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
%%
% -------------------------Feasibility Region---------------------------%
figure(6)   %feasbility region 
 hold on 
% [x1,x2,x3] = meshgrid(-0.5:0.01:0.5); 
for x1 = -10:0.5:10
    for x2 = -10:0.5:10
        for x3 = -10:0.5:10
        x = [x1; x2;x3]; 
        % solve the MPC problem
        [Uopt, fval, flag] = quadprog(H,L*x,Pc,qc+Sc*x);      
        if flag >= 1
            plot3(x1,x2,x3,'b.')
        end
        end 
        
    end
end
xlabel('x1'); ylabel('x2'); zlabel('x3'); 
title('Operating Region for N = 4');
print('-clipboard','-dmeta')
%%
% ------------------------------Feasibility Region plot----------------%
figure(6)   %feasbility region 
 hold on 
for x1 = -0.01:0.5:0.01
    for x2 = -10:0.5:10
        for x3 = -10:0.5:10
        x = [x1; x2;x3]; 
        % solve the MPC problem
        [Uopt, fval, flag] = quadprog(H,L*x,Pc,qc+Sc*x);      
        if flag >= 1
            plot3(x1,x2,x3,'b.')
        end
        end 
        
    end
end
xlabel('x1'); ylabel('x2'); zlabel('x3'); 
title('Operating Region for N = 5');
%%
%---------------------------States and Input response plot----------------%
% Cleanup arrays
xS = [x0 xS]; % Add initial condition to array
uS = [uS uS(:,end)]; % Replicate final column to show results properly when plotting with stairs

% Plot results
plotRng = 0:(size(xS,2)-1);
figure
plot(plotRng,xS);
title('States')
grid on
print('-clipboard','-dmeta')

figure
stairs(plotRng,uS);
title('Input Response') 
xlabel('Simulation Step')
grid on
print('-clipboard','-dmeta')
%% 
% -------------------------Feasibility Plot------------------------------% 
figure(7)
for x1 = -10:0.5:10
    for x2 = -10:0.5:10
        for x3 = -10:0.5:10
        x = [x1; x2;x3];
        
        % solve the MPC problem
        [Uopt, fval, flag] = quadprog(H,L*x,Pc,qc+Sc*x);
        
        if flag >= 1 % feasible
            plot3(x1,x2,x3,'b.','MarkerSize',5);
            hold on
        end
        end
 
    end
end
title('Region of Feasibility')
axis([-5 5 -5 5]);
grid on
xlabel('x_1'); ylabel('x_2'); zlabel('x_3')
print('-clipboard','-dmeta')

% -------------------------------End of Script--------------------------%