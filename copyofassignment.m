% Author: Awabullah Syed 
% Date: April 28th 2021 
% MPC contoller for Frequency control in a power system
% Description: the primary frequency control loop present in the system
% does not offer adequate control. therefore, the aim is to design a
% secondary freqency control loop that will adjust the reference power,
% delta_p^ref 
%------------------------Tuning MPC--------------------------------------%
x0 = [0;0;0.29];  %choose any value 
Ts = 0.1;
N = 5; % Increasing Horizon, reduces cost value, J 
Q_design = 1000;
% R_design = 10;
% ----------------------------System model paramters---------------------%
M_sys =10; 
D_sys = 0.8;
T_t = 0.5;
T_g = 0.2; 
R_sys = 0.1;
R = 0.25;
%-------------------System Matrices---------------------------------------% 
r = 0.3; %reference input (can be set to any reasonable value)
syms A_sys B_sys C_sys D
A_sys = [-D_sys/M_sys       1/M_sys       0; 
              0             -1/T_t        1/T_t;
         -1/(R_sys*T_g)     0             -1/T_g;]; %A, Dynamic Matrix 
B_sys = [0 ; 0 ; 1/T_g]; % B matrix of the system 
C_sys = [50 0 0]; D = 0;
sys = ss(A_sys,B_sys,C_sys,D_sys); %State space 
sysd = c2d(sys,Ts,'zoh'); %Discrete 
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;
%----------------------------Dimensions----------------------------------%
n = size(A,1);
m = size(B,2);
%--------------------------Cost Matrices---------------------------------%
Q = C'*C*Q_design;
check_ABQR(A,B,Q,R_sys); %if stable without terminal cost 
check_ABQR(A,B,Q,R); %if stable without terminal cost 
% Q = (C'*C)*Q; % Q^(1/2) 

%-----------------------------K -----------------------------------------%
%eig_con = [-1 -1.23 -5.5]; %Desired Eigenvalues for mode 2  
K = -dlqr(A,B,Q,R) %[0 0 0]); % Gain k 
P = dlyap((A+B*K)',Q+K'*R*K); % Lynaponouv (Q & R) - Tuning Parameter)
%--------------------------Constraints------------------------------------%
% input constant
umin = -0.5;   
umax = 0.5; 
Pu = [eye(m) ; -eye(m)];
Qu = [umax;-umin];
% State Constraint 
xmin = -0.5;
xmax = 0.5;
% px = [eye(n); eye(n)];
px = [C;-C];                      %[eye(3).*C; -eye(3).*C];
% qx = [xmax ; xmin];
qx = [xmax ; -xmin];            %[ones(3,1)*xmax; -ones(3,1)*xmin];
%---------------------- mode 1 extend------------------------------------%
Maux = [];  %Extended Mode-1 constraints based on K
for i =0:n-1    %N
    Maux = [Maux;(A+B*K)^(i)];
end 
Mm = kron(eye(n),[px; Pu*K]);       % N
pxf = Mm*Maux; 
qxf = kron(ones(3,1),[qx; Qu]); %qx; 
%-------------------MPC Prediction Matrices-------------------------------%
% prediction matrices
[F,G] = predict_mats(A,B,N);
[H, L, M] = cost_mats(F,G,Q,R,P);   %Cost Matrix 
[Pc,qc,Sc] = constraint_mats(F,G,Pu,Qu,px,qx,pxf,qxf);
%-------------------------Checking Stability----------------------------%
% optimal policy
S = -H\L*x0;
% optimal feedback law
KN = S(1,:);
% form closed-loop system
Phi = A+B*KN;
% stability check
rho = max(abs(eig(Phi)));
if rho >= 1
    display('System with terminal P is not stable')
else
    display('System with terminal P is stable')
end
%------------------Mode 2------------------------------------------------%
% Obtain a plot of x2(k+i|k) versus x1(k+i|k)
XOpt = (F + G*S)*x0;
% Xopt contains x(k+1|k), x(k+2|k), ..., x(k+N|k)
Xopt = [x0; XOpt];

% extract x1 and x2 elements
X1opt = Xopt(1:3:end,:);
X2opt = Xopt(2:3:end,:);
X3opt = Xopt(3:3:end,:);

% Xoptodd = Xopt(1:2:end,:);
% Xopteven = Xopt(2:2:end,:);

% x(k+N|k) is final 2 elements of Xopt
xN = Xopt(end-n+1:end,:);
% x(k+N+1|k) = A*x(k+N|k0) + B*K*x(k+N|k) = (A+B*K)*x(k+N|k)
% x(k+N+2|k) = (A+B*K)*x(k+N+1|k) = (A+B*K)^2 * x(k+N|k)
% x(k+N+i|k) = (A+B*K)^i * x(k+N|k)
% mode-2

N2 = 5;     % Horizon Length

[F2, G2] = predict_mats(A+B*K,B,N2);
% F2 is [(A+BK); (A+BK1)^2; ...; (A+BK)^N2]
% mode 2 state predictions
Xopt2 = [xN; F2*xN;];
% extract x1 and x2 elements
X1opt2 = Xopt2(1:3:end,:);
X2opt2 = Xopt2(2:3:end,:);
X3opt2 = Xopt2(3:3:end,:);

% Xopt2odd = Xopt2(1:2:end,:);
% Xopt2even = Xopt2(2:2:end,:);

figure 
plot3(X1opt,X2opt,X3opt,'bo-',X1opt2,X2opt2,X3opt2,'r*-');
legend({'mode-1','mode-2'}); 
xlabel('x1'); ylabel('x2'); zlabel('x3')
title ('MPC, Mode-1 | Mode-2')

% figure 
% plot(Xoptodd,Xopteven,'bo-'),hold on
% plot(Xopt2odd,Xopt2even,'r*-'); hold off 
% title ('States, Horizon = 3'); 
% xlabel('x1'); 
% ylabel('x2')
% legend({'mode 1','mode 2'},'Location','southeast')
%-----------------Constrained Control law u(k)Kn(x(k))-------------------%
nk = 15; %Simulation steps 
x = x0; 
xss = zeros(n,nk);
uss = zeros(m,nk);
opt = optimoptions('quadprog','ConstraintTolerance',1e-25);

for k = 1:nk+1     %( Addition of ys)
    % store x
    xss(:,k) = x; 
    % solve the MPC problem
    [Uopt, fval, flag] = quadprog(H,L*x,Pc,qc+Sc*x,[],[],[],[],[],opt);    
    if flag < 1
        disp(['Optimization infeasible at k = ' num2str(k)])
    break 
    end 
    % extract the first control
    u = Uopt(1:m); 
    % store u
    uss(:,k) = u; 
    % apply to the plant
    x = A*x + B*u;
    ys(:,k) = C*xss(:,k); 
end
figure(2)  %Output / States Changes 
plot((0:k-1),xss(:,1:k)'); hold on 
plot(0:k-1,repmat(0.01,1,k),'c--')
hold on
plot(0:k-1,repmat(-0.01,1,k),'c--')
hold off
% plot ((0:k-1),xmin);
% plot ((0:k-1),xmax); 
title('Change in State per time step (k)')
xlabel('Time step k'); 
ylabel('States')
legend('\Deltaw','\Deltap^m','\Deltap^v') ;
% legend('Output','xmin constraint','xmax constraint');
hold off

figure(3) % Input 
stairs((0:k-1),uss(:,1:k)'); hold on 
% constraints
plot(0:k-1,repmat(0.5,1,k),'c--')

plot(0:k-1,repmat(-0.5,1,k),'c--')
hold off
% plot((0:k-1),umin,'o')
% plot ((0:k-1),umax,'o')
title ('Input Step Changes')
xlabel('Time step k')
ylabel('Change in Reference Power ');
legend('Input','umax','umin');hold off

figure (4)  % Output 
stairs([0:k-1],ys(:,1:k)')
hold on
% constraints
plot(0:k-1,repmat(0.01,1,k),'c--')
hold on
plot(0:k-1,repmat(-0.01,1,k),'c--')
hold off
xlabel('Time Step (k) ')
ylabel('y (Hz)')
title('Output signal per time step')
%%
PxN = pxf;  %check these 
qxN = qxf;

% build MPC problem matrices
[F, G] = predict_mats(A,B,N);
[H, L, M] = cost_mats(F,G,Q,R,P);
[Pc,qc,Sc] = constraint_mats(F,G,Pu,Qu,px,qx,PxN,qxN);

figure(6)   %feasbility region 
 hold on 
% [x1,x2,x3] = meshgrid(-0.5:0.01:0.5); 
for x1 = -0.01:0.001:0.01
    for x2 = -10:0.1:10
        for x3 = -10:0.1:10
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
% surf(x1,x2,x3,'b.');
% A = [ 1 1 ; 0 1]; %lambda(A+BK) = 0
% B = [0.5;1];
% k = -[1 1.5];
% K = A+B*k
% eig(K)
%% function for checking properties
check_ABQR(A,B,Q,R);

function check_ABQR(A,B,Q,R)

if min(eig(R)) > 0
    display('R is positive definite')
elseif min(eig(R)) == 0
    display('R is positive semidefinite -- DLQR.M will fail!')
else
    display('R is not positive semidefinite  -- DLQR.M will fail!')
end

if min(eig(Q)) > 0
    display('Q is positive definite')
elseif min(eig(Q)) == 0
    display('Q is positive semidefinite')
    % check detectability
    [ABAR,BBAR,CBAR,T,K] = obsvf(A,B,sqrtm(Q));
    nobs = size(A,2) - sum(K);
    if max(abs(eig(ABAR(1:nobs,1:nobs)))) == 0
        display('(Q^0.5,A) is constructible')
    elseif max(abs(eig(ABAR(~K,~K)))) < 1
        display('(Q^0.5,A) is detectable')
    else
        display('(Q^0.5,A) is not detectable -- DLQR.M will fail!')
    end
    
    
else
    display('Q is not positive semidefinite -- DQLR.M will fail!')
end

r = rank(ctrb(A,B));
if r == 2
    display('(A,B) is reachable')
else
    display('(A,B) is not reachable')
    % check stabilizability
    [ABAR,BBAR,CBAR,T,K] = ctrbf(A,B,sqrtm(Q));
    nrea = size(A,2) - sum(K);
    if max(abs(eig(ABAR(1:nrea,1:nrea)))) == 0
        display('(A,B) is controllable')
    elseif max(abs(eig(ABAR(~K,~K)))) < 1
        display('(A,B) is stabilizable')
    else
        display('(A,B) is not stabilizable -- DLQR.M will fail!')
    end
end

end