clc
clear
load('Example2_TimeDependent.mat');

%% system dem
n = size(A1, 1);
m = size(B1, 2);
p = size(C1, 1);

%% generate a set of input/state data
T = 15;
X1_0 = [U1(:,1:T);
        U1(:,2:T+1);
        U1(:,3:T+2);
        Y1(1:T);
        Y1(2:T+1);
        Y1(3:T+2)];
X1_1 = [U1(:,2:T+1);
        U1(:,3:T+2);
        U1(:,4:T+3);
        Y1(2:T+1);
        Y1(3:T+2);
        Y1(4:T+3)];
% system 2
X2_0 = [U2(:,1:T);
        U2(:,2:T+1);
        U2(:,3:T+2);
        Y2(1:T);
        Y2(2:T+1);
        Y2(3:T+2)];
X2_1 = [U2(:,2:T+1);
        U2(:,3:T+2);
        U2(:,4:T+3);
        Y2(2:T+1);
        Y2(3:T+2);
        Y2(4:T+3)];     

% system 3
X3_0 = [U3(:,1:T);
        U3(:,2:T+1);
        U3(:,3:T+2);
        Y3(1:T);
        Y3(2:T+1);
        Y3(3:T+2)];
X3_1 = [U3(:,2:T+1);
        U3(:,3:T+2);
        U3(:,4:T+3);
        Y3(2:T+1);
        Y3(3:T+2);
        Y3(4:T+3)]; 
    
%% solve the matrix with dwell time
n_hat = n+n*m;
Delta = 2;
rho = sdpvar(1,1);
nu = sdpvar(1,1);
Q11 = sdpvar(T,n_hat);
Q22 = sdpvar(T,n_hat);
Q33 = sdpvar(T,n_hat);
Q12 = sdpvar(T,n_hat);
Q21 = sdpvar(T,n_hat);
Q12_1 = sdpvar(T,n_hat);
Q21_1 = sdpvar(T,n_hat);
Q13 = sdpvar(T,n_hat);
Q31 = sdpvar(T,n_hat);
Q13_1 = sdpvar(T,n_hat);
Q31_1 = sdpvar(T,n_hat);
Q32 = sdpvar(T,n_hat);
Q23 = sdpvar(T,n_hat);
Q32_1 = sdpvar(T,n_hat);
Q23_1 = sdpvar(T,n_hat);
S1 = sdpvar(m,p);
S2 = sdpvar(m,p);
S3 = sdpvar(m,p);
R1 = sdpvar(p,p);
R2 = sdpvar(p,p);
R3 = sdpvar(p,p);
P1 = sdpvar(n_hat,n_hat, 'symmetric');
P2 = sdpvar(n_hat,n_hat, 'symmetric');
P3 = sdpvar(n_hat,n_hat, 'symmetric');
con1 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con2 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con3 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con4 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con5 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con6 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con7 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con8 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con9 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
M1 = [P1-nu*eye(n_hat), X1_1*Q11; Q11'*X1_1', P1];
M2 = [P2-nu*eye(n_hat), X2_1*Q22; Q22'*X2_1', P2];
M3 = [P3-nu*eye(n_hat), X3_1*Q33; Q33'*X3_1', P3];
M4 = [P1, X1_1*Q12; Q12'*X1_1', P2];
M5 = [P2, X2_1*Q21; Q21'*X2_1', P1];
M6 = [P1, X1_1*Q13; Q13'*X1_1', P3];
M7 = [P3, X3_1*Q31; Q31'*X3_1', P1];
M8 = [P2, X2_1*Q23; Q23'*X2_1', P3];
M9 = [P3, X3_1*Q31; Q31'*X3_1', P2];

constraints = [con1==M1, con1>=0]; 
constraints = [constraints, con2==M2, con2>=0];
constraints = [constraints, con3==M3, con3>=0];
constraints = [constraints, con4==M4, con4>=0];
constraints = [constraints, con5==M5, con5>=0];
constraints = [constraints, con6==M6, con6>=0];
constraints = [constraints, con7==M7, con7>=0];
constraints = [constraints, con8==M8, con8>=0];
constraints = [constraints, con9==M9, con9>=0];
constraints = [constraints, eye(n_hat)>=P1>=rho*eye(n_hat),rho>=0,nu>=0];
constraints = [constraints, eye(n_hat)>=P2>=rho*eye(n_hat)];
constraints = [constraints, eye(n_hat)>=P3>=rho*eye(n_hat)];
constraints = [constraints, P1==X1_0*Q11];
constraints = [constraints, P2==X2_0*Q22];
constraints = [constraints, P3==X3_0*Q33];
constraints = [constraints, X1_1*Q12_1==X1_0*Q12, P2==X1_0*Q12_1];
constraints = [constraints, X2_1*Q21_1==X2_0*Q21, P1==X2_0*Q21_1];
constraints = [constraints, X1_1*Q13_1==X1_0*Q13, P3==X1_0*Q13_1];
constraints = [constraints, X3_1*Q31_1==X3_0*Q31, P1==X3_0*Q31_1];
constraints = [constraints, X2_1*Q23_1==X2_0*Q23, P3==X2_0*Q23_1];
constraints = [constraints, X3_1*Q32_1==X3_0*Q32, P2==X3_0*Q32_1];
constraints = [constraints, R1*Y1(n+1:T+n)==Y1(n+1:T+n)*Q11*X1_0];
constraints = [constraints, R2*Y2(n+1:T+n)==Y2(n+1:T+n)*Q22*X2_0];
constraints = [constraints, R3*Y3(n+1:T+n)==Y3(n+1:T+n)*Q33*X3_0];
constraints = [constraints, S1*Y1(n+1:T+n)==U1(:,(n+1):T+n)*Q11*X1_0];
constraints = [constraints, S2*Y2(n+1:T+n)==U2(:,(n+1):T+n)*Q22*X2_0];
constraints = [constraints, S3*Y3(n+1:T+n)==U3(:,(n+1):T+n)*Q33*X3_0];
objective = -nu;
solutions_out = {S1,R1,S2,R2,S3,R3,nu,rho,con1,con2,con3,con4,P1,P2};
% ops = sdpsettings('verbose',2, 'solver','mosek');
ops = sdpsettings('verbose',2);
controller = optimizer(constraints, objective,ops,[],solutions_out);
[solutions,diagnostics] = controller{[]};
if diagnostics == 1
    disp('The problem is infeasible');
end


S1 = solutions{1};
R1 = solutions{2};
S2 = solutions{3};
R2 = solutions{4};
S3 = solutions{5};
R3 = solutions{6};
K1 = S1*inv(R1);
K2 = S2*inv(R2);
K3 = S3*inv(R3);

%% switching control
x0 = [1;1;1];
N = 41;
x = x0;
u = [];
y = C1*x0;
role = [];
rn = randi(100,1,N);
for i = 1:N
    flag = 0;
    if i <= Delta
        role = [role, 1];
    else
        for j = 2 : Delta
            if role(i - 1) ~= role(i - j)
                role = [role, role(end)];
                flag = 1;
                break;
            end
        end
        if flag == 0
            if mod(i+rn(i),3)==1
                role = [role, 1];
            elseif mod(i+rn(i),3)==2
                role = [role, 2];
            else
                role = [role, 3];
            end
        end
    end
    
    if role(end) == 1
        u = [u K1*y(end)];
        x = [x, (A1+B1*K1*C1)*x(:,end)];
        y = [y, C1*x(:,end)];
    elseif role(end) == 2
        u = [u K2*y(end)];
        x = [x, (A2+B2*K2*C2)*x(:,end)];
        y = [y, C2*x(:,end)];
    else
        u = [u K3*y(end)];
        x = [x, (A3+B3*K3*C3)*x(:,end)];
        y = [y, C3*x(:,end)];        
    end
end

%% plot
figure
subplot(2,1,1)
stairs(0:N-1,role(1:end),'b','LineWidth',2);
axis([0 N-1 0.9 3.1])
ylabel('Switching signal')
grid on
% figure
% plot(1:N,x(1,1:end-1),'r','LineWidth',2);hold on;
% plot(1:N,x(2,1:end-1),'b','LineWidth',2);hold on;
% plot(1:N,x(3,1:end-1),'k','LineWidth',2);
% legend('x1','x2','x3');
% figure
subplot(2,1,2)
plot(0:N-1,y(1:end-1),'b','LineWidth',2);
grid on
ylabel('Output')
xlabel('k')
