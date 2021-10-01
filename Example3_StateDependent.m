clc
clear
load('Example3_StateDependent.mat');

%% system dem
n = size(A1, 1);
m = size(B1, 2);
p = size(C1, 1);

%% generate a set of input/output data
T = 10;
% % system 1
X1_0 = [U1(1:T);
        U1(2:T+1);
        Y1(1:T);
        Y1(2:T+1)];
X1_1 = [U1(2:T+1);
        U1(3:T+2);
        Y1(2:T+1);
        Y1(3:T+2)];
% system 2
X2_0 = [U2(1:T);
        U2(2:T+1);
        Y2(1:T);
        Y2(2:T+1)];
X2_1 = [U2(2:T+1);
        U2(3:T+2);
        Y2(2:T+1);
        Y2(3:T+2)];              

%% design controller
n_hat = 2*n;


theta1 = 1;
theta2 = 1;
mu1 = 0.1;
mu2 = 0.1;
tau1 = 0.1;
tau2 = 0.1;
iota1 = 0.1;
iota2 = 0.1;

rho = sdpvar(1,1);
nu = sdpvar(1,1);
P1 = sdpvar(n_hat,n_hat,'symmetric');
P1_1 = sdpvar(n_hat,n_hat,'symmetric');
W1 = sdpvar(n_hat,n_hat,'symmetric');
Q1 = sdpvar(T,n_hat);
Q12 = sdpvar(T,n_hat);
Q21 = sdpvar(T,n_hat);
R1 = sdpvar(p,p);
S1 = sdpvar(m,p);

P2 = sdpvar(n_hat,n_hat,'symmetric');
P2_1 = sdpvar(n_hat,n_hat,'symmetric');
W2 = sdpvar(n_hat,n_hat,'symmetric');
Q2 = sdpvar(T,n_hat);
R2 = sdpvar(p,p);
S2 = sdpvar(m,p);
W3 = sdpvar(n_hat,n_hat,'symmetric');

con1 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con2 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con3 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
con4 = sdpvar(2*n_hat, 2*n_hat, 'symmetric');
M1 = [X1_0*Q1-mu1*W1-nu*eye(n_hat), X1_1*Q1; Q1'*X1_1', X1_0*Q1];
M2 = [X2_0*Q2-mu2*W2-nu*eye(n_hat), X2_1*Q2; Q2'*X2_1', X2_0*Q2];
M3 = [-X1_0*Q1+iota1*W1, X1_1*Q12; Q12'*X1_1', X2_0*Q2-Q12'*X1_0'-X1_0*Q12+iota1*W2];
M4 = [-X2_0*Q2+iota2*W2, X2_1*Q21; Q21'*X2_1', X1_0*Q1-Q21'*X2_0'-X2_0*Q21+iota2*W1];

constraints = [con1==M1, con1>=0, rho>=0.0000001, nu>=0.0000001];
constraints = [constraints, con2==M2, con2>=0];
constraints = [constraints, con3==M3, con3<=0];
constraints = [constraints, con4==M4, con4<=0];
constraints = [constraints, X1_0*Q1-tau1*W1-rho*eye(n_hat)==P1, P1>=0];
constraints = [constraints, eye(n_hat)-X1_0*Q1-tau1*W1==P1_1, P1_1>=0];
constraints = [constraints, X2_0*Q2-tau2*W2-rho*eye(n_hat)==P2, P2>=0];
constraints = [constraints, eye(n_hat)-X2_0*Q2-tau2*W2==P2_1, P2_1>=0];
constraints = [constraints, R1*Y1(:,n+1:T+n)==Y1(:,n+1:T+n)*Q1*X1_0];
constraints = [constraints, R2*Y2(:,n+1:T+n)==Y2(:,n+1:T+n)*Q2*X2_0];
constraints = [constraints, S1*Y1(:,n+1:T+n)==U1(:,n+1:T+n)*Q1*X1_0];
constraints = [constraints, S2*Y2(:,n+1:T+n)==U2(:,n+1:T+n)*Q2*X2_0];
constraints = [constraints, W3 == theta1*W1+theta2*W2, W3>=0];

objective = -nu;
solutions_out = {S1,R1,S2,R2,W1,W2,nu,rho,con1,con2,con3,con4};
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
W1 = solutions{5};
W2 = solutions{6};
K1 = S1*inv(R1);
K2 = S2*inv(R2);

%% switching control
x0 = [1;-1];
N = 21;
u = [];
x = x0;
y = C1*x0;
role = [];
w = [];
for i = 1:N
    if i <= n
        role = [role, 2];
        u = [u K2*y(end)];
        x = [x, A2*x(:,end)+B2*u(end)];
        y = [y, C2*x(:,end)];
    else
        tmp = [u(end-1:end),y(end-2:end-1)];
        if tmp*W1*tmp'> tmp*W2*tmp'
            role = [role, 1];
            u = [u K1*y(end)];
            x = [x, (A1+B1*K1*C1)*x(:,end)];
            y = [y, C1*x(:,end)];
            w = [w,tmp*W1*tmp'];
        else 
            role = [role, 2];
            u = [u K2*y(end)];
            x = [x, (A2+B2*K2*C2)*x(:,end)];
            y = [y, C2*x(:,end)];
            w = [w,tmp*W2*tmp'];
        end        
    end

end

%% plot
figure
subplot(2,1,1)
stairs(0:N-1,role(1:end),'b','LineWidth',2);
axis([0 N-1 0.9 2.1])
ylabel('Switching signal')
grid on
% plot(0:N-1,x(1,1:end),'r','LineWidth',2);hold on;
% plot(0:N-1,x(2,1:end),'b','LineWidth',2);
% legend('x1','x2');
% figure
subplot(2,1,2)
plot(0:N-1,y(1:end-1),'b','LineWidth',2);
grid on
ylabel('Output')
xlabel('k')