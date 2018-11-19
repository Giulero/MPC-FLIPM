%% Take here the matrices for the new mpc problem

clear, clear all, clc

% syms  b k omega delta m1 m2 tau real

b = 1;
k= 10;
omega = 6;
delta = 0.02;
m1 = 5;
m2 = 1;
tau = 10;


M = (m1+m2)/(m1*m2);

A_2 = [   0    1;
       -k*M -b*M];
   
B_2 = [0;
       1/m2];

alpha = 1/(m1*m2)*(tau*(k - M*b^2) + b);
beta = k/(m1*m2)*(1 - M*b*tau);
   
C_a = [beta*m1 alpha*m2];
   
D_a = beta*tau/(m1*m2);

A = [0       1       0 0 0;
     omega^2 0 omega^2 0 0;
     0       0       0 0 0;
     0       0       0 0 0;
     0       0       0 0 0];
 
A(4:5,4:5) = A_2 - B_2/D_a*C_a;
A(4:5,1) = B_2/D_a*omega^2;
A(4:5,3) = B_2/D_a*omega^2;


B = [ 0 0 1 0 0]';

C = [0 0 0 1 0];

D = 0;

sys = ss(A,B,C,D);
d_sys = c2d(sys, delta);
A = d_sys.A;
B = d_sys.B;
C = d_sys.C;
D = d_sys.D;

N = 10;

MPC_matrix = zeros(N);
C_MPC_matrix = zeros(N,5);


Q = eye(5);
R = 1;
Q_t = eye(5*N);
R_t = eye(N);
S_t = zeros(5*N, N);
for i=1:N
    for j=1:i
        S_t((5*(i-1)+1):(5*(i-1)+5),j) = A^(i-j)*B;
    end
end
T_t = zeros(5*N, 5);
for i=1:N
    T_t((5*(i-1)+1):(5*(i-1)+5),:) = A^i
end
H = 2*(R_t+S_t'*Q_t*S_t);
F = 2*(T_t'*Q_t*S_t);
Y = 2*(Q+T_t'*Q_t*T_t);
size(H);
issymmetric(H)
[~,p] = chol(H)
