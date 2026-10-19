%% Inverted Pendulum State Space Model
m = 0.3;
l = 0.165;
I = m*l^2/3;
g = 9.81;

A=[ 0 1 0 0;
    0 0 0 0
    0 0 0 1
    0 0 m*g*l/(I+m*l^2) 0];
B=[0;1;0;m*l/(I+m*l^2)];
C=[ 1 0 0 0
    0 0 1 0];
%% Inverted Pendulum State Space Observer
TL = 0.005;

EW =[-1/TL -1/TL/1.01 -1/TL/1.02 -1/TL/1.03];

L = (place(A',C',EW))';
%% Inverted Pendulum State Space Controller
Q = diag([1000 0 1000 0]);

R = 1;

K = lqr(A,B,Q,R);

Cn=[1 0 0 0];

sys_ss = ss(A,B,Cn,0);

Nbar = rscale(sys_ss,K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Container Control

% m = 0.3;
% l = 0.165;
% I = m*l^2/3;
% g = 9.81;
% A=[ 0 1 0 0;
%     0 0 0 0
%     0 0 0 1
%     0 0 -m*g*l/(I+m*l^2) 0];
% B=[0;1;0;-m*l/(I+m*l^2)];
% C=[ 1 0 0 0
%     0 0 1 0];
% % Observer
% TL = 0.005;
% 
% EW =[-1/TL -1/TL/1.01 -1/TL/1.02 -1/TL/1.03];
% 
% L = (place(A',C',EW))';
% % Controller
% Q = diag([1000 0 1000 0]);
% 
% R = 1;
% 
% K = lqr(A,B,Q,R);
% 
% Cn=[1 0 0 0];
% 
% sys_ss = ss(A,B,Cn,0);
% 
% Nbar = rscale(sys_ss,K);