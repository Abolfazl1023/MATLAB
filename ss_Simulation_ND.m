% Sensor Noise & Disturbance Compensation Using Kalman Filter (KF) Gain & LQR Controller
% Author: Abolfazl Roostaeiparidari
% Email: roostaeiabolfazl@gmail.com
% Thanks to Steven L. Brunton - University of Washington

% State Matrices
Ix=1000;
A=[0 1
   0 0];
B=[0
    1/Ix];
C=[1 0];
D=zeros(size(C,1),size(B,2));
% Parameter Q
q1=Ix^2;
q2=Ix^2;
Q1=[q1 0
   0 q2];
dt=0.01; t=dt:dt:50;
R=0.001;
k1=lqr(A,B,Q1,R);

% System State space with Noise+Disturbance
Vd=0.02*eye(2); %Disturbance
Vn=0.2;
BF=[B Vd 0*B];
sysC=ss(A,BF,C,[0 0 0 Vn]);

% System True 
sysFullOutput=ss(A,BF,eye(2),zeros(2,size(BF,2)));

% Kalman filter design
% [KF, P, E]=lqe(A,Vd,C,Vd,Vn);
Kf=lqr(A',C',Vd,Vn);
KF=Kf';
sysKF=ss(A-KF*C,[B KF],eye(2),0*[B KF]);
%%
u=0*t;
u(100:120)=100;
u(1500:1520)=-100;

uD=0.1*randn(2,size(t,2));
uN=0.1*randn(size(t));
uAUG=[u; Vd*Vd*uD; uN];

[y,t]=lsim(sysC,uAUG,t);
plot(t,y,'b');
hold on

[Xtrue,t]=lsim(sysFullOutput,uAUG,t);
plot(t,Xtrue(:,1),'r','LineWidth',2);
hold on

[x,t]=lsim(sysKF,[u; y'],t);
plot(t,x(:,1),'k--','LineWidth',2);
title('Impulse Response')
legend('Faulty/Disturbed System','Correct System','Estimated System')
hold off




