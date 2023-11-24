Mp=0.127;
Lp=0.33655;
Jp=0.001198;
Jr=9.98E-04;
Lr=0.2159;
g=9.81;
Dr=0.0024;
Dp= 0.0024;
Mr = 0.257;
Kg=70;
kt= 0.007682;
km= 0.007677;
Rm=2.6;
zeta=0.7;
wn=4;

% State Space Representation
Jt = Jr*Jp + Mp*(Lp/2)^2*Jr + Jp*Mp*Lr^2;
A = [0 0 1 0;
     0 0 0 1;
     0 Mp^2*(Lp/2)^2*Lr*g/Jt -Dr*(Jp+Mp*(Lp/2)^2)/Jt -Mp*(Lp/2)*Lr*Dp/Jt;
     0 Mp*g*(Lp/2)*(Jr+Mp*Lr^2)/Jt -Mp*(Lp/2)*Lr*Dr/Jt -Dp*(Jr+Mp*Lr^2)/Jt];

B = [0; 0; (Jp+Mp*(Lp/2)^2)/Jt; Mp*(Lp/2)*Lr/Jt];
C = eye(2,4);
D = zeros(2,1);

% Add actuator dynamics


B = Kg * kt * B / Rm;
L=-Dr*(Jp+Mp*(Lp/2)^2)/Jt;
N=(Kg^2*kt*km/Rm)*B(3);
A(3,3)=L-N;

R=-Mp*(Lp/2)*Lr*Dr/Jt;
T= (Kg^2*kt*km/Rm)*B(4);
A(4,3) = R-T;



% Load into state-space system

rp_sys = ss(A,B,C,D);

pole=eig(A)


co = ctrb(rp_sys);
controllability = rank(co)


% Find desired poles
% Location of dominant poles along real-axis
sigma = zeta*wn;
% Location of dominant poles along img axis (damped natural freqency)
wd = wn*sqrt(1-zeta^2);
% Desired poles (-30 and -40 are given)
DP = [-sigma+1i*wd, -sigma-1i*wd, -30, -40];    
%
%% Find Tranformation Matrix W
% Characteristic equation: s^4 + a_4*s^3 + a_3*s^2 + a_2*s + a_1
% inverted pendulum characteristic equation: s^4 + 11.6*s^3 - 117.3*s^2 - 408.3*s
a = poly(A);
% 
% Companion matrices (Ac, Bc)
Ac = [  0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        -a(5) -a(4) -a(3) -a(2)];
%
Bc = [0; 0; 0; 1];
%
% Controllability
T = ctrb(A,B)
% Controllability of companion matrices
Tc = ctrb(Ac,Bc)
%
% Transformation matrices
W = T*inv(Tc)
%
%% Find Gain
% Desired polynomial from desired poles in vector "DP" (found in
% setup_rotpen.m)
a_des = poly(DP);
% Companion control gains
k1c = a_des(5) - a(5);
k2c = a_des(4) - a(4);
k3c = a_des(3) - a(3);
k4c = a_des(2) - a(2);
% Companion state-feedback control gain
Kc = [k1c k2c k3c k4c]; 
% The students would enter the gain manually from their pre-lab as
% Kc = [ 19200 8286 1725 64 ];
% Convert back from companion form
k= Kc*inv(W)
%
%% Closed-loop System Poles
% Find poles of closed-loop system. 
% Verify that they are the same as the desired poles.
cls_poles = eig(A-B*k)



states = {'theta' 'theta_dot' 'alpha' 'alpha_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

sysd = c2d(sys_ss,0.001);

%poles = eig(A);

%co = ctrb(sys_ss);
%controllability = rank(co)

Q = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
R = 0.4;
K = lqr(A,B,Q,R)

Ac = (A-B*K);
Bc = B;
Cc = C;
Dc = D;

states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:30;
r =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(AX,'FontName', 'Arial','FontSize',12)
set(get(AX(1),'Ylabel'),'String','rotary angle')
set(get(AX(2),'Ylabel'),'String','Pendulum angle')
title('Step Response')
xlabel('Time')

