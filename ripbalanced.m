%%%% Assigned the parameters into the system %%%%
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
kt= 0.0145;
km= 0.0145;
Rm=8.71;
eta_g = 0.90;
eta_m = 0.69;
K_ENC = 2 * pi / ( 4 * 1024 );
K_AMP = 1;

%%%% State Space Representation %%%%
Jt = Jr*Jp + Mp*(Lp/2)^2*Jr + Jp*Mp*Lr^2;
A = [0 0 1 0;
     0 0 0 1;
     0 Mp^2*(Lp/2)^2*Lr*g/Jt -Dr*(Jp+Mp*(Lp/2)^2)/Jt -Mp*(Lp/2)*Lr*Dp/Jt;
     0 Mp*g*(Lp/2)*(Jr+Mp*Lr^2)/Jt -Mp*(Lp/2)*Lr*Dr/Jt -Dp*(Jr+Mp*Lr^2)/Jt];
B = [0; 0; (Jp+Mp*(Lp/2)^2)/Jt; Mp*(Lp/2)*Lr/Jt];
C = eye(2,4);
D = zeros(2,1);

%%%% Add actuator dynamics %%%%

B = Kg * kt * B / Rm;
L=-Dr*(Jp+Mp*(Lp/2)^2)/Jt;
N=(Kg^2*kt*km/Rm)*B(3);
A(3,3)=L-N;

R=-Mp*(Lp/2)*Lr*Dr/Jt;
T= (Kg^2*kt*km/Rm)*B(4);
A(4,3) = R-T;

%%%% Load into state-space system %%%%

rp_sys = ss(A,B,C,D)

pole=eig(A)
step(rp_sys)
% TF=tf(rp_sys)
% [p,z] = pzmap(TF)
% [b,a] = ss2tf(A,B,C,D)
%%% Get the Transfer Function through the SS %%%
%[n,d]=ss2tf(A,B,C,D)
%G=tf(n,d)

co = ctrb(rp_sys);
controllability = rank(co)


% IMPORTANT: Make sure you run setup_rotpen.m first. You need the (A,B)
% state-space matrices.
%
%% Find Tranformation Matrix W
% Characteristic equation: s^4 + a_4*s^3 + a_3*s^2 + a_2*s + a_1
a = poly(A);
% 
% Companion matrices (Ac, Bc)
Ac = [  0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        -a(5) -a(4) -a(3) -a(2)];
%
Bc = [0; 0; 0; 1];
% T=ctrb(A,B)
% Tc=ctrb(Ac,Bc);
% w=T*inv(Tc)
% Kc=K/inv(w)


zeta=0.7;
wn=4;
%K = d_balance(A,B,zeta,wn,p3,p4)
    % Location of dominant poles along real-axis
    sigma = zeta*wn;
    % Location of dominant poles along img axis (damped natural freqency)
    wd = wn*sqrt(1-zeta^2);
    % Desired poles (-30 and -40 are given)
    DP = [-sigma+1i*wd, -sigma-1i*wd, -30, -40];
    % Find control gain using Matlab pole-placement command
    K = acker(A,B,DP)
    
    eig(A-B*K)
    
    
    states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

sysd = c2d(sys_ss,0.001);

poles = eig(A);

%co = ctrb(sys_ss);
%controllability = rank(co)

Q = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
R = 1;
klqr = lqr(A,B,Q,R)


Ac = (A-B*klqr);
Bc = B;
%Bc = (B-B*K);
Cc = C;
Dc = D;

states = {'alpha' 'alpha_dot' 'theta' 'theta_dot'};
inputs = {'vm'};
outputs = {'alpha'; 'theta'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs)
step(rp_sys)

t = 0:0.01:30;
r =0.2*ones(size(t));
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(AX,'FontName', 'Arial','FontSize',12)
set(get(AX(1),'Ylabel'),'String','rotary angle')
set(get(AX(2),'Ylabel'),'String','Pendulum angle')
title('Step Response')
xlabel('Time')

%Swing up Control
Ep=Mp*g*Lp;
epsilon = 12.0 * pi / 180;

    
