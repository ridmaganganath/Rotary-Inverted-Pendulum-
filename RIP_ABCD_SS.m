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

rp_sys = ss(A,B,C,D)
pole=eig(A)


co = ctrb(rp_sys);
controllability = rank(co)


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
    
    eig(A-B*K);
    
    
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
K = lqr(A,B,Q,R)

Ac = (A-B*K);
Bc = B;
Cc = C;
Dc = D;

states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'r'};
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




    
