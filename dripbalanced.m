Mp=0.127;
Lp=0.33655;
Jp=0.001198;
Jr=9.98E-04;
Lr=0.2159;
g=9.81;
Dr=0.0024;
Dp1= 0.0024;
Dp2= 0.0024;
Mr = 0.257;
Kg=70;
kt= 0.0145;
km= 0.0145;
Rm=8.71;
eta_g = 0.90;
eta_m = 0.69;
K_ENC = 2 * pi / ( 4 * 1024 );
K_AMP = 1;
K_IN2M = 0.0254;
Mh=0.1410;

% Pendulum Mass (with T-fitting)
    Mp1 = 0.127;
    % Pendulum Full Length (with T-fitting, from axis of rotation to tip)
    Lp1 = ( 13 + 1 / 4 ) * K_IN2M;
    % Distance from Pivot to Centre Of Gravity: calculated experimentally
    lp1 = ( 6 + 1 / 8 ) * K_IN2M;
    % Pendulum Moment of Inertia (kg.m^2) - approximation
    Jp1 = Mp * Lp1^2 / 12; 
    % Equivalent Viscous Damping Coefficient (N.m.s/rad)
    Bp1 = 0.0024;

    % Pendulum Mass (with T-fitting)
    Mp2 = 0.097;
    % Pendulum Full Length (with T-fitting, from axis of rotation to tip)
    Lp2 = (7 + 7 / 8) * K_IN2M;
    % Distance from Pivot to Centre Of Gravity: estimated
    lp2 = (6 + 7 / 16) * K_IN2M ;
    % Pendulum Moment of Inertia (kg.m^2)
    Jp2 = Mp1 * Lp2^2 / 12;
    % Equivalent Viscous Damping Coefficient (N.m.s/rad)
    Bp1 = 0.0024;

JT = -2*Mp1*lp1*Lr^2*Mh*Lp1+Mh*Lp1^2*Mp1*Lr^2+Mp1*lp1^2*Mh*Lr^2+Mh*Lp1^2*Jr+Mp1*lp1^2*Jr;

A = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];
    
 
A( 4, 2 ) = Lr*g*(Mp1^2*lp1^2+2*Mp1*lp1*Mh*Lp1+Mh*Lp1^2*Mp2+Mp1*lp1^2*Mp2+Mh^2*Lp1^2)/JT;
A( 4, 3 ) = -Mp1*lp1*Lr*Mp2*g*(-lp1+Lp1)/JT;
A( 4, 4 ) = -Dr*(Mh*Lp1^2+Mp1*lp1^2)/JT;
A( 4, 5 ) = -Lr*Dp1*(Mh*Lp1+Mp1*lp1)/JT;
A( 4, 6 ) = Lr*Dp2*(-Mp1*lp1^2+Mp1*lp1*lp2+Mh*Lp1*lp2+Mp1*lp1*Lp1)/lp2/JT;
A( 5, 1 ) = 0;
A( 5, 2 ) = g*(Lr^2*lp1*Mp1*Mp2+Lp1*Lr^2*Mh*Mp2+Lp1*Lr^2*Mh^2+Lr^2*lp1*Mp1^2+Jr*Lp1*Mh+Jr*lp1*Mp1+Lp1*Lr^2*Mh*Mp1+Lr^2*lp1*Mh*Mp1)/JT;
A( 5, 3 ) = -Mp2*g*(-Lr^2*lp1*Mp1+Lp1*Lr^2*Mp1+Jr*Lp1)/JT;
A( 5, 4 ) = -Dr*Lr*(Mh*Lp1+Mp1*lp1)/JT;
A( 5, 5 ) = -Dp1*(Jr+Mp1*Lr^2+Mh*Lr^2)/JT;
A( 5, 6 ) = Dp2*(Jr*lp2+Jr*Lp1-Lr^2*lp1*Mp1+Lp1*Lr^2*Mp1+Lr^2*lp2*Mh+Lr^2*lp2*Mp1)/lp2/JT;
A( 6, 1 ) = 0;
A( 6, 2 ) = -g/lp2*(Lp1*Lr^2*lp1*Mp1^2+Lp1*Lr^2*lp2*Mh*Mp1+Lr^2*lp1*lp2*Mh*Mp1+Lp1*Lr^2*lp2*Mh^2+Lr^2*lp1*lp2*Mp1^2-Mp1*lp1^2*Lr^2*Mp2+Lr^2*lp1*lp2*Mp1*Mp2+Lp1*Lr^2*lp2*Mh*Mp2+Jr*Lp1*lp2*Mh+Jr*lp1*lp2*Mp1-Mp1*lp1^2*Jr+Lp1*Jr*lp1*Mp1-Mp1*lp1^2*Mh*Lr^2+Lr^2*Lp1*Mp2*Mp1*lp1+Mp1*lp1*Lr^2*Mh*Lp1-Lr^2*Mp1^2*lp1^2)/JT;
A( 6, 3 ) = g/lp2*(Jr*Lp1*lp2*Mp2+Jr*Lp1^2*Mp2+Lp1^2*Lr^2*Mp1*Mp2+Mp1*lp1^2*Lr^2*Mp2+Lp1*Lr^2*lp2*Mp1*Mp2-Lr^2*lp1*lp2*Mp1*Mp2+Mp1*lp1^2*Jr+Mh*Lp1^2*Jr+Mh*Lp1^2*Mp1*Lr^2+Mp1*lp1^2*Mh*Lr^2-2*Lr^2*Lp1*Mp2*Mp1*lp1-2*Mp1*lp1*Lr^2*Mh*Lp1)/JT;
A( 6, 4 ) = Dr*Lr/lp2*(-Mp1*lp1^2+Mp1*lp1*lp2+Mh*Lp1*lp2+Mp1*lp1*Lp1)/JT;
A( 6, 5 ) = Dp1/lp2*(Jr*lp2+Jr*Lp1-Lr^2*lp1*Mp1+Lp1*Lr^2*Mp1+Lr^2*lp2*Mh+Lr^2*lp2*Mp1)/JT;
A( 6, 6 ) = -Dp2*(-2*Lr^2*Lp1*Mp2*Mp1*lp1-2*Mp1*lp1*Lr^2*Mh*Lp1+2*Lp1*Lr^2*lp2*Mp1*Mp2+Jr*Lp1^2*Mp2+Mp2*Jr*lp2^2+Mh*Lp1^2*Jr+Mp1*lp1^2*Jr+2*Jr*Lp1*lp2*Mp2-2*Lr^2*lp1*lp2*Mp1*Mp2+Mh*Lp1^2*Mp1*Lr^2+Mp1*lp1^2*Mh*Lr^2+Mp1*lp1^2*Lr^2*Mp2+Lp1^2*Lr^2*Mp1*Mp2+Mp2*Lr^2*lp2^2*Mh+Mp2*Lr^2*lp2^2*Mp1)/lp2^2/JT/Mp2;

B = [0;
    0;
    0;
    (Mh*Lp1^2+Mp1*lp1^2)/JT;
    Lr*(Mh*Lp1+Mp1*lp1)/JT;
    -Lr/lp2*(-Mp1*lp1^2+Mp1*lp1*lp2+Mh*Lp1*lp2+Mp1*lp1*Lp1)/JT];

C = eye(6,6);
D = zeros(6,1);

% Add actuator dynamics
B = Kg * kt * B / Rm;
A(4,4) = A(4,4) - Kg*km/Rm*B(4);
A(5,4) = A(5,4) - Kg*km/Rm*B(5);
A(6,4) = A(6,4) - Kg*km/Rm*B(6);


rp_sys= ss(A,B,C,D)
pole=eig(A)


pzmap(rp_sys)
step(rp_sys)
co = ctrb(rp_sys);
controllability = rank(co)


Ai=A;  
Bi=B;
Ai(7,1) = 1; Ai(7,7) = 0; 
Bi(7,1) = 0;


% %controller design
% 
%  

%Q = [1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0; 0 0 0 5 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 0.5];
Q = diag([ 1 1 1 5 1 1 0.5]);
R = 30;        
k = lqr(Ai, Bi, Q, R)
k


Ac=(Ai-Bi*k)
Bc=Bi
Cc=C
Dc=D
eig(Ac)

%pzmap(Ai-Bi*k)
%step(rp_sys)
%step(Ac)

%sys_cl = ss(Aci,Bci,Cc,Dc)



%LQRcontroller_sys=ss(Ai-Bi*K,Bi,C,D)

X0 = [-5, 1, -0.5, 0, 0, 0] * pi / 180

states = {'theta' 'theta_dot' 'alpha' 'alpha_dot' 'phi'  'phi_dot'  };
inputs = {'vm'};
outputs = {'x1'; 'theta'  };

%sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

% 
% t = 0:0.01:10;
% r =0.2*ones(size(t));
% [y,t,x]=lsim(rp_sys,r,t);
% [AX,H1,H2,H3] = plotyy(t,y(:,1),t,y(:,2),t,y(:,3),'plot');
% set(AX,'FontName', 'Arial','FontSize',12)
% set(get(AX(1),'Ylabel'),'String','Linear displacement')
% set(get(AX(2),'Ylabel'),'String','Tilt angular displacement')
% set(get(AX(3),'Ylabel'),'String','Tilt angular displacement')
% title('Step Response')
% xlabel('Time')

