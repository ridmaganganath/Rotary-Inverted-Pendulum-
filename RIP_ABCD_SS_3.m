
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
controllability = rank(co);


