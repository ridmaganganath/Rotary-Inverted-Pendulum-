
Mp=0.230;
Lp=0.6413;
Jp=0.0023;
Jr=0.0023;

Lr=0.2159;
g=9.81;
Dr=0.0024;
Dp=0.0024;

Kg=70;
kt= 0.00767;
km= 0.00767;
Rm=2.6;

A_32 = Mp^2*(Lp/2)^2*Lr*g/Jt;

A_33 = (-Dr*(Jp+Mp*(Lp/2)^2)/Jt)- Kg^2*kt*km/Rm*B(3);

A_34 = -Mp*(Lp/2)*Lr*Dp/Jt;

A_42 = Mp*g*(Lp/2)*(Jr+Mp*Lr^2)/Jt;

A_43 = (-Mp*(Lp/2)*Lr*Dr/Jt)- Kg^2*kt*km/Rm*B(4);

A_44 = Dp*(Jr+Mp*Lr^2)/Jt];

B_3 = (Jp+Mp*(Lp/2)^2)/Jt;

B_4 = Mp*(Lp/2)*Lr/Jt];

%----------------------State Matrices-----------------------------------%

A = [0   0    1   0
     0  0 0 1
     0   A_32  A_33     A_34
     0  A_42 A_43 A_44];

B = [0
     0
     B_3
     B_4];

C = [1 0 0 0
     0 1 0 0];

D = [0
     0];
B = Kg * kt * B / Rm;

states = {'theta' 'theta_dot' 'alpha' 'alpha_dot'};
inputs = {'u'};
outputs = {'x'; 'theta'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs)