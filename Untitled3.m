Mp=0.127;
Lp=0.337;
Jp=0.0012;
Jr=0.0020;
Lr=0.216;
g=9.81;
Dr=0.0024;
Dp=0.0024;

Kg=14;
Kt=7.6E-03;
Km=0.725;
Rm=2.6;

% State Space Representation

A = [0     0      1     0
     0     0      0     1_
     0   a_32   a_33  a_34
     0   a_42   a_43  a_44]