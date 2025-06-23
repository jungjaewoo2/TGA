// test acker.r
rfile acker
rfile c2d

F=[0,1;0,0];
G=[0;1];
H=[1,0];
J=0;
T=1;

A = c2d(F,G,T);
j=sqrt(-1);
pc =[0.78+0.18*j; 0.78-0.18*j];

K = acker(A.phi, A.gamma, pc)

pe =[0.2+0.2*j; 0.2-0.2*j];
L = acker(A.phi', H', pe)'

