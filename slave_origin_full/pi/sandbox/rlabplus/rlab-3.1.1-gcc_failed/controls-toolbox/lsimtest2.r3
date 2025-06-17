//
// This routine tests the lsim routine
//

rfile lsim
rfile rem

A=[0,1;-10000,-4];
B=[0;1];
C=[1,0;0,1];
D=[0;0];
t=(0.0:10.0:0.1)';

period=4.0;
u=(rem(t,period) >= period ./ 2);


J=lsim(A,B,C,D,u,t);

