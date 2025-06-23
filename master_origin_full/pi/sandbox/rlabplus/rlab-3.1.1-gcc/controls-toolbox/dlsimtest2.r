//
// This routine tests the dlsim routine
//

rfile dlsim

A=[0.6597,0.0053;-53.3507,0.1262];
B=[0;0.0053];
C=[1,0;0,1];
D=[0;0];
u=rand(4,1);


J=dlsim(A,B,C,D,u);

