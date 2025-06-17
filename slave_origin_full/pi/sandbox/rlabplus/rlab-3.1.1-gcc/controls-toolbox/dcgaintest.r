
// This rfile test the covar routine

rfile dcovar

A=[3,0,0;0,3.2,0;0,0,3.5];
B=[2;-1;-1];
C=[1,0,0];
D=[0];
W=0.5;

Adum=dcovar(A,B,C,D,W);
Adum.X
Adum.Y


