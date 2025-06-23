//
// This routine test initial()
//

rfile initial

A=[-0.5572,-0.7814;0.7814,0.0];
B=[1;0];
C=[1.9691,6.4493];
D=[0];
X0=[1,0];
T=0.0:20.0:0.1;


DD=initial(A,B,C,D,X0,T);


