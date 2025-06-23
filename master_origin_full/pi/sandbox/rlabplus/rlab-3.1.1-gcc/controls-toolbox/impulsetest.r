
rfile impulse

A=[-0.5562,-0.7814;0.7814,0];
B=[1;0];
C=[1.9691,6.4493];
D=[0];
//NTIM=512/4;
// old
//impulse(A,B,C,D,NTIM);
T=(0:20:0.1)';
impulse(A,B,C,D,1,T);

