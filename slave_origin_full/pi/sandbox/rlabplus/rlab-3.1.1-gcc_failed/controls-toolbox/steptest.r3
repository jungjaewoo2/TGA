//
// This routine tests the step routine
//

rfile step

// one output variable
A=[-0.5572,-0.7814;0.7814,0.0];
B=[1;0];
C=[1.9691,6.4493];
D=[0];
J=step(A,B,C,D,1);
pause();

// two output variables
C=[1.9691,6.4493;0,1];
D=[0;0];
J=step(A,B,C,D,1);

