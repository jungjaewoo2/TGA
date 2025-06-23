
// This routine tests the continuous lqr solution methods using examples
// from Laub's paper.

rfile banner
rfile lqr
rfile are2

banner("Continuous lqr QA");
// =============================================
//                     Test 1:
// =============================================

printf(" \n");
banner("Test1");

A=[0,1;0,0];
B=[0;1];
R=[1];
Q=[1,0;0,2];

P=[2,1;1,2];

P2=are2(A,B,Q,R);

// Check norm of the error
Pnorm=norm(P-P2);
printf("Norm of solution error = ");
sprintf(sdum,"%12.5e",Pnorm);
s="   "+sdum;
printf("%s\n",s);

// =============================================
//                     Test 2:
// =============================================

printf(" \n");
banner("Test2");

A=[4,3;-(9/2),-(7/2)];
B=[1;-1];
R=[1];
Q=[9,6;6,4];

// RLaB lqr routine solution
X=lqr(A,B,Q,R);
P=X.p

P2=are2(A,B,Q,R);

// Check norm of the error
Pnorm=norm(P-P2);
printf("Norm of solution error = ");
sprintf(sdum,"%12.5e",Pnorm);
s="   "+sdum;
printf("%s\n",s);

