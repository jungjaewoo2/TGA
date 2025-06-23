//-----------------------------------------------------------------------------
// ****************
// * ARE2 Testing *
// ****************
// Ref: MATLAB Control Toolbox Manual
// Ref: Laub, A. J., "A Schur Method for Solving Algebraic Riccati
//      Equations," IEEE Trans. AC-24, 1979, pp. 913-921.
//
// Version JBL 950106
//

rfile lqr are are2 banner

// =============================================
//                     Test 1:
// =============================================

printf("\n");
banner("Test1");

A=[0,1;0,0];
B=[0;1];
R=[1];
Q=[1,0;0,2];

P=[2,1;1,2];

Barg=B*inv(R)*B';
Carg=Q;
P2=are2(A,Barg,Carg);

// Check norm of the error
Pnorm=norm(P-P2);
printf ("ARE2 check (this matrix ought to be SMALL)\n");
P-P2
printf("%s","Norm of solution error = ");
sprintf(sdum,"%12.5e",Pnorm);
s="   "+sdum;
printf("%s\n",s);

// =============================================
//                     Test 2:
// =============================================

printf("\n");
banner("Test2");

A=[4,3;-(9/2),-(7/2)];
B=[1;-1];
R=[1];
Q=[9,6;6,4];

// RLaB lqr routine solution
X=lqr(A,B,Q,R);
P=X.s;

Barg=B*inv(R)*B';
Carg=Q;
P2=are2(A,Barg,Carg);

// Check norm of the error
Pnorm=norm(P-P2);
printf ("ARE2 check (this matrix ought to be SMALL)\n");
P-P2
printf("Norm of solution error = ");
sprintf(sdum,"%12.5e",Pnorm);
s="   "+sdum;
printf("%s\n",s);
printf("\n");


// =============================================
//                     Test 1:
// =============================================

printf("\n");
banner("Test1");

A=[0,1;0,0];
B=[0;1];
R=[1];
Q=[1,0;0,2];

P=[2,1;1,2];

Barg=B*inv(R)*B';
Carg=Q;
P2=are(A,Barg,Carg);

// Check norm of the error
Pnorm=norm(P-P2);
printf ("ARE check (this matrix ought to be SMALL)\n");
P-P2
printf("%s","Norm of solution error = ");
sprintf(sdum,"%12.5e",Pnorm);
s="   "+sdum;
printf("%s\n",s);

// =============================================
//                     Test 2:
// =============================================

printf("\n");
banner("Test2");

A=[4,3;-(9/2),-(7/2)];
B=[1;-1];
R=[1];
Q=[9,6;6,4];

// RLaB lqr routine solution
X=lqr(A,B,Q,R);
P=X.s;

Barg=B*inv(R)*B';
Carg=Q;
P2=are(A,Barg,Carg);

// Check norm of the error
Pnorm=norm(P-P2);
printf ("ARE check (this matrix ought to be SMALL)\n");
P-P2
printf("Norm of solution error = ");
sprintf(sdum,"%12.5e",Pnorm);
s="   "+sdum;
printf("%s\n",s);
printf("\n");



