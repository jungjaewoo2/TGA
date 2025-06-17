
// This routine tests the continuous lqr solution methods using examples
// from Laub's paper.

rfile banner
rfile lqr
rfile lqrmsf

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

// Exact Solution
Xsolution=[2,1;1,2];

// RLaB lqr routine solution
X=lqr(A,B,Q,R);

// RLaB lqrmsf routine solution
X2=lqrmsf(A,B,Q,R);

// Check norm of the error
Xnorm=norm(Xsolution-X.s);
X2norm=norm(Xsolution-X2.P);
printf("Norm of solution error = ");
sprintf(sdum,"%12.5e",Xnorm);
s="   "+sdum;
printf("%s\n",s);
X2norm

// =============================================
//                     Test 2:
// =============================================

printf(" \n");
banner("Test2");

A=[4,3;-(9/2),-(7/2)];
B=[1;-1];
R=[1];
Q=[9,6;6,4];

// Exact Solution
c=1.0+sqrt(2);
Xsolution=[9*c,6*c;6*c,4*c];

// RLaB lqr routine solution
X=lqr(A,B,Q,R);

// RLaB lqrmsf routine solution
X2=lqrmsf(A,B,Q,R);

// Check norm of the error
Xnorm=norm(Xsolution-X.s);
printf("Norm of solution error = ");
sprintf(sdum,"%12.5e",Xnorm);
s="   "+sdum;
printf("%s\n",s);
X2norm=norm(Xsolution-X2.P)

