
// This routine tests the continuous dlqr solution methods using examples
// from Laub's paper.

rfile banner
rfile daremsf

banner("Discrete lqr QA");
// =============================================
//                     Test 1:
// =============================================

printf(" \n");
banner("Test1");
printf(" \n");

A=[0.9512,0;0,0.9048];
B=[4.877,4.877;-1.1895,3.569];
R=[(1/3),0;0,3];
Q=[0.005,0;0,0.02];

// Exact Solution
Xsolution=[0.010459082320970,0.003224644477419;
           0.003224644477419,0.050397741135643];

// RLaB dlqr routine solution
X=daremsf(A,B,Q,R);

// Check norm of the error
Xnorm=norm(Xsolution-X.p);
printf("Norm of solution error = ");
sprintf(sdum,"%12.5e",Xnorm);
s="   "+sdum;
printf("%s\n",s);

