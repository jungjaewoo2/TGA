//-----------------------------------------------------------------------------
// ******************
// * APPEND Testing *
// ******************
//
// Version JBL 950106
//

rfile append

// Ref: MATLAB Control Toolbox Manual
A1=[0,1;2,3];
B1=[0;1];
C1=[1,0];
D1=[1];

A2=[0,1,0;-1,-3,-0.5;0,0,4];
B2=[0;1;1];
C2=[1,1,1];
D2=[0];

Dum=append(A1,B1,C1,D1,A2,B2,C2,D2);

Aex=[0,1,0,0,0;2,3,0,0,0;0,0,0,1,0;0,0,-1,-3,-0.5;
     0,0,0,0,4];
Bex=[0,0;1,0;0,0;0,1;0,1];
Cex=[1,0,0,0,0;0,0,1,1,1];
Dex=[1,0;0,0];

// Check results
// Appended A matrix
Enorm=norm(Dum.aa-Aex);
printf ("APPLEND: A matrix check (this matrix ought to be ZERO)\n");
Dum.aa-Aex
printf("Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);

// Appended B matrix
Enorm=norm(Dum.ba-Bex);
printf ("APPLEND: B matrix check (this matrix ought to be ZERO)\n");
Dum.ba-Bex
printf("%s","Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);

// Appended C matrix
Enorm=norm(Dum.ca-Cex);
printf ("APPLEND: C matrix check (this matrix ought to be ZERO)\n");
Dum.ca-Cex
printf("Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);

// Appended D matrix
Enorm=norm(Dum.da-Dex);
printf ("APPLEND: D matrix check (this matrix ought to be ZERO)\n");
Dum.da-Dex
printf("Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);


