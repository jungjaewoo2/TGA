//-----------------------------------------------------------------------------
// ********************
// * AUGSTATE Testing *
// ********************
//
// Version JBL 950106
//

rfile augstate

// Ref: MATLAB Control Toolbox Manual
A=[0,1;2,3];
B=[0;1];
C=[1,0];
D=[1];

Dum=augstate(A,B,C,D);

Aex=A;
Bex=B;
Cex=[1,0;1,0;0,1];
Dex=[1;0;0];

// Check results
// Applended A matrix
Enorm=norm(Dum.aa-Aex);
printf ("AUGSTATE: A matrix check (this matrix ought to be ZERO)\n");
Dum.aa-Aex
printf("Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);

// Applended B matrix
Enorm=norm(Dum.ba-Bex,1);
printf ("AUGSTATE: B matrix check (this matrix ought to be ZERO)\n");
Dum.ba-Bex
printf("Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);

// Applended C matrix
Enorm=norm(Dum.ca-Cex);
printf ("AUGSTATE: C matrix check (this matrix ought to be ZERO)\n");
Dum.ca-Cex
printf("Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);

// Applended D matrix
Enorm=norm(Dum.da-Dex,2);
printf ("AUGSTATE: D matrix check (this matrix ought to be ZERO)\n");
Dum.da-Dex
printf("Norm of difference = ");
sprintf(sdum,"%12.5e",Enorm);
s="   "+sdum;
printf("%s\n",s);


