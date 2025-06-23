
// This routine test the eigen sensitivity routines.

rfile eigsen
rfile eigsen2
rfile banner
rfile eignm

banner("Eigen Sensitivity Test.");
// Generate a random 4x4 matrix for A
A=5.*rand(4,4);
A=real(A);
//A=(A+A')/2;

// Generate a random 4x4 matrix for d(A)/dp matrix
DADP=2.*rand(4,4);

// Set B=I for eigsen (Junkins method)
B=eye(4,4);

// Set d(B)/dp = 0
DBDP = zeros(4,4);

// Compute eigenvalues, right and left eigenvectors of A
E=eig(A);
X=E.vec;
Lambda=E.val;
ET=eig(A');
Y=ET.vec;

// Clear unnecessary variables
clear(E,ET);

Adum=eignm(A,B);
X=Adum.Phi;
Y=Adum.Psi;
Lambda=Adum.Lambda;
// Compute Eigen sensitivity using Junkins' Method.
Es=eigsen(A,B,DADP,DBDP);
//clear(Adum);

// Compute Eigen sensitivity using Nelson's Method.
Es2=eigsen2(A,DADP);

// Compute norms of the errors
REnorm=norm(Es.REVTD-Es2.dX);

Lnorm=norm(Es.EVAD-Es2.de);

// Clean-up
clear(DADP,B,DBDP,X,Y,Lambda);

if (!(REnorm < 1.0e-12) ) {
    printf("FAILED Right Eigen Sensitivity test.\n");
    printf("The norm of the Right Eigenvector Sensitivity Error is");
    sprintf(sdum,"%12.5e",REnorm);
    s="   "+sdum;
    printf("%s\n",s);
} else {
    printf("Passed Right Eigenvector Sensitivity Check.\n");
    printf("The norm of the Right Eigenvector Sensitivity Error is");
    sprintf(sdum,"%12.5e",REnorm);
    s="   "+sdum;
    printf("%s\n",s);
}

printf(" \n");
if (!(Lnorm < 1.0e-12) ) {
    printf("Failed Eigenvalue sensitivity test.\n");
    printf("The norm of the Eigenvalue Sensitivity Error is");
    sprintf(sdum,"%12.5e",Lnorm);
    s="   "+sdum;
    printf("%s\n",s);
} else {
    printf("Passed Eigenvalue Sensitivity Test.\n");
    printf("The norm of the Eigenvalue Sensitivity Error is");
    sprintf(sdum,"%12.5e",Lnorm);
    s="   "+sdum;
    printf("%s\n",s);
}

clear(REnorm,Lnorm,s,sdum);

