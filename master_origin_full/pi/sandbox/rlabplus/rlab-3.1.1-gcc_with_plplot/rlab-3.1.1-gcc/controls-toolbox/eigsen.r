//---------------------------------------------------------------------------
//
// eigsen
//
// Syntax: D = eigsen(A,DADP)
//         where D.de = Gradient of the eigenvalues
//               D.dX = Gradient of the eigenvectors
//
// This routine computes the gradient of the eigenvalues of A and the
// gradient of the eigenvectors of A. It assumes that A has distinct
// eigenvalues. It uses Nelson's algorithm. The routine also normalizes the
// eigenvectors and normalizes the gradient of the eigenvectors.
//
// If the routine is called as,
//
//    D=eigsen(A,DADP)
//
// where A is the square matrix of interest and DADP is the gradient of
// A with respect to some scalar parameter p, it returns the gradient
// of the eigenvalues (D.de) and the gradient of the eigenvalues (D.dX).
//
// The results are returned in a list:
//
//     D.de = Gradient of the eigenvalues
//     D.dX = Gradient of the eigenvectors
//
// Note: The routine makes an effort to check if there are no repeated
//       eigenvalues. However, only a warning is printed and the routine
//       continues. It is up to the user to make sure there are no repeated
//       eigenvalues.
//
// Ref: Nelson, R.B. "Simplified Calculation of Eigenvector Derivatives,"
//      AIAA Journal, Vol. 14, No. 9, 1979, pp. 1201-1205.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//---------------------------------------------------------------------------

require eignm

eigsen = function(A,DADP)
{
   global(eps)

   if (nargs < 2) {
       error("EIGSEN: Wrong number of input arguments.");
   }

// Initial Error Traps
   if ( A.nr <= 1) {
       error("EIGSEN: A must be a matrix.");
       if (A.nr != A.nc) {
           error("EIGSEN: A is not square.");
       }
       if (DADP.nr != DADP.nc) {
           error("EIGSEN: DADP is not square.");
       }
       if ( (DADP.nr != A.nr) && (DADP.nc != A.nc) ) {
           error("EIGSEN: A and DADP are not compatible.");
       }
   }

// Initialize the eigenvector gradient matrix
   dXdp = zeros(A.nr,A.nc);

// Find eigenvalues and eigenvectors (both right and left) of A
   Adum=eig(A);
   D=conj(Adum.val');
   X=Adum.vec;

   AT=A';
   Adum=eig(AT);
   Y=Adum.vec;

// Check for repeated eigenvalues
   for (iloop in 1:(A.nr-1)) {
        arg=Adum.val[iloop]-Adum.val[iloop+1];
        if (abs(real(arg)) < (eps*10) ) {
            if (abs(imag(arg)) < (eps*10) ) {
                error("eigsen: Repeated eigenvalue in A");
            }
        }
   }

// Normalizing the right eigenvectors.
   for (iloop in 1:A.nr) {
       con=X[;iloop]'*X[;iloop];
       if (con > 0.0) {
           X[;iloop]=X[;iloop]/sqrt(con);
       }
   }

// Normalizing the left eigenvectors.
   for (iloop in 1:A.nr) {
       con=Y[;iloop]'*Y[;iloop];
       if (con > 0.0) {
           Y[;iloop]=Y[;iloop]/sqrt(con);
       }
   }

// Could just call eign to do the dirty work
   Adum=eignm(A,eye(A.nr,A.nc));
   D=Adum.lambda;
   X=Adum.phi;
   Y=Adum.psi;

// Initializing the right hand side matrix
   F=zeros(A.nr,1);

// Initialize the eigenvalue gradient vector.
   de=zeros(A.nr,1);

// Looping to find the particular solution in Nelson's
// algorithm

   for (i in 1:A.nr) {

// Compute the gradient of the eigenvalues.
       de[i]=Y[;i].'*DADP*X[;i];

// Finding the RHS of the modified problem of (A-lambda*I)
       if (abs(D[i]) > eps ) {
           F=X[;i]*(de[i])-DADP*X[;i];

// Now Finding the largest value of an element of the eigenvector
           abs(X[;i]).*abs(Y[;i])
           k=maxi(abs(X[;i]).*abs(Y[;i]))

// Setting up the modified set of equations (A-Lambda*I)
           Ap=A-(D[i]*eye(A.nr,A.nr));

           Ap[k;1:A.nr]=zeros(1,A.nr);
           Ap[1:A.nr;k]=zeros(A.nr,1);

           Ap[k;k]=1;

// Setting up the RHS of the problem to find V.
           F[k]=0;

// Solving for the vector V.
           V=solve(Ap,F);
           V[k]=0.0;

// Solving for the constant c.
           C=-real((X[;i]'*V));

// Calculating dX/dp.
           dXdp[;i]=V+C*X[;i];
       else
           D[i]
           error("eigsen: Bug in Looping.");
       }
   }

   return << de=de; dX=dXdp >>;
};

