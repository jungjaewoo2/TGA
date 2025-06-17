//-------------------------------------------------------------------------------
//
// eigsen
//
// Syntax: D = eigsen(A,B,DADP,DBDP)
//         where D.de   = Gradient of the Eigenvalues
//               D.dPhi = Gradient of the Right Eigenvectors
//               D.dPsi = Gradient of the Left Eigenvectors
//
// This routine computes the gradient of the eigenvalues, right eigenvectors,
// and left eigenvectors of a matrix with respect to a scalar parameter.It
// considers a general eigenvalue problem of the form:
//
//    Ax = (lambda) B x
//
// with the normalizations:
//
//    phi(i)'*B*phi(i) = 1
//    psi(i)'*A*phi(j) = krnecker delta(i,j)
//
// where phi is right eigenvector matrix, psi is the let eigenvector matrix,
// and where phi(i) is the ith right eigenvector and psi(i) is the ith left
// eigenvector, krnecker delta is the Kronecker Delta function.
//
// If the routine is called as,
//
//    D=eigsen2(A,B,DADP,DBDP)
//
// where DADP and DBDP are the gradients of A and B respectively, with respect
// to a scalar p, it returns the gradient of the eigenvalues (D.de), the
// gradient of the right eigenvectors (D.dPhi), and the gradient of the
// left eigenvectors (D.dPsi).
//
// The results are returned in a list:
//
//     D.de = Gradient of the Eigenvalues
//   D.dPhi = Gradient of the Right Eigenvectors
//   D.dPsi = Gradient of the Left Eigenvectors
//
// Note: The routine makes an effort to check if there are no repeated
//       eigenvalues. However, only a warning is printed and the routine
//       continues. It is up to the user to make sure there are no repeated
//       eigenvalues.
//
// Refs: (1) Lim, K.B., Junkins, J.L., and Wang, B.P.
//           "Re-examination of Eigenvector Derivatives"
//           Journal of Guidance, Control, and Dynamics, Vol. 10, No. 6, 1987,
//           pp. 581-587.
//
//       (2) Junkins and Kim, Dyn. & Ctrl of Structures, Introduction to Dynamics
//           and Control of Flexible Structures, AIAA, 1993.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//-------------------------------------------------------------------------------

require eignm

eigsen2 = function(A,B,DADP,DBDP)
{

   if (nargs < 4) {
       error("EIGSEN2: Wrong number of input arguments.");
   }

// Compute the spectral decomposition using the correct normalization
// by using the routine "eign"
   Adum=eignm(A,B);
   Lambda=Adum.lambda;

// Initialization
   n=max(size(Lambda));

// Initialize Results
   a=zeros(n,1);
   de=zeros(n,1);
   dPhi=zeros(n,n);
   dPsi=zeros(n,n);

// Begin Computations
   for (i in 1:n) {
        xi=Adum.phi[;i];
        yi=Adum.psi[;i];
        de[i]=yi.'*(DADP-Lambda[i]*DBDP)*xi;
        ssum=0;
        for (k in 1:n) {
             if (k != i) {
                 xk=Adum.phi[;k];
                 yk=Adum.psi[;k];
                 a[k]=yk.'*(DADP-Lambda[i]*DBDP)*xi/(Lambda[i]-Lambda[k]);
                 ssum=ssum+a[k]*xk.'*(B+B')*xi;
                 b[k]=yi.'*(DADP-Lambda[i]*DBDP)*xk/(Lambda[i]-Lambda[k]);
             }
        }
// Add in term when k=i
        a[i]=-(xi.'*DBDP*xi+ssum)/2;
        b[i]=-yi.'*DBDP*xi-a[i];
        for (k in 1:n) {
             dPhi[;i]=dPhi[;i]+a[k]*Adum.phi[;k];
             dPsi[;i]=dPsi[;i]+b[k]*Adum.psi[;k];
        }
   }

   return << de=de; dPhi=dPhi; dPsi=dPsi >>;
};

