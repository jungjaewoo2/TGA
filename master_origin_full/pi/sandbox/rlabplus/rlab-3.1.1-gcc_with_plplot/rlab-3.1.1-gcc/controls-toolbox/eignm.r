//---------------------------------------------------------------------
//
// eignm
//
// Syntax: D = eignm(A,B)
//         where D.lambda = eigenvalues
//               D.phi    = right eigenvectors
//               D.psi    = left eigenvectors
//
// This routine solves the generalized eigenvalues problem given by,
//
//    A*x = lambda*B*x
//
// It then orders and normalizes the eigenvectors according to the
// following normalizations:
//
//    phi(i)'*B*phi(i) = 1
//    psi(i)'*A*phi(j) = indrealnecker delta(i,j)
//
// where phi is the right eigenvector matrix, psi is the left eigenvector
// matrix and where phi(i) is the ith right eigenvector and psi(i) is the
// ith left eigenvector, indrealnecker delta is the Kronecker Delta function.
//
// The routine sorts them in the following order:
//
//   real eigenvalues (ascending order)
//   comples eigenvalues (ascending order according to real part)
//
// where "real" is defined as imag(lambda) < 1.0e-12
//
// The results are returned in a list:
//
//     D.lambda = eigenvalues
//     D.phi    = right eigenvectors
//     D.psi    = left eigenvectors
//
// Ref:  Junkins and Kim, Dyn. & Ctrl of Structures, Introduction to Dynamics
//       and Control of Flexible Structures, AIAA, 1993.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//---------------------------------------------------------------------

eignm = function(A,B)
{

   if (nargs < 2) {
       error("EIGNM: Wrong number of input arguments.");
   }

// Dimension Check
   if (A.nr != B.nr) {
       error("eign: A and B don't have the same number of rows.");
   }

// Initialize needed matrices
// Note: This RLaB from having to resize as elements are added.
   n=max(size(A));
   lambda=zeros(n,1);
   phi=zeros(n,n);
   psi=zeros(n,n);

// Solve Right Eigenvalues Problem.
   Adum=eig(A,B);
   LAMR=Adum.vec;
   PHI=Adum.val.';

// Solve Left Eigenvalue Problem.
   Adum=eig(A',B');
   LAML=Adum.vec;
   PSI=Adum.val.';
   clear(Adum);

// Initialize sorted storage
// Note: This keeps RLaB from having to resize as elements are added.
   indreal=zeros(n,1);
   indcompx=zeros(n,1);
   lamreal=zeros(n,1);
   lamcompx=zeros(n,1);

// Sort Right Eigenvectors and Eigenvalues
   numr=0;
   numc=0;
   for (i in 1:n) {
        if (abs(imag(PHI[i])) <= 1.e-12) {
            numr=numr+1;
            indreal[numr]=i;
            lamreal[numr]=PHI[i];
        else if (imag(PHI[i]) > 1.0e-12) {
            numc=numc+1;
            indcompx[numc]=i;
            lamcompx[numc]=PHI[i];
        }}
   }

// Sort Real and Complex Right Eigenvalues
   lamreal=lamreal[1:numr];
   lamreal=real(lamreal);
   lamcompx=lamcompx[1:numc];
   indreals=sort(lamreal).idx;
   indcompxs=sort(lamcompx).idx;

// Extract Sorted Right Eigenvectors and put them in output matrix
   ic=1;
   for (i in 1:(numr+numc)) {
        if (i <= numr) {
            phi[;i]=real(LAMR[;indreal[indreals[i]]]);
            lambda[i]=real(PHI[indreal[indreals[i]]]);
            ic=ic+1;
        else
            ip=i-numr;
            phi[;ic]=LAMR[;indcompx[indcompxs[ip]]];
            phi[;ic+1]=conj(phi[;ic]);
            lambda[ic]=PHI[indcompx[indcompxs[ip]]];
            lambda[ic+1]=conj(lambda[ic]);
            ic=ic+2;
        }
   }

// Initialize sorted storage
// Note: This RLaB from having to resize as elements are added.
   indreal=zeros(n,1);
   indcompx=zeros(n,1);
   lamreal=zeros(n,1);
   lamcompx=zeros(n,1);

// Sort Left Eigenvectors and "left" Eigenvalues.
   numr=0;
   numc=0;
   for (i in 1:n) {
        if (abs(imag(PSI[i])) <= 1.0e-12) {
            numr=numr+1;
            indreal[numr]=i;
            lamreal[numr]=PSI[i];
        else if (imag(PSI[i]) > 1.0e-12) {
            numc=numc+1;
            indcompx[numc]=i;
            lamcompx[numc]=PSI[i];
        }}
   }

// Sort Real and Complex Left Eigenvalues
   lamreal=lamreal[1:numr];
   lamreal=real(lamreal);
   lamcompx=lamcompx[1:numc];
   indreals=sort(lamreal).idx;
   indcompxs=sort(lamcompx).idx;

// Extract Sorted Left Eigenvectors and put them in output matrix
   ic=1;
   for (i in 1:numr+numc) {
       if (i <= numr) {
           psi[;i]=real(LAML[;indreal[indreals[i]]]);
           ic=ic+1;
       else
           ip=i-numr;
           psi[;ic]=LAML[;indcompx[indcompxs[ip]]];
           psi[;ic+1]=conj(psi[;ic]);
           ic=ic+2;
      }
   }

// Normalize the Sorted eigenvectors
   for (i in 1:n) {
        phii=phi[;i];
        psii=psi[;i];
        term1=conj(phii')*B*phii;
        phi[;i]=phi[;i]/sqrt(term1);
        term1=conj(psii')*B*phi[;i];
        psi[;i]=phi[;i]/term1;
   }

   return << lambda=lambda; phi=phi; psi=psi >>;
};

