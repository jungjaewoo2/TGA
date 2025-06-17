/* Copyright (C) 2006 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 harminv was written by Steven G. Johnson (stevenj@alum.mit.edu).
 */


struct harminv_data_struct
{
  const Complex *c;
  int n, K, J, nfreqs;
  double fmin, fmax;
  Complex *z;
  Complex *U0, *U1;
  Complex *G0, *G0_M, *D0; /* cached G, G_M, and D arrays for U0 */
  Complex *B, *u;  /* eigen-solutions of U1*B = u*U0*B */
  Complex *amps; /* mode amplitudes */
  double *errs; /* relative "error" estimates */
};

typedef struct
{
  int verbose;
  double fmin, fmax;
  int only_f_inrange;
  double err_thresh, rel_err_thresh, amp_thresh, rel_amp_thresh, Q_thresh;
  double min_err, max_amp;
  int num_ok;
} mode_ok_data;


typedef struct harminv_data_struct *harminv_data;
typedef int (*harminv_mode_ok_func)(harminv_data d, int k, void *);

#define CHECK(condition, message) do { \
     if (!(condition))  { \
          fprintf(stderr, "harminv: failure on line %d of " __FILE__ ": " \
      message "\n", __LINE__);  rerror("harminv\n"); \
} \
} while (0)

#define CHK_MALLOC(p, t, n) do {                              \
     size_t CHK_MALLOC_n_tmp = (n);                           \
     (p) = (t *) GC_malloc(sizeof(t) * CHK_MALLOC_n_tmp);        \
     CHECK((p) || CHK_MALLOC_n_tmp == 0, "out of memory!");   \
} while (0)

static char cn='N';
static char cv='V';
static char ct='T';


/**************************************************************************/

/* The harminv routines are designed to perform "harmonic inversion."
   That is, given a signal (a set of samples as a function of time),
   they decompose the signal into a finite number of (possibly
   exponentially-decaying) sinusoids.

   Essentially, because we assume that the signal has this form, we
   can determine the frequencies, decay constants, and amplitudes of
   the sinusoids much more accurately than we could via taking the FFT
   and looking at the peaks, for the same number of samples.

   We use a low-storage "filter diagonalization method" (FDM) for finding
   the sinusoids near a given frequency interval, described in:

   V. A. Mandelshtam and H. S. Taylor, "Harmonic inversion of time
   signals," J. Chem. Phys., vol. 107, no. 17, p. 6756-6769 (Nov. 1
   1997).  See also erratum, ibid, vol. 109, no. 10, p. 4128 (Sep. 8
   1998).

   with a refinement (for generate_U below) described in:

   Rongqing Chen and Hua Guo, "Efficient calculation of matrix
   elements in low storate filter diagonalization," J. Chem. Phys.,
   vol. 111, no. 2, p. 464-471(Jul. 8 1999).

   The seminal work (though less practical than the M&T algorithm) on
   this class of methods was done by:

   Michael R. Wall and Daniel Neuhauser, "Extraction, through
   filter-diagonalization, of general quantum eigenvalues or classical
   normal mode frequencies from a small number of residues or a
   short-time segment of a signal. I. Theory and application to a
   quantum-dynamics model," J. Chem. Phys., 102, no. 20, p. 8011-8022
   (May 22 1995).

   A more recent reference is:

   V. A. Mandelshtam, "On harmonic inversion of cross-correlation
   functions by the filter diagonalization method," J. Theoretical and
   Computational Chemistry 2 (4), 497-505 (2003).

*/

#define TWOPI 6.2831853071795864769252867665590057683943388

/* Computing powers by cumulative multiplication, below, is faster
   than calling complex_Pow_int repeatedly, but accumulates O(n) floating-point
   error.  As a compromise, we call complex_Pow_int every NPOW iterations, which
   accumulates only O(NPOW) error. */
#define NPOW 8

/* FIXME: instead of this, we should really do the first-order expansion
   of the U matrix in |z - z2|, below. */
#define SMALL (1e-12)
#define C_CLOSE(c1,c2) (cabs((c1) - (c2)) < SMALL)

/**************************************************************************/

int harminv_get_num_freqs(const harminv_data d)
{
  return d->nfreqs;
}


/* Initialize the JxJ2 matrix U = U_p(z,z2), as described in M&T.
   Also, if U1 != NULL, then set U1 = U_{p+1}(z,z2).  If z == z2, it
   must be the case that no two elements of z are the same and that
   J2 == J1; in this case the matrix U will be symmetric.  Note that c
   must be an array whose size n, is at least 2*K+p elements.  */
static void generate_U(Complex *U, Complex *U1,
           int p,
           const Complex *c, int n,
           int K,
           int J, int J2, const Complex *z, const Complex *z2,
           Complex **G0, Complex **G0_M, Complex **D0)
{
     int M = K - 1;
     int i, j, m;

/*
     temp. arrays for 1/z, z^(-m), z^(-M), the G function of C&G,
     and the diagonal elements D[i] = U(z[i],z[i]):
*/
     Complex *z_inv, *z_m, *z_M;
     Complex *G, *G_M, *D;
     Complex *z2_inv, *z2_m, *z2_M;
     Complex *G2, *G2_M;

     CHECK(U && c && z && z2, "invalid arguments to generate_U");
     CHECK(n >= 2*K + p, "too few coefficients in generate_U");
     CHECK(z != z2 || J == J2, "invalid sizes passed to generate_U");
     CHECK((!G0 && !G0_M && !D0) || (G0 && G0_M && D0),
     "G0/G0_M/D0 must be all-non-NULL/all-NULL");
     CHECK(!G0 || (!*G0 && !*G0_M && !*D0) || (*G0 && *G0_M && *D0),
     "*G0/*G0_M/*D0 must be all-non-NULL/all-NULL");

     /* Now, compute U according to eqs. 25-27 of Chen & Guo, but
        using the notation of eq. 25 of M&T.  This operation has
        complexity O(N*J + J*J).  At the same time, we can compute the
        matrix U1 as well by eqs. 29-30 of C&G, saving an extra pass
        over the input array. */

     /* first, set up some temporary arrays for caching things like
  z^m and 1/z, so we don't need to recompute them all the time. */
     CHK_MALLOC(z_inv, Complex, J);
     CHK_MALLOC(z_m, Complex, J);
     CHK_MALLOC(z_M, Complex, J);
     if (G0 && *G0)
     {
       G = *G0;
       G_M = *G0_M;
       D = *D0;
     }
     else
     {
       CHK_MALLOC(G, Complex, J);
       CHK_MALLOC(G_M, Complex, J);
       CHK_MALLOC(D, Complex, J);
       for (i = 0; i < J; ++i)
       { D[i] = G[i] = G_M[i] = 0; }
     }
     for (i = 0; i < J; ++i)
     {
       z_inv[i] = 1.0 / z[i];
       z_m[i] = 1;
       z_M[i] = complex_Pow_int(z[i], -M);
     }
     if (z2 != z)
     {
       CHK_MALLOC(z2_inv, Complex, J2);
       CHK_MALLOC(z2_m, Complex, J2);
       CHK_MALLOC(z2_M, Complex, J2);
       CHK_MALLOC(G2, Complex, J2);
       CHK_MALLOC(G2_M, Complex, J2);
       for (i = 0; i < J2; ++i)
       {
         z2_inv[i] = 1.0 / z2[i];
         z2_m[i] = 1;
         z2_M[i] = complex_Pow_int(z2[i], -M);
         G2[i] = G2_M[i] = 0;
       }
     }
     else
     {
       z2_inv = z2_m = z2_M = NULL;
       G2 = G2_M = NULL;
     }

     /* First, loop over the signal array (c), building up the
     spectral functions G and G_M (corresponding to G_p and
     G_{p+M+1} in C&G), as well as the diagonal matrix entries: */
     for (m = 0; m <= M; ++m)
     {
       Complex c1 = c[m + p], c2 = c[m + p + M + 1];
       double d = m + 1; /* M - fabs(M - m) + 1 */
       double d2 = M - m; /* M - fabs(M - (m + M + 1)) + 1 */

       if (!G0 || !*G0)
       {
         for (i = 0; i < J; ++i)
         {
           Complex x1 = z_m[i] * c1;
           Complex x2 = z_m[i] * c2;
           G[i] += x1;
           G_M[i] += x2;
           D[i] += x1 * d + x2 * d2 * z_M[i] * z_inv[i];
           if (m % NPOW == NPOW - 1)
             z_m[i] = complex_Pow_int(z_inv[i], m + 1);
           else
             z_m[i] *= z_inv[i];
         }
       }
       if (z2 != z)
         for (i = 0; i < J2; ++i)
       {
         G2[i] += z2_m[i] * c1;
         G2_M[i] += z2_m[i] * c2;
         if (m % NPOW == NPOW - 1)
           z2_m[i] = complex_Pow_int(z2_inv[i], m + 1);
         else
           z2_m[i] *= z2_inv[i];
       }
     }

     /* Compute U (or just the upper part if U is symmetric), via the
     formula from C&G; compute U1 at the same time as in C&G. */
     if (z2 != z)
     {
       for (i = 0; i < J; ++i)
         for (j = 0; j < J2; ++j)
       {
         if (C_CLOSE(z[i], z2[j]))
           U[i*J2 + j] = D[i];
         else
           U[i*J2 + j] = (z[i] * G2[j] - z2[j] * G[i] +
               z2_M[j] * G_M[i] - z_M[i] * G2_M[j])
            / (z[i] - z2[j]);
       }

       if (U1)
         for (i = 0; i < J; ++i)
           for (j = 0; j < J2; ++j)
       {
         if (C_CLOSE(z[i], z2[j]))
           U1[i*J2 + j] = z[i] * (D[i] - G[i]) +
               z_M[i] * G_M[i];
         else
           U1[i*J2 + j] = (z[i] * z2[j] * (G2[j] - G[i])
               + z2_M[j] * z[i] * G_M[i]
               - z_M[i] * z2[j] * G2_M[j])
           / (z[i] - z2[j]);
       }
     }
     else
     {  /* z == z2 */
       for (i = 0; i < J; ++i)
       {
         U[i*J + i] = D[i];
         for (j = i + 1; j < J; ++j)
         {
           U[i*J + j] = (z[i] * G[j] - z[j] * G[i] +
               z_M[j] * G_M[i] - z_M[i] * G_M[j])
       / (z[i] - z[j]);
         }
       }
       if (U1)
         for (i = 0; i < J; ++i)
       {
         U1[i*J + i] = z[i] * (D[i] - G[i]) + z_M[i] * G_M[i];
         for (j = i + 1; j < J; ++j) {
           U1[i*J + j] = (z[i] * z[j] * (G[j] - G[i])
               + z_M[j] * z[i] * G_M[i]
               - z_M[i] * z[j] * G_M[j])
            / (z[i] - z[j]);
         }
       }
     }

     /* finally, copy the upper to the lower triangle if U is symmetric: */
     if (z == z2)
     {
       for (i = 0; i < J; ++i)
         for (j = i + 1; j < J; ++j)
           U[j*J + i] = U[i*J + j];
       if (U1)
         for (i = 0; i < J; ++i)
           for (j = i + 1; j < J; ++j)
             U1[j*J + i] = U1[i*J + j];
     }

     GC_free(G2_M);
     GC_free(G2);
     GC_free(z2_M);
     GC_free(z2_m);
     GC_free(z2_inv);

     if (G0 && !*G0)
     {
       *G0 = G;
       *G0_M = G_M;
       *D0 = D;
     }
     else if (!G0)
     {
       GC_free(D);
       GC_free(G_M);
       GC_free(G);
     }

     GC_free(z_M);
     GC_free(z_m);
     GC_free(z_inv);
}

/**************************************************************************/

static void init_z(harminv_data d, int J, Complex *z)
{
  d->J = J;
  d->z = z;
  CHK_MALLOC(d->U0, Complex, J*J);
  CHK_MALLOC(d->U1, Complex, J*J);
  generate_U(d->U0, d->U1, 0, d->c, d->n, d->K, d->J, d->J, d->z, d->z,
             &d->G0, &d->G0_M, &d->D0);
}

/**************************************************************************/

harminv_data harminv_data_create(int n,
         const Complex *signal,
         double fmin, double fmax, int nf)
{
  int i;
  harminv_data d;

  CHECK(nf > 1, "# frequencies must > 1");
  CHECK(n > 0, "invalid number of data points");
  CHECK(signal, "invalid NULL signal array");
  CHECK(fmin < fmax, "should have fmin < fmax");

  CHK_MALLOC(d, struct harminv_data_struct, 1);
  d->c = signal;
  d->n = n;
  d->K = n/2 - 1;
  d->fmin = fmin;
  d->fmax = fmax;
  d->nfreqs = -1;  /* we haven't computed eigen-solutions yet */
  d->B = d->u = d->amps = d->U0 = d->U1 = (Complex *) NULL;
  d->G0 = d->G0_M = d->D0 = (Complex *) NULL;
  d->errs = (double *) NULL;

  CHK_MALLOC(d->z, Complex, nf);

  for (i = 0; i < nf; ++i)

    d->z[i] = cexp(-I * TWOPI * (fmin + i * ((fmax - fmin) / (nf - 1))));

  init_z(d, nf, d->z);

  return d;
}

/**************************************************************************/

void harminv_data_destroy(harminv_data d)
{
  if (d)
  {
    d->c = NULL; // signal is pointer from somewhere else
    GC_free(d->B);    GC_free(d->u);    GC_free(d->amps);
    GC_free(d->U0);   GC_free(d->U1);
    GC_free(d->G0);   GC_free(d->G0_M); GC_free(d->D0);
    GC_free(d->errs); GC_free(d->z);

    GC_free(d);
    d = NULL;
  }
}

/**************************************************************************/

/* Compute the symmetric dot product of x and y, both vectors of
   length n.  If they are column-vectors, this is: transpose(x) * y.
   (We could use the BLAS ZDOTU function for this, but calling Fortran
    functions, as opposed to subroutines, from C is problematic.) */
static Complex symmetric_dot(int n, Complex *x, Complex *y)
{
  Complex dot = 0;
  int i;
  for (i = 0; i < n; ++i)
    dot += x[i] * y[i];
  return dot;
}

/**************************************************************************/

/**************************************************************************/

/* Solve for the eigenvalues (v) and eigenvectors (rows of V) of the
   complex-symmetric n x n matrix A.  The eigenvectors are normalized
   to 1 according to the symmetric dot product (i.e. no complex
   conjugation). */
static void solve_eigenvects(int n, const Complex *A0, Complex *V, Complex *v)
{
  int lwork, info;
  Complex *work;
  double *rwork;
  Complex *A;

/*
  according to the ZGEEV documentation, the matrix A is overwritten,
  and we don't want to overwrite our input A0
*/

  CHK_MALLOC(A, Complex, n*n);
  {
    int n2 = n*n, one = 1;
    ZCOPY(&n2, A0, &one, A, &one);
  }

/*
  Unfortunately, LAPACK doesn't have a special solver for the
  complex-symmetric eigenproblem.  For now, just use the general
  non-symmetric solver, and realize that the left eigenvectors
  are the complex-conjugates of the right eigenvectors.
*/

  lwork = 4*n; /* minimum is 2*n; we'll be generous. */
  CHK_MALLOC(rwork, double, 2*n);
  CHK_MALLOC(work,  Complex,  lwork);

  ZGEEV(&cn, &cv, &n, A, &n, v, V, &n, V, &n, work, &lwork, rwork, &info);

  GC_free(work);
  GC_free(rwork);
  GC_free(A);

  CHECK(info >= 0, "invalid argument to ZGEEV");
  CHECK(info <= 0, "failed convergence in ZGEEV");

/*
  Finally, we need to fix the normalization of the eigenvectors,
  since LAPACK normalizes them under the ordinary dot product,
  i.e. with complex conjugation.  (In principle, do we also need
  to re-orthogonalize, for the case of degenerate eigenvalues?) */
  {
    int i, one = 1;
    for (i = 0; i < n; ++i)
    {
      Complex norm = 1.0 / csqrt(symmetric_dot(n, V+i*n, V+i*n));
      ZSCAL(&n, &norm, V+i*n, &one);
    }

  }
}

/**************************************************************************/

/* how conservative do we need to be for this? */
#define SINGULAR_THRESHOLD 1e-5

/*
   Solve the eigenvalue problem U1 b = u U0 b, where b is the eigenvector
   and u is the eigenvalue.  u = exp(iwt - at) then contains both the
   frequency and the decay constant.
*/
void harminv_solve_once(harminv_data d)
{
  int J, i, one=1;
  Complex zone = 1.0, zzero = 0.0;
  Complex *V0, *v0, *H1, *V1; /* for eigensolutions of U0 and U1 */
  double max_v0 = 0.0;

  J = d->J;
  CHK_MALLOC(V0, Complex, J*J);
  CHK_MALLOC(v0, Complex, J);

/*
  Unfortunately, U0 is very likely to be singular, so we must
  first extract the non-singular eigenvectors and only work in
  that sub-space.  See the Wall & Neuhauser paper.
*/

  solve_eigenvects(J, d->U0, V0, v0);

  /* find maximum |eigenvalue| */
  for (i = 0; i < J; ++i)
  {
    double v = cabs(v0[i]);
    if (v > max_v0)
      max_v0 = v;
  }

/*
  we must remove the singular components of U0, those
  that are less than some threshold times the maximum eigenvalue.
  Also, we need to scale the eigenvectors by 1/sqrt(eigenval).
*/
  d->nfreqs = J;
  for (i = 0; i < J; ++i)
  {
    if (cabs(v0[i]) < SINGULAR_THRESHOLD * max_v0)
    {
      v0[i] = 0; /* tag as a "hole" */
      d->nfreqs -= 1;
    }
    else
    {
      // not singular
      Complex s;
      int j;
      // move the eigenvector to the first "hole" left by
      // deleting singular eigenvalues:
      for (j = 0; j < i && v0[j] != 0.0; ++j);
      if (j < i)
      {
        ZCOPY(&J, V0 + i*J, &one, V0 + j*J, &one);
        v0[j] = v0[i];
        v0[i] = 0; /* tag as a "hole" */
      }
      s = 1.0 / csqrt(v0[j]);
      ZSCAL(&J, &s, V0 + j*J, &one);
    }
  }

  CHK_MALLOC(d->B, Complex, d->nfreqs * J);
  CHK_MALLOC(d->u, Complex, d->nfreqs);
  CHK_MALLOC(V1, Complex, d->nfreqs * d->nfreqs);
  CHK_MALLOC(H1, Complex, d->nfreqs * d->nfreqs);

  /* compute H1 = V0 * U1 * V0': */

  /* B = V0 * U1: */
  ZGEMM(&cn, &cn, &J, &d->nfreqs, &J,
        &zone, d->U1, &J, V0, &J, &zzero, d->B, &J);
  /* H1 = B * transpose(V0) */
  ZGEMM(&ct, &cn, &d->nfreqs, &d->nfreqs, &J,
        &zone, V0, &J, d->B, &J, &zzero, H1, &d->nfreqs);

  /* Finally, we can find the eigenvalues and eigenvectors: */
  solve_eigenvects(d->nfreqs, H1, V1, d->u);
  /* B = V1 * V0: */
  ZGEMM(&cn, &cn, &J, &d->nfreqs, &d->nfreqs,
        &zone, V0, &J, V1, &d->nfreqs, &zzero, d->B, &J);

  GC_free(H1);
  GC_free(V1);
  GC_free(v0);
  GC_free(V0);
}

/**************************************************************************/

/* After solving once, solve again using the solutions from last
   time as the input to the spectra estimator this time.

   Optionally, if mode_ok is not NULL, we can only use the solutions k from
   last time that pass ok(d, k, ok_d).  Currently, this is not recommended
   as it seems to make things worse.  */
void harminv_solve_again(harminv_data d, harminv_mode_ok_func ok, void *ok_d)
{
  int i, j;
  char *mode_ok = 0;
  CHECK(d->nfreqs >= 0, "haven't computed eigensolutions yet");

  if (!d->nfreqs) return; /* no eigensolutions to work with */

  if (ok)
  {
    CHK_MALLOC(mode_ok, char, d->nfreqs);
    ok(d, -1, ok_d); /* initialize */
    for (i = 0; i < d->nfreqs; ++i)
      mode_ok[i] = ok(d, i, ok_d);
  }

  GC_free(d->B);
  GC_free(d->U1); GC_free(d->U0); GC_free(d->G0); GC_free(d->G0_M); GC_free(d->D0);
  GC_free(d->z);
  GC_free(d->amps);
  GC_free(d->errs);
  d->B = d->U1 = d->U0 = d->z = d->amps = (Complex *) NULL;
  d->G0 = d->G0_M = d->D0 = (Complex *) NULL;
  d->errs = (double *) NULL;

  /* Spectral grid needs to be on the unit circle or system is unstable: */
  for (i = j = 0; i < d->nfreqs; ++i)
  {
    if (!ok || mode_ok[i])
      d->u[j++] = d->u[i] / cabs(d->u[i]);
  }
  d->nfreqs = j;

  if (ok)
  {
    ok(d, -2, ok_d); /* finish */
    GC_free(mode_ok);
  }

  d->u = (Complex *) GC_realloc(d->u, sizeof(Complex) * d->nfreqs);

  if (!d->nfreqs) return; /* no eigensolutions to work with */

  init_z(d, d->nfreqs, d->u);

  d->nfreqs = 0;
  d->B = d->u = NULL;

  harminv_solve_once(d);
}

/**************************************************************************/

/* Keep re-solving as long as spurious solutions are eliminated.

   Currently, it is recommended that you use harminv_solve (i.e. pass ok = 0);
   see harminv_solve_again. */
void harminv_solve_ok_modes(harminv_data d, harminv_mode_ok_func ok,void *ok_d)
{

  harminv_solve_once(d);

/*
  This is not in the papers, but seems to be a good idea:
  plug the u's back in as z's for another pass, and repeat
  as long as the number of eigenvalues decreases.   Effectively,
  this gives us more basis functions where the modes are.
*/
  {
    int prev_nf, cur_nf, nf_ok;
    cur_nf = harminv_get_num_freqs(d);
    do
    {
      prev_nf = cur_nf;
      harminv_solve_again(d, ok, ok_d);
      cur_nf = harminv_get_num_freqs(d);
      if (ok)
      {
        ok(d, -1, ok_d); /* initialize */
        for (nf_ok = 0; nf_ok < cur_nf
             && ok(d, nf_ok, ok_d); ++nf_ok);
        ok(d, -2, ok_d); /* finish */
      }
      else
        nf_ok = cur_nf;
    }
    while (cur_nf < prev_nf || nf_ok < cur_nf);
    /* FIXME: solve one more time for good measure? */
  }
}

/**************************************************************************/

void harminv_solve(harminv_data d)
{
  harminv_solve_ok_modes(d, NULL, NULL);
}

/**************************************************************************/

#define NORMSQR(c) (creal(c) * creal(c) + cimag(c) * cimag(c))

/* Returns an array (of size harminv_get_num_freqs(d)) of estimates
   for the |error| in the solution frequencies.  Solutions with
   errors much larger than the smallest error are likely to be spurious. */
double *harminv_compute_frequency_errors(harminv_data d)
{
     int i, J2, one = 1;
     Complex *U2, *U2b;
     double *freq_err;

     CHECK(d->nfreqs >= 0, "haven't computed eigensolutions yet");
     if (!d->nfreqs) return NULL;

     CHK_MALLOC(freq_err, double, d->nfreqs);

     J2 = d->J*d->J;
     CHK_MALLOC(U2, Complex, J2);
     generate_U(U2, NULL, 2, d->c, d->n, d->K, d->J, d->J, d->z, d->z,
    NULL, NULL, NULL);
     CHK_MALLOC(U2b, Complex, d->J);

     /* For each eigenstate, compute an estimate of the error, roughly
  as suggested in W&N, eq. (2.19). */

     for (i = 0; i < d->nfreqs; ++i) {
    Complex zone = 1.0, zzero = 0.0;

    /* compute U2b = U2 * B[i] */
    ZGEMV(&ct, &d->J, &d->J,
    &zone, U2, &d->J, d->B + i * d->J, &one,
    &zzero, U2b, &one);

    /* ideally, B[i] should satisfy U2 B[i] = u^2 U0 B[i].
       since B U0 B = 1, then we can get a second estimate
       for u by sqrt(B[i] U2 B[i]), and from this we compute
       the relative error in the (complex) frequency. */

    freq_err[i] =
         cabs(clog(csqrt(symmetric_dot(d->J, d->B + i * d->J, U2b))
       / d->u[i])) / cabs(clog(d->u[i]));
     }

     GC_free(U2b);
     GC_free(U2);
     return freq_err;
}

/**************************************************************************/

#define UNITY_THRESH 1e-4 /* FIXME? */

/* true if UNITY_THRESH < |u|^n < 1/UNITY_THRESH. */
static int u_near_unity(Complex u, int n)
{
     double nlgabsu = n * log(cabs(u));
     return (log(UNITY_THRESH) < nlgabsu && nlgabsu < -log(UNITY_THRESH));
}

/* Return an array (of size harminv_get_num_freqs(d)) of complex
   amplitudes of each sinusoid in the solution. */
Complex *harminv_compute_amplitudes(harminv_data d)
{
     int k, j;
     Complex *u;
     Complex *Uu;
     Complex *a; /* the amplitudes of the eigenfrequencies */
     int ku, nu;

     CHECK(d->nfreqs >= 0, "haven't computed eigensolutions yet");
     if (!d->nfreqs) return NULL;

     CHK_MALLOC(a, Complex, d->nfreqs);
     CHK_MALLOC(u, Complex, d->nfreqs);

     for (k = ku = 0; k < d->nfreqs; ++k)
    if (u_near_unity(d->u[k], d->n))
         u[ku++] = d->u[k];
     nu = ku;

     CHK_MALLOC(Uu, Complex, d->J * nu);
     generate_U(Uu, NULL, 0, d->c, d->n, d->K, d->J, nu, d->z, u,
    &d->G0, &d->G0_M, &d->D0);

     /* compute the amplitudes via eq. 27 of M&T, except when |u| is
  too small..in that case, the computation of Uu is unstable,
  and we use eq. 26 instead (which doesn't use half of the data,
  but doesn't blow up either): */
     for (k = ku = 0; k < d->nfreqs; ++k) {
    Complex asum = 0;
    if (u_near_unity(d->u[k], d->n)) { /* eq. 27 */
         for (j = 0; j < d->J; ++j)
        asum += d->B[k * d->J + j] * Uu[j * nu + ku];
         asum /= d->K;
         ku++;
    }
    else { /* eq. 26 */
         for (j = 0; j < d->J; ++j)
        asum += d->B[k * d->J + j] * d->G0[j];
    }
    a[k] = asum * asum;
     }

     GC_free(Uu);
     GC_free(u);
     return a;
}

/**************************************************************************/

double harminv_get_freq(harminv_data d, int k)
{
     CHECK(d->nfreqs >= 0, "haven't computed eigensolutions yet");
     CHECK(k >= 0 && k < d->nfreqs,
     "argument out of range in harminv_get_freq");
     return(-carg(d->u[k]) / TWOPI);
}

double harminv_get_decay(harminv_data d, int k)
{
     CHECK(d->nfreqs >= 0, "haven't computed eigensolutions yet");
     CHECK(k >= 0 && k < d->nfreqs,
     "argument out of range in harminv_get_decay");
     return(-log(cabs(d->u[k])));
}

double harminv_get_Q(harminv_data d, int k)
{
     CHECK(k >= 0 && k < d->nfreqs,
     "argument out of range in harminv_get_Q");
     return(TWOPI * fabs(harminv_get_freq(d, k))
      / (2 * harminv_get_decay(d, k)));
}

Complex harminv_get_omega(harminv_data d, int k)
{
     CHECK(d->nfreqs >= 0, "haven't computed eigensolutions yet");
     CHECK(k >= 0 && k < d->nfreqs,
     "argument out of range in harminv_get_omega");
     return(I * clog(d->u[k]));
}

Complex harminv_get_amplitude(harminv_data d, int k)
{
     CHECK(k >= 0 && k < d->nfreqs,
     "argument out of range in harminv_get_amplitude");
     if (!d->amps)
    d->amps = harminv_compute_amplitudes(d);
     return d->amps[k];
}

double harminv_get_freq_error(harminv_data d, int k)
{
     CHECK(k >= 0 && k < d->nfreqs,
     "argument out of range in harminv_get_freq_error");
     if (!d->errs)
    d->errs = harminv_compute_frequency_errors(d);
     return d->errs[k];
}

/**************************************************************************/

static int mode_ok(harminv_data d, int k, void *ok_d_)
{
  mode_ok_data *ok_d = (mode_ok_data *) ok_d_;
  double errk, ampk, f;
  int ok;

  if (k == -1) { /* initialize */
    int i;
    ok_d->num_ok = 0;
    if (!harminv_get_num_freqs(d))
      return 0;
    ok_d->min_err = harminv_get_freq_error(d, 0);;
    ok_d->max_amp = cabs( harminv_get_amplitude( d, 0));
    for (i = 1; i < harminv_get_num_freqs(d); ++i) {
      double err, amp;
      if ((err = harminv_get_freq_error(d, i)) < ok_d->min_err)
        ok_d->min_err = err;
      if ((amp = cabs(harminv_get_amplitude(d, i))) > ok_d->max_amp)
        ok_d->max_amp = amp;

    }
    return 0;
  }
  else if (k == -2)
  { /* finish */
    if (ok_d->verbose && harminv_get_num_freqs(d))
      printf("# harminv: %d/%d modes are ok: "
          "errs <= %e and %e * %e\n, "
          "amps >= %g, %e * %g, "
          "|Q| >= %g\n",
          ok_d->num_ok, harminv_get_num_freqs(d),
                                              ok_d->err_thresh, ok_d->rel_err_thresh, ok_d->min_err,
                                              ok_d->amp_thresh, ok_d->rel_amp_thresh, ok_d->max_amp,
                                              ok_d->Q_thresh);
    return 0;
  }

  f = fabs(harminv_get_freq(d, k));
  errk = harminv_get_freq_error(d, k);
  ampk = cabs(harminv_get_amplitude(d, k));

  ok = ((!ok_d->only_f_inrange || (f >= ok_d->fmin && f <= ok_d->fmax))
      && errk <= ok_d->err_thresh
      && errk <= ok_d->min_err * ok_d->rel_err_thresh
      && ampk >= ok_d->amp_thresh
      && ampk >= ok_d->rel_amp_thresh * ok_d->max_amp
      && fabs(harminv_get_Q(d,k)) >= ok_d->Q_thresh);

  ok_d->num_ok += ok;

  return ok;
}



