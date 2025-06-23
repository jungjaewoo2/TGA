static int  mem_alloc_TFR (type_TFR *tfr, double *ptr_freq_bins,
                    double *ptr_time_instants, Complex *ptr_tfr);
static void mem_free_TFR (type_TFR *tfr);
static int  mem_alloc_AF (type_AF *ambig_func, double *ptr_doppler_bins,
                   double *ptr_delay_bins,  Complex *ptr_af);
static void mem_free_AF (type_AF *ambig_func);
static void transpose (int N_line, int N_col, Complex *matrix);
static double powof (double x, double alpha);
static int po2 (int n);
static int fft (int siglen, int Nfft, Complex *signal);
static int ifft (int siglen, int Nfft, Complex *signal);
static double sinc (double x);
static int irem( double x, double y);
static int idx (int i_row, int j_col, int nb_row);
static double sqr (double x);
static int double_fftshift (double *vector_in, double *vector_out, int vector_length);
static int Complex_fftshift (Complex *vector_in, Complex *vector_out, int vector_length);


/*====================================================================*
 * this file contains several programs used in a signal processing    *
 * context. The programs are :                                        *
 * - fft       : Computes the fast Fourier transform of a complex     *
 *               signal                                               *
 * - po2       : Computes the next higher power of two                *
 * - idx       : In a matrix stored like a vector, computes the       *
 *               vector index corresponding to given line and column  *
 *               numbers in the matrix                                *
 * - ifft      : Computes the inverse fft                             *
 * - sqr       : square of a number                                   *
 * - powof     : Computes x to the power of alpha                     *
 * - Recover_Signal : In programs used inside the Matlab              *
 *                environement, recovers a signal from a matlab       *
 *                function                                            *
 * - gradient  : Computes the bidimensional gradient of a surface     *
 *               (called 'field of potential')                        *
 *====================================================================*/


/*====================================================================*
 * Name of the function : powof (double)                               *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     : 02 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *                                                                    *
 *  Computes x^alpha and takes the case x=0 into account              *
 *  !!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!              *
 *  x must be non-negative when alpha is not an interger              *
 *  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              *
 *====================================================================*/

 /*===================================================================*
 * Name of the function : transpose                 *
 * Author               : Manuel Davy                                 *
 * Date of creation     : 10 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM                                  *
 *                      *
 *                      *
 *====================================================================*
 * INPUT VARIABLES                          *
 * Name        | type  |              role                      *
 *             |       | kernel               *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                     *
 * Name        | type  |              role                        *
 *                      *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES                   *
 * Name        | type  |              role                        *
 *                      *
 *====================================================================*
 * SUBROUTINES USED HERE                      *
 *--------------------------------------------------------------------*
 * Name   |                                                       *
 * Action |                   *
 * Place  |                     *
 *====================================================================*/

static void
transpose (int N_line, int N_col, Complex *matrix)
{
  int            line, col, index;

  /* checks if the the matrix is not reduced to a single element */
  {
    if ((N_line > 1) && (N_col > 1))
      if (N_line == N_col)  /* the matrix is square */
    {
      /* requires an intermediary element */
      Complex        inter;
      int            index_1, index_2;

      for (line = 1; line < N_line; line++)
      {
        for (col = 0; col < line; col++)
        {
          index_1 = idx (line, col, N_line);  /* in the under triangle */
          index_2 = idx (col, line, N_col); /* in the upper triangle */

          inter = matrix[index_1];
          matrix[index_1] = matrix[index_2];
          matrix[index_2] = inter;
        }
      }
    }
    else
      /* the matrix is not square */
    {
      /* requires an intermediary matrix */
      Complex *inter;

      inter = (Complex *) GC_malloc (N_line * N_col * sizeof (Complex));

      /* recopy in a transpose matrix */
      for (line = 0; line < N_line; line++)
      {
        for (col = 0; col < N_col; col++)
        {
          inter[idx (col, line, N_col)] = matrix[idx (line, col, N_line)];
        }
      }
      /* recopy in the original matrix */
      for (index = 0; index < (N_line * N_col); index++)
      {
        matrix[index] = inter[index];
      }

      GC_free (inter);
    }
  }
}

static double
powof (double x, double alpha)
{
  double         resu;
  if (x >= 0) /* in this case, no problem */
  {
    if (x == 0)
      resu = 0.0;
    else
      resu = exp (alpha * log (x));
  }
  else /* there may be problems */
  {
    if (alpha == (int)(alpha)) /* if alpha is an integer */
    {
      /* then x^alpha is real-valued */
      resu = powof ( -x, alpha);
    /* and the sign of the results depends on
      wether alpha is ODD or EVEN */
      if (ISODD(alpha) == 1)
      {
        /* alpha is ODD */
        resu = -resu;
      }
    }
    else
    {
      if (ctftbx_debug)
        printf ("Attempt to compute x^alpha with x<0 : complex valued result\n");
      return 0;
    }
  }
  return resu;
}

// ****************************************************
//
// FFTPACK: find it and build it. link the library
// with -lfftpack
//
// ****************************************************
# define CFFTI dcffti_
# define CFFTF dcfftf_
# define CFFTB dcfftb_

extern int CFFTI ();
extern int CFFTF ();
extern int CFFTB ();

static int fft (int siglen, int Nfft, Complex *signal)
{
  int n, i;
  double  *wsave=0;
  Complex *wsignal=0;

  n = (int) powof (2.0, Nfft) + 1;
  wsave = (double *) GC_malloc ((4 * n + 15) * sizeof(double));
  CFFTI (&n, wsave);

  if (siglen >= n)
  {
    // only part of the signal is fft-ed
    CFFTF (&n, (double *) (&signal[0]), wsave);
  }
  else
  {
    wsignal = (Complex *) GC_malloc (n * sizeof(double));
    // copy signal to wsignal and pad the rest with zeros
    for (i=0; i<n; i++)
    {
      if (i<siglen)
      {
        wsignal[i] = signal[i];
      }
      else
      {
        wsignal[i] = 0.0;
      }
    }
    // fft the padded signal
    CFFTF (&n, (double *) (&wsignal[0]), wsave);
    for (i=0; i<siglen; i++)
    {
      signal[i] = wsignal[i];
    }
    GC_free (wsignal);
  }

  GC_free (wsave);

  return 0;
}

static int ifft (int siglen, int Nfft, Complex *signal)
{
  int n, i;
  double  *wsave=0;
  Complex *wsignal=0;

  n = (int) powof (2.0, Nfft) + 1;
  wsave = (double *) GC_malloc ((4 * n + 15) * sizeof(double));
  CFFTI (&n, wsave);

  if (siglen >= n)
  {
    // only part of the signal is fft-ed
    CFFTB (&n, (double *) (&signal[0]), wsave);
  }
  else
  {
    wsignal = (Complex *) GC_malloc (n * sizeof(double));
    // copy signal to wsignal and pad the rest with zeros
    for (i=0; i<n; i++)
    {
      if (i<siglen)
      {
        wsignal[i] = signal[i];
      }
      else
      {
        wsignal[i] = 0.0;
      }
    }
    // fft the padded signal
    CFFTB (&n, (double *) (&wsignal[0]), wsave);
    for (i=0; i<siglen; i++)
    {
      signal[i] = wsignal[i];
    }
    GC_free (wsignal);
  }

  GC_free (wsave);

  return 0;
}

/*====================================================================*
 * Name of the function : po2 (int)                                   *
 * Author               :                                             *
 * Date of creation     :                                             *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * Computes the next power of to of an integer n, i.e.                *
 * the smallest k such that m=2^k>=n                                  *
 *====================================================================*/
static int po2 (int n)
{
  int            nextpo2, m;

  m = 1;
  nextpo2 = 0;
  while (m < n)
    {
      ++nextpo2;
      m = m + m;
    }

  return (nextpo2);
}

/*====================================================================*
 * Name of the function : sinc (double)                               *
 * Author               : Emmanuel ROY                                *
 * Date of creation     : 22 - 06 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * Computes the sinc function of a number, defined by :               *
 *                                                                    *
 *   for any x in R,                                                  *
 *                                                                    *
 *             /                                                      *
 *             | 1 if x = 0                                           *
 *             |                                                      *
 *   sinc(x) = | sin (pi * x)                                         *
 *             | ------------ if x != 0                               *
 *             |   pi * x                                             *
 *             \                                                      *
 *                                                                    *
 *====================================================================*/
static double sinc (double x)
{
  double         resu;

  if (x == 0)
   	resu = 1.0;
  else
	resu = sin(pi*x)/(pi*x);

  return resu;
}


/*====================================================================*
 * Name of the function : irem (double x, double y)                   *
 * Author               : Emmanuel ROY                                *
 * Date of creation     : 22 - 06 -1999                               *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *                                                                    *
 * Computes the remainder after Euclidean division of two doubles     *
 *                                                                    *
 *====================================================================*/

static int irem( double x, double y)
{
 int result;

 if (y != 0)
 {
   result =  x-y* (int)(x/y);
 }
 else
 {
   result = 0;
   if (ctftbx_debug)
     printf("irem.c : attempt to divide by 0\n");
 }

 return result;
}

/*====================================================================*
 * Name of the function : idx (int i_row, int j_col, int nb_row)      *
 * Author               : Emmanuel ROY                                *
 * Date of creation     :                                             *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *                                                                    *
 * The matrices are stored as vectors, column by column               *
 * This program computes the vector index corresponding to the        *
 * specified line number (line), the column number (col) and the      *
 * number of lines   (nb_line)                                        *
 *====================================================================*/

static int
idx (int i_row, int j_col, int nb_row)
{
  if (i_row >= nb_row)
  {
    if (ctftbx_debug)
      printf("idx : incorrect row number\n");
    return 0;
  }
  else
  {
    return (i_row + (j_col * nb_row));
  }
}


/*====================================================================*
 * Name of the function : sqr (double)                                *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     : 02 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *								      *
 * Computes the square value of x                                     *
 *====================================================================*/

static double
sqr (double x)
{
  return (x * x);
}



/*====================================================================*
 * Name of the function : Recover_Signal (void)                       *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     : 02 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * In Matlab environment, recovers a signal given by a signo.m program*
 * the function signo.m must be :                                     *
 * signal=signo(N,class,number);                                      *
 * N     = length of the signal                                       *
 * class = class of the signal                                        *
 * number = number of the signal in its class                         *
 *====================================================================*/



/*====================================================================*
 * Name of the function : fftshift                                    *
 * Date of creation     : 02 - 06 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * swaps the first and second halves of a vector. Example             *
 * [1 2 3 4 5 6 7 8 9 10 ] becomes [6 7 8 9 10 1 2 3 4 5]             *
 * The parameters to pass are :                          	      *
 *   - the input vector		                            	      *
 *   - the output vector	                            	      *
 *   - its length                                                     *
 * if the length is odd, example [1 2 3 4 5] becomes [4 5 1 2 3]      *
 *====================================================================*/
static int
double_fftshift (double *vector_in, double *vector_out, int vector_length)
{
  double inter1, inter2;
  int i, half_length;


  /* computation of the half length in case of odd or even length */
  half_length = (int) (vector_length/2.0);


  /* case where the length is odd */
  if (ISODD(vector_length)==1)
  {
    inter2=vector_in[half_length];
    for (i=0; i<half_length; i++)
    {
      inter1 = vector_in[i];
      vector_out[i] = vector_in[half_length+i+1];
      vector_out[half_length + i ] = inter1;
    }
    vector_out[vector_length-1]=inter2;
  }
  /* case where the length is even */
  else
  {
    for (i=0; i<half_length; i++)
    {
      inter1 = vector_in[half_length + i ];
      vector_out[half_length + i] = vector_in[i];
      vector_out[i] = inter1;
    }
  }
  /* fftshifting of the vector */

  return 0;
}

static int
Complex_fftshift (Complex *vector_in, Complex *vector_out, int vector_length)
{
  Complex inter1, inter2;
  int i, half_length;


  /* computation of the half length in case of odd or even length */
  half_length = (int) (vector_length/2.0);


  /* case where the length is odd */
  if (ISODD(vector_length)==1)
  {
    inter2=vector_in[half_length];
    for (i=0; i<half_length; i++)
    {
      inter1 = vector_in[i];
      vector_out[i] = vector_in[half_length+i+1];
      vector_out[half_length + i ] = inter1;
    }
    vector_out[vector_length-1]=inter2;
  }
  /* case where the length is even */
  else
  {
    for (i=0; i<half_length; i++)
    {
      inter1 = vector_in[half_length + i ];
      vector_out[half_length + i] = vector_in[i];
      vector_out[i] = inter1;
    }
  }
  /* fftshifting of the vector */

  return 0;
}

/*====================================================================*
 * Name of the function : mem_alloc_signal                            *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * memory allocation for the type 'type_signal'                       *
 * the fields :                                                       *
 *   - length		                                 	      *
 *   - is_complex		                            	      *
 * must be previously initialized                                     *
 *====================================================================*/


/*====================================================================*
 * Name of the function : mem_alloc_AF                                *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * memory allocation for the type 'type_AF'                           *
 * the fields :                                                       *
 *   - N_doppler	                                 	      *
 *   - N_delay  	                                 	      *
 *   - is_complex		                            	      *
 * must be previously initialized                                     *
 *====================================================================*/

static int
mem_alloc_AF (type_AF *ambig_func, double *ptr_doppler_bins,
              double *ptr_delay_bins, Complex *ptr_af)
{

  /* some tests to make sure that all the fields are initialized */
  if ((*ambig_func).N_doppler <= 0)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_AF : AF.N_doppler incorrect\n");
    return 1;
  }

  if ((*ambig_func).N_delay <= 0)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_AF : AF.N_delay incorrect\n");
    return 1;
  }

  if (((*ambig_func).is_complex != TRUE) && ((*ambig_func).is_complex != FALSE))
  {
    if (ctftbx_debug)
      printf ("mem_alloc_AF : AF.is_complex incorrect\n");
    return 1;
  }


  /* memory allocation for the field 'doppler_bins' */
  if (ptr_doppler_bins == NULL)
  {
    (*ambig_func).doppler_bins = (double *) GC_malloc
        ((*ambig_func).N_doppler * (*ambig_func).N_delay * sizeof (double));

  }
  else
  {
    (*ambig_func).doppler_bins = ptr_doppler_bins;
  }

  if ((*ambig_func).doppler_bins == NULL)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_AF : memory allocation error\n");
    return 1;
  }

  /* memory allocation for the field 'delay_bins' */
  if (ptr_delay_bins == NULL)
  {
    (*ambig_func).delay_bins = (double *) GC_malloc ((*ambig_func).N_doppler
        * (*ambig_func).N_delay * sizeof (double));

  }
  else
  {
    (*ambig_func).delay_bins = ptr_delay_bins;
  }

  if ((*ambig_func).delay_bins == NULL)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_AF : memory allocation error\n");
    return 1;
  }


  /* memory allocation for the field 'af' */

  if (ptr_af == NULL)
  {
    (*ambig_func).af = (Complex *) GC_malloc ((*ambig_func).N_doppler *
        (*ambig_func).N_delay * sizeof (Complex));
  }
  else
  {
    (*ambig_func).af = ptr_af;
  }
}

/*====================================================================*
 * Name of the function : mem_free_AF                                 *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * free memory for variables of the type 'type_AF'                    *
 *====================================================================*/

static void
mem_free_AF (type_AF *ambig_func)
{
  GC_free ((*ambig_func).doppler_bins);
  GC_free ((*ambig_func).delay_bins);
  GC_free ((*ambig_func).af);

  // set pointers to zero to be on the safe side
  (*ambig_func).doppler_bins = 0;
  (*ambig_func).delay_bins = 0;
  (*ambig_func).af = 0;
}


/*====================================================================*
 * Name of the function : mem_alloc_TFR                               *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * memory allocation for the type 'type_TFR'                          *
 * the fields :                                                       *
 *   - N_freq   	                                 	      *
 *   - N_time   	                                 	      *
 *   - is_complex		                            	      *
 * must be previously initialized                                     *
 *====================================================================*/

static int
mem_alloc_TFR (type_TFR *tfr,
	       double *ptr_freq_bins, double *ptr_time_instants,
	       Complex *ptr_tfr)
{

  /* some tests to make sure that all the fields are initialized */
  if ((*tfr).N_freq <= 0)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_TFR : TFR.N_freq incorrect\n");
    return 1;
  }

  if ((*tfr).N_time <= 0)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_TFR : TFR.N_time incorrect\n");
    return 1;
  }

  /* memory allocation for the field 'freq_bins' */
  if (ptr_freq_bins == NULL)
  {
    (*tfr).freq_bins = (double *) GC_malloc ((*tfr).N_freq * (*tfr).N_time *
    sizeof (double));
  }
  else
  {
    (*tfr).freq_bins = ptr_freq_bins;
  }

  if ((*tfr).freq_bins == NULL)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_TFR : memory allocation error\n");
    return 1;
  }

  /* memory allocation for the field 'time_instants' */
  if (ptr_time_instants == NULL)
  {
    (*tfr).time_instants = (double *) GC_malloc ((*tfr).N_freq *
        (*tfr).N_time * sizeof (double));
  }
  else
  {
    (*tfr).time_instants = ptr_time_instants;
  }

  if ((*tfr).time_instants == NULL)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_TFR : memory allocation error\n");
    return 1;
  }


  /* memory allocation for the field 'real_part' */
  if (ptr_tfr == NULL)
  {
    (*tfr).tfr = (Complex *) GC_malloc ((*tfr).N_freq * (*tfr).N_time *
    sizeof (Complex));
  }
  else
  {
    (*tfr).tfr = ptr_tfr;
  }

  if ((*tfr).tfr == NULL)
  {
    if (ctftbx_debug)
      printf ("mem_alloc_TFR : memory allocation error\n");
    return 1;
  }

}

/*====================================================================*
 * Name of the function : mem_free_TFR                                *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * free memory for variables of the type 'type_TFR'                   *
 *====================================================================*/

static void
mem_free_TFR (type_TFR *tfr)
{
  GC_free ((*tfr).freq_bins);
  GC_free ((*tfr).time_instants);
  GC_free ((*tfr).tfr);

  // set pointers to zero to be on the safe side
  (*tfr).freq_bins = 0;
  (*tfr).time_instants = 0;
  (*tfr).tfr = 0;
}
