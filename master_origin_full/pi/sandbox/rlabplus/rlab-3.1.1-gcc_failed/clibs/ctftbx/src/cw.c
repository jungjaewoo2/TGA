/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRCW.C                           *
 *============================================================================*
 * Name of the function : cw.c (void)                                         *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Choi-Williams*
 * Distribution (CW) :                                                        *
 *                                                                            *
 *                 //                                                         *
 *                 ||     sqrt(sigma)     -mu^2 *sigma / (16 tau^2)           *
 *      CW(t,f) =  || -----------------  e                           .../... *
 *                 ||  2 sqrt(pi)|tau |                                       *
 *                //                                                          *
 *                                                                            *
 *                                                 -j2pi f tau                *
 *                     x(t+mu+tau/2)x*(t+mu-tau/2)e            dmu dtau       *
 *                                                                            *
 *                                                                            *
 * This function is real valued. Its computation requires a real or complex   *
 * signal, a vector containing time instants, the number of frequency bins, a *
 * time smoothing window, a frequency smoothing window and the kernel width.  *
 *                                                                            *
 *============================================================================*
 * INPUT VARIABLES                                                            *
 * Name                |              role                                    *
 * Signal              | The signal to analyze. No field modified             *
 *                     |                                                      *
 * WindowT             | Vector containing the points of the time moothing    *
 *                     | window                                               *
 * WindowT_Length      | Number of points of the time window (ODD number !)   *
 *                     |                                                      *
 * WindowF             | Vector containing the points of the frequency window *
 * WindowF_Length      | Number of points of the window (ODD number !)        *
 *                     |                                                      *
 * tfr                 | Matrix containing the resulting TFR (real)           *
 * tfr.time_instants   | positions of the smoothing window                    *
 * tfr.N_time          | length of '.time_instants' = number of cols.         *
 *                     | in the tfr matrix                                    *
 * tfr.N_freq          | number of frequency bins = number of rows in the tfr *
 *                     | matrix                                               *
 * tfr.is_complex      | must be set to FALSE (a CW tfr is real-valued)       *
 *                     |                                                      *
 * sigma               | the kernel width                                     *
 *                     |                                                      *
 *----------------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                           *
 * Name                |                role                                  *
 * tfr.real_part       | the output tfr matrix  (real_part)                   *
 * tfr.freq_bins       | vector of frequency bins (freqs where the tfr matrix *
 *                     | is computed)                                         *
 *----------------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                         *
 * Name                |                 role                                 *
 *                     |                                                      *
 * Nfft                | Next power of two to tfr.N_freq                      *
 * column, row         | variables of displacement in the matrices            *
 * time                | local time-instant variable to compute the tfr       *
 *                     |                                                      *
 * half_WindowT_Length | half-length of the time smoothing window             *
 *                     |                                                      *
 * half_WindowF_Length | half-length of the frequency smoothing window        *
 * normF               | normalization factor for the frequency window        *
 *                     |                                                      *
 * CWKernel            | variable to compute the Choi-Williams Kernel         *
 * normK               | normalization factor for the Kernel                  *
 * spreadfac           | factor = 16/sigma                                    *
 * index               | variable to locate position in the kernel matrix     *
 *                     |                                                      *
 * tau                 | time-lag variable                                    *
 * taumin              | local time-lag variable bounds. Used to take into    *
 * taumax              | accound the beginning and the end of the             *
 *                     | signal, where the window is cut                      *
 *                     |                                                      *
 * mu                  | time-smoothing variable                              *
 * mumin               | local time-smoothing variable bounds. Used to take   *
 * mumax               | into accound the beginning and the end of time       *
 *                     | smoothing procedure                                  *
 *                     |                                                      *
 * lacf_real           | real and imaginary parts of the local autocorrelation*
 * lacf_imag           | function of the signal                               *
 *                     |                                                      *
 * R1_real R1_imag     | used to compute real and imaginary parts of the time *
 * R2_real R2_imag     | smoothed-windowed local autocorrelation function     *
 *                     |                                                      *
 *============================================================================*
 * SUBROUTINES USED HERE                                                      *
 *----------------------------------------------------------------------------*
 * Name   | int idx(int i_row, int j_col, int nb_row)                         *
 * Action | computes the vector index for an element in a matrix given the row*
 *        | and column indices (i,j) and the total number of row              *
 * Place  | divers.c                                                          *
 *----------------------------------------------------------------------------*
 * Name   | void fft(int n, int m, double *x, double *y)                      *
 * Action | Computes the fft                                                  *
 * Place  | divers.c                                                          *
 *----------------------------------------------------------------------------*
 * Name   | int po2(int x)                                                    *
 * Action | Computes the next power of two of x                               *
 * Place  | divers.c                                                          *
 *============================================================================*/

int
ctftbx_cw (type_signal Signal, double *WindowT, int WindowT_Length,
           double *WindowF, int WindowF_Length, double sigma, type_TFR tfr)
{
  int            Nfft, column, row, time;
  int            half_WindowT_Length, half_WindowF_Length;
  int            taumin, taumax, tau;
  int            mumin, mumax, mu, index;
  Complex       *lacf=0;/* local autocorrelation function */
  double         normK, normF, spreadfac;
  Complex        R1, R2;
  double        *CWKernel=0;

 /*--------------------------------------------------------------------------*/
 /*                        Test the input variables                          */
 /*--------------------------------------------------------------------------*/


   if (tfr.is_complex == TRUE)
   {
     if (ctftbx_debug)
       printf ("cw.c : The tfr matrix must be real valued\n");
     return 1;
   }

   if (tfr.N_freq <= 0)
   {
     if (ctftbx_debug)
       printf ("cw.c : The field tfr.N_freq is not correctly set\n");
     return 1;
   }

   if (tfr.N_time <= 0)
   {
     if (ctftbx_debug)
       printf ("cw.c : The field tfr.N_time is not correctly set\n");
     return 1;
   }

   if (ISODD(WindowT_Length) == 0)
   {
     if (ctftbx_debug)
       printf ("cw.c : The time-window Length must be an ODD number\n");
     return 1;
   }

   if (ISODD(WindowF_Length) == 0)
   {
     if (ctftbx_debug)
       printf ("cw.c : The frequency-window Length must be an ODD number\n");
     return 1;
   }

 /*--------------------------------------------------------------------------*/
 /*                    Determines some internal constants                    */
 /*--------------------------------------------------------------------------*/

  half_WindowT_Length = (WindowT_Length - 1) / 2;

  half_WindowF_Length = (WindowF_Length - 1) / 2;
  normF=WindowF[half_WindowF_Length];

  /* normalization of the frequency smoothing window */
  for(row = 0; row < WindowF_Length; row++)
  {
    WindowF[row]=WindowF[row]/normF;
  }

 /*--------------------------------------------------------------------------*/
 /*          Memory allocation and computation of  the kernel                */
 /*--------------------------------------------------------------------------*/
  CWKernel = (double * ) GC_malloc (MIN(tfr.N_freq,half_WindowF_Length)
			                                    * WindowT_Length * sizeof(double));

  spreadfac = 16.0/sigma;
  taumax=MIN(tfr.N_freq,half_WindowF_Length);

  for(tau=1;tau<=taumax;tau++)
  {
    for(mu=-half_WindowT_Length;mu<=+half_WindowT_Length;mu++)
    {
      CWKernel[idx(tau-1,half_WindowT_Length+mu,taumax)] =
          exp(-1.0/(spreadfac * tau * tau) * mu * mu)
          * WindowT[half_WindowT_Length + mu];
    }
  }

 /*--------------------------------------------------------------------------*/
 /*           creation of the vector of frequency bins  (output)             */
 /*--------------------------------------------------------------------------*/
  Nfft = po2 (tfr.N_freq);

  for (row = 0; row < tfr.N_freq; row++)
  {
    tfr.freq_bins[row] = (double) (0.5 * row) / tfr.N_freq;
  }

 /*--------------------------------------------------------------------------*/
 /*     memory allocation and init. of the local autocorrelation fuction     */
 /*--------------------------------------------------------------------------*/
  lacf = (Complex *) GC_malloc (tfr.N_freq * sizeof (Complex));

 /* initialization of the intermediary vectors */
  for (row = 0; row < tfr.N_freq ; row++)
  {
    lacf[row] = 0.0;
  }

 /*--------------------------------------------------------------------------*/
 /*      computation of the fft for the local autocorrelation function       */
 /*--------------------------------------------------------------------------*/
  for (column = 0; column < tfr.N_time; column++)
  {
    /* time instants of interest to compute the tfr */
    time = ((int) tfr.time_instants[column]) - 1;

      /* maximum value of the delay in order to take the edges
    into account */
    taumax = MIN((time+half_WindowT_Length),
                  (Signal.length-time-1+half_WindowT_Length));
    taumax = MIN(taumax,(tfr.N_freq / 2 - 1));
    taumax = MIN(taumax, half_WindowF_Length);

    if (Signal.is_complex == TRUE)
    {
      lacf[0] = RE(Signal.signal[time]) *  RE(Signal.signal[time])
          +  IM(Signal.signal[time]) * IM(Signal.signal[time]);
    }
    else /* the signal is real-valued */
    {
      lacf[0] = RE(Signal.signal[time]) * RE(Signal.signal[time]);
    }

    /* The signal is windowed around the current time */
    for (tau = 1; tau <= taumax; tau++)
    {
      R1=0.0;
      R2=0.0;

      /* bound of mu in order to take into account the edges */
      mumin=MIN(half_WindowT_Length, (Signal.length-time-1-tau));
      mumax=MIN(half_WindowT_Length,time-tau);

      normK=0;
      for(mu=-mumin;mu<=mumax;mu++)
      {
        normK = normK +
            CWKernel[idx(tau-1,half_WindowT_Length+mu,
                         MIN(tfr.N_freq / 2,half_WindowF_Length))];
      }

      for(mu=-mumin;mu<=mumax;mu++)
      {
        /* case of complex valued signal */
        if (Signal.is_complex == TRUE)
        {
          index = idx(tau-1, half_WindowT_Length+mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R1 = R1 + Signal.signal[time + tau - mu] * conj(Signal.signal[time - tau - mu])
              * CWKernel[index]/normK;

          index = idx(tau-1,half_WindowT_Length-mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R2 = R2 + Signal.signal[time - tau - mu] * conj(Signal.signal[time + tau - mu])
              * CWKernel[index]/normK;
        }
        /* case of real-valued signal */
        else
        {
          index = idx(tau-1, half_WindowT_Length+mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R1 = R1 + RE(Signal.signal[time + tau - mu])
              * RE(Signal.signal[time - tau - mu])
              * CWKernel[index]/normK;

          index = idx(tau-1, half_WindowT_Length-mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R2 = R2 + RE(Signal.signal[time - tau - mu])
              * RE(Signal.signal[time + tau - mu])
              * CWKernel[index]/normK;
        }
      }

      lacf[tau] = R1 * WindowF[half_WindowF_Length + tau];

      lacf[tfr.N_freq - tau] = R2 * WindowF[half_WindowF_Length - tau];
    }



    tau=floor(tfr.N_freq/2);
    if ((time<=Signal.length-tau-1)&(time>=tau)&(tau<=half_WindowF_Length))
    {
      R1=0.0;
      R2=0.0;

      /* bound of mu in order to take into account the edges */
      mumin=MIN(half_WindowT_Length, (Signal.length-time-1-tau));
      mumax=MIN(half_WindowT_Length,time-tau);

      normK=0;
      for(mu=-mumin;mu<=mumax;mu++)
      {
        normK = normK + CWKernel[idx(tau-1,half_WindowT_Length+mu,
                                     MIN(tfr.N_freq / 2,half_WindowF_Length))];
      }

      for(mu=-mumin;mu<=mumax;mu++)
      {
        /* case of complex valued signal */
        if (Signal.is_complex == TRUE)
        {
          index = idx(tau-1, half_WindowT_Length+mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R1 = R1 + Signal.signal[time + tau - mu] * conj(Signal.signal[time - tau - mu])
              * CWKernel[index]/normK;

          index = idx(tau-1, half_WindowT_Length-mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R2 = R2 + Signal.signal[time - tau - mu] * conj(Signal.signal[time + tau - mu])
              * CWKernel[index]/normK;
        }
        /* case of real-valued signal */
        else
        {
          index = idx(tau-1, half_WindowT_Length+mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R1 = R1 + RE(Signal.signal[time + tau - mu])
              * RE(Signal.signal[time - tau - mu])
              * CWKernel[index]/normK;

          index = idx(tau-1, half_WindowT_Length-mu,
                      MIN(tfr.N_freq / 2,half_WindowF_Length));

          R2 = R2 + RE(Signal.signal[time - tau - mu])
              * RE(Signal.signal[time + tau - mu])
              * CWKernel[index]/normK;
        }
      }


      lacf[tau] = 0.5*(R1 * WindowF[half_WindowF_Length + tau]
          + R2 * WindowF[half_WindowF_Length - tau]);

    }

    /* fft of the local autocorrelation function lacf */
    fft (tfr.N_freq, Nfft, lacf);

    /* the fft is put in the tfr matrix  */
    for (row = 0; row < tfr.N_freq; row++)
    {
      tfr.tfr[idx (row,column,tfr.N_freq)] = RE(lacf[row]);
      lacf[row] = 0.0;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*                free the memory used in this program                      */
  /*--------------------------------------------------------------------------*/
  GC_free (lacf);
  GC_free (CWKernel);

  return 0;
}
