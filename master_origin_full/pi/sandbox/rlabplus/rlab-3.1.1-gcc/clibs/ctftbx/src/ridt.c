/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRRIDT.C                         *
 *============================================================================*
 * Name of the function : ridt.c (void)                                       *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Reduced      *
 * Interference Distribution with Triangular kernel (RIDT) :                  *
 *                                                                            *
 *                    /                                                       *
 *                    |                  -j2pi f tau                          *
 *        RIDT(t,f) = | h(tau) R(t,tau) e            dtau                     *
 *                    |                                                       *
 *                   /                                                        *
 *                                                                            *
 *                    / +|tau|/2                                              *
 *                    |            g(mu) /    2 |mu | \                       *
 * with    R(t,tau) = |          2 ----- |1 - -------  |    .../...           *
 *                    |            |tau| \     |tau | /                       *
 *                   / -|tau|/2                                               *
 *                                                                            *
 *                              x(t+mu+tau/2)x*(t+mu-tau/2) dmu               *
 *                                                                            *
 * This function is real valued. Its computation requires a real or complex   *
 * signal, a vector containing time instants, the number of frequency bins, a *
 * time smoothing window and a frequency smoothing window.                    *
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
 * tfr.is_complex      | must be set to FALSE (a RIDT tfr is real-valued)     *
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
 * normT               | normalization factor for the time window             *
 *                     |                                                      *
 * half_WindowF_Length | half-length of the frequency smoothing window        *
 * normF               | normalization factor for the frequency window        *
 *                     |                                                      *
 * RIDTKernel          | variable to compute the RIDT Kernel                  *
 * normK               | normalization factor for the Kernel                  *
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

int ctftbx_ridt (type_signal Signal, double *WindowT, int WindowT_Length,
                 double *WindowF, int WindowF_Length, type_TFR tfr)
{
  int            Mfft, Nfft, column, row, time;
  int            half_WindowT_Length, half_WindowF_Length;
  int            taumin, taumax, tau;
  int            mumin, mumax, mu;
  Complex       *lacf=0;		/* local autocorrelation function */
  double         normT, normF, normK;
  Complex        R1, R2;
  double        *RIDTKernel=0;

  /*--------------------------------------------------------------------*/
  /*                      Test the input variables                      */
  /*--------------------------------------------------------------------*/


  if (tfr.is_complex == TRUE)
  {
    if (ctftbx_debug)
      printf ("ridt.c : The tfr matrix must be real valued\n");
    return 1;
  }

  if (tfr.N_freq <= 0)
  {
    if (ctftbx_debug)
      printf ("ridt.c : The field tfr.N_freq is not correctly set\n");
    return 1;
  }

  if (tfr.N_time <= 0)
  {
    if (ctftbx_debug)
      printf ("ridt.c : The field tfr.N_time is not correctly set\n");
    return 1;
  }

  if (ISODD(WindowT_Length) == 0)
  {
    if (ctftbx_debug)
      printf ("ridt.c : The time-window Length must be an ODD number\n");
    return 1;
  }

  if (ISODD(WindowF_Length) == 0)
  {
    if (ctftbx_debug)
      printf ("ridt.c : The frequency-window Length must be an ODD number\n");
    return 1;
  }

  /*--------------------------------------------------------------------*/
  /*                    Determines some internal constants              */
  /*--------------------------------------------------------------------*/

  half_WindowT_Length = (WindowT_Length - 1) / 2;
  normT=WindowT[half_WindowT_Length];

  /* normalization of the time smoothing window */
  for(row = 0 ; row < WindowT_Length ; row++)
  {
    WindowT[row]=WindowT[row] / normT;
  }

  /* normalization of the frequency smoothing window */
  half_WindowF_Length = (WindowF_Length - 1) / 2;
  normF = WindowF[half_WindowF_Length];
  for(row = 0 ; row < WindowF_Length ; row++)
  {
    WindowF[row] = WindowF[row] / normF;
  }



  /*--------------------------------------------------------------------*/
  /*      Memory allocation and computation of  the kernel              */
  /*--------------------------------------------------------------------*/

  RIDTKernel = (double *) GC_malloc ( WindowT_Length * sizeof(double) );
  for(row = 0 ; row < WindowT_Length ; row++)
  {
    RIDTKernel[row]=0.0;
  }

  /*--------------------------------------------------------------------*/
  /*           creation of the vector of frequency bins  (output)       */
  /*--------------------------------------------------------------------*/
  Nfft = po2 (tfr.N_freq);

  for (row = 0; row < tfr.N_freq; row++)
  {
    tfr.freq_bins[row] = (double) (0.5 * row) / tfr.N_freq;
  }

  /*--------------------------------------------------------------------*/
  /*                memory allocation for the windowed signal           */
  /*--------------------------------------------------------------------*/
  lacf = (Complex *) GC_malloc (tfr.N_freq * sizeof (Complex));

  /* initialization of the intermediary vectors */
  for (row = 0; row < tfr.N_freq ; row++)
  {
    lacf[row] = 0.0;
  }

  /*--------------------------------------------------------------------*/
  /*      computation of the fft for the current windowed signal        */
  /*--------------------------------------------------------------------*/
  for (column = 0; column < tfr.N_time; column++)
  {

    /* time instants of interest to compute the tfr */
    time = ((int) tfr.time_instants[column]) - 1;

      /* maximum value of the delay in order to take the edges
    into account */
    taumax = MIN((time+half_WindowT_Length),(Signal.length-time-1+half_WindowT_Length));
    taumax = MIN(taumax,(tfr.N_freq / 2 - 1));
    taumax = MIN(taumax, half_WindowF_Length);

    /* initialization of the first local autocorrelation function */
    if (Signal.is_complex == TRUE)
    {
      lacf[0] = SQR(Signal.signal[time]) * WindowT[half_WindowT_Length];
    }
    else /* the signal is real-valued */
    {
      lacf[0] = RE(Signal.signal[time]) * RE(Signal.signal[time])
          * WindowT[half_WindowT_Length];
    }

    /* The signal is windowed around the current time */
    for (tau = 1; tau <= taumax; tau++)
    {
      R1=0.0;
      R2=0.0;

      /* bound of mu in order to take into account the edges */
      mumin=MIN(tau,half_WindowT_Length);
      mumin=MIN(mumin,(Signal.length-time-1-tau));

      mumax=MIN(tau,half_WindowT_Length);
      mumax=MIN(mumax,time-tau);

      normK = 0.0;
      for(mu = -mumin ; mu <= mumax ; mu++)
      {
        RIDTKernel[half_WindowT_Length + mu] =
            (1.0 - ABS((double)mu/(double)tau))
            *WindowT[half_WindowT_Length + mu];

        /* computation of the current kernel norm */
        normK = normK +  RIDTKernel[half_WindowT_Length + mu];
      }

      /* when the kernel is nearly null, no normalization */
      if (normK<EPS)
        normK=1.0;

      for(mu = -mumin ; mu <= mumax ; mu++)
      {
        /* case of complex valued signal */
        if (Signal.is_complex == TRUE)
        {
          R1 = R1 + Signal.signal[time+tau-mu] * conj(Signal.signal[time-tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;

          R2 = R2 + Signal.signal[time-tau-mu] * conj(Signal.signal[time+tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;
        }
        else
        {
          R1 = R1 + RE(Signal.signal[time+tau-mu]) * RE(Signal.signal[time-tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;

          R2 = R2 + RE(Signal.signal[time+tau-mu]) * RE(Signal.signal[time-tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;
        }

      }
      lacf[tau] = R1 * WindowF[half_WindowF_Length+tau];
      lacf[tfr.N_freq-tau] = R2 * WindowF[half_WindowF_Length-tau];
    }

    tau=floor(tfr.N_freq/2);
    if ((time<=Signal.length-tau-1)&(time>=tau)&(tau<=half_WindowF_Length))
    {
      R1=0.0;
      R2=0.0;

      /* bound of mu in order to take into account the edges */
      mumin=MIN(tau,half_WindowT_Length);
      mumin=MIN(mumin,(Signal.length-time-1-tau));

      mumax=MIN(tau,half_WindowT_Length);
      mumax=MIN(mumax,time-tau);


      normK = 0.0;
      for(mu = -mumin ; mu <= mumax ; mu++)
      {
        RIDTKernel[half_WindowT_Length + mu] =
            (1.0 - ABS((double)mu/(double)tau))
            *WindowT[half_WindowT_Length + mu];

        /* computation of the current kernel norm */
        normK = normK +  RIDTKernel[half_WindowT_Length + mu];
      }

      /* when the kernel is nearly null, no normalization */
      if (normK<EPS)
        normK=1.0;

      for(mu = -mumin ; mu <= mumax ; mu++)
      {
        /* case of complex valued signal */
        if (Signal.is_complex == TRUE)
        {
          R1 = R1 + Signal.signal[time+tau-mu] * conj(Signal.signal[time-tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;

          R2 = R2 + Signal.signal[time-tau-mu] * conj(Signal.signal[time+tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;
        }
        else
        {
          R1 = R1 + RE(Signal.signal[time+tau-mu]) * RE(Signal.signal[time-tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;

          R2 = R2 + RE(Signal.signal[time+tau-mu]) * RE(Signal.signal[time-tau-mu])
              * RIDTKernel[half_WindowT_Length+mu] / normK;
        }

      }
      lacf[tau] = 0.5*(R1 * WindowF[half_WindowF_Length+tau]
          +R2 * WindowF[half_WindowF_Length-tau]);
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

  /*--------------------------------------------------------------------*/
  /*                free the memory used in this program                */
  /*--------------------------------------------------------------------*/
  GC_free (lacf);
  GC_free (RIDTKernel);

  return 0;
}
