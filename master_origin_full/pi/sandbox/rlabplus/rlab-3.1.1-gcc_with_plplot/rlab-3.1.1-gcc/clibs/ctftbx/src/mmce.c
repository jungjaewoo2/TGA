/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRMMCE.C                         *
 *============================================================================*
 * Name of the function : mmce.c (void)                                       *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Minimu Mean  *
 * Cross-Entropy combination of spectrograms (MMCE) :                         *
 *                                                                            *
 *                         E                  __N                             *
 * MMCE(t,f) = ----------------------------   ||    |F(t,f;k)|^(2/N)          *
 *             ||  __N                    ||    k=1                           *
 *             ||  ||    |F(t,f;k)|^(2/N) ||                                  *
 *             ||    k=1                  ||1                                 *
 *                                                                            *
 * where || ||1 denotes the L1 norm, E the energy of the signal :             *
 *                                                                            *
 *         /               //                   ||         ||                 *
 *     E = | |x(t)|^2 dt = || MMCE(t,f) dt df = ||MMCE(t,f)||                 *
 *         /               //                   ||         ||1                *
 *                                                                            *
 * and F(t,f;k) the short-time Fourier transform of the signal, with the kieme*
 * analysis window h                                                          *
 *                                                                            *
 * This function is real valued. Its computation requires k frequency         *
 * smoothing windows, their displacement positions and the number of frequency*
 * bins to be computed with discrete variables.                               *
 *                                                                            *
 *============================================================================*
 * INPUT VARIABLES                                                            *
 * Name                |              role                                    *
 * Signal              | The signal to analyze. No field modified             *
 *                     |                                                      *
 * Window              | Matrix containing the points of the k frequency      *
 *                     | windows                                              *
 * Window_Length       | Number of points of the windows (ODD number !)       *
 * Window_col          | Number of analysis windows                           *
 *                     |                                                      *
 * tfr                 | Matrix containing the resulting TFR (real)           *
 * tfr.time_instants   | positions of the smoothing window                    *
 * tfr.N_time          | length of '.time_instants' = number of cols.         *
 *                     | in the tfr matrix                                    *
 * tfr.N_freq          | number of frequency bins = number of rows in the tfr *
 *                     | matrix                                               *
 * tfr.is_complex      | must be set to FALSE (a BJ tfr is real-valued)       *
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
 * half_Window_Length  | half-length of the frequency smoothing window        *
 * Nb_Window           | number of analysis windows h                         *
 * normh               | normalization factor for the frequency windows       *
 *                     |                                                      *
 * tau                 | time-lag variable                                    *
 * taumin              | local time-lag variable bounds. Used to take into    *
 * taumax              | accound the beginning and the end of the             *
 *                     | signal, where the window is cut                      *
 *                     |                                                      *
 * wind_sig_real       | real and imaginary parts of the windowed analysed    *
 * wind_sig_imag       | signal                                               *
 *                     |                                                      *
 * WIND_SIG_REAL       | matrix to store real and imaginary parts of the      *
 * WIND_SIG_IMAG       | signal windowed by the different analysis windows    *
 *                     |                                                      *
 *============================================================================*
 * SUBROUTINES USED HERE                                                      *
 *----------------------------------------------------------------------------*
 * Name   | int idx(int i_row, int j_col, int nb_row)                         *
 * Action | computes the vector index for an element in a matrix given the row*
 *        | and column indices (i,j) and the total number of row              *
 * Place  | divers.c                                                          *
 *----------------------------------------------------------------------------*
 * Name   | int irem( double x, double y)                                     *
 * Action | computes the remainder after Euclidean division of double         *
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
ctftbx_mmce (type_signal Signal, double *Window,
             int Window_Length, int Window_col, type_TFR tfr)
{
  int            Nfft, column, row, time;
  int            taumin, taumax, tau;
  int            half_Window_Length;
  Complex       *wind_sig=0; /* windowed signal */
  Complex       *WIND_SIG=0;
  double         normh, Nb_Window;
  int            i, j, hcol;
  double         WIN_ij, POW_Wind_sig, SUM_sig2, NORM_TFR, PROD_sig;
  char           STR[50];
 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/

  if (tfr.is_complex == TRUE)
  {
    if (ctftbx_debug)
      printf ("mmce.c : The tfr matrix must be real valued\n");
    return 1;
  }

  if (tfr.N_freq <= 0)
  {
    if (ctftbx_debug)
      printf ("mmce.c : The field tfr.N_freq is not correctly set\n");
    return 1;
  }

  if (tfr.N_time <= 0)
  {
    if (ctftbx_debug)
      printf ("mmce.c : The field tfr.N_time is not correctly set\n");
    return 1;
  }

 /*--------------------------------------------------------------------*/
 /*                   checks that the window length is odd             */
 /*--------------------------------------------------------------------*/

  if (ISODD(Window_Length) == 0)
  {
    if (ctftbx_debug)
      printf ("mmce.c : The window Length must be an ODD number\n");
    return 1;
  }

  if (Window_col < 2)
  {
    if (ctftbx_debug)
      printf ("mmce.c : The window must have at least 2 columns\n");
    return 1;
  }

  half_Window_Length = (Window_Length - 1) / 2;

  /* Normalization of all Windows */
  for(j=0;j<Window_col;j++)
  {
    normh=0;
    for(i=0;i<Window_Length;i++)
    {
      WIN_ij = Window[idx(i,j,Window_Length)];
      normh = normh + ABS( WIN_ij*WIN_ij );
    }
    normh = sqrt(normh)+EPS;

    for(i=0;i<Window_Length;i++)
    {
      Window[idx(i,j,Window_Length)] = Window[idx(i,j,Window_Length)]/normh;
    }
  }

  Nb_Window = Window_col;

  /*--------------------------------------------------------------------*/
  /*           creation of the vector of frequency bins  (output)       */
  /*--------------------------------------------------------------------*/
  Nfft = po2 (tfr.N_freq);

  for (row = 0; row < tfr.N_freq; row++)
  {
    tfr.freq_bins[row] = (double) row / tfr.N_freq;
  }

  /*--------------------------------------------------------------------*/
  /*                memory allocation for the windowed signal           */
  /*--------------------------------------------------------------------*/
  wind_sig = (Complex *) GC_malloc (tfr.N_freq * sizeof (Complex));

  /* initialization of the intermediary vectors */
  for (row = 0; row < tfr.N_freq; row++)
  {
    wind_sig[row] = 0.0;
  }


  WIND_SIG = (Complex *) GC_malloc (tfr.N_freq*Window_col * sizeof (Complex));

  /* initialization of the intermediary vectors */
  for (row = 0; row < (tfr.N_freq*Window_col); row++)
  {
    WIND_SIG[row] = 0.0;
  }
  /*--------------------------------------------------------------------*/
  /*      computation of the fft for the current windowed signal        */
  /*--------------------------------------------------------------------*/

  for (column = 0; column < tfr.N_time; column++)
  {

    /* time instants of interest to compute the tfr */
    time = ((int) tfr.time_instants[column]) - 1;


      /* the signal is multipied by the window between the instants
    time-taumin and time+taumax */
      /* when the window is wider than the number of desired frequencies (tfr.N_freq),
    the range is limited to tfr.N_freq */
    taumin = MIN (tfr.N_freq / 2, half_Window_Length);
    taumin = MIN (taumin, time);

    taumax = MIN ((tfr.N_freq / 2 - 1), half_Window_Length);
    taumax = MIN (taumax, (Signal.length - time - 1));

    /* The signal is windowed around the current time */
    for (tau = -taumin; tau <= taumax; tau++)
    {
      row = irem( (tfr.N_freq+tau),  tfr.N_freq ) ;

      for ( hcol=0 ; hcol<Window_col ; hcol++ )
      {
        WIND_SIG[idx(row,hcol,tfr.N_freq)] = Signal.signal[time + tau]
            * Window[idx((half_Window_Length + tau),hcol,Window_Length)];
      }
    }


    for ( hcol=0 ; hcol<Window_col ; hcol++ )
    {
      for (row = 0; row < tfr.N_freq; row++)
      {
        wind_sig[row] = WIND_SIG[idx(row,hcol,tfr.N_freq)];
      }

      /* fft of the windowed signal */
      fft (tfr.N_freq, Nfft, wind_sig);

      for (row = 0; row < tfr.N_freq; row++)
      {
        WIND_SIG[idx(row,hcol,tfr.N_freq)] =
            RE(wind_sig[row]) * RE(wind_sig[row]) + IM(wind_sig[row]) * IM(wind_sig[row]);
      }
    }


    for (row = 0; row < tfr.N_freq; row++)
    {
      PROD_sig = 1;
      for ( hcol=0 ; hcol<Window_col ; hcol++ )
      {
        PROD_sig = PROD_sig * RE(WIND_SIG[idx(row,hcol,tfr.N_freq)]);
        WIND_SIG[idx(row,hcol,tfr.N_freq)] = 0;
      }

      POW_Wind_sig = pow( PROD_sig,(1/Nb_Window));
      tfr.tfr[idx (row,column,tfr.N_freq)] = POW_Wind_sig;
      wind_sig[row] = 0.0;
    }



  }
  /*--------------------------------------------------------------------*/
  /*                free the memory used in this program                */
  /*--------------------------------------------------------------------*/
  GC_free (wind_sig);
  GC_free (WIND_SIG);

  SUM_sig2 = 0;
  NORM_TFR = 0;
  for(column=0;column<tfr.N_time;column++)
  {
    time = ((int) tfr.time_instants[column]) - 1;
    if (Signal.is_complex == TRUE)
    {
      SUM_sig2 = SUM_sig2 + RE(Signal.signal[time]) * RE(Signal.signal[time])
          + IM(Signal.signal[time]) * IM(Signal.signal[time]);
    }
    else
    {
      SUM_sig2 = SUM_sig2 + RE(Signal.signal[time]) * RE(Signal.signal[time]);
    }

    for(row=0; row<tfr.N_freq; row++)
    {
      NORM_TFR = NORM_TFR + RE(tfr.tfr[idx (row,column,tfr.N_freq)]);
    }
  }


  for(row=0; row<tfr.N_freq; row++)
  {
    for(column=0;column<tfr.N_time;column++)
    {
      tfr.tfr[idx (row,column,tfr.N_freq)] =  SUM_sig2/NORM_TFR
          * RE(tfr.tfr[idx (row,column,tfr.N_freq)]);
    }
  }

  return 0;
}
