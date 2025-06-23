/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRRI.C                           *
 *============================================================================*
 * Name of the function : ri.c (void)                                         *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Rihaczek     *
 * Time-Frequency Distribution (RI) :                                         *
 *                                                                            *
 *                                       -j2pi f t                            *
 *                RI(t,f) =  x(t) X*(f) e                                     *
 *                                                                            *
 *                                                                            *
 * This function is real valued. Its computation requires a real or complex   *
 * signal, a vector containing time instants and the number of frequency bins *
 *                                                                            *
 *============================================================================*
 * INPUT VARIABLES                                                            *
 * Name                |              role                                    *
 * Signal              | The signal to analyze. No field modified             *
 *                     |                                                      *
 * tfr                 | Matrix containing the resulting TFR (complex)        *
 * tfr.time_instants   | positions of the smoothing window                    *
 * tfr.N_time          | length of '.time_instants' = number of cols.         *
 *                     | in the tfr matrix                                    *
 * tfr.N_freq          | number of frequency bins = number of rows in the tfr *
 *                     | matrix                                               *
 * tfr.is_complex      | must be set to TRUE (a RI tfr is complex-valued)     *
 *                     |                                                      *
 *----------------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                           *
 * Name                |                role                                  *
 * tfr.real_part       | the output tfr matrix  (real_part)                   *
 * tfr_imag_part       | the output tfr matrix  (imag_part)                   *
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
 * tau                 | time-lag variable                                    *
 * taumin              | local time-lag variable bounds. Used to take into    *
 * taumax              | accound the beginning and the end of the             *
 *                     | signal, where the window is cut                      *
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

int ctftbx_ri (type_signal Signal, type_TFR tfr)
{
  int            Nfft, column, row, time;
  int            taumin, taumax, tau;
  Complex       *lacf=0; /* local autocorrelation function */

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/


  if (tfr.is_complex == FALSE)
  {
    if (ctftbx_debug)
      printf ("ri.c : The tfr matrix must be complex valued\n");
    return 1;
  }

  if (tfr.N_freq <= 0)
  {
    if (ctftbx_debug)
      printf ("ri.c : The field tfr.N_freq is not correctly set\n");
    return 1;
  }

  if (tfr.N_time <= 0)
  {
    if (ctftbx_debug)
      printf ("ri.c : The field tfr.N_time is not correctly set\n");
    return 1;
  }


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

    /* taumax and taumin enable the computation near the edges */
    taumin = MIN( (tfr.N_freq - time), (Signal.length - time - 1) );
    taumax = time;

    /* The signal is windowed around the current time */
    for (tau = -taumin; tau <= taumax; tau++)
    {
      row = irem( (tfr.N_freq+tau), tfr.N_freq ) ;

      if (Signal.is_complex)
        /* when the signal is complex valued */
      {
        RE(lacf[row]) =  RE(Signal.signal[time]) * RE(Signal.signal[time - tau])
            + IM(Signal.signal[time]) * IM(Signal.signal[time - tau]);

        IM(lacf[row]) =  IM(Signal.signal[time]) * RE(Signal.signal[time - tau])
            - RE(Signal.signal[time]) * IM(Signal.signal[time - tau]);
      }
      else
        /* when the signal is real valued */
      {
        lacf[row] =  RE(Signal.signal[time]) * RE(Signal.signal[time - tau]);
      }
    }


    /* fft of the local autocorrelation function lacf */
    fft (tfr.N_freq, Nfft, lacf);

    /* the fft is put in the tfr matrix  */
    for (row = 0; row < tfr.N_freq; row++)
    {
      tfr.tfr[idx (row,column,tfr.N_freq)] = lacf[row];

      lacf[row] = 0.0;
    }
  }
  /*--------------------------------------------------------------------*/
  /*                free the memory used in this program                */
  /*--------------------------------------------------------------------*/
  GC_free (lacf);

  return 0;
}
