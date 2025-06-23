

/* EXISTS AN INTERFACE PROGRAM TO MATLAB : DISTMEX.C */

/*=====================================================================
This file contains :
   - Renyi : computes the renyi info of a TFR
   - Jensen_inter_index : intermediary distance index
          necessary to compute the Jensen distance
   - distance : computes a distance between TFRs
=====================================================================*/





/*====================================================================*
 * Name of the function : Renyi (double)                              *
 * Author             : Manuel DAVY                                   *
 * Date of creation   : 20 - 01 - 1999                                *
 *------------------------------------------------------------------- *
 * THE ALGORITHM                                                      *
 * Computes the Renyi information, of exponenet alpha of a TFR        *
 *                                                                    *
 *====================================================================*
 * INPUT VARIABLES                                                    *
 * Name           |                   role                            *
 * TFR            | The TFR for which the Renyi information is        *
 *                | computed                                          *
 * .N_time        | Number of columns in the TFR matrix               *
 * .N_freq        | Number of rows in the TFR matrix                  *
 * alpha          | Exponent of the Renyi information formula         *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name           |                   role                            *
 * resu           | Value of the REnyi info of the TFR, computed with *
 *                | exponent alpha                                    *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                 *
 * Name           |                   role                            *
 * time, freq     | Row and column index in the matrix TFR            *
 *====================================================================*/


/* Computes the Renyi information of a TFR, for the order alpha */

static double ctftbx_renyi (type_TFR tfr, double alpha)
{
  double         resu, inter;		/* output result */
  int            time, freq;


  if (alpha != 1.0)
  {
    resu = 0;
    for (time = 0; time < tfr.N_time; time++)
    {
      for (freq = 0; freq < tfr.N_freq; freq++)
      {
        resu = resu + powof (RE(tfr.tfr[idx (freq, time, tfr.N_time)]), alpha);
      }
    }
    resu = (1.0 / (1.0 - alpha)) * log (resu) / log (2.0);  /* log in base 2 */
  }
  else /* Shanon entropy : the TFR must be >0, else take the absol. value */
  {
    resu = 0;
    for (time = 0; time < tfr.N_time; time++)
    {
      for (freq = 0; freq < tfr.N_freq; freq++)
      {
        inter =  ABS(RE(tfr.tfr[idx (freq, time, tfr.N_time)]));
        resu = resu + inter * log(inter)/log(2.0);
      }
    }
    resu = - resu;
  }
  return resu;
}




/*====================================================================*
 * Name of the function : Jensen_inter_index (double)                 *
 * Author               : Manuel DAVY                                 *
 * Date of creation     : 20 - 01 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM                                                      *
 * Computes an intermediairy index necessary for the Jensen distance  *
 * this index is the difference of Renyi info (of exponent alpha of   *
 * the TFR) of the geometrical average TFR and arithmetical average   *
 * TFR                                                                *
 *         ______________               TFR_1 + TFR_2                 *
 * Renyi(V TFR_1 * TFR_2,alpha)- Renyi(--------------- , alpha)       *
 *                                            2                       *
 *====================================================================*
 * INPUT VARIABLES                                                    *
 * Name           |                    role                           *
 * TFR_1, TFR_2   | TFR for which the intermediairy index is computed *
 * . N_time       | Number of columns (time instants) in  TFR_*       *
 * . N_freq       | Number of lines (frequency bins) in  TFR_*        *
 * alpha          | Exponent of the Renyi information formula         *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name           |                   role                            *
 * resu           | Result of the computation                         *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                 *
 * Name           |                  role                             *
 * TFR_inter      | Intermediary TFR containing the geometrical       *
 *                | Average of TFR_1 and TFR_2                        *
 * indice         | Index for the matrices TFR_1 and TFR_2            *
 *                | (stored as vectors)                               *
 *====================================================================*
 * SUBROUTINES USED HERE                                              *
 *--------------------------------------------------------------------*
 * Name   |                                                           *
 * Action |                                                           *
 * Place  |                                                           *
 *====================================================================*/

static double
ctftbx_Jensen_inter_index (type_TFR tfr_1, type_TFR tfr_2, double alpha)
{
  double         resu;		/* output result */
  type_TFR       TFR_inter;
  int            indice;

/* finds memory to store the geometrical average matrix */
  TFR_inter.N_time = tfr_1.N_time;
  TFR_inter.N_freq = tfr_1.N_freq;
  TFR_inter.is_complex = FALSE;

  TFR_inter.tfr = (Complex *) GC_malloc (TFR_inter.N_time *
      TFR_inter.N_freq * sizeof (Complex));


  /* computes the geometrical average of TFR_1 and TFR_2 */
  for (indice = 0; indice < tfr_1.N_time * tfr_1.N_freq; indice++)
    TFR_inter.tfr[indice] = sqrt (RE(tfr_1.tfr[indice])
        * RE(tfr_2.tfr[indice]));

  resu = ctftbx_renyi (TFR_inter, alpha) -
      (ctftbx_renyi (tfr_1, alpha) + ctftbx_renyi (tfr_2, alpha)) / 2;

  GC_free (TFR_inter.tfr);

  return resu;
}


/* EXISTS AN INTERFACE PROGRAM TO MATLAB : ---MEX.c                   *
 *====================================================================*
 * Name of the function : distance                                    *
 * Author               : Manuel DAVY                                 *
 * Date of creation     : 20 - 01 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM                                                      *
 * Given a distance measure, computes the distance between two        *
 * Time-Feequency Representations (TFRs). Some distances require a    *
 * normalization :                                                    *
 *                                                                    *
 *                     | TFR (t,f)|                                   *
 * TFR_norm(t,f) = ------------------------                           *
 *                 / /                                                *
 *                 | |  | TFR(t,f)| dtdf                              *
 *                 / /                                                *
 *                                                                    *
 *====================================================================*
 * INPUT VARIABLES                                                    *
 * Name           |                   role                            *
 *  first_TFR     | TFR  for which the distance is computed           *
 *   .N_time      |  Number of columns (resp. rows ) in the TFRs      *
 *   .N_freq      |                                                   *
 *  second_TFR    | TFR  for which the distance is computed           *
 *   .N_time      |  Number of columns (resp. rows ) in the TFRs      *
 *   .N_freq      |                                                   *
 *  name          | identifier for the distance measure to compute    *
 *  coef          | coefficient required to compute some distance     *
 *                | measures                                          *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name           |                   role                            *
 * dist           | Result of the computation : the distance between  *
 *                | TFRs                                              *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                 *
 * Name           |                   role                            *
 * time,freq      | index for columns (resp. lines) in the TFR        *
 *                | matrices                                          *
 * N_time,N_freq  | Number of columns (resp rows) in the TFR matrices *
 * first_sum      | Result of the  double integration  of the         *
 * second_sum     | absolute values of each TFR. Used to normalize    *
 *                | the TFRs                                          *
 * distan         | Value of the distance that is passed as a result  *
 *                | in '*dist'                                        *
 * inter          | Intermediary value in the computation of the TFRs *
 * tfr1_local     | An element of the matrix first_TFR (resp.         *
 * tfr2_local     | second_tfr) in line freq and column time,         *
 *                | normalized as explained before                    *
 * index          | Used to access an element in a TFR matrix (stored *
 *                | as a vector)                                      *
 *              Special for the Jensen divergence                     *
 * TFR_2_norm     | Second TFR normalized, to pass to                 *
 *                | Jensen_inter_index                                *
 * TFR_inter      | average normalized TFRs (first and second) to     *
 *                | to Jensen_inter_index                             *
 *====================================================================*
 * SUBROUTINES USED HERE                                              *
 *--------------------------------------------------------------------*
 * Name   |                                                           *
 * Action |                                                           *
 * Place  |                                                           *
 *====================================================================*/


static int
ctftbx_enum_distance (type_TFR first_TFR, type_TFR second_TFR,
                      enum ctftbx_dist name, double coef, double *dist)
{
  int            time, freq, N_time, N_freq;
  double         first_sum, second_sum, distan;
  double         inter, tfr1_local, tfr2_local;
  int            index;

  /* specific for the Jensen distance */
  type_TFR       TFR_2_norm, TFR_inter;

  if ((first_TFR.N_time != second_TFR.N_time) ||
       (first_TFR.N_freq != second_TFR.N_freq))
  {
    if (ctftbx_debug)
      printf ("distance.c : The two TFR do not have the same size\n");
    return 1;

  }

  N_time = first_TFR.N_time;
  N_freq = first_TFR.N_freq;

/*====================================================================
  Normalization of the TFRs (real-valued)
  ====================================================================*/

/* computation of the double integral
  of the absolute values of each TFR */
  if ((name != CTFTBX_LQ) && (name != CTFTBX_CORRELATION) && (name != CTFTBX_QUADRATIC))
  {
    first_sum = 0;
    second_sum = 0;
    for (time = 0; time < N_time; time++)
    {
      for (freq = 0; freq < N_freq; freq++)
      {
        index = idx (freq, time, N_freq);
        first_sum = first_sum + ABS (RE(first_TFR.tfr[index]));
        second_sum = second_sum + ABS (RE(second_TFR.tfr[index]));
      }
    }
  }

/*====================================================================
  computation of the distance
  ====================================================================*/
  switch (name)
  {

// ***************** LQ distance **********************************
//
//        -                                    -  1/coef
//        | / /                        coef    |
//     d= | | |  | TFR1(t,f)-TFR2(t,f)|    dtdf|
//        | / /                                |
//        -                                    -

    case CTFTBX_LQ:
      distan = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          tfr1_local = RE(first_TFR.tfr[index]);
          tfr2_local = RE(second_TFR.tfr[index]);
          inter = tfr1_local - tfr2_local;
          distan = distan + powof (ABS (inter), coef);
        }
      }
      distan = powof (distan, 1.0 / coef);
      break;

// ******************* Quadratic distance *************************
//
//     / /                        2
// d=  | |  | TFR1(t,f)-TFR2(t,f)| dtdf
//     / /

    case CTFTBX_QUADRATIC:
      distan = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          tfr1_local = RE(first_TFR.tfr[index]);
          tfr2_local = RE(second_TFR.tfr[index]);
          inter = tfr1_local - tfr2_local;
          distan = distan + inter * inter;
        }
      }

      break;

// ************************ Correlation distance ********************
//
//         / /
//  d= 1 - | | TFR1(t,f)*TFR2(t,f) dtdf
//         / /

    case CTFTBX_CORRELATION:
      distan = 0;
      first_sum = 0;
      second_sum = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index =  idx (time, freq, N_time);
          tfr1_local = RE(first_TFR.tfr[index]);
          tfr2_local = RE(second_TFR.tfr[index]);
          first_sum  = first_sum + tfr1_local * tfr1_local;
          second_sum = second_sum + tfr2_local * tfr2_local;

          inter = tfr1_local * tfr2_local;
          distan = distan + inter;
        }
      }
      distan = 1 - distan / (first_sum + second_sum);
      break;


// ************************* Kolmogorov distance *********************
//
//      / /
//  d=  | |  | TFR1_norm(t,f)-TFR2_norm(t,f)| dtdf
//      / /

    case CTFTBX_KOLMOGOROV:
      distan = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          tfr1_local = ABS (RE(first_TFR.tfr[index])) / first_sum;
          tfr2_local = ABS (RE(second_TFR.tfr[index])) / second_sum;
          inter = tfr1_local - tfr2_local;
          distan = distan + ABS (inter);
        }
      }
      break;

// *************************** Kullback distance *********************
//
//      / /                                    TFR1_norm(t,f)
//  d=  | |(TFR1_norm(t,f)-TFR2_norm(t,f))*log--------------- dtdf
//      / /                                    TFR2_norm(t,f)

    case CTFTBX_KULLBACK:
      distan = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          tfr1_local = ABS (RE(first_TFR.tfr[index])) / first_sum;
          tfr2_local = ABS (RE(second_TFR.tfr[index])) / second_sum;
          if ((tfr1_local != 0) && (tfr2_local != 0))
          {
            inter = (tfr1_local - tfr2_local) * log (tfr1_local / tfr2_local);
          }
          else
            /* the distance is not defined , then the null points do not count */
          {
            inter = 0;
          }
          distan = distan + ABS (inter);
        }
      }
      break;

// *********************** Chernoff distance *************************
//           -
//           | / /               coef                1/coef   |
//  d= -log  | | | TFR1_norm(t,f)     *TFR2_norm(t,f)     dtdf|
//           | / /                                            |
//           -                                                -

      case CTFTBX_CHERNOFF:    /*Chernoff distance */
        distan = 0;
        for (time = 0; time < N_time; time++)
        {
          for (freq = 0; freq < N_freq; freq++)
          {
            index = idx (freq, time, N_freq);
            tfr1_local = ABS (RE(first_TFR.tfr[index])) / first_sum;
            tfr2_local = ABS (RE(second_TFR.tfr[index])) / second_sum;
            inter = powof (tfr1_local, coef) * powof (tfr2_local, 1.0 - coef);
            distan = distan + inter;
          }
        }
        distan = -log (distan);
        break;

// ********************** Generalized Matusita distance **************
//
//       -                                                       -  1/coef
//      | / /               1/coef               1/coef  coef     |
//  d=  | | | |TFR1_norm(t,f)      -TFR2_norm(t,f)      |     dtdf|
//      | / /                                                     |
//       -                                                       -

    case CTFTBX_MATUSITA:
      distan = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          tfr1_local = ABS (RE(first_TFR.tfr[index])) / first_sum;
          tfr2_local = ABS (RE(second_TFR.tfr[index])) / second_sum;
          inter = powof (tfr1_local, 1.0 / coef) - powof (tfr2_local,
          1.0 / coef);
          distan = distan + powof (ABS (inter), coef);
        }
      }
      distan = powof (distan, 1.0 / coef);
      break;

      /********************** Normalized Lq distance **********************/
/*---------------------------------------------------------------
      -                                             -    1/coef
      | / /                                  coef    |
      d= | | |  | TFR1_norm(t,f)-TFR2_norm(t,f)|    dtdf|
      | / /                                          |
      -                                             -
      ---------------------------------------------------------------*/
    case CTFTBX_NLQ:
      distan = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          tfr1_local = ABS (RE(first_TFR.tfr[index])) / first_sum;
          tfr2_local = ABS (RE(second_TFR.tfr[index])) / second_sum;
          inter = tfr1_local - tfr2_local;
          distan = distan + powof (ABS (inter), coef);
        }
      }
      distan = powof (distan, 1.0 / coef);
      break;

// *************** q Normalized Log Spectral deviation **************
//
//      -                                                      -  1/coef
//     | / /                                            coef    |
//  d= | | |  | log(TFR1_norm(t,f))-log(TFR2_norm(t,f))|    dtdf|
//     | / /                                                    |
//      -                                                      -

    case CTFTBX_LSD:
      distan = 0;
      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          tfr1_local = ABS (RE(first_TFR.tfr[index])) / first_sum;
          tfr2_local = ABS (RE(second_TFR.tfr[index])) / second_sum;
          inter = log (tfr1_local) - log (tfr2_local);
          distan = distan + powof (ABS (inter), coef);
        }
      }
      distan = powof (distan, 1.0 / coef);
      break;


// *************************** Jensen Divergence ********************

    case CTFTBX_JENSEN:
      /* uses Renyi and Jensen_inter_index */
      TFR_inter.N_time = N_time;
      TFR_inter.N_freq = N_freq;
      TFR_2_norm.N_time = N_time;
      TFR_2_norm.N_freq = N_freq;
      TFR_inter.is_complex = FALSE;
      TFR_2_norm.is_complex = FALSE;

      TFR_inter.tfr = (Complex *) GC_malloc (TFR_inter.N_time *
          TFR_inter.N_freq * sizeof (Complex));
      TFR_2_norm.tfr = (Complex *) GC_malloc (TFR_2_norm.N_time *
          TFR_2_norm.N_freq * sizeof (Complex));


      for (time = 0; time < N_time; time++)
      {
        for (freq = 0; freq < N_freq; freq++)
        {
          index = idx (freq, time, N_freq);
          TFR_2_norm.tfr[index] = ABS (RE(second_TFR.tfr[index]))
              / second_sum;
          TFR_inter.tfr[index]  = (ABS (RE(first_TFR.tfr[index]))
              / first_sum
              + RE(TFR_2_norm.tfr[index])) / 2;
        }
      }
      distan = ctftbx_Jensen_inter_index (TFR_inter, TFR_2_norm, coef);

      GC_free (TFR_inter.tfr);
      GC_free (TFR_2_norm.tfr);

      break;

  }


/*====================================================================
  the final distance is stored in the output variable
  ====================================================================*/
  *dist = distan;

  return 1;
}

int
ctftbx_distance (type_TFR first_TFR, type_TFR second_TFR,
                 char *name_char, double coef, double *dist)
{
  int j;
  enum ctftbx_dist name = CTFTBX_UNKNOWN_DIST;

  if (!name_char)
    return 2;

  switch (name_char[0])
  {
    case 'L':
    case 'l':
      if (name_char[1]=='Q' || name_char[1]=='q')
        name = CTFTBX_LQ;
      else if (name_char[1]=='S' || name_char[1]=='s')
        name = CTFTBX_LSD;
      break;

    case 'Q':
    case 'q':
      name = CTFTBX_QUADRATIC;
      break;

    case 'C':
    case 'c':
      if (name_char[1]=='O' || name_char[1]=='o')
        name = CTFTBX_CORRELATION;
      else if (name_char[1]=='H' || name_char[1]=='h')
        name = CTFTBX_CHERNOFF;
      break;

    case 'K':
    case 'k':
      if (name_char[1]=='O' || name_char[1]=='o')
        name = CTFTBX_KOLMOGOROV;
      else if (name_char[1]=='U' || name_char[1]=='u')
        name = CTFTBX_KULLBACK;
      break;

    case 'M':
    case 'm':
      name = CTFTBX_MATUSITA;
      break;

    case 'N':
    case 'n':
      name = CTFTBX_NLQ;
      break;

    case 'J':
    case 'j':
      name = CTFTBX_JENSEN;
      break;
  }

  if (name == CTFTBX_UNKNOWN_DIST)
    return 2;

  j = ctftbx_enum_distance (first_TFR, second_TFR, name, coef, dist);

  return j;
}


