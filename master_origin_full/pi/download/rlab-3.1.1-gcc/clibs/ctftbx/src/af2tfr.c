/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 2008  Marijan Kostrun

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   See the file ./COPYING
   ********************************************************************** */

/* EXISTS AN INTERFACE PROGRAM TO MATLAB : TFRMEX.C		      *
 *====================================================================*
 * Name of the function :  af2tfr (void)        		      *
 * Author               :  Manuel DAVY - IRCYN                        *
 * Date of creation     :  20 - 01 - 1999                             *
 *--------------------------------------------------------------------*
 * THE ALGORITHM             				              *
 * Given the narrow band ambiguity function (AF) of a signal, and a   *
 * kernel, computes the TFR obtained by :                             *
 * 1 - masking the AF with the kernel                                 *
 * 2 - fft and ifft along the lines and columns		              *
 *                                                                    *
 *                                                                    *
 * The resulting Time-frequency representation (TFR) has as many      *
 * lines as the AF matrix has columns and as many oolumns as the AF   *
 * has lines                                                          *
 *                                                                    *
 *                                                                    *
 * The correct algorithm is summed up in this matlab command line :   *
 * TFR = real(fft(flipud(ifft(fftshift(AF.*kernel))')))               *
 *                                                                    *
 *====================================================================*
 * INPUT VARIABLES   					              *
 * Name        | type  |              role               	      *
 * AF.N_doppler| int   | Nombre de lignes du noyau i.e.               *
 *             |       | nombre de doppler bins                       *
 * N_Delay     | int   | Nombre de colonnes du noyau ie               *
 *             |       | nombre de delay bins                         *
 * AF_real     |double*| matrix (stored as a vector) containing the   *
 *             |       | real part of the signal's Ambiguity Function *
 * AF_imag     |double*| matrix (stored as a vector) containing the   *
 *             |       | imaginary part of the signal's Ambiguity     *
 *             |       | Function                                     *
 * kernel      |double*| matrix (stored as a vector) containing the   *
 *             |       | TFR kernel                                   *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES    						      *
 * Name        | type  |              role                 	      *
 * tfr         |double*| matrix (stored as a vector) containing the   *
 *             |       | TFR after computations			      *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES 						      *
 * Name        | type  |              role                	      *
 * delay       | int   |           loop variables used to             *
 * doppler     | int   |       specify an element in a matrix         *
 * time        | int   |       (the AF, the kernel. or any            *
 * freq        | int   |        intermediary matrix )                 *
 * N_time      | int   | Number of time bins (= N_Doppler)            *
 * po2_doppler | int   | result of po2(N_Doppler)                     *
 * po2_delay   | int   | result of po2(N_delay)                       *
 * imag        |double*| imaginary part of the TFR and the            *
 *             |       | intermediary matrices between (doppler,delay)*
 *             |       | and (time,frequency). The same size as 'tfr' *
 * inter       |double*| an intermediary matrix used to transpose     *
 *             |       | and flip up-down 'tfr' and 'imag'            *
 *====================================================================*
 * Name   | int idx(int line, int row, int nb_row)                    *
 * Action | computes the vector index for an element in a matrix given*
 *        | the line and column indices and the number of lines       *
 * Place  | divers.c                                                  *
 *--------------------------------------------------------------------*
 * Name   | void fft(int n, int m, double *x, double *y)              *
 * Action | Computes the fft                                          *
 * Place  | divers.c                                                  *
 *--------------------------------------------------------------------*
 * Name   | void ifft(int n, int m, double *x, double *y)             *
 * Action |Computes the inverse fft                                   *
 * Place  | divers.c                                                  *
 *--------------------------------------------------------------------*
 * Name   | int po2(int x)                                            *
 * Action | Computes the next power of two of x                       *
 * Place  | divers.c                                                  *
 *====================================================================*/

int
ctftbx_af2tfr (type_AF ambif, type_AF kernel, type_TFR tfr)

{
  int            delay, doppler, time, freq;
  int            index_1, index_2;
  int            po2_doppler, po2_delay;
  double         inter;


  /* at the moment programmed only for real-valued kernels */
  if (kernel.is_complex == TRUE)
  {
    if (ctftbx_debug)
      printf ("af2tfr: complex-valued kernels not implemented yet \n");
    return 1;
  }

  /* checks the sizes of the AF and the kernel */
  if ((ambif.N_doppler != kernel.N_doppler) || (ambif.N_delay != kernel.N_delay))
  {
    if (ctftbx_debug)
      printf ("af2tfr.c : kernel size different from AF size \n");
    return 1;
  }


  /*================================================================
     computation of the masked ambiguity function  and fft shifting
     this operation consists in inverting the four quadrans of the
     masked AF :
          -         |         -
          |       _ | _        |
          |      |\ | /|       |
          |        \|/         |
          |---------\----------|
          |       / |\         |
          |     |_  | _|       |
          |         |          |
          -                   -
   ============================================================== */

  for (delay = 0; delay < (ambif.N_delay - 1.0) / 2.0; delay++)
  {
    for (doppler = 0; doppler < (ambif.N_doppler - 1.0) / 2.0; doppler++)
    {
      /* First quadran */
      index_1 = idx (doppler, delay, ambif.N_doppler);
      index_2 = idx (doppler + ambif.N_doppler / 2, delay +
          ambif.N_delay / 2, ambif.N_doppler);

      tfr.tfr[index_1] = ambif.af[index_2] * RE(kernel.af[index_2]);

      /* Second quadran */
      index_1 = idx (doppler + ambif.N_doppler / 2, delay, ambif.N_doppler);
      index_2 = idx (doppler, delay + ambif.N_delay / 2, ambif.N_doppler);

      tfr.tfr[index_1] = ambif.af[index_2] * RE(kernel.af[index_2]);

      /* Third quadran */
      index_1 = idx (doppler, delay + ambif.N_delay / 2, ambif.N_doppler);
      index_2 = idx (doppler + ambif.N_doppler / 2, delay, ambif.N_doppler);

      tfr.tfr[index_1] = ambif.af[index_2] * RE(kernel.af[index_2]);

      /* Fouth quadran */
      index_1 = idx (doppler + ambif.N_doppler / 2, delay +
          ambif.N_delay / 2, ambif.N_doppler);
      index_2 = idx (doppler, delay, ambif.N_doppler);

      tfr.tfr[index_1] = ambif.af[index_2] * RE(kernel.af[index_2]);
    }
  }

  /*------------------------------------------------------*/
  /*  for the fft shifted matrix : transformation of the  */
  /*   doppler to the time : inverse fft for the columns  */
  /*------------------------------------------------------*/
  po2_doppler = po2 (ambif.N_doppler);
  for (delay = 0; delay < ambif.N_delay; delay++)
  {
    index_1 = idx (0, delay, ambif.N_doppler);
    ifft (ambif.N_doppler, po2_doppler, (Complex *) &tfr.tfr[index_1]);
  }


  /*------------------------------------------------------*/
  /* in order to transpose and conjugate the result,      */
  /* recopy of the  matrix in the inter matrix            */
  /*------------------------------------------------------*/

  transpose (tfr.N_time, ambif.N_delay, tfr.tfr);

  /* flips up - down and conjugate */
  for (delay = 0; delay <= ambif.N_delay / 2.0; delay++)  /* the first lines */
  {
    for (time = 0; time < tfr.N_time; time++)   /* all the columns */
    {
      index_1 = idx (delay, time, ambif.N_delay);
      index_2 = idx (ambif.N_delay - 1 - delay, time, ambif.N_delay);

      inter = RE(tfr.tfr[index_1]);
      RE(tfr.tfr[index_1]) = RE(tfr.tfr[index_2]);
      RE(tfr.tfr[index_2]) = inter;

      inter = -IM(tfr.tfr[index_1]);
      IM(tfr.tfr[index_1]) = -IM(tfr.tfr[index_2]);
      IM(tfr.tfr[index_2]) = inter;

    }
  }

  /*------------------------------------------------------*/
  /* transformation of the delay to the frequency :       */
  /*   fft over the columns                               */
  /*------------------------------------------------------*/
  po2_delay = po2 (ambif.N_delay);
  for (time = 0; time < tfr.N_time; time++)
  {
    /* computes the index of the first element of the  column no 'time' */
    index_1 = idx (0, time, ambif.N_delay);
    fft (ambif.N_delay, po2_delay, (Complex *) &tfr.tfr[index_1]);
  }
  tfr.is_complex = 0;

  return 0;
}
