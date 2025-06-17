/* tftb.h: time-frequency toolbox */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
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

/* Signal structure */
typedef struct SIG
{
    int            length;          // Length of the signal in points
    double         sample_freq;     // Sample frequency of the signal
    double        *time_instants;   // instants of sample for the signal
    unsigned char  is_complex;      // TRUE if there exists an imag part
    Complex       *signal;          // real part of the signal

}
type_signal, kernel_type_signal;

typedef struct Time_freq_rep
{
    int            N_freq;        // number of freq bins in the TFR matrix
    int            N_time;        // number of time_bins in the TFR matrix
    double        *freq_bins;     // fqs for each line of the matrix
    double        *time_instants; // instant for each column of the TFR
    unsigned char  is_complex;    // TRUE if there exists an imag part
    Complex       *tfr;           // tfr
}
type_TFR;

typedef struct Ambi_func
{
    int            N_doppler;     // number of doppler bins in the AF
    int            N_delay;	      // number of delay bins in the AF matrix
    double        *doppler_bins;  // doppler bin for each line of the AF: x-coordinate
    double        *delay_bins;    // delay bin for each column of the AF: y-coordinate
    unsigned char  is_complex;    // TRUE if there exists an imag part
    Complex       *af;            // AF
}
type_AF, kernel_type_AF;


