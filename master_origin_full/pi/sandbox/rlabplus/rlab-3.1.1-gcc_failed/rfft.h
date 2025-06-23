/* rfft.h */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle

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

#ifndef RLAB_RFFT_H
#define RLAB_RFFT_H

#include "rlab.h"
#include "bltin.h"

#include <stdio.h>

extern void class_fft_init (void);

extern Ent *Real (int nargs, Datum args[]);
extern Ent *Imag (int nargs, Datum args[]);
extern Ent *Conj (int nargs, Datum args[]);
extern Ent *FFT  (int nargs, Datum args[]);
extern Ent *IFFT (int nargs, Datum args[]);
extern Ent *Filter (int nargs, Datum args[]);

//
// ctftbx library
//
extern Ent *ent_ctftbx_window   (int nargs, Datum args[]);
extern Ent *ent_ctftbx_stft     (int nargs, Datum args[]);
extern Ent *ent_ctftbx_fftshift (int nargs, Datum args[]);
// tfd:
extern Ent *ent_ctftbx_tfd_butter (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_bj (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_cw (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_grd (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_mh (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_page (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_pmh (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_ppage (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_pwv (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_ri (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_spwv (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_wv (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfd_zam (int nargs, Datum args[]);
// tfrid
extern Ent *ent_ctftbx_tfrid_bessel (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfrid_binom (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfrid_tri (int nargs, Datum args[]);
extern Ent *ent_ctftbx_tfrid_hann (int nargs, Datum args[]);
// af, kernel
extern Ent *ent_ctftbx_af (int nargs, Datum args[]);
extern Ent *ent_ctftbx_kernel (int nargs, Datum args[]);
extern Ent *ent_ctftbx_af2tfr (int nargs, Datum args[]);

#endif /* RLAB_RFFT_H */
