/* libtftbx.h: time-frequency toolbox include file for rlab */

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

//
// time-frequency structures
//
#include "tftb.h"

//
// functions provided by the library
//
extern int double_fftshift  (double *vector_in, double *vector_out, int vector_length);
extern int Complex_fftshift (Complex *vector_in, Complex *vector_out, int vector_length);

extern int ctftbx_af2tfr(type_AF ambif, type_AF kernel, type_TFR tfr);
extern int ctftbx_af    (type_signal Signal, type_AF AF);
extern int ctftbx_bj    (type_signal Signal,  double *WindowT, int WindowT_Length,
                         double *WindowF, int WindowF_Length,  type_TFR tfr);
extern int ctftbx_bud   (type_signal Signal, double *WindowT, int WindowT_Length,
                         double *WindowF, int WindowF_Length, double sigma, type_TFR tfr);
extern int ctftbx_cw    (type_signal Signal, double *WindowT, int WindowT_Length,
                         double *WindowF, int WindowF_Length, double sigma, type_TFR tfr);
extern int ctftbx_distance (type_TFR first_TFR, type_TFR second_TFR,
                            enum ctftbx_dist name, double coef, double *dist);
extern int ctftbx_grd (type_signal Signal, double *WindowT, int WindowT_Length,
                       double *WindowF, int WindowF_Length, double rs, double MoverN,
                       type_TFR tfr);
extern int ctftbx_hough (type_TFR tfr, double nb_theta,double  nb_rho,
                         double* transfo_hough, double* rho_vect, double* theta_vect);
extern int ctftbx_kernel (enum ctftbx_kernel_shape kernel_type,
                          double *parameters, int nb_param, kernel_type_AF ker);
extern int ctftbx_mh (type_signal Signal, type_TFR tfr);
extern int ctftbx_mhs (type_signal Signal, double *WindowG, int WindowG_Length,
                       double *WindowH, int WindowH_Length, type_TFR tfr );
extern int ctftbx_mmce (type_signal Signal, double *Window,
                        int Window_Length, int Window_col, type_TFR tfr);
extern int ctftbx_page(type_signal Signal, type_TFR tfr);
extern int ctftbx_pmh (type_signal Signal, double *Window, int Window_Length, type_TFR tfr);
extern int ctftbx_ppage (type_signal Signal, double *Window, int Window_Length, type_TFR tfr);
extern int ctftbx_pwv (type_signal Signal, double *Window, int Window_Length, type_TFR tfr);
extern int ctftbx_ri (type_signal Signal, type_TFR tfr);
extern int ctftbx_ridb (type_signal Signal, double *WindowT, int WindowT_Length,
                        double *WindowF, int WindowF_Length, type_TFR tfr);
extern int ctftbx_ridbn (type_signal Signal, double *WindowT, int WindowT_Length,
                  double *WindowF, int WindowF_Length, type_TFR tfr);
extern int ctftbx_ridh (type_signal Signal, double *WindowT, int WindowT_Length,
                        double *WindowF, int WindowF_Length, type_TFR tfr);
extern int ctftbx_ridt (type_signal Signal, double *WindowT, int WindowT_Length,
                        double *WindowF, int WindowF_Length, type_TFR tfr);
extern int ctftbx_spwv (type_signal Signal, double *WindowT, int WindowT_Length,
                        double *WindowF, int WindowF_Length, type_TFR tfr);
extern int ctftbx_stft (type_signal Signal, double *Window, int Window_Length,
                        type_TFR tfr, double *norm_vector);
extern int ctftbx_wv (type_signal Signal, type_TFR tfr);
extern int ctftbx_zam (type_signal Signal, double *WindowT, int WindowT_Length,
                       double *WindowF, int WindowF_Length, type_TFR tfr);









