// init.c
/* init.c: Perform initialization of the main symbol-table. i.e.
   install builtin functions, predefined variables etc. */

/* This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1993, 1994, 1995, 1996  Ian R. Searle
   This file is a part of rlabplus
   Copyright (C) 2005-2017 Marijan Kostrun

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
// names of some solvers, and their parameters
//
#include "rlab_solver_parameters_names.h"

#include "rlab.h"
#include "symbol.h"
#include "util.h"
#include "print.h"
#include "btreef1.h"

#include <stdio.h>

/* main.c */
extern Ent *Eval (int nargs, Ent * args);

//
//rlabplus extensions, Marijan Kostrun (C) 2005
//
    // mmio
extern Ent *ent_readmm          (int nargs, Datum args[]);
extern Ent *ent_writemm         (int nargs, Datum args[]);
    // ode toolkit
extern Ent *ent_odebv (int nargs, Datum args[]);
extern Ent *ent_odesl_eign (int nargs, Datum args[]);
extern Ent *ent_odesl_coup (int nargs, Datum args[]);
extern Ent *ent_odaei (int nargs, Datum args[]);
    // signal processing kit: gcv, stl
extern Ent *ent_gcv_init (int nargs, Datum args[]);
extern Ent *ent_gcv_val (int nargs, Datum args[]);
extern Ent *ent_stl2 (int nargs, Datum args[]);
extern Ent *ent_harminv(int nargs, Datum args[]);
    // rlabplus extensions: chaos.c
extern Ent *ent_recurrence_map (int nargs, Datum args[]);
extern Ent *ent_xcorr (int nargs, Datum args[]);
extern Ent *ent_runvar (int nargs, Datum args[]);
extern Ent *ent_avgmutinfo (int nargs, Datum args[]);
extern Ent *ent_falsenn (int nargs, Datum args[]);
extern Ent *ent_poincare (int nargs, Datum args[]);

#ifdef HAVE_CRYPTO_HASH
    // rlabplus extension: crypto hash
extern Ent *ent_openssl_hash (int nargs, Datum args[]);
#endif
extern Ent *ent_crc32 (int nargs, Datum args[]);

//
// rlabplus PUBLIC (GPL and such) extensions:
//
// rlabplus extensions: GNU scientific library
#include "rlabplus_gsl.h"
// rlabplus extensions: CLAWPACK
#include "rlabplus_claw.h"
// rlabplus extensions: PDECOL
#include "rlabplus_pde.h"
//  rlabplus extensions: ODRPACK
#include "rlabplus_odr.h"
//  rlabplus extensions: ARPACK
#include "rlabplus_arpack.h"
//  rlabplus extensions: pzeros
#include "rlabplus_pz.h"
//  rlabplus extensions: dloess
#include "rlabplus_dloess.h"
//  rlabplus extensions: dierckx
#include "rlabplus_dierckx.h"
//  rlabplus extension: strings, time, display, misc
#include "getline.h"
// asa
extern Ent * ent_asa_min (int nargs, Datum args[]);

// gnuplot
#include "r_gnuplot.h"
#define RLAB_PLOT_SUPPORT "none"

#ifdef HAVE_RLAB_PLPLOT
#include "r_plplot.h"
#undef  RLAB_PLOT_SUPPORT
#define RLAB_PLOT_SUPPORT "plplot"
#endif

#ifdef HAVE_RLAB_PGPLOT
#include "r_pgplot.h"
#undef  RLAB_PLOT_SUPPORT
#define RLAB_PLOT_SUPPORT "pgplot"
#endif


//
// end of rlabplus extensions
//


#include "bltin.h"
#include "bltin1.h"
#include "bltin2.h"
#include "bltin3.h"
#include "bltin4.h"
#include "bltin5.h"
#include "rfft.h"
#include "rfileio.h"
#include "rdl.h"

//
// sparse matrix operations
//
#include "msrf2.h"
#include "mscf2.h"

// sep - ask me why.
#include "r_sep.h"

// more on curl
#include "r_curl.h"

#include "r_sys.h"

#ifdef HAVE_GLPK
#include "rlabplus_glpk.h"
#endif

Bltin rlab_bltin[] = {
  {BLTIN, "isvector", ent_IsVec},
  {BLTIN, "isstring", ent_IsString},
  {BLTIN, "ismatrix", ent_IsMat},
  {BLTIN, "isscalar", ent_IsScalar},
  {BLTIN, "isnumber", ent_IsNumber},
  {BLTIN, "isempty", ent_IsEmpty},
#ifdef HAVE_GLPK
  {BLTIN, "_glpk_read",  ent_glpk_read_file},
  {BLTIN, "_glpk_write", ent_glpk_write_file},
  {BLTIN, "_glpk_solve", ent_glpk_solve_lp},
#endif
  {BLTIN, "all", All},
  {BLTIN, "minmax", MinMax},
  {BLTIN, "max", Max},
  {BLTIN, "maxi", MaxI},
  {BLTIN, "maxi2", Max2I},
  {BLTIN, "min", Min},
  {BLTIN, "mini", MinI},
  {BLTIN, "mini2", Min2I},
  {BLTIN, "abs", Abs},
  {BLTIN, "entinfo", EntInfo},
  {BLTIN, "fsymtab", FixSymTable},
  {BLTIN, "any", Any},
  {BLTIN, RLAB_MEMBER_CLASS, Group},
  {BLTIN, "error", Error},
  {BLTIN, "stop", Stop},
  {BLTIN, "eval", Eval},
  {BLTIN, "exist", Exist},
  {BLTIN, "members", Members},
  {BLTIN, "printf", Printf},
  {BLTIN, "fprintf", FPrintf},
  {BLTIN, "sprintf", SPrintf},
  {BLTIN, "open", Open},
  {BLTIN, "close", Close},
  {BLTIN, "lof", ListOpenFiles},
  {BLTIN, "getline", Getline},
  {BLTIN, "scanf", Scanf},
  {BLTIN, "argv", Argv},
  {BLTIN, "isinteractive", Interactive},
  {BLTIN, "basename", Basename},
  {BLTIN, "assign", ent_Assign},
  {BLTIN, "unpack", ent_Unpack},
  {BLTIN, "curr_file", ent_Filename}, 
#ifdef HAVE_HDF5_SO
  {BLTIN, "h5isfile", ent_hdf5_isfile},
  {BLTIN, "h5ls", ent_h5ls},
  {BLTIN, "h5mv", ent_h5mv},
  {BLTIN, "h5cp", ent_h5cp},
  {BLTIN, "h5ln", ent_h5ln},
#endif
  {BLTIN, "size", Size},
  {BLTIN, "length", Length},
  {BLTIN, "type", Type},
  {BLTIN, "load", Load},
#if 0
  {BLTIN, "debug", Debug},
#endif
  {BLTIN, "mod", Mod},
  {BLTIN, "reshape", Reshape},
  {BLTIN, "resize", ent_reshape_resize},
  {BLTIN, "int", Int},
  {BLTIN, "ceil", Ceil},
  {BLTIN, "floor", Floor},
  {BLTIN, "round", Round},
  {BLTIN, "issymm", IsSymm},
  {BLTIN, "mnorm", MNorm},
  {BLTIN, "rownorm", RowNorm},
  {BLTIN, "zeros", Zeros},
  {BLTIN, "ones", Ones},
  {BLTIN, "inf", Inf},
  {BLTIN, "nan", Nan},
  {BLTIN, "isinf", IsInf},
  {BLTIN, "isnan", IsNan},
  {BLTIN, "sin", Sin},
  {BLTIN, "cos", Cos},
  {BLTIN, "tan", Tan},
  {BLTIN, "asin", ASin},
  {BLTIN, "acos", ACos},
  {BLTIN, "atan", ATan},
  {BLTIN, "sqrt", Sqrt},
  {BLTIN, "log", Log},
  {BLTIN, "log10", Log10},
  {BLTIN, "exp", Exp},
  {BLTIN, "givens", ent_blas_givens},
  {BLTIN, "mexp", Mexp},
  {BLTIN, "mpow", Mpow},
  {BLTIN, "diag", Diag},
  {BLTIN, "clear", Clear},
  {BLTIN, "eig", Eig},
  {BLTIN, "sum", Sum},
  {BLTIN, "cumsum", CumSum},
  {BLTIN, "real", Real},
  {BLTIN, "imag", Imag},
  {BLTIN, "conj", Conj},
  {BLTIN, "fft", FFT},
  {BLTIN, "ifft", IFFT},
  {BLTIN, "window", ent_ctftbx_window},
  {BLTIN, "stft", ent_ctftbx_stft},
  {BLTIN, "fftshift", ent_ctftbx_fftshift},
  {BLTIN, "tframbfun", ent_ctftbx_af},
  {BLTIN, "tfrkernel", ent_ctftbx_kernel},
  {BLTIN, "af2tfr", ent_ctftbx_af2tfr},
  {BLTIN, "mkfilter", ent_mkfilter},
  {BLTIN, "find", Find},
  {BLTIN, "findvec", FindVector},
  {BLTIN, "sort", Sort},
  {BLTIN, "unique", VectorSet},
  {BLTIN, "union", VectorUnion},
  {BLTIN, "intersect", VectorIntersect},
  {BLTIN, "setdiff", VectorComplement},
  {BLTIN, "merge", Merge},
  {BLTIN, "compact", Compact},
  {BLTIN, "strtod", Strtod},
  {BLTIN, "strtol", Strtol},
  {BLTIN, "factor", Factor},
  {BLTIN, "backsub", Backsub},
  {BLTIN, "vpnorm", PNorm},
  {BLTIN, "sign", Sign},
  {BLTIN, "svd", Svd},
  {BLTIN, "readm", ReadM},
  {BLTIN, "writem", WriteM},
  {BLTIN, "readb", ReadB},
  {BLTIN, "writeb", WriteB},
  {BLTIN, "read_ascii", ReadASCII},
  {BLTIN, "write_ascii", WriteASCII},
  {BLTIN, "chomp", ent_Chomp},
  {BLTIN, "flipud", ent_FlipUD},
  {BLTIN, "fliplr", ent_FlipLR},
  {BLTIN, "shiftu", ent_ShiftU},
  {BLTIN, "shiftd", ent_ShiftD},
  {BLTIN, "shiftl", ent_ShiftL},
  {BLTIN, "shiftr", ent_ShiftR},
  {BLTIN, "rot90", ent_Rot90},
  {BLTIN, "rotx", ent_rotate_x},
  {BLTIN, "roty", ent_rotate_y},
  {BLTIN, "rotz", ent_rotate_z},
// some internal functions
  {BLTIN, "copy", ent_copyentity},
//
// rlabplus extensions: BEGIN
//
  {BLTIN, "asamin",   ent_asa_min},
  {BLTIN, "diffevol", ent_gsl_diffevol},
#ifdef HAVE_CRYPTO_HASH
  {BLTIN, "hash", ent_openssl_hash},
#endif
  {BLTIN, "crc32", ent_crc32},
  // string et c.
  {BLTIN, "grep", Grep},
  {BLTIN, "lstrip", Lstrip},
  {BLTIN, "rstrip", Rstrip},
  {BLTIN, "strip", Strip},
  {BLTIN, "strlen", Strlen},
  {BLTIN, "strindex", ent_string_index},
  {BLTIN, "strsplt", Strsplt},
  {BLTIN, "findstr", Findstr},
  {BLTIN, "strindexr", ent_string_index_pattern},
  {BLTIN, "substr", ent_string_substr},
  {BLTIN, "reads", ReadS},
  {BLTIN, "spinner", ent_spinner},
  {BLTIN, "smiley", ent_smiley},
  {BLTIN, "tolower", ent_string_tolower},
  {BLTIN, "toupper", ent_string_toupper},
  {BLTIN, "capitalize", ent_string_capitalize},
  {BLTIN, "num2str", ent_string_text},
  {BLTIN, "blank", ent_string_blank},
  {BLTIN, "ascii", ent_string_ascii},
  {BLTIN, "char", ent_string_char},
  {BLTIN, "gsub", ent_string_substitute},
  {BLTIN, "ifelse", ent_ifelse},
  // time functions
  {BLTIN, "gmtime", ent_gmtime},
  {BLTIN, "clock", ent_Clock},
  {BLTIN, "seconds", ent_Seconds},
  {BLTIN, "time2dstr", ent_strftime},
  {BLTIN, "etime", ent_Etime},
  {BLTIN, "dayofweek", ent_dayofweek},
  {BLTIN, "dayofyear", ent_dayofyear},
  {BLTIN, "dstr2time", ent_dstr2time},
  // file and locking functions
  {BLTIN, "stat", ent_Stat},
  {BLTIN, "getpid", ent_Getpid},
  {BLTIN, "fnctl_lock", ent_Fnctl},
  // entropy
  {BLTIN, "_iseed", ent_iseed},
  // system
  {BLTIN, "lsdir", ent_LsDir},
  // terminal
  {BLTIN, "openpty", ent_openpty_test},
  {BLTIN, "clrscr", ent_termcap_clrscr},
  {BLTIN, "mvcrsr", ent_termcap_mvcrsr},
  {BLTIN, "clrpos", ent_termcap_clrpos},
  {BLTIN, "colors", ent_termcap_colors},
  {BLTIN, "showkey", ent_showkey},
  // mmio
  {BLTIN, "readmm",  ent_readmm},
  {BLTIN, "writemm", ent_writemm},
  // tisean et c.
  {BLTIN, "rmap", ent_recurrence_map},
  {BLTIN, "xcorr", ent_xcorr},
  {BLTIN, "ami", ent_avgmutinfo},
  {BLTIN, "falsenn", ent_falsenn},
  {BLTIN, "poincare", ent_poincare},
  {BLTIN, "runvar", ent_runvar},
  // signal et c.
  {BLTIN, "bsplinefit", ent_bspline_fit},
  {BLTIN, "bsplineval", ent_bspline_val},
  {BLTIN, "bsplinefit2", ent_bspline_fit2},
  {BLTIN, "bsplineval2", ent_bspline_val2},
  {BLTIN, "gcvsplinefit", ent_gcv_init},
  {BLTIN, "gcvsplineval", ent_gcv_val},
  {BLTIN, "stl", ent_stl2},
  {BLTIN, "harminv", ent_harminv},
  // ode
  {BLTIN, "odeiv", ent_gsl_odeiv},
  {BLTIN, "odebv", ent_odebv},
  {BLTIN, "odaei", ent_odaei},
  // pde
  {BLTIN, "claw1", ent_claw1_basic},
  {BLTIN, "pdecol", ent_pdecol},
  {BLTIN, "bacol", ent_bacol},
  // The GSL and related functions
  {BLTIN, "dwt", ent_gsl_dwt},
  {BLTIN, "idwt", ent_gsl_idwt},
  {BLTIN, "dwt2", ent_gsl_dwt2},
  {BLTIN, "idwt2", ent_gsl_idwt2},
  {BLTIN, "levinu", ent_gsl_levinu},
  {BLTIN, "siman", ent_gsl_siman},
  {BLTIN, "cspleval", ent_cubicspline_eval},
  {BLTIN, "csplinterp", ent_cubicspline_init},
  {BLTIN, "linterp", ent_lininterp_eval},
  {BLTIN, "linterp1", ent_lininterp1_eval},
  {BLTIN, "harmsum", ent_harmsum_eval},
  {BLTIN, "ndiffs", ent_ndiffs},
  {BLTIN, "ndiv", ent_div},
  {BLTIN, "lsfit", ent_gsl_lsfit},
  {BLTIN, "polyfit", ent_polyfit},
  {BLTIN, "odrfit", ent_odr_basic},
  {BLTIN, "findmin", ent_minimize_1d_bracket},
  {BLTIN, "findmins", ent_gsl_findmins},
  {BLTIN, "conmins", ent_conmax},
  {BLTIN, "findroot", ent_findroot},
  {BLTIN, "findroots", ent_gsl_findroots},
  {BLTIN, "crvtrack", ent_ctrack},
  {BLTIN, "mcnint", ent_nintegrate_mc},
  {BLTIN, "tchebyfit", ent_tchebyshev_init},
  {BLTIN, "tchebyval", ent_tchebyshev_eval},
  {BLTIN, "tchebyder", ent_tchebyshev_diff},
  {BLTIN, "tchebyint", ent_tchebyshev_int},
  {BLTIN, "ndiff", ent_ndiff},
  {BLTIN, "covar", ent_covariance},
  {BLTIN, "mean", ent_mean},
  {BLTIN, "var", ent_variance},
  {BLTIN, "absdev", ent_absdev},
  {BLTIN, "skew", ent_skew},
  {BLTIN, "kurtosis", ent_kurtosis},
  {BLTIN, "median", ent_median},
  {BLTIN, "quantile", ent_quantile},
  {BLTIN, "invperm", ent_invperm},
  {BLTIN, "revperm", ent_revperm},
  {BLTIN, "prevperm", ent_prevperm},
  {BLTIN, "nextperm", ent_nextperm},
  {BLTIN, "validperm", ent_isaperm},
  {BLTIN, "mulperm", ent_mulperm},
  {BLTIN, "prevcomb", ent_prevcomb},
  {BLTIN, "nextcomb", ent_nextcomb},
  {BLTIN, "validcomb", ent_iscomb},
  {BLTIN, "ncomb", ent_ncomb},
  {BLTIN, "clip", ent_gsl_clip},
  {BLTIN, "hist", ent_histogram_create},
  {BLTIN, "hist2", ent_histogram2d_create},
  {BLTIN, "nintegrate", ent_nintegrate},
  {BLTIN, "nintegrates", ent_nintegrate_qagp},
  {BLTIN, "cauchypval", ent_nintegrate_qawc},
  {BLTIN, "nintqaws", ent_nintegrate_qaws},
  {BLTIN, "nintsin", ent_nintegrate_qawo_sin},
  {BLTIN, "nintcos", ent_nintegrate_qawo_cos},
  {BLTIN, "nintmd", ent_nintegrate_genz_hypercube},
  {BLTIN, "nintsimplex", ent_nintegrate_genz_simplex},
  {BLTIN, "nintconv", ent_convolution},
  {BLTIN, "nintlt", ent_laplacetransform},
  {BLTIN, "invlt",  ent_invlt},
  {BLTIN, "AiryAi", ent_gsl_sf_airy_Ai},
  {BLTIN, "AiryBi", ent_gsl_sf_airy_Bi},
  {BLTIN, "AiryAiPrime", gsl_sf_airy_Ai_deriv},
  {BLTIN, "AiryBiPrime", gsl_sf_airy_Bi_deriv},
  {BLTIN, "AiryZeroAi", ent_gsl_sf_airy_zero_Ai},
  {BLTIN, "AiryZeroBi", ent_gsl_sf_airy_zero_Bi},
  {BLTIN, "AiryZeroAiPrime", ent_gsl_sf_airy_zero_Ai_deriv},
  {BLTIN, "AiryZeroBiPrime", ent_gsl_sf_airy_zero_Bi_deriv},
  {BLTIN, "BesselI", ent_gsl_sf_bessel_I},
  {BLTIN, "BesselJ", ent_gsl_sf_bessel_J},
  {BLTIN, "BesselK", ent_gsl_sf_bessel_K},
  {BLTIN, "BesselY", ent_gsl_sf_bessel_Y},
  {BLTIN, "Besseli", ent_gsl_sf_bessel_i},
  {BLTIN, "Besselj", ent_gsl_sf_bessel_j},
  {BLTIN, "Besselk", ent_gsl_sf_bessel_k},
  {BLTIN, "Bessely", ent_gsl_sf_bessel_y},
  {BLTIN, "BesselZeroJ", ent_gsl_sf_bessel_zeroJ},
  {BLTIN, "Clausen", ent_gsl_sf_clausen},
  {BLTIN, "Hydrogen", ent_gsl_sf_hydrogenicR},
  {BLTIN, "Dawson", ent_gsl_sf_dawson},
  {BLTIN, "Debye", ent_gsl_sf_debye},
  {BLTIN, "EllipticE", ent_gsl_sf_elliptic_E},
  {BLTIN, "EllipticF", ent_gsl_sf_ellint_F},
  {BLTIN, "EllipticK", ent_gsl_sf_ellint_Kcomp},
  {BLTIN, "EllipticP", ent_gsl_sf_ellint_P},
  {BLTIN, "EllipticD", ent_gsl_sf_ellint_D},
  {BLTIN, "EllipticRC", ent_gsl_sf_ellint_RC},
  {BLTIN, "EllipticRD", ent_gsl_sf_ellint_RD},
  {BLTIN, "EllipticRF", ent_gsl_sf_ellint_RF},
  {BLTIN, "EllipticRJ", ent_gsl_sf_ellint_RJ},
  {BLTIN, "EllipticJacobi", ent_gsl_sf_elliptic_jacobi},
  {BLTIN, "Erf", ent_gsl_sf_erf},
  {BLTIN, "InvErf", ent_gsl_sf_inverf},
  {BLTIN, "Erfc", ent_gsl_sf_erfc},
  {BLTIN, "ErfQ", ent_gsl_sf_erfq},
  {BLTIN, "ErfZ", ent_gsl_sf_erfz},
  {BLTIN, "ExpIntegralE", ent_gsl_sf_expint},
  {BLTIN, "ExpIntegralEi", ent_gsl_sf_expint_ei},
  {BLTIN, "SinhIntegral", ent_gsl_sf_Shi},
  {BLTIN, "CoshIntegral", ent_gsl_sf_Chi},
  {BLTIN, "SinIntegral", ent_gsl_sf_Si},
  {BLTIN, "CosIntegral", ent_gsl_sf_Ci},
  {BLTIN, "AtanIntegral", ent_gsl_sf_atanint},
  {BLTIN, "FermiDiracIntegral", ent_gsl_sf_fermidiracint},
  {BLTIN, "LogGamma", ent_gsl_sf_loggamma},
  {BLTIN, "Gamma", ent_gsl_sf_gamma},
  {BLTIN, "GammaRegularized", ent_gsl_sf_gamma_inc_Q},
  {BLTIN, "GammaRegularizedC", ent_gsl_sf_gamma_inc_P},
  {BLTIN, "RecGamma", ent_gsl_sf_gammainv},
  {BLTIN, "LogBeta", ent_gsl_sf_lnbeta},
  {BLTIN, "Beta", ent_gsl_sf_beta},
  {BLTIN, "BetaRegularized", ent_gsl_sf_beta_inc},
  {BLTIN, "Pochhammer", ent_gsl_sf_poch},
  {BLTIN, "GegenbauerC", ent_gsl_sf_gegenpoly},
  {BLTIN, "Hypergeometric0F1", ent_gsl_sf_hyperg_0F1},
  {BLTIN, "Hypergeometric1F1", ent_gsl_sf_hypergeom_1F1},
  {BLTIN, "HypergeometricU", ent_gsl_sf_hypergeom_U},
  {BLTIN, "Hypergeometric2F0", ent_gsl_sf_hyperg_2F0},
  {BLTIN, "Hypergeometric2F1", ent_gsl_sf_hyperg_2F1},
  {BLTIN, "LaguerreL", ent_gsl_sf_laguerre},
  {BLTIN, "ProductLog", ent_gsl_sf_lambert_W0},
  {BLTIN, "LegendreP", ent_gsl_sf_legendreP},
  {BLTIN, "LegendreQ", ent_gsl_sf_legendreQ},
  {BLTIN, "LegendreSphericalP", ent_gsl_sf_legendre_sphPlm},
  {BLTIN, "LegendreConicalP", ent_gsl_sf_conicalP},
  {BLTIN, "LegendreH3d", ent_gsl_sf_legendreH3d},
  {BLTIN, "Digamma", ent_gsl_sf_digamma},
  {BLTIN, "PolyGamma", ent_gsl_sf_psi_n},
  {BLTIN, "SynchrotronF", ent_gsl_sf_synchrotron},
  {BLTIN, "TransportF", ent_gsl_sf_transport},
  {BLTIN, "Zeta", ent_gsl_sf_zeta},
  {BLTIN, "ErlingN", ent_sf_erling},
  {BLTIN, "irng", ent_gsl_irng},
  {BLTIN, "irand", ent_gsl_irand},
  {BLTIN, "irng_state", ent_gsl_irng_state},
  {BLTIN, "rng", ent_gsl_rng},
  {BLTIN, "rand", ent_gsl_rand},
  {BLTIN, "randomize", ent_gsl_randomize},
  {BLTIN, "uniform", ent_gsl_urand},
  {BLTIN, "gaussian", ent_gsl_gauss},
  {BLTIN, "shuffle", ent_gsl_shuffle},
  {BLTIN, "sample", ent_gsl_sample},
  {BLTIN, "drng", ent_gsl_setdrng},
  {BLTIN, "drand", ent_gsl_drand},
  {BLTIN, "hrng", ent_gsl_sethrng},
  {BLTIN, "hrand", ent_gsl_hrand},
  {BLTIN, "protect", btree_protect},
  {BLTIN, "release", btree_release},
  {BLTIN, "isprot", btree_IsConst},
  // math
  {BLTIN, "polyder", ent_poly_eval_diff},
  {BLTIN, "polyint", ent_poly_eval_int},
  {BLTIN, "polyval", ent_poly_eval},
  {BLTIN, "cosh",    ent_Cosh},
  {BLTIN, "sinh",    ent_Sinh},
  {BLTIN, "tanh",    ent_Tanh},
  {BLTIN, "acosh",   ent_Acosh},
  {BLTIN, "asinh",   ent_Asinh},
  {BLTIN, "atanh",   ent_Atanh},
  // 3d centered stuff
  {BLTIN, "rot3d",   ent_3d_rotate},
  {BLTIN, "bwlabel",   ent_bwlabel},
  {BLTIN, "bwblobs",   ent_bwblobs},
  // pzeros
  {BLTIN, "polyroots", ent_pzeros},
  {BLTIN, "nseries", ent_nseries},
  // sparse
  {BLTIN, "spsolve", SpSolve_BF},
  {BLTIN, "spwrite", SpWrite},
  // isolated eigenvalues
#ifdef HAVE_ARPACK
  {BLTIN, "eigs", Eigs  },
#endif
  // gnuplot interface
  {BLTIN, "gnustart", ent_gnuplot_init },
  {BLTIN, "gnuclose", ent_gnuplot_close },
  {BLTIN, "gnuwin", ent_gnuplot_win },
  {BLTIN, "_gnuwins", ent_gnuplot_gnuwins},
  {BLTIN, "gnucmd", ent_gnuplot_cmd },
  {BLTIN, "gnuprint", ent_gnuplot_print },
  // data processing
  {BLTIN, "locextri", localExtremaeIdx },
  {BLTIN, "isabsmono", isAbsMonotone },
  {BLTIN, "isrelmono", isRelMonotone },
  {BLTIN, "distance",    ent_distance_matrix},
  {BLTIN, "rank",    ent_rank},
  {BLTIN, "sammon",  ent_sprannlib_sammon},
  {BLTIN, "ffnneval", ent_ffnn_eval},
  {BLTIN, "samann",    ent_ffnn_samann},
  {BLTIN, "ffnn",    ent_ffnn},
  // bitwise
  {BLTIN, "bitshiftl", ent_bit_shift_left},
  {BLTIN, "bitshiftr", ent_bit_shift_right},
  {BLTIN, "bytesplit", ent_byte_split},
  {BLTIN, "bytejoin", ent_byte_join},
  {BLTIN, "bitsplit", ent_bit_split},
  // misc, yes?
  {BLTIN, "conrec", ent_Contour},
  //
  // rlabplus extensions: END
  //
  {BLTIN, "solve", Solve},
  {BLTIN, "hess", Hess},
  {BLTIN, "balance", Balance},
  {BLTIN, "qr", QR},
  {BLTIN, "filter", Filter},
  {BLTIN, "sylv", Sylv},
  {BLTIN, "sizeof", Sizeof},
  {BLTIN, "system", System},
  {BLTIN, "fork", Fork},
  {BLTIN, "prod", Prod},
  {BLTIN, "cumprod", CumProd},
  {BLTIN, "frexp", Frexp},
  {BLTIN, "ldexp", Ldexp},
  {BLTIN, "tic", Tic},
  {BLTIN, "toc", Toc},
  {BLTIN, "format", Format},
  {BLTIN, "diary", Diary},
#ifdef HAVE_GETENV
  {BLTIN, "getenv", Getenv},
#endif
  {BLTIN, "cd", Cd},
#ifdef HAVE_PUTENV
  {BLTIN, "putenv", Putenv},
#endif
  {BLTIN, "sleep", Sleep},
//  {BLTIN, "tmpnam", Tmpnam},
  {BLTIN, "rcond", Rcond},
  {BLTIN, "schur", Schur},
#ifdef HAVE_SO			/* Don't try this if shared objects are not supported. */
  {BLTIN, "dlopen", DLopen},
#endif /* HAVE_SO */
  {BLTIN, "finite", Finite},
  {BLTIN, "fread", Fread},
  {BLTIN, "fwrite", Fwrite},
  {BLTIN, "fseek", Fseek},
  {BLTIN, "atan2", ATan2},
  {BLTIN, "chol", Chol},
  {BLTIN, "det", Det},
  {BLTIN, "sparse", Sparse},
  {BLTIN, "full", Dense},
  {BLTIN, "spfactor", SpFactor_BF},
  {BLTIN, "spconvert", Spconvert},
  {BLTIN, "sporder", SpOrder},
  {BLTIN, "readgraph", ReadGraph},
#ifdef HAVE_LOGB
  {BLTIN, "logb", Logb},
#endif /* HAVE_LOGB */
#ifdef HAVE_RLAB_PLPLOT
  {BLTIN, "_plsetopts", ent_plplot_plsetopt},
  {BLTIN, "_plprintf", ent_plplot_printf},
  {BLTIN, "_plprintf3", ent_plplot3_printf},
  {BLTIN, "_plchr", ent_plplot_plchr},
  {BLTIN, "_plgcol0", ent_plplot_plgcol0},
  {BLTIN, "_plmkstrm", ent_plplot_plmkstrm},
  {BLTIN, "_pllegend", ent_plplot_pllegend},
  {BLTIN, "_plfill", ent_plplot_plfill},
  {BLTIN, "_plprint", _plot_plprint},
  {BLTIN, "_replot", _plot_plreplot},
  {BLTIN, "_plsstrm", _plot_plsstrm},
  {BLTIN, "_plssub", _plot_plssub},
  {BLTIN, "_plinit", _plot_plinit},
  {BLTIN, "_plstart", _plot_plstart},
  {BLTIN, "_plsfile", _plot_plsfile},
  {BLTIN, "_plenv", _plot_plenv},
  {BLTIN, "_plline", _plot_plline},
  {BLTIN, "_plline3", _plot_plline3},
  {BLTIN, "_plpoly3", _plot_plpoly3},
  {BLTIN, "_plend", _plot_plend},
  {BLTIN, "_plend1", _plot_plend1},
  {BLTIN, "_pllab", _plot_pllab},
  {BLTIN, "_plcol", _plot_plcol},
  {BLTIN, "_plscolbg", _plot_plscolbg},
  {BLTIN, "_plscol0", _plot_plscol0},
  {BLTIN, "_pllsty", _plot_pllsty},
  {BLTIN, "_plclr", _plot_plclr},
  {BLTIN, "_plpoin", _plot_plpoin},
  {BLTIN, "_plpoin3", _plot_plpoin3},
  {BLTIN, "_plhist", _plot_plhist},
  {BLTIN, "_plspage", _plot_plspage},
  {BLTIN, "_pladv", _plot_pladv},
  {BLTIN, "_plgra", _plot_plgra},
  {BLTIN, "_pltext", _plot_pltext},
  {BLTIN, "_plflush", _plot_plflush},
  {BLTIN, "_plbox", _plot_plbox},
  {BLTIN, "_plbox3", _plot_plbox3},
  {BLTIN, "_plvsta", _plot_plvsta},
  {BLTIN, "_plvasp", _plot_plvasp},
  {BLTIN, "_plwind", _plot_plwind},
  {BLTIN, "_plot3d", _plot_plot3d},
  {BLTIN, "_plmesh", _plot_plmesh},
  {BLTIN, "_plw3d", _plot_plw3d},
  {BLTIN, "_plmtex", _plot_plmtex},
  {BLTIN, "_plspause", _plot_plspause},
  {BLTIN, "_plwid", _plot_plwid},
  {BLTIN, "_plptex", _plot_plptex},
  {BLTIN, "_plfontld", _plot_plfontld},
  {BLTIN, "_plfont", _plot_plfont},
  {BLTIN, "_plsori", _plot_plsori},
  {BLTIN, "_plscolor", _plot_plscolor},
  {BLTIN, "_plcont", _plot_plcont},
  {BLTIN, "_plvpor", _plot_plvpor},
  {BLTIN, "_plvpas", _plot_plvpas},
  {BLTIN, "_plerry", _plot_plerry},
  {BLTIN, "_plssym", _plot_plssym},
  {BLTIN, "_plpsty", _plot_plpsty},
  {BLTIN, "_plprec", ent_plplot_plprec},
  {BLTIN, "_plstring", ent_plplot_plstring},
  {BLTIN, "_plstring3", ent_plplot_plstring3},
  {BLTIN, "_plstyl", ent_plplot_plstyl},
  {BLTIN, "_plspal1", ent_plplot_plspal1},
  {BLTIN, "_plspal0", ent_plplot_plspal0},
  {BLTIN, "_plsmin", ent_plplot_plsmin},
  {BLTIN, "_plsmaj", ent_plplot_plsmaj},
  {BLTIN, "_plshades", ent_plplot_plshades},
#endif /* HAVE_RLAB_PLPLOT */
#ifdef HAVE_RLAB_PGPLOT
  {BLTIN, "_pgdrawgrad", ent_pgplot_draw_grad_plane},
  {BLTIN, "_pgprintf", ent_pgplot_printf},
  {BLTIN, "_pgqcr", ent_pgplot_pgqcr},
  {BLTIN, "_pgcolor", _pg_color},
  {BLTIN, "pglen", _pg_cpglen},
  {BLTIN, "pgqwin", _pg_pgqwin},
  {BLTIN, "pgqvp", _pg_cpgqvp},
  {BLTIN, "pgvsiz", _pg_cpgvsiz},
  {BLTIN, "pgbeg", _pg_cpgbeg},
  {BLTIN, "pgopen", _pg_cpgopen},
  {BLTIN, "pgsvp", _pg_cpgsvp},
  {BLTIN, "pgvstd", _pg_cpgvstd},
  {BLTIN, "pgswin", _pg_cpgswin},
  {BLTIN, "pgwnad", _pg_cpgwnad},
  {BLTIN, "pgslct", _pg_cpgslct},
  {BLTIN, "pgsubp", _pg_cpgsubp},
  {BLTIN, "pgpanl", _pg_cpgpanl},
  {BLTIN, "pgpap", _pg_cpgpap},
  {BLTIN, "pgenv", _pg_cpgenv},
  {BLTIN, "pgline", _pg_cpgline},
  {BLTIN, "pgpt", _pg_cpgpt},
  {BLTIN, "pgimag", _pg_cpgimag},
  {BLTIN, "pggray", _pg_cpggray},
  {BLTIN, "pgwedg", _pg_cpgwedg},
  {BLTIN, "pgconb", _pg_cpgconb},
  {BLTIN, "pgconl", _pg_cpgconl},
  {BLTIN, "pgcont", _pg_cpgcont},
  {BLTIN, "pgconf", _pg_cpgconf},
  {BLTIN, "pgbox", _pg_cpgbox},
  {BLTIN, "pgtbox", _pg_cpgtbox},
  {BLTIN, "pgstbg", _pg_cpgstbg},
  {BLTIN, "pgsls", _pg_cpgsls},
  {BLTIN, "pgslw", _pg_cpgslw},
  {BLTIN, "pgend", _pg_cpgend},
  {BLTIN, "pgclos", _pg_cpgclos},
  {BLTIN, "pgpage", _pg_cpgpage},
  {BLTIN, "pgeras", _pg_cpgeras},
  {BLTIN, "pgetxt", _pg_cpgetxt},
  {BLTIN, "pgtext", _pg_cpgtext},
  {BLTIN, "pgptxt", _pg_cpgptxt},
  {BLTIN, "pgmtxt", _pg_cpgmtxt},
  {BLTIN, "pgscf", _pg_cpgscf},
  {BLTIN, "pgsch", _pg_cpgsch},
  {BLTIN, "pgask", _pg_cpgask},
  {BLTIN, "pglab", _pg_cpglab},
  {BLTIN, "pgupdt", _pg_cpgupdt},
  {BLTIN, "pgsci", _pg_cpgsci},
  {BLTIN, "pgscir", _pg_cpgscir},
  {BLTIN, "pgqci", _pg_cpgqci},
  {BLTIN, "pgqcir", _pg_cpgqcir},
  {BLTIN, "pgbbuf", _pg_cpgbbuf},
  {BLTIN, "pgqcol", _pg_cpgqcol},
  {BLTIN, "pgmove", _pg_cpgmove},
  {BLTIN, "pgdraw", _pg_cpgdraw},
  {BLTIN, "pgsfs", _pg_cpgsfs},
  {BLTIN, "pgshs", _pg_cpgshs},
  {BLTIN, "pgpoly", _pg_cpgpoly},
  {BLTIN, "pgbin", _pg_cpgbin},
  {BLTIN, "pghist", _pg_cpghist},
  {BLTIN, "pgrnge", _pg_cpgrnge},
  {BLTIN, "pgsave", _pg_cpgsave},
  {BLTIN, "pgunsa", _pg_cpgunsa},
  {BLTIN, "pgrect", _pg_cpgrect},
#endif /* HAVE_RLAB_PGPLOT */
  // MDE's are here
  {BLTIN, "cell", Cell},
  {BLTIN, "vcell", VCell},
  {BLTIN, "hcell", HCell},
  {BLTIN, "iscell", IsCell},
  {0, 0, 0},
};

#include "gsl_pdf.h"

//
// odesl(), list of sturm-liouville solvers
//
Bltin rlab_odesl_bltin[] = {
  {BLTIN, "eign", ent_odesl_eign},
  {BLTIN, "coup", ent_odesl_coup},
  {0, 0, 0},
};


//
// sparams(), auxiliary list of functions for sparse solvers
//
Bltin rlab_sparams_bltin[] = {
  {BLTIN, "realsolv", ent_sparse_realsolv},
  {BLTIN, "compsolv", ent_sparse_compsolv},
  {BLTIN, "iterator", ent_sparse_iter},
  {BLTIN, "precond", ent_sparse_prec},
  {BLTIN, "tol", ent_sparse_tol},
  {BLTIN, "dpt", ent_sparse_dpt},
  {0, 0, 0},
};

//
// clawparams(), auxiliary list of functions for  claw1  solver
//
Bltin rlab_clawparams_bltin[] = {
  {BLTIN, "bc", ent_claw1par_bc},
  {BLTIN, "method", ent_claw1par_method},
  {BLTIN, "cour", ent_claw1par_courmaxit},
  {BLTIN, "b4step1", ent_claw1par_b4s1},
  {BLTIN, "tsrc", ent_claw1par_tsrc},
  {BLTIN, "lim", ent_claw1par_lim},
  {0, 0, 0},
};

//
// loess list
//
Bltin rlab_loess_bltin[] = {
  // dloess
  {BLTIN, "setup",   ent_loess_init},
  {BLTIN, "main",    ent_loess},
  {BLTIN, "predict", ent_loess_predict},
  {BLTIN, "anova",   ent_loess_anova},
  {0, 0, 0},
};

//
// classify
//
Bltin rlab_classify_bltin[] = {
  {BLTIN, "parzen",  ent_classify_parzen },
  {BLTIN, "knn",     ent_classify_knn},
//   {BLTIN, "maha",    ent_sprannlib_classify_maha},
  {BLTIN, "fisher",  ent_sprannlib_classify_fisher},
  {BLTIN, "fisherq", ent_sprannlib_classify_fisherq},
  {0, 0, 0}
};

//
// cluster
//
Bltin rlab_cluster_bltin[] = {
  {BLTIN, "iso",    ent_cluster_isodata},
  {BLTIN, "pao",    ent_cluster_pao},
  {BLTIN, "knn",    ent_cluster_nn},
  {0, 0, 0},
};

#ifdef HAVE_RLABPLUS_SO

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>

//
// GSL physical constants
//
Doubleconst gsl_mks_double[] = {
  {"c", GSL_CONST_MKSA_SPEED_OF_LIGHT},
  {"Na", GSL_CONST_NUM_AVOGADRO},
  {"G", GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT},
  {"h", GSL_CONST_MKSA_PLANCKS_CONSTANT_H},
  {"hbar", GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR},
  {"au", GSL_CONST_MKSA_ASTRONOMICAL_UNIT},
  {"ly", GSL_CONST_MKSA_LIGHT_YEAR},
  {"parsec", GSL_CONST_MKSA_PARSEC},
  {"g", GSL_CONST_MKSA_GRAV_ACCEL},
  {"eV", GSL_CONST_MKSA_ELECTRON_VOLT},
  {"me", GSL_CONST_MKSA_MASS_ELECTRON},
  {"mmu", GSL_CONST_MKSA_MASS_MUON},
  {"mp", GSL_CONST_MKSA_MASS_PROTON},
  {"mn", GSL_CONST_MKSA_MASS_NEUTRON},
  {"Ry", GSL_CONST_MKSA_RYDBERG},
  {"kb", GSL_CONST_MKSA_BOLTZMANN},
  {"mub", GSL_CONST_MKSA_BOHR_MAGNETON},
  {"mun", GSL_CONST_MKSA_NUCLEAR_MAGNETON},
  {"mue", GSL_CONST_MKSA_ELECTRON_MAGNETIC_MOMENT},
  {"mup", GSL_CONST_MKSA_PROTON_MAGNETIC_MOMENT},
  {"R0", GSL_CONST_MKSA_MOLAR_GAS},
  {"V0", GSL_CONST_MKSA_STANDARD_GAS_VOLUME},
  {"min", GSL_CONST_MKSA_MINUTE},
  {"hour", GSL_CONST_MKSA_HOUR},
  {"day", GSL_CONST_MKSA_DAY},
  {"week", GSL_CONST_MKSA_WEEK},
  {"inch", GSL_CONST_MKSA_INCH},
  {"ft", GSL_CONST_MKSA_FOOT},
  {"yd", GSL_CONST_MKSA_YARD},
  {"mile", GSL_CONST_MKSA_MILE},
  {"nmi", GSL_CONST_MKSA_NAUTICAL_MILE},
  {"fathom", GSL_CONST_MKSA_FATHOM},
  {"mil", GSL_CONST_MKSA_MIL},
  {"pt", GSL_CONST_MKSA_POINT},
  {"texpt", GSL_CONST_MKSA_TEXPOINT},
  {"micron", GSL_CONST_MKSA_MICRON},
  {"a", GSL_CONST_MKSA_ANGSTROM},
  {"ha", GSL_CONST_MKSA_HECTARE},
  {"acre", GSL_CONST_MKSA_ACRE},
  {"barn", GSL_CONST_MKSA_BARN},
  {"L", GSL_CONST_MKSA_LITER},
  {"gal", GSL_CONST_MKSA_US_GALLON},
  {"qt", GSL_CONST_MKSA_QUART},
  {"pt", GSL_CONST_MKSA_PINT},
  {"cup", GSL_CONST_MKSA_CUP},
  {"floz", GSL_CONST_MKSA_FLUID_OUNCE},
  {"Tsp", GSL_CONST_MKSA_TABLESPOON},
  {"tsp", GSL_CONST_MKSA_TEASPOON},
  {"cgal", GSL_CONST_MKSA_CANADIAN_GALLON},
  {"ukgal", GSL_CONST_MKSA_UK_GALLON},
  {"mph", GSL_CONST_MKSA_MILES_PER_HOUR},
  {"kmh", GSL_CONST_MKSA_KILOMETERS_PER_HOUR},
  {"knot", GSL_CONST_MKSA_KNOT},
  {"lb", GSL_CONST_MKSA_POUND_MASS},
  {"oz", GSL_CONST_MKSA_OUNCE_MASS},
  {"uston", GSL_CONST_MKSA_TON},
  {"ton", GSL_CONST_MKSA_METRIC_TON},
  {"ukton", GSL_CONST_MKSA_UK_TON},
  {"toz", GSL_CONST_MKSA_TROY_OUNCE},
  {"ct", GSL_CONST_MKSA_CARAT},
  {"uam", GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS},
  {"gf", GSL_CONST_MKSA_GRAM_FORCE},
  {"pf", GSL_CONST_MKSA_POUND_FORCE},
  {"kpf", GSL_CONST_MKSA_KILOPOUND_FORCE},
  {"pal", GSL_CONST_MKSA_POUNDAL},
  {"cal", GSL_CONST_MKSA_CALORIE},
  {"kcal", 1e3 * GSL_CONST_MKSA_CALORIE},
  {"tntton", 1e9 * GSL_CONST_MKSA_CALORIE},
  {"btu", GSL_CONST_MKSA_BTU},
  {"therm", GSL_CONST_MKSA_THERM},
  {"hp", GSL_CONST_MKSA_HORSEPOWER},
  {"bar", GSL_CONST_MKSA_BAR},
  {"atm", GSL_CONST_MKSA_STD_ATMOSPHERE},
  {"torr", GSL_CONST_MKSA_TORR},
  {"mHg", GSL_CONST_MKSA_METER_OF_MERCURY},
  {"inHg", GSL_CONST_MKSA_INCH_OF_MERCURY},
  {"inH2O", GSL_CONST_MKSA_INCH_OF_WATER},
  {"psi", GSL_CONST_MKSA_PSI},
  {"poise", GSL_CONST_MKSA_POISE},
  {"stokes", GSL_CONST_MKSA_STOKES},
  {"F", GSL_CONST_MKSA_FARADAY},
  {"e", GSL_CONST_MKSA_ELECTRON_CHARGE},
  {"gauss", GSL_CONST_MKSA_GAUSS},
  {"stilb", GSL_CONST_MKSA_STILB},
  {"lumen", GSL_CONST_MKSA_LUMEN},
  {"lux", GSL_CONST_MKSA_LUX},
  {"phot", GSL_CONST_MKSA_PHOT},
  {"ftcan", GSL_CONST_MKSA_FOOTCANDLE},
  {"lam", GSL_CONST_MKSA_LAMBERT},
  {"ftlam", GSL_CONST_MKSA_FOOTLAMBERT},
  {"curie", GSL_CONST_MKSA_CURIE},
  {"roe", GSL_CONST_MKSA_ROENTGEN},
  {"rad", GSL_CONST_MKSA_RAD},
  {"mSun", GSL_CONST_MKSA_SOLAR_MASS},
  {"a0", GSL_CONST_MKSA_BOHR_RADIUS},
  {"N", GSL_CONST_MKSA_NEWTON},
  {"dyne", GSL_CONST_MKSA_DYNE},
  {"J", GSL_CONST_MKSA_JOULE},
  {"erg", GSL_CONST_MKSA_ERG},
  {"kSB", GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT},
  {"tcs", GSL_CONST_MKSA_THOMSON_CROSS_SECTION},
  {"eps0", GSL_CONST_MKSA_VACUUM_PERMITTIVITY},
  {"mu0", GSL_CONST_MKSA_VACUUM_PERMEABILITY},
  {"alpha",GSL_CONST_NUM_FINE_STRUCTURE},
  // more units:
  {"mm", 1e-3},
  {"cm", 1e-2},
  {"dm", 1e-1},
  {"km", 1e3},
  {"mm2", 1e-6},
  {"cm2", 1e-4},
  {"dm2", 1e-2},
  {"km2", 1e6},
  {"mm3", 1e-9},
  {"cm3", 1e-6},
  {"dm3", 1e-3},
  {"km3", 1e9},
  {0, 0},
};

#include <gsl/gsl_math.h>
Doubleconst gsl_math_double[] = {
  {"e", M_E},
  {"log2e", M_LOG2E},
  {"log10e", M_LOG10E},
  {"sqrt2", M_SQRT2},
  {"sqrt2r", M_SQRT1_2},
  {"sqrt3", M_SQRT3},
  {"pi", M_PI},
  {"pihalf", M_PI_2},
  {"piquarter", M_PI_4},
  {"sqrtpi", M_SQRTPI},
  {"rpi", M_1_PI},
  {"tworpi", M_2_PI},
  {"ln10", M_LN10},
  {"ln2", M_LN2},
  {"lnpi", M_LNPI},
  {"euler", M_EULER},
  {"deg_per_rad", 180.0 / M_PI},
  {"rad_per_deg", M_PI / 180.0},
  {0, 0},
};

#endif  /* HAVE_RLABPLUS_SO */

// time-frequency distributions:
Bltin rlab_ctftbx_tfd_bltin[] = {
  {BLTIN, "butter", ent_ctftbx_tfd_butter},
  {BLTIN, "bj", ent_ctftbx_tfd_bj},
  {BLTIN, "cw", ent_ctftbx_tfd_cw},
  {BLTIN, "gr", ent_ctftbx_tfd_grd},
  {BLTIN, "mh", ent_ctftbx_tfd_mh},
  {BLTIN, "page", ent_ctftbx_tfd_page},
  {BLTIN, "pmh", ent_ctftbx_tfd_pmh},
  {BLTIN, "ppage", ent_ctftbx_tfd_ppage},
  {BLTIN, "pwv", ent_ctftbx_tfd_pwv},
  {BLTIN, "ri", ent_ctftbx_tfd_ri},
  {BLTIN, "spwv", ent_ctftbx_tfd_spwv},
  {BLTIN, "wv", ent_ctftbx_tfd_wv},
  {BLTIN, "zam", ent_ctftbx_tfd_zam},
  {0, 0, 0},
};
// time-frequency reduced interference distributions:
Bltin rlab_ctftbx_tfrid_bltin[] = {
  {BLTIN, "bessel", ent_ctftbx_tfrid_bessel},
  {BLTIN, "binom", ent_ctftbx_tfrid_binom},
  {BLTIN, "tri", ent_ctftbx_tfrid_tri},
  {BLTIN, "hann", ent_ctftbx_tfrid_hann},
  {0, 0, 0},
};

//
//  rlabplus extensions: python
//
#ifdef HAVE_PYTHON
#include "rlabplus_python.h"
#endif

//
// rlabplus extensions: java
//
#ifdef HAVE_JVM
#include "rjava.h"
#endif

//
// rlabplus extensions: ngspice
//
#ifdef HAVE_NGSPICE
#include "rngspice.h"
#endif

//
// rlabplus extensions: gphoto2
//
#ifdef HAVE_GPHOTO2
#include "rgphoto2.h"
#endif

//
// rlabplus extensions: imagemagick
//
#ifdef HAVE_IMAGEMAGICK  /* Set in config.h */
#include "r_imagemagick.h"
#endif  /* HAVE_IMAGEMAGICK */

/* **************************************************************
 * Install built-ins in symbol table.
 * Also, setup some global variables that define RLaB's configuration.
 * ************************************************************** */
void
init_symbol_table ()
{
  int i;
  Ent *ent;
  Btree *config_list=0;
  Btree *sp_config_list=0;
  Btree *pdf_config_list=0, *cpf_config_list=0, *odesl_config_list=0;
  Btree *claw_config_list=0;
  Btree *loess_config_list=0;
  Btree *ctftbx_tfd_config_list=0, *ctftbx_tfrid_config_list=0;
  MDS *config_string=0;

  symbol_table_create ();	/* Create RLaB global symbol table */

  //
  // Add sys list of functions for providing system wide parameters for different solvers
  //
  // sys = <<>>
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  Btree * sys_config_list = btree_Create ();
  ent_data (ent) = sys_config_list;
  install (0, "sys", ent);
  // sys.op = <<>>
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  Btree * sys_op_config_list = btree_Create ();
  ent_data (ent) = sys_op_config_list;
  install (sys_config_list, "op", ent);
  // sys.op.zero = <<>>
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  Btree * sys_op_zero_config_list = btree_Create ();
  ent_data (ent) = sys_op_zero_config_list;
  install (sys_op_config_list, "zero", ent);
  // sys.op.zero = rel(), abs()
  for (i = 0; rlab_sys_op_zero[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_sys_op_zero[i];
    install (sys_op_zero_config_list, rlab_sys_op_zero[i].name, ent);
  }
  for (i = 0; rlab_sys_op[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_sys_op[i];
    install (sys_op_config_list, rlab_sys_op[i].name, ent);
  }

#ifdef HAVE_PYTHON
  //
  // Add python list of functions
  //
  Btree *python_config_list=0;

  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  python_config_list = btree_Create ();
  python_config_list->isconst = 1;
  ent_data (ent) = python_config_list;
  install (0, (RLAB_NAME_EMBED_PYTHON), ent);
  for (i = 0; rlab_python_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_python_bltin[i];
    install (python_config_list, rlab_python_bltin[i].name, ent);
  }
#endif


#ifdef HAVE_JVM
  //
  // Add java list of functions
  //
  Btree *java_config_list=0;

  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  java_config_list = btree_Create ();
  java_config_list->isconst = 1;
  ent_data (ent) = java_config_list;
  install (0, (RLAB_NAME_EMBED_JVM), ent);
  for (i = 0; rlab_java_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_java_bltin[i];
    install (java_config_list, rlab_java_bltin[i].name, ent);
  }
#endif


#ifdef HAVE_NGSPICE

  //
  // Add ngspice list of functions
  //
  Btree *ngspice_config_list=0;

  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  ngspice_config_list = btree_Create ();
  ngspice_config_list->isconst = 1;
  ent_data (ent) = ngspice_config_list;
  install (0, (RLAB_NAME_EMBED_SPICE), ent);
  for (i = 0; rlab_ngspice_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_ngspice_bltin[i];
    install (ngspice_config_list, rlab_ngspice_bltin[i].name, ent);
  }
#endif


#ifdef HAVE_GPHOTO2
  //
  // Add gphoto2 list of functions
  //
  Btree *gphoto2_config_list=0;

  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  gphoto2_config_list = btree_Create ();
  gphoto2_config_list->isconst = 1;
  ent_data (ent) = gphoto2_config_list;
  install (0, (RLAB_NAME_EMBED_GPHOTO2), ent);
  for (i = 0; rlab_gphoto2_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_gphoto2_bltin[i];
    install (gphoto2_config_list, rlab_gphoto2_bltin[i].name, ent);
  }
#endif


#ifdef HAVE_IMAGEMAGICK
  //
  // Add imagemagick list of functions
  //
  Btree *im_config_list=0;

  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  im_config_list = btree_Create ();
  im_config_list->isconst = 1;
  ent_data (ent) = im_config_list;
  install (0, (RLAB_NAME_EMBED_IM), ent);
  for (i = 0; rlab_im_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_im_bltin[i];
    install (im_config_list, rlab_im_bltin[i].name, ent);
  }
#endif

  //
  // Add loess list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  loess_config_list = btree_Create ();
  loess_config_list->isconst = 1;
  ent_data (ent) = loess_config_list;
  install (0, ("loess"), ent);
  for (i = 0; rlab_loess_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_loess_bltin[i];
    install (loess_config_list, rlab_loess_bltin[i].name, ent);
  }
  //
  // Add  clawparams  list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  claw_config_list = btree_Create ();
  claw_config_list->isconst = 1;
  ent_data (ent) = claw_config_list;
  install (0, ("clawparams"), ent);
  for (i = 0; rlab_clawparams_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_clawparams_bltin[i];
    install (claw_config_list, rlab_clawparams_bltin[i].name, ent);
  }
  //
  // Add pdf list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  pdf_config_list = btree_Create ();
  pdf_config_list->isconst = 1;
  ent_data (ent) = pdf_config_list;
  install (0, RLAB_PDF_BLTIN, ent);
  for (i = 0; rlab_pdf_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_pdf_bltin[i];
    install (pdf_config_list, rlab_pdf_bltin[i].name, ent);
  }
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  cpf_config_list = btree_Create ();
  cpf_config_list->isconst = 1;
  ent_data (ent) = cpf_config_list;
  install (0, RLAB_CDF_BLTIN, ent);
  for (i = 0; rlab_cpf_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_cpf_bltin[i];
    install (cpf_config_list, rlab_cpf_bltin[i].name, ent);
  }
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  Btree * cpfinv_config_list = btree_Create ();
  cpfinv_config_list->isconst = 1;
  ent_data (ent) = cpfinv_config_list;
  install (0, RLAB_CDFINV_BLTIN, ent);
  for (i = 0; rlab_cpfinv_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_cpfinv_bltin[i];
    install (cpfinv_config_list, rlab_cpfinv_bltin[i].name, ent);
  }

  Btree *cluster_config_list;
  //
  // Add cluster list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  cluster_config_list = btree_Create ();
  cluster_config_list->isconst = 1;
  ent_data (ent) = cluster_config_list;
  install (0, "cluster", ent);
  for (i = 0; rlab_cluster_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_cluster_bltin[i];
    install (cluster_config_list, rlab_cluster_bltin[i].name, ent);
  }

  Btree *classify_config_list;
  //
  // Add classify list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  classify_config_list = btree_Create ();
  classify_config_list->isconst = 1;
  ent_data (ent) = classify_config_list;
  install (0, cpstr ("classify"), ent);
  for (i = 0; rlab_classify_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_classify_bltin[i];
    install (classify_config_list, rlab_classify_bltin[i].name, ent);
  }

//   Btree *ffnn_config_list;
//   //
//   // Add ffnn list of functions
//   //
//   ent = ent_Create ();
//   ent_SetType (ent, BTREE);
//   ffnn_config_list = btree_Create ();
//   ffnn_config_list->isconst = 1;
//   ent_data (ent) = ffnn_config_list;
//   install (0, cpstr ("ffnn"), ent);
//   for (i = 0; rlab_ffnn_bltin[i].name; i++)
//   {
//     ent = ent_Create ();
//     ent_SetType (ent, BLTIN);
//     ent_data (ent) = &rlab_ffnn_bltin[i];
//     install (ffnn_config_list, rlab_ffnn_bltin[i].name, ent);
//   }


#ifdef HAVE_RLABPLUS_SO

  Btree *math_const_list, *mks_const_list;

  //
  // Add mks units to RLaB: see bltin1.c, code.c for how to protect any other
  // variable from deletion or changing content
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  mks_const_list = btree_Create ();
  mks_const_list->isconst = 1;
  ent_data (ent) = mks_const_list;
  install (0, ("mks"), ent);
  for (i = 0; gsl_mks_double[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, MATRIX_DENSE_REAL);
    ent_data (ent) = mdr_CreateScalar (gsl_mks_double[i].val);
    install (mks_const_list, gsl_mks_double[i].name, ent);
  }
  //
  // Add math constants to RLaB
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  math_const_list = btree_Create ();
  math_const_list->isconst = 1;
  ent_data (ent) = math_const_list;
  install (0, ("const"), ent);
  for (i = 0; gsl_math_double[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, MATRIX_DENSE_REAL);
    ent_data (ent) = mdr_CreateScalar (gsl_math_double[i].val);
    install (math_const_list, gsl_math_double[i].name, ent);
  }

#endif /* HAVE_RLABPLUS_SO */

  //
  // Add curl list of functions
  //
  // curl = <<>>
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  Btree * curl_config_list = btree_Create ();
  ent_data (ent) = curl_config_list;
  install (0, "curl", ent);
  for (i = 0; rlab_curl_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_curl_bltin[i];
    install (curl_config_list, rlab_curl_bltin[i].name, ent);
  }
  // curl.get = <<>>
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  Btree * curl_get_config_list = btree_Create ();
  ent_data (ent) = curl_get_config_list;
  install (curl_config_list, "get", ent);
  for (i = 0; rlab_curl_get_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_curl_get_bltin[i];
    install (curl_get_config_list, rlab_curl_get_bltin[i].name, ent);
  }

  //
  // Add sep list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  Btree * sep_config_list = btree_Create ();
  sep_config_list->isconst = 1;
  ent_data (ent) = sep_config_list;
  install (0, "sep", ent);
  for (i = 0; rlab_sep_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_sep_bltin[i];
    install (sep_config_list, rlab_sep_bltin[i].name, ent);
  }

  //
  // Add odesl() list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  odesl_config_list = btree_Create ();
  odesl_config_list->isconst = 1;
  ent_data (ent) = odesl_config_list;
  install (0, ("odesl"), ent);
  for (i = 0; rlab_odesl_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_odesl_bltin[i];
    install (odesl_config_list, rlab_odesl_bltin[i].name, ent);
  }

  //
  // Add sparams() list of functions
  //
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  sp_config_list = btree_Create ();
  sp_config_list->isconst = 1;
  ent_data (ent) = sp_config_list;
  install (0, ("sparams"), ent);
  for (i = 0; rlab_sparams_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_sparams_bltin[i];
    install (sp_config_list, rlab_sparams_bltin[i].name, ent);
  }

  //
  // time-frequency toolbox
  //
  //Add tfd() list of functions
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  ctftbx_tfd_config_list = btree_Create ();
  ctftbx_tfd_config_list->isconst = 1;
  ent_data (ent) = ctftbx_tfd_config_list;
  install (0, ("tfrd"), ent);
  for (i = 0; rlab_ctftbx_tfd_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_ctftbx_tfd_bltin[i];
    install (ctftbx_tfd_config_list, rlab_ctftbx_tfd_bltin[i].name, ent);
  }
  //Add tfd() list of functions
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  ctftbx_tfrid_config_list = btree_Create ();
  ctftbx_tfrid_config_list->isconst = 1;
  ent_data (ent) = ctftbx_tfrid_config_list;
  install (0, ("tfrid"), ent);
  for (i = 0; rlab_ctftbx_tfrid_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_ctftbx_tfrid_bltin[i];
    install (ctftbx_tfrid_config_list, rlab_ctftbx_tfrid_bltin[i].name, ent);
  }

  //
  // load all RLaB built-in functions
  //
  for (i = 0; rlab_bltin[i].name; i++)
  {
    ent = ent_Create ();
    ent_SetType (ent, BLTIN);
    ent_data (ent) = &rlab_bltin[i];
    install (0, rlab_bltin[i].name, ent);
  }

  //
  // Define, within rlab's global workspace, rlab's configuration.
  // Obviously, these are configuration details that are determined
  // at compile time.
  //
  // Create the var: RLAB_SYSVAR_CONFIG = "_rlab_config"
  ent = ent_Create ();
  ent_SetType (ent, BTREE);
  config_list = btree_Create ();
  config_list->isconst = 1;
  ent_data (ent) = config_list;
  install (0, RLAB_SYSVAR_CONFIG, ent);

  // Install version of rlab
  install (config_list, RLAB_SYSVAR_CONFIG_VERSION, ent_Create_Rlab_String(version_string));
  // Install plot support
  install (config_list, RLAB_SYSVAR_CONFIG_PLOT, ent_Create_Rlab_String(RLAB_PLOT_SUPPORT));

  // Install the various configuration information
  ent = ent_Create ();
  ent_SetType (ent, MATRIX_DENSE_STRING);
#if defined (HAVE_GC)
  config_string = mds_CreateScalar ( "gc" );
#else
  config_string = mds_CreateScalar ( "none" );
#endif
  ent_data (ent) = config_string;
  install (config_list, ("garbage_collector"), ent);

  ent = ent_Create ();
  ent_SetType (ent, MATRIX_DENSE_STRING);
#if defined (HAVE_METIS3)
  config_string = mds_CreateScalar ( "metis" );
#else
  config_string = mds_CreateScalar ( "none" );
#endif
  ent_data (ent) = config_string;
  install (config_list, ("graph_support"), ent);


  ent = ent_Create ();
  ent_SetType (ent, MATRIX_DENSE_STRING);
#if defined (HAVE_READLINE)
  config_string = mds_CreateScalar ( "readline" );
#else
  config_string = mds_CreateScalar ( "none" );
#endif
  ent_data (ent) = config_string;
  install (config_list, ("command_line_support"), ent);

  ent = ent_Create ();
  ent_SetType (ent, MATRIX_DENSE_STRING);
#if defined (WE_ARE_BIG_ENDIAN)
  config_string = mds_CreateScalar ( "big" );
#else
  config_string = mds_CreateScalar ( "little" );
#endif
  ent_data (ent) = config_string;
  install (config_list, ("endianess"), ent);
}
