// rlabplus extensions: GNU scientific library
// gsl_odeiv.c
extern Ent *  ent_gsl_odeiv(int nargs, Datum args[]);
// gsl_levin.c
extern Ent *  ent_gsl_levinu(int nargs, Datum args[]);
// gsl_spline.c
extern Ent *  ent_cubicspline_init(int nargs, Datum args[]);
extern Ent *  ent_cubicspline_eval(int nargs, Datum args[]);
extern Ent *  ent_lininterp_eval(int nargs, Datum args[]);
extern Ent *  ent_lininterp1_eval(int nargs, Datum args[]);
extern Ent *  ent_harmsum_eval(int nargs, Datum args[]);
// gsl_ndiff_md.c
extern Ent *  ent_ndiffs(int nargs, Datum args[]);
extern Ent *  ent_div(int nargs, Datum args[]);
// gsl_leastsquares.c
extern Ent *  ent_gsl_lsfit(int nargs, Datum args[]);
extern Ent *  ent_polyfit(int nargs, Datum args[]);
// gsl_minimize.c
extern Ent *  ent_minimize_1d_bracket(int nargs, Datum args[]);
extern Ent *  ent_gsl_findmins(int nargs, Datum args[]);
extern Ent *  ent_conmax(int nargs, Datum args[]);
// gsl_findroot.c
extern Ent *  ent_findroot(int nargs, Datum args[]);
extern Ent *  ent_gsl_findroots(int nargs, Datum args[]);
extern Ent *  ent_ctrack(int nargs, Datum args[]);
// gsl_diffint.c
extern Ent *  ent_tchebyshev_init(int nargs, Datum args[]);
extern Ent *  ent_tchebyshev_eval(int nargs, Datum args[]);
extern Ent *  ent_tchebyshev_diff(int nargs, Datum args[]);
extern Ent *  ent_tchebyshev_int(int nargs, Datum args[]);
extern Ent *  ent_ndiff(int nargs, Datum args[]);
// gsl_stat.c
extern Ent *  ent_covariance(int nargs, Datum args[]);
extern Ent *  ent_mean(int nargs, Datum args[]);
extern Ent *  ent_variance(int nargs, Datum args[]);
extern Ent *  ent_absdev(int nargs, Datum args[]);
extern Ent *  ent_skew(int nargs, Datum args[]);
extern Ent *  ent_kurtosis(int nargs, Datum args[]);
extern Ent *  ent_wmean(int nargs, Datum args[]);
extern Ent *  ent_wvariance(int nargs, Datum args[]);
extern Ent *  ent_wvariance_fm(int nargs, Datum args[]);
extern Ent *  ent_wabsdev(int nargs, Datum args[]);
extern Ent *  ent_wskew(int nargs, Datum args[]);
extern Ent *  ent_wkurtosis(int nargs, Datum args[]);
extern Ent *  ent_median(int nargs, Datum args[]);
extern Ent *  ent_quantile(int nargs, Datum args[]);
extern Ent *  ent_invperm(int nargs, Datum args[]);
extern Ent *  ent_revperm(int nargs, Datum args[]);
extern Ent *  ent_nextperm(int nargs, Datum args[]);
extern Ent *  ent_prevperm(int nargs, Datum args[]);
extern Ent *  ent_isaperm(int nargs, Datum args[]);
extern Ent *  ent_mulperm(int nargs, Datum args[]);
extern Ent *  ent_prevcomb(int nargs, Datum args[]);
extern Ent *  ent_nextcomb(int nargs, Datum args[]);
extern Ent *  ent_iscomb(int nargs, Datum args[]);
extern Ent *  ent_ncomb(int nargs, Datum args[]);
extern Ent *  ent_histogram_create(int nargs, Datum args[]);
extern Ent *  ent_histogram2d_create(int nargs, Datum args[]);

// gsl_nintegrate.c
extern Ent *  ent_nintegrate(int nargs, Datum args[]);
extern Ent *  ent_nintegrate_qagp(int nargs, Datum args[]);
extern Ent *  ent_nintegrate_qawc(int nargs, Datum args[]);
extern Ent *  ent_nintegrate_qaws(int nargs, Datum args[]);
extern Ent *  ent_nintegrate_qawo_sin(int nargs, Datum args[]);
extern Ent *  ent_nintegrate_qawo_cos(int nargs, Datum args[]);
extern Ent *  ent_convolution        (int nargs, Datum args[]);
extern Ent *  ent_laplacetransform   (int nargs, Datum args[]);
extern Ent *  ent_nintegrate_genz_hypercube (int nargs, Datum args[]);
extern Ent *  ent_nintegrate_genz_simplex (int nargs, Datum args[]);
extern Ent *  ent_invlt (int nargs, Datum args[]);

// gsl_specfunc.c
extern Ent *  ent_gsl_sf_airy_ai (int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_airy_bi(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_airy_dai (int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_airy_dbi(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_airy_zeroai (int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_airy_zerobi(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_airy_zerodai (int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_airy_zerodbi (int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_I(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_J(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_K(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_Y(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_i(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_j(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_k(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_y(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_bessel_zeroJ(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_clausen(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_hydrogen(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_dawson(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_debye(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_E(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_F(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_K(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_P(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_D(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_RC(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_RD(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_RF(int nargs, Datum args[]);
extern Ent * ent_gsl_sf_elliptic_RJ(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_elliptic_jacobi(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_erf(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_erfc(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_erfq(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_erfz(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_expint(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_expint_ei(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_shi(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_chi(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_si(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_ci(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_atanint(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_fermidiracint(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_loggamma(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_gamma(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_gammareg(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_gammaregc(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_recgamma(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_logbeta(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_beta(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_betareg(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_poch(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_gegen(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_hyperg_0F1(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_hyperg_1F1(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_hyperg_U(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_hyperg_2F0(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_hyperg_2F1(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_laguerre(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_lambert_W0(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_legendreP(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_legendreQ(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_legendre_sph_Plm(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_legendre_con_Plm(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_legendre_h3d(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_psi(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_polygamma(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_synchrotron(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_transport(int nargs, Datum args[]);
extern Ent *  ent_gsl_sf_zeta(int nargs, Datum args[]);
extern Ent *  ent_sf_erling(int nargs, Datum args[]);
// gsl_rng.c, gsl_rng_nintegrate_mc.c, gsl_rng_siman.c
extern Ent *  ent_gsl_irng       (int nargs, Datum args[]);
extern Ent *  ent_gsl_irand      (int nargs, Datum args[]);
extern Ent *  ent_gsl_irng_state (int nargs, Datum args[]);
extern Ent *  ent_gsl_rng       (int nargs, Datum args[]);
extern Ent *  ent_gsl_rand      (int nargs, Datum args[]);
extern Ent *  ent_gsl_urand      (int nargs, Datum args[]);
extern Ent *  ent_gsl_gauss      (int nargs, Datum args[]);
extern Ent *  ent_gsl_randomize (int nargs, Datum args[]);
extern Ent *  ent_gsl_rand2(int nargs, Datum args[]);
extern Ent *  ent_gsl_shuffle(int nargs, Datum args[]);
extern Ent *  ent_gsl_sample(int nargs, Datum args[]);
extern Ent *  ent_gsl_setdrng(int nargs, Datum args[]);
extern Ent *  ent_gsl_drand  (int nargs, Datum args[]);
extern Ent *  ent_gsl_sethrng(int nargs, Datum args[]);
extern Ent *  ent_gsl_hrand  (int nargs, Datum args[]);
extern Ent *  ent_nintegrate_mc(int nargs, Datum args[]);
extern Ent *  ent_gsl_siman(int nargs, Datum args[]);
extern Ent *  ent_gsl_diffevol(int nargs, Datum args[]);
// gsl_pdf.c functions
//  - non gsl
extern Ent *  ent_misc_pdf_driftdiffusion (int nargs, Datum args[]);
//  - gsl
extern Ent *  ent_gsl_pdf_poisson (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_geometric (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_log (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_binomial (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_negbinomial (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_hypergeom (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_multinom (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_dirichlet (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_normal (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_exp (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_laplace (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_cauchy (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_rayleigh (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_chisq (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_t (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_logistic (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_normaltail (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_exppow (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_rayleightail(int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_gamma (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_lognormal (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_F (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_beta (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_pareto (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_weibull (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_uniform (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_gumbel1 (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_gumbel2 (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_normalbiv (int nargs, Datum args[]);

extern Ent *  ent_gsl_cpf_poisson (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_geometric (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_binomial (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_negbinomial (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_hypergeom (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_normal (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_exp (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_laplace (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_cauchy (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_rayleigh (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_chisq (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_t (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_logistic (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_exppow (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_gamma (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_lognormal (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_F (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_beta (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_pareto (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_weibull (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_uniform (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_gumbel1 (int nargs, Datum args[]);
extern Ent *  ent_gsl_cpf_gumbel2 (int nargs, Datum args[]);


// gsl_wavelet.c functions
extern Ent *  ent_gsl_dwt (int nargs, Datum args[]);
extern Ent *  ent_gsl_idwt (int nargs, Datum args[]);
extern Ent *  ent_gsl_dwt2 (int nargs, Datum args[]);
extern Ent *  ent_gsl_idwt2 (int nargs, Datum args[]);

// gnunet stuff
extern Ent *  ent_cluster_isodata (int nargs, Datum args[]);
extern Ent *  ent_cluster_pao     (int nargs, Datum args[]);
extern Ent *  ent_cluster_nn      (int nargs, Datum args[]);
extern Ent *  ent_rank            (int nargs, Datum args[]);
extern Ent *  ent_distance_matrix (int nargs, Datum args[]);
extern Ent *  ent_classify_parzen (int nargs, Datum args[]);
extern Ent *  ent_sprannlib_classify_fisher (int nargs, Datum args[]);
extern Ent *  ent_sprannlib_classify_fisherq(int nargs, Datum args[]);
extern Ent *  ent_classify_knn    (int nargs, Datum args[]);
extern Ent *  ent_sprannlib_sammon (int nargs, Datum args[]);
extern Ent *  ent_ffnn_samann (int nargs, Datum args[]);
extern Ent *  ent_ffnn_eval (int nargs, Datum args[]);
extern Ent *  ent_ffnn (int nargs, Datum args[]);



