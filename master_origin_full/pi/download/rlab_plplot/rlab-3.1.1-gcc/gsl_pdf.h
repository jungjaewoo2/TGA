//
// pdf(), list of probability distribution functions (GSL)
//
// 2 parameters
extern Ent *  ent_gsl_ran_logistic_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_tdist_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_chisq_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_rayleigh_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_cauchy_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_exponential_pdf(int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_poisson_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_geometric_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_logarithmic_pdf(int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_laplace_pdf (int nargs, Datum args[]);
// 3 parameters
extern Ent *  ent_gsl_ran_flat_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_gumbel1_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_gumbel2_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_weibull_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_pareto_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_beta_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_fdist_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_lognormal_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_gamma_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_rayleigh_tail_pdf(int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_exppow_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_gaussian_tail_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_binomial_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_negative_binomial_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_pascal_pdf (int nargs, Datum args[]);
// 4 parameters
extern Ent *  ent_gsl_ran_hypergeometric_pdf(int nargs, Datum args[]);
// arrays of parameters
extern Ent *  ent_gsl_ran_dirichlet_pdf (int nargs, Datum args[]);
extern Ent *  ent_gsl_ran_multinomial_pdf (int nargs, Datum args[]);

extern Ent *  ent_rlab2gsl_ran_gaussian_pdf (int nargs, Datum args[]);
extern Ent *  ent_rlab2gsl_pdf_normal (int nargs, Datum args[]);
extern Ent *  ent_gsl_pdf_normalbiv (int nargs, Datum args[]);

//  - non gsl
extern Ent *  ent_misc_pdf_driftdiffusion (int nargs, Datum args[]);

#define RLAB_PDF_BLTIN "pdf"
Bltin rlab_pdf_bltin[] = {
  // 2-parameter distributions
  {BLTIN, "logistic", ent_gsl_ran_logistic_pdf},
  {BLTIN, "t", ent_gsl_ran_tdist_pdf},
  {BLTIN, "chisq", ent_gsl_ran_chisq_pdf},
  {BLTIN, "rayleigh", ent_gsl_ran_rayleigh_pdf},
  {BLTIN, "cauchy", ent_gsl_ran_cauchy_pdf},
  {BLTIN, "exp", ent_gsl_ran_exponential_pdf},
  {BLTIN, "poisson", ent_gsl_ran_poisson_pdf},
  {BLTIN, "geometric", ent_gsl_ran_geometric_pdf},
  {BLTIN, "log", ent_gsl_ran_logarithmic_pdf},
  {BLTIN, "laplace", ent_gsl_ran_laplace_pdf},
  // 3-parameter distributions
  {BLTIN, "uniform", ent_gsl_ran_flat_pdf},
  {BLTIN, "gumbel1", ent_gsl_ran_gumbel1_pdf},
  {BLTIN, "gumbel2", ent_gsl_ran_gumbel2_pdf},
  {BLTIN, "weibull", ent_gsl_ran_weibull_pdf},
  {BLTIN, "pareto", ent_gsl_ran_pareto_pdf},
  {BLTIN, "beta", ent_gsl_ran_beta_pdf},
  {BLTIN, "F", ent_gsl_ran_fdist_pdf},
  {BLTIN, "lognormal", ent_gsl_ran_lognormal_pdf},
  {BLTIN, "gamma", ent_gsl_ran_gamma_pdf},
  {BLTIN, "rayleightail", ent_gsl_ran_rayleigh_tail_pdf},
  {BLTIN, "exppow", ent_gsl_ran_exppow_pdf},
  {BLTIN, "normaltail", ent_gsl_ran_gaussian_tail_pdf},
  {BLTIN, "binomial", ent_gsl_ran_binomial_pdf},
  {BLTIN, "negbinomial", ent_gsl_ran_negative_binomial_pdf},
  {BLTIN, "pascal", ent_gsl_ran_pascal_pdf},
  // 4-parameters
  {BLTIN, "hypergeom", ent_gsl_ran_hypergeometric_pdf},
  //
  {BLTIN, "dirichlet", ent_gsl_ran_dirichlet_pdf},
  {BLTIN, "multinom", ent_gsl_ran_multinomial_pdf},
  //
  {BLTIN, "normal", ent_rlab2gsl_pdf_normal},
  {BLTIN, "gaussian", ent_rlab2gsl_pdf_normal},
  {BLTIN, "normalbiv", ent_gsl_pdf_normalbiv},
  {BLTIN, "dd", ent_misc_pdf_driftdiffusion},
  {0, 0, 0}
};


extern Ent *  ent_gsl_cdf_logistic_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_tdist_P (int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_chisq_P (int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_poisson_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_geometric_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_binomial_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_negative_binomial_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_hypergeometric_P(int nargs, Datum args[]);
//extern Ent *  ent_gsl_cpf_normal (int nargs, Datum args[]);
extern Ent *  ent_rlab2gsl_cdf_normal_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_exponential_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_laplace_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_cauchy_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_rayleigh_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_exppow_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_gamma_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_lognormal_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_fdist_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_beta_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_pareto_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_weibull_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_flat_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_gumbel1_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_gumbel2_P(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_pascal_P(int nargs, Datum args[]);

#define RLAB_CDF_BLTIN "cdf"
Bltin rlab_cpf_bltin[] = {
  {BLTIN, "logistic", ent_gsl_cdf_logistic_P},
  {BLTIN, "t", ent_gsl_cdf_tdist_P},
  {BLTIN, "chisq", ent_gsl_cdf_chisq_P},
  {BLTIN, "poisson", ent_gsl_cdf_poisson_P},
  {BLTIN, "geometric", ent_gsl_cdf_geometric_P},
  {BLTIN, "binomial", ent_gsl_cdf_binomial_P},
  {BLTIN, "negbinomial", ent_gsl_cdf_negative_binomial_P},
  {BLTIN, "hypergeom", ent_gsl_cdf_hypergeometric_P},
  {BLTIN, "normal", ent_rlab2gsl_cdf_normal_P},
  {BLTIN, "gaussian", ent_rlab2gsl_cdf_normal_P},
//   {BLTIN, "normal", ent_gsl_cpf_normal},
//   {BLTIN, "gaussian", ent_gsl_cpf_normal},
  {BLTIN, "exp", ent_gsl_cdf_exponential_P},
  {BLTIN, "laplace", ent_gsl_cdf_laplace_P},
  {BLTIN, "cauchy", ent_gsl_cdf_cauchy_P},
  {BLTIN, "rayleigh", ent_gsl_cdf_rayleigh_P},
  {BLTIN, "exppow", ent_gsl_cdf_exppow_P},
  {BLTIN, "gamma", ent_gsl_cdf_gamma_P},
  {BLTIN, "lognormal", ent_gsl_cdf_lognormal_P},
  {BLTIN, "F", ent_gsl_cdf_fdist_P},
  {BLTIN, "beta", ent_gsl_cdf_beta_P},
  {BLTIN, "pareto", ent_gsl_cdf_pareto_P},
  {BLTIN, "weibull", ent_gsl_cdf_weibull_P},
  {BLTIN, "uniform", ent_gsl_cdf_flat_P},
  {BLTIN, "gumbel1", ent_gsl_cdf_gumbel1_P},
  {BLTIN, "gumbel2", ent_gsl_cdf_gumbel2_P},
  {BLTIN, "pascal", ent_gsl_cdf_pascal_P},
  {0, 0, 0}
};


extern Ent *  ent_gsl_cdf_gamma_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_rayleigh_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_cauchy_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_exponential_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_laplace_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_flat_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_lognormal_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_chisq_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_fdist_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_tdist_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_beta_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_logistic_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_pareto_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_weibull_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_gumbel1_Pinv(int nargs, Datum args[]);
extern Ent *  ent_gsl_cdf_gumbel2_Pinv(int nargs, Datum args[]);
extern Ent *  ent_rlab2gsl_cdf_normal_invP(int nargs, Datum args[]);
#define RLAB_CDFINV_BLTIN "invcdf"
Bltin rlab_cpfinv_bltin[] = {
  {BLTIN, "normal", ent_rlab2gsl_cdf_normal_invP},
  {BLTIN, "gaussian", ent_rlab2gsl_cdf_normal_invP},
  {BLTIN, "exp", ent_gsl_cdf_exponential_Pinv},
  {BLTIN, "laplace", ent_gsl_cdf_laplace_Pinv},
  {BLTIN, "cauchy", ent_gsl_cdf_cauchy_Pinv},
  {BLTIN, "rayleigh", ent_gsl_cdf_rayleigh_Pinv},
  {BLTIN, "rayleigh", ent_gsl_cdf_gamma_Pinv},
  {BLTIN, "uniform", ent_gsl_cdf_flat_Pinv},
  {BLTIN, "lognormal", ent_gsl_cdf_lognormal_Pinv},
  {BLTIN, "chisq", ent_gsl_cdf_chisq_Pinv},
  {BLTIN, "F", ent_gsl_cdf_fdist_Pinv},
  {BLTIN, "t", ent_gsl_cdf_tdist_Pinv},
  {BLTIN, "beta", ent_gsl_cdf_beta_Pinv},
  {BLTIN, "logistic", ent_gsl_cdf_logistic_Pinv},
  {BLTIN, "pareto", ent_gsl_cdf_pareto_Pinv},
  {BLTIN, "weibull", ent_gsl_cdf_weibull_Pinv},
  {BLTIN, "gumbel1", ent_gsl_cdf_gumbel1_Pinv},
  {BLTIN, "gumbel2", ent_gsl_cdf_gumbel2_Pinv},
  {0,0,0}
};













