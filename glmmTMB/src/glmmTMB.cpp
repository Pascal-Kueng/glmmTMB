#define EIGEN_DONT_PARALLELIZE // see https://github.com/kaskr/adcomp/issues/390
#define DISABLE_AD if(isDouble<Type>::value)
#define NAN std::numeric_limits<double>::quiet_NaN()
// set correlation matrix for output, unless corr output has been disabled
// (make this a function instead?
#define SET_COR    if (term.fullCor == 1) { \
      term.corr = nldens.cov(); \
    } else { \
      term.corr.resize(1,1); \
      term.corr(0,0) = NAN; \
    }

#include <TMB.hpp>
#include <R_ext/Error.h>
#include "init.h"
#include "distrib.h"
#include "cordistrib.h"

// don't need to include omp.h; we get it via TMB.hpp

namespace glmmtmb{
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
bool notFinite(Type x) {
  return (!R_FINITE(asDouble(x)));
}
}

enum valid_family {
  gaussian_family = 0,
  binomial_family = 100,
  betabinomial_family =101,
  beta_family =200,
  ordbeta_family = 201,
  Gamma_family =300,
  poisson_family =400,
  truncated_poisson_family =401,
  genpois_family =402,
  compois_family =403,
  truncated_genpois_family =404,
  truncated_compois_family =405,
  nbinom1_family =500,
  nbinom2_family =501,
  nbinom12_family =502,
  truncated_nbinom1_family =550,
  truncated_nbinom2_family =551,
  t_family =600,
  tweedie_family = 700,
  lognormal_family = 800,
  skewnormal_family = 900,
  bell_family = 1000
};

// capitalize Family so this doesn't get picked up by the 'enum' scraper
bool trunc_Family(int family) {
  return (family == truncated_poisson_family ||
	  family == truncated_genpois_family ||
	  family == truncated_compois_family ||
	  family == truncated_nbinom1_family ||
	  family == truncated_nbinom2_family);
  
}

enum valid_link {
  log_link                 = 0,
  logit_link               = 1,
  probit_link              = 2,
  inverse_link             = 3,
  cloglog_link             = 4,
  identity_link            = 5,
  sqrt_link                = 6,
  lambertW_link            = 7
};

enum valid_covStruct {
  diag_covstruct = 0,
  us_covstruct   = 1,
  cs_covstruct   = 2,
  ar1_covstruct  = 3,
  ou_covstruct   = 4,
  exp_covstruct = 5,
  gau_covstruct = 6,
  mat_covstruct = 7,
  toep_covstruct = 8,
  rr_covstruct = 9,
  homdiag_covstruct = 10,
  propto_covstruct = 11,
  // should perhaps be next to homdiag but don't want to mess
  //  up interpretation of stored fits ...
  hetar1_covstruct = 12,
  homcs_covstruct = 13,
  homtoep_covstruct = 14,
  equalto_covstruct = 15,
  separable_covstruct = 16
};

enum separable_density_kind {
  // These values must match `.sep_density_kind_code` in R/utils_covstruct.R.
  // They are deliberately separate from valid_covStruct: many covariance
  // structures can share the same separable density kind.  For example, homcs
  // and us are both "dense correlation" margins in the current prototype.
  dense_corr_sep = 1,
  ar1_sep = 2,
  diag_sep = 3,
  spatial_sep = 4,
  toep_sep = 5
};

enum separable_dispatch {
  // Evaluator codes for separable covariance structures.
  corr_corr_dispatch = 1
};

enum separable_scale_mode {
  // These values must match `.sep_scale_mode_code` in R/utils_covstruct.R.
  margin_sep_scale = 1,
  global_sep_scale = 2,
  product_sep_scale = 3,
  selected_product_sep_scale = 4
};

// should probably be named just 'predictCode';
// originally for enabling z-i prediction
// 'corrected' = mean prediction incorporates z-i effects
// 'uncorrected' = mean not accounting for z-i
// 'prob' = zero-inflation on back-transformed (probability) scale
// 'disp' = report value of dispersion parameter
enum valid_ziPredictCode {
  corrected_zipredictcode = 0,
  uncorrected_zipredictcode = 1,
  prob_zipredictcode = 2,
  disp_zipredictcode = 3
};


enum valid_simCode {
  zero_simcode = 0,
  fix_simcode = 1,
  random_simcode = 2
};
  
// codes for prior distributions
enum valid_prior {
  // real-valued
  normal_prior = 0,
  t_prior = 1,
  cauchy_prior = 2,
  // non-negative
  gamma_prior = 10,
  // (0,1), e.g. for zi prob
  beta_prior = 20,
  // correlations
  lkj_prior = 30
};

// codes for parameter (vec) to apply prior to
enum valid_vprior {
  beta_vprior = 0,
  betazi_vprior = 1,
  betadisp_vprior = 2,
  theta_vprior = 10,
  thetazi_vprior = 20,
  psi_vprior = 30
};

template<class Type>
Type inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = exp(eta);
    break;
  case identity_link:
    ans = eta;
    break;
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  case cloglog_link:
    ans = Type(1) - exp(-exp(eta));
    break;
  case inverse_link:
    ans = Type(1) / eta;
    break;
  case sqrt_link:
    ans = eta*eta; // pow(eta, Type(2)) doesn't work ... ?
    break;
  case lambertW_link:
    // for Bell distribution: mean = theta*exp(theta), theta 
    ans = exp(eta)*exp(exp(eta));
    break;

    // TODO: Implement remaining links
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

template<class Type>
Type linkfun(Type mu, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = log(mu);
    break;
  case identity_link:
    ans = mu;
    break;
  case logit_link:
    ans = logit(mu);
    break;
  case probit_link:
    ans = qnorm(mu);
    break;
  case cloglog_link:
    ans = log(-log(Type(1)-mu));
    break;
  case inverse_link:
    ans = Type(1) / mu;
    break;
  case sqrt_link:
    ans = sqrt(mu);
    break;
    // TODO: Implement remaining links
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

/* logit transformed inverse_linkfun without losing too much
   accuracy */
template<class Type>
Type logit_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case logit_link:
    ans = eta;
    break;
  case probit_link:
    ans = glmmtmb::logit_pnorm(eta);
    break;
  case cloglog_link:
    ans = glmmtmb::logit_invcloglog(eta);
    break;
  default:
    ans = logit( inverse_linkfun(eta, link) );
  } // End switch
  return ans;
}

/* log transformed inverse_linkfun without losing too much accuracy */
template<class Type>
Type log_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = eta;
    break;
  case logit_link:
    ans = -logspace_add(Type(0), -eta);
    break;
  default:
    ans = log( inverse_linkfun(eta, link) );
  } // End switch
  return ans;
}

/* log transformed (1-inverse_linkfun) without losing too much accuracy */
template<class Type>
Type log1m_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case log_link:
    ans = logspace_sub(Type(0), eta);
    break;
  case logit_link:
    ans = -logspace_add(Type(0), eta);
    break;
  default:
    ans = logspace_sub(Type(0), log( inverse_linkfun(eta, link) ));
  } // End switch
  return ans;
}

/* log-prob of non-zero value in conditional distribution  */
template<class Type>
Type calc_log_nzprob(Type mu, Type phi, Type eta, Type etadisp, int family,
		     int link) {
  Type ans, s1, s2;
  switch (family) {
  case truncated_nbinom1_family:
    s2 = logspace_add( Type(0), etadisp);      // log(1. + phi(i)
    ans = logspace_sub( Type(0), -mu / phi * s2 ); // 1-prob(0)
    break;
  case truncated_nbinom2_family:
    // s1 is repeated computation from main loop ...
    s1 = log_inverse_linkfun(eta, link);          // log(mu)
    // s2 := log( 1. + mu(i) / phi(i) )
    s2 = logspace_add( Type(0), s1 - etadisp );
    ans = logspace_sub( Type(0), -phi * s2 );
    break;
  case truncated_poisson_family:
    ans = logspace_sub(Type(0), -mu);  // log(1-exp(-mu(i))) = P(x>0)
    break;
  case truncated_genpois_family:
    s1 = mu / sqrt(phi); //theta
    s2 = Type(1) - Type(1)/sqrt(phi); //lambda
    ans = logspace_sub(Type(0), -s1);
    break;
  case truncated_compois_family:
    ans = logspace_sub(Type(0), dcompois2(Type(0), mu, 1/phi, true));
    break;
  default: ans = Type(0);
  }
  return ans;
}

template <class Type>
struct per_term_info {
  // Input from R
  int blockCode;     // Code that defines structure
  int blockSize;     // Size of one block
  int blockReps;     // Repeat block number of times
  int blockNumTheta; // Parameter count per block
  int simCode;       // Simulation code (zero, fixed, or draw new random deviate?)
  int fullCor;       // Compute/store full correlation matrix?
  matrix<Type> dist;
  vector<Type> times;// For ar1 case
  // Metadata for separable covariance structures.
  //
  // The random effects still arrive in the usual glmmTMB representation:
  // one flat vector per grouping level.  For a two-dimensional separable term,
  // the flat vector is conceptually an array with dimensions:
  //
  //   sepDims(0) = number of columns in the first product margin
  //   sepDims(1) = number of columns in the second product margin
  //
  // The storage order follows the same convention as sepgrid()/numFactor():
  //
  //   flat index = coord1 + sepDims(0) * coord2
  //
  // `sepCodes` stores the covariance-structure code for each margin in this
  // same order.  `sepDensityKinds` is a smaller C++ margin-kind code:
  //
  //   1 = dense correlation margin (currently cs, homcs, or us)
  //   2 = AR(1) correlation margin (currently ar1 or hetar1)
  //   3 = diagonal correlation margin (currently diag or homdiag)
  //
  // `sepDispatch` selects the C++ evaluator.  The current backend uses one
  // evaluator for any two margins that can be represented by correlation
  // matrices.
  //
  // `sepScaleMode` and `sepScaleSpec` describe how cell standard deviations are
  // built: one selected margin, one global scale, all scale-capable margins, or
  // an explicit product of selected margins.
  vector<int> sepDims;
  vector<int> sepCodes;
  vector<int> sepDensityKinds;
  vector<int> sepDispatch;
  vector<int> sepScaleMode;
  vector<int> sepScaleSpec;
  vector<int> sepDistStarts;
  vector<Type> sepDists;
  // Report output
  matrix<Type> corr;
  vector<Type> sd;
  matrix<Type> fact_load; // For rr case
};

template <class Type>
struct terms_t : vector<per_term_info<Type> > {
  terms_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int blockCode = (int) REAL(getListElement(y, "blockCode", &isNumericScalar))[0];
      int blockSize = (int) REAL(getListElement(y, "blockSize", &isNumericScalar))[0];
      int blockReps = (int) REAL(getListElement(y, "blockReps", &isNumericScalar))[0];
      int blockNumTheta = (int) REAL(getListElement(y, "blockNumTheta", &isNumericScalar))[0];
      int simCode = (int) REAL(getListElement(y, "simCode", &isNumericScalar))[0];
      int fullCor = (int) REAL(getListElement(y, "fullCor", &isNumericScalar))[0];
      (*this)(i).blockCode = blockCode;
      (*this)(i).blockSize = blockSize;
      (*this)(i).blockReps = blockReps;
      (*this)(i).blockNumTheta = blockNumTheta;
      (*this)(i).simCode = simCode;
      (*this)(i).fullCor = fullCor;
      // Optionally, pass time vector:
      SEXP t = getListElement(y, "times");
      if(!Rf_isNull(t)){
	RObjectTestExpectedType(t, &Rf_isNumeric, "times");
	(*this)(i).times = asVector<Type>(t);
      }
      // Optionally, pass distance matrix:
      SEXP d = getListElement(y, "dist");
      if(!Rf_isNull(d)){
	RObjectTestExpectedType(d, &Rf_isMatrix, "dist");
	(*this)(i).dist = asMatrix<Type>(d);
      }
      // Optionally, pass separable margin metadata:
      //
      // These fields are absent for all existing covariance structures.  Keeping
      // them optional means the rest of the template keeps using the same
      // `per_term_info` object without special wrappers.
      SEXP sdims = getListElement(y, "sepDims");
      if(!Rf_isNull(sdims)){
	RObjectTestExpectedType(sdims, &Rf_isNumeric, "sepDims");
	(*this)(i).sepDims = asVector<int>(sdims);
      }
      SEXP scodes = getListElement(y, "sepCodes");
      if(!Rf_isNull(scodes)){
	RObjectTestExpectedType(scodes, &Rf_isNumeric, "sepCodes");
	(*this)(i).sepCodes = asVector<int>(scodes);
      }
      SEXP skinds = getListElement(y, "sepDensityKinds");
      if(!Rf_isNull(skinds)){
	RObjectTestExpectedType(skinds, &Rf_isNumeric, "sepDensityKinds");
	(*this)(i).sepDensityKinds = asVector<int>(skinds);
      }
      SEXP sdispatch = getListElement(y, "sepDispatch");
      if(!Rf_isNull(sdispatch)){
	RObjectTestExpectedType(sdispatch, &Rf_isNumeric, "sepDispatch");
	(*this)(i).sepDispatch = asVector<int>(sdispatch);
      }
      SEXP sscalemode = getListElement(y, "sepScaleMode");
      if(!Rf_isNull(sscalemode)){
	RObjectTestExpectedType(sscalemode, &Rf_isNumeric, "sepScaleMode");
	(*this)(i).sepScaleMode = asVector<int>(sscalemode);
      }
      SEXP sscalespec = getListElement(y, "sepScaleSpec");
      if(!Rf_isNull(sscalespec)){
	RObjectTestExpectedType(sscalespec, &Rf_isNumeric, "sepScaleSpec");
	(*this)(i).sepScaleSpec = asVector<int>(sscalespec);
      }
      SEXP sdiststarts = getListElement(y, "sepDistStarts");
      if(!Rf_isNull(sdiststarts)){
	RObjectTestExpectedType(sdiststarts, &Rf_isNumeric, "sepDistStarts");
	(*this)(i).sepDistStarts = asVector<int>(sdiststarts);
      }
      SEXP sdists = getListElement(y, "sepDists");
      if(!Rf_isNull(sdists)){
	RObjectTestExpectedType(sdists, &Rf_isNumeric, "sepDists");
	(*this)(i).sepDists = asVector<Type>(sdists);
      }
    }
  }
};

template <class Type>
Type fill_separable_array(array<Type>& z, array<Type> &U,
			  const vector<Type>& cell_sd, int g) {
  Type logscale = 0;
  for (int k = 0; k < z.size(); k++) {
    z(k) = U(k, g) / cell_sd(k);
    logscale += log(cell_sd(k));
  }
  return logscale;
}

template <class Type, class Density0, class Density1>
Type separable_2d_nll(array<Type> &U, const vector<Type>& cell_sd,
		      Density0 density0, Density1 density1,
		      per_term_info<Type>& term) {
  // Evaluate a two-margin separable block without forming the full Kronecker
  // covariance matrix. TMB's SEPARABLE arguments are supplied in reverse
  // array-dimension order.
  int n0 = term.sepDims(0);
  int n1 = term.sepDims(1);
  vector<int> dim(2);
  dim << n0, n1;
  Type ans = 0;

  for (int g = 0; g < term.blockReps; g++) {
    array<Type> z(dim);
    Type logscale = fill_separable_array(z, U, cell_sd, g);
    ans += density::SEPARABLE(density1, density0)(z) + logscale;
  }
  return ans;
}

template <class Type, class Density0, class Density1, class Density2>
Type separable_3d_nll(array<Type> &U, const vector<Type>& cell_sd,
		      Density0 density0, Density1 density1, Density2 density2,
		      per_term_info<Type>& term) {
  int n0 = term.sepDims(0);
  int n1 = term.sepDims(1);
  int n2 = term.sepDims(2);
  vector<int> dim(3);
  dim << n0, n1, n2;
  Type ans = 0;

  for (int g = 0; g < term.blockReps; g++) {
    array<Type> z(dim);
    Type logscale = fill_separable_array(z, U, cell_sd, g);
    ans += density::SEPARABLE(density2,
			      density::SEPARABLE(density1, density0))(z) +
      logscale;
  }
  return ans;
}

template <class Type, class Density0, class Density1, class Density2,
	  class Density3>
Type separable_4d_nll(array<Type> &U, const vector<Type>& cell_sd,
		      Density0 density0, Density1 density1, Density2 density2,
		      Density3 density3, per_term_info<Type>& term) {
  int n0 = term.sepDims(0);
  int n1 = term.sepDims(1);
  int n2 = term.sepDims(2);
  int n3 = term.sepDims(3);
  vector<int> dim(4);
  dim << n0, n1, n2, n3;
  Type ans = 0;

  for (int g = 0; g < term.blockReps; g++) {
    array<Type> z(dim);
    Type logscale = fill_separable_array(z, U, cell_sd, g);
    ans += density::SEPARABLE(
      density3,
      density::SEPARABLE(density2,
			 density::SEPARABLE(density1, density0)))(z) +
      logscale;
  }
  return ans;
}

bool is_dense_corr_margin(int code) {
  return code == cs_covstruct || code == homcs_covstruct || code == us_covstruct;
}

bool is_diag_margin(int code) {
  return code == diag_covstruct || code == homdiag_covstruct;
}

bool is_toep_margin(int code) {
  return code == toep_covstruct || code == homtoep_covstruct;
}

bool is_spatial_margin(int code) {
  return code == ou_covstruct || code == exp_covstruct ||
    code == gau_covstruct || code == mat_covstruct;
}

template <class Type>
matrix<Type> compound_symmetry_corr(int n, Type corr_transf) {
  Type a = Type(1) / (Type(n) - Type(1));
  Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
  matrix<Type> corr(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      corr(i, j) = (i == j ? Type(1) : rho);
  return corr;
}

template <class Type>
matrix<Type> identity_corr(int n) {
  matrix<Type> corr(n, n);
  corr.setZero();
  for (int i = 0; i < n; i++) corr(i, i) = Type(1);
  return corr;
}

template <class Type>
void parse_separable_dense_margin(int code, int n, const vector<Type>& theta,
				  int& theta_pos, bool scale_here,
				  vector<Type>& margin_sd,
				  matrix<Type>& corr,
				  vector<Type>& us_corr_params) {
  if (!is_dense_corr_margin(code))
    error("unsupported dense margin for separable covariance structure");

  margin_sd.resize(n);
  margin_sd.fill(Type(1));
  if (scale_here) {
    if (code == homcs_covstruct) {
      margin_sd.fill(exp(theta(theta_pos++)));
    } else {
      vector<Type> logsd = theta.segment(theta_pos, n);
      theta_pos += n;
      margin_sd = exp(logsd);
    }
  }

  if (code == cs_covstruct || code == homcs_covstruct) {
    corr = compound_symmetry_corr(n, theta(theta_pos++));
  } else {
    int n_corr = n * (n - 1) / 2;
    us_corr_params = theta.segment(theta_pos, n_corr);
    theta_pos += n_corr;
  }
}

template <class Type>
void parse_separable_diag_margin(int code, int n, const vector<Type>& theta,
				 int& theta_pos, bool scale_here,
				 vector<Type>& margin_sd,
				 matrix<Type>& corr) {
  if (!is_diag_margin(code))
    error("unsupported diagonal margin for separable covariance structure");

  margin_sd.resize(n);
  margin_sd.fill(Type(1));
  if (scale_here) {
    if (code == homdiag_covstruct) {
      margin_sd.fill(exp(theta(theta_pos++)));
    } else {
      vector<Type> logsd = theta.segment(theta_pos, n);
      theta_pos += n;
      margin_sd = exp(logsd);
    }
  }
  corr = identity_corr<Type>(n);
}

template <class Type>
void parse_separable_toep_margin(int code, int n, const vector<Type>& theta,
				 int& theta_pos, bool scale_here,
				 vector<Type>& margin_sd,
				 matrix<Type>& corr) {
  if (!is_toep_margin(code))
    error("unsupported Toeplitz margin for separable covariance structure");

  margin_sd.resize(n);
  margin_sd.fill(Type(1));
  if (scale_here) {
    if (code == homtoep_covstruct) {
      margin_sd.fill(exp(theta(theta_pos++)));
    } else {
      vector<Type> logsd = theta.segment(theta_pos, n);
      theta_pos += n;
      margin_sd = exp(logsd);
    }
  }

  vector<Type> corr_params = theta.segment(theta_pos, n - 1);
  theta_pos += n - 1;
  corr_params = corr_params / sqrt(Type(1) + corr_params * corr_params);
  corr.resize(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      corr(i, j) = (i == j ? Type(1) :
		    corr_params((i > j ? i - j : j - i) - 1));
}

template <class Type>
matrix<Type> separable_margin_dist(per_term_info<Type>& term, int m) {
  int n = term.sepDims(m);
  if (term.sepDistStarts.size() != term.sepDims.size() ||
      term.sepDistStarts(m) < 0)
    error("separable spatial margin is missing distance metadata");
  int start = term.sepDistStarts(m);
  if (start + n * n > term.sepDists.size())
    error("separable spatial margin has invalid distance metadata");

  matrix<Type> dist(n, n);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      dist(i, j) = term.sepDists(start + i + n * j);
  return dist;
}

template <class Type>
void parse_separable_spatial_margin(int code, int n, const vector<Type>& theta,
				    int& theta_pos, bool scale_here,
				    const matrix<Type>& dist,
				    vector<Type>& margin_sd,
				    matrix<Type>& corr) {
  if (!is_spatial_margin(code))
    error("unsupported spatial margin for separable covariance structure");
  if (scale_here)
    error("spatial separable margins cannot carry scale parameters");

  margin_sd.resize(n);
  margin_sd.fill(Type(1));
  corr.resize(n, n);

  Type theta0 = theta(theta_pos++);
  Type theta1 = Type(0);
  if (code == mat_covstruct) theta1 = theta(theta_pos++);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      switch (code) {
      case ou_covstruct:
	corr(i, j) = (i == j ? Type(1) : exp(-exp(theta0) * dist(i, j)));
	break;
      case exp_covstruct:
	corr(i, j) = (i == j ? Type(1) : exp(-dist(i, j) * exp(-theta0)));
	break;
      case gau_covstruct:
	corr(i, j) = (i == j ? Type(1) : exp(-pow(dist(i, j), 2) *
					     exp(Type(-2) * theta0)));
	break;
      case mat_covstruct:
	corr(i, j) = (i == j ? Type(1) :
		      matern(dist(i, j), exp(theta0), exp(theta1)));
	break;
      default:
	error("unsupported spatial margin for separable covariance structure");
      }
    }
  }
}

template <class Type>
bool separable_margin_has_scale(per_term_info<Type>& term, int m) {
  for (int i = 0; i < term.sepScaleSpec.size(); i++)
    if (term.sepScaleSpec(i) == m) return true;
  return false;
}

template <class Type>
void check_separable_metadata(per_term_info<Type>& term) {
  if (term.sepDims.size() < 2 || term.sepCodes.size() != term.sepDims.size() ||
      term.sepDensityKinds.size() != term.sepDims.size() ||
      term.sepDispatch.size() != 1 ||
      term.sepScaleMode.size() != 1)
    error("separable covariance structure is missing margin metadata");
  int n = 1;
  for (int i = 0; i < term.sepDims.size(); i++)
    n *= term.sepDims(i);
  if (n != term.blockSize)
    error("separable dimensions do not match block size");
  int mode = term.sepScaleMode(0);
  if (mode < margin_sep_scale || mode > selected_product_sep_scale)
    error("unknown separable scale mode");
  if (mode == global_sep_scale && term.sepScaleSpec.size() != 0)
    error("separable global scale should not specify scale margins");
  if (mode != global_sep_scale && term.sepScaleSpec.size() == 0)
    error("separable scale mode is missing selected margins");
  for (int i = 0; i < term.sepScaleSpec.size(); i++)
    if (term.sepScaleSpec(i) < 0 ||
	term.sepScaleSpec(i) >= term.sepDims.size())
      error("separable margin scale index is out of range");
}

template <class Type>
vector<Type> separable_cell_sd(Type global_sd, const vector<Type>& sd0,
			       const vector<Type>& sd1) {
  int n0 = sd0.size();
  int n1 = sd1.size();
  vector<Type> cell_sd(n0 * n1);
  for (int i1 = 0; i1 < n1; i1++) {
    for (int i0 = 0; i0 < n0; i0++) {
      int k = i0 + n0 * i1;
      cell_sd(k) = global_sd * sd0(i0) * sd1(i1);
    }
  }
  return cell_sd;
}

template <class Type>
vector<Type> separable_cell_sd(Type global_sd,
			       const vector<vector<Type> >& margin_sd,
			       const vector<int>& dims) {
  int n = 1;
  for (int i = 0; i < dims.size(); i++) n *= dims(i);
  vector<Type> cell_sd(n);
  for (int k = 0; k < n; k++) {
    int rest = k;
    cell_sd(k) = global_sd;
    for (int m = 0; m < dims.size(); m++) {
      int coord = rest % dims(m);
      rest /= dims(m);
      cell_sd(k) *= margin_sd(m)(coord);
    }
  }
  return cell_sd;
}

template <class Type>
matrix<Type> ar1_corr(int n, Type phi) {
  matrix<Type> corr(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      corr(i, j) = pow(phi, abs(i - j));
  return corr;
}

template <class Type>
Type parse_ar1_phi(const vector<Type>& theta, int& theta_pos) {
  Type corr_transf = theta(theta_pos++);
  // Same transform used by existing glmmTMB AR1.
  return corr_transf / sqrt(Type(1) + pow(corr_transf, 2));
}

template <class Type>
void parse_separable_ar1_margin(int code, int n, const vector<Type>& theta,
				int& theta_pos, bool scale_here,
				vector<Type>& margin_sd, Type& phi) {
  if (code != ar1_covstruct && code != hetar1_covstruct)
    error("unsupported AR1 margin for separable covariance structure");

  margin_sd.resize(n);
  margin_sd.fill(Type(1));
  if (scale_here) {
    if (code != hetar1_covstruct)
      error("homogeneous AR1 separable margins cannot carry scale parameters");
    vector<Type> logsd = theta.segment(theta_pos, n);
    theta_pos += n;
    margin_sd = exp(logsd);
  }
  phi = parse_ar1_phi(theta, theta_pos);
}

template <class Type>
void report_separable_2d(const vector<Type>& cell_sd,
			 const matrix<Type>& corr0,
			 const matrix<Type>& corr1,
			 per_term_info<Type>& term) {
  // Build report objects for VarCorr().
  int n0 = term.sepDims(0);
  int n1 = term.sepDims(1);
  term.sd = cell_sd;
  if (term.fullCor == 1) {
    // Construct the full correlation only for reporting/testing.
    int n = n0 * n1;
    term.corr.resize(n, n);
    for (int b1 = 0; b1 < n1; b1++) {
      for (int a1 = 0; a1 < n0; a1++) {
        int k1 = a1 + n0 * b1;
        for (int b2 = 0; b2 < n1; b2++) {
          for (int a2 = 0; a2 < n0; a2++) {
            int k2 = a2 + n0 * b2;
            term.corr(k1, k2) = corr0(a1, a2) * corr1(b1, b2);
          }
        }
      }
    }
  } else {
    // Follow the existing compact-report convention for structured covariance
    // terms that do not store the whole correlation matrix.
    term.corr.resize(1,1);
    term.corr(0,0) = NAN;
  }
}

template <class Type>
void report_separable_dense(const vector<Type>& cell_sd,
			    const matrix<Type>& corr,
			    per_term_info<Type>& term) {
  term.sd = cell_sd;
  if (term.fullCor == 1) {
    term.corr = corr;
  } else {
    term.corr.resize(1,1);
    term.corr(0,0) = NAN;
  }
}

template <class Type, class Density0, class Density1>
Type eval_separable_2d(array<Type> &U, const vector<Type>& cell_sd,
		       Density0 density0, Density1 density1,
		       const matrix<Type>& corr0, const matrix<Type>& corr1,
		       per_term_info<Type>& term) {
  Type ans = separable_2d_nll(U, cell_sd, density0, density1, term);
  DISABLE_AD {
    report_separable_2d(cell_sd, corr0, corr1, term);
  }
  return ans;
}

template <class Type, class Density0, class Density1, class Density2>
Type eval_separable_3d(array<Type> &U, const vector<Type>& cell_sd,
		       Density0 density0, Density1 density1, Density2 density2,
		       const matrix<Type>& corr, per_term_info<Type>& term) {
  Type ans = separable_3d_nll(U, cell_sd, density0, density1, density2, term);
  DISABLE_AD {
    report_separable_dense(cell_sd, corr, term);
  }
  return ans;
}

template <class Type, class Density0, class Density1, class Density2,
	  class Density3>
Type eval_separable_4d(array<Type> &U, const vector<Type>& cell_sd,
		       Density0 density0, Density1 density1, Density2 density2,
		       Density3 density3, const matrix<Type>& corr,
		       per_term_info<Type>& term) {
  Type ans = separable_4d_nll(U, cell_sd, density0, density1, density2,
			      density3, term);
  DISABLE_AD {
    report_separable_dense(cell_sd, corr, term);
  }
  return ans;
}

template <class Type>
Type eval_separable_dense(array<Type> &U, const vector<Type>& cell_sd,
			  const matrix<Type>& corr,
			  per_term_info<Type>& term) {
  density::MVNORM_t<Type> density(corr);
  Type ans = 0;
  for (int g = 0; g < term.blockReps; g++) {
    vector<Type> z = U.col(g);
    Type logscale = 0;
    for (int k = 0; k < z.size(); k++) {
      z(k) /= cell_sd(k);
      logscale += log(cell_sd(k));
    }
    ans += density(z) + logscale;
  }
  DISABLE_AD {
    report_separable_dense(cell_sd, corr, term);
  }
  return ans;
}

bool is_separable_corr_matrix_kind(int kind) {
  return kind == dense_corr_sep || kind == ar1_sep || kind == diag_sep ||
    kind == spatial_sep || kind == toep_sep;
}

template <class Type>
struct sep_corr_margin_pars {
  vector<Type> sd;
  matrix<Type> corr;
};

template <class Type>
sep_corr_margin_pars<Type> parse_separable_corr_margin(int code, int kind,
						       int n,
						       const vector<Type>& theta,
						       int& theta_pos,
						       bool scale_here,
						       per_term_info<Type>& term,
						       int m) {
  sep_corr_margin_pars<Type> out;

  if (kind == diag_sep) {
    parse_separable_diag_margin(code, n, theta, theta_pos, scale_here,
				out.sd, out.corr);
  } else if (kind == ar1_sep) {
    Type phi = Type(0);
    parse_separable_ar1_margin(code, n, theta, theta_pos, scale_here,
			       out.sd, phi);
    out.corr = ar1_corr(n, phi);
  } else if (kind == dense_corr_sep) {
    vector<Type> us_corr_params(0);
    parse_separable_dense_margin(code, n, theta, theta_pos, scale_here,
				 out.sd, out.corr, us_corr_params);
    if (code == us_covstruct) {
      density::UNSTRUCTURED_CORR_t<Type> us_density(us_corr_params);
      out.corr = us_density.cov();
    }
  } else if (kind == toep_sep) {
    parse_separable_toep_margin(code, n, theta, theta_pos, scale_here,
				out.sd, out.corr);
  } else if (kind == spatial_sep) {
    matrix<Type> dist = separable_margin_dist(term, m);
    parse_separable_spatial_margin(code, n, theta, theta_pos, scale_here,
				   dist, out.sd, out.corr);
  } else {
    error("unsupported separable margin for correlation-matrix evaluator");
  }

  return out;
}

template <class Type>
struct sep_corr_corr_pars {
  vector<Type> cell_sd;
  matrix<Type> corr0;
  matrix<Type> corr1;
  matrix<Type> corr2;
  matrix<Type> corr3;
  matrix<Type> corr;
};

template <class Type>
matrix<Type> kronecker_corr(const matrix<Type>& left,
			    const matrix<Type>& right) {
  matrix<Type> out(left.rows() * right.rows(), left.cols() * right.cols());
  for (int i = 0; i < left.rows(); i++) {
    for (int j = 0; j < left.cols(); j++) {
      for (int k = 0; k < right.rows(); k++) {
	for (int l = 0; l < right.cols(); l++) {
	  out(i * right.rows() + k, j * right.cols() + l) =
	    left(i, j) * right(k, l);
	}
      }
    }
  }
  return out;
}

template <class Type>
sep_corr_corr_pars<Type> parse_separable_corr_corr(const vector<Type>& theta,
						   per_term_info<Type>& term) {
  sep_corr_corr_pars<Type> out;
  int n_margin = term.sepDims.size();
  for (int m = 0; m < n_margin; m++) {
    int kind = term.sepDensityKinds(m);
    if (!is_separable_corr_matrix_kind(kind))
      error("unsupported separable margin for correlation-matrix evaluator");
  }

  Type global_sd = Type(1);
  int theta_pos = 0;
  if (term.sepScaleMode(0) == global_sep_scale)
    global_sd = exp(theta(theta_pos++));

  vector<vector<Type> > margin_sd(n_margin);
  for (int m = 0; m < n_margin; m++) {
    sep_corr_margin_pars<Type> margin =
      parse_separable_corr_margin(term.sepCodes(m), term.sepDensityKinds(m),
				  term.sepDims(m), theta, theta_pos,
				  separable_margin_has_scale(term, m),
				  term, m);
    margin_sd(m) = margin.sd;
    if (m == 0) {
      out.corr0 = margin.corr;
      out.corr = margin.corr;
    } else if (m == 1) {
      out.corr1 = margin.corr;
      out.corr = kronecker_corr(margin.corr, out.corr);
    } else if (m == 2) {
      out.corr2 = margin.corr;
      out.corr = kronecker_corr(margin.corr, out.corr);
    } else if (m == 3) {
      out.corr3 = margin.corr;
      out.corr = kronecker_corr(margin.corr, out.corr);
    } else {
      out.corr = kronecker_corr(margin.corr, out.corr);
    }
  }
  out.cell_sd = separable_cell_sd(global_sd, margin_sd, term.sepDims);

  if (theta_pos != theta.size())
    error("separable covariance theta parsing mismatch");

  return out;
}


// compute log-likelihood of b (conditional modes) conditional on theta (var/cov)
//  for a specified random-effects term 
template <class Type>
Type termwise_nll(array<Type> &U, vector<Type> theta, per_term_info<Type>& term, bool do_simulate = false) {
  Type ans = 0;
  if (term.blockCode == diag_covstruct) {
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    for(int i = 0; i < term.blockReps; i++) {
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
      if (do_simulate) {
	// FIXME this can be abstracted as a more general function (zero and fix are always the same,
	// random_simcode will differ by type) ¿can we make 'covariance structure' a class, with nll and simulate methods??
        switch(term.simCode) {
	case fix_simcode:
          // do nothing, leave U values as is
	  break;
	case zero_simcode:
	  for (int j=0; j < U.rows(); j++) {
	    U(j,i) = Type(0);
	  };
	  break;
	case random_simcode:
          U.col(i) = rnorm(Type(0), sd);
          break;
	default: error ("unknown simcode");
	} // simcode
      } // simulate
    } // loop over blocks
    term.sd = sd; // For report
  }
  else if (term.blockCode == homdiag_covstruct) {
    // case: homdiag_covstruct
    Type sd = exp(theta(0));
    for(int i = 0; i < term.blockReps; i++){
          for (int j = 0; j < U.rows(); j++) {
	    ans -= dnorm(Type(U(j,i)), Type(0), sd, true);
	    if (do_simulate) {
	      if (term.simCode != random_simcode) {
		Rf_error("simcode not yet implemented for homdiag cov struct");
	      }
	      U(j,i) = rnorm(Type(0), sd);
	    }	      
	  }
    }
    int n = term.blockSize;
    vector<Type> sdvec(n);
    for(int i = 0; i < n; i++) {
      sdvec(i) = sd;
    }
    
    term.sd = sdvec; // For report
  }

  else if (term.blockCode == us_covstruct){
    // case: us_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    vector<Type> sd = exp(logsd);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        switch(term.simCode) {
        case fix_simcode:
          // do nothing, leave U values as is
          break;
        case zero_simcode:
          for (int j=0; j < U.rows(); j++) {
            U(j,i) = Type(0);
          };
          break;
        case random_simcode:
          U.col(i) = sd * nldens.simulate();
          break;
        default: error ("unknown simcode");
        } // simcode
        
      }
    }

    // FIXME: DRY/make this into a function or macro
    SET_COR;
    
    term.sd = sd;             // For report
  }
  else if (term.blockCode == cs_covstruct || term.blockCode == homcs_covstruct) {
    // case: cs_covstruct or homcs_covstruct
    int n = term.blockSize;
    vector<Type> logsd(n);
    for (int i = 0; i < n; i++) {
      if (term.blockCode == cs_covstruct) {
        logsd(i) = theta(i);
      } else {
        logsd(i) = theta(0);
      }
    }
    Type corr_transf = theta.tail(1)(0);
    vector<Type> sd = exp(logsd);
    Type a = Type(1) / (Type(n) - Type(1));
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        corr(i,j) = (i==j ? Type(1) : rho);
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        if (term.simCode != random_simcode) {
          Rf_error("simcode not yet implemented for cs cov struct");
        }
        U.col(i) = sd * nldens.simulate();
      }
    }
    SET_COR;
    term.sd = sd;             // For report
  }
  else if (term.blockCode == toep_covstruct || term.blockCode == homtoep_covstruct) {
    // case: toep_covstruct or homtoep_covstruct
    int n = term.blockSize;
    vector<Type> logsd(n);
    for (int i = 0; i < n; i++) {
      if (term.blockCode == toep_covstruct) {
        logsd(i) = theta(i);
      } else {
        logsd(i) = theta(0);
      }
    }
    vector<Type> sd = exp(logsd);
    vector<Type> corr_params = theta.tail(n-1);
    corr_params = corr_params / sqrt(Type(1.0) + corr_params * corr_params);
    matrix<Type> corr(n,n);
    for (int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        corr(i,j) = (i==j ? Type(1) :
                     corr_params( (i > j ? i-j : j-i) - 1 ) );
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        if (term.simCode != random_simcode) {
          Rf_error("simcode not yet implemented for toep and homtoep cov struct");
        }
        U.col(i) = sd * nldens.simulate();
      }
    }
    SET_COR;
    term.sd = sd;             // For report
  }
  else if (term.blockCode == ar1_covstruct ||
	   term.blockCode == hetar1_covstruct) {
    // case: ar1_covstruct
    //  * NOTE: Valid parameter space is phi in [-1, 1]
    //  * NOTE: 'times' not used as we assume unit distance between consecutive time points.
    int n = term.blockSize;
    vector<Type> logsd = theta.head(term.blockNumTheta-1);
    Type corr_transf = theta(term.blockNumTheta-1);
    Type phi = corr_transf / sqrt(1.0 + pow(corr_transf, 2));
    vector<Type> sd = exp(logsd);
    
    for(int j = 0; j < term.blockReps; j++){
      // Initialize
      if (term.blockCode == hetar1_covstruct) {
        ans -= dnorm(U(0, j) / sd(0), Type(0), Type(1), true);
        ans += logsd(0);
      } else { // ar1_covstruct
        ans -= dnorm(U(0, j), Type(0), sd(0), true);
      }
      if (do_simulate) {
        switch(term.simCode) {
          case fix_simcode:
            break;
          case zero_simcode:
            U(0,j) = Type(0);
            break;
          case random_simcode:
            U(0, j) = rnorm(Type(0), sd(0));
            break;
        }
      }
      for(int i=1; i<n; i++) {

        if (term.blockCode == hetar1_covstruct) {
          ans -= dnorm(U(i, j) / sd(i), phi * U(i-1, j) / sd(i-1), sqrt(1 - phi*phi), true);
          ans += logsd(i);
        } else { // ar1_covstruct
          ans -= dnorm(U(i, j), phi * U(i-1, j), sd(0) * sqrt(1 - phi*phi), true);
        }
      
        if (do_simulate) {
	        switch(term.simCode) {
            case fix_simcode:
              // do nothing, leave U values as is
	            break;
            case zero_simcode:
              for (int i=0; i < U.rows(); i++) {
                U(i, j) = Type(0);
	            };
              break;
            case random_simcode:
              if (term.blockCode == hetar1_covstruct) {
                  U(i, j) = sd(i) * rnorm(phi * U(i-1, j) / sd(i-1), sqrt(1 - phi*phi));
              } else {
                  U(i, j) = rnorm(phi * U(i-1, j), sd(0) * sqrt(1 - phi*phi));
              }
              break;
            default: error ("unknown simcode");
	        } // term.simCode
      	} // do_simulate
      } // loop over lags
    } // loop over blocks
    DISABLE_AD { // Disable AD for this part
      if (term.fullCor == 0) {
        // report *only* phi in the corr struct
        term.corr.resize(1,1);
        term.corr(0,0) = phi;
      } else {
        term.corr.resize(n,n);
        for (int i=0; i<n; i++) {
          for(int j=0; j<n; j++){
            term.corr(i,j) = pow(phi, abs(i-j));
          }
        }
      }
      if (term.blockCode == hetar1_covstruct) {
        term.sd.resize(n);
        for(int i=0; i<n; i++){
          term.sd(i) = sd(i);
        }
      } else { // ar1_covstruct
        term.sd.resize(1);
        term.sd(0) = sd(0);
      }
    } // DISABLE_AD 
  } // [het]ar1_covstruct
  else if (term.blockCode == ou_covstruct){
    // case: ou_covstruct
    //  * NOTE: this is the continuous time version of ar1.
    //          One-step correlation must be non-negative
    //  * NOTE: 'times' assumed sorted !
    int n = term.times.size();
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type sd = exp(logsd);
    for(int i = 0; i < term.blockReps; i++){
      for(int j=1; j<n; j++){
        Type rho = exp(-exp(corr_transf) * (term.times(j) - term.times(j-1)));
        ans -= dnorm(U(j, i), rho * U(j-1, i), sd * sqrt(1 - rho*rho), true);
      }
      ans -= dnorm(U(0, i), Type(0), sd, true);   // Initialize
      if (do_simulate) {
        switch(term.simCode) {
        case fix_simcode:
          // do nothing, leave U values as is
          break;
        case zero_simcode:
          for (int j=0; j < U.rows(); j++) {
            U(j,i) = Type(0);
          };
          break;
        case random_simcode:
          U(0, i) = rnorm(Type(0), sd);
          for(int j=1; j<n; j++){
            Type rho = exp(-exp(corr_transf) * (term.times(j) - term.times(j-1)));
            U(j, i) = rnorm( rho * U(j-1, i), sd * sqrt(1 - rho*rho));
          }
          break;
        default: error ("unknown simcode");
        }
      } // do_simulate

      
    }
    // Report only the sd and the decay parameter
    DISABLE_AD { // Disable AD for this part
      if (term.fullCor==1) {
        term.corr.resize(n,n);
        term.sd.resize(n);
        for(int i=0; i<n; i++) {
          term.sd(i) = sd;
          for(int j=0; j<n; j++){
            term.corr(i,j) =
              exp(-exp(corr_transf) * CppAD::abs(term.times(i) - term.times(j)));
          }
        }
      } else {
        term.corr.resize(1,1);
        term.sd.resize(1);
        term.sd(0) = sd;
        term.corr(0,0) = exp(corr_transf);
      }
    }
  } // OU covstruct
  
  // Spatial correlation structures
  else if (term.blockCode == exp_covstruct ||
           term.blockCode == gau_covstruct ||
           term.blockCode == mat_covstruct){
    int n = term.blockSize;
    matrix<Type> dist = term.dist;
    if(! ( dist.cols() == n && dist.rows() == n ) )
      error ("Dimension of distance matrix must equal blocksize.");
    // First parameter is sd
    Type sd = exp( theta(0) );
    // Setup correlation matrix
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        switch (term.blockCode) {
        case exp_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -dist(i,j) * exp(-theta(1)) ) );
          break;
        case gau_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -pow(dist(i,j),2) * exp(-2. * theta(1)) ) );
          break;
        case mat_covstruct:
          corr(i,j) = (i==j ? Type(1) : matern( dist(i,j),
                                                exp(theta(1)) /* range */,
                                                exp(theta(2)) /* smoothness */) );
          break;
        default:
          error("Not implemented");
        }
      }
    }
    density::MVNORM_t<Type> nldens(corr);
    density::SCALE_t<density::MVNORM_t<Type> > scnldens = density::SCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
	if (term.simCode != random_simcode) {
	  Rf_error("simcode not yet implemented for spatial cov structs");
	}
        U.col(i) = sd * nldens.simulate();
      }
    }
    if (term.fullCor==1) term.corr = corr;   // For report
    term.sd.resize(n);  // For report
    term.sd.fill(sd);
  }
  else if (term.blockCode == rr_covstruct){
    // case: reduced rank

    // computing log-likelihood based on *spherical* (iid N(0,1)) random effects
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), 1, true).sum();
      if (do_simulate) {
        switch(term.simCode) {
        case fix_simcode:
          // do nothing, leave U values as is
          break;
        case zero_simcode:
          for (int j=0; j < U.rows(); j++) {
            U(j,i) = Type(0);
          };
          break;
        case random_simcode:
          U.col(i) = rnorm(U.rows(), Type(0), Type(1));
          break;
        default: error ("unknown simcode");
        } // simcode
      }
    }

    // now construct the factor matrix and convert the spherical random
    //  effects back to the 'data scale', and *replace them* in the U matrix

    // constructing the factor loadings matrix
    int p = term.blockSize;
    int nt = theta.size();
    int rank = (2*p + 1 -  (int)sqrt(pow(2.0*p + 1, 2) - 8*nt) ) / 2 ;
    matrix<Type> Lambda(p, rank);
    vector<Type> lam_diag = theta.head(rank);
    vector<Type> lam_lower = theta.tail(nt - rank);
    for (int j = 0; j < rank; j++){
      for (int i = 0; i < p; i++){
        if (j > i)
          Lambda(i, j) = 0;
        else if(i == j)
          Lambda(i, j) = lam_diag(j);
        else
          Lambda(i, j) = lam_lower(j*p - (j + 1)*j/2 + i - 1 - j); //Fills by column
      }
    }

    // transforming u to b by multiplying by the loadings matrix
    // if simcode is 'fixed', do **not** multiply by loadings matrix
    // (i.e. user is assumed to be passing the 'non-spherical' latent variables)
    if (term.simCode != fix_simcode) {
      for(int i = 0; i < term.blockReps; i++){
        vector<Type> usub = U.col(i).segment(0, rank);
        U.col(i) = Lambda * usub;
      }
    }

    // computing the correlation matrix and std devs
    // (the same D^(-1/2) L L^T D^(-1/2) transformation that we use for correlations
    term.fact_load = Lambda;
    if (term.fullCor==1) {
      DISABLE_AD {
        term.corr = Lambda * Lambda.transpose();
        term.sd = term.corr.diagonal().array().sqrt();
        term.corr.array() /= (term.sd.matrix() * term.sd.matrix().transpose()).array();
      }
    }
  }
  else if (term.blockCode == propto_covstruct){
    // case: propto_covstruct
    int n = term.blockSize;
    Type loglambda = theta( theta.size() - 1);
    vector<Type> logsd = theta.head(n);
    vector<Type> sd =  exp(logsd + loglambda/2) ;
    vector<Type> corr_transf = theta.segment(n, theta.size() - n - 1);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    SET_COR;
    term.sd = sd;             // For report
  }
  else if (term.blockCode == equalto_covstruct){
    // case: equalto_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> sd =  exp(logsd);
    vector<Type> corr_transf = theta.segment(n, theta.size() - n);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == separable_covstruct) {
    // R validates the separable margins and supplies an evaluator code
    // (`sepDispatch`) plus coordinate-order metadata (`sepDensityKinds`,
    // `sepScaleMode`, `sepScaleSpec`).  Products up to four margins use
    // nested TMB SEPARABLE calls; longer products use a dense fallback.
    if (do_simulate)
      error("simulation is not yet implemented for separable covariance structures");

    check_separable_metadata(term);
    switch (term.sepDispatch(0)) {
    case corr_corr_dispatch: {
      sep_corr_corr_pars<Type> sep = parse_separable_corr_corr(theta, term);
      if (term.sepDims.size() == 2) {
	density::MVNORM_t<Type> density0(sep.corr0);
	density::MVNORM_t<Type> density1(sep.corr1);
	ans += eval_separable_2d(U, sep.cell_sd, density0, density1,
				 sep.corr0, sep.corr1, term);
      } else if (term.sepDims.size() == 3) {
	density::MVNORM_t<Type> density0(sep.corr0);
	density::MVNORM_t<Type> density1(sep.corr1);
	density::MVNORM_t<Type> density2(sep.corr2);
	ans += eval_separable_3d(U, sep.cell_sd, density0, density1,
				 density2, sep.corr, term);
      } else if (term.sepDims.size() == 4) {
	density::MVNORM_t<Type> density0(sep.corr0);
	density::MVNORM_t<Type> density1(sep.corr1);
	density::MVNORM_t<Type> density2(sep.corr2);
	density::MVNORM_t<Type> density3(sep.corr3);
	ans += eval_separable_4d(U, sep.cell_sd, density0, density1,
				 density2, density3, sep.corr, term);
      } else {
	ans += eval_separable_dense(U, sep.cell_sd, sep.corr, term);
      }
      break;
    }
    default:
      error("unsupported separable covariance dispatch");
    }
  }
  else error("covStruct not implemented!");
  return ans;
}

template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
		  vector<per_term_info<Type> >& terms,
                  bool do_simulate = false) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    // Note: 'blockNumTheta=0' ==> Same parameters as previous term.
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i), do_simulate);
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_MATRIX(X);
  bool sparseX = X.rows()==0 && X.cols()==0;
  DATA_SPARSE_MATRIX(Z);
  DATA_MATRIX(Xzi);
  bool sparseXzi = Xzi.rows()==0 && Xzi.cols()==0;
  DATA_SPARSE_MATRIX(Zzi);
  DATA_MATRIX(Xdisp);
  bool sparseXdisp = Xdisp.rows()==0 && Xdisp.cols()==0;
  DATA_SPARSE_MATRIX(Zdisp);
  DATA_VECTOR(yobs);
  DATA_VECTOR(size); //only used in binomial
  DATA_VECTOR(weights);
  DATA_VECTOR(offset);
  DATA_VECTOR(zioffset);
  DATA_VECTOR(dispoffset);

  // Define covariance structure for the conditional model
  DATA_STRUCT(terms, terms_t);

  // Define covariance structure for the zero inflation
  DATA_STRUCT(termszi, terms_t);

  // Define covariance structure for the dispersion model
  DATA_STRUCT(termsdisp, terms_t);
  
  // Parameters related to design matrices
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(betazi);
  PARAMETER_VECTOR(betadisp);
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(bzi);
  PARAMETER_VECTOR(bdisp);
  
  // Joint vector of covariance parameters
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(thetazi);
  PARAMETER_VECTOR(thetadisp);
  
  // Extra family specific parameters
  // tweedie, t, ordbetareg, nbinom12, skewnormal;
  // see .extraParamFamilies in R code for expected length
  PARAMETER_VECTOR(psi);

  DATA_INTEGER(family);
  DATA_INTEGER(link);

  // Flags
  DATA_INTEGER(ziPredictCode);
  bool zi_flag = (betazi.size() > 0);
  // 0 = no prediction; 1 = predictions on link scale; 2 = predictions on
  // data scale; 3 = predictions of latent variables (b)
  DATA_INTEGER(doPredict);
  DATA_IVECTOR(whichPredict);
  // One-Step-Ahead (OSA) residuals
  DATA_VECTOR_INDICATOR(keep, yobs);

  // Prior info

  DATA_IVECTOR(prior_distrib);    // specify distribution
  DATA_IVECTOR(prior_whichpar);   // specify parameter
  DATA_IVECTOR(prior_elstart);    // starting element index
  DATA_IVECTOR(prior_elend);      // ending element index
  DATA_IVECTOR(prior_npar);       // number of parameters (based on prior distrib)
  DATA_VECTOR(prior_params);      // specify parameters (concatenated)
  
  // Joint negative log-likelihood
  Type jnll=0;

  // Random effects
  PARALLEL_REGION jnll += allterms_nll(b, theta, terms, this->do_simulate);
  PARALLEL_REGION jnll += allterms_nll(bzi, thetazi, termszi, this->do_simulate);
  PARALLEL_REGION jnll += allterms_nll(bdisp, thetadisp, termsdisp, this->do_simulate);

  // int nz = 0, nnz=0;
  // for (int i=0; i<b.size(); i++) {
  //   if (b(i)==0) nz++; else nnz++;
  // }
  // printf("b nz = %d, nnz = %d\n", nz, nnz);
  
  // Linear predictor
  vector<Type> eta = Z * b + offset;
  if (!sparseX) {
    eta += X*beta;
  } else {
    DATA_SPARSE_MATRIX(XS);
    eta += XS*beta;
  }
  vector<Type> etazi = Zzi * bzi + zioffset;
  if (!sparseXzi) {
    etazi += Xzi*betazi;
  } else {
    DATA_SPARSE_MATRIX(XziS);
    etazi += XziS*betazi;
  }
  vector<Type> etadisp = Zdisp * bdisp + dispoffset;
  if (!sparseXdisp) {
    etadisp += Xdisp*betadisp;
  } else {
    DATA_SPARSE_MATRIX(XdispS);
    etadisp += XdispS*betadisp;
  }

  // Apply link
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);
  vector<Type> pz = invlogit(etazi);
  vector<Type> phi = exp(etadisp);
  vector<Type> log_nzprob(eta.size());
  if (!trunc_Family(family)) {
    log_nzprob.setZero();
  } else {
    for (int i = 0; i < log_nzprob.size(); i++) {
      log_nzprob(i) =  calc_log_nzprob(mu(i), phi(i), eta(i), etadisp(i),
				       family, link);
    }
  }


// "zero-truncated" likelihood: ignore zeros in positive distributions
// exact zero: use for positive distributions (Gamma, beta)
#define zt_lik_zero(x,loglik_exp) (zi_flag && (x == Type(0)) ? -INFINITY : loglik_exp)
// close to zero: use for count data (cf binomial()$initialize)
#define zt_lik_nearzero(x,loglik_exp) (zi_flag && (x < Type(0.001)) ? -INFINITY : loglik_exp)

  // Observation likelihood
  Type s1, s2, s3;
  Type tmp_loglik;

  for (int i=0; i < yobs.size(); i++) PARALLEL_REGION {
    if ( !glmmtmb::isNA(yobs(i)) ) {
      switch (family) {
      case gaussian_family:
        tmp_loglik = dnorm(yobs(i), mu(i), phi(i), true);
        SIMULATE{yobs(i) = rnorm(mu(i), phi(i));}
        break;
      case skewnormal_family:
        s1 = mu(i);
        s2 = phi(i);
        s3 = psi(0);
        tmp_loglik = glmmtmb::dskewnorm(yobs(i), s1, s2, s3, true);
        SIMULATE{yobs(i) = glmmtmb::rskewnorm(s1, s2, s3);}
        break;
      case poisson_family:
        tmp_loglik = dpois(yobs(i), mu(i), true);
        SIMULATE{yobs(i) = rpois(mu(i));}
        break;
      case binomial_family:
        s1 = logit_inverse_linkfun(eta(i), link); // logit(p)
        tmp_loglik = dbinom_robust(yobs(i), size(i), s1, true);
        SIMULATE{yobs(i) = rbinom(size(i), mu(i));}
        break;
      case Gamma_family:
        s1 = phi(i);           // shape
        s2 = mu(i) / phi(i);   // scale
        tmp_loglik = zt_lik_zero(yobs(i),dgamma(yobs(i), s1, s2, true));
        SIMULATE{yobs(i) = rgamma(s1, s2);}
        break;
      case beta_family:
        // parameterization after Ferrari and Cribari-Neto 2004, betareg package
        s1 = mu(i)*phi(i);
        s2 = (Type(1)-mu(i))*phi(i);
        tmp_loglik = zt_lik_zero(yobs(i),dbeta(yobs(i), s1, s2, true));
        SIMULATE{yobs(i) = rbeta(s1, s2);}
        break;
      case ordbeta_family:
	// https://github.com/saudiwin/ordbetareg_pack/blob/master/R/modeling.R#L565-L573
	if (yobs(i) == 0.0) {
	  tmp_loglik = log1m_inverse_linkfun(eta(i) - psi(0), logit_link);
	  // std::cout << "zero " << asDouble(eta(i)) << " " << asDouble(psi(0)) << " " << asDouble(tmp_loglik) << std::endl;
	} else if (yobs(i) == 1.0) {
	  tmp_loglik = log_inverse_linkfun(eta(i) - psi(1), logit_link);
	  // std::cout << "one " << asDouble(eta(i)) << " " << asDouble(psi(1)) << " " << asDouble(tmp_loglik) << std::endl;
	} else {
	  s1 = mu(i)*phi(i);
	  s2 = (Type(1)-mu(i))*phi(i);
	  s3 = logspace_sub(log_inverse_linkfun(eta(i) - psi(0), logit_link),
			    log_inverse_linkfun(eta(i) - psi(1), logit_link));
	  tmp_loglik = s3 + dbeta(yobs(i), s1, s2, true);

	  // std::cout << "middle " << asDouble(eta(i)) << " " << asDouble(psi(0)) << " " << asDouble(psi(1)) << " " << asDouble(s3) << " " << asDouble(tmp_loglik) << " " << asDouble(s1) << " " << asDouble(s2) << " " << asDouble(mu(i)) << " " << asDouble(phi(i)) << std::endl;
	}
	SIMULATE{
	  s3 = invlogit(psi(0) - eta(i));
	  if (runif(Type(0), Type(1)) < s3) {
	    yobs(i) = 0;
	  } else {
	    s3 = invlogit(eta(i) - psi(1));
	    if (runif(Type(0), Type(1)) < s3) {
	      yobs(i) = 1;
	    } else {
	      s1 = mu(i)*phi(i);
	      s2 = (Type(1)-mu(i))*phi(i);
	      yobs(i) = rbeta(s1, s2);
	    }
	  }
	}
	break;
      case betabinomial_family:
        // Transform to logit scale independent of link
        s3 = logit_inverse_linkfun(eta(i), link); // logit(p)
        // Was: s1 = mu(i) * phi(i);
        s1 = log_inverse_linkfun( s3, logit_link) + log(phi(i)); // s1 = log(mu*phi)
        // Was: s2 = (Type(1) - mu(i)) * phi(i);
        s2 = log_inverse_linkfun(-s3, logit_link) + log(phi(i)); // s2 = log((1-mu)*phi)
        tmp_loglik = glmmtmb::dbetabinom_robust(yobs(i), s1, s2, size(i), true);
        SIMULATE {
          yobs(i) = rbinom(size(i), rbeta(exp(s1), exp(s2)) );
        }
        break;
      case nbinom1_family:
      case truncated_nbinom1_family:
        // Was:
        //   s1 = mu(i);
        //   s2 = mu(i) * (Type(1)+phi(i));  // (1+phi) guarantees that var >= mu
        //   tmp_loglik = dnbinom2(yobs(i), s1, s2, true);
        s1 = log_inverse_linkfun(eta(i), link);          // log(mu)
        s2 = s1 + etadisp(i) ;                              // log(var - mu)
        tmp_loglik = dnbinom_robust(yobs(i), s1, s2, true);
	if (family != truncated_nbinom1_family) {
	SIMULATE {
	  s1 = mu(i);
	  s2 = mu(i) * (Type(1)+phi(i));  // (1+phi) guarantees that var >= mu
	  yobs(i) = rnbinom2(s1, s2);
	  }
	} else {
          tmp_loglik -= log_nzprob(i);
	  tmp_loglik = zt_lik_nearzero(yobs(i), tmp_loglik);
          SIMULATE{
            s1 = mu(i)/phi(i); //sz
	    yobs(i) = glmmtmb::rtruncated_nbinom(asDouble(s1), 0, asDouble(mu(i)));
          }
        }
        break;
      case nbinom2_family:
      case truncated_nbinom2_family:
        s1 = log_inverse_linkfun(eta(i), link);          // log(mu)
        s2 = 2. * s1 - etadisp(i) ;                         // log(var - mu)
        tmp_loglik = dnbinom_robust(yobs(i), s1, s2, true);
        SIMULATE {
          s1 = mu(i);
          s2 = mu(i) * (Type(1) + mu(i) / phi(i));
          yobs(i) = rnbinom2(s1, s2);
        }
        if (family == truncated_nbinom2_family) {
          tmp_loglik -= log_nzprob(i);
          tmp_loglik = zt_lik_nearzero( yobs(i), tmp_loglik);
          SIMULATE{
		  yobs(i) = glmmtmb::rtruncated_nbinom(asDouble(phi(i)), 0, asDouble(mu(i)));
          }
        }
        break;
      case nbinom12_family:
	s1 = log_inverse_linkfun(eta(i), link);          // log(mu)
	// log(var - mu) = log(mu) + log(phi + mu/psi)
	s2 = s1 + logspace_add(etadisp(i), s1 - psi(0));
	tmp_loglik = dnbinom_robust(yobs(i), s1, s2, true);
	SIMULATE{
	  s1 = mu(i);
	  s2 = mu(i) * (Type(1)+phi(i) + mu(i)/exp(psi(0)));
	  yobs(i) = rnbinom2(s1, s2);
	}
	break;
      case truncated_poisson_family:
        tmp_loglik = dpois(yobs(i), mu(i), true) - log_nzprob(i);
        tmp_loglik = zt_lik_nearzero(yobs(i), tmp_loglik);
        SIMULATE{
		yobs(i) = glmmtmb::rtruncated_poisson(0, asDouble(mu(i)));
        }
        break;
     case genpois_family:
        s1 = mu(i) / sqrt(phi(i)); //theta
        s2 = Type(1) - Type(1)/sqrt(phi(i)); //lambda
        tmp_loglik = glmmtmb::dgenpois(yobs(i), s1, s2, true);
        SIMULATE{yobs(i)=glmmtmb::rgenpois(mu(i) / sqrt(phi(i)), Type(1) - Type(1)/sqrt(phi(i)));}
        break;
      case truncated_genpois_family:
        s1 = mu(i) / sqrt(phi(i)); //theta
        s2 = Type(1) - Type(1)/sqrt(phi(i)); //lambda
        tmp_loglik = zt_lik_nearzero(yobs(i),
		    glmmtmb::dgenpois(yobs(i), s1, s2, true) - log_nzprob(i));
        SIMULATE{yobs(i)=glmmtmb::rtruncated_genpois(mu(i) / sqrt(phi(i)), Type(1) - Type(1)/sqrt(phi(i)));}
        break;
      case compois_family:
        s1 = mu(i); //mean
        s2 = 1/phi(i); //nu
        tmp_loglik = dcompois2(yobs(i), s1, s2, true);
        SIMULATE{yobs(i)=rcompois2(mu(i), 1/phi(i));}
        break;
      case truncated_compois_family:
        s1 = mu(i); //mean
        s2 = 1/phi(i); //nu
        log_nzprob(i) = logspace_sub(Type(0), dcompois2(Type(0), s1, s2, true));
        tmp_loglik = zt_lik_nearzero(yobs(i),
			    dcompois2(yobs(i), s1, s2, true) - log_nzprob(i));
        SIMULATE{yobs(i)=glmmtmb::rtruncated_compois2(mu(i), 1/phi(i));}
        break;
      case tweedie_family:
        s1 = mu(i);  // mean
        s2 = phi(i); // phi
        s3 = invlogit(psi(0)) + Type(1); // p, 1<p<2
        tmp_loglik = dtweedie(yobs(i), s1, s2, s3, true);
        SIMULATE {
          yobs(i) = glmmtmb::rtweedie(s1, s2, s3);
        }
	break;
      case lognormal_family:
	// parameterized in terms of mean and SD on *data* scale, i.e.
	// mu = exp(logmu + logsd^2/2)
	// sd = sqrt((exp(logsd^2)-1)*exp(2*logmu + logsd^2)) = mu*sqrt(exp(logsd^2)-1)
	// 1+(sd/mu)^2 = exp(logsd^2)
	// logvar = log(1+(sd/mu)^2)
	// logsd = sqrt(logvar)
	// logmu = log(mu)-logvar/2
	// logvar via logspace_add() [log1p not compatible with CppAD]
        s1 = logspace_add(2*(log(phi(i))-log(mu(i))), Type(0)); // log-scale var
        s2 = log(mu(i)) - s1/2; //log-scale mean
        s3 = sqrt(s1);          //log-scale sd
	tmp_loglik = zt_lik_zero(yobs(i),
			 dnorm(log(yobs(i)), s2, s3, true) - log(yobs(i)));
	SIMULATE{
	  yobs(i) = exp(rnorm(s2, s3));
	}  // untested
        break;
      case t_family:
        s1 = (yobs(i) - mu(i))/phi(i);
	s2 = exp(psi(0));
	// since resid was scaled above, density needs to be divided by log(sd) = log(var)/2 = etadisp(i)/2
	tmp_loglik = dt(s1, s2, true) - etadisp(i);
	SIMULATE{
	  yobs(i) = mu(i)+phi(i)*rt(s2);
	}  // untested
	break;
      case bell_family:
	// unfortunately need to back-transform from mu to underlying theta via Lambert W ...

	// see https://stackoverflow.com/questions/92396/why-cant-variables-be-declared-in-a-switch-statement for {} 
	{
	  Type btheta;
	  btheta = glmmtmb::LambertW(mu(i));
	  tmp_loglik = glmmtmb::dbell(yobs(i), btheta, true);
	  SIMULATE{yobs(i) = glmmtmb::rbell(btheta);}
	  break;
	}
      default:
        error("Family not implemented!");
      } // End switch

      // Add zero inflation
      if(zi_flag){
        Type logit_pz = etazi(i) ;
        Type log_pz   = -logspace_add( Type(0) , -logit_pz );
        Type log_1mpz = -logspace_add( Type(0) ,  logit_pz );
        if(yobs(i) == Type(0)){
          // Was:
          //   tmp_loglik = log( pz(i) + (1.0 - pz(i)) * exp(tmp_loglik) );
          tmp_loglik = logspace_add( log_pz, log_1mpz + tmp_loglik );
        } else {
          // Was:
          //   tmp_loglik += log( 1.0 - pz(i) );
          tmp_loglik += log_1mpz ;
        }
        SIMULATE{yobs(i) = yobs(i)*rbinom(Type(1), Type(1)-pz(i));}
      }
      tmp_loglik *= weights(i);

      // Add up
      jnll -= keep(i) * tmp_loglik;
    } // if !is.na(obs)
    } // loop over observations

    // Add priors
    vector<Type> parvec;
    int np = prior_distrib.size();
    Type parval, logpriorval;
    int par_ind = 0; // parameter index
    for (int i = 0; i < np; i++) {
      switch(prior_whichpar[i]) {
      case beta_vprior: parvec = beta; break;
      case betazi_vprior: parvec = betazi; break;
      case betadisp_vprior: parvec = betadisp; break;
      case theta_vprior: parvec = theta; break;
      case thetazi_vprior: parvec = thetazi; break;
      case psi_vprior: parvec = psi; break;
      }

      // need an if-clause here for multivariate distrib (lkj/corr parameters
      // otherwise go element-by-element
      if (prior_distrib[i] == lkj_prior) {
	vector<Type> corpars = parvec.segment(prior_elstart[i] ,
					      prior_elend[i]-prior_elstart[i]+1);
	jnll -= glmmtmb::dlkj(corpars, prior_params[par_ind], true);
      } else {
	for (int j = prior_elstart[i]; j <= prior_elend[i]; j++) { // <= is on purpose here
	  if (j >= parvec.size()) {
	    // FIXME: should also check upstream ...
	    error("Bad prior index!");
	  };

	parval = parvec[j];
	switch(prior_distrib[i]) {
	case normal_prior:
	  s1 = prior_params[par_ind];           // mean
	  s2 = prior_params[par_ind+1];         // sd
	  logpriorval = dnorm(parval, s1, s2, true);
	break;
	case gamma_prior:
	  s1 = prior_params[par_ind+1];           // shape
	  s2 = prior_params[par_ind] / prior_params[par_ind+1];   // scale
	  logpriorval = dgamma(exp(parval), s1, s2, true);
	break;
	case t_prior:
	  s1 = prior_params[par_ind];           // mean
	  s2 = prior_params[par_ind+1];         // sd
	  s3 = prior_params[par_ind+2];         // df
	  parval = (parval - s1)/s2;   // scale value
	  // see note at t_family about adjusting density for scaling
	  logpriorval = dt(parval, s3, true) - log(s2);
	  break;
	case cauchy_prior:
	  s1 = prior_params[par_ind];           // loc
	  s2 = prior_params[par_ind+1];         // scale
	  logpriorval = glmmtmb::dcauchy(parval, s1, s2, true);
	break;
	default:
	  error("Prior distribution not implemented!");
	}

	jnll -= logpriorval;
	
	} // loop over elements
      }
      // step forward in prior-parameter vector
      par_ind += prior_npar[i];
      
    } // loop over priors
    
    

  // Report / ADreport / Simulate Report
  vector<matrix<Type> > corr(terms.size());
  vector<vector<Type> > sd(terms.size());
  for(int i=0; i<terms.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(terms(i).blockNumTheta > 0){
      corr(i) = terms(i).corr;
      sd(i) = terms(i).sd;
    }
  }
  vector<matrix<Type> > corrzi(termszi.size());
  vector<vector<Type> > sdzi(termszi.size());
  for(int i=0; i<termszi.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(termszi(i).blockNumTheta > 0){
      corrzi(i) = termszi(i).corr;
      sdzi(i) = termszi(i).sd;
    }
  }
  vector<matrix<Type> > corrdisp(termsdisp.size());
  vector<vector<Type> > sddisp(termsdisp.size());
  for(int i=0; i<termsdisp.size(); i++){
  	// NOTE: Dummy terms reported as empty
  	if(termsdisp(i).blockNumTheta > 0){
  		corrdisp(i) = termsdisp(i).corr;
  		sddisp(i) = termsdisp(i).sd;
  	}
  }
  vector<matrix<Type> > fact_load(terms.size());
  for(int i=0; i<terms.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(terms(i).blockNumTheta > 0){
      fact_load(i) = terms(i).fact_load;
    }
  }

  REPORT(corr);
  REPORT(sd);
  REPORT(corrzi);
  REPORT(sdzi);
  REPORT(corrdisp);
  REPORT(sddisp);
  REPORT(fact_load);
  REPORT(b);
  REPORT(bzi);
  SIMULATE {
    REPORT(yobs);
  }
  // For predict
  if(ziPredictCode == disp_zipredictcode) {
    // predict dispersion
    // zi irrelevant; just reusing variable
    switch(family){
    case Gamma_family:
      mu = 1/sqrt(phi);
      break;
    default:
      mu = phi;
    }
  } else {
    if (trunc_Family(family)) {
      // convert from mean of *un-truncated* to mean of *truncated* distribution
      mu /= exp(log_nzprob);
    }
    if (zi_flag) {
      switch(ziPredictCode){
      case corrected_zipredictcode:
	mu *= (Type(1) - pz); // Account for zi in prediction
	break;
      case uncorrected_zipredictcode:
	//mu = mu; // Predict mean of 'family' //commented out for clang 7.0.0. with no effect
	break;
      case prob_zipredictcode:
	mu = pz;     // Predicted zi probability
	eta = etazi; // want to return linear pred for zi
	break;
      default:
	error("Invalid 'ziPredictCode'");
      }
    }}

  whichPredict -= 1; // R-index -> C-index
  vector<Type> mu_predict = mu(whichPredict);
  vector<Type> eta_predict = eta(whichPredict);

  DATA_FACTOR(aggregate);
  if (aggregate.size() > 0) {
    if (aggregate.size() != mu_predict.size())
      Rf_error("'aggregate' wrong size");
    vector<Type> tmp(NLEVELS(aggregate));
    tmp.setZero();
    for (int i=0; i<aggregate.size(); i++) {
      tmp[aggregate[i]] += mu_predict[i];
    }
    mu_predict = tmp;
    for (int i=0; i<tmp.size(); i++) {
      tmp[i] = linkfun(tmp[i], link);
    }
    eta_predict = tmp;
  }

  REPORT(mu_predict);
  REPORT(eta_predict);
  // ADREPORT expensive for long vectors - only needed by predict() method
  if (doPredict==1) {
	  ADREPORT(mu_predict);
  } else if (doPredict == 2) {
	  ADREPORT(eta_predict);
  } else if (doPredict == 3) {
           ADREPORT(b);
	   ADREPORT(bzi);
	   ADREPORT(bdisp);
  }
  return jnll;
}
