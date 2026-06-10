#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <RcppGSL.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include<trunclst.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

#define crossprod(x) symmatu(x.t() * x)
#define tcrossprod(x) symmatu(x * x.t())
#define ssq(x) dot(x,x)

// [[Rcpp::export]]
double log_sum_exp(const arma::vec &x)
{
  arma::uword max_ind=x.index_max();
  double maxVal= x(max_ind);
  
  double sum_exp=0.0;
  
  for (unsigned  i = 0; i < x.n_elem ; ++i){
    if(i !=max_ind)
      sum_exp += exp(  (x(i) - maxVal));
  }
  return log1p(sum_exp)+maxVal ;
}

// [[Rcpp::export]]
double surv_fn_lognorm(const double x, const unsigned nu,const  double ss_survtime,const  double survtime,
                       const double df0,const double a0,const double mu0,const double beta0,const unsigned n){//nu=1 implies failure time
  ///don't need x_aug argument here!!
  double df_post=df0+ n , alpha_post=a0+ n/2.0;
  double tmp=(   n *df0 )/df_post, ss_j, survtime_j, mean_j, mu_post, beta_post, sigma_post;
  
  if(n){
    survtime_j = survtime ;
    mean_j= survtime_j/ n  ;
    ss_j= ss_survtime - n *gsl_pow_2(mean_j);
  } else ss_j=  survtime_j = mean_j= 0.0;
  
  
  mu_post = (df0*mu0 + survtime_j) / df_post;
  beta_post=beta0+  (ss_j   + tmp* gsl_pow_2(mean_j - mu0) )/2.0; //ss_j differently defined from the rest of the code
  sigma_post =sqrt( beta_post * (df_post+1)/(df_post *alpha_post) );
  
  // Rcpp::Rcout<<"nu= " << nu << std::endl;
  // Rcpp::Rcout<<"mu0= "<<mu0 << " beta0= "<<beta0<<endl;
  // Rcpp::Rcout<<"mean_j= "<<mean_j<<" survtime_j= "<<survtime_j << endl;
  // Rcpp::Rcout << "ss_survtime = " << ss_survtime << endl;
  // Rcpp::Rcout<< " ss_j= " << ss_j<< " gsl_pow_2(mean_j - mu0)= " <<gsl_pow_2(mean_j - mu0)<<" n="<<n << " tmp="<<tmp<< " df_post="<<df_post<< " alpha_post="<<alpha_post<<
  // " beta_post= "<<beta_post<<" sigma_post= "<<sigma_post<<endl;
  // Rcpp::Rcout << "ss_survtime= " << ss_survtime<<endl;
  // Rcpp::Rcout<<"survtime_j= "<<survtime_j<<" ss_j= "<<ss_j<<  " df = "<< 2*df_post<< " mu_post = "<< mu_post <<" sigma_post = "<<sigma_post<< endl;
  // Rcpp::Rcout<<"x= "<<x <<endl;
  
  return nu ?  ( d_lst( x,  2.0*alpha_post, mu_post, sigma_post, 1) ) : ( p_lst(x, 2.0*alpha_post, mu_post, sigma_post, 0, 1 ) );
}

// [[Rcpp::export]]
double post_t_dens(const double x, const  double ss_survtime,const  double survtime,
                   const double df0,const double a0,const double mu0,const double beta0,const unsigned n){
  ///ss_survtime is ssq and survtime is just the sum
  double df_post=df0+ n , alpha_post=a0+ n/2.0;
  double tmp=(   n *df0 )/df_post, ss_j, survtime_j, mean_j, mu_post, beta_post, sigma_post;
  
  if(n){
    survtime_j = survtime ;
    mean_j= survtime_j/ n  ;
    ss_j= ss_survtime -n *gsl_pow_2(mean_j);
  } else ss_j=  survtime_j = mean_j= 0.0;
  
  
  mu_post = (df0*mu0 + survtime_j) / df_post;
  
  beta_post=beta0+  (ss_j   + tmp* gsl_pow_2(mean_j - mu0) )/2.0; //ss_j differently defined from the rest of the code
  sigma_post =sqrt( beta_post * (df_post+1)/(df_post *alpha_post) );
  
  // Rcpp::Rcout << "2 * alpha_post " << alpha_post << std::endl;
  // Rcpp::Rcout << "survtime_j " << survtime << std::endl;
  // Rcpp::Rcout << "mean_j " << mean_j << std::endl;
  // Rcpp::Rcout << "n " << n << std::endl;
  // Rcpp::Rcout << "ss_survtime " << ss_survtime << std::endl;
  // Rcpp::Rcout << "ss_j " << ss_j << std::endl;
  // Rcpp::Rcout << "mu_post " << mu_post << std::endl;
  // Rcpp::Rcout << "beta_post " << beta_post << std::endl;
  // Rcpp::Rcout << "sigma_post " << sigma_post << std::endl;
  
  // gsl_sf_lnpoch(a,b)=log (gamma(a+b)/gamma(a)  )
  
  
  /*double df_final=2.0*alpha_post;
   double denss=gsl_sf_lnpoch(alpha_post,.5) - (M_LNPI/2 + log( sqrt(df_final)* sigma_post ) +(alpha_post+.5)* log1p( gsl_pow_2( (x-mu_post)/sigma_post )/df_final  ) );
   Rcpp::Rcout<<"mu_post="<<mu_post<<" ss_j="<<ss_j <<  " manual dens="<<denss<<" RcppDIST dens="<<d_lst( x,  2.0*alpha_post, mu_post, sigma_post, 1)<<endl;*/
  
  return d_lst( x,  2.0*alpha_post, mu_post, sigma_post, 1);
}

arma::vec sim_lognorm_params( const  double ss_j,const  double survtime_j,
                              const double df0,const double a0,const double mu0,const double beta0,const unsigned n){
  double df_post=df0+ n, alpha_post=a0+ n/2.0;
  double tmp=( n*df0 )/df_post;
  double mean_j=  n ? (survtime_j/ n) : 0.0;
  double mu_post = (df0*mu0 + survtime_j) / df_post;
  double beta_post=beta0+  (ss_j -  ((double) n)  * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
  double sigma_post =sqrt( beta_post  / (alpha_post * df_post ) );
  
  arma::vec ret(2);
  // Rcpp::Rcout<<"Flag lognormal_param 0 alpha_post= "<<alpha_post<<" beta_post= "<< beta_post<<" data beta= "<<(ss_j -  ((double) n)  * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0<<endl;
  ret(0)=r_lst(2*alpha_post,  mu_post, sigma_post );
  ret(1)= 1/ randg( distr_param(alpha_post, 1/beta_post) );
  return ret;
}

/******* updated functions required for HMC on main location scale controlling parameters**********/
double logpi_lognorm_extend(const arma::mat &nj_val,
                            const arma::mat &survtime, const arma::mat &ss_survtime,
                            const arma::vec &params,
                            const double a0, const double df0,
                            const double mu_m, const double mu_v,
                            const double b_m, const double b_v,
                            const unsigned t){
  //mu_v and b_v are prior variances, not stds!!!
  double mu0=params(0),  beta0=exp(params(1)), log_b0=params(1);
  double sum= - (gsl_pow_2(mu0-mu_m)/mu_v + gsl_pow_2( log_b0 - b_m) / b_v) /2; //normal prior
  //log_normpdf(mu0,mu_m, sqrt(mu_v)) + log_normpdf(log_b0 , b_m, sqrt(b_v));
  
  for (unsigned j = 0; j < nj_val.n_rows; ++j)
  {
    if (nj_val(j, t))  // only use this treatment’s column
    {
      double df_post = df0 + nj_val(j, t);
      double alpha_post = a0 + nj_val(j, t) / 2.0;
      double tmp = (nj_val(j, t) * df0) / df_post;
      double ss_j = ss_survtime(j, t);
      double survtime_j = survtime(j, t);
      double mean_j = survtime_j / nj_val(j, t);
      double beta_post = beta0 + (ss_j - nj_val(j, t) * gsl_pow_2(mean_j)
                                    + tmp * gsl_pow_2(mean_j - mu0)) / 2.0;
      sum += (a0 * log_b0 - alpha_post * log(beta_post));
    }
  }
  return sum;
}


arma::vec delpi_lognorm_extend(const arma::mat &nj_val,
                               const arma::mat &survtime, const arma::mat &ss_survtime,
                               const arma::vec &params,
                               const double a0, const double df0,
                               const double mu_m, const double mu_v,
                               const double b_m, const double b_v,
                               const unsigned t){
  //mu_v and b_v are prior variances, not sds!!!
  double mu0=params(0),  beta0=exp(params(1));
  double del_m=-(mu0-mu_m)/mu_v, del_b= -( params(1) - b_m) / b_v;//derivative of normal prior
  
  for (unsigned j = 0; j < nj_val.n_rows; ++j)
  {
    if (nj_val(j, t))
    {
      double df_post = df0 + nj_val(j, t);
      double alpha_post = a0 + nj_val(j, t) / 2.0;
      double tmp = (nj_val(j, t) * df0) / df_post;
      double ss_j = ss_survtime(j, t);
      double survtime_j = survtime(j, t);
      double mean_j = survtime_j / nj_val(j, t);
      double beta_post = beta0 + (ss_j - nj_val(j, t) * gsl_pow_2(mean_j)
                                    + tmp * gsl_pow_2(mean_j - mu0)) / 2.0;
      double common_part = alpha_post / beta_post;
      del_m += (common_part * tmp * (mean_j - mu0));
      del_b += (a0 - common_part * beta0);
    }
  }
  arma::vec ret(2);
  ret(0)=del_m; ret(1)= del_b;
  return ret;
}

void leapfrog_lognorm_hyper_extend(const unsigned nstep,
                                   const double delta, arma::vec &v_old,
                                   arma::vec &p_lam, arma::vec &params,
                                   const arma::mat &nj_val,
                                   const arma::mat &survtime, const arma::mat &ss_survtime,
                                   const double a0, const double df0,
                                   const double mu_m, const double mu_v,
                                   const double b_m, const double b_v,
                                   const unsigned t){
  for(unsigned i=0;i<nstep;++i){
    // Rcpp::Rcout<<"flag -1 leap i="<<i<<endl;
    params+=(delta)*(p_lam -(delta/2)*v_old);
    // Rcpp::Rcout<<"flag 0 leap i="<<i<<endl;
    params(0)=std::clamp(params(0), -1e3, 1e5);
    params(1)=std::clamp(params(1), -1e1, 1e1);
    
    arma::vec v_new = -delpi_lognorm_extend(nj_val, survtime, ss_survtime,
                                            params, a0, df0, mu_m, mu_v, b_m, b_v, t);
    // Rcpp::Rcout<<"flag 1.5 leap i="<<i<<endl;
    p_lam-=(delta/2)* ( v_old+v_new);
    v_old=v_new;
  }
}


void update_lognorm_hyper_extend(const arma::mat &nj_val,
                                 const arma::mat &survtime, const arma::mat &ss_survtime,
                                 const double a0,const  double df0,
                                 const double mu_m, const  double mu_v,
                                 const double b_m, const double b_v,
                                 const arma::vec &del_range, const unsigned nleapfrog,
                                 arma::vec &params, unsigned &acceptance,
                                 const unsigned t){
  double ll_old = -logpi_lognorm_extend(nj_val, survtime, ss_survtime,
                                        params, a0, df0, mu_m, mu_v, b_m, b_v, t);
  arma::vec v_old = -delpi_lognorm_extend(nj_val, survtime, ss_survtime,
                                          params, a0, df0, mu_m, mu_v, b_m, b_v, t);
  
  arma::vec p_lam(2, fill::randn);
  double kin_energy= dot(p_lam, p_lam);
  
  arma::vec params_new=params, v_new=v_old;
  
  unsigned pois_draw=(unsigned) R::rpois(nleapfrog);
  unsigned nstep=GSL_MAX_INT(1,pois_draw);
  double delta= arma::randu(distr_param(del_range(0),del_range(1)));
  // Rcpp::Rcout<<"Flag X params -.5!!"<<endl;
  leapfrog_lognorm_hyper_extend(nstep, delta, v_new, p_lam, params_new,
                                nj_val, survtime, ss_survtime,
                                a0, df0, mu_m, mu_v, b_m, b_v, t);
  // params_new.print("params_new");
  // Rcpp::Rcout<<"Flag X params 0!!"<<endl;
  double ll_new = -logpi_lognorm_extend(nj_val, survtime, ss_survtime,
                                        params_new, a0, df0, mu_m, mu_v, b_m, b_v, t);
  
  //get H_ll_new
  double H_new= ll_new+  ssq(p_lam) /2, H_old=ll_old+kin_energy/2;
  
  if(log(randu())< -(H_new-H_old) ){
    params=params_new;
    /*ll_old=ll_new;
     v_old=v_new;*/
    ++acceptance;
  }
  // Rcpp::Rcout << "Acceptance: " << acceptance << std::endl;
  /*Rcpp::Rcout<<"i= "<<i<< "Acceptance rate="<<(((double)acceptance)/ ((double) i))<<" alpha= "<<exp(l_alpha)<<endl;
   Rcpp::Rcout<<"Acceptance rate="<<(acceptance/n_mc)<<endl;
   return exp(alpha_vec);*/
}
/*****************************************************************/


/*******functions required for HMC on Dirichlet mixture parameters**********/
double ll_alpha(const double l_alpha, const unsigned K, const arma::vec &nj_vec,
                const double mu_alp, const double sig_alp){
  const double alpha=exp(l_alpha);
  const unsigned Kn=nj_vec.n_elem;
  if(K<Kn)
    Rcpp::stop("Error in ll_alpha: # Mixtures < # Occupied clusters! ");
  
  const double al_by_k= alpha/K;
  
  double ret= lgamma(alpha) - lgamma(alpha+ accu(nj_vec)) + accu(lgamma(nj_vec+al_by_k))
    - Kn * lgamma(al_by_k) + log_normpdf(l_alpha , mu_alp, sig_alp );
  
  return ret;
}

double del_ll_alpha(const double l_alpha, const unsigned K, const arma::vec &nj_vec,
                    const double mu_alp, const double sig_alp){
  const double alpha=exp(l_alpha);
  const unsigned Kn=nj_vec.n_elem;
  if(K<Kn)
    Rcpp::stop("Error in ll_alpha: # Mixtures < # Occupied clusters! ");
  
  const double al_by_k= alpha/((double) K);
  
  vec tmp=nj_vec+ al_by_k;
  tmp.transform( [](double val) { return ( gsl_sf_psi(val) ); } );
  
  // Rcpp::Rcout<<"alpha="<<alpha<<" al_by_k="<<al_by_k<<" alpha+ accu(nj_vec)="<<alpha+ accu(nj_vec) <<endl;
  double ret= (gsl_sf_psi(alpha) - gsl_sf_psi(alpha+ accu(nj_vec)) + (accu(tmp) -  Kn*gsl_sf_psi(al_by_k) )/ ((double) K) )*alpha
  - (l_alpha - mu_alp)/gsl_pow_2(sig_alp) ;
  
  return ret;
}

void leapfrog_dir_alpha(const unsigned nstep,const double delta, double &v_old,  double &p_lam, double &l_alpha,
                        const unsigned K, const arma::vec &nj_vec, const double mu_alp, const double sig_alp){
  for(unsigned i=0;i<nstep;++i){
    l_alpha+=(delta)*(p_lam -(delta/2)*v_old);
    l_alpha = std::clamp(l_alpha, -1e1, 1e1);
    
    double v_new=-del_ll_alpha(l_alpha, K, nj_vec, mu_alp, sig_alp);
    p_lam-=(delta/2)* ( v_old+v_new);
    v_old=v_new;
  }
}

void update_alpha(const unsigned K, const arma::vec &nj_vec, const double mu_alp, const double sig_alp,
                  const arma::vec &del_range_alp, const unsigned nleapfrog_alp, double &l_alpha, unsigned &acceptance){
  /*double mu_alp, sig_alp;
   sig_alp=log1p(ps_hyper(1) /gsl_pow_2(ps_hyper(0)));mu_alp=log(ps_hyper(0))-sig_alp/2;sig_alp=sqrt(sig_alp);*/
  
  double ll_old=-ll_alpha(l_alpha, K, nj_vec, mu_alp, sig_alp), v_old=-del_ll_alpha(l_alpha, K, nj_vec, mu_alp, sig_alp);
  
  double p_lam=randn();
  double kin_energy= gsl_pow_2(p_lam);
  
  double l_alpha_new=l_alpha, v_new=v_old;
  
  unsigned pois_draw=(unsigned) R::rpois(nleapfrog_alp); //randomly generating no. of leapfrog steps
  unsigned nstep=GSL_MAX_INT(1,pois_draw);
  // nstep=GSL_MIN_INT(nstep,leapmax);
  double delta= arma::randu(distr_param(del_range_alp(0),del_range_alp(1)));//randomly generating \delta t
  // Rcpp::Rcout<<"Flag X params -.5!!"<<endl;
  leapfrog_dir_alpha(nstep, delta, v_new,  p_lam, l_alpha_new, K, nj_vec, mu_alp, sig_alp);
  // params_new.print("params_new");
  // Rcpp::Rcout<<"Flag X params 0!!"<<endl;
  double ll_new=- ll_alpha(l_alpha_new, K, nj_vec, mu_alp, sig_alp); //log-likelihood at the new proposed value
  
  //get H_ll_new
  double H_new= ll_new+  gsl_pow_2(p_lam) /2, H_old=ll_old+kin_energy/2;
  
  if(log(randu())< -(H_new-H_old) ){
    l_alpha=l_alpha_new;
    /*ll_old=ll_new;
     v_old=v_new;*/
    ++acceptance;
  }
  /*Rcpp::Rcout<<"i= "<<i<< "Acceptance rate="<<(((double)acceptance)/ ((double) i))<<" alpha= "<<exp(l_alpha)<<endl;
   Rcpp::Rcout<<"Acceptance rate="<<(acceptance/n_mc)<<endl;
   return exp(alpha_vec);*/
}
/***************************************************/

/*****************************************************************************/
// Build a 3D field of level-count vectors from the Stage-1 noccu_old list
// - noccu_old: R list of length n_mc_old; each item is a (nmix × k) list (dim attr)
// - ncat:      length-k vector with #levels for each categorical covariate
// Returns: F(j,c,s) = uvec(ncat(c)) of counts at draw s
// Build field<arma::uvec>(nmix, k, n_mc_old) from flat noccu_old
arma::field<arma::uvec> build_noccu_field(const Rcpp::List& noccu_old,
                                          unsigned nmix, unsigned k,
                                          const arma::uvec& ncat) {
  const unsigned n_mc_old = noccu_old.size();
  field<arma::uvec> F(nmix, k, n_mc_old);
  
  for (unsigned s = 0; s < n_mc_old; ++s) {
    Rcpp::List Ls = noccu_old[s];                 // length should be nmix * k
    if ((unsigned)Ls.size() != nmix * k) {
      Rcpp::stop("noccu_old[[%d]] has %d elements; expected %d (= nmix*k).",
                 (int)s, (int)Ls.size(), (int)(nmix * k));
    }
    
    for (unsigned c = 0; c < k; ++c) {
      const unsigned L = ncat(c);
      for (unsigned j = 0; j < nmix; ++j) {
        const unsigned idx = j + nmix * c;        // <-- flat index
        
        // Accept INT or REAL; coerce to integer vector of counts
        Rcpp::RObject elem = Ls[idx];
        Rcpp::IntegerVector v;
        switch (TYPEOF(elem)) {
        case INTSXP:  v = Rcpp::as<Rcpp::IntegerVector>(elem); break;
        case REALSXP: v = Rcpp::as<Rcpp::IntegerVector>(Rcpp::as<Rcpp::NumericVector>(elem)); break;
        default:
          Rcpp::stop("noccu_old[[%d]][idx=%d] has unsupported SEXP type %d.",
                     (int)s, (int)idx, (int)TYPEOF(elem));
        }
        
        arma::uvec out(L, arma::fill::zeros);
        const unsigned m = std::min<unsigned>(v.size(), L);
        for (unsigned u = 0; u < m; ++u) {
          int x = v[u];
          out(u) = (Rcpp::IntegerVector::is_na(x) || x < 0) ? 0u : (unsigned)x;
        }
        F(j, c, s) = std::move(out);
      }
    }
  }
  return F;
}

/*****************************************************************************/

/*****************************************************************************/
// for initial cluster allocations in Stage 2


// [[Rcpp::export]]
Rcpp::List stage2_cov_loglik_precompute(
    const Rcpp::IntegerVector& init_ids,   // length R, 0-based draw indices
    const arma::mat& eta_cat,              // [N x k_cat], 0-based codes; may contain NA (but we won't scan for NA)
    const arma::mat& eta_cont,             // [N x k_cont], scaled; may contain NA (but we won't scan for NA)
    Rcpp::List non_na_obs1,                // length N; each is uvec of non-missing cat positions (0-based)
    Rcpp::List non_na_obs1_cont,           // length N; each is uvec of non-missing cont positions (0-based)
    const arma::cube& nobs_cube,           // [K x k_cat x M]
    const arma::cube& nj_x_cube,           // [K x k_cont x M]
    const arma::cube& sum_x_cube,          // [K x k_cont x M]
    const arma::cube& ss_x_cube,           // [K x k_cont x M]
    const Rcpp::List& noccu_old,           // length M; each is flat list of length K*k_cat of level-count vectors
    const arma::uvec& ncat,                // length k_cat
    const double df_x,                     // continuous predictive hyper
    const double alpha_x,
    const double mu_x,
    const double beta_x
){
  // ---- dimensions (all 0-based logic) ----
  const unsigned N      = eta_cat.n_rows;               // # current cohort rows for init
  const unsigned k_cat  = eta_cat.n_cols;               // # categorical covariates
  const unsigned k_cont = eta_cont.n_cols;              // # continuous covariates
  const unsigned K      = nobs_cube.n_rows;             // # clusters
  const unsigned M      = nobs_cube.n_slices;           // # Stage-1 draws
  const unsigned R      = init_ids.size();              // # init draws used
  
  // sanity checks
  if ( (k_cat > 0 && (nobs_cube.n_cols != k_cat)) ||
       (k_cont> 0 && (nj_x_cube.n_cols != k_cont || sum_x_cube.n_cols != k_cont || ss_x_cube.n_cols != k_cont)) ||
       (nobs_cube.n_rows != K) || (nj_x_cube.n_rows != K) ||
       (sum_x_cube.n_rows != K) || (ss_x_cube.n_rows != K) )
  {
    Rcpp::stop("Dimension mismatch among K/k_cat/k_cont across input cubes.");
  }
  
  // ---- unpack non-missing index lists (0-based) ----
  arma::field<arma::uvec> non_na_obs(N), non_na_obs_cont(N);
  for (unsigned i = 0; i < N; ++i){
    non_na_obs(i)      = Rcpp::as<arma::uvec>(non_na_obs1[i]);
    non_na_obs_cont(i) = Rcpp::as<arma::uvec>(non_na_obs1_cont[i]);
  }
  
  // ---- build categorical level-count field: F(j, c, m) -> uvec(ncat(c)) ----
  arma::field<arma::uvec> F;
  if (k_cat > 0) {
    if (ncat.n_elem != k_cat) {
      Rcpp::stop("ncat length (%d) must equal k_cat (%d).", (int)ncat.n_elem, (int)k_cat);
    }
    F = build_noccu_field(noccu_old, K, k_cat, ncat); // dims (K, k_cat, M)
  }
  
  // ---- output cubes: [R x K x N] to match R's cov_loglik_*[rr, j, i] ----
  arma::cube cov_loglik_cont(R, K, N, arma::fill::zeros);
  arma::cube cov_loglik_cat (R, K, N, arma::fill::zeros);
  
  // ---- main loops (all indices are 0-based) ----
  for (unsigned rr = 0; rr < R; ++rr) {
    int m0 = init_ids[rr];     // 0-based Stage-1 draw index
    if (m0 < 0 || (unsigned)m0 >= M) {
      Rcpp::stop("init_ids[%d]=%d out of range [0, %d).", (int)rr, (int)m0, (int)M);
    }
    
    // continuous block
    if (k_cont > 0) {
      for (unsigned j = 0; j < K; ++j) {
        for (unsigned i = 0; i < N; ++i) {
          double li = 0.0;
          const arma::uvec& nz = non_na_obs_cont(i);   // non-missing positions in eta_cont(i, .)
          for (auto c : nz) {
            if (c >= k_cont) continue;                 // guard (shouldn't happen)
            double x_i = eta_cont(i, c);
            // if x_i is NA/Inf, skip (defensive)
            if (!std::isfinite(x_i)) continue;
            
            double nC  = nj_x_cube(j, c, (unsigned)m0);
            double sC  = sum_x_cube(j, c, (unsigned)m0);
            double ssC = ss_x_cube(j, c, (unsigned)m0);
            
            li += post_t_dens(
              x_i,
              ssC,
              sC,
              df_x,
              alpha_x,
              mu_x,
              beta_x,
              (unsigned)std::max(0.0, nC)
            );
          }
          cov_loglik_cont(rr, j, i) = li;
        }
      }
    }
    
    // categorical block (Dirichlet-α=1 predictive)
    if (k_cat > 0) {
      for (unsigned j = 0; j < K; ++j) {
        for (unsigned i = 0; i < N; ++i) {
          double lj = 0.0;
          const arma::uvec& nz = non_na_obs(i);        // non-missing positions in eta_cat(i, .)
          for (auto v : nz) {
            if (v >= k_cat) continue;                  // guard
            double val = eta_cat(i, v);
            if (!std::isfinite(val) || val < 0.0) continue;  // defensive
            unsigned lev = static_cast<unsigned>(val);        // 0-based level code
            
            double nobsJ = nobs_cube(j, v, (unsigned)m0);     // total obs for (j,v) at draw m0
            const arma::uvec& cnt = F(j, v, (unsigned)m0);    // level counts vector
            const unsigned Lv = cnt.n_elem;
            
            double nlev = (lev < Lv) ? (double)cnt(lev) : 0.0;
            lj += ( std::log(nlev + 1.0) - std::log(nobsJ + (double)Lv) );
          }
          cov_loglik_cat(rr, j, i) = lj;
        }
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("cov_loglik_cont") = cov_loglik_cont,
    Rcpp::Named("cov_loglik_cat")  = cov_loglik_cat
  );
}


// [[Rcpp::export]]
arma::imat stage2_init_alloc_cov_only(
    const arma::cube& cov_loglik_cat,   // [R x K x N] log-lik from categorical X (per draw, cluster, subject)
    const arma::cube& cov_loglik_cont,  // [R x K x N] log-lik from continuous  X (per draw, cluster, subject)
    const arma::uvec& init_ids,         // [R] 0-based draw indices (kept for interface parity; must match rows)
    const unsigned    nmix,             // K
    const unsigned    rng_seed = 500    // GSL seed
){
  // ----- dimension checks -----
  const arma::uword R = cov_loglik_cat.n_rows;
  const arma::uword K = cov_loglik_cat.n_cols;
  const arma::uword N = cov_loglik_cat.n_slices;
  
  if (cov_loglik_cont.n_rows   != R ||
      cov_loglik_cont.n_cols   != K ||
      cov_loglik_cont.n_slices != N) {
    Rcpp::stop("cov_loglik_cat and cov_loglik_cont must have identical dims [R x K x N].");
  }
  if (K != nmix) {
    Rcpp::stop("nmix (%u) != K (%u).", (unsigned)nmix, (unsigned)K);
  }
  if (init_ids.n_elem != R) {
    Rcpp::stop("init_ids length (%u) must equal R (%u), the # of rows in cov_loglik_*.", 
               (unsigned)init_ids.n_elem, (unsigned)R);
  }
  
  // ----- RNG -----
  const gsl_rng_type * T_rng;
  gsl_rng * r;
  gsl_rng_env_setup();
  T_rng = gsl_rng_default;
  r = gsl_rng_alloc(T_rng);
  gsl_rng_set(r, rng_seed);
  
  arma::imat Z_init(R, N, arma::fill::zeros);   // 0-based labels
  arma::vec  logp(K, arma::fill::zeros);
  arma::vec  pr(K,    arma::fill::zeros);
  
  // ----- main loops -----
  for (arma::uword i = 0; i < N; ++i) {
    for (arma::uword rr = 0; rr < R; ++rr) {
      // score per cluster (covariates only)
      for (arma::uword j = 0; j < K; ++j) {
        const double lxc = cov_loglik_cat(rr, j, i);
        const double lcc = cov_loglik_cont(rr, j, i);
        logp(j) = lxc + lcc;     // no Y term
      }
      
      // stable normalize
      const double lp_max = logp.max();
      pr = arma::exp(logp - lp_max);
      const double denom = arma::accu(pr);
      if (denom <= 0.0 || !arma::is_finite(denom)) {
        pr.fill(1.0 / double(K));   // uniform fallback if degenerate
      } else {
        pr /= denom;
      }
      
      // sample one-hot via GSL multinomial, then argmax
      arma::uvec one_hot(K, arma::fill::zeros);
      gsl_ran_multinomial(r, (size_t)K, 1, pr.memptr(), one_hot.memptr());
      Z_init(rr, i) = one_hot.index_max();   // 0-based cluster label
    }
  }
  
  gsl_rng_free(r);
  return Z_init;
}

/*****************************************************************************/

/************************** compatibility ************************************/

double logpi_lognorm(const arma::vec &nj_val1, const arma::vec &nj_val2, 
                     const arma::vec &survtime1, const arma::vec &ss_survtime1,
                     const arma::vec &survtime2, const arma::vec &ss_survtime2,
                     const arma::vec &params,
                     const double a0, const double df0, 
                     const double mu_m, const double mu_v,
                     const double b_m, const double b_v){
  //mu_v and b_v ae prior variances, not sds!!!
  double mu0=params(0),  beta0=exp(params(1)), log_b0=params(1);
  double sum= - (gsl_pow_2(mu0-mu_m)/mu_v + gsl_pow_2( log_b0 - b_m) / b_v) /2; //normal prior
  //log_normpdf(mu0,mu_m, sqrt(mu_v)) + log_normpdf(log_b0 , b_m, sqrt(b_v));
  
  
  // Rcpp::Rcout<<"sum= "<<sum<<endl;
  for(unsigned j=0; j< nj_val2.n_elem; ++j){//Iterate over the clusters
    ///notations mostly follow Wikipedia conjugate prior NIG
    if(nj_val2(j)){
      double df_post=df0+ nj_val2(j), alpha_post=a0+ nj_val2(j)/2.0; 
      double tmp=( nj_val2(j)*df0 )/df_post;
      double ss_j= ss_survtime2(j) , survtime_j = survtime2(j) ;
      double mean_j= survtime_j/nj_val2(j) ;
      double beta_post=beta0+  (ss_j -  nj_val2(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
      sum+= (a0 *log_b0 -alpha_post *  log( beta_post));
      /*Rcpp::Rcout<<"df_post= "<<df_post<<" tmp="<<tmp<<" beta_post= "<<beta_post<<" ss_j= "<<ss_j<<" nj_val2(j) * gsl_pow_2(mean_j )="<<nj_val2(j) * gsl_pow_2(mean_j )<<  endl;
       Rcpp::Rcout<<"density at j= "<< j<<"is "<<(mu0 *log_b0 -alpha_post *  log( beta_post)) <<endl;*/
      
      if(nj_val1(j)){
        df_post=df0+ nj_val1(j); alpha_post=a0+ nj_val1(j) /2.0;
        tmp=( nj_val1(j) *df0 )/df_post;
        ss_j= ss_survtime1(j), survtime_j = survtime1(j);
        mean_j= survtime_j/nj_val1(j);
        beta_post=beta0+  (ss_j -  nj_val1(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
        sum+= (a0 *log_b0 -alpha_post *  log( beta_post));
      }
    }
  }
  return sum;
}

arma::vec delpi_lognorm(const arma::vec &nj_val1, const arma::vec &nj_val2, 
                        const arma::vec &survtime1, const arma::vec &ss_survtime1,
                        const arma::vec &survtime2, const arma::vec &ss_survtime2,
                        const arma::vec &params,
                        double a0, double df0, 
                        double mu_m, double mu_v,
                        double b_m, double b_v){
  //mu_v and b_v are prior variances, not sds!!!
  double mu0=params(0),  beta0=exp(params(1));
  double del_m=-(mu0-mu_m)/mu_v, del_b= -( params(1) - b_m) / b_v;//derivative of normal prior
  
  for(unsigned j=0; j<nj_val2.n_elem; ++j){//Iterate over the clusters
    ///notations mostly follow Wikipedia conjugate prior NIG
    ////UPDATE for the experimental arm
    if(nj_val2(j)){
      double df_post=df0+ nj_val2(j), alpha_post=a0+ nj_val2(j)/2.0;
      double tmp=( nj_val2(j)*df0 )/df_post;
      double ss_j= ss_survtime2(j) , survtime_j = survtime2(j) ;
      double mean_j= survtime_j/nj_val2(j) ;
      double beta_post=beta0+  (ss_j -   nj_val2(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
      double common_part=alpha_post / beta_post ;
      del_m+= (common_part * tmp*(mean_j-mu0));
      del_b+= (a0 -common_part*beta0) ;
      
      if(nj_val1(j)){
        df_post=df0+ nj_val1(j); alpha_post=a0+ nj_val1(j) /2.0;
        tmp=( nj_val1(j) *df0 )/df_post;
        ss_j= ss_survtime1(j), survtime_j = survtime1(j);
        mean_j= survtime_j/nj_val1(j);
        beta_post=beta0+  (ss_j -  nj_val1(j) * gsl_pow_2(mean_j )  + tmp* gsl_pow_2(mean_j - mu0) )/2.0;
        common_part=alpha_post / beta_post ;
        del_m+= (common_part * tmp*(mean_j-mu0));
        del_b+= (a0 -common_part*beta0) ;
      }
    }
  }
  arma::vec ret(2);
  ret(0)=del_m; ret(1)= del_b;
  return ret;
}


void leapfrog_lognorm_hyper(const unsigned nstep,const double delta, arma::vec &v_old,  arma::vec &p_lam, arma::vec &params, 
                            const arma::vec &nj_val1, const arma::vec &nj_val2, 
                            const arma::vec &survtime1, const arma::vec &ss_survtime1,
                            const arma::vec &survtime2, const arma::vec &ss_survtime2,
                            const double a0,const  double df0, 
                            const double mu_m, const  double mu_v,
                            const double b_m, const double b_v){
  for(unsigned i=0;i<nstep;++i){
    // Rcpp::Rcout<<"flag -1 leap i="<<i<<endl;
    params+=(delta)*(p_lam -(delta/2)*v_old);
    // Rcpp::Rcout<<"flag 0 leap i="<<i<<endl;
    params(0)=std::clamp(params(0), -1e3, 1e5);
    params(1)=std::clamp(params(1), -1e1, 1e1);
    
    arma::vec v_new=-delpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                                   params, a0, df0, mu_m, mu_v, b_m, b_v);
    // Rcpp::Rcout<<"flag 1.5 leap i="<<i<<endl;
    p_lam-=(delta/2)* ( v_old+v_new);
    v_old=v_new;
  }
}

void update_lognorm_hyper(const arma::vec &nj_val1, const arma::vec &nj_val2, 
                          const arma::vec &survtime1, const arma::vec &ss_survtime1,
                          const arma::vec &survtime2, const arma::vec &ss_survtime2,
                          const double a0,const  double df0, 
                          const double mu_m, const  double mu_v,
                          const double b_m, const double b_v,
                          const arma::vec &del_range, const unsigned nleapfrog, arma::vec &params, unsigned &acceptance){
  double ll_old=-logpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                               params, a0, df0, mu_m, mu_v, b_m, b_v);
  arma::vec v_old=-delpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                                 params, a0, df0, mu_m, mu_v, b_m, b_v);
  
  arma::vec p_lam(2, fill::randn);
  double kin_energy= dot(p_lam, p_lam);
  
  arma::vec params_new=params, v_new=v_old;
  
  unsigned pois_draw=(unsigned) R::rpois(nleapfrog);
  unsigned nstep=GSL_MAX_INT(1,pois_draw);
  double delta= R::runif(del_range(0),del_range(1) );
  // Rcpp::Rcout<<"Flag X params -.5!!"<<endl;
  leapfrog_lognorm_hyper(nstep, delta, v_new, p_lam, params_new, 
                         nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                         a0, df0, mu_m, mu_v, b_m, b_v);
  // params_new.print("params_new");
  // Rcpp::Rcout<<"Flag X params 0!!"<<endl;
  double ll_new=- logpi_lognorm(nj_val1, nj_val2, survtime1, ss_survtime1, survtime2, ss_survtime2,
                                params_new, a0, df0, mu_m, mu_v, b_m, b_v); //log-likelihood
  
  //get H_ll_new
  double H_new= ll_new+  ssq(p_lam) /2, H_old=ll_old+kin_energy/2;
  
  if(log(randu())< -(H_new-H_old) ){
    params=params_new;
    /*ll_old=ll_new;
     v_old=v_new;*/
    ++acceptance;
  }
  /*Rcpp::Rcout<<"i= "<<i<< "Acceptance rate="<<(((double)acceptance)/ ((double) i))<<" alpha= "<<exp(l_alpha)<<endl;
   Rcpp::Rcout<<"Acceptance rate="<<(acceptance/n_mc)<<endl;
   return exp(alpha_vec);*/
}

/*****************************************************************************/



// [[Rcpp::export]]
Rcpp::List background_MCMC_storage(const arma::uvec &dat_index,
                                   const arma::uvec &trt_index,
                                   arma::vec &st, arma::uvec &nu,
                                   arma::uvec &del,
                                   arma::imat &eta, arma::mat &eta_cont,
                                   Rcpp::List non_na_obs1,
                                   Rcpp::List non_na_obs1_cont,
                                   const unsigned nmix, arma::uvec ncat,
                                   const double a0, const double df0,
                                   const arma::vec &mu_m_t,const arma::vec &b_v_t,
                                   const arma::vec &b_m_t, const arma::vec &mu_v_t,
                                   const arma::vec &del_range_lognorm_ref,
                                   const unsigned nleapfrog_lognorm_ref,
                                   const arma::vec &del_range_lognorm_oth,
                                   const unsigned nleapfrog_lognorm_oth,
                                   const arma::vec &alpha_hyper,
                                   const arma::vec &del_range_alp1,
                                   const unsigned nleapfrog_alp1,
                                   const arma::vec &del_range_alp2,
                                   const unsigned nleapfrog_alp2,
                                   const int nrun,const int burn,
                                   const int thin)
{
  const unsigned n = trt_index.n_elem;
  const unsigned n_trt = max(trt_index) + 1;
  const unsigned num_cohort = max(dat_index) + 1;
  
  // acceptance per dataset (for the covariates alpha_s HMC)
  arma::uvec acceptance_alph(num_cohort, arma::fill::zeros);
  double sig_alp = std::sqrt(std::log1p(alpha_hyper(1) / gsl_pow_2(alpha_hyper(0))));
  double mu_alp  = std::log(alpha_hyper(0)) - 0.5*sig_alp*sig_alp;
  
  // acceptance per treatment (for the lognormal hyper HMC)
  arma::uvec acceptance_y(n_trt, arma::fill::zeros);
  // double mu0 = mu_m;
  // double beta0 = std::exp(b_v/2.0 + b_m);
  // arma::vec current_params(2); current_params(0)=mu0; current_params(1)=std::log(beta0);
  arma::vec mu0_t(n_trt);
  arma::vec logb0_t(n_trt);
  for (unsigned t = 0; t < n_trt; ++t) {
    mu0_t(t)   = mu_m_t(t);      // prior mean for μ0,t
    logb0_t(t) = b_m_t(t) + 0.5 * b_v_t(t);       // prior mean for log β0,t  (not log E[β])
  }
  
  arma::mat params_t(n_trt, 2);
  for (unsigned t = 0; t < n_trt; ++t) {
    params_t(t, 0) = mu0_t(t);
    params_t(t, 1) = logb0_t(t);
  }
  
  const unsigned k      = eta.n_cols;
  const unsigned k_cont = eta_cont.n_cols;
  
  arma::mat eta_cont_sq = arma::square(eta_cont);
  const double max_st = 5*st.max();
  
  // unpack non-missing index lists
  arma::field<arma::uvec> non_na_obs(n), non_na_obs_cont(n);
  for (unsigned i=0; i<n; ++i){
    non_na_obs(i)      = Rcpp::as<arma::uvec>(non_na_obs1[i]);
    non_na_obs_cont(i) = Rcpp::as<arma::uvec>(non_na_obs1_cont[i]);
  }
  
  // storage
  const unsigned n_mc = std::floor((nrun+burn)/thin);
  arma::umat alloc_var_mat(n_mc, n);
  
  arma::field<arma::mat> weightsSyn(num_cohort > 0 ? num_cohort-1 : 0);
  for (unsigned s=0; s+1<num_cohort; ++s){
    arma::uvec idx = arma::find(dat_index == (s+1));
    weightsSyn(s).set_size(n_mc, idx.n_elem);
  }
  
  arma::cube pimat(n_mc, nmix, num_cohort, arma::fill::none);
  arma::cube lognormal_mu(n_mc, nmix, n_trt, arma::fill::none);
  arma::cube lognormal_sig(n_mc, nmix, n_trt, arma::fill::none);
  arma::mat  dir_alpha_mat(n_mc, num_cohort, arma::fill::none);
  arma::cube  hyperparams(n_mc, 2, n_trt, arma::fill::none);
  arma::mat  unifmat(n_mc, n, arma::fill::none);
  
  arma::cube nj_val_mat(nmix, num_cohort, n_mc, arma::fill::none);
  
  Rcpp::List noccu_list(n_mc);
  arma::ucube nobs_mat(nmix, k, n_mc, arma::fill::zeros);
  arma::ucube nj_x_mat(nmix, k_cont, n_mc, arma::fill::zeros);
  arma::cube  sum_j_x_mat(nmix, k_cont, n_mc, arma::fill::zeros);
  arma::cube  ss_j_x_mat(nmix , k_cont, n_mc, arma::fill::zeros);
  arma::mat   nj_val_all_mat(nmix, n_mc, arma::fill::none);
  arma::cube  nj_val_shared_mat(nmix, n_trt, n_mc, arma::fill::none);
  arma::cube  survtime_mat(nmix, n_trt, n_mc, arma::fill::none);
  arma::cube  ss_survtime_mat(nmix, n_trt, n_mc, arma::fill::none);
  
  // working buffers
  arma::uvec d(nmix, arma::fill::zeros);
  arma::vec  probs(nmix, arma::fill::zeros);
  arma::vec  log_probs(nmix, arma::fill::value(arma::datum::log_min));
  arma::vec  alpha_vec(k, arma::fill::ones);
  
  // categorical atoms
  arma::field<arma::uvec> noccu(nmix, k);
  for (unsigned j=0;j<nmix;++j){
    for (unsigned v=0; v<k; ++v){
      noccu(j,v).set_size(ncat(v));
      noccu(j,v).zeros();
    }
  }
  
  // continuous atoms
  arma::umat nj_x(nmix, k_cont, arma::fill::zeros);
  arma::mat  sum_j_x(nmix, k_cont, arma::fill::zeros);
  arma::mat  ss_j_x(nmix ,k_cont, arma::fill::zeros);
  
  arma::field<arma::umat> inds_eq_j_shared(nmix, n_trt);
  arma::field<arma::umat> inds_eq_j(nmix, num_cohort);
  arma::field<arma::uvec> inds_eq_j_all(nmix);
  arma::vec  nj_val_all(nmix, arma::fill::zeros);
  
  // initialize cluster stats from current labels
  for (unsigned j=0;j<nmix;++j){
    inds_eq_j_all(j) = arma::find(del == j);
    nj_val_all(j)    = inds_eq_j_all(j).n_elem;
    
    for (auto i : inds_eq_j_all(j)){
      for (auto v : non_na_obs(i)) ++noccu(j,v)(eta(i,v));
      for (auto c : non_na_obs_cont(i)){
        ++nj_x(j,c);
        sum_j_x(j,c) += eta_cont(i,c);
        ss_j_x(j,c)  += eta_cont_sq(i,c);
      }
    }
  }
  
  // cohort-wise occupancies
  arma::mat nj_val(nmix, num_cohort, arma::fill::zeros);
  for (unsigned j=0;j<nmix;++j){
    for (unsigned s=0; s<num_cohort; ++s){
      arma::uvec idx_s = arma::find(dat_index == s);
      arma::uvec in_j  = arma::find(del == j);
      inds_eq_j(j,s)   = arma::intersect(idx_s, in_j);
      nj_val(j,s)      = inds_eq_j(j,s).n_elem;
    }
  }
  
  // continuous-covariate prior (simulation-friendly default)
  const double df_x=1.0, alpha_x = k_cont + 30.0;
  const double beta_x = 1;
  const double mu_x   = 0;
  
  // response stats per (cluster, treatment)
  arma::mat survtime(nmix, n_trt, arma::fill::zeros);
  arma::mat ss_survtime(nmix, n_trt, arma::fill::zeros);
  arma::mat nj_val_shared(nmix, n_trt, arma::fill::zeros);
  for (unsigned j=0;j<nmix;++j){
    for (unsigned t=0;t<n_trt;++t){
      arma::uvec idx_t = arma::find(trt_index == t);
      arma::uvec in_j  = arma::find(del == j);
      inds_eq_j_shared(j,t) = arma::intersect(idx_t, in_j);
      nj_val_shared(j,t)    = inds_eq_j_shared(j,t).n_elem;
      if (nj_val_shared(j,t)){
        survtime(j,t)   = arma::sum(st(inds_eq_j_shared(j,t)));
        ss_survtime(j,t)= ssq(st(inds_eq_j_shared(j,t)));
      }
    }
  }
  
  // GSL rng (for Dirichlet/multinomial)
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, 500);
  
  const arma::vec  st_original = st;
  const arma::uvec censored_indices = arma::find(nu == 0);
  
  // nobs per (cluster, categorical var)
  arma::umat nobs(nmix, k, arma::fill::zeros);
  for (unsigned j=0;j<nmix;++j)
    for (unsigned v=0; v<k; ++v)
      nobs(j,v) = arma::sum(noccu(j,v));
  
  // per-cohort non-empty
  arma::field<arma::uvec> non_empty_clusters(num_cohort);
  for (unsigned s=0;s<num_cohort;++s)
    non_empty_clusters(s) = arma::find(nj_val.col(s));
  
  arma::vec tmpp = nj_val_all;
  arma::uvec occu_hist = arma::find(tmpp);
  unsigned nmix1 = occu_hist.n_elem;
  
  arma::vec alpha(num_cohort, arma::fill::ones);
  arma::vec dir_prec(num_cohort, arma::fill::ones);
  dir_prec = alpha / nmix;                 // precision per cohort
  arma::vec prob_empty = arma::log(dir_prec);
  
  arma::vec unif(n, arma::fill::zeros);
  
  // ======================= MCMC ==========================
  for (unsigned i=0; i < (unsigned)(nrun + burn); ++i){
    
    // --- augment censored Y (cluster-wise)
    for (unsigned j=0;j<nmix;++j){
      for (unsigned t=0;t<n_trt;++t){
        const unsigned njt = nj_val_shared(j,t);
        const double mu0_t_cur  = params_t(t,0);
        const double beta0_t_cur = std::exp(params_t(t,1));
        if (!njt) continue;
        
        const double df_post    = df0 + njt - 1.0;
        const double alpha_post = a0  + (njt - 1.0)/2.0;
        const double tmp_shrink = ( (njt - 1.0) * df0 ) / df_post;
        
        for (auto ii : inds_eq_j_shared(j,t)){
          if (nu(ii)) continue; // only censored
          double ss_j, sum_j, mean_j;
          if (njt > 1){
            ss_j  = ss_survtime(j,t) - gsl_pow_2(st(ii));
            sum_j = survtime(j,t)    - st(ii);
            mean_j= sum_j / (njt - 1.0);
          } else {
            ss_j = sum_j = mean_j = 0.0;
          }
          double mu_post   = (df0*mu0_t_cur + sum_j) / df_post;
          double beta_post = beta0_t_cur + ( ss_j - (njt - 1.0)*gsl_pow_2(mean_j)
                                               + tmp_shrink * gsl_pow_2(mean_j - mu0_t_cur) )/2.0;
          double sigma_post= std::sqrt( beta_post * (df_post+1.0)/(df_post*alpha_post) );
          st(ii) = r_trunclst(2.0*alpha_post, mu_post, sigma_post, st_original(ii), std::numeric_limits<double>::max());
          st(ii) = GSL_MIN_DBL(st(ii), max_st );
          // update stats
          ss_survtime(j,t) = ss_j + gsl_pow_2(st(ii));
          survtime(j,t)    = sum_j + st(ii);
        }
      }
    }
    
    
    // --- update allocations
    for (unsigned ii=0; ii<n; ++ii){
      const unsigned cur_j = del(ii);
      const unsigned s     = dat_index(ii);
      const unsigned t     = trt_index(ii);
      const double st_sq   = gsl_pow_2(st(ii));
      const double mu0_t_cur  = params_t(t,0);
      const double beta0_t_cur = std::exp(params_t(t,1));
      
      
      // empty-component predictive (common to all)
      double dens_empty = surv_fn_lognorm(st_original(ii), nu(ii),
                                          0.0, 0.0, df0, a0, mu0_t_cur, beta0_t_cur, 0);
      for (auto c : non_na_obs_cont(ii))
        dens_empty += post_t_dens(eta_cont(ii,c), 0.0, 0.0, df_x, alpha_x, mu_x, beta_x, 0);
      
      log_probs.fill(arma::datum::log_min);
      
      // score all j != cur_j
      for (unsigned j=0;j<nmix;++j){
        if (j == cur_j) continue;
        
        // Rcpp::Rcout << "j = " << j << std::endl;
        double dens, lp_clust;
        if (nj_val_all(j) == 0){                       // empty atom
          dens     = dens_empty;
          lp_clust = std::log(alpha(s)/nmix);
        } else {
          dens = surv_fn_lognorm(st_original(ii), nu(ii),
                                 ss_survtime(j,t), survtime(j,t),
                                 df0, a0, mu0_t_cur, beta0_t_cur, nj_val_shared(j,t));
          // Rcpp::Rcout << "st_original(ii): " << st_original(ii) << std::endl;
          // Rcpp::Rcout << "ss_survtime(j, t): " << ss_survtime(j, t) << std::endl;
          // Rcpp::Rcout << "survtime(j, t): " << survtime(j, t) << std::endl;
          // Rcpp::Rcout << "st(ii): " << st(ii) << std::endl;
          // Rcpp::Rcout << "st_sq: " << st_sq << std::endl;
          // Rcpp::Rcout << "nj_val_all(j): " << nj_val_all(j) << std::endl;
          // Rcpp::Rcout << "nj_val_shared(j, t): " << nj_val_shared(j, t) << std::endl;
          // Rcpp::Rcout << "Y part: " << dens << std::endl;
          // categorical X
          for (auto v : non_na_obs(ii))
          {
            dens += ( std::log( noccu(j,v)(eta(ii,v)) + alpha_vec(v) )
                        -  std::log( nobs(j,v) + ncat(v)*alpha_vec(v) ) );
            // Rcpp::Rcout << "X part categorical: " << dens << std::endl;
          }
          
          // continuous X
          for (auto c : non_na_obs_cont(ii))
          {
            dens += post_t_dens(eta_cont(ii,c),
                                ss_j_x(j,c), sum_j_x(j,c),
                                df_x, alpha_x, mu_x, beta_x, nj_x(j,c));
            // Rcpp::Rcout << "X part continuous: " << dens << std::endl;
          }
          
          
          lp_clust = std::log(nj_val(j,s) + alpha(s)/nmix);
        }
        log_probs(j) = dens + lp_clust;
      }
      
      
      // score cur_j (leave-one-out)
      double dens_cur, lp_cur;
      // Rcpp::Rcout << "j = " << cur_j << std::endl;
      if (nj_val_all(cur_j) == 1){
        dens_cur = dens_empty;
        lp_cur   = std::log(alpha(s)/nmix);
      } else {
        dens_cur = surv_fn_lognorm(st_original(ii), nu(ii),
                                   ss_survtime(cur_j,t) - st_sq,
                                   survtime(cur_j,t)    - st(ii),
                                   df0, a0, mu0_t_cur, beta0_t_cur, nj_val_shared(cur_j,t) - 1);
        // Rcpp::Rcout << "st_original(ii): " << st_original(ii) << std::endl;
        // Rcpp::Rcout << "ss_survtime(cur_j, t): " << ss_survtime(cur_j, t) << std::endl;
        // Rcpp::Rcout << "survtime(cur_j, t): " << survtime(cur_j, t) << std::endl;
        // Rcpp::Rcout << "st(ii): " << st(ii) << std::endl;
        // Rcpp::Rcout << "st_sq: " << st_sq << std::endl;
        // Rcpp::Rcout << "nj_val_all(cur_j): " << nj_val_all(cur_j) << std::endl;
        // Rcpp::Rcout << "nj_val_shared(cur_j, t): " << nj_val_shared(cur_j, t) << std::endl;
        // Rcpp::Rcout << "Y part: " << dens_cur << std::endl;
        for (auto v : non_na_obs(ii))
        {
          dens_cur += ( std::log( noccu(cur_j,v)(eta(ii,v)) - 1 + alpha_vec(v) )
                          -  std::log( nobs(cur_j,v) - 1 + ncat(v)*alpha_vec(v) ) );
          // Rcpp::Rcout << "X part categorical: " << dens_cur << std::endl;
        }
        
        for (auto c : non_na_obs_cont(ii))
        {
          dens_cur += post_t_dens(eta_cont(ii,c),
                                  ss_j_x(cur_j,c) - eta_cont_sq(ii,c),
                                  sum_j_x(cur_j,c)- eta_cont(ii,c),
                                  df_x, alpha_x, mu_x, beta_x, nj_x(cur_j,c) - 1);
          // Rcpp::Rcout << "X part continuous: " << dens_cur << std::endl;
        }
        
        lp_cur = std::log(nj_val(cur_j,s) - 1 + alpha(s)/nmix);
      }
      log_probs(cur_j) = dens_cur + lp_cur;
      
      // log_probs.print("take a closer look log: ");
      
      // normalize and sample
      const double log_DEN = log_sum_exp(log_probs);
      probs = arma::exp(log_probs - log_DEN);
      
      // probs.print("take a closer look: ");
      
      gsl_ran_multinomial(r, nmix, 1, probs.memptr(), d.memptr());
      arma::uvec dd = arma::find(d==1, 1, "first");
      const unsigned new_j = dd(0);
      
      if (new_j != cur_j){
        // decrement
        --nj_val(cur_j, s);
        --nj_val_all(cur_j);
        --nj_val_shared(cur_j, t);
        survtime(cur_j, t)    -= st(ii);
        ss_survtime(cur_j, t) -= st_sq;
        arma::uvec tmp_current_ind={cur_j}, tmp_ii={ii};
        nj_x(tmp_current_ind, non_na_obs_cont(ii))-=1;
        
        ss_j_x(tmp_current_ind,non_na_obs_cont(ii))-= eta_cont_sq(tmp_ii,non_na_obs_cont(ii));
        
        sum_j_x(tmp_current_ind,non_na_obs_cont(ii))-= eta_cont(tmp_ii,non_na_obs_cont(ii));
        
        // increment
        ++nj_val(new_j, s);
        ++nj_val_all(new_j);
        ++nj_val_shared(new_j, t);
        survtime(new_j, t)    += st(ii);
        ss_survtime(new_j, t) += st_sq;
        arma::uvec tmp_del_ii={new_j};
        
        nj_x(tmp_del_ii,non_na_obs_cont(ii))+=1;
        
        ss_j_x(tmp_del_ii,non_na_obs_cont(ii)) += eta_cont_sq(tmp_ii,non_na_obs_cont(ii));
        
        sum_j_x(tmp_del_ii,non_na_obs_cont(ii))+= eta_cont(tmp_ii,non_na_obs_cont(ii));
        
        del(ii) = new_j;
      }
      
      tmpp = nj_val_all; // refresh
      
      // Rcpp::Rcout << "Subject ii: " << ii << std::endl;
      
    } // end i over subjects
    
    // refresh supports
    for (unsigned s=0;s<num_cohort;++s)
      non_empty_clusters(s) = arma::find(nj_val.col(s));
    occu_hist = arma::find(tmpp);
    nmix1     = occu_hist.n_elem;
    
    // --- HEAVY CONSISTENCY CHECK (do once per sweep) ---
    for (unsigned j = 0; j < nmix; ++j) {
      // (A) all-subjects view
      if (nj_val_all(j)) {
        inds_eq_j_all(j) = arma::find(del == j);
      } else {
        inds_eq_j_all(j).reset();
      }
      if (inds_eq_j_all(j).n_elem != nj_val_all(j)) {
        Rcpp::Rcout << "Mismatch: inds_eq_j_all size=" << inds_eq_j_all(j).n_elem
                    << " vs nj_val_all(" << j << ")=" << nj_val_all(j) << std::endl;
        Rcpp::stop("inds_eq_j_all vs nj_val_all mismatch");
      }
      
      // (B) per-cohort view
      for (unsigned s = 0; s < num_cohort; ++s) {
        if (nj_val(j, s)) {
          arma::uvec sidx = arma::find(dat_index == s);
          arma::uvec didx = arma::find(del == j);
          inds_eq_j(j, s) = arma::intersect(sidx, didx);
        } else {
          inds_eq_j(j, s).reset();
        }
        if (inds_eq_j(j, s).n_elem != nj_val(j, s)) {
          Rcpp::Rcout << "Mismatch: inds_eq_j(" << j << "," << s << ") size="
                      << inds_eq_j(j, s).n_elem
                      << " vs nj_val(" << j << "," << s << ")="
                      << nj_val(j, s) << std::endl;
          Rcpp::stop("inds_eq_j(.,.) vs nj_val(.,.) mismatch");
        }
      }
      
      // (C) per-treatment view
      for (unsigned t = 0; t < n_trt; ++t) {
        if (nj_val_shared(j, t)) {
          arma::uvec tidx = arma::find(trt_index == t);
          arma::uvec didx = arma::find(del == j);
          inds_eq_j_shared(j, t) = arma::intersect(tidx, didx);
        } else {
          inds_eq_j_shared(j, t).reset();
        }
        if (inds_eq_j_shared(j, t).n_elem != nj_val_shared(j, t)) {
          Rcpp::Rcout << "Mismatch: inds_eq_j_shared(" << j << "," << t << ") size="
                      << inds_eq_j_shared(j, t).n_elem
                      << " vs nj_val_shared(" << j << "," << t << ")="
                      << nj_val_shared(j, t) << std::endl;
          Rcpp::stop("inds_eq_j_shared vs nj_val_shared mismatch");
        }
      }
    }
    // --- END HEAVY CONSISTENCY CHECK ---
    
    
    
    // mixture weights
    arma::mat pi(nmix, num_cohort, arma::fill::zeros);
    
    // sample outcome params per (j,t)
    arma::mat mu(nmix, n_trt, arma::fill::zeros), sig(nmix, n_trt, arma::fill::zeros);
    for (unsigned j=0;j<nmix;++j){
      for (unsigned t=0;t<n_trt;++t){
        const double mu0_t_cur  = params_t(t,0);
        const double beta0_t_cur = std::exp(params_t(t,1));
        const unsigned nC = static_cast<unsigned>(nj_val_shared(j,t));
        const double sC   = survtime(j,t);
        const double ssC  = ss_survtime(j,t);
        arma::vec draw = sim_lognorm_params(ssC, sC, df0, a0, mu0_t_cur, beta0_t_cur, nC);
        mu(j,t)  = draw(0);
        sig(j,t) = draw(1);
      }
    }
    
    // RWD cohorts: Dirichlet over nmix with nj_val + prior
    for (unsigned s=1;s<num_cohort;++s){
      arma::vec dir_alpha_s = nj_val.col(s) + (alpha(s)/nmix);
      arma::vec out(nmix);
      gsl_ran_dirichlet(r, nmix, dir_alpha_s.memptr(), out.memptr());
      pi.col(s) = out;
    }
    
    // current cohort: Dirichlet over ALL nmix
    {
      arma::vec dir_alpha0 = nj_val.col(0) + (alpha(0)/nmix);
      arma::vec out(nmix);
      gsl_ran_dirichlet(r, nmix, dir_alpha0.memptr(), out.memptr());
      pi.col(0) = out;
    }
    
    // synthetic weights (unchanged structure, using current dir_alpha0)
    arma::field<arma::vec> wght_xSyn_new(num_cohort > 0 ? num_cohort-1 : 0);
    for (unsigned s=1; s<num_cohort; ++s){
      arma::uvec one = {s};
      arma::vec mean_pi1(nmix, arma::fill::zeros);
      arma::vec wghtSyn(nmix,   arma::fill::zeros);
      arma::uvec ne = non_empty_clusters(s);
      if (ne.n_elem){
        arma::vec base = nj_val(ne, one) + (alpha(0)/nmix);
        mean_pi1(ne) = arma::normalise(base - (alpha(0)/nmix), 1.0);
        wghtSyn(occu_hist) = mean_pi1(occu_hist) / nj_val(occu_hist, one);
      }
      wght_xSyn_new(s-1) = wghtSyn( del( arma::find(dat_index == s) ) );
    }
    
    // PIT uniforms
    for (unsigned t=0; t<n_trt; ++t){
      arma::uvec tmp  = arma::find(trt_index == t);   // rows in arm t
      if (tmp.is_empty()) continue;
      
      arma::uvec tmp1 = del.elem(tmp);                // clusters for those rows
      
      arma::vec mu_t  = mu.col(t);                    // nmix x 1
      arma::vec sig_t = sig.col(t);                   // nmix x 1 (variance)
      
      arma::vec mu_i  = mu_t.elem(tmp1);              // per-row means
      arma::vec sd_i  = arma::sqrt(sig_t.elem(tmp1)); // per-row SDs
      
      unif.elem(tmp) = arma::normcdf(st_original.elem(tmp), mu_i, sd_i);
    }
    
    unif(censored_indices) += arma::randu(censored_indices.n_elem) % (1.0 - unif(censored_indices));
    
    // HMC updates
    for (unsigned t = 0; t < n_trt; ++t) {
      
      // choose tuning based on whether t is the reference arm
      const arma::vec &del_rng_t = (t == 0) ? del_range_lognorm_ref
      : del_range_lognorm_oth;
      const unsigned nlf_t       = (t == 0) ? nleapfrog_lognorm_ref
      : nleapfrog_lognorm_oth;
      
      // Build per-t tallies: each is nmix × 1 (columns = that treatment only)
      arma::mat NJ_t  = nj_val_shared.col(t);
      arma::mat S_t   = survtime.col(t);
      arma::mat SS_t  = ss_survtime.col(t);
      
      arma::vec cp = params_t.row(t).t();      // [μ0,t, log β0,t]
      
      update_lognorm_hyper_extend(NJ_t, S_t, SS_t,
                                  a0, df0, mu_m_t(t), mu_v_t(t), b_m_t(t), b_v_t(t),
                                  del_rng_t, nlf_t,
                                  cp, acceptance_y(t), t);
      params_t(t,0) = cp(0);
      params_t(t,1) = cp(1);
    }
    
    // alpha_0 (current)
    {
      double l_alpha = std::log(alpha(0));
      arma::uvec one = {0};
      arma::uvec ne  = non_empty_clusters(0);
      update_alpha(nmix1, nj_val(ne, one), mu_alp, sig_alp,
                   del_range_alp1, nleapfrog_alp1, l_alpha, acceptance_alph(0));
      alpha(0)   = std::exp(l_alpha);
      dir_prec(0)= alpha(0)/nmix;
      prob_empty(0) = std::log(dir_prec(0));
    }
    // alpha_s, s>=1
    for (unsigned s=1;s<num_cohort;++s){
      double l_alpha = std::log(alpha(s));
      arma::uvec one = {s};
      arma::uvec ne  = non_empty_clusters(s);
      update_alpha(nmix, nj_val(ne, one), mu_alp, sig_alp,
                   del_range_alp2, nleapfrog_alp2, l_alpha, acceptance_alph(s));
      alpha(s)   = std::exp(l_alpha);
      dir_prec(s)= alpha(s)/nmix;
      prob_empty(s) = std::log(dir_prec(s));
    }
    
    
    if ( ((i + 1) % thin) == 0 ) {
      double denom = static_cast<double>(i + 1);  // avoids /0 and forces floating division
      Rcpp::Rcout << "MCMC iteration: " << i + 1 << std::endl;
      arma::vec acc_rates_y = arma::conv_to<arma::vec>::from(acceptance_y) / denom;
      acc_rates_y.t().print("Treatment Acceptance (lognormal hyper):");
      
      arma::vec acc_rates = arma::conv_to<arma::vec>::from(acceptance_alph) / denom;
      acc_rates.t().print("Dataset Acceptance (alpha):");
    }
    
    
    // thinning
    const unsigned rem = (i+1) % thin;
    if (rem == 0){
      const unsigned q = (i+1)/thin;
      alloc_var_mat.row(q-1) = del.t();
      dir_alpha_mat.row(q-1) = alpha.t();
      for (unsigned s=0;s<num_cohort;++s) pimat.slice(s).row(q-1) = pi.col(s).t();
      for (unsigned s=0; s+1<num_cohort; ++s) weightsSyn(s).row(q-1) = wght_xSyn_new(s).t();
      for (unsigned t=0;t<n_trt;++t){
        lognormal_mu.slice(t).row(q-1)  = mu.col(t).t();
        lognormal_sig.slice(t).row(q-1) = sig.col(t).t();
        hyperparams(q-1, 0, t) = params_t(t, 0);
        hyperparams(q-1, 1, t) = std::exp(params_t(t, 1));
      }
      unifmat.row(q-1)     = unif.t();
      
      noccu_list[q-1]          = Rcpp::wrap(noccu);
      nj_x_mat.slice(q-1)      = nj_x;
      sum_j_x_mat.slice(q-1)   = sum_j_x;
      ss_j_x_mat.slice(q-1)    = ss_j_x;
      nj_val_all_mat.col(q-1)  = nj_val_all;
      nj_val_mat.slice(q-1)    = nj_val;
      nj_val_shared_mat.slice(q-1)= nj_val_shared;
      survtime_mat.slice(q-1)  = survtime;
      ss_survtime_mat.slice(q-1)= ss_survtime;
      nobs_mat.slice(q-1)      = nobs;
    }
    // Rcpp::Rcout << "MCMC iteration: " << i+1 << std::endl;
  } // end MCMC
  
  arma::vec acceptance_all(num_cohort + n_trt);
  acceptance_all.subvec(0, num_cohort-1) =
    arma::conv_to<arma::vec>::from(acceptance_alph)/((double)(nrun+burn));
  acceptance_all.subvec(num_cohort, (num_cohort + n_trt - 1)) = 
    arma::conv_to<arma::vec>::from(acceptance_y)/((double)(nrun+burn));
  
  arma::ivec MCMC_iter = {nrun, burn, thin};
  
  gsl_rng_free(r);
  
  return Rcpp::List::create(
    Rcpp::Named("picube")                = pimat,
    Rcpp::Named("Weights_Syn")           = weightsSyn,
    Rcpp::Named("Lognormal_Mu_Cube")     = lognormal_mu,
    Rcpp::Named("Lognormal_Sig_Cube")    = lognormal_sig,
    Rcpp::Named("Unifs")                 = unifmat,
    Rcpp::Named("Lognormal_hyperparams") = hyperparams,
    Rcpp::Named("Dirichlet_params")      = dir_alpha_mat,
    Rcpp::Named("Acceptance_rates")      = acceptance_all,
    Rcpp::Named("Allocation_variables")  = alloc_var_mat,
    Rcpp::Named("noccu_list")            = noccu_list,
    Rcpp::Named("nj_x_cube")             = nj_x_mat,
    Rcpp::Named("sum_j_x_cube")          = sum_j_x_mat,
    Rcpp::Named("ss_j_x_cube")           = ss_j_x_mat,
    Rcpp::Named("nj_val_all_mat")        = nj_val_all_mat,
    Rcpp::Named("nj_val_cube")           = nj_val_mat,
    Rcpp::Named("nj_val_shared_cube")    = nj_val_shared_mat,
    Rcpp::Named("survtime_cube")         = survtime_mat,
    Rcpp::Named("ss_survtime_cube")      = ss_survtime_mat,
    Rcpp::Named("nobs_cube")             = nobs_mat,
    Rcpp::Named("MCMC_iter")             = MCMC_iter
  );
}

// [[Rcpp::export]]
Rcpp::List common_atoms_cat_lognormal_shared_approx(const arma::uvec &dat_index,
                                                    const arma::uvec &trt_index,
                                                    arma::vec &st, arma::uvec &nu,
                                                    arma::uvec &del,
                                                    arma::imat &eta, arma::mat &eta_cont,
                                                    Rcpp::List non_na_obs1,
                                                    Rcpp::List non_na_obs1_cont,
                                                    arma::imat &eta_pred, arma::mat &eta_cont_pred,
                                                    Rcpp::List non_na_obs1_pred,
                                                    Rcpp::List non_na_obs1_cont_pred,
                                                    const unsigned nmix, arma::uvec ncat,
                                                    const double a0, const double df0,
                                                    const arma::vec &mu_m_t,const arma::vec &b_v_t,
                                                    const arma::vec &b_m_t, const arma::vec &mu_v_t,
                                                    const arma::vec &del_range_lognorm_ref,
                                                    const unsigned nleapfrog_lognorm_ref,
                                                    const arma::vec &del_range_lognorm_oth,
                                                    const unsigned nleapfrog_lognorm_oth,
                                                    const arma::vec &alpha_hyper,
                                                    const arma::vec &del_range_alp1,
                                                    const unsigned nleapfrog_alp1,
                                                    const arma::vec &del_range_alp2,
                                                    const unsigned nleapfrog_alp2,
                                                    const int nrun,const int burn,
                                                    const int thin,
                                                    const arma::uvec &dat_index_old,
                                                    const arma::umat &trt_convert,
                                                    const arma::mat nj_val_all_old,   // [nmix x M]
                                                    const arma::cube  nj_val_old,
                                                    const arma::cube nj_val_shared_old,
                                                    Rcpp::List noccu_old,
                                                    const arma::cube nobs_old,
                                                    const arma::cube ss_j_x_old,
                                                    const arma::cube sum_j_x_old,
                                                    const arma::cube nj_x_old,
                                                    const arma::cube ss_survtime_old,
                                                    const arma::cube survtime_old,
                                                    const arma::cube lognormal_mu_old,
                                                    const arma::cube lognormal_sig_old,
                                                    const arma::mat dirichlet_alpha_mat_old,
                                                    const arma::cube pimat_old,
                                                    const bool freeze_control,
                                                    const arma::vec &mu0_control_draws,
                                                    const arma::vec &beta0_control_draws)
{
  
  // -------- basic sizes --------
  const unsigned n_curr = dat_index.n_elem - dat_index_old.n_elem;   // current sample size
  const unsigned n_trt  = max(trt_index) + 1;                        // # shared trts (new index space)
  const unsigned S      = max(dat_index) + 1;                        // # cohorts (0 = current)
  const unsigned n_pred = eta_pred.n_rows;
  
  if(n_pred == 0){
    eta_pred.reset();
    eta_cont_pred.reset();
  }
  
  // ====== (NEW) Precompute historical-occupied cluster sets per Stage-1 draw ======
  // nj_val_all_old is assumed nmix x M, columns correspond to Stage-1 draw index i
  const unsigned M_hist = nj_val_all_old.n_cols;
  if (M_hist == 0u) Rcpp::stop("nj_val_all_old has zero columns.");
  arma::field<arma::uvec> occu_hist(M_hist);
  for (unsigned i = 0; i < M_hist; ++i) {
    occu_hist(i) = arma::find(nj_val_all_old.col(i) > 0.0); // clusters with historical support at draw i
    if (occu_hist(i).is_empty()) {
      Rcpp::stop("Historical occupancy set is empty for draw i=%u; unexpected under your assumptions.", i);
    }
  }
  // ================================================================================
  
  // -------- hyper for alpha from Stage 1 (empirical) --------
  arma::vec log_alpha_hist_all = arma::log(arma::vectorise(dirichlet_alpha_mat_old));
  double mu_alp  = arma::median(log_alpha_hist_all);
  double sig_alp = arma::stddev(log_alpha_hist_all);
  
  // -------- log-normal hyper start --------
  arma::vec  mu0_t(n_trt);
  arma::vec  logb0_t(n_trt);
  // Rcpp::Rcout << "n_trt: " << n_trt << std::endl;
  // 
  // mu_m_t.print("mu_m_t: ");
  // b_m_t.print("b_m_t: ");
  // b_v_t.print("b_v_t: ");
  // mu_v_t.print("mu_v_t: ");
  // 
  // 
  
  for (unsigned t = 0; t < n_trt; ++t) {
    mu0_t(t)   = mu_m_t(t);
    // You said you're switching to b_m_t + 0.5*b_v_t for the init on log β0,t
    logb0_t(t) = b_m_t(t) + 0.5 * b_v_t(t);
  }
  
  // live parameters that HMC updates each sweep
  arma::mat params_t(n_trt, 2); // [μ0,t, log β0,t]
  for (unsigned t = 0; t < n_trt; ++t) {
    params_t(t, 0) = mu0_t(t);
    params_t(t, 1) = logb0_t(t);
  }
  
  // -------- covariate dims --------
  const unsigned k      = eta.n_cols;
  const unsigned k_cont = eta_cont.n_cols;
  
  arma::mat eta_cont_sq = arma::square(eta_cont);
  const double max_st   = 5.0 * max(st);
  
  // -------- handy accessors for non-NA indices --------
  field<arma::uvec> non_na_obs(n_curr), non_na_obs_cont(n_curr);
  unsigned cat_na_count=0;
  for(unsigned i=0;i<n_curr;++i){
    non_na_obs(i)      = Rcpp::as<arma::uvec>(non_na_obs1[i]);
    non_na_obs_cont(i) = Rcpp::as<arma::uvec>(non_na_obs1_cont[i]);
    cat_na_count += (k - non_na_obs(i).n_elem);
  }
  bool cat_na = (cat_na_count == eta.n_rows * eta.n_cols);
  
  // -------- storage (unchanged shape) --------
  const unsigned n_mc = std::floor((nrun + burn)/thin);
  arma::umat  alloc_var_mat(n_mc, n_curr);
  arma::cube  pimat(n_mc, nmix, S);
  arma::cube  lognormal_mu(n_mc, nmix, n_trt);
  arma::cube  lognormal_sig(n_mc, nmix, n_trt);
  arma::mat   dir_alpha_mat(n_mc, S);
  arma::cube  hyperparams(n_mc, 2, n_trt);
  
  arma::umat  alloc_var_mat_pred(n_mc, n_pred);
  arma::cube  pimat_pred(n_mc, nmix, n_pred, arma::fill::zeros);
  arma::cube  y_pred_cube(n_mc, n_pred, n_trt);
  
  // -------- running tallies for current data --------
  field<arma::uvec> inds_eq_j_all(nmix);
  field<arma::uvec> inds_eq_j_shared(nmix, n_trt);
  arma::mat   nj_val_curr(nmix, S, fill::zeros);
  arma::mat   nj_val_shared_curr(nmix, n_trt, fill::zeros);
  arma::mat   survtime_curr(nmix, n_trt, fill::zeros);
  arma::mat   ss_survtime_curr(nmix, n_trt, fill::zeros);
  
  // cat atoms
  field<arma::uvec> noccu(nmix, k);
  for(unsigned j=0;j<nmix;++j){
    for(unsigned c=0;c<k;++c){
      noccu(j,c).set_size(ncat(c));
      noccu(j,c).zeros();
    }
  }
  arma::umat nobs_curr(nmix, k, fill::zeros);
  
  // cont atoms
  arma::umat nj_x_curr(nmix, k_cont, fill::zeros);
  arma::mat  sum_j_x_curr(nmix, k_cont, fill::zeros);
  arma::mat  ss_j_x_curr(nmix, k_cont, fill::zeros);
  
  // initialize current tallies from initial del
  for(unsigned j=0;j<nmix;++j){
    inds_eq_j_all(j) = find(del == j);
    if(inds_eq_j_all(j).n_elem){
      // per-treatment
      for(unsigned t=0;t<n_trt;++t){
        arma::uvec id_t = find(trt_index.head(n_curr) == t);
        arma::uvec id_j = inds_eq_j_all(j);
        inds_eq_j_shared(j,t) = intersect(id_t, id_j);
        nj_val_shared_curr(j,t) = inds_eq_j_shared(j,t).n_elem;
        if(nj_val_shared_curr(j,t)){
          survtime_curr(j,t)    = sum(st(inds_eq_j_shared(j,t)));
          ss_survtime_curr(j,t) = dot(st(inds_eq_j_shared(j,t)), st(inds_eq_j_shared(j,t)));
        }
      }
      // per-cohort (only s=0 exists in current block)
      nj_val_curr(j,0) = inds_eq_j_all(j).n_elem;
      
      // covariates
      for(auto ii : inds_eq_j_all(j)){
        for(auto c : non_na_obs(ii)){
          ++noccu(j,c)( eta(ii,c) );
          ++nobs_curr(j,c);
        }
        for(auto cc : non_na_obs_cont(ii)){
          ++nj_x_curr(j,cc);
          sum_j_x_curr(j,cc)  += eta_cont(ii,cc);
          ss_j_x_curr(j,cc)   += eta_cont_sq(ii,cc);
        }
      }
    }
  }
  
  // -------- RNG for Dirichlet / Multinomial --------
  const gsl_rng_type *T; gsl_rng *r;
  gsl_rng_env_setup(); T = gsl_rng_default; r = gsl_rng_alloc(T); gsl_rng_set(r, 500);
  
  // -------- constants for covariate priors --------
  const double df_x    = 1.0;
  const double alpha_x = k_cont + 30.0;
  const double beta_x  = 1;
  const double mu_x    = 0;
  
  // -------- helper lambdas --------
  auto map_t_old = [&](unsigned t)->int {
    arma::uvec new_idx = trt_convert.col(1);
    arma::uvec where   = find(new_idx == t, 1, "first");
    if(where.is_empty()) return -1;
    return (int)trt_convert(where(0), 0); // mapped old t
  };
  
  auto combine_resp = [&](unsigned j, unsigned t, unsigned i_draw,
                          double &ssC, double &sC, unsigned &nC){
    ssC = ss_survtime_curr(j,t);
    sC  = survtime_curr(j,t);
    nC  = (unsigned) nj_val_shared_curr(j,t);
    int t_old = map_t_old(t);
    if(t_old >= 0){
      ssC += ss_survtime_old.slice(i_draw)(j, (unsigned)t_old);
      sC  += survtime_old.slice(i_draw)(j, (unsigned)t_old);
      nC  += (unsigned) nj_val_shared_old.slice(i_draw)(j, (unsigned)t_old);
    }
  };
  
  auto combine_cont = [&](unsigned j, unsigned c, unsigned i_draw,
                          double &ssC, double &sC, unsigned &nC){
    ssC = ss_j_x_curr(j,c);
    sC  = sum_j_x_curr(j,c);
    nC  = (unsigned) nj_x_curr(j,c);
    ssC += ss_j_x_old.slice(i_draw)(j,c);
    sC  += sum_j_x_old.slice(i_draw)(j,c);
    nC  += (unsigned) nj_x_old.slice(i_draw)(j,c);
  };
  
  // convert noccu_old into a field
  arma::field<arma::uvec> noccu_hist = build_noccu_field(noccu_old, nmix, k, ncat);
  
  // helpers now become trivial and fast:
  auto combine_cat_counts = [&](unsigned j, unsigned c, unsigned lev, unsigned i_draw)->unsigned {
    unsigned n_cur = (unsigned) noccu(j,c)(lev);       // current tallies
    unsigned n_old = (lev < noccu_hist(j,c,i_draw).n_elem) ? noccu_hist(j,c,i_draw)(lev) : 0u;
    return n_cur + n_old;
  };
  
  auto combined_nobs_cat = [&](unsigned j, unsigned c, unsigned i_draw)->unsigned {
    return (unsigned)nobs_curr(j,c) + arma::accu(noccu_hist(j,c,i_draw));
  };
  
  // -------- working vectors --------
  arma::vec alpha(S, fill::value(1.0));
  arma::vec dir_prec = alpha / nmix; // CURRENT cohort precision uses nmix
  
  const arma::vec st_original = st;
  const arma::uvec censored_indices = find(nu == 0);
  
  arma::uvec acceptance_y(n_trt, fill::zeros);
  arma::uvec acceptance_alph(S, fill::zeros);
  
  // guard once (outside loop is fine too):
  if (freeze_control) {
    if (mu0_control_draws.n_elem != M_hist || beta0_control_draws.n_elem != M_hist) {
      Rcpp::stop("Freeze-control: draw vector length mismatch.");
    }
  }
  
  // ========================= MCMC =========================
  for(int i=0; i<(nrun+burn); ++i){
    
    if (freeze_control) {
      // Pin control (t=0) to Stage-1 per-draw values for draw i
      const double mu0_ctl   = mu0_control_draws((unsigned)i);
      const double beta0_ctl = beta0_control_draws((unsigned)i);
      if (!std::isfinite(mu0_ctl) || !std::isfinite(beta0_ctl) || beta0_ctl <= 0.0)
        Rcpp::stop("Freeze-control: non-finite/invalid μ0 or β0 at draw %d.", i);
      
      params_t(0, 0) = mu0_ctl;                // μ0, t=0
      params_t(0, 1) = std::log(beta0_ctl);    // store log β0, t=0
    }
    
    
    // ---- 1) augment censored (current data only) ----
    for(unsigned j=0;j<nmix;++j){
      for(unsigned t=0;t<n_trt;++t){
        const double mu0_t_cur = params_t(t, 0);
        const double beta0_t_cur = std::exp(params_t(t, 1));
        unsigned njt = (unsigned)nj_val_shared_curr(j,t);
        if(njt){
          double df_post = df0 + (njt - 1);
          double alpha_post = a0 + ( (double)njt - 1.0 )/2.0;
          double tmp = ((njt-1)*df0)/df_post;
          
          for(auto ii: inds_eq_j_shared(j,t)){
            if(!nu(ii)){
              double ss_j=0.0, s_j=0.0; unsigned n_j=0;
              // leave this censored obs out of the current stats
              double st2 = st(ii)*st(ii);
              ss_j = ss_survtime_curr(j,t) - st2;
              s_j  = survtime_curr(j,t)     - st(ii);
              n_j  = njt - 1;
              
              double mean_j = (n_j ? s_j/n_j : 0.0);
              double mu_post = (df0*mu0_t_cur + s_j)/df_post;
              double beta_post = beta0_t_cur + (ss_j - n_j*mean_j*mean_j + tmp*(mean_j - mu0_t_cur)*(mean_j - mu0_t_cur))/2.0;
              double sigma_post = std::sqrt( beta_post * (df_post+1.0)/(df_post*alpha_post) );
              
              st(ii) = r_trunclst(2*alpha_post, mu_post, sigma_post, st_original(ii), std::numeric_limits<double>::max());
              st(ii) = std::min(st(ii), max_st);
              // restore
              ss_survtime_curr(j,t) = ss_j + st(ii)*st(ii);
              survtime_curr(j,t)    = s_j + st(ii);
            }
          }
        }
      }
    }
    
    // Rcpp::Rcout << "Censored observations augmented!" << std::endl;
    
    // ---- 2) update allocations (RESTRICT to historically-occupied clusters) ----
    // identify allowed set for this draw i
    if ((unsigned)i >= M_hist) Rcpp::stop("Historical arrays shorter than (nrun+burn): i=%d >= M_hist=%u", i, M_hist);
    const arma::uvec &A = occu_hist((unsigned)i);  // allowed cluster indices for Stage-1 draw i
    const unsigned A_size = A.n_elem;
    
    for(unsigned ii=0; ii<n_curr; ++ii){
      // Rcpp::Rcout << "Patient ii: " << ii << std::endl;
      unsigned j_cur = del(ii);
      unsigned s  = dat_index(ii);   // should be 0 for current
      unsigned t  = trt_index(ii);
      const double mu0_t_cur   = params_t(t, 0);
      const double beta0_t_cur = std::exp(params_t(t, 1));
      
      arma::vec log_probs(nmix, fill::value(-std::numeric_limits<double>::infinity()));
      
      // empty-predictive pieces for a brand-new cluster (combined empty)
      double log_pdf_empty_y  = surv_fn_lognorm(st_original(ii), nu(ii), 0.0, 0.0, df0, a0, mu0_t_cur, beta0_t_cur, 0);
      double log_pdf_empty_xc = 0.0;
      for(auto c: non_na_obs(ii)){
        (void)c; // placeholder; kept to mirror old structure
      }
      for(auto c: non_na_obs_cont(ii)){
        log_pdf_empty_xc += post_t_dens(eta_cont(ii,c), 0.0, 0.0, df_x, alpha_x, mu_x, beta_x, 0);
      }
      
      double st_sq = st(ii)*st(ii);
      
      // Rcpp::Rcout << "Empty values assigned!" << std::endl;
      
      
      // ONLY iterate over allowed clusters A
      for(auto j : A){
        // cluster prior mass (current cohort only)
        double cluster_prob = std::log( (double)nj_val_curr(j,0) + (j==j_cur ? -1.0 : 0.0) + alpha(0)/nmix );
        
        // outcome likelihood using COMBINED stats
        double ssC, sC; unsigned nC;
        combine_resp(j, t, (unsigned)i, ssC, sC, nC);
        
        if(j==j_cur){
          if(nj_val_shared_curr(j,t) > 0){
            ssC -= st_sq;
            sC  -= st(ii);
            nC  -= 1;
          }
        }
        
        double log_like_y = (nC==0)
          ? log_pdf_empty_y
        : surv_fn_lognorm(st_original(ii), nu(ii), ssC, sC, df0, a0, mu0_t_cur, beta0_t_cur, nC);
        
        // categorical covariates
        double log_like_xcat = 0.0;
        for(auto c: non_na_obs(ii)){
          unsigned lev = (unsigned)eta(ii,c);
          unsigned nlev = combine_cat_counts(j, c, lev, (unsigned)i);
          unsigned nobsC = combined_nobs_cat(j, c, (unsigned)i);
          if(nobsC==0){
            log_like_xcat += std::log( 1.0 / ( (double) ncat(c) ) );
          }else{
            log_like_xcat += ( std::log( (double)nlev + 1.0 )
                                 - std::log( (double)nobsC + (double)ncat(c) ) );
          }
        }
        
        // continuous covariates
        double log_like_xcont = 0.0;
        for(auto c: non_na_obs_cont(ii)){
          double ssXC, sXC; unsigned nXC;
          combine_cont(j, c, (unsigned)i, ssXC, sXC, nXC);
          if(j==j_cur){
            ssXC -= eta_cont_sq(ii,c);
            sXC  -= eta_cont(ii,c);
            nXC  -= 1;
          }
          if(nXC==0){
            log_like_xcont += post_t_dens(eta_cont(ii,c), 0.0, 0.0, df_x, alpha_x, mu_x, beta_x, 0);
          }else{
            log_like_xcont += post_t_dens(eta_cont(ii,c), ssXC, sXC, df_x, alpha_x, mu_x, beta_x, nXC);
          }
        }
        
        log_probs(j) = cluster_prob + log_like_y + log_like_xcat + log_like_xcont;
      }
      
      // Rcpp::Rcout << "Iterated over allowed clusters!" << std::endl;
      
      // normalize
      double log_den = log_sum_exp(log_probs);
      arma::vec probs = arma::exp(log_probs - log_den);
      
      // log_probs.print("log probability: ");
      // probs.print("probability: ");
      
      // sample new label over full nmix (disallowed remain prob=0)
      arma::uvec d(nmix, fill::zeros);
      gsl_ran_multinomial(r, nmix, 1, probs.begin(), d.begin());
      unsigned j_new = (unsigned) d.index_max();
      
      if(j_new != j_cur){
        // update tallies (current) — remove from old
        --nj_val_curr(j_cur,0);
        --nj_val_shared_curr(j_cur,t);
        survtime_curr(j_cur,t)    -= st(ii);
        ss_survtime_curr(j_cur,t) -= st_sq;
        for(auto c: non_na_obs(ii)){
          --noccu(j_cur,c)( eta(ii,c) );
          --nobs_curr(j_cur,c);
        }
        for(auto c: non_na_obs_cont(ii)){
          --nj_x_curr(j_cur,c);
          sum_j_x_curr(j_cur,c)  -= eta_cont(ii,c);
          ss_j_x_curr(j_cur,c)   -= eta_cont_sq(ii,c);
        }
        
        // add to new
        ++nj_val_curr(j_new,0);
        ++nj_val_shared_curr(j_new,t);
        survtime_curr(j_new,t)    += st(ii);
        ss_survtime_curr(j_new,t) += st_sq;
        for(auto c: non_na_obs(ii)){
          ++noccu(j_new,c)( eta(ii,c) );
          ++nobs_curr(j_new,c);
        }
        for(auto c: non_na_obs_cont(ii)){
          ++nj_x_curr(j_new,c);
          sum_j_x_curr(j_new,c)  += eta_cont(ii,c);
          ss_j_x_curr(j_new,c)   += eta_cont_sq(ii,c);
        }
        del(ii) = j_new;
      }
    }
    
    // Rcpp::Rcout << "Starting heavy consistency check!" << std::endl;
    // ================= HEAVY CONSISTENCY CHECK (Stage 2; current cohort only) =================
    {
      // (A) membership vs nj_val_curr(.,0)
      for (unsigned j = 0; j < nmix; ++j) {
        arma::uvec allj = arma::find( del.head(n_curr) == j );
        if (allj.n_elem != (unsigned)nj_val_curr(j,0)) {
          Rcpp::Rcout << "[A] j=" << j
                      << " | count(del==j)=" << allj.n_elem
                      << " vs nj_val_curr(j,0)=" << nj_val_curr(j,0) << std::endl;
          Rcpp::stop("Stage2 check: del vs nj_val_curr(.,0) mismatch");
        }
        inds_eq_j_all(j) = allj;
      }
      
      // (B) per-treatment membership vs nj_val_shared_curr
      for (unsigned j = 0; j < nmix; ++j) {
        arma::uvec allj = inds_eq_j_all(j);
        for (unsigned t = 0; t < n_trt; ++t) {
          arma::uvec tidx = arma::find( trt_index.head(n_curr) == t );
          arma::uvec didx = arma::intersect(tidx, allj);
          
          if (didx.n_elem != (unsigned)nj_val_shared_curr(j,t)) {
            Rcpp::Rcout << "[B] j=" << j << " t=" << t
                        << " | size(intersect(trt=t,del==j))=" << didx.n_elem
                        << " vs nj_val_shared_curr=" << nj_val_shared_curr(j,t) << std::endl;
            Rcpp::stop("Stage2 check: inds_eq_j_shared vs nj_val_shared_curr mismatch");
          }
          inds_eq_j_shared(j,t) = didx;
          
          // (C) response tallies (current only)
          double s_sum = 0.0, ss_sum = 0.0;
          for (auto ii : didx) { s_sum += st(ii); ss_sum += st(ii)*st(ii); }
          if (std::fabs(s_sum  - survtime_curr(j,t))    > 1e-8 ||
              std::fabs(ss_sum - ss_survtime_curr(j,t)) > 1e-8) {
            Rcpp::Rcout << "[C] resp mismatch j=" << j << " t=" << t
                        << " | survtime: have=" << survtime_curr(j,t) << " recalc=" << s_sum
                        << " | ss: have=" << ss_survtime_curr(j,t) << " recalc=" << ss_sum << std::endl;
            Rcpp::stop("Stage2 check: response tallies (current) mismatch");
          }
        }
      }
      // (D) categorical covariates tallies
      for (unsigned j = 0; j < nmix; ++j) {
        const arma::uvec idx = inds_eq_j_all(j);
        for (unsigned v = 0; v < k; ++v) {
          std::vector<unsigned> cnt( ncat(v), 0U );
          unsigned nobs_re = 0U;
          
          for (auto ii : idx) {
            double val = eta(ii,v);
            if (std::isfinite(val)) {
              unsigned lev = (unsigned)val;
              if (lev < cnt.size()) { ++cnt[lev]; ++nobs_re; }
            }
          }
          
          if ((unsigned)nobs_curr(j,v) != nobs_re) {
            Rcpp::Rcout << "[D] cat nobs mismatch j=" << j << " v=" << v
                        << " | have=" << (unsigned)nobs_curr(j,v)
                        << " recalc=" << nobs_re << std::endl;
            Rcpp::stop("Stage2 check: nobs_curr mismatch (categorical)");
          }
          for (unsigned lev = 0; lev < cnt.size(); ++lev) {
            if ((unsigned)noccu(j,v)(lev) != cnt[lev]) {
              Rcpp::Rcout << "[D] cat level mismatch j=" << j << " v=" << v << " lev=" << lev
                          << " | have=" << (unsigned)noccu(j,v)(lev)
                          << " recalc=" << cnt[lev] << std::endl;
              Rcpp::stop("Stage2 check: noccu counts mismatch");
            }
          }
        }
      }
      
      // (E) continuous covariates tallies
      for (unsigned j = 0; j < nmix; ++j) {
        const arma::uvec idx = inds_eq_j_all(j);
        for (unsigned c = 0; c < k_cont; ++c) {
          unsigned n_re = 0U; double s_re = 0.0, ss_re = 0.0;
          for (auto ii : idx) {
            double x = eta_cont(ii,c);
            if (std::isfinite(x)) { ++n_re; s_re += x; ss_re += x*x; }
          }
          
          if ((unsigned)nj_x_curr(j,c) != n_re ||
              std::fabs(s_re  - sum_j_x_curr(j,c)) > 1e-8 ||
              std::fabs(ss_re - ss_j_x_curr(j,c))  > 1e-8) {
            Rcpp::Rcout << "[E] cont mismatch j=" << j << " c=" << c
                        << " | n: have=" << (unsigned)nj_x_curr(j,c) << " recalc=" << n_re
                        << " | sum: have=" << sum_j_x_curr(j,c) << " recalc=" << s_re
                        << " | ss: have=" << ss_j_x_curr(j,c)  << " recalc=" << ss_re << std::endl;
            Rcpp::stop("Stage2 check: continuous tallies mismatch");
          }
        }
      }
    }
    // ================= END HEAVY CONSISTENCY CHECK =================
    // Rcpp::Rcout << "Ending heavy consistency check!" << std::endl;
    
    // ---- 3) draw mu/sig per (j,t) using combined stats ----
    arma::mat mu_draw(nmix, n_trt, fill::zeros), sig_draw(nmix, n_trt, fill::zeros);
    for(unsigned j=0;j<nmix;++j){
      for(unsigned t=0;t<n_trt;++t){
        double ssC, sC; unsigned nC;
        combine_resp(j,t,(unsigned)i, ssC, sC, nC);
        const double mu0_t_cur   = params_t(t, 0);
        const double beta0_t_cur = std::exp(params_t(t, 1));
        arma::vec v = sim_lognorm_params(ssC, sC, df0, a0, mu0_t_cur, beta0_t_cur, nC);
        mu_draw(j,t)  = v(0);
        sig_draw(j,t) = v(1);
      }
    }
    
    // Rcpp::Rcout << "Draw mu/sig(j,t) using combined stats!" << std::endl;
    
    // ---- 4) mixture weights: current cohort with α0/nmix; others from Stage 1 ----
    arma::mat pi(nmix, S, arma::fill::zeros);
    {
      if (A_size == 0u) Rcpp::stop("|A| is zero.");
      
      // Build alpha vector only on allowed clusters
      arma::vec alpha_A(A_size);
      for (unsigned idx = 0; idx < A_size; ++idx) {
        unsigned j = A(idx);
        alpha_A(idx) = nj_val_curr(j,0) + alpha(0)/A_size; // strictly > 0
      }
      
      // Sample Dirichlet on A
      arma::vec tmpA(A_size);
      gsl_ran_dirichlet(r, A_size, alpha_A.begin(), tmpA.begin());
      
      // Scatter back into full-length vector (zeros elsewhere)
      arma::vec tmp(nmix, arma::fill::zeros);
      for (unsigned idx = 0; idx < A_size; ++idx) tmp( A(idx) ) = tmpA(idx);
      
      pi.col(0) = tmp;
      
      // historical cohorts as before
      for (unsigned s = 1; s < S; ++s) pi.col(s) = pimat_old.slice(s-1).row(i).t();
    }
    
    // Rcpp::Rcout << "mixture weights: current cohort with α0/nmix; others from Stage 1!" << std::endl;
    
    // ---- 5) HMC update (mu{0,t},beta{0,t}) using combined outcome tallies ----
    // build combined tallies for HMC (per treatment)
    arma::mat NJ_comb(nmix, n_trt, arma::fill::zeros),
    S_comb (nmix, n_trt, arma::fill::zeros),
    SS_comb(nmix, n_trt, arma::fill::zeros);
    
    for (unsigned j = 0; j < nmix; ++j) {
      for (unsigned t = 0; t < n_trt; ++t) {
        double ssC, sC; unsigned nC;
        combine_resp(j, t, (unsigned)i, ssC, sC, nC);
        NJ_comb(j,t) = nC;
        S_comb (j,t) = sC;
        SS_comb(j,t) = ssC;
      }
    }
    
    // Rcpp::Rcout << "build combined tallies!" << std::endl;
    
    for (unsigned t = 0; t < n_trt; ++t) {
      if (freeze_control && t == 0u) {
        // do NOT update control; params_t(0,.) already set for this draw
        continue;
      }
      
      arma::vec cp = params_t.row(t).t();  // [μ0,t, log β0,t]
      const arma::vec &del_rng_t = (t == 0) ? del_range_lognorm_ref : del_range_lognorm_oth;
      const unsigned nlf_t       = (t == 0) ? nleapfrog_lognorm_ref : nleapfrog_lognorm_oth;
      
      update_lognorm_hyper_extend(
        NJ_comb, S_comb, SS_comb,
        a0, df0,
        mu_m_t(t), mu_v_t(t),
        b_m_t(t),  b_v_t(t),
        del_rng_t, nlf_t,
        cp, acceptance_y(t), t
      );
      params_t(t, 0) = cp(0);
      params_t(t, 1) = cp(1);
    }
    
    // Rcpp::Rcout << "HMC update of treatment-based lognormal parameters!" << std::endl;
    
    // ---- 6) HMC update of α0 (current cohort only) ----
    {
      double l_alpha = std::log(alpha(0));
      arma::uvec ind0 = {0};
      update_alpha(A_size, nj_val_curr(A, ind0), mu_alp, sig_alp,
                   del_range_alp1, nleapfrog_alp1, l_alpha, acceptance_alph(0));
      alpha(0)   = std::exp(l_alpha);
      dir_prec(0)= alpha(0)/A_size;
    }
    
    // Rcpp::Rcout << "HMC update of alpha_0 (current cohort only)!" << std::endl;
    
    // ---- 7) progress printing (unchanged) ----
    if (((i + 1) % thin) == 0) {
      double denom = static_cast<double>(i + 1);
      Rcpp::Rcout << "MCMC iteration: " << i + 1 << std::endl;
      arma::vec acc_y = arma::conv_to<arma::vec>::from(acceptance_y) / denom;
      acc_y.t().print("Treatment Acceptance (lognormal hyper):");
      
      arma::vec acc_alph = arma::conv_to<arma::vec>::from(acceptance_alph) / denom;
      acc_alph.t().print("Dataset Acceptance (alpha):");
    }
    
    // ---- 8) thinning storage ----
    unsigned remainder = (i+1);
    unsigned q = (unsigned)std::floor(remainder/thin);
    remainder -= q*thin;
    if(remainder==0){
      alloc_var_mat.row(q-1)     = del.t();
      dir_alpha_mat.row(q-1)     = alpha.t();
      pimat.slice(0).row(q-1)    = pi.col(0).t();
      for(unsigned s=1;s<S;++s)  pimat.slice(s).row(q-1) = pimat_old.slice(s-1).row(i);
      for(unsigned t=0;t<n_trt;++t){
        lognormal_mu.slice(t).row(q-1)  = mu_draw.col(t).t();
        lognormal_sig.slice(t).row(q-1) = sig_draw.col(t).t();
        hyperparams(q-1, 0, t) = params_t(t, 0);                  // μ0,t
        hyperparams(q-1, 1, t) = std::exp(params_t(t, 1));        // β0,t  (stored on natural scale like Stage 1)
      }
    }
    
    // ---- 9) prediction branch (RESTRICT to historical A for this draw) ----
    if(n_pred){
      field<arma::uvec> non_na_obs_pred(n_pred), non_na_obs_cont_pred(n_pred);
      arma::mat pr_mat(nmix, n_pred, arma::fill::zeros);
      unsigned cat_na_count_pred = 0;
      for(unsigned p=0;p<n_pred;++p){
        non_na_obs_pred(p)      = Rcpp::as<arma::uvec>(non_na_obs1_pred[p]);
        non_na_obs_cont_pred(p) = Rcpp::as<arma::uvec>(non_na_obs1_cont_pred[p]);
        cat_na_count_pred += (k - non_na_obs_pred(p).n_elem);
      }
      arma::uvec del_pred(n_pred);
      arma::mat  y_pred(n_pred, n_trt, fill::zeros);
      
      for(unsigned p=0;p<n_pred;++p){
        arma::vec logp(nmix, fill::value(-std::numeric_limits<double>::infinity()));
        arma::vec cluster_log_prior = arma::log(nj_val_curr.col(0) + alpha(0)/nmix);
        
        // only over allowed set A for draw i
        for(auto j : A){
          double lxc=0.0, lcc=0.0;
          
          for(auto c: non_na_obs_cont_pred(p)){
            double ssXC, sXC; unsigned nXC;
            combine_cont(j, c, (unsigned)i, ssXC, sXC, nXC);
            if(nXC==0) lcc += post_t_dens(eta_cont_pred(p,c), 0.0, 0.0, df_x, alpha_x, mu_x, beta_x, 0);
            else       lcc += post_t_dens(eta_cont_pred(p,c), ssXC, sXC, df_x, alpha_x, mu_x, beta_x, nXC);
          }
          for(auto c: non_na_obs_pred(p)){
            unsigned lev  = (unsigned)eta_pred(p,c);
            unsigned nlev = combine_cat_counts(j, c, lev, (unsigned)i);
            unsigned nobsC= combined_nobs_cat(j, c, (unsigned)i);
            if(nobsC==0) lxc += std::log(1.0 / (double)ncat(c));
            else         lxc += ( std::log((double)nlev + 1.0) - std::log((double)nobsC + (double)ncat(c)) );
          }
          
          logp(j) = cluster_log_prior(j) + lcc + lxc;
        }
        
        double logZ = log_sum_exp(logp);
        arma::vec pr = arma::exp(logp - logZ);
        pr_mat.col(p) = pr;
        arma::uvec d_pred(nmix, fill::zeros);
        gsl_ran_multinomial(r, nmix, 1, pr.begin(), d_pred.begin());
        unsigned jstar = (unsigned)d_pred.index_max();
        del_pred(p) = jstar;
        
        for(unsigned t=0;t<n_trt;++t){
          double ssC, sC; unsigned nC;
          combine_resp(jstar, t, (unsigned)i, ssC, sC, nC);
          const double mu0_t_cur   = params_t(t, 0);
          const double beta0_t_cur = std::exp(params_t(t, 1));
          arma::vec v = sim_lognorm_params(ssC, sC, df0, a0, mu0_t_cur, beta0_t_cur, nC);
          y_pred(p,t) = R::rnorm(v(0), std::sqrt(v(1)));
        }
      }
      
      unsigned remp=(i+1), qp=(unsigned)std::floor(remp/thin); remp-=qp*thin;
      if(remp==0){
        alloc_var_mat_pred.row(qp-1) = del_pred.t();
        for(unsigned t=0;t<n_trt;++t) y_pred_cube.slice(t).row(qp-1) = y_pred.col(t).t();
        for(unsigned p = 0; p < n_pred; ++p){
          pimat_pred.slice(p).row(qp-1) = pr_mat.col(p).t();
        }
      }
    } // end prediction
    // Rcpp::Rcout << "MCMC iteration i: " << i + 1 << std::endl;
  }   // end MCMC
  
  arma::vec acceptance_all(S + n_trt, arma::fill::zeros);
  acceptance_all.subvec(0, S-1) = 
    arma::conv_to<arma::vec>::from(acceptance_alph) / ((double)(nrun+burn));
  acceptance_all.subvec(S, S + n_trt - 1) =
    arma::conv_to<arma::vec>::from(acceptance_y) / ((double)(nrun+burn));
  
  gsl_rng_free(r);
  
  if(n_pred){
    return Rcpp::List::create(
      Rcpp::Named("picube")                 = pimat,
      Rcpp::Named("Lognormal_Mu_Cube")      = lognormal_mu,
      Rcpp::Named("Lognormal_Sig_Cube")     = lognormal_sig,
      Rcpp::Named("Lognormal_hyperparams")  = hyperparams,
      Rcpp::Named("Dirichlet_params")       = dir_alpha_mat,
      Rcpp::Named("Acceptance_rates")       = acceptance_all,
      Rcpp::Named("Allocation_variables")   = alloc_var_mat,
      Rcpp::Named("New_Categorical_Covariates") = eta_pred,
      Rcpp::Named("New_Continuous_Covariates")  = eta_cont_pred,
      Rcpp::Named("Predicted_Allocation_variables") = alloc_var_mat_pred,
      Rcpp::Named("Predicted_Y")            = arma::exp(y_pred_cube)
    );
  } else {
    return Rcpp::List::create(
      Rcpp::Named("picube")                 = pimat,
      Rcpp::Named("Lognormal_Mu_Cube")      = lognormal_mu,
      Rcpp::Named("Lognormal_Sig_Cube")     = lognormal_sig,
      Rcpp::Named("Lognormal_hyperparams")  = hyperparams,
      Rcpp::Named("Dirichlet_params")       = dir_alpha_mat,
      Rcpp::Named("Acceptance_rates")       = acceptance_all,
      Rcpp::Named("Allocation_variables")   = alloc_var_mat
    );
  }
}
