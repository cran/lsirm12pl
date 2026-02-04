#include "progress.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

static inline double inv_logit_stable(const double x) {
  if (x >= 0.0) {
    return 1.0 / (1.0 + std::exp(-x));
  }
  const double ex = std::exp(x);
  return ex / (1.0 + ex);
}

template <typename RowLike>
static inline double log_grm_prob_y0based(const int y, const RowLike &beta_i,
                                          const double eta) {
  // Graded response model with latent-space distance effect:
  // P(Y >= k) = logistic(eta + beta_{i,k})  for k = 1,...,K-1
  // where beta_i has length (K-1) and corresponds to thresholds k=1..K-1.
  // y is 0-based category in {0,...,K-1}.
  const int Kminus1 = (int)beta_i.n_elem;
  const int K = Kminus1 + 1;

  // Use larger epsilon for better numerical stability
  // log(1e-10) ≈ -23, which is safer than log(1e-16) ≈ -37
  const double eps = 1e-10;

  auto p_ge = [&](const int k) -> double {
    // k in 1..K-1
    return inv_logit_stable(eta + beta_i((arma::uword)(k - 1)));
  };

  double p = 0.0;
  if (y <= 0) {
    p = 1.0 - p_ge(1);
  } else if (y >= (K - 1)) {
    p = p_ge(K - 1);
  } else {
    // y in 1..K-2
    p = p_ge(y) - p_ge(y + 1);
  }

  // Clamp probability to avoid log(0) or log(negative)
  if (!std::isfinite(p) || p <= eps)
    p = eps;
  if (p >= 1.0 - eps)
    p = 1.0 - eps;
  return std::log(p);
}

static inline double log_prior_gamma_lognormal(const double gamma,
                                               const double meanlog,
                                               const double sdlog) {
  return R::dlnorm(gamma, meanlog, sdlog, 1);
}

struct AdaptControlGrm {
  bool use_adapt = false;
  int adapt_interval = 100;
  double adapt_rate = 1.0;
  double decay_rate = 0.5;
  double jump_min = 1e-4;
  double jump_max = 10.0;
  double target_accept_beta = 0.44;
  double target_accept_theta = 0.44;
  double target_accept_gamma = 0.44;
  double target_accept_alpha = 0.44;
  double target_accept_zw = 0.234;
};

static inline double clamp_double(const double x, const double lo,
                                  const double hi) {
  if (x < lo)
    return lo;
  if (x > hi)
    return hi;
  return x;
}

static inline bool list_has(const Rcpp::List &x, const char *name) {
  return x.containsElementNamed(name) &&
         !Rcpp::as<Rcpp::RObject>(x[name]).isNULL();
}

static inline double list_double(const Rcpp::List &x, const char *name,
                                 const double def) {
  if (!list_has(x, name))
    return def;
  return Rcpp::as<double>(x[name]);
}

static inline int list_int(const Rcpp::List &x, const char *name,
                           const int def) {
  if (!list_has(x, name))
    return def;
  return Rcpp::as<int>(x[name]);
}

static inline bool list_bool(const Rcpp::List &x, const char *name,
                             const bool def) {
  if (!list_has(x, name))
    return def;
  return Rcpp::as<bool>(x[name]);
}

static inline AdaptControlGrm
parse_adapt_grm(const Rcpp::Nullable<Rcpp::List> &adapt, const int ndim,
                const bool fixed_gamma, const bool use_alpha) {
  AdaptControlGrm ac;
  if (adapt.isNull())
    return ac;

  const Rcpp::List a(adapt);
  ac.use_adapt = list_bool(a, "use_adapt", ac.use_adapt);
  ac.adapt_interval = list_int(a, "adapt_interval", ac.adapt_interval);
  ac.adapt_rate = list_double(a, "adapt_rate", ac.adapt_rate);
  ac.decay_rate = list_double(a, "decay_rate", ac.decay_rate);
  ac.jump_min = list_double(a, "jump_min", ac.jump_min);
  ac.jump_max = list_double(a, "jump_max", ac.jump_max);

  const double target_scalar = list_double(a, "target_accept", 0.44);
  ac.target_accept_beta = list_double(a, "target_accept_beta", target_scalar);
  ac.target_accept_theta = list_double(a, "target_accept_theta", target_scalar);
  ac.target_accept_gamma = list_double(a, "target_accept_gamma", target_scalar);
  ac.target_accept_alpha = list_double(a, "target_accept_alpha", target_scalar);
  ac.target_accept_zw = list_double(a, "target_accept_zw", ac.target_accept_zw);

  if (!(ac.adapt_interval > 0))
    ac.use_adapt = false;
  ac.adapt_rate = (ac.adapt_rate > 0.0) ? ac.adapt_rate : 1.0;
  if (!(ac.jump_min > 0.0))
    ac.jump_min = 1e-8;
  if (!(ac.jump_max > ac.jump_min))
    ac.jump_max = ac.jump_min * 1e4;
  ac.target_accept_beta = clamp_double(ac.target_accept_beta, 1e-4, 0.999);
  ac.target_accept_theta = clamp_double(ac.target_accept_theta, 1e-4, 0.999);
  ac.target_accept_gamma = clamp_double(ac.target_accept_gamma, 1e-4, 0.999);
  ac.target_accept_alpha = clamp_double(ac.target_accept_alpha, 1e-4, 0.999);
  ac.target_accept_zw = clamp_double(ac.target_accept_zw, 1e-4, 0.999);

  if (fixed_gamma)
    ac.target_accept_gamma = NA_REAL;
  if (!use_alpha)
    ac.target_accept_alpha = NA_REAL;
  return ac;
}

static Rcpp::List lsirmgrm_internal(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta, const double pr_a_theta,
    const double pr_b_theta, const double pr_mean_theta, double pr_sd_theta,
    const double pr_spike_mean, const double pr_spike_sd,
    const double pr_slab_mean, const double pr_slab_sd, const double pr_beta_a,
    const double pr_beta_b, const double jump_gamma, const double missing,
    const bool verbose, const bool fix_theta_sd, const bool fixed_gamma,
    const bool spike_slab, const bool use_alpha, const double jump_alpha,
    const double pr_mean_alpha, const double pr_sd_alpha,
    const bool fix_alpha_1, const bool missing_mar,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {

  int i, j, k, c, count, accept;
  double num, den, old_like_beta, new_like_beta, old_like_theta, new_like_theta;
  double old_like_z, new_like_z, old_like_w, new_like_w, old_like_gamma,
      new_like_gamma;
  double old_like_alpha = 0.0, new_like_alpha = 0.0;
  double ratio, un, post_a, post_b, dist_temp, dist_new_temp;
  double pr_mean_z = 0.0, pr_sd_z = 1.0, pr_mean_w = 0.0, pr_sd_w = 1.0, mle;

  const int nsample = data.n_rows;
  const int nitem = data.n_cols;
  const int K = ncat;
  const int Kminus1 = K - 1;

  if (K < 2) {
    Rcpp::stop("ncat must be >= 2.");
  }

  const AdaptControlGrm adapt_ctrl =
      parse_adapt_grm(adapt, ndim, fixed_gamma, use_alpha);
  const bool do_adapt = adapt_ctrl.use_adapt && (nburn > 0);

  double cur_jump_beta = clamp_double((jump_beta > 0.0 ? jump_beta : 1e-8),
                                      adapt_ctrl.jump_min, adapt_ctrl.jump_max);
  double cur_jump_theta =
      clamp_double((jump_theta > 0.0 ? jump_theta : 1e-8), adapt_ctrl.jump_min,
                   adapt_ctrl.jump_max);
  double cur_jump_z = clamp_double((jump_z > 0.0 ? jump_z : 1e-8),
                                   adapt_ctrl.jump_min, adapt_ctrl.jump_max);
  double cur_jump_w = clamp_double((jump_w > 0.0 ? jump_w : 1e-8),
                                   adapt_ctrl.jump_min, adapt_ctrl.jump_max);
  double cur_jump_gamma =
      clamp_double((jump_gamma > 0.0 ? jump_gamma : 1e-8), adapt_ctrl.jump_min,
                   adapt_ctrl.jump_max);
  double cur_jump_alpha =
      clamp_double((jump_alpha > 0.0 ? jump_alpha : 1e-8), adapt_ctrl.jump_min,
                   adapt_ctrl.jump_max);

  double log_jump_beta = std::log(cur_jump_beta);
  double log_jump_theta = std::log(cur_jump_theta);
  double log_jump_z = std::log(cur_jump_z);
  double log_jump_w = std::log(cur_jump_w);
  double log_jump_gamma = std::log(cur_jump_gamma);
  double log_jump_alpha = std::log(cur_jump_alpha);

  double win_acc_beta = 0.0, win_prop_beta = 0.0;
  double win_acc_theta = 0.0, win_prop_theta = 0.0;
  double win_acc_gamma = 0.0, win_prop_gamma = 0.0;
  double win_acc_alpha = 0.0, win_prop_alpha = 0.0;
  double win_acc_z = 0.0, win_prop_z = 0.0;
  double win_acc_w = 0.0, win_prop_w = 0.0;

  double burn_acc_beta = 0.0, burn_prop_beta = 0.0;
  double burn_acc_theta = 0.0, burn_prop_theta = 0.0;
  double burn_acc_gamma = 0.0, burn_prop_gamma = 0.0;
  double burn_acc_alpha = 0.0, burn_prop_alpha = 0.0;
  double burn_acc_z = 0.0, burn_prop_z = 0.0;
  double burn_acc_w = 0.0, burn_prop_w = 0.0;

  // Acceptance rate in the most recent *completed* adaptation window.
  double lastwin_accept_beta = NA_REAL;
  double lastwin_accept_theta = NA_REAL;
  double lastwin_accept_gamma = NA_REAL;
  double lastwin_accept_alpha = NA_REAL;
  double lastwin_accept_z = NA_REAL;
  double lastwin_accept_w = NA_REAL;

  int adapt_step = 0;

  // detect if data is 1-based (Likert 1..K) or 0-based (0..K-1)
  const double data_min = data.min();
  const int data_is_one_based = (data_min >= 1.0);

  // Precompute integer responses (0-based) and missing mask to avoid repeated
  // casts/clamping.
  arma::imat yint(nsample, nitem, fill::value(-1));
  int nmissing = 0;
  arma::ivec miss_row;
  arma::ivec miss_col;
  for (int rr = 0; rr < nsample; rr++) {
    for (int cc = 0; cc < nitem; cc++) {
      if (data(rr, cc) == missing)
        continue;
      int y = (int)(data(rr, cc));
      if (data_is_one_based)
        y -= 1;
      if (y < 0)
        y = 0;
      if (y > (K - 1))
        y = K - 1;
      yint(rr, cc) = y;
    }
  }

  if (missing_mar) {
    // record missing positions for data augmentation
    nmissing = 0;
    for (int rr = 0; rr < nsample; rr++) {
      for (int cc = 0; cc < nitem; cc++) {
        if (data(rr, cc) == missing)
          nmissing++;
      }
    }
    miss_row = arma::ivec(nmissing);
    miss_col = arma::ivec(nmissing);
    int mi = 0;
    for (int rr = 0; rr < nsample; rr++) {
      for (int cc = 0; cc < nitem; cc++) {
        if (data(rr, cc) == missing) {
          miss_row(mi) = rr;
          miss_col(mi) = cc;
          mi++;
        }
      }
    }
  }

  arma::dmat oldbeta(nitem, Kminus1, fill::randu);
  oldbeta = oldbeta * 4.0 - 2.0;
  // enforce descending order per item
  for (i = 0; i < nitem; i++) {
    arma::rowvec bi = oldbeta.row(i);
    bi = sort(bi, "descend");
    oldbeta.row(i) = bi;
  }

  arma::dvec oldtheta(nsample, fill::randu);
  oldtheta = oldtheta * 4.0 - 2.0;
  arma::dvec newtheta = oldtheta;

  arma::dmat oldz(nsample, ndim, fill::randu);
  oldz = oldz * 2.0 - 1.0;
  arma::dmat newz = oldz;

  arma::dmat oldw(nitem, ndim, fill::randu);
  oldw = oldw * 2.0 - 1.0;
  arma::dmat neww = oldw;

  arma::dvec oldalpha;
  arma::dvec newalpha;
  if (use_alpha) {
    oldalpha = arma::dvec(nitem, fill::randu);
    oldalpha = oldalpha + 0.5;
    if (fix_alpha_1)
      oldalpha(0) = 1.0;
    newalpha = oldalpha;
  }

  double oldgamma = 1.0, newgamma = 1.0;
  if (fixed_gamma) {
    oldgamma = 1.0;
    newgamma = 1.0;
  }

  int pi_gamma = 1;
  double xi = 0.5;
  if (spike_slab) {
    // match lsirm1pl_ss initialization
    pi_gamma = (int)R::rbinom(1.0, 0.5);
    xi = R::rbeta(pr_beta_a + (double)pi_gamma,
                  pr_beta_b + 1.0 - (double)pi_gamma);
  }

  const int nmcmc = (niter - nburn) / nthin;
  arma::dcube samp_beta(nmcmc, nitem, Kminus1, fill::zeros);
  arma::dmat samp_theta(nmcmc, nsample, fill::zeros);
  arma::dcube samp_z(nmcmc, nsample, ndim, fill::zeros);
  arma::dcube samp_w(nmcmc, nitem, ndim, fill::zeros);
  arma::dvec samp_sd_theta(nmcmc, fill::zeros);
  arma::dvec sample_mle(nmcmc, fill::zeros);
  arma::dvec samp_gamma(nmcmc, fill::zeros);
  arma::dmat samp_alpha;
  if (use_alpha) {
    samp_alpha = arma::dmat(nmcmc, nitem, fill::zeros);
  }
  arma::ivec samp_pi;
  arma::dvec samp_xi;
  if (spike_slab) {
    samp_pi = arma::ivec(nmcmc, fill::zeros);
    samp_xi = arma::dvec(nmcmc, fill::zeros);
  }

  arma::dmat samp_impute;
  arma::dvec impute_col;
  if (missing_mar) {
    samp_impute = arma::dmat(nmcmc, nmissing, fill::zeros);
    impute_col = arma::dvec(nmissing, fill::zeros);
  }

  arma::dmat accept_beta(nitem, Kminus1, fill::zeros);
  arma::dvec accept_theta(nsample, fill::zeros);
  arma::dvec accept_z(nsample, fill::zeros);
  arma::dvec accept_w(nitem, fill::zeros);
  double accept_gamma = 0.0;
  arma::dvec accept_alpha;
  if (use_alpha) {
    accept_alpha = arma::dvec(nitem, fill::zeros);
  }

  accept = count = 0;

  // Cache distances and update incrementally on accepted z/w moves.
  arma::dmat dist(nsample, nitem, fill::zeros);
  arma::dvec old_dist_k(nitem, fill::zeros);
  arma::dvec new_dist_k(nitem, fill::zeros);
  arma::dvec old_dist_i(nsample, fill::zeros);
  arma::dvec new_dist_i(nsample, fill::zeros);

  // Initial distance matrix
  for (i = 0; i < nitem; i++) {
    for (k = 0; k < nsample; k++) {
      dist_temp = 0.0;
      for (j = 0; j < ndim; j++) {
        const double dkj = oldz(k, j) - oldw(i, j);
        dist_temp += dkj * dkj;
      }
      dist(k, i) = std::sqrt(dist_temp);
    }
  }

  // helper for sampling an ordinal category from GRM probabilities (0-based)
  auto sample_grm_y = [&](const arma::rowvec &beta_i,
                          const double eta_i) -> int {
    arma::vec logp(K, fill::zeros);
    double maxlp = -std::numeric_limits<double>::infinity();
    for (int yy = 0; yy < K; yy++) {
      logp(yy) = log_grm_prob_y0based(yy, beta_i, eta_i);
      if (logp(yy) > maxlp)
        maxlp = logp(yy);
    }
    arma::vec p(K, fill::zeros);
    double s = 0.0;
    for (int yy = 0; yy < K; yy++) {
      p(yy) = std::exp(logp(yy) - maxlp);
      s += p(yy);
    }
    if (!(s > 0.0) || !std::isfinite(s)) {
      // fallback: uniform
      return (int)std::floor(R::runif(0.0, (double)K));
    }
    p /= s;
    const double u = R::runif(0.0, 1.0);
    double cdf = 0.0;
    for (int yy = 0; yy < K; yy++) {
      cdf += p(yy);
      if (u <= cdf)
        return yy;
    }
    return K - 1;
  };

  // Reused temporaries to reduce loop-local declarations/allocations.
  arma::rowvec beta_old_i;
  arma::rowvec beta_new_i;
  double eta = 0.0;
  double eta_new = 0.0;
  double eta_old = 0.0;
  int y = -1;

  for (int iter = 0; iter < niter; iter++) {
    if (iter % 10 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Imputation step (MAR data augmentation)
    if (missing_mar && nmissing > 0) {
      for (int mi = 0; mi < nmissing; mi++) {
        const int rr = miss_row(mi);
        const int cc = miss_col(mi);
        const double theta_part =
            use_alpha ? (oldalpha(cc) * oldtheta(rr)) : oldtheta(rr);
        const double eta_imp = theta_part - oldgamma * dist(rr, cc);
        const arma::rowvec beta_i = oldbeta.row(cc);
        const int ynew = sample_grm_y(beta_i, eta_imp);
        yint(rr, cc) = ynew;
        // store on the original data scale (0..K-1 or 1..K)
        impute_col(mi) = (double)(ynew + (data_is_one_based ? 1 : 0));
      }
    }

    // beta(thresholds) update
    // Individual update with neighbor constraint check
    for (i = 0; i < nitem; i++) {
      beta_old_i = oldbeta.row(i);
      
      // Update each threshold individually
      for (c = 0; c < Kminus1; c++) {
        const double beta_proposal = R::rnorm(beta_old_i((uword)c), cur_jump_beta);

        if (do_adapt && iter < nburn) {
          win_prop_beta += 1.0;
          burn_prop_beta += 1.0;
        }

        // Check ordering constraint with neighbors
        // For descending order: beta[c-1] > beta[c] > beta[c+1]
        bool accept_constraint = true;
        if (c > 0 && !(beta_old_i((uword)(c - 1)) > beta_proposal)) {
          accept_constraint = false;
        }
        if (c < Kminus1 - 1 && !(beta_proposal > beta_old_i((uword)(c + 1)))) {
          accept_constraint = false;
        }
        
        // Skip MH step if constraint violated
        if (!accept_constraint) {
          continue;
        }

        // MH accept/reject step
        // Compute likelihood with proposed value
        beta_new_i = beta_old_i;
        beta_new_i((uword)c) = beta_proposal;
        
        old_like_beta = new_like_beta = 0.0;
        for (k = 0; k < nsample; k++) {
          y = yint(k, i);
          if (y < 0)
            continue;

          const double theta_part =
              use_alpha ? (oldalpha(i) * oldtheta(k)) : oldtheta(k);
          eta = theta_part - oldgamma * dist(k, i);
          new_like_beta += log_grm_prob_y0based(y, beta_new_i, eta);
          old_like_beta += log_grm_prob_y0based(y, beta_old_i, eta);
        }

        num = new_like_beta + R::dnorm4(beta_proposal, pr_mean_beta, pr_sd_beta, 1);
        den = old_like_beta + R::dnorm4(beta_old_i((uword)c), pr_mean_beta, pr_sd_beta, 1);
        ratio = num - den;

        if (ratio > 0.0)
          accept = 1;
        else {
          un = R::runif(0, 1);
          if (std::log(un) < ratio)
            accept = 1;
          else
            accept = 0;
        }

        if (accept == 1) {
          beta_old_i((uword)c) = beta_proposal;
          if (!do_adapt || iter >= nburn) {
            accept_beta(i, c) += 1.0 / ((niter - nburn) * 1.0);
          }
          if (do_adapt && iter < nburn) {
            win_acc_beta += 1.0;
            burn_acc_beta += 1.0;
          }
        }
      }
      
      // Save updated thresholds back to oldbeta
      oldbeta.row(i) = beta_old_i;
    }

    // theta update
    for (k = 0; k < nsample; k++) {
      newtheta(k) = R::rnorm(oldtheta(k), cur_jump_theta);
      if (do_adapt && iter < nburn) {
        win_prop_theta += 1.0;
        burn_prop_theta += 1.0;
      }
      old_like_theta = new_like_theta = 0.0;

      for (i = 0; i < nitem; i++) {
        y = yint(k, i);
        if (y < 0)
          continue;

        const double alpha_i = use_alpha ? oldalpha(i) : 1.0;
        eta_new = alpha_i * newtheta(k) - oldgamma * dist(k, i);
        eta_old = alpha_i * oldtheta(k) - oldgamma * dist(k, i);
        new_like_theta += log_grm_prob_y0based(y, oldbeta.row(i), eta_new);
        old_like_theta += log_grm_prob_y0based(y, oldbeta.row(i), eta_old);
      }

      num = new_like_theta +
            R::dnorm4(newtheta(k), pr_mean_theta, pr_sd_theta, 1);
      den = old_like_theta +
            R::dnorm4(oldtheta(k), pr_mean_theta, pr_sd_theta, 1);
      ratio = num - den;

      if (ratio > 0.0)
        accept = 1;
      else {
        un = R::runif(0, 1);
        if (std::log(un) < ratio)
          accept = 1;
        else
          accept = 0;
      }

      if (accept == 1) {
        oldtheta(k) = newtheta(k);
        if (!do_adapt || iter >= nburn) {
          accept_theta(k) += 1.0 / ((niter - nburn) * 1.0);
        }
        if (do_adapt && iter < nburn) {
          win_acc_theta += 1.0;
          burn_acc_theta += 1.0;
        }
      } else {
        newtheta(k) = oldtheta(k);
      }
    }

    // gamma update
    if (!fixed_gamma) {
      newgamma = R::rlnorm(std::log(oldgamma), cur_jump_gamma);
      if (do_adapt && iter < nburn) {
        win_prop_gamma += 1.0;
        burn_prop_gamma += 1.0;
      }
      old_like_gamma = new_like_gamma = 0.0;

      for (k = 0; k < nsample; k++) {
        for (i = 0; i < nitem; i++) {
          y = yint(k, i);
          if (y < 0)
            continue;

          const double theta_part =
              use_alpha ? (oldalpha(i) * oldtheta(k)) : oldtheta(k);
          eta_new = theta_part - newgamma * dist(k, i);
          eta_old = theta_part - oldgamma * dist(k, i);
          new_like_gamma += log_grm_prob_y0based(y, oldbeta.row(i), eta_new);
          old_like_gamma += log_grm_prob_y0based(y, oldbeta.row(i), eta_old);
        }
      }

      // proposal + prior
      num = new_like_gamma +
            R::dlnorm(oldgamma, std::log(newgamma), cur_jump_gamma, 1);
      den = old_like_gamma +
            R::dlnorm(newgamma, std::log(oldgamma), cur_jump_gamma, 1);
      if (spike_slab) {
        if (pi_gamma == 1) {
          num += log_prior_gamma_lognormal(newgamma, pr_slab_mean, pr_slab_sd);
          den += log_prior_gamma_lognormal(oldgamma, pr_slab_mean, pr_slab_sd);
        } else {
          num +=
              log_prior_gamma_lognormal(newgamma, pr_spike_mean, pr_spike_sd);
          den +=
              log_prior_gamma_lognormal(oldgamma, pr_spike_mean, pr_spike_sd);
        }
      } else {
        num += log_prior_gamma_lognormal(newgamma, pr_slab_mean, pr_slab_sd);
        den += log_prior_gamma_lognormal(oldgamma, pr_slab_mean, pr_slab_sd);
      }
      ratio = num - den;

      if (ratio > 0.0)
        accept = 1;
      else {
        un = R::runif(0, 1);
        if (std::log(un) < ratio)
          accept = 1;
        else
          accept = 0;
      }

      if (accept == 1) {
        oldgamma = newgamma;
        if (!do_adapt || iter >= nburn) {
          accept_gamma += 1.0 / ((niter - nburn) * 1.0);
        }
        if (do_adapt && iter < nburn) {
          win_acc_gamma += 1.0;
          burn_acc_gamma += 1.0;
        }
      } else {
        newgamma = oldgamma;
      }

      if (spike_slab) {
        // pi update (indicator)
        double log_d_slab = R::dlnorm(oldgamma, pr_slab_mean, pr_slab_sd, 1);
        double log_d_spike = R::dlnorm(oldgamma, pr_spike_mean, pr_spike_sd, 1);

        double log_num_pi, log_den_pi;
        if (xi <= 0.0) {
          log_num_pi = -1e300; // Effectively -Inf
        } else {
          log_num_pi = std::log(xi) + log_d_slab;
        }

        double log_term_spike;
        if (xi >= 1.0) {
          log_term_spike = -1e300; // Effectively -Inf
        } else {
          log_term_spike = std::log(1.0 - xi) + log_d_spike;
        }

        double max_val = std::max(log_term_spike, log_num_pi);
        log_den_pi = max_val + std::log(std::exp(log_term_spike - max_val) +
                                        std::exp(log_num_pi - max_val));

        double prob_pi = std::exp(log_num_pi - log_den_pi);
        if (std::isnan(prob_pi))
          prob_pi = 0.0;

        pi_gamma = (int)R::rbinom(1.0, prob_pi);

        // xi update
        xi = R::rbeta(pr_beta_a + (double)pi_gamma,
                      pr_beta_b + 1.0 - (double)pi_gamma);
      }
    }

    // alpha update (2PL)
    if (use_alpha) {
      for (i = (fix_alpha_1 ? 1 : 0); i < nitem; i++) {
        newalpha(i) = R::rlnorm(std::log(oldalpha(i)), cur_jump_alpha);
        if (do_adapt && iter < nburn) {
          win_prop_alpha += 1.0;
          burn_prop_alpha += 1.0;
        }
        old_like_alpha = new_like_alpha = 0.0;

        for (k = 0; k < nsample; k++) {
          y = yint(k, i);
          if (y < 0)
            continue;

          const double eta_new_a =
              newalpha(i) * oldtheta(k) - oldgamma * dist(k, i);
          const double eta_old_a =
              oldalpha(i) * oldtheta(k) - oldgamma * dist(k, i);
          new_like_alpha += log_grm_prob_y0based(y, oldbeta.row(i), eta_new_a);
          old_like_alpha += log_grm_prob_y0based(y, oldbeta.row(i), eta_old_a);
        }

        num = new_like_alpha +
              R::dlnorm(oldalpha(i), std::log(newalpha(i)), cur_jump_alpha, 1) +
              R::dlnorm(newalpha(i), pr_mean_alpha, pr_sd_alpha, 1);
        den = old_like_alpha +
              R::dlnorm(newalpha(i), std::log(oldalpha(i)), cur_jump_alpha, 1) +
              R::dlnorm(oldalpha(i), pr_mean_alpha, pr_sd_alpha, 1);
        ratio = num - den;

        if (ratio > 0.0)
          accept = 1;
        else {
          un = R::runif(0, 1);
          if (std::log(un) < ratio)
            accept = 1;
          else
            accept = 0;
        }

        if (accept == 1) {
          oldalpha(i) = newalpha(i);
          if (!do_adapt || iter >= nburn) {
            accept_alpha(i) += 1.0 / ((niter - nburn) * 1.0);
          }
          if (do_adapt && iter < nburn) {
            win_acc_alpha += 1.0;
            burn_acc_alpha += 1.0;
          }
        } else {
          newalpha(i) = oldalpha(i);
        }
      }
    }

    // z update
    for (k = 0; k < nsample; k++) {
      for (j = 0; j < ndim; j++)
        newz(k, j) = R::rnorm(oldz(k, j), cur_jump_z);
      if (do_adapt && iter < nburn) {
        win_prop_z += 1.0;
        burn_prop_z += 1.0;
      }
      old_like_z = new_like_z = 0.0;

      old_dist_k = dist.row(k).t();
      for (i = 0; i < nitem; i++) {
        dist_new_temp = 0.0;
        for (j = 0; j < ndim; j++) {
          const double dkj = newz(k, j) - oldw(i, j);
          dist_new_temp += dkj * dkj;
        }
        new_dist_k(i) = std::sqrt(dist_new_temp);
      }

      for (i = 0; i < nitem; i++) {
        y = yint(k, i);
        if (y < 0)
          continue;

        const double theta_part =
            use_alpha ? (oldalpha(i) * oldtheta(k)) : oldtheta(k);
        eta_new = theta_part - oldgamma * new_dist_k(i);
        eta_old = theta_part - oldgamma * old_dist_k(i);
        new_like_z += log_grm_prob_y0based(y, oldbeta.row(i), eta_new);
        old_like_z += log_grm_prob_y0based(y, oldbeta.row(i), eta_old);
      }

      num = den = 0.0;
      for (j = 0; j < ndim; j++) {
        num += R::dnorm4(newz(k, j), pr_mean_z, pr_sd_z, 1);
        den += R::dnorm4(oldz(k, j), pr_mean_z, pr_sd_z, 1);
      }
      num += new_like_z;
      den += old_like_z;
      ratio = num - den;

      if (ratio > 0.0)
        accept = 1;
      else {
        un = R::runif(0, 1);
        if (std::log(un) < ratio)
          accept = 1;
        else
          accept = 0;
      }

      if (accept == 1) {
        for (j = 0; j < ndim; j++)
          oldz(k, j) = newz(k, j);
        dist.row(k) = new_dist_k.t();
        if (!do_adapt || iter >= nburn) {
          accept_z(k) += 1.0 / ((niter - nburn) * 1.0);
        }
        if (do_adapt && iter < nburn) {
          win_acc_z += 1.0;
          burn_acc_z += 1.0;
        }
      } else {
        for (j = 0; j < ndim; j++)
          newz(k, j) = oldz(k, j);
      }
    }

    // w update
    for (i = 0; i < nitem; i++) {
      for (j = 0; j < ndim; j++)
        neww(i, j) = R::rnorm(oldw(i, j), cur_jump_w);
      if (do_adapt && iter < nburn) {
        win_prop_w += 1.0;
        burn_prop_w += 1.0;
      }
      old_like_w = new_like_w = 0.0;

      old_dist_i = dist.col(i);
      for (k = 0; k < nsample; k++) {
        dist_new_temp = 0.0;
        for (j = 0; j < ndim; j++) {
          const double dkj = oldz(k, j) - neww(i, j);
          dist_new_temp += dkj * dkj;
        }
        new_dist_i(k) = std::sqrt(dist_new_temp);
      }

      for (k = 0; k < nsample; k++) {
        y = yint(k, i);
        if (y < 0)
          continue;

        const double theta_part =
            use_alpha ? (oldalpha(i) * oldtheta(k)) : oldtheta(k);
        eta_new = theta_part - oldgamma * new_dist_i(k);
        eta_old = theta_part - oldgamma * old_dist_i(k);
        new_like_w += log_grm_prob_y0based(y, oldbeta.row(i), eta_new);
        old_like_w += log_grm_prob_y0based(y, oldbeta.row(i), eta_old);
      }

      num = den = 0.0;
      for (j = 0; j < ndim; j++) {
        num += R::dnorm4(neww(i, j), pr_mean_w, pr_sd_w, 1);
        den += R::dnorm4(oldw(i, j), pr_mean_w, pr_sd_w, 1);
      }
      num += new_like_w;
      den += old_like_w;
      ratio = num - den;

      if (ratio > 0.0)
        accept = 1;
      else {
        un = R::runif(0, 1);
        if (std::log(un) < ratio)
          accept = 1;
        else
          accept = 0;
      }

      if (accept == 1) {
        for (j = 0; j < ndim; j++)
          oldw(i, j) = neww(i, j);
        dist.col(i) = new_dist_i;
        if (!do_adapt || iter >= nburn) {
          accept_w(i) += 1.0 / ((niter - nburn) * 1.0);
        }
        if (do_adapt && iter < nburn) {
          win_acc_w += 1.0;
          burn_acc_w += 1.0;
        }
      } else {
        for (j = 0; j < ndim; j++)
          neww(i, j) = oldw(i, j);
      }
    }

    // sigma_theta update (Gibbs)
    if (!fix_theta_sd) {
      post_a = 2 * pr_a_theta + nsample;
      post_b = pr_b_theta;
      for (j = 0; j < nsample; j++)
        post_b += std::pow((oldtheta(j) - pr_mean_theta), 2.0) / 2.0;
      pr_sd_theta = std::sqrt(2 * post_b * (1.0 / R::rchisq(post_a)));
    }

    // Adapt proposal SDs during burn-in only (Robbins--Monro on log scale)
    if (do_adapt && iter < nburn &&
        ((iter + 1) % adapt_ctrl.adapt_interval == 0)) {
      const double eta_rm =
          adapt_ctrl.adapt_rate /
          std::pow((double)adapt_step + 1.0, adapt_ctrl.decay_rate);

      // record acceptance rates for diagnostics at the end of this window
      lastwin_accept_beta =
          (win_prop_beta > 0.0) ? (win_acc_beta / win_prop_beta) : NA_REAL;
      lastwin_accept_theta =
          (win_prop_theta > 0.0) ? (win_acc_theta / win_prop_theta) : NA_REAL;
      lastwin_accept_gamma = (!fixed_gamma && win_prop_gamma > 0.0)
                                 ? (win_acc_gamma / win_prop_gamma)
                                 : NA_REAL;
      lastwin_accept_alpha = (use_alpha && win_prop_alpha > 0.0)
                                 ? (win_acc_alpha / win_prop_alpha)
                                 : NA_REAL;
      lastwin_accept_z =
          (win_prop_z > 0.0) ? (win_acc_z / win_prop_z) : NA_REAL;
      lastwin_accept_w =
          (win_prop_w > 0.0) ? (win_acc_w / win_prop_w) : NA_REAL;

      auto update_jump = [&](double &log_jump, double &cur_jump,
                             const double acc, const double prop,
                             const double target) {
        if (!(prop > 0.0) || !std::isfinite(target))
          return;
        const double acc_rate = acc / prop;
        log_jump += eta_rm * (acc_rate - target);
        cur_jump = clamp_double(std::exp(log_jump), adapt_ctrl.jump_min,
                                adapt_ctrl.jump_max);
        log_jump = std::log(cur_jump);
      };

      update_jump(log_jump_beta, cur_jump_beta, win_acc_beta, win_prop_beta,
                  adapt_ctrl.target_accept_beta);
      update_jump(log_jump_theta, cur_jump_theta, win_acc_theta, win_prop_theta,
                  adapt_ctrl.target_accept_theta);
      if (!fixed_gamma)
        update_jump(log_jump_gamma, cur_jump_gamma, win_acc_gamma,
                    win_prop_gamma, adapt_ctrl.target_accept_gamma);
      if (use_alpha)
        update_jump(log_jump_alpha, cur_jump_alpha, win_acc_alpha,
                    win_prop_alpha, adapt_ctrl.target_accept_alpha);
      update_jump(log_jump_z, cur_jump_z, win_acc_z, win_prop_z,
                  adapt_ctrl.target_accept_zw);
      update_jump(log_jump_w, cur_jump_w, win_acc_w, win_prop_w,
                  adapt_ctrl.target_accept_zw);

      win_acc_beta = win_prop_beta = 0.0;
      win_acc_theta = win_prop_theta = 0.0;
      win_acc_gamma = win_prop_gamma = 0.0;
      win_acc_alpha = win_prop_alpha = 0.0;
      win_acc_z = win_prop_z = 0.0;
      win_acc_w = win_prop_w = 0.0;
      adapt_step++;
    }

    if (iter >= nburn && iter % nthin == 0) {
      for (i = 0; i < nitem; i++) {
        for (c = 0; c < Kminus1; c++) {
          samp_beta(count, i, c) = oldbeta(i, c);
        }
      }
      for (k = 0; k < nsample; k++)
        samp_theta(count, k) = oldtheta(k);
      if (use_alpha) {
        for (i = 0; i < nitem; i++)
          samp_alpha(count, i) = oldalpha(i);
      }
      for (i = 0; i < nitem; i++) {
        for (j = 0; j < ndim; j++)
          samp_w(count, i, j) = oldw(i, j);
      }
      for (k = 0; k < nsample; k++) {
        for (j = 0; j < ndim; j++)
          samp_z(count, k, j) = oldz(k, j);
      }

      samp_gamma(count) = oldgamma;
      samp_sd_theta(count) = pr_sd_theta;
      if (spike_slab) {
        samp_pi(count) = pi_gamma;
        samp_xi(count) = xi;
      }

      if (missing_mar && nmissing > 0) {
        for (int mi = 0; mi < nmissing; mi++)
          samp_impute(count, mi) = impute_col(mi);
      }

      mle = 0.0;
      for (i = 0; i < nitem; i++) {
        for (c = 0; c < Kminus1; c++) {
          mle += R::dnorm4(oldbeta(i, c), pr_mean_beta, pr_sd_beta, 1);
        }
      }
      for (k = 0; k < nsample; k++)
        mle += R::dnorm4(oldtheta(k), pr_mean_theta, pr_sd_theta, 1);
      if (spike_slab) {
        if (pi_gamma == 1)
          mle += log_prior_gamma_lognormal(oldgamma, pr_slab_mean, pr_slab_sd);
        else
          mle +=
              log_prior_gamma_lognormal(oldgamma, pr_spike_mean, pr_spike_sd);
      } else {
        mle += log_prior_gamma_lognormal(oldgamma, pr_slab_mean, pr_slab_sd);
      }
      for (i = 0; i < nitem; i++)
        for (j = 0; j < ndim; j++)
          mle += R::dnorm4(oldw(i, j), pr_mean_w, pr_sd_w, 1);
      for (k = 0; k < nsample; k++)
        for (j = 0; j < ndim; j++)
          mle += R::dnorm4(oldz(k, j), pr_mean_z, pr_sd_z, 1);

      for (k = 0; k < nsample; k++) {
        for (i = 0; i < nitem; i++) {
          y = yint(k, i);
          if (y < 0)
            continue;

          const double theta_part =
              use_alpha ? (oldalpha(i) * oldtheta(k)) : oldtheta(k);
          eta = theta_part - oldgamma * dist(k, i);
          mle += log_grm_prob_y0based(y, oldbeta.row(i), eta);
        }
      }

      sample_mle(count) = mle;
      count++;
    }

    if (verbose) {
      int percent = 0;
      if (iter % nprint == 0 || iter == (niter - 1)) {
        percent = (iter * 100) / niter;
        Rprintf("Iteration: %.5u %3d%% ", iter, percent);
        for (c = 0; c < std::min(3, Kminus1); c++) {
          Rprintf("% .3f ", oldbeta(0, c));
        }
        if (Kminus1 > 3)
          Rprintf("... ");
        if (use_alpha) {
          Rprintf(" a1=%.3f ", oldalpha(0));
        }
        Rprintf(" %.3f ", oldgamma);
        Rprintf(" %.3f\n", pr_sd_theta);
      }
    } else {
      progressbar(iter + 1, niter);
    }
  } // for iter

  Rcpp::List output;
  output["beta"] = samp_beta;
  output["theta"] = samp_theta;
  output["z"] = samp_z;
  output["w"] = samp_w;
  output["gamma"] = samp_gamma;
  if (use_alpha) {
    output["alpha"] = samp_alpha;
    output["accept_alpha"] = accept_alpha;
  }
  output["sigma_theta"] = samp_sd_theta;
  output["map"] = sample_mle;
  output["accept_beta"] = accept_beta;
  output["accept_theta"] = accept_theta;
  output["accept_z"] = accept_z;
  output["accept_w"] = accept_w;
  output["accept_gamma"] = accept_gamma;

  Rcpp::List tuning;
  tuning["use_adapt"] = do_adapt;
  tuning["adapt_interval"] = adapt_ctrl.adapt_interval;
  tuning["adapt_rate"] = adapt_ctrl.adapt_rate;
  tuning["jump_min"] = adapt_ctrl.jump_min;
  tuning["jump_max"] = adapt_ctrl.jump_max;
  tuning["target_accept_beta"] = adapt_ctrl.target_accept_beta;
  tuning["target_accept_theta"] = adapt_ctrl.target_accept_theta;
  tuning["target_accept_gamma"] = adapt_ctrl.target_accept_gamma;
  tuning["target_accept_alpha"] = adapt_ctrl.target_accept_alpha;
  tuning["target_accept_zw"] = adapt_ctrl.target_accept_zw;
  tuning["jump_beta_init"] = jump_beta;
  tuning["jump_theta_init"] = jump_theta;
  tuning["jump_gamma_init"] = jump_gamma;
  tuning["jump_alpha_init"] = jump_alpha;
  tuning["jump_z_init"] = jump_z;
  tuning["jump_w_init"] = jump_w;
  tuning["jump_beta_final"] = cur_jump_beta;
  tuning["jump_theta_final"] = cur_jump_theta;
  tuning["jump_gamma_final"] = cur_jump_gamma;
  tuning["jump_alpha_final"] = cur_jump_alpha;
  tuning["jump_z_final"] = cur_jump_z;
  tuning["jump_w_final"] = cur_jump_w;
  tuning["accept_beta_burn"] =
      (burn_prop_beta > 0.0) ? (burn_acc_beta / burn_prop_beta) : NA_REAL;
  tuning["accept_theta_burn"] =
      (burn_prop_theta > 0.0) ? (burn_acc_theta / burn_prop_theta) : NA_REAL;
  tuning["accept_gamma_burn"] = (!fixed_gamma && burn_prop_gamma > 0.0)
                                    ? (burn_acc_gamma / burn_prop_gamma)
                                    : NA_REAL;
  tuning["accept_alpha_burn"] = (use_alpha && burn_prop_alpha > 0.0)
                                    ? (burn_acc_alpha / burn_prop_alpha)
                                    : NA_REAL;
  tuning["accept_z_burn"] =
      (burn_prop_z > 0.0) ? (burn_acc_z / burn_prop_z) : NA_REAL;
  tuning["accept_w_burn"] =
      (burn_prop_w > 0.0) ? (burn_acc_w / burn_prop_w) : NA_REAL;

  // Acceptance in the most recent *completed* adaptation window.
  tuning["accept_beta_lastwin"] = lastwin_accept_beta;
  tuning["accept_theta_lastwin"] = lastwin_accept_theta;
  tuning["accept_gamma_lastwin"] = lastwin_accept_gamma;
  tuning["accept_alpha_lastwin"] = lastwin_accept_alpha;
  tuning["accept_z_lastwin"] = lastwin_accept_z;
  tuning["accept_w_lastwin"] = lastwin_accept_w;

  tuning["n_adapt_steps"] = adapt_step;
  output["tuning"] = tuning;
  if (spike_slab) {
    output["pi"] = samp_pi;
    output["xi"] = samp_xi;
  }
  if (missing_mar) {
    output["impute"] = samp_impute;
  }

  return output;
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta, const double pr_a_theta,
    const double pr_b_theta, const double pr_mean_theta, double pr_sd_theta,
    const double pr_mean_gamma, const double pr_sd_gamma,
    const double jump_gamma, const double missing, const bool verbose,
    const bool fix_theta_sd, Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, 0.0, 1.0, pr_mean_gamma, pr_sd_gamma, 1.0,
      1.0, jump_gamma, missing, verbose, fix_theta_sd, false, false, false, 1.0,
      0.5, 1.0, true, false, adapt);
} // function end

// [[Rcpp::export]]
Rcpp::List lsirmgrm_fixed_gamma_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta, const double pr_a_theta,
    const double pr_b_theta, const double pr_mean_theta, double pr_sd_theta,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(data, ndim, ncat, niter, nburn, nthin, nprint,
                           jump_beta, jump_theta, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_a_theta, pr_b_theta, pr_mean_theta,
                           pr_sd_theta, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.025,
                           missing, verbose, fix_theta_sd, true, false, false,
                           1.0, 0.5, 1.0, true, false, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm_ss_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_beta_a, const double pr_beta_b,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, pr_spike_mean, pr_spike_sd, pr_slab_mean,
      pr_slab_sd, pr_beta_a, pr_beta_b, jump_gamma, missing, verbose,
      fix_theta_sd, false, true, false, 1.0, 0.5, 1.0, true, false, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm_mar_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta, const double pr_a_theta,
    const double pr_b_theta, const double pr_mean_theta, double pr_sd_theta,
    const double pr_mean_gamma, const double pr_sd_gamma,
    const double jump_gamma, const double missing, const bool verbose,
    const bool fix_theta_sd, Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, 0.0, 1.0, pr_mean_gamma, pr_sd_gamma, 1.0,
      1.0, jump_gamma, missing, verbose, fix_theta_sd, false, false, false, 1.0,
      0.5, 1.0, true, true, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm_fixed_gamma_mar_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta, const double pr_a_theta,
    const double pr_b_theta, const double pr_mean_theta, double pr_sd_theta,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(data, ndim, ncat, niter, nburn, nthin, nprint,
                           jump_beta, jump_theta, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_a_theta, pr_b_theta, pr_mean_theta,
                           pr_sd_theta, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.025,
                           missing, verbose, fix_theta_sd, true, false, false,
                           1.0, 0.5, 1.0, true, true, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm_ss_mar_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_beta_a, const double pr_beta_b,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, pr_spike_mean, pr_spike_sd, pr_slab_mean,
      pr_slab_sd, pr_beta_a, pr_beta_b, jump_gamma, missing, verbose,
      fix_theta_sd, false, true, false, 1.0, 0.5, 1.0, true, true, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm2pl_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_alpha, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_mean_gamma,
    const double pr_sd_gamma, const double jump_gamma,
    const double pr_mean_alpha, const double pr_sd_alpha, const double missing,
    const bool verbose, const bool fix_theta_sd, const bool fix_alpha_1,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, 0.0, 1.0, pr_mean_gamma, pr_sd_gamma, 1.0,
      1.0, jump_gamma, missing, verbose, fix_theta_sd, false, false, true,
      jump_alpha, pr_mean_alpha, pr_sd_alpha, fix_alpha_1, false, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm2pl_fixed_gamma_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_alpha, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_mean_alpha,
    const double pr_sd_alpha, const double missing, const bool verbose,
    const bool fix_theta_sd, const bool fix_alpha_1,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.025, missing,
      verbose, fix_theta_sd, true, false, true, jump_alpha, pr_mean_alpha,
      pr_sd_alpha, fix_alpha_1, false, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm2pl_ss_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_alpha, const double jump_gamma,
    const double jump_z, const double jump_w, const double pr_mean_beta,
    const double pr_sd_beta, const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_beta_a, const double pr_beta_b,
    const double pr_mean_alpha, const double pr_sd_alpha, const double missing,
    const bool verbose, const bool fix_theta_sd, const bool fix_alpha_1,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, pr_spike_mean, pr_spike_sd, pr_slab_mean,
      pr_slab_sd, pr_beta_a, pr_beta_b, jump_gamma, missing, verbose,
      fix_theta_sd, false, true, true, jump_alpha, pr_mean_alpha, pr_sd_alpha,
      fix_alpha_1, false, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm2pl_mar_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_alpha, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_mean_gamma,
    const double pr_sd_gamma, const double jump_gamma,
    const double pr_mean_alpha, const double pr_sd_alpha, const double missing,
    const bool verbose, const bool fix_theta_sd, const bool fix_alpha_1,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, 0.0, 1.0, pr_mean_gamma, pr_sd_gamma, 1.0,
      1.0, jump_gamma, missing, verbose, fix_theta_sd, false, false, true,
      jump_alpha, pr_mean_alpha, pr_sd_alpha, fix_alpha_1, true, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm2pl_fixed_gamma_mar_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_alpha, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_mean_alpha,
    const double pr_sd_alpha, const double missing, const bool verbose,
    const bool fix_theta_sd, const bool fix_alpha_1,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, 0.0, 1.0, 0.5, 1.0, 1.0, 1.0, 0.025, missing,
      verbose, fix_theta_sd, true, false, true, jump_alpha, pr_mean_alpha,
      pr_sd_alpha, fix_alpha_1, true, adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirmgrm2pl_ss_mar_cpp(
    arma::mat data, const int ndim, const int ncat, const int niter,
    const int nburn, const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_alpha, const double jump_gamma,
    const double jump_z, const double jump_w, const double pr_mean_beta,
    const double pr_sd_beta, const double pr_a_theta, const double pr_b_theta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_beta_a, const double pr_beta_b,
    const double pr_mean_alpha, const double pr_sd_alpha, const double missing,
    const bool verbose, const bool fix_theta_sd, const bool fix_alpha_1,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirmgrm_internal(
      data, ndim, ncat, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_a_theta, pr_b_theta,
      pr_mean_theta, pr_sd_theta, pr_spike_mean, pr_spike_sd, pr_slab_mean,
      pr_slab_sd, pr_beta_a, pr_beta_b, jump_gamma, missing, verbose,
      fix_theta_sd, false, true, true, jump_alpha, pr_mean_alpha, pr_sd_alpha,
      fix_alpha_1, true, adapt);
}
