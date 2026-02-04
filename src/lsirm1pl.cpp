#include "progress.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;

namespace {

inline double clamp_double(const double x, const double lo, const double hi) {
  if (x < lo)
    return lo;
  if (x > hi)
    return hi;
  return x;
}

inline bool list_has(const Rcpp::List &x, const char *name) {
  return x.containsElementNamed(name);
}

inline bool list_get_bool(const Rcpp::List &x, const char *name,
                          const bool fallback) {
  if (!list_has(x, name))
    return fallback;
  return Rcpp::as<bool>(x[name]);
}

inline int list_get_int(const Rcpp::List &x, const char *name,
                        const int fallback) {
  if (!list_has(x, name))
    return fallback;
  return Rcpp::as<int>(x[name]);
}

inline double list_get_double(const Rcpp::List &x, const char *name,
                              const double fallback) {
  if (!list_has(x, name))
    return fallback;
  return Rcpp::as<double>(x[name]);
}

struct AdaptControlLsirm {
  bool use_adapt = false;
  int adapt_interval = 100;
  double adapt_rate = 1.0;
  double decay_rate = 0.5;
  double jump_min = 1e-6;
  double jump_max = 10.0;
  double target_accept_beta = 0.44;
  double target_accept_theta = 0.44;
  double target_accept_gamma = 0.44;
  double target_accept_zw = 0.234;
  double target_accept_alpha = 0.44;
};

inline AdaptControlLsirm
parse_adapt_lsirm(const Rcpp::Nullable<Rcpp::List> &adapt, const int ndim) {
  AdaptControlLsirm out;
  if (adapt.isNull())
    return out;

  const Rcpp::List a(adapt);
  out.use_adapt = list_get_bool(a, "use_adapt", false);
  out.adapt_interval =
      std::max(1, list_get_int(a, "adapt_interval", out.adapt_interval));
  out.adapt_rate = list_get_double(a, "adapt_rate", out.adapt_rate);
  out.decay_rate = list_get_double(a, "decay_rate", out.decay_rate);
  out.jump_min = list_get_double(a, "jump_min", out.jump_min);
  out.jump_max = list_get_double(a, "jump_max", out.jump_max);

  const double target_scalar = list_get_double(a, "target_accept", 0.44);
  out.target_accept_beta =
      list_get_double(a, "target_accept_beta", target_scalar);
  out.target_accept_theta =
      list_get_double(a, "target_accept_theta", target_scalar);
  out.target_accept_gamma =
      list_get_double(a, "target_accept_gamma", target_scalar);
  out.target_accept_alpha =
      list_get_double(a, "target_accept_alpha", target_scalar);
  out.target_accept_zw =
      list_get_double(a, "target_accept_zw", out.target_accept_zw);

  return out;
}

inline void adapt_log_jump(double &log_jump, const double acc_rate,
                           const double target, const AdaptControlLsirm &ctrl,
                           const int adapt_step) {
  const double step_size =
      ctrl.adapt_rate /
      std::pow(static_cast<double>(std::max(1, adapt_step)), ctrl.decay_rate);
  log_jump += step_size * (acc_rate - target);
  const double new_jump =
      clamp_double(std::exp(log_jump), ctrl.jump_min, ctrl.jump_max);
  log_jump = std::log(new_jump);
}

inline double safe_rate(const double acc, const double prop) {
  if (prop <= 0.0)
    return NA_REAL;
  return acc / prop;
}

} // namespace

static Rcpp::List lsirm1pl_internal(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_mean_gamma,
    const double pr_sd_gamma, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_a_theta, const double pr_b_theta,
    const double pr_beta_a, const double pr_beta_b, const double pr_a_eps,
    const double pr_b_eps, const double missing, const bool verbose,
    const bool fix_theta_sd, const bool fixed_gamma, const bool mar,
    const bool spike_slab, const bool normal_model,
    Rcpp::Nullable<Rcpp::List> adapt) {

  const AdaptControlLsirm adapt_ctrl = parse_adapt_lsirm(adapt, ndim);
  const bool do_adapt = adapt_ctrl.use_adapt && (nburn > 0);

  double cur_jump_beta = jump_beta;
  double cur_jump_theta = jump_theta;
  double cur_jump_z = jump_z;
  double cur_jump_w = jump_w;
  double cur_jump_gamma = jump_gamma;

  double log_jump_beta = std::log(cur_jump_beta);
  double log_jump_theta = std::log(cur_jump_theta);
  double log_jump_z = std::log(cur_jump_z);
  double log_jump_w = std::log(cur_jump_w);
  double log_jump_gamma = std::log(cur_jump_gamma);

  int adapt_step = 0;

  double win_prop_beta = 0.0, win_acc_beta = 0.0;
  double win_prop_theta = 0.0, win_acc_theta = 0.0;
  double win_prop_z = 0.0, win_acc_z = 0.0;
  double win_prop_w = 0.0, win_acc_w = 0.0;
  double win_prop_gamma = 0.0, win_acc_gamma = 0.0;

  double burn_prop_beta = 0.0, burn_acc_beta = 0.0;
  double burn_prop_theta = 0.0, burn_acc_theta = 0.0;
  double burn_prop_z = 0.0, burn_acc_z = 0.0;
  double burn_prop_w = 0.0, burn_acc_w = 0.0;
  double burn_prop_gamma = 0.0, burn_acc_gamma = 0.0;

  double lastwin_accept_beta = NA_REAL;
  double lastwin_accept_theta = NA_REAL;
  double lastwin_accept_z = NA_REAL;
  double lastwin_accept_w = NA_REAL;
  double lastwin_accept_gamma = NA_REAL;

  int i, j, k, count, accept;
  double num, den, old_like_beta, new_like_beta, old_like_theta, new_like_theta;
  double old_like_z, new_like_z, old_like_w, new_like_w, old_like_gamma,
      new_like_gamma;
  double ratio, un, dist_temp;
  double pr_mean_z = 0.0, pr_sd_z = 1.0, pr_mean_w = 0.0, pr_sd_w = 1.0, mle;

  const int nsample = data.n_rows;
  const int nitem = data.n_cols;

  arma::dvec oldbeta(nitem, fill::randu);
  oldbeta = oldbeta * 4.0 - 2.0;
  arma::dvec newbeta = oldbeta;

  arma::dvec oldtheta(nsample, fill::randu);
  oldtheta = oldtheta * 4.0 - 2.0;
  arma::dvec newtheta = oldtheta;

  arma::dmat oldz(nsample, ndim, fill::randu);
  oldz = oldz * 2.0 - 1.0;
  arma::dmat newz = oldz;

  arma::dmat oldw(nitem, ndim, fill::randu);
  oldw = oldw * 2.0 - 1.0;
  arma::dmat neww = oldw;

  double oldgamma = 1.0, newgamma = 1.0;
  int pi_gamma = 1;
  double xi = 0.5;

  if (spike_slab) {
    pi_gamma = (int)R::rbinom(1.0, 0.5);
    xi = R::rbeta(pr_beta_a + (double)pi_gamma,
                  pr_beta_b + 1.0 - (double)pi_gamma);
  }

  // For Normal model
  double pr_sd = 1.0;
  double post_a_sigma, post_b_sigma;

  // MAR
  arma::dvec impute_col;
  int nmissing = 0;
  arma::ivec miss_row, miss_col;
  if (mar) {
    for (int rr = 0; rr < nsample; rr++) {
      for (int cc = 0; cc < nitem; cc++) {
        if (data(rr, cc) == missing)
          nmissing++;
      }
    }
    miss_row = arma::ivec(nmissing);
    miss_col = arma::ivec(nmissing);
    impute_col = arma::dvec(nmissing);
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

  const int nmcmc = (niter - nburn) / nthin;
  arma::dmat samp_beta(nmcmc, nitem, fill::zeros);
  arma::dmat samp_theta(nmcmc, nsample, fill::zeros);
  arma::dcube samp_z(nmcmc, nsample, ndim, fill::zeros);
  arma::dcube samp_w(nmcmc, nitem, ndim, fill::zeros);
  arma::dvec samp_sd_theta(nmcmc, fill::zeros);
  arma::dvec samp_mle(nmcmc, fill::zeros);
  arma::dvec samp_gamma(nmcmc, fill::zeros);
  arma::dvec samp_sd(nmcmc, fill::zeros); // Only used if normal_model
  arma::dmat samp_impute;
  if (mar)
    samp_impute = arma::dmat(nmcmc, nmissing, fill::zeros);
  arma::dvec samp_pi;
  arma::dvec samp_xi;
  if (spike_slab) {
    samp_pi = arma::dvec(nmcmc, fill::zeros);
    samp_xi = arma::dvec(nmcmc, fill::zeros);
  }

  arma::dvec accept_beta(nitem, fill::zeros);
  arma::dvec accept_theta(nsample, fill::zeros);
  arma::dvec accept_z(nsample, fill::zeros);
  arma::dvec accept_w(nitem, fill::zeros);
  double accept_gamma = 0;

  accept = count = 0;

  arma::dmat dist(nsample, nitem, fill::zeros);
  arma::dvec old_dist_k(nitem, fill::zeros);
  arma::dvec new_dist_k(nitem, fill::zeros);
  arma::dvec old_dist_i(nsample, fill::zeros);
  arma::dvec new_dist_i(nsample, fill::zeros);

  for (int iter = 0; iter < niter; iter++) {
    if (iter % 10 == 0)
      Rcpp::checkUserInterrupt();

    // Distance calculation
    dist.fill(0.0);
    for (i = 0; i < nitem; i++) {
      for (k = 0; k < nsample; k++) {
        dist_temp = 0.0;
        for (j = 0; j < ndim; j++)
          dist_temp += std::pow((oldz(k, j) - oldw(i, j)), 2.0);
        dist(k, i) = std::sqrt(dist_temp);
      }
    }

    // Imputation
    if (mar) {
      for (int mi = 0; mi < nmissing; mi++) {
        int rr = miss_row(mi);
        int cc = miss_col(mi);
        double eta = oldbeta(cc) + oldtheta(rr) - oldgamma * dist(rr, cc);
        if (normal_model) {
          data(rr, cc) = R::rnorm(eta, pr_sd);
        } else {
          double prob = 1.0 / (1.0 + std::exp(-eta));
          data(rr, cc) = (R::runif(0, 1) < prob) ? 1.0 : 0.0;
        }
        impute_col(mi) = data(rr, cc);
      }
    }

    // Beta Update
    for (i = 0; i < nitem; i++) {
      newbeta(i) = R::rnorm(oldbeta(i), cur_jump_beta);
      if (do_adapt && iter < nburn) {
        win_prop_beta += 1.0;
        burn_prop_beta += 1.0;
      }

      old_like_beta = new_like_beta = 0.0;
      for (k = 0; k < nsample; k++) {
        if (data(k, i) != missing) {
          if (normal_model) {
            new_like_beta += -std::pow(data(k, i) - newbeta(i) - oldtheta(k) +
                                           oldgamma * dist(k, i),
                                       2.0) /
                             (2.0 * std::pow(pr_sd, 2.0));
            old_like_beta += -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                           oldgamma * dist(k, i),
                                       2.0) /
                             (2.0 * std::pow(pr_sd, 2.0));
          } else {
            double term_new = newbeta(i) + oldtheta(k) - oldgamma * dist(k, i);
            double term_old = oldbeta(i) + oldtheta(k) - oldgamma * dist(k, i);
            if (data(k, i) == 1.0) {
              new_like_beta += -std::log(1.0 + std::exp(-term_new));
              old_like_beta += -std::log(1.0 + std::exp(-term_old));
            } else {
              new_like_beta += -std::log(1.0 + std::exp(term_new));
              old_like_beta += -std::log(1.0 + std::exp(term_old));
            }
          }
        }
      }
      num = new_like_beta + R::dnorm4(newbeta(i), pr_mean_beta, pr_sd_beta, 1);
      den = old_like_beta + R::dnorm4(oldbeta(i), pr_mean_beta, pr_sd_beta, 1);
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
        oldbeta(i) = newbeta(i);
        if (!do_adapt || iter >= nburn) {
          accept_beta(i) += 1.0 / (niter - nburn);
        }
        if (do_adapt && iter < nburn) {
          win_acc_beta += 1.0;
          burn_acc_beta += 1.0;
        }
      }
    }

    // Theta Update
    for (k = 0; k < nsample; k++) {
      newtheta(k) = R::rnorm(oldtheta(k), cur_jump_theta);
      if (do_adapt && iter < nburn) {
        win_prop_theta += 1.0;
        burn_prop_theta += 1.0;
      }

      old_like_theta = new_like_theta = 0.0;
      for (i = 0; i < nitem; i++) {
        if (data(k, i) != missing) {
          if (normal_model) {
            new_like_theta += -std::pow(data(k, i) - oldbeta(i) - newtheta(k) +
                                            oldgamma * dist(k, i),
                                        2.0) /
                              (2.0 * std::pow(pr_sd, 2.0));
            old_like_theta += -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                            oldgamma * dist(k, i),
                                        2.0) /
                              (2.0 * std::pow(pr_sd, 2.0));
          } else {
            double term_new = oldbeta(i) + newtheta(k) - oldgamma * dist(k, i);
            double term_old = oldbeta(i) + oldtheta(k) - oldgamma * dist(k, i);
            if (data(k, i) == 1.0) {
              new_like_theta += -std::log(1.0 + std::exp(-term_new));
              old_like_theta += -std::log(1.0 + std::exp(-term_old));
            } else {
              new_like_theta += -std::log(1.0 + std::exp(term_new));
              old_like_theta += -std::log(1.0 + std::exp(term_old));
            }
          }
        }
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
          accept_theta(k) += 1.0 / (niter - nburn);
        }
        if (do_adapt && iter < nburn) {
          win_acc_theta += 1.0;
          burn_acc_theta += 1.0;
        }
      }
    }

    // Gamma Update
    if (!fixed_gamma) {
      newgamma = R::rlnorm(std::log(oldgamma), cur_jump_gamma);
      if (do_adapt && iter < nburn) {
        win_prop_gamma += 1.0;
        burn_prop_gamma += 1.0;
      }

      old_like_gamma = new_like_gamma = 0.0;
      for (k = 0; k < nsample; k++) {
        for (i = 0; i < nitem; i++) {
          if (data(k, i) != missing) {
            if (normal_model) {
              new_like_gamma +=
                  -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                newgamma * dist(k, i),
                            2.0) /
                  (2.0 * std::pow(pr_sd, 2.0));
              old_like_gamma +=
                  -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                oldgamma * dist(k, i),
                            2.0) /
                  (2.0 * std::pow(pr_sd, 2.0));
            } else {
              double term_new =
                  oldbeta(i) + oldtheta(k) - newgamma * dist(k, i);
              double term_old =
                  oldbeta(i) + oldtheta(k) - oldgamma * dist(k, i);
              if (data(k, i) == 1.0) {
                new_like_gamma += -std::log(1.0 + std::exp(-term_new));
                old_like_gamma += -std::log(1.0 + std::exp(-term_old));
              } else {
                new_like_gamma += -std::log(1.0 + std::exp(term_new));
                old_like_gamma += -std::log(1.0 + std::exp(term_old));
              }
            }
          }
        }
      }

      num = new_like_gamma +
            R::dlnorm(oldgamma, std::log(newgamma), cur_jump_gamma, 1);
      den = old_like_gamma +
            R::dlnorm(newgamma, std::log(oldgamma), cur_jump_gamma, 1);

      if (spike_slab) {
        if (pi_gamma == 1) {
          num += R::dlnorm(newgamma, pr_slab_mean, pr_slab_sd, 1);
          den += R::dlnorm(oldgamma, pr_slab_mean, pr_slab_sd, 1);
        } else {
          num += R::dlnorm(newgamma, pr_spike_mean, pr_spike_sd, 1);
          den += R::dlnorm(oldgamma, pr_spike_mean, pr_spike_sd, 1);
        }
      } else {
        num += R::dlnorm(newgamma, pr_mean_gamma, pr_sd_gamma, 1);
        den += R::dlnorm(oldgamma, pr_mean_gamma, pr_sd_gamma, 1);
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
        accept_gamma += 1.0 / niter;
        if (do_adapt && iter < nburn) {
          win_acc_gamma += 1.0;
          burn_acc_gamma += 1.0;
        }
      } else {
        newgamma = oldgamma;
      }

      if (spike_slab) {
        if (spike_slab) {
          double log_d_slab = R::dlnorm(oldgamma, pr_slab_mean, pr_slab_sd, 1);
          double log_d_spike =
              R::dlnorm(oldgamma, pr_spike_mean, pr_spike_sd, 1);

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

          // log_den_pi = log(exp(log_term_spike) + exp(log_num_pi))
          double max_val = std::max(log_term_spike, log_num_pi);
          log_den_pi = max_val + std::log(std::exp(log_term_spike - max_val) +
                                          std::exp(log_num_pi - max_val));

          double prob_pi = std::exp(log_num_pi - log_den_pi);
          if (std::isnan(prob_pi))
            prob_pi = 0.0; // Fail-safe

          pi_gamma = (int)R::rbinom(1.0, prob_pi);
          xi = R::rbeta(pr_beta_a + (double)pi_gamma,
                        pr_beta_b + 1.0 - (double)pi_gamma);
        }
      }
    }

    // Z Update
    for (k = 0; k < nsample; k++) {
      for (j = 0; j < ndim; j++)
        newz(k, j) = R::rnorm(oldz(k, j), cur_jump_z);
      if (do_adapt && iter < nburn) {
        win_prop_z += 1.0;
        burn_prop_z += 1.0;
      }

      old_like_z = new_like_z = 0.0;

      for (i = 0; i < nitem; i++) {
        double dist_new_sq = 0.0;
        double dist_old_sq = 0.0;
        for (j = 0; j < ndim; j++) {
          dist_new_sq += std::pow(newz(k, j) - oldw(i, j), 2.0);
          dist_old_sq += std::pow(oldz(k, j) - oldw(i, j), 2.0);
        }
        new_dist_k(i) = std::sqrt(dist_new_sq);
        old_dist_k(i) = std::sqrt(dist_old_sq);

        if (data(k, i) != missing) {
          if (normal_model) {
            new_like_z += -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                        oldgamma * new_dist_k(i),
                                    2.0) /
                          (2.0 * std::pow(pr_sd, 2.0));
            old_like_z += -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                        oldgamma * old_dist_k(i),
                                    2.0) /
                          (2.0 * std::pow(pr_sd, 2.0));
          } else {
            double term_new =
                oldbeta(i) + oldtheta(k) - oldgamma * new_dist_k(i);
            double term_old =
                oldbeta(i) + oldtheta(k) - oldgamma * old_dist_k(i);
            if (data(k, i) == 1.0) {
              new_like_z += -std::log(1.0 + std::exp(-term_new));
              old_like_z += -std::log(1.0 + std::exp(-term_old));
            } else {
              new_like_z += -std::log(1.0 + std::exp(term_new));
              old_like_z += -std::log(1.0 + std::exp(term_old));
            }
          }
        }
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
          accept_z(k) += 1.0 / (niter - nburn);
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

    // W Update
    for (i = 0; i < nitem; i++) {
      for (j = 0; j < ndim; j++)
        neww(i, j) = R::rnorm(oldw(i, j), cur_jump_w);
      if (do_adapt && iter < nburn) {
        win_prop_w += 1.0;
        burn_prop_w += 1.0;
      }

      old_like_w = new_like_w = 0.0;
      for (k = 0; k < nsample; k++) {
        double dist_new_sq = 0.0;
        double dist_old_sq = 0.0;
        for (j = 0; j < ndim; j++) {
          dist_new_sq += std::pow(oldz(k, j) - neww(i, j), 2.0);
          dist_old_sq += std::pow(oldz(k, j) - oldw(i, j), 2.0);
        }
        new_dist_i(k) = std::sqrt(dist_new_sq);
        old_dist_i(k) = std::sqrt(dist_old_sq);

        if (data(k, i) != missing) {
          if (normal_model) {
            new_like_w += -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                        oldgamma * new_dist_i(k),
                                    2.0) /
                          (2.0 * std::pow(pr_sd, 2.0));
            old_like_w += -std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                        oldgamma * old_dist_i(k),
                                    2.0) /
                          (2.0 * std::pow(pr_sd, 2.0));
          } else {
            double term_new =
                oldbeta(i) + oldtheta(k) - oldgamma * new_dist_i(k);
            double term_old =
                oldbeta(i) + oldtheta(k) - oldgamma * old_dist_i(k);
            if (data(k, i) == 1.0) {
              new_like_w += -std::log(1.0 + std::exp(-term_new));
              old_like_w += -std::log(1.0 + std::exp(-term_old));
            } else {
              new_like_w += -std::log(1.0 + std::exp(term_new));
              old_like_w += -std::log(1.0 + std::exp(term_old));
            }
          }
        }
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
          accept_w(i) += 1.0 / (niter - nburn);
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

    // Sigma Theta Update
    if (!fix_theta_sd) {
      double post_a = 2.0 * pr_a_theta + nsample;
      double post_b = pr_b_theta;
      for (j = 0; j < nsample; j++)
        post_b += std::pow((oldtheta(j) - pr_mean_theta), 2.0) / 2.0;
      pr_sd_theta = std::sqrt(2.0 * post_b * (1.0 / R::rchisq(post_a)));
    }

    // Sigma Update (Normal)
    if (normal_model) {
      post_a_sigma = 2.0 * pr_a_eps + nsample * nitem;
      post_b_sigma = pr_b_eps;
      for (k = 0; k < nsample; k++) {
        for (i = 0; i < nitem; i++) {
          post_b_sigma += std::pow(data(k, i) - oldbeta(i) - oldtheta(k) +
                                       oldgamma * dist(k, i),
                                   2.0) /
                          2.0;
        }
      }
      pr_sd = std::sqrt(2.0 * post_b_sigma * (1.0 / R::rchisq(post_a_sigma)));
    }

    // Adapt
    if (do_adapt && (iter + 1) % adapt_ctrl.adapt_interval == 0) {
      adapt_step++;

      double rate;
      rate = safe_rate(win_acc_beta, win_prop_beta);
      if (!Rcpp::NumericVector::is_na(rate)) {
        adapt_log_jump(log_jump_beta, rate, adapt_ctrl.target_accept_beta,
                       adapt_ctrl, adapt_step);
        cur_jump_beta = std::exp(log_jump_beta);
        lastwin_accept_beta = rate;
      }

      rate = safe_rate(win_acc_theta, win_prop_theta);
      if (!Rcpp::NumericVector::is_na(rate)) {
        adapt_log_jump(log_jump_theta, rate, adapt_ctrl.target_accept_theta,
                       adapt_ctrl, adapt_step);
        cur_jump_theta = std::exp(log_jump_theta);
        lastwin_accept_theta = rate;
      }

      if (!fixed_gamma) {
        rate = safe_rate(win_acc_gamma, win_prop_gamma);
        if (!Rcpp::NumericVector::is_na(rate)) {
          adapt_log_jump(log_jump_gamma, rate, adapt_ctrl.target_accept_gamma,
                         adapt_ctrl, adapt_step);
          cur_jump_gamma = std::exp(log_jump_gamma);
          lastwin_accept_gamma = rate;
        }
      }

      rate = safe_rate(win_acc_z, win_prop_z);
      if (!Rcpp::NumericVector::is_na(rate)) {
        adapt_log_jump(log_jump_z, rate, adapt_ctrl.target_accept_zw,
                       adapt_ctrl, adapt_step);
        cur_jump_z = std::exp(log_jump_z);
        lastwin_accept_z = rate;
      }

      rate = safe_rate(win_acc_w, win_prop_w);
      if (!Rcpp::NumericVector::is_na(rate)) {
        adapt_log_jump(log_jump_w, rate, adapt_ctrl.target_accept_zw,
                       adapt_ctrl, adapt_step);
        cur_jump_w = std::exp(log_jump_w);
        lastwin_accept_w = rate;
      }

      win_prop_beta = win_acc_beta = 0.0;
      win_prop_theta = win_acc_theta = 0.0;
      win_prop_z = win_acc_z = 0.0;
      win_prop_w = win_acc_w = 0.0;
      win_prop_gamma = win_acc_gamma = 0.0;
    }

    // Save
    if (iter >= nburn && iter % nthin == 0) {
      for (i = 0; i < nitem; i++)
        samp_beta(count, i) = oldbeta(i);
      for (k = 0; k < nsample; k++)
        samp_theta(count, k) = oldtheta(k);
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
      if (normal_model)
        samp_sd(count) = pr_sd;
      if (spike_slab) {
        samp_pi(count) = pi_gamma; // Store indicator
        samp_xi(count) = xi;       // Store xi (probability)
      }
      if (mar) {
        for (int mi = 0; mi < nmissing; mi++)
          samp_impute(count, mi) = impute_col(mi);
      }

      // MLE
      mle = 0.0;
      for (i = 0; i < nitem; i++)
        mle += R::dnorm4(oldbeta(i), pr_mean_beta, pr_sd_beta, 1);
      for (k = 0; k < nsample; k++)
        mle += R::dnorm4(oldtheta(k), pr_mean_theta, pr_sd_theta, 1);
      if (!fixed_gamma && !spike_slab)
        mle += R::dlnorm(oldgamma, pr_mean_gamma, pr_sd_gamma, 1);
      // Spike slab prior not added to MLE for simplcity in original code,
      // following pattern

      for (i = 0; i < nitem; i++)
        for (j = 0; j < ndim; j++)
          mle += R::dnorm4(oldw(i, j), pr_mean_w, pr_sd_w, 1);
      for (k = 0; k < nsample; k++)
        for (j = 0; j < ndim; j++)
          mle += R::dnorm4(oldz(k, j), pr_mean_z, pr_sd_z, 1);

      for (k = 0; k < nsample; k++) {
        for (i = 0; i < nitem; i++) {
          if (data(k, i) != missing) {
            double term = oldbeta(i) + oldtheta(k) - oldgamma * dist(k, i);
            if (normal_model) {
              mle += R::dnorm4(data(k, i), term, pr_sd, 1);
            } else if (data(k, i) == 1.0) {
              mle += -std::log(1.0 + std::exp(-term));
            } else {
              mle += -std::log(1.0 + std::exp(term));
            }
          }
        }
      }
      samp_mle(count) = mle;
      count++;
    }

    if (verbose) {
      if (iter % nprint == 0) {
        Rprintf("Iteration: %.5u %3d%% ", iter, (iter * 100) / niter);
        Rprintf("Gamma: %.3f ", oldgamma);
        if (normal_model)
          Rprintf("Sigma: %.3f ", pr_sd);
        Rprintf("\n");
      }
    } else {
      progressbar(iter + 1, niter);
    }
  }

  Rcpp::List output;
  output["beta"] = samp_beta;
  output["theta"] = samp_theta;
  output["z"] = samp_z;
  output["w"] = samp_w;
  output["gamma"] = samp_gamma;
  output["sigma_theta"] = samp_sd_theta;
  if (normal_model)
    output["sigma"] = samp_sd;
  output["map"] = samp_mle;
  output["accept_beta"] = accept_beta;
  output["accept_theta"] = accept_theta;
  output["accept_z"] = accept_z;
  output["accept_w"] = accept_w;
  output["accept_gamma"] = fixed_gamma ? NA_REAL : accept_gamma;
  if (mar)
    output["impute"] = samp_impute;
  if (spike_slab) {
    output["pi"] = samp_pi;
    output["xi"] = samp_xi;
  }

  Rcpp::List tuning;
  tuning["use_adapt"] = adapt_ctrl.use_adapt;
  tuning["adapt_interval"] = adapt_ctrl.adapt_interval;
  tuning["adapt_rate"] = adapt_ctrl.adapt_rate;
  tuning["target_accept_beta"] = adapt_ctrl.target_accept_beta;
  tuning["target_accept_theta"] = adapt_ctrl.target_accept_theta;
  tuning["target_accept_zw"] = adapt_ctrl.target_accept_zw;
  tuning["target_accept_gamma"] = adapt_ctrl.target_accept_gamma;
  tuning["target_accept_alpha"] = NA_REAL; // 1PL doesn't use alpha

  tuning["jump_beta_init"] = jump_beta;
  tuning["jump_theta_init"] = jump_theta;
  tuning["jump_z_init"] = jump_z;
  tuning["jump_w_init"] = jump_w;
  tuning["jump_gamma_init"] = fixed_gamma ? NA_REAL : jump_gamma;
  tuning["jump_alpha_init"] = NA_REAL;

  tuning["jump_beta_final"] = cur_jump_beta;
  tuning["jump_theta_final"] = cur_jump_theta;
  tuning["jump_z_final"] = cur_jump_z;
  tuning["jump_w_final"] = cur_jump_w;
  tuning["jump_gamma_final"] = fixed_gamma ? NA_REAL : cur_jump_gamma;
  tuning["jump_alpha_final"] = NA_REAL;

  tuning["accept_beta_burn"] = safe_rate(burn_acc_beta, burn_prop_beta);
  tuning["accept_theta_burn"] = safe_rate(burn_acc_theta, burn_prop_theta);
  tuning["accept_z_burn"] = safe_rate(burn_acc_z, burn_prop_z);
  tuning["accept_w_burn"] = safe_rate(burn_acc_w, burn_prop_w);
  tuning["accept_gamma_burn"] =
      fixed_gamma ? NA_REAL : safe_rate(burn_acc_gamma, burn_prop_gamma);
  tuning["accept_alpha_burn"] = NA_REAL;

  tuning["accept_beta_lastwin"] = lastwin_accept_beta;

  tuning["accept_theta_lastwin"] = lastwin_accept_theta;
  tuning["accept_z_lastwin"] = lastwin_accept_z;
  tuning["accept_w_lastwin"] = lastwin_accept_w;
  tuning["accept_gamma_lastwin"] = fixed_gamma ? NA_REAL : lastwin_accept_gamma;
  tuning["accept_alpha_lastwin"] = NA_REAL;

  output["tuning"] = tuning;

  return output;
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_cpp(arma::mat data, const int ndim, const int niter,
                        const int nburn, const int nthin, const int nprint,
                        const double jump_beta, const double jump_theta,
                        const double jump_gamma, const double jump_z,
                        const double jump_w, const double pr_mean_beta,
                        const double pr_sd_beta, const double pr_mean_theta,
                        double pr_sd_theta, const double pr_mean_gamma,
                        const double pr_sd_gamma, const double pr_a_theta,
                        const double pr_b_theta, const bool verbose,
                        const bool fix_theta_sd,
                        Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, pr_mean_gamma, pr_sd_gamma, 0.0, 1.0, 0.0, 1.0, // Spike slab
      pr_a_theta, pr_b_theta, 1.0, 1.0, // Beta priors
      1.0, 1.0,                         // Eps priors
      99.0, verbose, fix_theta_sd, false, false, false,
      false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_mcar_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_mean_gamma,
    const double pr_sd_gamma, const double pr_a_theta, const double pr_b_theta,
    double pr_sd_theta, const double missing, const bool verbose,
    const bool fix_theta_sd, const bool fixed_gamma,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, pr_mean_gamma, pr_sd_gamma, 0.0, 1.0, 0.0, 1.0, // Spike slab
      pr_a_theta, pr_b_theta, 1.0, 1.0, // Beta priors
      1.0, 1.0,                         // Eps priors
      missing, verbose, fix_theta_sd, fixed_gamma, false, false,
      false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_mcar_ss_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_a_theta, const double pr_b_theta,
    const double pr_beta_a, const double pr_beta_b, double pr_sd_theta,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, 0.0, 1.0, // Gamma mean/sd unused
      pr_spike_mean, pr_spike_sd, pr_slab_mean, pr_slab_sd, pr_a_theta,
      pr_b_theta, pr_beta_a, pr_beta_b, 1.0, 1.0, missing, verbose,
      fix_theta_sd, false, false, true, false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_mar_cpp(arma::mat data, const int ndim, const int niter,
                            const int nburn, const int nthin, const int nprint,
                            const double jump_beta, const double jump_theta,
                            const double jump_gamma, const double jump_z,
                            const double jump_w, const double pr_mean_beta,
                            const double pr_sd_beta, const double pr_mean_theta,
                            const double pr_mean_gamma,
                            const double pr_sd_gamma, const double pr_a_theta,
                            const double pr_b_theta, double pr_sd_theta,
                            const double missing, const bool verbose,
                            const bool fix_theta_sd, const bool fixed_gamma,
                            Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, pr_mean_gamma, pr_sd_gamma, 0.0, 1.0, 0.0, 1.0, pr_a_theta,
      pr_b_theta, 1.0, 1.0, 1.0, 1.0, missing, verbose, fix_theta_sd,
      fixed_gamma, true, false, false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_mar_ss_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_a_theta, const double pr_b_theta,
    const double pr_beta_a, const double pr_beta_b, double pr_sd_theta,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, 0.0, 1.0, // Gamma mean/sd unused
      pr_spike_mean, pr_spike_sd, pr_slab_mean, pr_slab_sd, pr_a_theta,
      pr_b_theta, pr_beta_a, pr_beta_b, 1.0, 1.0, missing, verbose,
      fix_theta_sd, false, true, true, false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_mean_gamma,
    const double pr_sd_gamma, const double pr_a_theta, const double pr_b_theta,
    double pr_sd_theta, const double pr_a_eps, const double pr_b_eps,
    const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, pr_mean_gamma, pr_sd_gamma, 0.0, 1.0, 0.0, 1.0, pr_a_theta,
      pr_b_theta, 1.0, 1.0, pr_a_eps, pr_b_eps, 99.0, verbose, fix_theta_sd,
      false, false, false, true, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_ss_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_a_theta, const double pr_b_theta,
    double pr_sd_theta, const double pr_beta_a, const double pr_beta_b,
    const double pr_a_eps, const double pr_b_eps, const bool verbose,
    const bool fix_theta_sd, Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(data, ndim, niter, nburn, nthin, nprint, jump_beta,
                           jump_theta, jump_gamma, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0, 1.0,
                           pr_spike_mean, pr_spike_sd, pr_slab_mean, pr_slab_sd,
                           pr_a_theta, pr_b_theta, pr_beta_a, pr_beta_b,
                           pr_a_eps, pr_b_eps, 99.0, verbose, fix_theta_sd,
                           false, false, true, true, // mar, spike_slab, normal
                           adapt);
}
// [[Rcpp::export]]
Rcpp::List lsirm1pl_ss_cpp(arma::mat data, const int ndim, const int niter,
                           const int nburn, const int nthin, const int nprint,
                           const double jump_beta, const double jump_theta,
                           const double jump_gamma, const double jump_z,
                           const double jump_w, const double pr_mean_beta,
                           const double pr_sd_beta, const double pr_mean_theta,
                           const double pr_spike_mean, const double pr_spike_sd,
                           const double pr_slab_mean, const double pr_slab_sd,
                           const double pr_a_theta, const double pr_b_theta,
                           const double pr_beta_a, const double pr_beta_b,
                           double pr_sd_theta, const bool verbose,
                           const bool fix_theta_sd,
                           Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(data, ndim, niter, nburn, nthin, nprint, jump_beta,
                           jump_theta, jump_gamma, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0, 1.0,
                           pr_spike_mean, pr_spike_sd, pr_slab_mean, pr_slab_sd,
                           pr_a_theta, pr_b_theta, pr_beta_a, pr_beta_b, 1.0,
                           1.0, 99.0, verbose, fix_theta_sd, false, false, true,
                           false, // mar, spike_slab, normal
                           adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_mar_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_mean_gamma,
    const double pr_sd_gamma, const double pr_a_theta, const double pr_b_theta,
    double pr_sd_theta, const double pr_a_eps, const double pr_b_eps,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(data, ndim, niter, nburn, nthin, nprint, jump_beta,
                           jump_theta, jump_gamma, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_mean_theta, pr_sd_theta,
                           pr_mean_gamma, pr_sd_gamma, 0.0, 1.0, 0.0, 1.0,
                           pr_a_theta, pr_b_theta, 1.0, 1.0, pr_a_eps, pr_b_eps,
                           missing, verbose, fix_theta_sd, false, true, false,
                           true, // mar, spike_slab, normal
                           adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_mcar_ss_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_a_theta, const double pr_b_theta,
    double pr_sd_theta, const double pr_beta_a, const double pr_beta_b,
    const double pr_a_eps, const double pr_b_eps, const double missing,
    const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, 0.0, 1.0, pr_spike_mean, pr_spike_sd, pr_slab_mean,
      pr_slab_sd, pr_a_theta, pr_b_theta, pr_beta_a, pr_beta_b, pr_a_eps,
      pr_b_eps, missing, verbose, fix_theta_sd, false, false, true,
      true, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_mar_ss_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_spike_mean,
    const double pr_spike_sd, const double pr_slab_mean,
    const double pr_slab_sd, const double pr_a_theta, const double pr_b_theta,
    double pr_sd_theta, const double pr_beta_a, const double pr_beta_b,
    const double pr_a_eps, const double pr_b_eps, const double missing,
    const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      jump_gamma, jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta,
      pr_sd_theta, 0.0, 1.0, pr_spike_mean, pr_spike_sd, pr_slab_mean,
      pr_slab_sd, pr_a_theta, pr_b_theta, pr_beta_a, pr_beta_b, pr_a_eps,
      pr_b_eps, missing, verbose, fix_theta_sd, false, true, true,
      true, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_mcar_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_gamma, const double jump_z,
    const double jump_w, const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, const double pr_mean_gamma,
    const double pr_sd_gamma, const double pr_a_theta, const double pr_b_theta,
    double pr_sd_theta, const double pr_a_eps, const double pr_b_eps,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(data, ndim, niter, nburn, nthin, nprint, jump_beta,
                           jump_theta, jump_gamma, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_mean_theta, pr_sd_theta,
                           pr_mean_gamma, pr_sd_gamma, 0.0, 1.0, 0.0, 1.0,
                           pr_a_theta, pr_b_theta, 1.0, 1.0, pr_a_eps, pr_b_eps,
                           missing, verbose, fix_theta_sd, false, false, false,
                           true, // mar, spike_slab, normal
                           adapt);
}
// [[Rcpp::export]]
Rcpp::List lsirm1pl_fixed_gamma_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_a_theta,
    const double pr_b_theta, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta,
      0.0, // jump_gamma
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0,
      1.0,                                        // gamma priors
      0.0, 1.0, 0.0, 1.0,                         // spike slab
      pr_a_theta, pr_b_theta, 1.0, 1.0, 1.0, 1.0, // eps priors
      99.0, verbose, fix_theta_sd,
      true,                // fixed_gamma
      false, false, false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_fixed_gamma_mar_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_a_theta,
    const double pr_b_theta, const double missing, const bool verbose,
    const bool fix_theta_sd, Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta, 0.0,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0,
      1.0, 0.0, 1.0, 0.0, 1.0, pr_a_theta, pr_b_theta, 1.0, 1.0, 1.0, 1.0,
      missing, verbose, fix_theta_sd, true, true, false,
      false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_fixed_gamma_mcar_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_a_theta,
    const double pr_b_theta, const double missing, const bool verbose,
    const bool fix_theta_sd, Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta, 0.0,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0,
      1.0, 0.0, 1.0, 0.0, 1.0, pr_a_theta, pr_b_theta, 1.0, 1.0, 1.0, 1.0,
      missing, verbose, fix_theta_sd, true, false, false,
      false, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_fixed_gamma_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_a_theta,
    const double pr_b_theta, const double pr_a_eps, const double pr_b_eps,
    const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(
      data, ndim, niter, nburn, nthin, nprint, jump_beta, jump_theta, 0.0,
      jump_z, jump_w, pr_mean_beta, pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0,
      1.0, 0.0, 1.0, 0.0, 1.0, pr_a_theta, pr_b_theta, 1.0, 1.0, pr_a_eps,
      pr_b_eps, 99.0, verbose, fix_theta_sd, true, false, false,
      true, // mar, spike_slab, normal
      adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_fixed_gamma_mar_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_a_theta,
    const double pr_b_theta, const double pr_a_eps, const double pr_b_eps,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(data, ndim, niter, nburn, nthin, nprint, jump_beta,
                           jump_theta, 0.0, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0, 1.0,
                           0.0, 1.0, 0.0, 1.0, pr_a_theta, pr_b_theta, 1.0, 1.0,
                           pr_a_eps, pr_b_eps, missing, verbose, fix_theta_sd,
                           true, true, false, true, // mar, spike_slab, normal
                           adapt);
}

// [[Rcpp::export]]
Rcpp::List lsirm1pl_normal_fixed_gamma_mcar_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const int nprint, const double jump_beta,
    const double jump_theta, const double jump_z, const double jump_w,
    const double pr_mean_beta, const double pr_sd_beta,
    const double pr_mean_theta, double pr_sd_theta, const double pr_a_theta,
    const double pr_b_theta, const double pr_a_eps, const double pr_b_eps,
    const double missing, const bool verbose, const bool fix_theta_sd,
    Rcpp::Nullable<Rcpp::List> adapt = R_NilValue) {
  return lsirm1pl_internal(data, ndim, niter, nburn, nthin, nprint, jump_beta,
                           jump_theta, 0.0, jump_z, jump_w, pr_mean_beta,
                           pr_sd_beta, pr_mean_theta, pr_sd_theta, 0.0, 1.0,
                           0.0, 1.0, 0.0, 1.0, pr_a_theta, pr_b_theta, 1.0, 1.0,
                           pr_a_eps, pr_b_eps, missing, verbose, fix_theta_sd,
                           true, false, false, true, // mar, spike_slab, normal
                           adapt);
}
