% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsirm1pl_mcar_ss.R
\name{lsirm1pl_mcar_ss}
\alias{lsirm1pl_mcar_ss}
\title{1pl LSIRM model with model selection approach for missing completely at random data.}
\usage{
lsirm1pl_mcar_ss(
  data,
  ndim = 2,
  niter = 15000,
  nburn = 2500,
  nthin = 5,
  nprint = 500,
  jump_beta = 0.4,
  jump_theta = 1,
  jump_gamma = 1,
  jump_z = 0.5,
  jump_w = 0.5,
  pr_mean_beta = 0,
  pr_sd_beta = 1,
  pr_mean_theta = 0,
  pr_spike_mean = -3,
  pr_spike_sd = 1,
  pr_slab_mean = 0.5,
  pr_slab_sd = 1,
  pr_a_theta = 0.001,
  pr_b_theta = 0.001,
  pr_xi_a = 1,
  pr_xi_b = 1,
  missing = 99
)
}
\arguments{
\item{data}{Matrix; binary item response matrix to be analyzed. Each row is assumed to be respondent and its column values are assumed to be response to the corresponding item.}

\item{ndim}{Numeric; dimension of latent space. default value is 2.}

\item{niter}{Numeric; number of iterations to run MCMC sampling. default value is 15000.}

\item{nburn}{Numeric; number of initial, pre-thinning, MCMC iterations to discard. default value is 2500.}

\item{nthin}{Numeric;number of thinning, MCMC iterations to discard. default value is 5.}

\item{nprint}{Numeric; MCMC samples is displayed during execution of MCMC chain for each \code{nprint}. default value is 500.}

\item{jump_beta}{Numeric; jumping rule of the proposal density for beta. default value is 0.4.}

\item{jump_theta}{Numeric; jumping rule of the proposal density for theta. default value is 1.0.}

\item{jump_gamma}{Numeric; jumping rule of the proposal density for gamma. default value is 1.0.}

\item{jump_z}{Numeric; jumping rule of the proposal density for z. default value is 0.5.}

\item{jump_w}{Numeric; jumping rule of the proposal density for w. default value is 0.5.}

\item{pr_mean_beta}{Numeric; mean of normal prior for beta. default value is 0.}

\item{pr_sd_beta}{Numeric; standard deviation of normal prior for beta. default value is 1.0}

\item{pr_mean_theta}{Numeric; mean of normal prior for theta. default value is 0.}

\item{pr_spike_mean}{Numeric; mean of spike prior for log gamma default value is -3.}

\item{pr_spike_sd}{Numeric; standard deviation of spike prior for log gamma default value is 1.}

\item{pr_slab_mean}{Numeric; mean of spike prior for log gamma default value is 0.5.}

\item{pr_slab_sd}{Numeric; standard deviation of spike prior for log gamma default value is 1.}

\item{pr_a_theta}{Numeric; shape parameter of inverse gamma prior for variance of theta. default value is 0.001.}

\item{pr_b_theta}{Numeric; scale parameter of inverse gamma prior for variance of theta. default value is 0.001.}

\item{pr_xi_a}{Numeric; first shape parameter of beta prior for latent variable xi. default value is 1.}

\item{pr_xi_b}{Numeric; second shape parameter of beta prior for latent variable xi. default value is 1.}

\item{missing}{Numeric; a number to replace missing values. default value is 99.}
}
\value{
\code{lsirm1pl_mcar_ss} returns an object of  list containing the following components:
 \item{beta_estimate}{posterior estimation of beta.}
 \item{theta_estimate}{posterior estimation of theta.}
 \item{sigma_theta_estimate}{posterior estimation of standard deviation of theta.}
 \item{gamma_estimate}{posterior estimation of gamma.}
 \item{z_estimate}{posterior estimation of z.}
 \item{w_estimate}{posterior estimation of w.}
 \item{pi_estimate}{posterior estimation of phi. inclusion probability of gamma. if estimation of phi is less than 0.5, choose Rasch model with gamma = 0, otherwise latent space model with gamma > 0. }
 \item{beta}{posterior samples of beta.}
 \item{theta}{posterior samples of theta.}
 \item{theta_sd}{posterior samples of standard deviation of theta.}
 \item{gamma}{posterior samples of gamma.}
 \item{z}{posterior samples of z. The output is 3-dimensional matrix with last axis represent the dimension of latent space.}
 \item{w}{posterior samples of w. The output is 3-dimensional matrix with last axis represent the dimension of latent space.}
 \item{pi}{posterior samples of phi which is indicator of spike and slab prior. If phi is 1, log gamma follows the slab prior, otherwise follows the spike prior. }
 \item{accept_beta}{accept ratio of beta.}
 \item{accept_theta}{accept ratio of theta.}
 \item{accept_w}{accept ratio of w.}
 \item{accept_z}{accept ratio of z.}
 \item{accept_gamma}{accept ratio of gamma.}
}
\description{
\link{lsirm1pl_mcar_ss} is used to fit LSIRM model with model selection approach based on spike-and-slab priors in incomplete data assumed to be missing completely at random. 
\link{lsirm1pl_mcar_ss} factorizes item response matrix into column-wise item effect, row-wise respondent effect and further embeds interaction effect in a latent space, while ignoring the missing element under the assumption of missing completely at random. The resulting latent space provides an interaction map that represents interactions between respondents and items.
}
\details{
\code{lsirm1pl_mcar_ss} models the probability of correct response by respondent \eqn{j} to item \eqn{i} with item effect \eqn{\beta_i}, respondent effect \eqn{\theta_j} and the distance between latent position \eqn{w_i} of item \eqn{i} and latent position \eqn{z_j} of respondent \eqn{j} in the shared metric space, with \eqn{\gamma} represents the weight of the distance term: \deqn{logit(P(Y_{j,i} = 1 |\theta_j,\beta_i,\gamma,z_j,w_i))=\theta_j+\beta_i-\gamma||z_j-w_i||} Under the assumption of missing completely at random, the model ignores the missing element in doing inference. For the details of missing completely at random assumption and data augmentation, see References. \code{lsirm1pl_mcar_ss} model include model selection approach based on spike-and-slab priors for log gamma. For detail of spike-and-slab priors, see References.
}
\examples{
\donttest{
# generate example item response matrix
data     <- matrix(rbinom(500, size = 1, prob = 0.5),ncol=10,nrow=50)

# generate example missing indicator matrix
missing_mat     <- matrix(rbinom(500, size = 1, prob = 0.2),ncol=10,nrow=50)

# make missing value with missing indicator matrix
data[missing_mat==1] <- 99 

lsirm_result <- lsirm1pl_mcar_ss(data)

# The code following can achieve the same result.
lsirm_result <- lsirm1pl(data, spikenslab = TRUE, fixed_gamma = FALSE, missing_data = 'mcar')
}

}
\references{
Little, R. J., & Rubin, D. B. (2019). Statistical analysis with missing data (Vol. 793). John Wiley & Sons.
Ishwaran, H., & Rao, J. S. (2005). Spike and slab variable selection: frequentist and Bayesian strategies. The Annals of Statistics, 33(2), 730-773.
}