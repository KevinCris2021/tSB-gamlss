// ============================================================================
// tsb_core.cpp  |  Numerical kernel for the t-Simplex-Binomial distribution
// Consolidated candidate version (GAMLSS / Rcpp), 2026-11-11
//
// Parametrization: Z | Y ~ Binomial(m,Y); Y ~ tS(mu,phi,nu)
// sigma (in GAMLSS) = phi. Parameters phi and nu are strictly positive.
// In particular, nu -> infinity yields the Simplex-Binomial limit.
//
// Rcpp exports:
//  - tSB_moments_cpp : log-PMF, score (3), observed Hessian (6)
//  - tSB_information_cpp : expected information (6), mass, E[score]
//
// Derivatives use Fisher/Louis identities rather than finite
// differences. Gauss-Legendre nodes and weights are explicit arguments.
// ============================================================================

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>

using namespace Rcpp;

namespace {

struct Quadrature {
  const NumericVector& nodes;
  const NumericVector& weights;
  std::vector<double> log_weights;
  std::vector<double> log_y1y;

  Quadrature(const NumericVector& n, const NumericVector& w) :
    nodes(n), weights(w), log_weights(n.size()), log_y1y(n.size()) {
    if (n.size() < 2 || n.size() != w.size()) {
      Rcpp::stop("The quadrature rule has invalid dimensions.");
    }
    double mass = 0.0;
    for (R_xlen_t j = 0; j < n.size(); ++j) {
      if (!std::isfinite(n[j]) || n[j] <= 0.0 || n[j] >= 1.0 ||
          !std::isfinite(w[j]) || w[j] <= 0.0) {
        Rcpp::stop("The quadrature nodes or weights are invalid.");
      }
      log_weights[j] = std::log(w[j]);
      log_y1y[j] = std::log(n[j]) + std::log1p(-n[j]);
      mass += w[j];
    }
    if (std::fabs(mass - 1.0) > 1e-8) {
      Rcpp::stop("Quadrature weights must sum to one.");
    }
  }
};

void check_parameters(double z, double m, double mu, double phi, double nu) {
  if (!std::isfinite(z) || !std::isfinite(m) ||
      !std::isfinite(mu) || !std::isfinite(phi) || !std::isfinite(nu) ||
      m < 1.0 || m != std::floor(m) ||
      z < 0.0 || z > m || z != std::floor(z) ||
      mu <= 0.0 || mu >= 1.0 || phi <= 0.0 || nu <= 0.0) {
    Rcpp::stop("Count or tSB parameters outside their domain.");
  }
}

struct Result {
  double log_prob;
  std::array<double, 3> score;
  // order: mu-mu, mu-phi, mu-nu, phi-phi, phi-nu, nu-nu
  std::array<double, 6> hessian;
};

Result evaluate_one(double z, double m, double mu, double phi, double nu,
                    const Quadrature& gl) {
  check_parameters(z, m, mu, phi, nu);
  const int K = static_cast<int>(gl.nodes.size());
  const double LOG2PI = std::log(2.0 * std::acos(-1.0));
  const double t = nu + 1.0;
  const double x = nu + 2.0;
  const double k = nu + 2.5;
  const double g = 1.0 + 1.5/t;
  const double phiinv = 1.0/phi;
  const double a = -2.0/mu + 2.0/(1.0-mu);
  const double a_prime = 2.0/(mu*mu) + 2.0/((1.0-mu)*(1.0-mu));
  const double mu_prod = mu*(1.0-mu);

  // Normalizing constant and derivatives with respect to nu.
  // For large nu, avoid lgamma(k)-lgamma(x), digamma, and trigamma.
  double logC, C1, C2;
  if (nu >= 1e4) {
    const double invx = 1.0/x;
    const double invt = 1.0/t;
    const double invx2 = invx*invx;
    const double invx3 = invx2*invx;
    const double invx4 = invx2*invx2;
    const double invx5 = invx4*invx;
    const double invx6 = invx3*invx3;
    const double invx7 = invx6*invx;
    logC = -0.5*(LOG2PI + std::log(phi)) +
      0.5*std::log1p(invt) - invx/8.0 + invx3/192.0 - invx5/640.0;
    C1 = -0.5*invt*invx + invx2/8.0 - invx4/64.0 + invx6/128.0;
    // This form avoids forming t^2 and x^2 as huge numbers.
    C2 = 0.5*(invt*invt*invx + invt*invx2) -
      invx3/4.0 + invx5/16.0 - 3.0*invx7/64.0;
  } else {
    logC = std::lgamma(k)-std::lgamma(x) -
      0.5*(LOG2PI + std::log(phi) + std::log(t));
    C1 = R::digamma(k)-R::digamma(x)-0.5/t;
    C2 = R::trigamma(k)-R::trigamma(x)+0.5/(t*t);
  }

  std::vector<double> lw(K);
  std::vector<std::array<double, 3>> sc(K);
  std::vector<std::array<double, 6>> hc(K);
  double maximum = -std::numeric_limits<double>::infinity();

  for (int j = 0; j < K; ++j) {
    const double u = gl.nodes[j];
    const double r = u - mu;
    const double c = 1.0/(u*(1.0-u)*mu_prod*mu_prod);
    const double d = c*r*r;
    const double dmu = c*(-2.0*r + a*r*r);
    const double dmumu = c*(2.0-4.0*a*r+(a*a+a_prime)*r*r);
    // Compute q without forming 2*phi*t if it would overflow.
    const double q = (d/(2.0*phi))/t;
    const double v = 1.0+q, v2=v*v;
    if (!std::isfinite(q) || !std::isfinite(v2)) {
      Rcpp::stop("Parameters/nodes too extreme for the tSB kernel.");
    }
    lw[j] = gl.log_weights[j] + R::dbinom(z,m,u,true) +
      logC - 1.5*gl.log_y1y[j] - k*std::log1p(q);
    if (!std::isfinite(lw[j])) {
      Rcpp::stop("Non-finite tSB integrand: check parameters or GL rule.");
    }
    maximum=std::max(maximum,lw[j]);

    sc[j][0] = -g*dmu/(2.0*phi*v);
    sc[j][1] = (-0.5+k*q/v)*phiinv;
    double dif;
    if (q < 1e-3) {
      // q/(1+q)-log(1+q) evaluated without cancellation.
      dif=q*q*(-0.5+q*(2.0/3.0+q*(-3.0/4.0+q*(4.0/5.0-5.0*q/6.0))));
    } else {
      dif=q/v-std::log1p(q);
    }
    sc[j][2]=C1+dif+1.5*q/(t*v);

    hc[j][0]=-g*dmumu/(2.0*phi*v)+
      (g/t)*(dmu/(2.0*phi))*(dmu/(2.0*phi))/v2;
    hc[j][1]=g*dmu/(2.0*phi*phi*v2);
    hc[j][2]=dmu*(1.5/t-q)/(2.0*phi*t*v2);
    hc[j][3]=(0.5-k*q*(q+2.0)/v2)/(phi*phi);
    hc[j][4]=q*(q-1.5/t)/(phi*v2);
    hc[j][5]=C2+q*(q-1.5*(2.0+q)/t)/(t*v2);
  }

  double sum_exp=0.0;
  for (int j=0;j<K;++j) sum_exp+=std::exp(lw[j]-maximum);
  if (!std::isfinite(sum_exp) || sum_exp<=0.0)
    Rcpp::stop("Could not normalize the tSB integrand.");

  Result out;
  out.log_prob=maximum+std::log(sum_exp);
  out.score.fill(0.0);
  out.hessian.fill(0.0);
  std::array<double,6> Ess={{0,0,0,0,0,0}};
  for (int j=0;j<K;++j) {
    const double p=std::exp(lw[j]-out.log_prob);
    for (int cidx=0;cidx<3;++cidx) out.score[cidx]+=p*sc[j][cidx];
    for (int hidx=0;hidx<6;++hidx) out.hessian[hidx]+=p*hc[j][hidx];
    const double s0=sc[j][0],s1=sc[j][1],s2=sc[j][2];
    Ess[0]+=p*s0*s0; Ess[1]+=p*s0*s1; Ess[2]+=p*s0*s2;
    Ess[3]+=p*s1*s1; Ess[4]+=p*s1*s2; Ess[5]+=p*s2*s2;
  }
  const double s0=out.score[0],s1=out.score[1],s2=out.score[2];
  const std::array<double,6> outer={{s0*s0,s0*s1,s0*s2,s1*s1,s1*s2,s2*s2}};
  for (int j=0;j<6;++j) {
    out.hessian[j]+=Ess[j]-outer[j];
    if (!std::isfinite(out.hessian[j]))
      Rcpp::stop("Non-finite tSB Hessian.");
  }
  for (int j=0;j<3;++j)
    if (!std::isfinite(out.score[j])) Rcpp::stop("Non-finite tSB score.");
  return out;
}

const char* hessian_names[6] =
  {"mumu", "muphi", "munu", "phiphi", "phinu", "nunu"};

} // namespace

// [[Rcpp::export]]
Rcpp::List tSB_moments_cpp(
    Rcpp::NumericVector z,
    Rcpp::NumericVector m,
    Rcpp::NumericVector mu,
    Rcpp::NumericVector phi,
    Rcpp::NumericVector nu,
    Rcpp::NumericVector nodes,
    Rcpp::NumericVector weights) {
  const R_xlen_t n = z.size();
  if (m.size() != n || mu.size() != n || phi.size() != n || nu.size() != n)
    Rcpp::stop("All five observation vectors must have the same length.");
  Quadrature gl(nodes, weights);
  NumericVector log_prob(n);
  NumericMatrix scores(n, 3), hess(n, 6);
  colnames(scores) = CharacterVector::create("mu", "phi", "nu");
  colnames(hess) = CharacterVector::create(
    "mumu", "muphi", "munu", "phiphi", "phinu", "nunu");

  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 500 == 0) Rcpp::checkUserInterrupt();
    const Result r = evaluate_one(z[i], m[i], mu[i], phi[i], nu[i], gl);
    log_prob[i] = r.log_prob;
    for (int j = 0; j < 3; ++j) scores(i,j) = r.score[j];
    for (int j = 0; j < 6; ++j) hess(i,j) = r.hessian[j];
  }
  return List::create(_["log_prob"] = log_prob,
                      _["score"] = scores,
                      _["hessian"] = hess);
}

// [[Rcpp::export]]
Rcpp::List tSB_information_cpp(
    Rcpp::NumericVector m,
    Rcpp::NumericVector mu,
    Rcpp::NumericVector phi,
    Rcpp::NumericVector nu,
    Rcpp::NumericVector nodes,
    Rcpp::NumericVector weights) {
  const R_xlen_t n = m.size();
  if (mu.size() != n || phi.size() != n || nu.size() != n)
    Rcpp::stop("m, mu, phi, and nu must have the same length.");
  Quadrature gl(nodes, weights);
  NumericMatrix info(n, 6), escore(n, 3);
  NumericVector mass(n);
  colnames(info) = CharacterVector::create(
    "mumu", "muphi", "munu", "phiphi", "phinu", "nunu");
  colnames(escore) = CharacterVector::create("mu", "phi", "nu");

  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    if (!std::isfinite(m[i]) || m[i] < 1.0 || m[i] != std::floor(m[i]) ||
        m[i] > static_cast<double>(std::numeric_limits<int>::max()-1))
      Rcpp::stop("Invalid or excessively large number of trials.");
    const int mi = static_cast<int>(m[i]);
    double I[6] = {0,0,0,0,0,0};
    double s[3] = {0,0,0};
    double probability_mass = 0.0;
    for (int z = 0; z <= mi; ++z) {
      const Result r = evaluate_one(z, m[i], mu[i], phi[i], nu[i], gl);
      const double p = std::exp(r.log_prob);
      probability_mass += p;
      const double s0 = r.score[0], s1 = r.score[1], s2 = r.score[2];
      s[0] += p*s0; s[1] += p*s1; s[2] += p*s2;
      I[0] += p*s0*s0; I[1] += p*s0*s1; I[2] += p*s0*s2;
      I[3] += p*s1*s1; I[4] += p*s1*s2; I[5] += p*s2*s2;
    }
    if (!std::isfinite(probability_mass) ||
        std::fabs(probability_mass - 1.0) > 1e-6)
      Rcpp::stop("Quadrature not normalized: increase the number of GL nodes.");
    mass[i] = probability_mass;
    for (int j = 0; j < 6; ++j) info(i,j) = I[j];
    for (int j = 0; j < 3; ++j) escore(i,j) = s[j];
  }
  return List::create(_["information"] = info,
                      _["expected_score"] = escore,
                      _["mass"] = mass);
}
