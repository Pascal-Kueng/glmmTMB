#ifndef GLMMTMB_KRON_H
#define GLMMTMB_KRON_H

#include <limits>
#include <vector>

namespace glmmtmb {
namespace kron {

template <class Type>
using margin_type = density::GMRF_t<Type>;

/*
 * Gaussian distribution with covariance
 *
 *   scale^2 * (Sigma[M-1] (x) ... (x) Sigma[0]).
 *
 * `dim[0]` is the fastest-running array dimension, as in R and TMB arrays.
 * Each marginal GMRF supplies its precision action, log determinant, and
 * covariance square root.  The full Kronecker matrix is never constructed.
 */
template <class Type>
class GaussianProduct {
public:
  typedef glmmtmb::kron::margin_type<Type> margin_distribution;

private:
  vector<int> dim_;
  int size_;
  Type scale_;
  std::vector<margin_distribution> margin_;
  Type base_normalizer_;

  static int checked_size(const vector<int>& dim, int n_margin) {
    if (dim.size() == 0)
      error("Kronecker Gaussian needs at least one margin");
    if (dim.size() != n_margin)
      error("Kronecker Gaussian dimensions and margins do not match");

    long long n = 1;
    for (int m = 0; m < dim.size(); ++m) {
      if (dim(m) <= 0)
        error("Kronecker Gaussian dimensions must be positive");
      n *= dim(m);
      if (n > std::numeric_limits<int>::max())
        error("Kronecker Gaussian is too large for a TMB array");
    }
    return static_cast<int>(n);
  }

  void initialize_normalizer() {
    base_normalizer_ = Type(0);
    for (int m = 0; m < dim_.size(); ++m) {
      vector<Type> zero(dim_(m));
      zero.setZero();
      base_normalizer_ += Type(size_ / dim_(m)) * margin_[m](zero);
    }
    base_normalizer_ -= Type(dim_.size() - 1) * Type(size_) *
      Type(log(sqrt(2.0 * M_PI)));
  }

public:
  GaussianProduct(const vector<int>& dim,
                  const std::vector<margin_distribution>& margin,
                  Type scale = Type(1))
    : dim_(dim), size_(checked_size(dim, margin.size())), scale_(scale),
      margin_(margin) {
    initialize_normalizer();
  }

  Type operator()(vector<Type> x) {
    if (x.size() != size_)
      error("Kronecker Gaussian vector has the wrong length");
    vector<Type> z = x / scale_;
    array<Type> z_array(z, dim_);
    array<Type> qz = z_array;
    /* GMRF::jacobian acts on the last dimension. Rotating after each
       reverse-order application visits every margin and restores the order. */
    for (int m = dim_.size() - 1; m >= 0; --m) {
      qz = margin_[m].jacobian(qz);
      qz = qz.rotate(1);
    }
    return base_normalizer_ + Type(size_) * log(scale_) +
      Type(0.5) * (z_array * qz).sum();
  }

  /* Transform an iid standard-normal vector by a covariance square root. */
  vector<Type> sqrt_cov_scale(vector<Type> u) {
    if (u.size() != size_)
      error("Kronecker Gaussian vector has the wrong length");
    array<Type> x(u, dim_);
    for (int m = dim_.size() - 1; m >= 0; --m) {
      const int d = dim_(m);
      const int n_line = size_ / d;
      vector<Type> line(d);
      for (int i = 0; i < n_line; ++i) {
        for (int j = 0; j < d; ++j)
          line(j) = x(i + n_line * j);
        line = margin_[m].sqrt_cov_scale(line);
        for (int j = 0; j < d; ++j)
          x(i + n_line * j) = line(j);
      }
      x = x.rotate(1);
    }
    return scale_ * x.vec();
  }

  vector<Type> simulate() {
    vector<Type> u(size_);
    density::rnorm_fill(u);
    return density::zero_derivatives(sqrt_cov_scale(u));
  }
};

} // namespace kron
} // namespace glmmtmb

#endif
