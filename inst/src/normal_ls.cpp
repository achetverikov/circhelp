#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y);
  DATA_MATRIX(X_mu);
  DATA_MATRIX(X_sigma);
  DATA_VECTOR(weights);

  DATA_MATRIX(D_mu);
  DATA_SCALAR(lambda_mu);

  PARAMETER_VECTOR(beta_mu);
  PARAMETER_VECTOR(beta_sigma);

  vector<Type> mu = X_mu * beta_mu;
  vector<Type> eta_sigma = X_sigma * beta_sigma;

  Type nll = 0.0;

  for (int i = 0; i < y.size(); i++) {
    if (weights(i) > Type(0)) {
      Type sigma = exp(eta_sigma(i));
      nll -= weights(i) * dnorm(y(i), mu(i), sigma, true);
    }
  }

  if (lambda_mu > Type(0) && D_mu.rows() > 0) {
    vector<Type> Db = D_mu * beta_mu;
    nll += Type(0.5) * lambda_mu * (Db * Db).sum();
  }

  REPORT(mu);
  REPORT(eta_sigma);

  ADREPORT(mu);
  ADREPORT(eta_sigma);

  return nll;
}
