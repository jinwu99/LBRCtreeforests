#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::export]]
List vardiCpp(NumericMatrix dat, NumericVector w,
                           double eps = 1e-7, int max_iter = 100) {
  int n = dat.nrow();

  NumericVector y(n);
  IntegerVector d(n);
  for(int i = 0; i < n; i++) {
    y[i] = dat(i, 1);        // R’s dat[,2]
    d[i] = static_cast<int>(dat(i, 2)); // R’s dat[,3]
  }

  // build unique, sorted event times t
  std::vector<double> tvec(y.begin(), y.end());
  std::sort(tvec.begin(), tvec.end());
  tvec.erase(std::unique(tvec.begin(), tvec.end()), tvec.end());
  int k = tvec.size();
  NumericVector t(k);
  for(int j = 0; j < k; j++) t[j] = tvec[j];

  // initialize EM variables
  NumericVector p(k, 1.0/k), p_old(k), denom(n), ratio(n), q(k), S(k);
  // at‐risk indicator matrix M[i][j] = (y[i] <= t[j])
  std::vector< std::vector<char> > M(n, std::vector<char>(k));
  for(int i = 0; i < n; i++)
    for(int j = 0; j < k; j++)
      M[i][j] = (y[i] <= t[j]);

  // count events at each unique time
  NumericVector events(k, 0.0);
  for(int i = 0; i < n; i++) if(d[i] == 1) {
    for(int j = 0; j < k; j++) {
      if(y[i] == t[j]) { events[j] += w[i]; break; }
    }
  }

  // total weight
  double sum_w = std::accumulate(w.begin(), w.end(), 0.0);

  // EM loop
  for(int iter = 0; iter < max_iter; iter++) {
    p_old = clone(p);

    // a) compute denom[i] = sum_j (p[j]/t[j]) * M[i][j] * w[i]
    for(int i = 0; i < n; i++) {
      double sum = 0;
      for(int j = 0; j < k; j++)
        sum += (p[j]/t[j]) * M[i][j];
      denom[i] = sum;

      // handle zero‐denom just in case
      if(denom[i] <= 0) denom[i] = 1e-16;
    }

    // b) ratio[i] = w[i] * (1 - d[i]) / denom[i]
    for(int i = 0; i < n; i++) {
      ratio[i] = w[i] * (1 - d[i]) / denom[i];
    }

    // c) update p[j]
    for(int j = 0; j < k; j++) {
      double sum = 0;
      for(int i = 0; i < n; i++)
        sum += ratio[i] * M[i][j];
      p[j] = (events[j] + (p[j]/t[j]) * sum) / sum_w;
    }

    // d) check convergence
    double err = 0;
    for(int j = 0; j < k; j++)
      err += std::abs(p[j] - p_old[j]);
    if(err <= eps) break;
  }

  // form the survival curve S_pred
  for(int j = 0; j < k; j++)
    q[j] = p[j] / t[j];
  double sumq = std::accumulate(q.begin(), q.end(), 0.0);
  double cum = 0;
  for(int j = 0; j < k; j++) {
    cum    += q[j] / sumq;
    S[j]    = 1.0 - cum;
    if(S[j] < 0) S[j] = 0;
  }

  // Return both t, p, and S
  return List::create(Named("t") = t,
                      Named("p") = q,
                      Named("S") = S);
}
