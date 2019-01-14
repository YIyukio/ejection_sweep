functions{
  vector npcurve(real b0, real b1, real b2, real b3, vector speed) {
    return (b0 + b1 * (1-exp(-b2 / b1 * speed))) .* exp(-b3/b1 * speed);
  }
  vector rpcurve(real b0, real b1, vector speed) {
    return (b0 + b1 * speed);
  }
  real npcurve_real(real b0, real b1, real b2, real b3, real speed) {
    return (b0 + b1 * (1-exp(-b2 / b1 * speed))) .* exp(-b3/b1 * speed);
  }
  real rpcurve_real(real b0, real b1, real speed) {
    return (b0 + b1 * speed);
  }
}

data {
  int<lower=0> N;                 // データの数
  int<lower=0> M;                 // 期待値に使う結果の数
  vector[N] speed;                // 観測した速度
  vector[N] net;                  // 観測した総一次生産量
  vector[N] resp;                  // 観測した総一次生産量
  vector[M] speed_predict;
  vector[6] prior_log_sigma;      // 各パラメータの事前分散
  vector[6] prior_log_mu;         // 各パラメータの事前期待値
}

transformed data {
  vector[N] net100;
  vector[N] resp100;
  net100 = net*100;
  resp100 = resp*100;
}

parameters {
  real<lower=0> offset_np;
  real offset_rp;
  real<lower=0> pmax;
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma_np;
  real<lower=0> sigma_rp;
  real slope;
}

transformed parameters {
  vector[N] fitted_np100;
  vector[N] fitted_rp100;
  fitted_np100 = npcurve(offset_np, pmax, alpha, beta, speed);
  fitted_rp100 = rpcurve(offset_rp, slope, speed);
}

model {
  offset_rp ~ normal(prior_log_mu[1], prior_log_sigma[1]);
  slope ~     normal(prior_log_mu[2], prior_log_sigma[2]);

  offset_np ~ normal(prior_log_mu[3], prior_log_sigma[3]);
  pmax   ~    normal(prior_log_mu[4], prior_log_sigma[4]);
  alpha  ~ lognormal(prior_log_mu[5], prior_log_sigma[5]);
  beta   ~ lognormal(prior_log_mu[6], prior_log_sigma[6]);
  sigma_np ~ cauchy(0, 2.5);
  sigma_rp ~ cauchy(0, 2.5);

  net100 ~ normal(fitted_np100, sigma_np);
  resp100 ~ normal(fitted_rp100, sigma_rp);
}

generated quantities {
  vector[M] yfit_np_100;
  vector[M] yfit_np;
  vector[M] ypred_np;

  vector[M] yfit_rp_100;
  vector[M] yfit_rp;
  vector[M] ypred_rp;

  vector[M] yfit_gp;

  vector[2*N] log_lik;

  yfit_np_100 = npcurve(offset_np, pmax, alpha, beta, speed_predict);
  yfit_rp_100 = rpcurve(offset_rp, slope, speed_predict);

  yfit_np = yfit_np_100/100;
  yfit_rp = yfit_rp_100/100;

  for(m in 1:M){
    ypred_np[m] = normal_rng(yfit_np_100[m], sigma_np);
    ypred_rp[m] = normal_rng(yfit_rp_100[m], sigma_rp);
  }
  ypred_np = ypred_np / 100;
  ypred_rp = ypred_rp / 100;
  for(m in 1:M) {
    yfit_gp[m] = yfit_np[m] + yfit_rp[m];
  }

  for(n in 1:N) {
    log_lik[n] = normal_lpdf(net100[n] | npcurve_real(offset_np, pmax, alpha, beta, speed[n]), sigma_np);
  }
  for(n in (N+1):(2*N)) {
    log_lik[n] = normal_lpdf(resp100[n-N] | rpcurve_real(offset_rp,slope, speed[n-N]), sigma_rp);
  }

}
