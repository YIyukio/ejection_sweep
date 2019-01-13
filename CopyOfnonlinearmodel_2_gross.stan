functions{
  vector Arrhenius(real EA, vector kelvin, real kelvin_mean){
    real gas_constant;
    gas_constant = 8.617/100; // boltzman constant
    return exp(- EA / gas_constant * (1 ./kelvin - 1/kelvin_mean));
  }
  
  vector modelcurve(real b0, real b1, real b2, real b3, vector speed) {
    return (b0 + b1 * (1-exp(-b2 / b1 * speed)) .* exp(-b3 / b1 * speed));
    //return (b0 + b1 * speed ./ (b2 + speed) - b3 * speed);
  }
}
data {
  int<lower=0> N;                 // データの数
  int<lower=0> M;                 // 期待値に使う結果の数
  vector[N] speed;                // 観測した速度
  vector[N] gross;                  // 観測した総一次生産量
  vector[4] prior_log_sigma;      // 各パラメータの事前分散
  vector[4] prior_log_mu;         // 各パラメータの事前期待値
}
transformed data {
  vector[N] gross10;
  vector[M] speed_predict;
  real speed_increment;
  gross10 = gross*10;
  speed_increment = max(speed) / M;
  speed_predict[1] = 0;
  
  for(m in 2:M) {
    speed_predict[m] = speed_predict[m-1] + speed_increment;
  }
}
parameters {
  real<lower=0> offset;
  real<lower=0> pmax;
  real<lower=0> alpha;
  real<lower=0> beta;
  
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] fitted;
  fitted = modelcurve(offset, pmax, alpha, beta, speed);
}
model {
  offset ~ lognormal(prior_log_mu[1], prior_log_sigma[1]);
  pmax   ~ lognormal(prior_log_mu[2], prior_log_sigma[2]);
  alpha  ~ lognormal(prior_log_mu[3], prior_log_sigma[3]);
  beta   ~ lognormal(prior_log_mu[4], prior_log_sigma[4]);
  sigma ~ cauchy(0, 2.5);
  gross10 ~ lognormal(log(fitted), sigma);
}
generated quantities {
  vector[M] expected10;
  vector[M] expected;
  vector[N] log_lik;
  expected10 = modelcurve(offset, pmax, alpha, beta, speed_predict);
  expected = expected10/10;
  for(n in 1:N) {
    log_lik[n] = lognormal_lpdf(log(gross10[n]) | log(fitted[n]), sigma);
  }
}


