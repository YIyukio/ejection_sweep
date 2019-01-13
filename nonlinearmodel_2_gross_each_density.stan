functions{
  real modelcurve(real b0, real b1, real b2, real b3, real speed) {
    return (b0 + b1 * (1-exp(-b2 / b1 * speed)) * exp(-b3 / b1 * speed));
  }
}

data {
  int<lower=0> N;                 // データの数
  int<lower=0> M;                 // 期待値に使う結果の数
  int idx[N];                     // Densityのインデックス
  vector[N] speed;                // 観測した速度
  vector[N] gross;                  // 観測した総一次生産量
  vector[M] speed_predict;
  vector[4] prior_log_sigma;      // 各パラメータの事前分散
  vector[4] prior_log_mu;         // 各パラメータの事前期待値
}

transformed data {
  vector[N] gross100;
  gross100 = gross*100;
}

parameters {
  vector<lower=0>[3] offset;
  vector<lower=0>[3] pmax;
  vector<lower=0>[3] alpha;
  vector<lower=0>[3] beta;
  vector<lower=0>[3] sigma;
}

model {
  offset ~ lognormal(prior_log_mu[1], prior_log_sigma[1]);
  pmax   ~ lognormal(prior_log_mu[2], prior_log_sigma[2]);
  alpha  ~ lognormal(prior_log_mu[3], prior_log_sigma[3]);
  beta   ~ lognormal(prior_log_mu[4], prior_log_sigma[4]);
  sigma ~ cauchy(0, 2.5);
  
  for(i in 1:N){
  gross100[i] ~ lognormal(log(modelcurve(offset[idx[i]], pmax[idx[i]], alpha[idx[i]], beta[idx[i]], speed[i])), sigma[idx[i]]);  
  }
}

generated quantities {
  vector[M] yfit_01_100;
  vector[M] yfit_01;
  vector[M] yfit_02_100;
  vector[M] yfit_02;
  vector[M] yfit_03_100;
  vector[M] yfit_03;
  //vector[3] sigma;
  vector[M] ypred_01;
  vector[M] ypred_02;
  vector[M] ypred_03;
  
  for(n in 1:M){
  yfit_01_100[n] = modelcurve(offset[1], pmax[1], alpha[1], beta[1], speed_predict[n]);
  yfit_02_100[n] = modelcurve(offset[2], pmax[2], alpha[2], beta[2], speed_predict[n]);
  yfit_03_100[n] = modelcurve(offset[3], pmax[3], alpha[3], beta[3], speed_predict[n]);
  }
  
  yfit_01 = yfit_01_100/100;
  yfit_02 = yfit_02_100/100;
  yfit_03 = yfit_03_100/100;
  //sigma = sigma100/100;
  
  for(n in 1:M){
  ypred_01[n] = lognormal_rng(log(yfit_01_100[n]), sigma[1]);
  ypred_02[n] = lognormal_rng(log(yfit_02_100[n]), sigma[2]);
  ypred_03[n] = lognormal_rng(log(yfit_03_100[n]), sigma[3]);
  }
  ypred_01 = ypred_01 / 100;
  ypred_02 = ypred_02 / 100;
  ypred_03 = ypred_03 / 100;
}
