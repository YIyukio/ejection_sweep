library(tidyverse)
library(lubridate)
library(zoo)

# Declare the function --------------------------
# xには流速の残差が代入される。差分法による微分値を算出する。
gradient = function(x){
  n = length(x)
  j = 2:(n-1)
  c(x[2] - x[1],
    sapply(j, function(j, z){(z[j+1] - z[j-1])*0.5}, z = x),
    x[n] - x[n-1])
}
#advdata には、流速の生データが代入される。残りの引数はそれぞれの閾値を示す。
# 閾値は、論文より引用して決めるが、便宜上以下の値を使用。
# 返り値は、閾値から外れたデータを NA に変化したもの。
snr_corr_filter = function(advdata, snr.threshold=2, corr.threshold=60) {
  # A filter for the SNR and CORR.
  require(tidyverse)
  advdata  =
    advdata %>%
    mutate_at(vars(starts_with("snr")), funs(na=. < snr.threshold)) %>%
    mutate_at(vars(starts_with("corr")), funs(na=. < corr.threshold)) %>%
    mutate(xvel  = ifelse(snr1_na, NA, xvel)) %>%
    mutate(yvel  = ifelse(snr2_na, NA, yvel)) %>%
    mutate(zvel  = ifelse(snr3_na, NA, zvel)) %>%
    mutate(z2vel = ifelse(snr4_na, NA, z2vel)) %>%
    mutate(xvel  = ifelse(corr1_na, NA, xvel)) %>%
    mutate(yvel  = ifelse(corr2_na, NA, yvel)) %>%
    mutate(zvel  = ifelse(corr3_na, NA, zvel)) %>%
    mutate(z2vel = ifelse(corr4_na, NA, z2vel)) %>%
    mutate_at(vars(contains("vel")), funs(out = zoo::na.locf(., na.rm=F))) %>%
    rename(u = xvel_out, v = yvel_out, w = zvel_out, w2 = z2vel_out) %>%
    mutate(tau = 1:length(xvel))
}

interpolate = function(value) {
  # Conducts a cubic spline interpolation of removed spikes.
  x = which(!is.na(value))
  y = value[x]
  xi = seq_len(length(value))
  spline(x,y, xout=xi, method="fmm")$y
}
ellipse_radius=function(ab, theta, alpha=0, x0=0, y0=0) {
  # Use the polar equation of the ellipse to determine the radius.
  x=x0+ab[1]*cos(theta)*cos(alpha)-ab[2]*sin(theta)*sin(alpha)
  y=y0+ab[1]*cos(theta)*sin(alpha)+ab[2]*sin(theta)*cos(alpha)
  radius=ab[1]*ab[2]/(sqrt(ab[1]^2*sin(theta-alpha)^2+ab[2]^2*cos(theta-alpha)^2))
  data_frame(x,y,radius)
}
distance_radius=function(x,y, ab, alpha=0) {
  # Determine the distance to the point and the distance (radius)
  # to the edge of the ellipse
  theta = atan(y/x)
  distance = sqrt(x^2+y^2)
  theta[x<0] = theta[x<0]+pi
  radius = ellipse_radius(ab, theta, alpha=alpha)
  data_frame(theta, distance, radius=radius$radius)
}
make_ellipse = function(ab,alpha=0,x0=0,y0=0,size=500) {
  theta=seq(0,2*pi,length=size)
  x=x0+ab[1]*cos(theta)*cos(alpha)-ab[2]*sin(theta)*sin(alpha)
  y=y0+ab[1]*cos(theta)*sin(alpha)+ab[2]*sin(theta)*cos(alpha)
  data_frame(x,y)
}
# 残差と残差の微分値(一階、二階)の相関図は、楕円形となる。
# 楕円作図用関数
plotdata=function(f,ft,ftt,f_ft, ft_ftt, f_ftt,f_ft_bad, ft_ftt_bad, f_ftt_bad, theta) {
  require(ggplot2)
  ell_f_ft = make_ellipse(f_ft,0)
  ell_ft_ftt = make_ellipse(ft_ftt,0)
  ell_f_ftt = make_ellipse(f_ftt,theta)
  plot1=ggplot(data.frame(f,ft)) +
    geom_point(aes(f,ft, color=f_ft_bad)) +
    geom_path(aes(x,y), data=ell_f_ft)
  plot2=ggplot(data.frame(ft,ftt)) +
    geom_point(aes(ft,ftt, color=ft_ftt_bad)) +
    geom_path(aes(x,y), data=ell_ft_ftt)
  plot3=ggplot(data.frame(f,ftt)) +
    geom_point(aes(f,ftt, color=f_ftt_bad)) +
    geom_path(aes(x,y), data=ell_f_ftt)
  gridExtra::grid.arrange(plot1,plot2,plot3, ncol=2)
}
MAD = function(x) {
  # See Wahl (2003) Discussion of "Despiking Acoustic Doppler Velocimeter Data"
  # by Derek G. Goring and Vladimir I. Nikora, Journal of Hydraulic Engineering, 129(6), 484-487.
  # and Goring and Nikora (2003) Closure to "Depiking Acoustic Doppler
  # Velocimeter Data" by Derek G. Goring and Vladimir I. Nikora,
  # Journal of Hydraulic Engineering, 129(6), 487-488.
  # This function calculates the median absolute deviation. The value 1.483
  # makes it similar to the standard deviation.
  sigma = median(abs(x - median(x, na.rm=T)), na.rm=T)
  k = 1.483 # for normal distribution
  k*sigma
}
despike=function(f, good_data, figures=TRUE) {
  ft = gradient(f)
  ftt = gradient(ft)
  theta = atan(sum(f*ftt)/sum(f^2))
  fsigma = MAD(f)
  ftsigma = MAD(ft)
  fttsigma = MAD(ftt)
  nobs = length(f)
  lambda = sqrt(2*log(nobs)) # Universal threshold
  # Set the minor and major axes for the ellipse
  f_ft = c(lambda*fsigma, lambda*ftsigma)
  f_ftt = c(lambda*fsigma, lambda*fttsigma)
  ft_ftt = c(lambda*ftsigma, lambda*fttsigma)
  # Determine the distance and radius of the ellipse for each point
  f_ft_tr = distance_radius(f, ft, f_ft)
  f_ftt_tr = distance_radius(f, ftt, f_ftt, alpha=theta)
  ft_ftt_tr = distance_radius(ft, ftt, ft_ftt)
  # Determine of the point is beyond the ellipse
  f_ft_bad = apply(f_ft_tr [, c("radius","distance")], 1, diff) > 0
  f_ftt_bad = apply(f_ftt_tr [, c("radius","distance")], 1, diff) > 0
  ft_ftt_bad = apply(ft_ftt_tr[, c("radius","distance")], 1, diff) > 0
  # Flag the spikes to remove
  makeNA = unique(sort(c(which(f_ft_bad), which(f_ftt_bad), which(ft_ftt_bad))))
  f_tmp=f
  makeNA = makeNA[!(makeNA %in% good_data)]
  f[(makeNA)] = NA
  f[1] = ifelse(is.na(f[1]), 0, f[1]) # 初期値が NA の場合エラーを返すので NA の場合は便宜上 0 を代入
  # f_new = interpolate(f)
  f_new = zoo::na.locf(f,na.rm=FALSE)
  if(is.na(f_new[1])) f_new[1] = f_new[2]
  if(figures) {
    plotdata(f_tmp, ft, ftt,
             f_ft, ft_ftt, f_ftt,
             f_ft_bad, ft_ftt_bad, f_ftt_bad, theta)
  }
  mat = matrix(c(f,ft,ftt), ncol=3)
  mat2 = matrix(c(f_ft, f_ftt, ft_ftt), ncol=2, byrow=T)
  list(f = f_new, N = length(makeNA), navalues = makeNA, mat = mat, mat2 = mat2, theta=theta)
}

nk_despike = function(f, just_despike=TRUE, quiet=TRUE) {
  f_mean = mean(f, na.rm=T)
  i=1
  Nold=100
  deltaN = 100
  lambda = sqrt(2*log(length(f)))
  C1 = 1.483
  C2 = 1.483
  f = f - f_mean
  fsigma = MAD(f)
  um = C2 * fsigma * lambda
  good_data = which(f >= -C1*fsigma, f <= C1*fsigma)
  remove = which(abs(f) > um)
  f[remove] = NA
  f[1] = ifelse(is.na(f[1]), 0, f[1]) # 初期値が NA の場合エラーを返すので NA の場合は便宜上 0 を代入
  f = zoo::na.locf(f, na.rm=FALSE)
  f = f + f_mean
  number_of_spikes = rep(0, 10)
  while(deltaN > 0 & i < 10) {
    f = f - f_mean
    out=despike(f, good_data, figure=F)
    # length(out$f)
    # dim(out$mat)
    number_of_spikes[i] = out$N
    Nnew = out$N
    f=out$f + f_mean
    # f_mean = mean(f, na.rm=T)
    f_mean = median(f, na.rm=T)
    deltaN = abs(Nold - Nnew)
    Nold = Nnew
    if(!quiet) {cat(paste0(i, ": Number of spikes detected: ", Nnew, "\n"))}
    i = i+1
  }
  spikes_removed = number_of_spikes[1]-number_of_spikes[i-1]
  if(!quiet) {
    cat(paste0("Spikes identified: ", number_of_spikes[1],
               "\nProportion of spikes: ", round((number_of_spikes[1] / length(f)),4)*100,
               "% of data.\n",
               "Spikes removed: ", spikes_removed))}
  if(just_despike) {
    f
  } else {
    list(f=f, out=out)
  }
}



cnames_vectrino=c("time", "counter", "status",
                  "xvel", "yvel", "zvel", "z2vel",
                  "amp1", "amp2", "amp3", "amp4",
                  "snr1", "snr2", "snr3", "snr4",
                  "corr1", "corr2", "corr3", "corr4")
ct = "nncnnnnnnnnnnnnnnnn"
# vector = lapply(vector_names, read_table , col_names=cnames_vectrino)
# vector_test1 = vector("list", length(vector_names))
# vector_test2 = vector("list", length(vector_names))

dset = data_frame(fnames = dir("Data/Velocity/", pattern = "dat", full =T))

dset =
  dset %>%
  mutate(data = map(fnames, read_table, col_names = cnames_vectrino, col_types = ct))

dset2 = dset %>%
  mutate(data = map(data, snr_corr_filter)) %>%
  mutate(data = map(data, function(X) {
    X %>%
      mutate_at(vars(u,v,w,w2), funs(out = nk_despike)) %>%
      select(time, tau, xvel:z2vel, u:w2, u_out:w2_out)
  })) %>%
  mutate(fnames = basename(fnames)) %>%
  separate(fnames, into = c("Hz", "Depth", "Distance", "Density"), sep = "_") %>%
  mutate(Depth = as.numeric(Depth), Hz = as.numeric(Hz))
