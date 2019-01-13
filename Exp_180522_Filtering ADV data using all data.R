#' ---
#' title: "Filtering ADV data"
#' author: "Inoue Yukio"
#' date: "2018 May 24"
#' output: html_document
#' ---


# library(tidyr)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(gtable)
library(scales)
library(grid)
library(scales)
library(stringr)
library(xtable)
library(dplyr)

# Read packages ----------------------
library(tidyverse)
library(lubridate)
library(zoo)

# Declare the function --------------------------
gradient = function(x){
  n = length(x)
  j = 2:(n-1)
  c(x[2] - x[1],
    sapply(j, function(j, z){(z[j+1] - z[j-1])*0.5}, z = x), #functionで幅を持たせて変位を計算している
    x[n] - x[n-1])
}

snr_corr_filter = function(advdata, snr.threshold=2, corr.threshold=30) {  
  # A filter for the SNR and CORR.
  require(tidyverse)
  snr = advdata %>% select(dplyr::starts_with("snr")) %>% mutate_all(funs(na=. < snr.threshold))
  corr = advdata %>% select(dplyr::starts_with("corr")) %>% mutate_all(funs(na=. < corr.threshold))
  tmp = as.matrix(advdata %>% select(dplyr::contains("vel")))
  tmp[as.matrix(snr %>% select(dplyr::ends_with("na")))] =NA
  tmp[as.matrix(corr %>% select(dplyr::ends_with("na")))] =NA
  tmp
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
  data.frame(x,y,radius)
}
distance_radius=function(x,y, ab, alpha=0) {
  # Determine the distance to the point and the distance (radius)
  # to the edge of the ellipse
  theta = atan(y/x)
  distance=sqrt(x^2+y^2)
  theta[x<0]=theta[x<0]+pi
  radius = ellipse_radius(ab, theta, alpha=alpha)
  data.frame(theta, distance, radius=radius$radius)
}

make_ellipse = function(ab,alpha=0,x0=0,y0=0,size=500) {
  theta=seq(0,2*pi,length=size)
  x=x0+ab[1]*cos(theta)*cos(alpha)-ab[2]*sin(theta)*sin(alpha)
  y=y0+ab[1]*cos(theta)*sin(alpha)+ab[2]*sin(theta)*cos(alpha)
  data.frame(x,y)
}
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
  f[1] = ifelse(is.na(f[1]), 0, f[1]) #初期値がNAの場合エラーを返すのでNAの場合は便宜上0を代入
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
  f[1] = ifelse(is.na(f[1]), 0, f[1]) #初期値がNAの場合エラーを返すのでNAの場合は便宜上0を代入
  f = zoo::na.locf(f, na.rm=FALSE)
  f = f + f_mean
  number_of_spikes = rep(0, 10)
  while(deltaN > 0 & i < 10) {
    f = f - f_mean
    out=despike(f, good_data, figure=F)
    length(out$f)
    dim(out$mat)
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

#Read all datas---------------------------------------------------------------
#vectorino data
vector_names = dir("../回流水槽光合成実験5/Data/Velocity/", pattern = "dat", full = T)

cnames_vectrino=c("time", "counter", "status",
                  "xvel", "yvel", "zvel", "z2vel",
                  "amp1", "amp2", "amp3", "amp4",
                  "snr1", "snr2", "snr3", "snr4",
                  "corr1", "corr2", "corr3", "corr4")

vector = lapply(vector_names, read_table , col_names=cnames_vectrino)

vector_test1 = vector("list", length(vector_names))
vector_test2 = vector("list", length(vector_names))
vector_test_result = vector("list", length(vector_names))

for(i in 1:length(vector_names)){
  # Try on a small subset of data
  # Trial run on the first 10000 data points
  vector_test1[[i]] = vector[[i]] %>% slice(1000:11000) #testのため一部抽出
  
  # Remove spikes
  # First clean up the data based on the SNR and CORR
  x = names(vector_test1[[i]] %>% select(dplyr::contains("vel")))
  vector_test1[[i]][,x]=snr_corr_filter(vector_test1[[i]])
  
  # Next run nk_despike() to clean it all up.
  vector_test2[[i]] = vector_test1[[i]] %>% mutate(tau = 1:length(xvel),
                                                   u = zoo::na.locf(xvel,na.rm=F), #NAを直近のNA以外の数値に置き換える
                                                   v = zoo::na.locf(yvel,na.rm=F),
                                                   w = zoo::na.locf(zvel,na.rm=F))
  
  vector_test_result[[i]] = vector_test2[[i]] %>% mutate_at(vars(u,v,w), nk_despike)
}

basename(vector_names)
names(vector_test_result) = paste("", str_extract(vector_names, "[0-9]{1,2}_[0-9]{1,2}_[0-9]{1,3}_[HLM]"), sep = "")
# write.csv(vector_test_result, "../回流水槽光合成実験5/vector_test_result.csv")

vector_test_result2 = vector_test_result %>%
  bind_rows(.id = "location") %>%
  mutate(location2 = location) %>%
  tidyr::separate(location2, into=c("Hz", "Depth", "Distance", "Density"), sep = '_') %>%
  select(Hz, Depth,
         Distance, Density,
         time,
         xvel, yvel, zvel, z2vel,
         tau, u, v, w) %>%
  mutate(Depth = as.numeric(Depth),
         Hz = as.numeric(Hz)) %>% 
  mutate(Distance = as.numeric(Distance)) %>% 
  arrange(Hz, Distance, Depth, Density) %>%
  tbl_df()

vector_all = vector_test_result2 %>%
  mutate(Hz = ifelse(Hz == 5, 0.5, Hz)) 

########
# str(vector_all)
# range(na.omit(vector_all[vector_all$Hz == 4 &
#                  vector_all$Depth == 6 &
#                  vector_all$Distance == 41 &
#                  vector_all$Density == "H" ,]$time))
# 
# 
# a =  na.omit(vector_all[vector_all$Hz == 4 &
#                           vector_all$Depth == 6 &
#                           vector_all$Distance == 41 &
#                           vector_all$Density == "H" ,]$time)
# 
# range(a[!is.infinite(a)])
# 
# 
# sum(is.na(vector_all))
# nrow(vector_all)*ncol(vector_all)
# nrow(vector_all)
# 
# na_count = vector_all %>% group_by(Hz, Depth, Distance, Density) %>% 
#   summarise_at(vars(xvel, yvel, zvel, u, v, w), funs(sum(is.na(.)), length(.)))
# tail(data.frame(na_count))
# nrow(data.frame(na_count))
# unique(vector_all$type)

###################
vector_mean=vector_all %>%
  group_by(Hz, Depth, Distance, Density) %>%
  summarise_at(vars(xvel, yvel, zvel, u, v, w),funs(mean, sd)) %>% 
  data.frame()

vector_zansa = vector_all %>%
  left_join(vector_mean, by=c("Depth", "Distance", "Hz", "Density")) %>%
  mutate(zansaX = u - u_mean, zansaY = v - v_mean, zansaZ = w - w_mean, stress = zansaX * zansaZ)

vector_zansa=data.frame(vector_zansa)

vector_zansa_mean = vector_zansa %>% group_by(Depth, Distance, Hz, Density) %>%
  summarise(zansaX_mean = mean(zansaX^2), zansaY_mean = mean(zansaY^2), zansaZ_mean = mean(zansaZ^2),
            stress_mean = mean(stress),TKE = 0.5*(zansaX_mean + zansaY_mean + zansaZ_mean))
vector_zansa_mean=data.frame(vector_zansa_mean)

upstream_u_mean = vector_mean %>%
  filter(Distance == 26) %>%
  select(Hz, Depth, Density, upstream_u_mean = u_mean) %>% 
  group_by(Hz, Density) %>% 
  summarise(upstream_u_mean = mean(upstream_u_mean))

vector_mean = vector_mean %>% left_join(upstream_u_mean, by = c("Hz", "Density"))
vector_zansa_mean = vector_zansa_mean %>% left_join(upstream_u_mean, by = c("Hz", "Density"))

vector_mean2 = vector_mean %>%
  mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
         Depth2 = recode(Depth2, "16" = "above", "6" = "in")) 

vector_zansa_mean2 = vector_zansa_mean %>%
  mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
         Depth2 = recode(Depth2, "16" = "above", "6" = "in")) 

write.csv(vector_all, "../回流水槽光合成実験5/Modified_Data/vector_all.csv")
write.csv(vector_zansa, "../回流水槽光合成実験5/Modified_Data/vector_zansa.csv")
write.csv(vector_zansa_mean, "../回流水槽光合成実験5/Modified_Data/vector_zansa_mean.csv")
write.csv(vector_mean, "../回流水槽光合成実験5/Modified_Data/vector_mean.csv")
write.csv(vector_zansa_mean2, "../回流水槽光合成実験5/Modified_Data/vector_zansa_mean2.csv")
write.csv(vector_mean2, "../回流水槽光合成実験5/Modified_Data/vector_mean2.csv")

# Plot of the Vectrino data------------------------------------------------
p1 = ggplot(vector_all) +
  # geom_line(aes(tau, xvel)) +
  geom_line(aes(tau, u))+
  facet_grid(Density + Depth ~ Distance + Hz)
p1

