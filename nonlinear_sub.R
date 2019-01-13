#' ---
#' title: "ヨレモク群落の流速に対する生産量の変動に関する非線形モデル"
#' author: "Inoue yukio"
#' date: "2018 June 16"
#' output: html_document
#' ---
library(tidyverse)
library(lubridate)
library(ggplot2)
library(gridExtra)
library(gtable)
library(dplyr)
library(scales)
library(grid)
library(scales)
library(stringr)
library(xtable)
library(marelac) # for estimating Saturation concentration of Oxygen in water  
library(rstan)

#データの読み込み-------------------------------------------------------------------
fnames = dir("../回流水槽光合成実験_まとめ/Modified_Data/", pattern = "Photosynthesis_WW_k.csv", full = T)
Photosynthesis_WW_k = read_csv(fnames)
Photosynthesis_WW_k %>% print(n = Inf)

######################################################################################################################################
all_data = all_data %>% filter(!Hz == 0) %>% mutate(Photosynthesis_net_k_WW_alph = Photosynthesis_net_k_WW * 10000)

#nls関数による推定
# Model1
model1 = function(theta, speed){
  theta[1] + theta[2] * (1-exp(-theta[3] / theta[2] * speed)) * exp(-theta[4] / theta[2] * speed)
} 

m1out = nls(Photosynthesis_net_k_WW_alph ~ model1(theta, upstream_u_mean), data = all_data, 
            start = list(theta = c(10, 5, -9, 12)), trace = T)

summary(m1out)

theta1 = coef(m1out)

x = seq(min(all_data$upstream_u_mean),max(all_data$upstream_u_mean), length = 50)

plot(x = all_data$upstream_u_mean, y =all_data$Photosynthesis_net_k_WW_alph, xlab = "Velocity", ylab = "Photosynthesis_net")
lines(x = x, y = model1(theta1, x), col = "blue")


#Model2
model2 = function(theta, speed){
  theta[1]*speed / (theta[2] + speed) + theta[3]
} 

m1out = nls(Photosynthesis_net_k_WW_alph ~ model2(theta, upstream_u_mean), data = all_data, 
            start = list(theta = c(1,1,1)), trace = T)

summary(m1out)

theta1 = coef(m1out)

x = seq(min(all_data$upstream_u_mean),max(all_data$upstream_u_mean), length = 50)

plot(x = all_data$upstream_u_mean, y =all_data$Photosynthesis_net_k_WW_alph, xlab = "Velocity", ylab = "Photosynthesis_net")
lines(x = x, y = model2(theta1, x), col = "blue")

#Model3
model3 = function(theta, speed, Tr){
  theta[1]*speed / (theta[2] + speed) + theta[3]*speed + theta[4]
} 

m1out = nls(Photosynthesis_net_k_WW_alph ~ model3(theta, upstream_u_mean), data = all_data, 
            start = list(theta = c(1,1,1,1)), trace = T)

summary(m1out)

theta1 = coef(m1out)

x = seq(min(all_data$upstream_u_mean),max(all_data$upstream_u_mean), length = 50)

plot(x = all_data$upstream_u_mean, y =all_data$Photosynthesis_net_k_WW_alph, xlab = "Velocity", ylab = "Photosynthesis_net")
lines(x = x, y = model3(theta1, x), col = "blue")

######################################################################################################################################

#Stanに渡す引数を準備------------------------------------------------------------------------
M = length(unique(df3$Date))
#以下引数
St1 = filter(df3, Station == "St.C")
St2 = filter(df3, Station == "St.R")
St3 = filter(df3, Station == "St.S")

S1_height = (St1 %>% select(Height))[,1]
S2_height = (St2 %>% select(Height))[,1]
S3_height = (St3 %>% select(Height))[,1]

K1 = length(S1_height) #観測した月の数で割る
K2 = length(S2_height) 
K3 = length(S3_height)

Date = unique(df3$Date)/max(unique(df3$Date))
Date1 = rep(Date, each = K1/M)
Date2 = rep(Date, each = K2/M)
Date3 = rep(Date, each = K3/M)

Temp = unique(Temp$Temperature)                                                     #arrheniusmodelで使用
Temp1 = rep(Temp, each = K1/M)
Temp2 = rep(Temp, each = K2/M)
Temp3 = rep(Temp, each = K3/M)

zizen_theta1_mu =50
zizen_theta2_mu =0.5 　　　　　　　　　　　　　　　　　　　
zizen_theta3_mu =10
zizen_theta4_mu =10
zizen_theta1_tau =10^2
zizen_theta2_tau =10　　　　　　　　　　　　　　　　　　　 
zizen_theta3_tau =10^2
zizen_theta4_tau =10^2

nmodel = 50
modelDate = seq(0, 1, length = nmodel)

a = data.frame(Temp = Temp, Date = Date)                                            #arrheniusmodelで使用
modelTemp = data.frame(approx(a$Date, a$Temp, method = "linear", n = nmodel))[,2]
S = 2.5 #modelerror:対数正規分布の分布のSgima

#Stanmodel実行----------------------------------------------------------------------------------------
stanmodel = stan_model(file = "nonlinearmodel_log_theta4plus.stan")

rng_seed　= 1111
nchains = 4
iter = 3000
ncores = 4 
stanout = sampling(stanmodel, 
                    chains = nchains, 
                    iter = iter, 
                    cores = ncores,
                    seed = rng_seed)

# stanout2 = sampling(stanmodel, 
#                    chains = nchains, cores = nchains,
#                    iter = iter, 
#                    seed = rng_seed)
# 
# p1 = stan_hist(stanout2, pars = c("S1_theta1_mu","S2_theta1_mu","S3_theta1_mu"), ncol =1)
# p1 + xlim(0, 100)

stanout@model_pars
pars =  c( "S1_theta1_mu","S1_theta3_mu","S1_theta4_mu","S2_theta1_mu","S2_theta3_mu","S2_theta4_mu","S3_theta1_mu","S3_theta3_mu","S3_theta4_mu","kotei_theta2_mu")
pars1 =  c("S1_yhat","S1_ypred", "S1_yhatD")
pars2 =  c("S2_yhat","S2_ypred", "S2_yhatD")
pars3 =  c("S3_yhat","S3_ypred", "S3_yhatD")
# pairs(stanout, pars = pars1)
traceplot(stanout, pars = c(pars))
print(stanout, prob=c(0.025, 0.5, 0.975), par = pars)

#期待値yhat
S1_yhat = summary(stanout, pars = pars1[1], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S2_yhat = summary(stanout, pars = pars2[1], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S3_yhat = summary(stanout, pars = pars3[1], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S1_yhat = as.data.frame(S1_yhat)
S2_yhat = as.data.frame(S2_yhat)
S3_yhat = as.data.frame(S3_yhat)
yhat = S1_yhat %>% bind_rows(S2_yhat, S3_yhat) %>% 
  mutate(Station = rep(c("St.C", "St.R", "St.S"), each = nmodel),
         Date = rep(seq(min(df3$Date), max(df3$Date), length = nmodel), 3),
         Date2 = rep(seq(min(df3$Date2), max(df3$Date2), length = nmodel), 3))


#予測値ypred
S1_ypred = summary(stanout, pars = pars1[2], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S2_ypred = summary(stanout, pars = pars2[2], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S3_ypred = summary(stanout, pars = pars3[2], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S1_ypred = as.data.frame(S1_ypred)
S2_ypred = as.data.frame(S2_ypred)
S3_ypred = as.data.frame(S3_ypred)
ypred = S1_ypred %>% bind_rows(S2_ypred, S3_ypred) %>% 
  mutate(Station = rep(c("St.C", "St.R", "St.S"), each = nmodel),
         Date = rep(seq(min(df3$Date), max(df3$Date), length = nmodel), 3),
         Date2 = rep(seq(min(df3$Date2), max(df3$Date2), length = nmodel), 3))

#微分yhatD
S1_yhatD = summary(stanout, pars = pars1[3], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S2_yhatD = summary(stanout, pars = pars2[3], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S3_yhatD = summary(stanout, pars = pars3[3], probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
S1_yhatD = as.data.frame(S1_yhatD)
S2_yhatD = as.data.frame(S2_yhatD)
S3_yhatD = as.data.frame(S3_yhatD)
yhatD = S1_yhatD %>% bind_rows(S2_yhatD, S3_yhatD) %>% 
  mutate(Station = rep(c("St.C", "St.R", "St.S"), each = nmodel),
         Date = rep(seq(min(df3$Date), max(df3$Date), length = nmodel), 3),
         Date2 = rep(seq(min(df3$Date2), max(df3$Date2), length = nmodel), 3),
         mean = mean / max(df3$Date),
         `2.5%` = `2.5%` / max(df3$Date),
         `50%` = `50%` / max(df3$Date),
         `97.5%` = `97.5%` / max(df3$Date)) #単位がcm/(day/max(day))からcm/dayに変換

str(yhat)

#作図yhat,ypred,yhatD
FONTSIZE = 10
xscale = scale_x_datetime(
  breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1",
                          "2015-11-1","2015-12-1","2016-1-1","2016-2-1","2016-3-1","2016-4-1","2016-5-1","2016-6-1","2016-7-1","2016-8-1","2016-9-1","2016-10-1"))),
  labels = c("15-4","5","6","7","8","9","10",
             "11","12","16-1","2","3","4","5","6", "7", "8", "9","10"),
  limits = c(as.POSIXct("2015-4-1"), NA))

fg_yhat = ggplot() +
  geom_line(aes(x = Date2, y = Height, group = No), data = df3, alpha = 0.5) +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = "95%"),data = yhat, alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = "期待値"), data = yhat, size = 2) +
  geom_line(aes(x = Date2, y = `50%`, color = "中央値"), data = yhat, size = 2) +
  labs(x = "", y = "ノコギリモクの生長 (yhat) [cm]", fill = "95％信用区間", color = "") +
  xscale+
  theme_bw(FONTSIZE)+
  theme(
    legend.position = c(0.1,0.75), 
    legend.background = element_blank(), 
    legend.key = element_rect(fill=NA, color=NA), 
    legend.key.size=unit(2, "lines")) + 
  facet_grid(. ~ Station)
fg_yhat

fg_ypred = ggplot() +
  geom_line(aes(x = Date2, y = Height, group = No), data = df3, alpha = 0.25) +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = "95%"),data = ypred, alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = "log分布の平均値"), data = ypred, size = 2) +
  geom_line(aes(x = Date2, y = `50%`, color = "log分布の中央値"), data = ypred, size = 2) +
  labs(x = "", y = "ノコギリモクの生長 (ypred) [cm]", fill = "95％予測区間", color = "") +
  xscale+
  lims(y = c(0,150))+
  theme_bw(FONTSIZE)+
  theme(
    legend.position = c(0.1,0.75), 
    legend.background = element_blank(), 
    legend.key = element_rect(fill=NA, color=NA), 
    legend.key.size=unit(2, "lines")) + 
  facet_grid(. ~ Station)
fg_ypred

fg_yhatD = ggplot() +
  # geom_line(aes(x = Date, y = Height, group = No), data = df3, alpha = 0.25) +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = Station),data = yhatD, alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = Station), data = yhatD, size = 2) +
  # geom_line(aes(x = Date, y = `50%`, color = Station), data = yhatD, size = 2) +
  labs(x = "月", y = "生長速度 (dyhat/dx) [cm / day]", fill = "95％信用区間", color = "") +
  xscale+
  lims(y = c(-0.1, 0.2))+
  theme_bw(FONTSIZE)+
  theme(
    legend.position = c(1.3,0.5), 
    legend.background = element_blank(), 
    legend.key = element_rect(fill=NA, color=NA), 
    legend.key.size=unit(2, "lines")) #+
# facet_grid(. ~ Station)
fg_yhatD
# パラメーターの比較----------------------------------------------------------------------------------------------------
#ヒストグラムで各パラメーターの確認
stan_hist(stanout, par=pars, ncol=3, bins = 30)

#パラメーターの抽出
PM = extract(stanout, par=pars)
str(PM) #listの形式

#theta1比較（最大値）
par(mfrow = c(2,2))
hist(PM$S2_theta1_mu - PM$S1_theta1_mu, col="skyblue", border=NA, breaks = 20, main = "最大値 (theta1) St.R-St.C")  
abline(v=0, lty = 2)

hist(PM$S3_theta1_mu - PM$S1_theta1_mu, col="skyblue", border=NA, breaks = 20, main = "最大値 (theta1) St.S-St.C")  
abline(v=0, lty = 2)

# theta3比較(生長速度？)
hist(PM$S2_theta3_mu - PM$S1_theta3_mu, col="skyblue", border=NA, breaks = 20, main = "生長速度 (theta3) St.R-St.C")  
abline(v=0, lty = 2)

hist(PM$S3_theta3_mu - PM$S1_theta3_mu, col="skyblue", border=NA, breaks = 20, main = "生長速度 (theta3) St.S-St.C")  
abline(v=0, lty = 2)

#確率で比較
theta1R_C = PM$S2_theta1_mu - PM$S1_theta1_mu
theta1S_C = PM$S3_theta1_mu - PM$S1_theta1_mu
theta3R_C = PM$S2_theta3_mu - PM$S1_theta3_mu
theta3S_C = PM$S3_theta3_mu - PM$S1_theta3_mu
theta4R_C = PM$S2_theta4_mu - PM$S1_theta4_mu
theta4S_C = PM$S3_theta4_mu - PM$S1_theta4_mu

a = c("theta1R_C",
      sum(theta1R_C>0)/length(theta1R_C),
      mean(PM$S2_theta1_mu / PM$S1_theta1_mu),
      sd(PM$S2_theta1_mu / PM$S1_theta1_mu))

b = c("theta1S_C",
      sum(theta1S_C>0)/length(theta1S_C),
      mean(PM$S3_theta1_mu / PM$S1_theta1_mu),
      sd(PM$S3_theta1_mu / PM$S1_theta1_mu))

c = c("theta3R_C",
      sum(theta3R_C>0)/length(theta3R_C),
      mean(PM$S2_theta3_mu / PM$S1_theta3_mu),
      sd(PM$S2_theta3_mu / PM$S1_theta3_mu))

d = c("theta3S_C",
      sum(theta3S_C>0)/length(theta3S_C),
      mean(PM$S3_theta3_mu / PM$S1_theta3_mu),
      sd(PM$S3_theta3_mu / PM$S1_theta3_mu))

e = c("theta4R_C",
      sum(theta4R_C>0)/length(theta4R_C),
      mean(PM$S2_theta4_mu / PM$S1_theta4_mu),
      sd(PM$S2_theta4_mu / PM$S1_theta4_mu))

f = c("theta4S_C",
      sum(theta4S_C>0)/length(theta4S_C),
      mean(PM$S3_theta4_mu / PM$S1_theta4_mu),
      sd(PM$S3_theta4_mu / PM$S1_theta4_mu))


PM_hikaku = data.frame(matrix(c(a,b,c,d,e,f),ncol = 4, byrow = T))
colnames(PM_hikaku) =c("type", "確率", "Cのx倍の平均", "Cのx倍の標準偏差")
PM_hikaku                                   

#作図 パラメーターのヒストグラム
C=as.data.frame(PM$S1_theta1_mu)
C= C %>% mutate(Station="St.C") %>% select(Station, "theta1"=`PM$S1_theta1_mu`)

R=as.data.frame(PM$S2_theta1_mu)
R= R %>% mutate(Station="St.R") %>% select(Station, "theta1"=`PM$S2_theta1_mu`)

S=as.data.frame(PM$S3_theta1_mu)
S= S %>% mutate(Station="St.S") %>% select(Station, "theta1"=`PM$S3_theta1_mu`)

theta1_all= C %>% bind_rows(R, S) 
theta1_all=as.data.frame(theta1_all)

ylabel="頻度"
xlabel=expression("パラメーター"~theta[1])

theta1_hist = ggplot()+
  geom_histogram(aes(theta1), theta1_all, bins = 30)+
  labs(y=ylabel, x=xlabel)+
  facet_grid(Station~.)+
  theme_bw(FONTSIZE)+
  labs(title ="モデル2のパラメーターの分布(成熟個体データ使用)")+
  theme(
    legend.position = c(0.1,0.75), 
    legend.background = element_blank(), 
    legend.key = element_rect(fill=NA, color=NA), 
    legend.key.size=unit(2, "lines"))  
theta1_hist                                      #PM_hikakuでは100％Rが高いが図で見るとかぶってる!?


C=as.data.frame(PM$S1_theta3_mu)
C= C %>% mutate(Station="St.C") %>% select(Station, "theta3"=`PM$S1_theta3_mu`)

R=as.data.frame(PM$S2_theta3_mu)
R= R %>% mutate(Station="St.R") %>% select(Station, "theta3"=`PM$S2_theta3_mu`)

S=as.data.frame(PM$S3_theta3_mu)
S= S %>% mutate(Station="St.S") %>% select(Station, "theta3"=`PM$S3_theta3_mu`)

theta3_all= C %>% bind_rows(R, S) 
theta3_all=as.data.frame(theta3_all)

ylabel=""
xlabel=expression("パラメーター"~theta[3])

theta3_hist = ggplot()+
  geom_histogram(aes(theta3), theta3_all, bins = 30)+
  labs(y=ylabel, x=xlabel)+
  facet_grid(Station~.)+
  theme_bw(FONTSIZE)+
  theme(
    legend.position = c(0.1,0.75), 
    legend.background = element_blank(), 
    legend.key = element_rect(fill=NA, color=NA), 
    legend.key.size=unit(2, "lines"))  
theta3_hist                          #PM_hikakuでは100％Rが高いが図で見るとかぶってる!?

C=as.data.frame(PM$S1_theta4_mu)
C= C %>% mutate(Station="St.C") %>% select(Station, "theta4"=`PM$S1_theta4_mu`)

R=as.data.frame(PM$S2_theta4_mu)
R= R %>% mutate(Station="St.R") %>% select(Station, "theta4"=`PM$S2_theta4_mu`)

S=as.data.frame(PM$S3_theta4_mu)
S= S %>% mutate(Station="St.S") %>% select(Station, "theta4"=`PM$S3_theta4_mu`)

theta4_all= C %>% bind_rows(R, S) 
theta4_all=as.data.frame(theta4_all)

ylabel="頻度"
xlabel=expression("パラメーター"~theta*4)

theta4_hist = ggplot()+
  geom_histogram(aes(theta4), theta4_all, bins = 30)+
  labs(y=ylabel, x=xlabel)+
  facet_grid(Station~.)+
  theme_bw(FONTSIZE)+
  theme(
    legend.position = c(0.1,0.75), 
    legend.background = element_blank(), 
    legend.key = element_rect(fill=NA, color=NA), 
    legend.key.size=unit(2, "lines"))  
theta4_hist                                      #PM_hikakuでは100％Rが高いが図で見るとかぶってる!?

###############################################################################################################################################
#得られた曲線を微分する(手動)
a1 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][1]
b1 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][10]
c1 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][2]
d1 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][3]

a2 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][4]
b2 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][10]
c2 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][5]
d2 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][6]

a3 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][7]
b3 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][10]
c3 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][8]
d3 = summary(stanout, pars = pars, probs = c(0.025, 0.5, 0.975))$summary[, c("mean")][9]

f = expression((theta1/ (1 + (theta1/theta2-1) * exp(-theta3* x))) * exp(-theta4*x)) #用いるもとのモデル
D = D(f,"x") #微分式の確認

tmp = deriv(f, c("x"), function(theta1,theta2,theta3,theta4,x){} )
y1 = tmp(a1,b1,c1,d1,mx)
y2 = tmp(a2,b2,c2,d2,mx)
y3 = tmp(a3,b3,c3,d3,mx)

dy1 = attributes(y1)$gradient
dy2 = attributes(y2)$gradient
dy3 = attributes(y3)$gradient

plot(x = df3$Date, y =df3$Height, xlab = "Date", ylab = "Height(cm)", ylim = c(-30,70))
lines(x = x, y = y1)
lines(x = x, y = y2, col = "red")
lines(x = x, y = y3, col = "blue")

lines(x = x, y = dy1)
lines(x = x, y = dy2, col = "red") 
lines(x = x, y = dy3, col = "blue")

# 図の集約################################################################################################################################################
fg1  = ggplotGrob(fg_yhat)
fg2  = ggplotGrob(fg_ypred)
fg3  = ggplotGrob(fg_yhatD)

maxwidth = unit.pmax(fg1[["widths"]],
                     fg2[["widths"]])

fg1[["widths"]]=maxwidth
fg2[["widths"]]=maxwidth
fg3[["widths"]]=maxwidth

fg_growth = rbind(fg1,fg2,fg3, size="first")


width = 210/25.4 #mmをinchに変換、A4の大きさに設定
height = 297/25.4

cairo_pdf(file = "fg_growth_model2_kisyo.pdf",
          width = width,
          height  = height,
          family = "Times New Roman",
          fallback_resolution = 1200)

grid.draw(fg_growth)
dev.off()

fg4  = ggplotGrob(theta1_hist)
fg5  = ggplotGrob(theta3_hist)
fg6  = ggplotGrob(theta4_hist)

maxwidth = unit.pmax(fg4[["widths"]],
                     fg5[["widths"]],
                     fg6[["widths"]])

fg4[["widths"]]=maxwidth
fg5[["widths"]]=maxwidth
fg6[["widths"]]=maxwidth

fg_theta = cbind(fg4,fg5,fg6, size="first")


width = 297/25.4 #mmをinchに変換、A4の大きさに設定
height = 210/25.4

cairo_pdf(file = "fg_theta_model2_kisyo.pdf",
          width = width,
          height  = height,
          family = "Times New Roman",
          fallback_resolution = 1200)

grid.draw(fg_theta)
dev.off()


#161120報告会用作図yhat,ypred,yhatD------------------------------------------------------------------------
# FONTSIZE = 18
# xscale = scale_x_datetime(
#       breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1","2015-11-1","2015-12-1",
#                               "2016-1-1","2016-2-1","2016-3-1", "2016-4-1", "2016-5-1", "2016-6-1", "2016-7-1", "2016-8-1", "2016-9-1", "2016-10-1"))),
#       labels = c("4","","6","","8","","10","","12","","2","", "4", "", "6", "", "8", "", "10"),
#       limits = c(as.POSIXct("2015-4-1"), as.POSIXct("2016-10-1")))
# 
# xscale2 =  scale_x_datetime(
#   breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1","2015-11-1","2015-12-1",
#                           "2016-1-1","2016-2-1","2016-3-1", "2016-4-1", "2016-5-1", "2016-6-1", "2016-7-1", "2016-8-1", "2016-9-1", "2016-10-1"))),
#   labels = c("4","5","6","7","8","9","10","11","12","1","2","3", "4", "5", "6", "7", "8", "9", "10"),
#   limits = c(as.POSIXct("2015-4-1"), as.POSIXct("2016-10-1")))
# 
# ylabel = expression("Total length [cm]")
# # fg_yhat = ggplot() +
# #   geom_line(aes(x = Date2, y = Height, group = No), data = df3, alpha = 0.5) +
# #   geom_line(aes(x = Date2, y = mean, color = "期待値"), data = yhat, size = 2) +
# #   geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = "95%"),data = yhat, alpha = 0.15) +
# #   # geom_line(aes(x = Date2, y = `50%`, color = "中央値"), data = yhat, size = 2) +
# #   labs(x = "月", y = ylabel, fill = "", color = "") +
# #   xscale+
# #   scale_fill_discrete(labels = c("95%信用区間"))+
# #   theme_gray(FONTSIZE)+
# #   theme(
# #     legend.position = c(0.1,0.75),
# #     legend.background = element_blank(),
# #     legend.key = element_rect(fill=NA, color=NA),
# #     legend.box.just = "left")+
# #     # legend.key.size=unit(2, "lines")) +
# #   facet_grid(. ~ Station)
# # fg_yhat
# 
# df3 = df3 %>% filter(Station %in% c("St.C", "St.R"))
# df3$Station = factor(df3$Station,
#                        levels = c("St.R","St.C"), labels = c("リン添加区", "対照区"))
# 
# yhat = yhat %>% filter(Station %in% c("St.C", "St.R"))
# yhat$Station = factor(yhat$Station,
#                       levels = c("St.R","St.C"), labels = c("リン添加区", "対照区"))
# 
# ypred = ypred %>% filter(Station %in% c("St.C", "St.R"))
# ypred$Station = factor(ypred$Station,
#                        levels = c("St.R","St.C"), labels = c("リン添加区", "対照区"))
# 
# 
# fg_ypred = ggplot() +
#   geom_line(aes(x = Date2, y = Height, group = No), data = df3, alpha = 0.5) +
#   geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = "95%"),data = ypred, alpha = 0.15) +
#   # geom_line(aes(x = Date2, y = mean, color = "log分布の平均値"), data = ypred, size = 2) +
#   # geom_line(aes(x = Date2, y = `50%`, color = "log分布の中央値"), data = ypred, size = 2) +
#   geom_line(aes(x = Date2, y = mean, color = "期待値"), data = yhat, size = 2) +
#   labs(x = "Month", y = ylabel, fill = "", color = "") +
#     xscale+
#     scale_y_continuous(breaks = c(seq(0, 100, by = 25)), limits = c(0,100))+
#     scale_fill_discrete(labels = c("95%予測区間"))+
#   theme_classic(FONTSIZE)+
#   theme(
#     panel.border = element_rect(color = "black", fill = NA),
#       legend.position = c(0.1,0.75),
#       legend.background = element_blank(),
#       legend.key = element_rect(fill=NA, color=NA),
#       legend.box.just = "left")+
#     facet_grid(. ~ Station)
# fg_ypred
# 
# yhatD = yhatD[yhatD$Station %in% c("St.C", "St.R"),]
# yhatD$Station = factor(yhatD$Station,
#                        levels = c("St.R","St.C"), labels = c("リン添加区", "対照区"))
# 
# 
# ylabel = expression("Growth rate"~"["~cm~day^-1~"]")
# fg_yhatD = ggplot() +
#   geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = Station),data = yhatD, alpha = 0.15) +
#   geom_line(aes(x = Date2, y = mean, color = Station), data = yhatD, size = 2) +
#   labs(x = "Month", y = ylabel) +
#   xscale2+
#   lims(y = c(-0.12, 0.2))+
#   theme_classic(FONTSIZE)+
#   theme(
#     panel.border = element_rect(color = "black", fill = NA),
#     legend.position = c(0.23,0.9),
#     legend.background = element_blank(),
#     legend.key = element_rect(fill=NA, color=NA),
#     legend.title = element_blank())
# fg_yhatD
# 
# 
# library(gridExtra)
# library(gtable)
# library(grid)
# 
# gp1 = ggplotGrob(fg_ypred)
# a = grep("guide-box", gp1$layout$name)
# b = gp1[["grobs"]][[a]]
# c1 = b[["grobs"]][[1]]
# c2 = b[["grobs"]][[2]]
# 
# gtable_show_layout(b)
# grid.draw(b)
# 
# gtable_show_layout(c1)
# grid.draw(c1)
# 
# gtable_show_layout(c2)
# grid.draw(c2)
# 
# c1$heights[1:3] = unit(0,'mm') #凡例の上のスペースを0にする
# c2$heights[1:3] = unit(0,'mm')
# b[["grobs"]][[1]] = c1
# b[["grobs"]][[2]] = c2
# b$heights[2] = unit(8,'mm')
# b$heights[4] = unit(8,'mm')
# gp1[["grobs"]][[a]] = b
# 
# 
# width = 240
# height= 140
# 
# png(file = "fgypred.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp1)
# dev.off()
# 
# width = 130
# height= 130
# 
# png(file = "fgyhatD.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(fg_yhatD)
# dev.off()

#161227修論用作図yhat,ypred,yhatD------------------------------------------------------------------------
FONTSIZE = 15

df3$Station = factor(df3$Station,
                     levels = c("St.C", "St.R", "St.S"), labels =  c("St.C", "St.R", "St.S"))

yhat$Station = factor(yhat$Station,
                      levels = c("St.R", "St.S", "St.C"), labels = c("St.R", "St.S", "St.C"))

ypred$Station = factor(ypred$Station,
                       levels = c("St.C", "St.R", "St.S"), labels =  c("St.C", "St.R", "St.S"))

yhatD$Station = factor(yhatD$Station,
                       levels = c("St.R", "St.S", "St.C"), labels = c("St. R", "St. S", "St. C"))

ylabel = expression("ノコギリモクの全長 [cm]")
xlabel = "月"
fg_ypred_all = ggplot() +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = Station),data = yhat, alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = Station), data = yhat, size = 2) +
  labs(x = xlabel, y = ylabel, fill = "", color = "") +
  scale_x_datetime(
    breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1","2015-11-1","2015-12-1",
                            "2016-1-1","2016-2-1","2016-3-1", "2016-4-1", "2016-5-1", "2016-6-1", "2016-7-1", "2016-8-1", "2016-9-1", "2016-10-1"))),
    labels = c("4\nH 27","5","6","7","8","9","10","11","12","1\nH 28","2","3", "4", "5", "6", "7", "8", "9", "10"),
    limits = c(as.POSIXct("2015-4-1"), as.POSIXct("2016-10-1"))) +
  scale_y_continuous(breaks = c(seq(0, 50, by = 10)), limits = c(0,50))+
  scale_color_manual(breaks=c("St. R", "St. S", "St. C"), labels = c("リントル試験区", "点滴灌水チューブ試験区", "筏内対照区"), values = c("#F8766D", "#7CAE00", "#00BFB4"))+
  scale_fill_manual(breaks=c("St. R", "St. S", "St. C"), labels = c("リントル試験区", "点滴灌水チューブ試験区", "筏内対照区"), values = c("#F8766D", "#7CAE00", "#00BFB4"))+
  theme_classic(FONTSIZE)+
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title.y =element_text(size = FONTSIZE - 6),
    axis.title.x =element_text(size = FONTSIZE - 5),
    legend.position = "right",
    legend.background = element_blank(),
    legend.justification =  c(0,0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = FONTSIZE - 5),
    legend.key = element_rect(fill=NA, color=NA),
    legend.key.size=unit(0.5, "lines"))


ylabel = expression("ノコギリモクの全長 [cm]")
xlabel = "月"
fg_ypred_C=ggplot() +
  geom_line(aes(x = Date2, y = Height, group = No), data = df3[df3$Station == "St.C",], alpha = 0.5) +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = "95%"),data = yhat[yhat$Station == "St.C",], fill = "#00BFB4", alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = "期待値"), data = yhat[yhat$Station == "St.C",], color = "#00BFB4", size = 2) +
  labs(x = xlabel, y = ylabel, fill = "", color = "") +
  scale_x_datetime(
    breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1","2015-11-1","2015-12-1",
                            "2016-1-1","2016-2-1","2016-3-1", "2016-4-1", "2016-5-1", "2016-6-1", "2016-7-1", "2016-8-1", "2016-9-1", "2016-10-1"))),
    labels = c("4\nH 27","5","6","7","8","9","10","11","12","1\nH 28","2","3", "4", "5", "6", "7", "8", "9", "10"),
    limits = c(as.POSIXct("2015-4-1"), as.POSIXct("2016-10-1"))) +
  scale_y_continuous(breaks = c(seq(0, 100, by = 25)), limits = c(0,100))+
  scale_fill_discrete(labels = c("95%予測区間"))+
  annotate("text", x = as.POSIXct("2015-4-1"), y = Inf, label = "筏内対照区", vjust = 1.8, hjust =0)+
  theme_classic(FONTSIZE)+
  theme(
    axis.title.y =element_blank(),
    axis.title.x =element_text(size = FONTSIZE - 5),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position ="none")


# ylabel = expression("Total Length [cm]")


fg_ypred_R = ggplot() +
  geom_line(aes(x = Date2, y = Height, group = No), data = df3[df3$Station == "St.R",],alpha = 0.5) +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = "95%"),data = yhat[yhat$Station == "St.R",], fill = "#F8766D", alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = "期待値"), data = yhat[yhat$Station == "St.R",], color = "#F8766D", size = 2) +
  labs(x = xlabel, y = ylabel, fill = "", color = "") +
  scale_x_datetime(
    breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1","2015-11-1","2015-12-1",
                            "2016-1-1","2016-2-1","2016-3-1", "2016-4-1", "2016-5-1", "2016-6-1", "2016-7-1", "2016-8-1", "2016-9-1", "2016-10-1"))),
    labels = c("4\nH 27","5","6","7","8","9","10","11","12","1\nH 28","2","3", "4", "5", "6", "7", "8", "9", "10"),
    limits = c(as.POSIXct("2015-4-1"), as.POSIXct("2016-10-1"))) +
  scale_y_continuous(breaks = c(seq(0, 100, by = 25)), limits = c(0,100))+
  scale_fill_discrete(labels = c("95%予測区間"))+
  annotate("text", x = as.POSIXct("2015-4-1"), y = Inf, label = "リントル試験区", vjust = 1.8, hjust =0)+
  theme_classic(FONTSIZE)+
  theme(
    axis.title.y =element_blank(),
    axis.title.x =element_blank(),
    axis.text.x =element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_blank(),
    legend.position =c(0.25,0.65),
    legend.title = element_blank(),
    legend.text = element_text(size = FONTSIZE - 5),
    legend.key = element_rect(fill=NA, color=NA),
    legend.key.size = unit(0.5, "lines"),
    legend.box.just = "left")


fg_ypred_S = ggplot() +
  geom_line(aes(x = Date2, y = Height, group = No), data = df3[df3$Station == "St.S",], alpha = 0.5) +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = "95%"),data = yhat[yhat$Station == "St.S",], fill = "#7CAE00", alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = "期待値"), data = yhat[yhat$Station == "St.S",], color = "#7CAE00", size = 2) +
  labs(x = xlabel, y = ylabel, fill = "", color = "") +
  scale_x_datetime(
    breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1","2015-11-1","2015-12-1",
                            "2016-1-1","2016-2-1","2016-3-1", "2016-4-1", "2016-5-1", "2016-6-1", "2016-7-1", "2016-8-1", "2016-9-1", "2016-10-1"))),
    labels = c("4\nH 27","5","6","7","8","9","10","11","12","1\nH 28","2","3", "4", "5", "6", "7", "8", "9", "10"),
    limits = c(as.POSIXct("2015-4-1"), as.POSIXct("2016-10-1"))) +
  scale_y_continuous(breaks = c(seq(0, 100, by = 25)), limits = c(0,100))+
  scale_fill_discrete(labels = c("95%予測区間"))+
  annotate("text", x = as.POSIXct("2015-4-1"), y = Inf, label = "点滴灌水チューブ試験区", vjust = 1.8, hjust =0)+
  theme_classic(FONTSIZE)+
  theme(
    axis.title.y = element_text(size = FONTSIZE - 5),
    axis.title.x =element_blank(),
    axis.text.x =element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position ="none")

fg_ypred_C
fg_ypred_R
fg_ypred_S


ylabel = expression("ノコギリモクの生長速度"~"["~cm~day^-1~"]")
fg_yhatD = ggplot() +
  geom_ribbon(aes(x = Date2, y = mean, ymin = `2.5%`, ymax = `97.5%`, fill = Station),data = yhatD, alpha = 0.15) +
  geom_line(aes(x = Date2, y = mean, color = Station), data = yhatD, size = 2) +
  labs(x = xlabel, y = ylabel) +
  scale_x_datetime(
    breaks = c(as.POSIXct(c("2015-4-1","2015-5-1","2015-6-1","2015-7-1","2015-8-1","2015-9-1","2015-10-1","2015-11-1","2015-12-1",
                            "2016-1-1","2016-2-1","2016-3-1", "2016-4-1", "2016-5-1", "2016-6-1", "2016-7-1", "2016-8-1", "2016-9-1", "2016-10-1"))),
    labels = c("4\nH 27","5","6","7","8","9","10","11","12","1\nH 28","2","3", "4", "5", "6", "7", "8", "9", "10"),
    limits = c(as.POSIXct("2015-4-1"), as.POSIXct("2016-10-1"))) +
  scale_color_manual(breaks=c("St. R", "St. S", "St. C"), labels = c("リントル試験区", "点滴灌水チューブ試験区", "筏内対照区"), values = c("#F8766D", "#7CAE00", "#00BFB4"))+
  scale_fill_manual(breaks=c("St. R", "St. S", "St. C"), labels = c("リントル試験区", "点滴灌水チューブ試験区", "筏内対照区"), values = c("#F8766D", "#7CAE00", "#00BFB4"))+
  lims(y = c(-0.12, 0.2))+
  theme_classic(FONTSIZE)+
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title.y =element_text(size = FONTSIZE - 6),
    axis.title.x =element_text(size = FONTSIZE - 5),
    legend.background = element_blank(),
    legend.justification =  c(0,0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = FONTSIZE - 5),
    legend.key = element_rect(fill=NA, color=NA), 
    legend.key.size=unit(0.5, "lines"))  
fg_yhatD


library(gridExtra)
library(gtable)
library(grid)

gp1 = ggplotGrob(fg_ypred_R)
gp2 = ggplotGrob(fg_ypred_S)
gp3 = ggplotGrob(fg_ypred_C)
gp4 = ggplotGrob(fg_yhatD)
gp5 = ggplotGrob(fg_ypred_all)

maxwidth = unit.pmax(gp1[["widths"]], gp2[["widths"]], gp3[["widths"]], gp4[["widths"]])

gp1[["widths"]] = maxwidth
gp2[["widths"]] = maxwidth
gp3[["widths"]] = maxwidth
gp4[["widths"]] = maxwidth


a = grep("guide-box", gp4$layout$name)
b = gp4[["grobs"]][[a[1]]]

c = b[["grobs"]][[1]]
c$heights[3] = unit(0,'mm') #凡例の上のスペースを0にする
c$heights[4:7] = c$heights[4:5]*1.6　#凡例の行間を1.1倍にする
b[["grobs"]][[1]] = c
gp4[["grobs"]][[a[1]]] = b


gp123 = rbind(gp1, gp2, gp3, size="first")


width = 160
height= 200

png(file = "fgypred_all_kisyo.png",
    type = "cairo",
    width = width,
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
grid.draw(gp123)
dev.off()


width = 160
height= 70

png(file = "fgyhatD_all_kisyo.png",
    type = "cairo",
    width = width,
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
grid.draw(gp4)
dev.off()

width = 160
height= 160

png(file = "fg_ypred_all.png",
    type = "cairo",
    width = width,
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
grid.draw(gp5)
dev.off()


# 数値まとめ-----------------------------------------------------------------------------------------------------------------------------

R = yhat[yhat$Station == "St.R",]
S = yhat[yhat$Station == "St.S",]
C = yhat[yhat$Station == "St.C",]

R %>% filter(mean == max(mean))
S %>% filter(mean == max(mean))
C %>% filter(mean == max(mean))

RD = yhatD[yhat$Station == "St.R",]
SD = yhatD[yhat$Station == "St.S",]
CD = yhatD[yhat$Station == "St.C",]

RD %>% filter(mean == max(mean))
SD %>% filter(mean == max(mean))
CD %>% filter(mean == max(mean))

RD %>% filter(mean == min(mean))
SD %>% filter(mean == min(mean))
CD1 = CD %>% filter(Date >= 300)
CD1 %>% filter(mean == min(mean))
