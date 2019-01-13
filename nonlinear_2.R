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
fname = dir("../回流水槽光合成実験5/Modified_Data/", pattern = "all_data", full = T)
all_data = read_csv(fname)[-1] %>% 
  mutate(Density = as.factor(Density),
         Density = factor(Density, levels = c("High", "Middle", "Low")))

#Stan----------
Gross_data = all_data %>%
  filter(Hz != 0, type == "GP")

N = nrow(Gross_data)                                      #データの数
speed = Gross_data$upstream_u_mean*100                        #観測した速度
gross = Gross_data$Slope                                  #観測した総一次生産量 

#密度条件ごとに変化しないパラメーター
M = 50                                                  #期待値に使う結果の数
prior_log_sigma = as.vector(c(1,1,1,1))               #各パラメータの事前分散
prior_log_mu = as.vector(c(1,1,1,1))


#使用するスタンモデルのコンパイル
stanmodel = stan_model(file = "nonlinearmodel_2_gross.stan")

#Stanの実行
rng_seed　= 1111
nchains = 4
iter = 5000
ncores = 4 

stanout = sampling(stanmodel, 
                   chains = nchains, 
                   iter = iter, 
                   cores = ncores,
                   seed = rng_seed,
                   control = list(adapt_delta = 0.9))


pars =  c( "offset", "pmax", "alpha", "beta", "sigma")
pars_fitted =  c( "fitted")
pars_expected =  c( "expected")
traceplot(stanout, pars = c(pars))
traceplot(stanout, pars = c(pars_fitted))
traceplot(stanout, pars = c(pars_expected))
print(stanout, prob=c(0.025, 0.5, 0.975), par = pars)


#期待値expexted
expected = summary(stanout, pars = pars_expected, probs = c(0.025, 0.5, 0.975))$summary[, c("mean","2.5%", "50%","97.5%")]
expected = as.data.frame(expected)
expected = expected %>% mutate(speed = seq(0, max(speed), length = M))


#作図
ylabel_photo = expression("Gross photosynthesis rate"~"["~mu*g~O[2]~gww^-1~min^-1~"]")
xlabel = expression("Upstream water velocity"~"["~cm~sec^-1~"]")


colnames(Gross_data)
Gross_data = Gross_data %>% mutate(label = factor(Distance_individuals, labels = c("3 cm", "6 cm", "8 cm")))
Fontsize = 19
#横軸上流流速で固定
fg1 = ggplot()+
  geom_line(aes(x = speed, y = mean, color = "期待値"), data = expected) +
  geom_ribbon(aes(x = speed, ymin = `2.5%`, ymax = `97.5%`, fill = "95%信用区間"), data = expected, alpha = 0.50) +
  geom_point(aes(x = upstream_u_mean*100, y = Slope, 
                 color = label), data = Gross_data)+
  labs(x = xlabel, y = ylabel_photo)+
  guides(color= guide_legend(override.aes = 
                               list(shape=c(19,19,19,NA), 
                                    linetype = c(0, 0, 0, 1),
                                    size = c(3,3,3,2)))) +
  ylim(0, NA)+
  # scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm", "Model"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

#保存-----
width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_gross_model_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg1
dev.off()


#Stan each density----------------------------------------------------------------------------------------
Gross_data = all_data %>%
  filter(Hz != 0, type == "GP")

N = nrow(Gross_data)#データの数
idx = as.integer(Gross_data$Density)
speed = Gross_data$upstream_u_mean*100                        #観測した速度
gross = Gross_data$Slope                                  #観測した総一次生産量 

#密度条件ごとに変化しないパラメーター
M = 50                                                  #期待値に使う結果の数
speed_increment = max(speed) / M
speed_predict = speed_increment * seq(0,(M-1))
prior_log_sigma = as.vector(c(1,1,1,1))               #各パラメータの事前分散
prior_log_mu =    as.vector(c(1,1,1,1))

#使用するスタンモデルのコンパイル
stanmodel_each = stan_model(file = "nonlinearmodel_2_gross_each_density.stan")

#Stanの実行
rng_seed　= 1111
nchains = 4
iter = 5000
ncores = 4 

stanout = sampling(stanmodel_each, 
                   chains = nchains, 
                   iter = iter, 
                   cores = ncores,
                   seed = rng_seed,
                   control = list(adapt_delta = 0.9))
# http://mc-stan.org/loo/reference/compare.html
pars =  c( "offset", "pmax", "alpha", "beta", "sigma")
pars_ypred_High =  c( "ypred_01")
pars_ypred_Middle =  c( "ypred_02")
pars_ypred_Low =  c( "ypred_03")

pars_yfit_High =  c( "yfit_01")
pars_yfit_Middle =  c( "yfit_02")
pars_yfit_Low =  c( "yfit_03")

traceplot(stanout, pars = c(pars_ypred_High))
traceplot(stanout, pars = c(pars_yfit_High))

print(stanout, prob=c(0.025, 0.5, 0.975), par = pars_ypred_High)

#予測値
ypred_High = summary(stanout,   pars = pars_ypred_High,   probs = c(0.10, 0.5, 0.90))$summary[, c("mean","10%", "50%","90%")] %>% as.tibble()
ypred_Middle = summary(stanout, pars = pars_ypred_Middle, probs = c(0.10, 0.5, 0.90))$summary[, c("mean","10%", "50%","90%")] %>% as.tibble()
ypred_Low = summary(stanout,    pars = pars_ypred_Low,    probs = c(0.10, 0.5, 0.90))$summary[, c("mean","10%", "50%","90%")]  %>% as.tibble()
ypred = list(High = ypred_High, Middle = ypred_Middle, Low = ypred_Low) %>% bind_rows(.id = "Density")
ypred = ypred %>% mutate(speed = rep(speed_predict, times= 3))
ypred = ypred %>% rename(low = `10%`,
                         med = `50%`,
                         hig = `90%`)
#作図
ylabel_photo = expression("Gross photosynthesis rate"~"["~mu*g~O[2]~gww^-1~min^-1~"]")
xlabel = expression("Upstream water velocity"~"["~cm~sec^-1~"]")

colnames(Gross_data)
Fontsize = 19

#横軸上流流速で固定
fg1 = 
  ggplot()+
  geom_point(aes(x = upstream_u_mean*100, y = Slope, color = Density), data = Gross_data)+
  geom_ribbon(aes(x = speed, ymin = low, ymax = hig, fill = Density), data = ypred, alpha = 0.50) +
  geom_line(aes(x = speed, y = mean, group = Density), data = ypred) +
  labs(x = xlabel, y = ylabel_photo)+
  # guides(color= guide_legend(override.aes = 
  #                              list(shape=c(19,19,19,NA), 
  #                                   linetype = c(0, 0, 0, 1),
  #                                   size = c(3,3,3,2)))) +
  ylim(0, NA)+
  # scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm", "Model"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL
fg1

N =43
N = 50
ypred_01_43 = as.matrix(stanout, pars = "ypred_01")[,N]
ypred_02_43 = as.matrix(stanout, pars = "ypred_02")[,N]
ypred_03_43 = as.matrix(stanout, pars = "ypred_03")[,N]

test = data_frame(ypred = c(ypred_01_43, ypred_02_43, ypred_03_43),
           density = rep(c(1,2,3), each = length(ypred_01_43)),
           N = rep(seq(1:length(ypred_01_43)), 3)) %>% 
  mutate(density = factor(density, label = c("High", "Middle", "Low")))

ggplot(test) +
  geom_histogram(aes(x=ypred)) +
  facet_wrap("density", ncol = 1)

test %>% 
  spread(key="density", value="ypred") %>% 
  mutate(HM = High - Middle,
         HL = High - Low,
         ML = Middle - Low) %>% 
  select(HM, HL, ML) %>% 
  gather() %>% 
  ggplot() +
  geom_histogram(aes(x=value)) +
  geom_vline(xintercept = 0) +
  facet_wrap("key", ncol = 1)

test %>% 
  spread(key="density", value="ypred") %>% 
  mutate(HM = High - Middle,
         HL = High - Low,
         ML = Middle - Low) %>% 
  select(HM, HL, ML) %>% 
  gather() %>% 
  group_by(key) %>% 
  summarise(Pval = sum(value>=0)/length(value))

#保存-----
width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_gross_model_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg1
dev.off()



# パラメーターの比較----------------------------------------------------------------------------------------------------
#ヒストグラムで各パラメーターの確認
stan_hist(stanout, par=pars_expected, bins = 30)

