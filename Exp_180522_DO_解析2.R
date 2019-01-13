#' ---
#' title: "回流水槽光合成実験DO_Flux補正方法を変更_ヨレモク"
#' author: "Inoue yukio"
#' date: "2018 May 23"
#' output: html_document
#' ---

library(tidyverse)
library(lubridate)
library(gridExtra)
library(gtable)
library(scales)
library(grid)
library(scales)
library(stringr)
library(xtable)
library(marelac) # for estimating Saturation concentration of Oxygen in water  

## データの読み込み------------------------------------------
fname = dir("../回流水槽光合成実験5/Data/DO/", pattern = "csv", full = T)
ID = read.csv(fname[1], header = T)
Trident = read.csv(fname[2], header = T, skip = 5,fileEncoding = "sjis" )

Trident = Trident %>%
  tidyr::unite(DATE, TIME, col = "Time", sep = " ") %>% 
  select(Time, ID = 3, Barometer = 6, DO = 8, Temp = 9, Sal = 12) %>% 
  mutate(Time = as.POSIXct(Time), Date = as.POSIXct(round(Time, "mins")))

fg_DO_Trident = ggplot(Trident)+
  geom_point(aes(x = Time, y = DO)) +
  facet_wrap("ID")
fg_DO_Trident

#データに実験条件の情報を追加:IDのファイルを結合---------------------------------------------------------------------------------
Time_ID = ID %>% 
  select(ID, Hz, Density, Light, starttime = Exp_start_time, endtime = Exp_end_time) %>%
  mutate_at(.vars = c("starttime", "endtime"), .funs = as.POSIXct) %>% 
  arrange(starttime, endtime, Light) %>% 
  select(ID, Hz, Density, Light, starttime, endtime)

Time_ID2 = Time_ID
Time_ID = Time_ID %>%
  filter(!Density %in% c("High") | !ID %in% c(5,6,1,2,9,10)) %>% 
  filter(!Density %in% c("Low") | !ID %in% c()) %>% 
  filter(!Density %in% c("Middle") | !ID %in% c()) 
  
Trident$Hz = NA
Trident$Light = NA
Trident$Density = NA
Trident$ID = NA

removetime = 900 #s 初めの切り取り時間
for (i in 1:nrow(Time_ID)) {

  Trident[Trident$Time >= Time_ID$starttime[i]+ removetime &
                       Trident$Time <= Time_ID$endtime[i], "Hz"] = Time_ID$Hz[i]
  Trident[Trident$Time >= Time_ID$starttime[i]+ removetime &
                       Trident$Time <= Time_ID$endtime[i], "Light"] = Time_ID$Light[i]
  Trident[Trident$Time >= Time_ID$starttime[i]+ removetime &
            Trident$Time <= Time_ID$endtime[i], "Density"] = as.character(Time_ID$Density[i])
  Trident[Trident$Time >= Time_ID$starttime[i]+ removetime &
            Trident$Time <= Time_ID$endtime[i], "ID"] = Time_ID$ID[i]
}

Trident = Trident %>%
  filter(!Hz == is.na(Hz)) %>% 
  group_by(ID, Hz, Light, Density) %>%
  arrange(Density, Hz, Light, Time) %>% 
  mutate(Elapsed_time = as.integer((Date - min(Date)))/60,
         N = row_number()) %>% 
  filter(Elapsed_time <= 40) %>%                                                  #データは1.5時間以下のみ使用
  select(ID, Hz, Light, Elapsed_time, Density, N, DO, Temp, Barometer, Sal)


#WWの追加
WW_sum = data.frame(Density = c("High", "Middle", "Low"),
                    WW_sum = c(453.4, 238.5, 153.8))

Trident = Trident %>% full_join(WW_sum, by = c("Density")) %>% 
  na.omit()

Trident$Density = factor(Trident$Density, levels =  c("Low", "Middle", "High"))

#glmへのあてはめ---------------------------------------------

Glmmodel = function(y, x){
  summary(glm(y ~ x, family = Gamma(link = "log")))
}　　　　　　　　　　　　　　　　　　　　　　　　

Slope = function(x){　　　　　　　　　　　　　　　
  coef(x)[2]
} 

Slope_error = function(x){　　　　　　　　　　　　　　　
  coef(x)[4]
} 


#DO
DO_analysis = Trident %>%
  group_by(ID, Hz, Light, WW_sum, Density) %>%
  nest() %>% 
  mutate(model = map(data, ~ glm(DO~Elapsed_time, data = ., family = Gamma(link = "log")))) %>% 
  mutate(summary = map(model, ~summary(.))) %>% 
  mutate(slope = map_dbl(model, ~Slope(.)),
         stderr = map_dbl(model, ~Slope_error(summary(.))),
         out = map(model, ~coef(summary(.)))) 

Exp1_DO_model = Trident %>% 
  group_by(ID, Hz, Density, WW_sum, Light) %>% 
  summarise(DO_Slope = Slope(Glmmodel(DO, Elapsed_time)), DO_Slope_error = Slope_error(Glmmodel(DO, Elapsed_time))) %>% 
  mutate(DO_Slope_WW = DO_Slope / WW_sum, DO_Slope_error_WW = DO_Slope_error / WW_sum)

ggplot()+
  geom_pointrange(aes(x = ID, y = DO_Slope, ymin = DO_Slope - DO_Slope_error, ymax = DO_Slope + DO_Slope_error, color = as.factor(Hz)), Exp1_DO_model)+
  facet_grid(.~ Light + Density)

ggplot()+
  geom_pointrange(aes(x = Hz, y = DO_Slope, ymin = DO_Slope - DO_Slope_error, ymax = DO_Slope + DO_Slope_error, color = as.factor(ID)), Exp1_DO_model)+
  facet_grid(.~ Light + Density)

ggplot()+
  geom_pointrange(aes(x = ID, y = DO_Slope_WW, ymin = DO_Slope_WW - DO_Slope_error_WW, ymax = DO_Slope_WW + DO_Slope_error_WW, color = as.factor(Hz)), Exp1_DO_model)+
  facet_grid(.~ Light + Density)

ggplot()+
  geom_pointrange(aes(x = Hz, y = DO_Slope_WW, ymin = DO_Slope_WW - DO_Slope_error_WW, ymax = DO_Slope_WW + DO_Slope_error_WW, color = ID), Exp1_DO_model)+
  facet_grid(.~ Light + Density)


#Temp
Temp_analysis=Trident %>%
  group_by(ID, Hz, Light, WW_sum, Density) %>%
  nest() %>% 
  mutate(model = map(data, ~glm(Temp~Elapsed_time, data = .))) %>% 
  mutate(summary = map(model, ~summary(.))) %>% 
  mutate(slope = map_dbl(model, ~Slope(.)),
         stderr = map_dbl(model, ~Slope_error(summary(.))),
         out = map(model, ~coef(summary(.)))) 

Exp1_Temp_model = Trident %>% 
  group_by(ID, Hz, Density, Light) %>% 
  summarise(Temp_Slope = Slope(Glmmodel(Temp, Elapsed_time)), Temp_Slope_error = Slope_error(Glmmodel(Temp, Elapsed_time)))

ggplot()+
  geom_pointrange(aes(x = ID, y = Temp_Slope, ymin = Temp_Slope - Temp_Slope_error, ymax = Temp_Slope + Temp_Slope_error, color = as.factor(Hz)), Exp1_Temp_model)+
  facet_grid(. ~ Light + Density)

ggplot()+
  geom_pointrange(aes(x = Hz, y = Temp_Slope, ymin = Temp_Slope - Temp_Slope_error, ymax = Temp_Slope + Temp_Slope_error, color = ID), Exp1_Temp_model)+
  facet_grid(.~ Light + Density)

#大気-海水間の酸素フラックスの算出-------------------------------------------------------------

Trident_flux = read.csv("../回流水槽光合成実験3/Modified_data/DO_air_water_Flux_data.csv") 

#データを用いてKの算出（ｋが負になると酸素の放出、正は吸収）-----------------------------
library(nlstools)
model = function (k , C, t, DO_saturation){
  C * exp(-k * t) + DO_saturation
}

#0.5Hzの場合： /h

DO = Trident_flux %>% filter(Hz == 0.5, ID == "ECSER02") %>% pull(DO)
DO_saturation = Trident_flux %>% filter(Hz == 0.5, ID == "ECSER02") %>% pull(DO_saturation)
t = Trident_flux %>% filter(Hz == 0.5, ID == "ECSER02") %>% pull(Elapsed_time)
df1 = data_frame(obs=DO, DO_saturation, t)
nlsout05 = nls(obs~model(k0,C0,t, DO_saturation), data = df1, start = list(k0=0.11, C0=-4))
Mdata05 = df1 %>% mutate(DO_Flux = (DO_saturation - obs)*coefficients(nlsout05)[[1]])

#4Hzの場合： /h

DO = Trident_flux %>% filter(Hz == 4) %>% pull(DO)
DO_saturation = Trident_flux %>% filter(Hz == 4) %>% pull(DO_saturation)
t = Trident_flux %>% filter(Hz == 4) %>% pull(Elapsed_time)
df1 = data_frame(obs=DO, DO_saturation, t)
nlsout4 = nls(obs~model(k0,C0,t, DO_saturation), data = df1, start = list(k0=0.11, C0=-4))
Mdata4 = df1 %>% mutate(DO_Flux = (DO_saturation - obs)*coefficients(nlsout4)[[1]])
df1 = df1 %>% mutate(dDO = DO_saturation  - obs)


#6Hzの場合： /h

DO = Trident_flux %>% filter(Hz == 6, ID == "ECSER02") %>% pull(DO)
DO_saturation = Trident_flux %>% filter(Hz == 6, ID == "ECSER02") %>% pull(DO_saturation)
t = Trident_flux %>% filter(Hz == 6, ID == "ECSER02") %>% pull(Elapsed_time)
df1 = data_frame(obs=DO, DO_saturation, t)
nlsout6 = nls(obs~model(k0,C0,t, DO_saturation), data = df1, start = list(k0=0.11, C0=-4))
Mdata6 = df1 %>% mutate(DO_Flux = (DO_saturation - obs)*coefficients(nlsout6)[[1]])
df1 = df1 %>% mutate(dDO = DO_saturation  - obs)


#gamによる推定方法（6Hz）
# ggplot(df1) +
#   geom_point(aes(x=t, y=dDO)) +
#   geom_smooth(aes(x=t, y=dDO),
#               method = "gam",
#               formula = y~s(x))
# 
# gamout = mgcv::gam(dDO ~ s(t, k = 40), data = df1)
# summary(gamout)
# mdata = data_frame(t=df1$t)
# ddo_predict = predict(gamout, mdata, type = "response")
# eps = 1e-7
# ddo_0 = predict(gamout, mdata, type = "lpmatrix")
# mdata_offset = mdata + eps
# ddo_1 = predict(gamout, mdata_offset, type = "lpmatrix")
# ddo_finite_difference = (ddo_1 - ddo_0)/eps
# ddo_finite_difference = ddo_finite_difference %*% coef(gamout) %>% as.numeric()
# mdata = mdata %>% mutate(slope = ddo_finite_difference)
# df1 = full_join(df1, mdata)
# 
# ggplot(df1) +
#   geom_point(aes(x=t, y=dDO)) +
#   geom_smooth(aes(x=t, y=dDO),
#               method = "gam",
#               formula = y~s(x)) +
#   geom_line(aes(x=t, y=slope), data = df1)
# 
# ggplot(df1) +
#   geom_point(aes(x=dDO, y = -slope))

#8Hzの場合： /h

DO = Trident_flux %>% filter(Hz == 8) %>% pull(DO)
DO_saturation = Trident_flux %>% filter(Hz == 8) %>% pull(DO_saturation)
t = Trident_flux %>% filter(Hz == 8) %>% pull(Elapsed_time)
df1 = data_frame(obs=DO, DO_saturation, t)
nlsout8 = nls(obs~model(k0,C0,t, DO_saturation), data = df1, start = list(k0=0.11, C0=-4))
Mdata8 = df1 %>% mutate(DO_Flux = (DO_saturation - obs)*coefficients(nlsout8)[[1]])
df1 = df1 %>% mutate(dDO = DO_saturation  - obs)

#20Hzの場合： /h

DO = Trident_flux %>% filter(Hz == 20) %>% pull(DO)
DO_saturation = Trident_flux %>% filter(Hz == 20) %>% pull(DO_saturation)
t = Trident_flux %>% filter(Hz == 20) %>% pull(Elapsed_time)
df1 = data_frame(obs=DO, DO_saturation, t)
nlsout20 = nls(obs~model(k0,C0,t, DO_saturation), data = df1, start = list(k0=0.11, C0=-4))
Mdata20 = df1 %>% mutate(DO_Flux = (DO_saturation - obs)*coefficients(nlsout20)[[1]])
df1 = df1 %>% mutate(dDO = DO_saturation  - obs)


#作図
fg_Hz05 = ggplot()+
  geom_line(aes(x = t, y = DO_Flux), Mdata05, color = "blue")+
  scale_y_continuous(limits = c(0,0.65))+
  NULL

fg_Hz4 = ggplot()+
  geom_line(aes(x = t, y = DO_Flux), Mdata4, color = "blue")+
  scale_y_continuous(limits = c(0, 0.65))+
  NULL

fg_Hz6 = ggplot()+
  geom_line(aes(x = t, y = DO_Flux), Mdata6, color = "blue")+
  scale_y_continuous(limits = c(0, 0.65))+
  NULL

fg_Hz8 = ggplot()+
  geom_line(aes(x = t, y = DO_Flux), Mdata8, color = "blue")+
  scale_y_continuous(limits = c(0, 0.65))+
  NULL

fg_Hz20 = ggplot()+
  geom_line(aes(x = t, y = DO_Flux), Mdata20, color = "blue")+
  scale_y_continuous(limits = c(0, 0.65))+
  NULL

gridExtra::grid.arrange(fg_Hz05, fg_Hz4, fg_Hz6, fg_Hz8, fg_Hz20)

#算出した酸素フラックス速度定数kを持ちいて、実験中の酸素フラックスを見積もる:使用するモデルはパターン2のDO_saturation を変動モデル-----------------------------------------------------------------------------
#現在（2018-03-21現在）ではすべてのHzのおいてｋを算出していないため、Hzとｋは直線回帰で近似できると仮定して、ｋを補足する

k = data_frame(k = c(coefficients(nlsout05)[[1]],NA, NA, coefficients(nlsout4)[[1]], coefficients(nlsout6)[[1]], coefficients(nlsout8)[[1]], NA, NA, coefficients(nlsout20)[[1]]),
               Hz = c(0.5, 1, 2, 4, 6, 8, 10, 15, 20))

#Hz6の推定が悪いので、NAにする
k = k %>% mutate(k = ifelse(Hz == 6, NA, k))

#一次式で近似
k_pre = data_frame(Hz = k$Hz, k = predict(lm(k ~ Hz, k), data = data_frame(Hz = k$Hz),newdata=data.frame(Hz = k$Hz)))

ggplot()+
  geom_point(aes(x = Hz, y = k), data = k)+
  geom_line(aes(x = Hz, y = k), data = k_pre)

#kを用いてDO濃度を補正
DO_offset_Flux = Trident %>% 
  mutate(DO_saturation = gas_O2sat(S = Sal, t = Temp),
         DO_diff = DO_saturation - DO) %>% 
  left_join(k_pre, by = "Hz") %>%
  mutate(DO_Flux_min = DO_diff * k) %>% 
  group_by(ID, Hz, Light, Density) %>% 
  mutate(DO_Flux_min_cum = cumsum(DO_Flux_min),
         DO_offset_Flux = DO - DO_Flux_min_cum) %>%
  return()

ggplot()+
  # geom_point(aes(x = Elapsed_time, y = DO_Flux_min, color = "Flux"), data = DO_offset_Flux)+
  # geom_point(aes(x = Elapsed_time, y = DO_offset_Flux, color = "offset"), data = DO_offset_Flux)+
  # geom_smooth(aes(x = Elapsed_time, y = DO_offset_Flux, color = "offset_line"), method = "glm", data = DO_offset_Flux)+
  # geom_point(aes(x = Elapsed_time, y = DO_saturation, color = "Saturation"), data = DO_offset_Flux)+
  # geom_point(aes(x = Elapsed_time, y = DO, color = "raw"), data = DO_offset_Flux)+
  # geom_smooth(aes(x = Elapsed_time, y = DO, color = "raw_line"), method = "glm", data = DO_offset_Flux)+
  geom_point(aes(x = Elapsed_time, y = DO_diff, color = "diff"), data = DO_offset_Flux)+
  geom_point(aes(x = Elapsed_time, y = DO_Flux_min_cum, color = "cum"), data = DO_offset_Flux)+
  # geom_line(aes(x = Elapsed_time, y = k, color = "k"), data = DO_offset_Flux)+
  facet_grid(Density ~ Light + Hz)


xlabel = expression("Elapsed time"~"["~min~"]")
ylabel = expression("Dissolved Oxygen"~"["~ml~L^-1~"]")
  
fg1 = ggplot(DO_offset_Flux %>% filter(Hz == 2, Light == T, Density == "High"))+
  geom_point(aes(x = Elapsed_time, y = DO_offset_Flux))+
  geom_smooth(aes(x = Elapsed_time, y = DO_offset_Flux), method = "lm", formula = y ~ x)+
  labs(x = xlabel, y = ylabel)

width = 100
height = 100   
png(file = "fg_Slope_DO.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg1
dev.off()

#DO_Slopを再計算
Exp1_DO_offset_Flux_model = DO_offset_Flux %>% 
  group_by(ID, Hz, Density, WW_sum, Light) %>% 
  summarise(DO_Slope = Slope(Glmmodel(DO_offset_Flux, Elapsed_time)),
            DO_Slope_error = Slope_error(Glmmodel(DO_offset_Flux, Elapsed_time)),
            Temp_mean = mean(Temp)) %>% 
  mutate(DO_Slope_WW = DO_Slope / WW_sum, DO_Slope_error_WW = DO_Slope_error / WW_sum)

ggplot()+
  geom_pointrange(aes(x = ID, y = DO_Slope, ymin = DO_Slope - DO_Slope_error, ymax = DO_Slope + DO_Slope_error, color = as.factor(Hz)), Exp1_DO_offset_Flux_model)+
  facet_grid(.~ Light + Density)

ggplot()+
  geom_pointrange(aes(x = Hz, y = DO_Slope, ymin = DO_Slope - DO_Slope_error, ymax = DO_Slope + DO_Slope_error, color = ID), Exp1_DO_offset_Flux_model)+
  facet_grid(.~ Light + Density)

ggplot()+
  geom_pointrange(aes(x = ID, y = DO_Slope_WW, ymin = DO_Slope_WW - DO_Slope_error_WW, ymax = DO_Slope_WW + DO_Slope_error_WW, color = as.factor(Hz)), Exp1_DO_offset_Flux_model)+
  facet_grid(.~ Light + Density)

ggplot()+
  geom_pointrange(aes(x = Hz, y = DO_Slope_WW, ymin = DO_Slope_WW - DO_Slope_error_WW, ymax = DO_Slope_WW + DO_Slope_error_WW, color = ID), Exp1_DO_offset_Flux_model)+
  facet_grid(.~ Light + Density)

#回流水槽の容量
V = 400 #L, 容量は0.4トン、海水の比重計算が必要かも

#kによる酸素フラックスの補正と総光合成速度の算出-------------------------------------------------------------------------------------
#総光合成速度の算出:kによる補正なし
Photosynthesis = 
  Exp1_DO_model %>% 
  ungroup() %>% 
  select(ID, Hz, Density, WW_sum, Light, DO_Slope) %>% 
  tidyr::spread(key = Light, value = DO_Slope) %>% 
  rename(Respiration = 'FALSE', Photosynthesis_net = 'TRUE') %>% 
  mutate(Respiration = Respiration,
         Photosynthesis_gross = Photosynthesis_net - Respiration,
         Respiration_WW = Respiration / WW_sum,
         Photosynthesis_net_WW = Photosynthesis_net / WW_sum,
         Photosynthesis_gross_WW = Photosynthesis_gross / WW_sum) %>% 
  mutate_at(.vars = c("Respiration", "Photosynthesis_gross", "Photosynthesis_net", "Respiration_WW", "Photosynthesis_net_WW", "Photosynthesis_gross_WW"), .funs = function(x){x*V})

fg_Photosynthesis = ggplot(Photosynthesis)+
  geom_point(aes(x = Hz, y = Photosynthesis_net_WW, color = "Photosynthesis_net"))+
  geom_smooth(aes(x = Hz, y = Photosynthesis_net_WW, color = "Photosynthesis_net"),  method = "loess", span = 1, se = F)+
  geom_point(aes(x = Hz, y = Photosynthesis_gross_WW, color = "Photosynthesis_gross"))+
  geom_smooth(aes(x = Hz, y = Photosynthesis_gross_WW, color = "Photosynthesis_gross"), method = "loess", span = 1, se = F)+
  geom_point(aes(x = Hz, y = Respiration_WW, color = "Respiration"))+
  geom_smooth(aes(x = Hz, y = Respiration_WW, color = "Respiration"), method = "loess", span = 1, se = F)+
  facet_wrap("Density")

# ggplot(Photosynthesis)+
#   geom_point(aes(x = Hz, y = Photosynthesis_net_WW, color = "Photosynthesis_net"))+
#   geom_smooth(aes(x = Hz, y = Photosynthesis_net_WW, color = "Photosynthesis_net"), method = "loess", span = 1, se = F)+
#   geom_point(aes(x = Hz, y = Photosynthesis_gross_WW, color = "Photosynthesis_gross"))+
#   geom_smooth(aes(x = Hz, y = Photosynthesis_gross_WW, color = "Photosynthesis_gross"), method = "loess", span = 1, se = F)+
#   geom_point(aes(x = Hz, y = Respiration_WW, color = "Respiration"))+
#   geom_smooth(aes(x = Hz, y = Respiration_WW, color = "Respiration"), method = "loess", span = 1, se = F)+
#   facet_wrap("Density")

# GNN
# Exp1_DO_model %>% select(ID:DO_Slope)%>% 
#   inner_join(DO_Flux %>% select(-DO_Flux_mean, -DO_Flux_sd), by = c("ID", "Hz", "Light")) %>% 
#   tidyr::unnest(data) %>% 
#   dplyr::mutate(DO_Flux_TRUE = DO_Slope  - DO_Flux) %>% 
#   ggplot() +
#   geom_point(aes(x=Elapsed_time, y=DO_Flux_TRUE,color=ID)) +
#   geom_hline(yintercept = 0) +
#   facet_grid(Hz~Light, scales = "free")





#総光合成速度の算出:kによる補正あり
bolzmann = 8.617e-5

RP = Exp1_DO_offset_Flux_model %>% 
  ungroup() %>%
  filter(Light == FALSE) %>% 
  mutate(K0 = mean(Temp_mean) + 273.15,
         K = Temp_mean + 273.15, 
         invK = (1 / K0 * bolzmann - 1 / K * bolzmann), 
         Respration_k_WW = -DO_Slope_WW,
         slresp = scale(log(Respration_k_WW), scale = F)) %>% 
  return()

model = lm(slresp ~ invK, data = RP) %>% summary()
Ea = model$coefficients %>% 
  as_tibble(rownames = "coefficient") %>% 
  filter(coefficient == "invK") %>% 
  pull(Estimate)

Ea_sub = 0.65

RP %>%
  mutate(Respration_k_WW_temp = Respration_k_WW * exp(-Ea_sub * invK)) %>%
  ggplot()+
  # geom_point(aes(x = K, y = Respration_k_WW, color = "Respration_k_WW"))+
  # geom_point(aes(x = K, y = Respration_k_WW_temp, color = "Respration_k_WW_temp"))
  geom_point(aes(x = Hz, y = -Respration_k_WW, color = "Respration_k_WW", shape = Density), size = 3)+
  geom_smooth(aes(x = Hz, y = -Respration_k_WW, color = "Respration_k_WW"), method = "loess", span = 1)+
  geom_point(aes(x = Hz, y = -Respration_k_WW_temp, color = "Respration_k_WW_temp", shape = Density), size = 3)+
  geom_smooth(aes(x = Hz, y =-Respration_k_WW_temp, color = "Respration_k_WW_temp"), method = "loess", span = 1)


NP = Exp1_DO_offset_Flux_model %>% 
  ungroup() %>%
  filter(Light == TRUE) %>% 
  mutate(K0 = mean(Temp_mean) + 273.15,
         K = Temp_mean + 273.15, 
         invK = (1 / K0 * bolzmann - 1 / K * bolzmann), 
         Photosynthesis_k_WW = DO_Slope_WW,
         slresp = scale(log(Photosynthesis_k_WW), scale = F)) %>% 
  return()

model = lm(slresp ~ invK, data = NP) %>% summary()
Ea = model$coefficients %>% 
  as_tibble(rownames = "coefficient") %>% 
  filter(coefficient == "invK") %>% 
  pull(Estimate)

NP %>%
  mutate(Photosynthesis_k_WW_temp = Photosynthesis_k_WW * exp(-Ea_sub * invK)) %>%
  ggplot()+
  # geom_point(aes(x = K, y = Photosynthesis_k_WW, color = "Photosynthesis_k_WW"))+
  # geom_point(aes(x = K, y = Photosynthesis_k_WW_temp, color = "Photosynthesis_k_WW_temp"))
  geom_point(aes(x = Hz, y = Photosynthesis_k_WW, color = "Photosynthesis_k_WW", shape = Density), size = 3)+
  geom_smooth(aes(x = Hz, y = Photosynthesis_k_WW, color = "Photosynthesis_k_WW"), method = "loess", span = 1, se = T)+
  geom_point(aes(x = Hz, y = Photosynthesis_k_WW_temp, color = "Photosynthesis_k_WW_temp", shape = Density), size = 3)+
  geom_smooth(aes(x = Hz, y =Photosynthesis_k_WW_temp, color = "Photosynthesis_k_WW_temp"), method = "loess", span = 1, se = F)







model = lm(slresp ~ invK, data = RP)
summary(model)

model = lm(slresp ~ invK, data = NP)
summary(model)

ggplot(NP)+
  geom_point(aes(x = K, y = Photosynthesis_k_WW, color = as.factor(Hz), shape = Density), size = 5)

ggplot(NP)+
  geom_point(aes(x = invK, y = Photosynthesis_k_WW, color = as.factor(Hz), shape = Density), size = 5)

ggplot(RP)+
  geom_point(aes(x = K, y = Respration_k_WW, color = as.factor(Hz), shape = Density), size = 5)

ggplot(RP)+
  geom_point(aes(x = invK, y = Respration_k_WW, color = as.factor(Hz), shape = Density), size = 5)




Photosynthesis_k = Exp1_DO_offset_Flux_model %>% 
  ungroup() %>% 
  select(Hz, Density, WW_sum, Light, DO_Slope) %>%
  tidyr::spread(key = Light, value = DO_Slope) %>% 
  rename(Respiration_k = 'FALSE', Photosynthesis_net_k = 'TRUE') %>% 
  mutate(Respiration_k = Respiration_k,
         Photosynthesis_gross_k = Photosynthesis_net_k - Respiration_k,
         Respiration_k_WW = Respiration_k / WW_sum,
         Photosynthesis_net_k_WW = Photosynthesis_net_k / WW_sum,
         Photosynthesis_gross_k_WW = Photosynthesis_gross_k / WW_sum) %>% 
  mutate_at(.vars = c("Respiration_k", "Photosynthesis_gross_k", "Photosynthesis_net_k", "Respiration_k_WW", "Photosynthesis_net_k_WW", "Photosynthesis_gross_k_WW"), .funs = function(x){x*V})





fg_Photosynthesis_k = ggplot(Photosynthesis_k)+
  geom_point(aes(x = Hz, y = Photosynthesis_net_k_WW, color = "Photosynthesis_net_k"))+
  geom_smooth(aes(x = Hz, y = Photosynthesis_net_k_WW, color = "Photosynthesis_net_k"),  method = "loess", span = 1, se = F)+
  geom_point(aes(x = Hz, y = Photosynthesis_gross_k_WW, color = "Photosynthesis_gross_k"))+
  geom_smooth(aes(x = Hz, y = Photosynthesis_gross_k_WW, color = "Photosynthesis_gross_k"), method = "loess", span = 1, se = F)+
  geom_point(aes(x = Hz, y = Respiration_k_WW, color = "Respiration_k"))+
  geom_smooth(aes(x = Hz, y = Respiration_k_WW, color = "Respiration_k"), method = "loess", span = 1, se = F)+
  facet_wrap("Density")

ggplot(Photosynthesis_k %>% filter(!Hz %in% c() | !Density %in% c("High", "Middle", "Low")))+
  # geom_point(aes(x = Hz, y = Photosynthesis_net_k_WW, color = Density, shape = "Photosynthesis_net_k"), size = 3)+
  # geom_smooth(aes(x = Hz, y = Photosynthesis_net_k_WW, color = Density, linetype = "Photosynthesis_net_k"), method = "loess", span = 1, se = F)+
  geom_point(aes(x = Hz, y = Photosynthesis_gross_k_WW, color = Density, shape = "Photosynthesis_gross_k"), size = 3)+
  geom_line(aes(x = Hz, y = Photosynthesis_gross_k_WW, color = Density, shape = "Photosynthesis_gross_k"), size = 3)+
  # geom_smooth(aes(x = Hz, y = Photosynthesis_gross_k_WW, color = Density, linetype = "Photosynthesis_gross_k"), method = "loess", span = 1, se = F)+
  geom_line(aes(x = Hz, y = Respiration_k_WW, color = Density, shape = "Respiration_k"), size = 3)+
  # geom_smooth(aes(x = Hz, y = Respiration_k_WW, color = Density, linetype = "Respiration_k"), method = "loess", span = 1, se = F)+
  # facet_wrap("Density")
  NULL

ggplot(Photosynthesis_k %>% filter(!Hz %in% c(10) | !Density %in% c("High", "Middle", "Low")))+
  # geom_line(aes(x = as.integer(Density), y = Photosynthesis_gross_k_WW, color = as.factor(Hz), shape = "Photosynthesis_gross_k"), size = 3)+
  geom_line(aes(x = as.integer(Density), y = Photosynthesis_net_k_WW, color = as.factor(Hz), shape = "Photosynthesis_net_k"), size = 3)+
  # geom_line(aes(x = as.integer(Density), y = Respiration_k_WW, color = as.factor(Hz), shape = "Respiration_k"), size = 3)+
  facet_wrap("Hz")
  NULL

gridExtra::grid.arrange(fg_Photosynthesis, fg_Photosynthesis_k)


#データの保存-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# write.csv(Photosynthesis_k, "../回流水槽光合成実験5/Modified_Data/DO_slope.csv")

