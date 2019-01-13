#' ---
#' title: "回流水槽光合成実験_解析まとめ"
#' author: "Inoue yukio"
#' date: "2018 June 13"
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

#データの読み込み-------------------------------------------------------------------------------------------------------
fnames = dir("../回流水槽光合成実験5/Modified_Data/", pattern = "csv", full.names = T)[-1]
list_alldata = lapply(fnames, read.csv, header = T)
Object_name = sub("\\.", "", paste("", str_extract(basename(fnames), "[a-zA-Z_1-9]*\\.")))


for(i in 1:length(fnames)){
  eval(parse(text = paste(Object_name[i],
                          "=list_alldata[[", i, "]]", sep = "")))                                               #リストの各データをオブジェクトに入れる
}

Object_name #オブジェクト名を確認

#データの編集------------------------------------------------------------------------------------------------------------
#NAの確認
# sum(is.na(Area))
# sum(is.na(DO_slope))
# sum(is.na(light_mean))
# sum(is.na(vector_mean))
# sum(is.na(vector_zansa_mean))
# sum(is.na(vector_mean2))
# sum(is.na(vector_zansa_mean2))

#ジョイント
Vel_upstream = vector_mean %>% 
  filter(Distance == 26) %>% 
  select(Hz, Density, 
         upstream_u_mean = u_mean, 
         upstream_v_mean = v_mean, 
         upstream_w_mean = w_mean) %>% 
  group_by(Hz,Density) %>% 
  summarise_all(.funs = function(x){mean(x) * 100}) %>% 
  ungroup() %>% 
  mutate(Density = recode(Density, 
                          "H" = "High",
                          "L" = "Low",
                          "M" = "Middle"))

light_mean = light_mean %>% 
  ungroup() %>% 
  mutate(Density = recode(Density, 
                          "H" = "High",
                          "L" = "Low",
                          "M" = "Middle"))

vector_mean = vector_mean %>% 
  mutate(Density = recode(Density, 
                          "H" = "High",
                          "L" = "Low",
                          "M" = "Middle"))

vector_zansa_mean = vector_zansa_mean %>% 
  mutate(Density = recode(Density, 
                          "H" = "High",
                          "L" = "Low",
                          "M" = "Middle"))

vector_mean2 = vector_mean2 %>% 
  mutate(Density = recode(Density, 
                          "H" = "High",
                          "L" = "Low",
                          "M" = "Middle"))

vector_zansa_mean2 = vector_zansa_mean2 %>% 
  mutate(Density = recode(Density, 
                          "H" = "High",
                          "L" = "Low",
                          "M" = "Middle"))

all_data = light_mean %>% 
  full_join(Vel_upstream, by = c("Hz", "Density")) %>% 
  full_join(DO_slope, by = c("Hz", "Density")) %>%
  select(-X.x, -X.y) %>% 
  mutate(upstream_u_mean = ifelse(is.na(upstream_u_mean), 0, upstream_u_mean))

#滞留時間の計算---------------------------------------------------------------------------------------------------------------------------------------------
#流速の変化により藻場が変形するため、藻場の体積は変化する。また、密度によって藻場内の海藻の体積も変化する。この2点を滞留時間の計算に用いる。

#藻場の高さ
Seaweedbed_hight = data.frame(Density = rep(unique(all_data$Density), each = length(unique(all_data$Hz))),
                              Hz = rep(unique(all_data$Hz), length(unique(all_data$Density))),
                              Hight = c(12,12,12,12,12,11,10,9,8,
                                        12,12,12,11,11,10,9,7,6,
                                        12,12,12,12,11,10,8,7,6)) %>% 
  mutate(Hight = Hight + 1)
  
#上流下流の流速
vel_updown = vector_mean %>%
  filter(Distance %in% c(90, 26), Depth == 6) %>% 
  select(Hz, Depth, Distance, Density, u_mean) %>% 
  spread(key = Distance, value = u_mean) %>% 
  rename(u_26 = `26`, u_90 = `90`)

#藻場の縁の流速
vel_aroud = vector_mean %>% 
  inner_join(Seaweedbed_hight, by = c("Hz", "Density", "Depth" = "Hight")) %>% 
  select(Hz, Depth, Distance, Density, w_mean) %>% 
  spread(key = Distance, value = w_mean) %>% 
  rename(w_26 = `26`, w_41 = `41`, w_71 = `71`, w_90 = `90`) %>%
  full_join(vel_updown, by = c("Hz", "Density")) %>%
  select(-Depth.y) %>%
  rename(Depth = Depth.x) %>% 
  return()                                         #Depthは藻場の高さ


#藻場上部からの流入あり
Length = (90-26)/100 #m 藻場の長さ
Width = 0.30 #m 藻場の幅
Seaweed_volume = data.frame(Density = c("High", "Low", "Middle"), 
                            Seaweed_volume = c(630, 175, 227)/1000/1000) #ml -> m^3 海藻の体積

Tr_in = vel_aroud %>% 
  left_join(Seaweed_volume) %>% 
  mutate(Depth = Depth - 1,
         type = "in",                       #藻場に流入している量から計算
         w_26 = ifelse(w_26 < 0, -w_26, 0),
         w_41 = ifelse(w_41 < 0, -w_41, 0),
         w_71 = ifelse(w_71 < 0, -w_71, 0),
         w_90 = ifelse(w_90 < 0, -w_90, 0),
         u_90 = ifelse(u_90 < 0, -u_90, 0),
         u_90 = ifelse(is.na(u_90), 0, u_90), #u_90がNAとなっている場所があるので、他の流速を鑑みて0を代入
         u_26 = ifelse(u_26 < 0, 0, u_26),
         Total_bed_volume = Depth/100 * Length * Width,
         Net_bed_volume = Total_bed_volume - Seaweed_volume,
         # Tr_in = Net_bed_volume / ((w_26+w_41+w_71+w_90)/4 * Length*Width + u_26 * Depth/100*Width + u_90 * Depth/100*Width)) %>% 
         Tr_in = Net_bed_volume / ((w_41+w_71)/2 * Length*Width + u_26 * Depth/100*Width + u_90 * Depth/100*Width)) %>%
         return()

all_data = all_data %>% 
  full_join(Tr_in, by = c("Hz", "Density")) 

# write.csv(all_data, file = "../回流水槽光合成実験5/Modified_Data/all_data.csv")

## 海藻の角度データの読み込み------------------------------------------
fname = dir("../画像解析/Data/", pattern = "csv", full = T)
angle = lapply(fname, read.csv, header = T)

#データの編集
names(angle) = basename(fname)

all_data = angle %>%
  bind_rows(.id = "filenames") %>% 
  tidyr::separate(filenames, into = c("Density", "a"), sep = "_") %>% 
  select(Density, Angle) %>% 
  mutate(Hz = rep(c(0,0.5,1,2,4,6,8,10,15), 3),
         Density = recode(Density, H = "High", M = "Middle", L = "Low")) %>% 
  group_by(Density) %>% 
  mutate(angle_rate = Angle / max(Angle)) %>% 
  right_join(all_data, by = c("Density", "Hz"))

#角度の変化速度の算出---------------------------------------------
library(mgcv)

x = all_data
eps = 1/1000
angle_differ = all_data %>%
  group_by(Density) %>% 
  nest() %>% 
  mutate(model = map(data, function(x){gam(angle_rate ~ s(upstream_u_mean, k = 6), data = x)}),
         summary = map(model, function(x){summary(x)}),
         df1 = map(data, function(x){data.frame(upstream_u_mean = seq(min(na.omit(x$upstream_u_mean)), max(na.omit(x$upstream_u_mean)), length = 50))}),
         df1 = map2(df1, model, function(x, y){data.frame(upstream_u_mean = x$upstream_u_mean, y1 = predict(y, newdata = x))}),
         df2 = map(df1, function(x){x %>% select(upstream_u_mean) %>% mutate(upstream_u_mean = upstream_u_mean + eps)}),
         df2 = map2(df2, model, function(x,y){data.frame(y2 = predict(y, newdata = x))}),
         differ = map2(df1, df2, function(x, y){bind_cols(x, y) %>% mutate(angle_differ = (y2-y1) / eps)})) %>% 
  select(1, 7) %>%
  unnest() %>% 
  select(1,2, angle_fit = y1, 5) %>% 
  return()

## 海藻の面積データの読み込み------------------------------------------
fname = dir("../回流水槽光合成実験5/Data/Area/", pattern = "area.csv", full = T)
Area = lapply(fname, read.csv, header = T)

#データの編集
names(Area) = basename(fname)

all_data = Area %>%
  bind_rows(.id = "filenames") %>% 
  tidyr::separate(filenames, into = c("Density", "a"), sep = "_") %>% 
  select(Density, Area) %>% 
  mutate(Hz = rep(c(0,0.5,1,2,4,6,8,10,15), 3),
         Density = recode(Density, H = "High", M = "Middle", L = "Low")) %>% 
  group_by(Density) %>% 
  mutate(Area_rate = Area / max(Area))  %>% 
  right_join(all_data, by = c("Density", "Hz"))

#面積の変化速度の算出---------------------------------------------
library(mgcv)

eps = 1/1000
Area_differ = all_data %>%
  group_by(Density) %>% 
  nest() %>% 
  mutate(model = map(data, function(x){gam(Area_rate ~ s(upstream_u_mean, k = 6), data = x)}),
         summary = map(model, function(x){summary(x)}),
         df1 = map(data, function(x){data.frame(upstream_u_mean = seq(min(na.omit(x$upstream_u_mean)), max(na.omit(x$upstream_u_mean)), length = 50))}),
         df1 = map2(df1, model, function(x, y){data.frame(upstream_u_mean = x$upstream_u_mean, y1 = predict(y, newdata = x))}),
         df2 = map(df1, function(x){x %>% select(upstream_u_mean) %>% mutate(upstream_u_mean = upstream_u_mean + eps)}),
         df2 = map2(df2, model, function(x,y){data.frame(y2 = predict(y, newdata = x))}),
         differ = map2(df1, df2, function(x, y){bind_cols(x, y) %>% mutate(Area_differ = (y2-y1) / eps)})) %>% 
  select(1, 7) %>%
  unnest() %>% 
  select(1,2, Area_fit = y1, 5) %>% 
  return()

## 海藻の面積データ（cm^2）と体積密度の算出------------------------------------------
fname = dir("../回流水槽光合成実験5/Data/Area/", pattern = "area_cm.csv", full = T)
Area_cm = lapply(fname, read.csv, header = T)

#データの編集
names(Area_cm) = basename(fname)

all_data =Area_cm %>%
  bind_rows(.id = "filenames") %>% 
  tidyr::separate(filenames, into = c("Density", "a", "b"), sep = "_") %>% 
  select(Density, Area_cm = "Area") %>% 
  mutate(Hz = rep(c(0,0.5,1,2,4,6,8,10,15), 3),
         Density = recode(Density, H = "High", M = "Middle", L = "Low")) %>% 
  group_by(Density) %>% 
  mutate(Area_cm_rate = 1-(Area_cm / max(Area_cm)))  %>% 
  right_join(all_data, by = c("Density", "Hz")) %>% 
  mutate(Volume_density = Seaweed_volume / (Area_cm/10000 * Width),
         Seawatar_Volume = (Area_cm/10000 * Width) - Seaweed_volume,
         Seawatar_Volume_biomass = ((Area_cm/10000 * Width) - Seaweed_volume) / WW_sum) 


# #光合成速度をGAMで近似-----------------------------------------------------------------
# model1 = gam(Photosynthesis_gross_k_WW ~ s(upstream_u_mean, k = 8), data = all_data %>% filter(Density == "High"))
# summary(model1)
# # df1 = data.frame(upstream_u_mean = seq(min(na.omit(all_data$upstream_u_mean)), max(na.omit(all_data$upstream_u_mean)), length = 50))
# df1 = data.frame(upstream_u_mean = na.omit(all_data$upstream_u_mean))
# Photosynthesis_fit = data.frame(upstream_u_mean = df1$upstream_u_mean, y1 = predict(model1, newdata = df1))
# 
# fit_data = bind_cols(Photosynthesis_fit, Angle_differ_fit) %>% select(1, Photosynthesis_fit = 2, "Angle_differ_fit" = 6)
# 
# ggplot(fit_data)+
#   geom_point(aes(x = Angle_differ_fit, y = Photosynthesis_fit, color = upstream_u_mean), size = 3)+
#   NULL

#流動パラメーターの追加------------------------------------------------------------
all_data_velosity = all_data %>%
  full_join(vector_mean2, by = c("Hz", "Density")) %>% 
  select(Hz, Distance, upstream_u_mean = upstream_u_mean.x, Density, Depth2,
         Photosynthesis_gross_k_WW, Respiration_k_WW,
         u_mean, v_mean, w_mean, u_sd, v_sd, w_sd) %>% 
  mutate_at(.vars = vars("u_mean","v_mean","w_mean" , contains("sd")), .funs = function(x){x * 100}) %>%
  return()

all_data =
  vector_mean2 %>%
  mutate_at(.vars = vars("u_mean","v_mean","w_mean" , contains("sd")), .funs = function(x){x * 100}) %>%
  filter(Distance %in% c(41, 71), Depth2 %in% c("in")) %>%
  mutate(Velosity = sqrt((u_mean^2 + v_mean^2 + w_mean^2)/3)) %>% 
  group_by(Hz, Density) %>%
  summarise(Velosity_mean = mean(Velosity), Velosity_sd = sd(Velosity)) %>% 
  full_join(all_data, by = c("Hz", "Density")) %>% 
  # mutate(Seawatar_Flux_biomass = Seawatar_Volume_biomass * Velosity_mean) %>%  #意味のある数値か？
  return()

all_data_zansa = vector_zansa_mean2 %>%
  full_join(all_data, by = c("Hz", "Density")) %>% 
  select(Hz, Distance, upstream_u_mean = upstream_u_mean.x, Density, Depth2,
         Photosynthesis_gross_k_WW, Respiration_k_WW,
         stress_mean, TKE) %>% 
  return()

all_data = vector_zansa_mean2 %>%
  filter(Distance %in% c(41, 71), Depth2 %in% c("in")) %>%
  group_by(Hz, Density) %>%
  summarise(TKE_mean = mean(TKE), TKE_sd = sd(TKE),
            stress_mean_mean = mean(stress_mean), stress_mean_sd = sd(stress_mean)) %>%
  full_join(all_data, by = c("Hz", "Density")) %>% 
  return()

#滞留時間を藻場内の平均流速を用いて再度算出----------------------------------------------------
all_data = all_data %>% mutate(Tr_mean = Total_bed_volume / (Velosity_mean/100* Depth * Width))
  
#レイノルズ数の計算-------------------------------------
#Re = vL/u v:velocity, L:物体の大きさ
#計算に用いる物体のスケールLは藻場の高さ（m）, 海藻間の距離(m)とする。
u = 0.94*10^-6 #m/s 海水２５度の時の動粘度係数(Bilger and Atkinson 1992, Kays and Crawford 1993, Thomas and Atkinson 1997)

all_data = all_data %>% 
  mutate(Re_height = ((upstream_u_mean/100) * (Depth/100))/u,
         Distance_individuals = recode(Density, "High" = 3, "Middle" = 6, "Low" = 8),
         Re_distance = ((upstream_u_mean/100) * (Distance_individuals/100))/u,
         Re_height_distance = ((upstream_u_mean/100) * ((Distance_individuals+Depth)/100))/u) %>% 
  return()

#Densityの順位を変える
all_data$Density = factor(all_data$Density, levels = c("High", "Middle", "Low"))

#作図--------------------------------------------------------------------------------------------------------------------------------------------------------
ylabel_photo = expression("Gross photosynthesis rate"~"["~mu*g~O[2]~gww^-1~min^-1~"]")
xlabel = expression("Upstream water velocity"~"["~cm~sec^-1~"]")

Fontsize = 17
#横軸上流流速で固定
fg1 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Photosynthesis_gross_k_WW*1000, color = Density), size = 3)+
  geom_smooth(aes(x = upstream_u_mean, y = Photosynthesis_gross_k_WW*1000), method = "gam", formula = y~s(x))+
  labs(x = xlabel, y = ylabel_photo)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_gross_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg1
dev.off()

fg2 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Photosynthesis_net_k_WW, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Photosynthesis_net_k_WW, color = Density))+
  facet_grid(. ~ Density)+
  NULL

fg3 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Respiration_k_WW, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Respiration_k_WW, color = Density))+
  facet_grid(. ~ Density)+
  NULL

label_vel_bed = expression("Water velocity"~"["~cm~sec^-1~"]")
fg4 = ggplot(all_data)+
  geom_pointrange(aes(x = upstream_u_mean, y = Velosity_mean,
                      ymin = Velosity_mean - Velosity_sd/sqrt(4), ymax = Velosity_mean + Velosity_sd/sqrt(4), color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Velosity_mean, color = Density))+
  labs(x = xlabel, y = ylabel_vel_bed)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_vel_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg4
dev.off()


ylabel_tke = expression("Turbulent kinetic energy"~"["~cm^2~sec^-2~"]")
fg5 =  ggplot(all_data)+
  geom_pointrange(aes(x = upstream_u_mean, y = TKE_mean,
                      ymin = TKE_mean - TKE_sd/sqrt(4), ymax = TKE_mean + TKE_sd/sqrt(4), color = Density))+
  geom_line(aes(x = upstream_u_mean, y = TKE_mean, color = Density))+
  labs(x = xlabel, y = ylabel_tke)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_tke_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg5
dev.off()


fg6 = ggplot(all_data)+
  geom_pointrange(aes(x = upstream_u_mean, y = stress_mean_mean,
                 ymin = stress_mean_mean - stress_mean_sd, ymax = stress_mean_mean + stress_mean_sd, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = stress_mean_mean, color = Density))+
  facet_grid(. ~ Density)


ylabel_Tr = expression("Residence time"~"["~sec~"]")
fg7 = ggplot(all_data)+
  geom_line(aes(x = upstream_u_mean, y = Tr_mean, color = Density))+
  geom_point(aes(x = upstream_u_mean, y = Tr_mean, color = Density), size = 3)+
  labs(x = xlabel, y = ylabel_Tr)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  scale_y_log10()+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL
width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_Tr_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg7
dev.off()

fg8 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Area_cm, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Area_cm, color = Density))+
  facet_grid(. ~ Density)+
  NULL

ylabel_area = expression("Change of rate of seaweed bed area")
fg9 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Area_cm_rate, color = Density), size = 3)+
  geom_line(aes(x = upstream_u_mean, y = Area_cm_rate, color = Density))+
  labs(x = xlabel, y = ylabel_area)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_Arearate_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600, 
    family = "Times New Roman")
fg9
dev.off()

fg10 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Angle, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Angle, color = Density))+
  facet_grid(. ~ Density)+
  NULL

fg11 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = angle_rate, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = angle_rate, color = Density))+
  facet_grid(. ~ Density)+
  NULL

fg12 = ggplot(all_data)+
  geom_pointrange(aes(x = upstream_u_mean, y = Light_mean, ymax = Light_mean + Light_sd, ymin = Light_mean - Light_sd, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Light_mean, color = Density))+
  facet_grid(. ~ Density)+
  NULL

fg13 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Volume_density, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Volume_density, color = Density))+
  facet_grid(. ~ Density)+
  NULL

fg14 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Seawatar_Volume_biomass, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = Seawatar_Volume_biomass, color = Density))+
  facet_grid(. ~ Density)+
  NULL

ylabel_Reynolds = expression("Reynolds number")
fg15 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Re_height, color = Density), size = 3)+
  labs(x = xlabel, y = ylabel_Reynolds)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_Reynolds_h_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600, 
    family = "Times New Roman")
fg15
dev.off()

ylabel_Reynolds = expression("Reynolds number")
fg16 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Re_distance, color = Density), size = 3)+
  labs(x = xlabel, y = ylabel_Reynolds)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_Reynolds_dis_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600, 
    family = "Times New Roman")
fg16
dev.off()

ylabel_Reynolds = expression("Reynolds number")
fg17 = ggplot(all_data)+
  geom_point(aes(x = upstream_u_mean, y = Re_height_distance, color = Density), size = 3)+
  labs(x = xlabel, y = ylabel_Reynolds)+
  ylim(0, NA)+
  scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

width = 150
height = 150   
png(file = "~/回流水槽光合成実験5/fg/fg_Reynolds_h_dis_upstream.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600, 
    family = "Times New Roman")
fg17
dev.off()

fg1
fg2
fg3
fg4
fg5
fg6
fg7
fg8
fg9
fg10
fg11
fg12
fg13
fg14
fg15
fg16
fg17

#相関図(縦軸総光合成速度)
ggplot(all_data)+
  geom_point(aes(x = Velosity_mean, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  labs(y = ylabel)+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = TKE_mean, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  labs(y = ylabel)+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = stress_mean_mean, y = Photosynthesis_gross_k_WW*1000, color = Density))+
  # scale_x_log10()+
  labs(y = ylabel)+
  # facet_grid(. ~ Density)+
  NULL

width = 160
height = 170   
fg_photo_tr = ggplot(all_data)+
  geom_point(aes(x = Tr_mean, y = Photosynthesis_gross_k_WW*1000,
                 color = upstream_u_mean), size = 3)+
  scale_x_log10()+
  labs(x = ylabel_Tr, y = ylabel_photo)+
  scale_color_gradient2(name = expression("Upstream water velocity"~"["~cm~sec^-1~"]  "),
                        high = "red", mid = "gray50",
                        low = "green" , midpoint = mean(all_data$upstream_u_mean))+
  ylim(0, NA)+
  # scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

png(file = "~/回流水槽光合成実験5/fg/fg_photo_tr.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg_photo_tr
dev.off()

ggplot(all_data)+
  geom_point(aes(x = Area_cm, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  labs(y = ylabel)+
  # facet_grid(. ~ Density)+
  NULL

fg_photo_area = ggplot(all_data)+
  geom_point(aes(x = Area_cm_rate, y = Photosynthesis_gross_k_WW*1000, color = upstream_u_mean), size = 3)+
  # scale_x_log10()+
  labs(x = ylabel_area, y = ylabel_photo)+
  scale_color_gradient2(name = expression("Upstream water velocity"~"["~cm~sec^-1~"]  "),
                        high = "red", mid = "gray50",
                        low = "green" , midpoint = mean(all_data$upstream_u_mean))+
  ylim(0, NA)+
  # scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

png(file = "~/回流水槽光合成実験5/fg/fg_photo_arearate.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg_photo_area
dev.off()

fg_photo_reynols_dis = ggplot(all_data)+
  geom_point(aes(x = Re_distance, y = Photosynthesis_gross_k_WW*1000, color = upstream_u_mean), size = 3)+
  # scale_x_log10()+
  labs(x = ylabel_Reynolds, y = ylabel_photo)+
  scale_color_gradient2(name = expression("Upstream water velocity"~"["~cm~sec^-1~"]  "),
                        high = "red", mid = "gray50",
                        low = "green" , midpoint = mean(all_data$upstream_u_mean))+
  ylim(0, NA)+
  # scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

fg_photo_reynols_h = ggplot(all_data)+
  geom_point(aes(x = Re_height, y = Photosynthesis_gross_k_WW*1000, color = upstream_u_mean), size = 3)+
  # scale_x_log10()+
  labs(x = ylabel_Reynolds, y = ylabel_photo)+
  scale_color_gradient2(name = expression("Upstream water velocity"~"["~cm~sec^-1~"]  "),
                        high = "red", mid = "gray50",
                        low = "green" , midpoint = mean(all_data$upstream_u_mean))+
  ylim(0, NA)+
  # scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

fg_photo_reynols_h_dis = ggplot(all_data)+
  geom_point(aes(x = Re_height_distance, y = Photosynthesis_gross_k_WW*1000, color = upstream_u_mean), size = 3)+
  # scale_x_log10()+
  labs(x = ylabel_Reynolds, y = ylabel_photo)+
  scale_color_gradient2(name = expression("Upstream water velocity"~"["~cm~sec^-1~"]  "),
                        high = "red", mid = "gray50",
                        low = "green" , midpoint = mean(all_data$upstream_u_mean))+
  ylim(0, NA)+
  # scale_color_discrete(labels = c("3 cm", "6 cm", "8 cm"))+
  # facet_grid(. ~ Density)+
  theme_gray(Fontsize)+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "gray60"),
        panel.background = element_blank())+
  NULL

png(file = "~/回流水槽光合成実験5/fg/fg_photo_reynols_dis.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg_photo_reynols_dis
dev.off()

png(file = "~/回流水槽光合成実験5/fg/fg_photo_reynols_h.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg_photo_reynols_h
dev.off()


png(file = "~/回流水槽光合成実験5/fg/fg_photo_reynols_h_dis.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg_photo_reynols_h_dis
dev.off()

ggplot(all_data)+
  geom_point(aes(x = Angle, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = angle_rate, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Light_mean, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Volume_density, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Seawatar_Volume_biomass, y = Photosynthesis_gross_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

#相関図(縦軸呼吸速度)
ggplot(all_data)+
  geom_point(aes(x = Velosity_mean, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = TKE_mean, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = stress_mean_mean, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Tr_in, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Area_cm, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Area_cm_rate, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Angle, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = angle_rate, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Light_mean, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Volume_density, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Seawatar_Volume_biomass, y = Respiration_k_WW, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

#その他の相関
ggplot(all_data)+
  geom_point(aes(x = Velosity_mean, y = TKE_mean, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = Velosity_mean, y = Area_cm_rate, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL

ggplot(all_data)+
  geom_point(aes(x = TKE_mean, y = Area_cm_rate, color = Density))+
  # scale_x_log10()+
  # facet_grid(. ~ Density)+
  NULL





#all_dataの保存--------------------------------------
write_csv(all_data, path = "~/回流水槽光合成実験5/Modified_Data/all_data.csv")




