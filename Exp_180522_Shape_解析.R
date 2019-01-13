library(tidyverse)
library(lubridate)
library(gridExtra)
library(gtable)
library(scales)
library(grid)
library(scales)
library(stringr)
library(xtable)
library(marelac) 

## 海藻の面積データの読み込み------------------------------------------
fname = dir("../回流水槽光合成実験5/Data/Area/", pattern = "csv", full = T)
Area = lapply(fname, read.csv, header = T)

#データの編集
names(Area) = basename(fname)

Area = Area %>%
  bind_rows(.id = "filenames") %>% 
  tidyr::separate(filenames, into = c("Density", "a"), sep = "_") %>% 
  select(Density, Area) %>% 
  mutate(Hz = rep(c(0,0.5,1,2,4,6,8,10,15), 3),
         Density = recode(Density, H = "High", M = "Middle", L = "Low")) %>% 
  group_by(Density) %>% 
  mutate(Area_rate = Area / max(Area)) 

ggplot(Area)+
  geom_point(aes(x = Hz, y = Area, color = Density)) +
  # geom_line(aes(x = Hz, y = Area_rate, color = Density)) +
  geom_smooth(aes(x = Hz, y = Area, color = Density), method = "gam", formula = y ~s(x, k = 7))+
  NULL

#面積の変化速度の算出---------------------------------------------
library(mgcv)
head(Area)

eps = 1/1000
Area_differ = Area %>%
  group_by(Density) %>% 
  nest() %>% 
  mutate(model = map(data, function(x){gam(Area ~ s(Hz, k = 7), data = x)}),
         summary = map(model, function(x){summary(x)}),
         df1 = map(data, function(x){data.frame(Hz = seq(min(na.omit(x$Hz)), max(na.omit(x$Hz)), length = 50))}),
         df1 = map2(df1, model, function(x, y){data.frame(Hz = x$Hz, y1 = predict(y, newdata = x))}),
         df2 = map(df1, function(x){x %>% select(Hz) %>% mutate(Hz = Hz + eps)}),
         df2 = map2(df2, model, function(x,y){data.frame(y2 = predict(y, newdata = x))}),
         differ = map2(df1, df2, function(x, y){bind_cols(x, y) %>% mutate(Area_differ = (y2-y1) / eps)})) %>% 
  select(1, 7) %>%
  unnest() %>% 
  select(1,2, Area_fit = y1, 5) %>% 
  return()

ggplot()+
  geom_line(aes(x = Hz, y = Area_differ/100, color = "Area_differ*100"), data = Area_differ)+
  geom_line(aes(x = Hz, y = Area_fit/1000, color = "Area_fit*1000"), data = Area_differ)+
  facet_grid(.~ Density)+ 
  NULL

# write.csv(Area, "../回流水槽光合成実験5/Modified_Data/Area.csv")
