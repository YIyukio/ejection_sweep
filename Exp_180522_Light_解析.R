#' ---
#' title: "回流水槽光合成実験Light"
#' author: "Inoue yukio"
#' date: "2018 Jun 8"
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

#データの読み込み---------------------------------------------------------------------------------------------
fnames = dir("../回流水槽光合成実験5/Data/Light/", pattern = "csv", full = T)
listdata_light = lapply(fnames[-7], read.csv, header = T, skip = 1)
light_ID = read.csv(fnames[7], header = T)
names(listdata_light) = paste("", str_extract(basename(fnames[-7]), "Inoue_[a-z]{2,4}_[0-9]{1,6}_[LMH]"), sep = "")
listdata_light = lapply(listdata_light, function(x){x %>% select(Time = 2, Light_lux = 4)})

#データの編集--------------------------------------------------------------------------------------------------
light = listdata_light %>%
  bind_rows(.id = "filename") %>% 
  tidyr::separate(col = filename, into = c("Inoue", "Position", "a", "Density"), sep = "_") %>%
  select(Position, Time, Light_lux, Density) %>% 
  mutate(Time = as.POSIXct(parse_date_time(Time, orders = "%m/%d/%y %I:%M:%S %p", tz = "JST", locale = "ja_JP.utf8")))

light_ID = light_ID %>% 
  mutate(starttime = as.POSIXct(as.character(Starttime), tz = "JST")) %>% 
  select(-3)

light$Hz = NA
light$Density = NA

ggplot(light)+
  geom_point(aes(x = Time, y = Light_lux))

removetime = 30 #s 初めの切り取り時間
extracttime = removetime + 20 #s 使用する時間
for (i in 1:nrow(light_ID)) {
  
  light[light$Time >= light_ID$starttime[i] + removetime &
          light$Time <= light_ID$starttime[i] + extracttime, "Hz"] = light_ID$Hz[i]
  
  light[light$Time >= light_ID$starttime[i] + removetime &
          light$Time <= light_ID$starttime[i] + extracttime, "Density"] = paste(light_ID$Density[i])
}

light = light %>%
  na.omit() %>% 
  group_by(Hz, Density, Position) %>% 
  mutate(Elapsed_time = as.integer((Time - min(Time))))

ggplot()+
  geom_point(aes(x = Elapsed_time, y = Light_lux), light)+
  facet_grid(Density ~  Position + Hz)

#平均値、ばらつきを計算-----------------------------------------------------------------------------------------------
light_mean = light %>% group_by(Density, Hz) %>% 
  summarise(Light_mean = mean(Light_lux), Light_sd = sd(Light_lux))

ggplot()+
  geom_pointrange(aes(x = Hz, y = Light_mean, ymin = Light_mean - Light_sd,
                      ymax = Light_mean + Light_sd), light_mean)+
  geom_line(aes(x = Hz, y = Light_mean), light_mean)+
  facet_wrap("Density")

#データの保存---------------------------------------------------------------------------------------------------------
write.csv(light, "../回流水槽光合成実験5/Modified_Data/light_all.csv")
write.csv(light_mean, "../回流水槽光合成実験5/Modified_Data/light_mean.csv")






  

