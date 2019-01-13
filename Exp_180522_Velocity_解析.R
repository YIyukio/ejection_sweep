#' ---
#' title: "ADVの解析"
#' author: "Inoue Yukio"
#' date: "2018 May 24"
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


#downstrean
#再読み込み----------------------------------------
fname = dir("../回流水槽光合成実験5/Modified_Data/", pattern = "vector", full.names = T)
vector_all = read.csv(fname[1], header = T)[,-1]
vector_mean = read.csv(fname[3], header = T)[,-1]
vector_zansa = read.csv(fname[6], header = T)[,-1]
vector_zansa_mean = read.csv(fname[5], header = T)[,-1]

head(vector_all)

xlabel = expression("Elapsed time"~"["~sec~"]")
ylabel = expression("Water velocity"~"["~cm~sec^-1~"]")

fg1 = ggplot(vector_all %>% filter(Hz == 2, Distance == 41, Depth == 6, Density == "H", time <= 2))+
  geom_line(aes(x = time, y = u))+
  labs(x = xlabel, y = ylabel)

width = 100
height = 100   
png(file = "~/回流水槽光合成実験5/fg/fg_velosity_u.png",
    type = "cairo",
    width = width, 
    height = height,
    units = "mm",
    res = 600,
    family = "Times New Roman")
fg1
dev.off()

#Initial ajustment length (Xd)の算出----------------------------
#Chen et al.,2013より引用
# w(x) = ws

#まずは作図
head(vector_zansa_mean)

ggplot(vector_mean %>%
         mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         filter(Depth2 == "edge"))+
  geom_point(aes(x = Distance, y = w_mean, color = Density))+
  geom_line(aes(x = Distance, y = w_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")


ggplot(vector_mean %>%
         mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         filter(Depth2 == "above"))+
  geom_point(aes(x = Distance, y = w_mean, color = Density))+
  geom_line(aes(x = Distance, y = w_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")


ggplot(vector_mean %>%
         mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         filter(Depth2 == "in"))+
  geom_point(aes(x = Distance, y = w_mean, color = Density))+
  geom_line(aes(x = Distance, y = w_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")

###
ggplot(vector_mean %>%
         mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         filter(Depth2 == "in"))+
  geom_point(aes(x = Distance, y = u_mean, color = Density))+
  geom_line(aes(x = Distance, y = u_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")

ggplot(vector_mean %>%
         mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         filter(Depth2 == "edge"))+
  geom_point(aes(x = Distance, y = u_mean, color = Density))+
  geom_line(aes(x = Distance, y = u_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")

ggplot(vector_mean %>%
         mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         filter(Depth2 == "above"))+
  geom_point(aes(x = Distance, y = u_mean, color = Density))+
  geom_line(aes(x = Distance, y = u_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")


###
ggplot(vector_zansa_mean %>% 
         filter(Depth2 == "edge"))+
  geom_point(aes(x = Distance, y = -stress_mean, color = Density))+
  geom_line(aes(x = Distance, y = -stress_mean, color = Density))+
  facet_wrap("Hz")

ggplot(vector_zansa_mean %>% 
         filter(Depth2 == "in"))+
  geom_point(aes(x = Distance, y = -stress_mean, color = Density))+
  geom_line(aes(x = Distance, y = -stress_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")

ggplot(vector_zansa_mean %>% 
         filter(Depth2 == "above"))+
  geom_point(aes(x = Distance, y = -stress_mean, color = Density))+
  geom_line(aes(x = Distance, y = -stress_mean, color = Density))+
  facet_wrap("Hz",scales = "free_y")


###
ggplot(vector_zansa_mean %>% 
         filter(Depth2 == "edge"))+
  geom_point(aes(x = Distance, y = TKE, color = Density))+
  geom_line(aes(x = Distance, y = TKE, color = Density))+
  facet_wrap("Hz")

ggplot(vector_zansa_mean %>% 
         filter(Depth2 == "in"))+
  geom_point(aes(x = Distance, y = TKE, color = Density))+
  geom_line(aes(x = Distance, y = TKE, color = Density))+
  facet_wrap("Hz")

ggplot(vector_zansa_mean %>% 
         filter(Depth2 == "above"))+
  geom_point(aes(x = Distance, y = TKE, color = Density))+
  geom_line(aes(x = Distance, y = TKE, color = Density))+
  facet_wrap("Hz")

#流速平均
X = ggplot()+
  geom_line(aes(y = u_mean, x = Depth, color = as.factor(Hz)),
            alpha = 0.4, 
            vector_mean)+
  geom_pointrange(aes(y = u_mean, ymin = u_mean - u_sd,
                      ymax = u_mean + u_sd,
                      x = Depth, color = as.factor(Hz)), size=0.2, vector_mean)+
  coord_flip()+
  facet_grid(Density ~ Distance)

Y = ggplot()+
  geom_line(aes(y = v_mean, x = Depth, color = as.factor(Hz)),
            alpha = 0.4, 
            vector_mean)+
  geom_pointrange(aes(y = v_mean, ymin = v_mean - v_sd,
                      ymax = v_mean + v_sd,
                      x = Depth, color = as.factor(Hz)), size=0.2, vector_mean)+
  coord_flip()+
  facet_grid(Density ~ Distance)

Z = ggplot()+
  geom_line(aes(y = w_mean, x = Depth, color = as.factor(Hz)),
            alpha = 0.4, 
            vector_mean)+
  geom_pointrange(aes(y = w_mean, ymin = w_mean - w_sd,
                      ymax = w_mean + w_sd,
                      x = Depth, color = as.factor(Hz)), size=0.2, vector_mean)+
  coord_flip()+
  facet_grid(Density ~ Distance)

X
Y
Z

#藻場内の流速
head(vector_mean)
upstream_u_mean = vector_mean %>%
  filter(Distance == 26) %>%
  select(Hz, Depth, Density, upstream_u_mean = u_mean) %>% 
  group_by(Hz, Density) %>% 
  summarise(upstream_u_mean = mean(upstream_u_mean))

vector_mean = vector_mean %>% left_join(upstream_u_mean, by = c("Hz", "Density"))

ggplot(vector_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                              Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(upstream_u_mean, Depth, Depth2, Density, Distance, u_mean))+
  geom_point(aes(x = upstream_u_mean, y = u_mean, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = u_mean, color = Density))+
  facet_grid(Depth2 ~ Distance)

ggplot(vector_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                              Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(upstream_u_mean, Depth, Depth2, Density, Distance, v_mean))+
  geom_point(aes(x = upstream_u_mean, y = v_mean, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = v_mean, color = Density))+
  facet_grid(Depth2 ~ Distance)

ggplot(vector_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                              Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(upstream_u_mean, Depth, Depth2, Density, Distance, w_mean))+
  geom_point(aes(x = upstream_u_mean, y = w_mean, color = Density))+
  geom_line(aes(x = upstream_u_mean, y = w_mean, color = Density))+
  facet_grid(Depth2 ~ Distance)

ggplot(vector_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                              Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(upstream_u_mean, Depth, Depth2, Density, Distance, u_mean, v_mean, w_mean))+
  geom_point(aes(x = upstream_u_mean, y = abs(w_mean)+abs(u_mean)+abs(v_mean), color = Density))+
  geom_line(aes(x = upstream_u_mean, y = abs(w_mean)+abs(u_mean)+abs(v_mean), color = Density))+
  facet_grid(Depth2 ~ Distance)

#残差2乗平均
zansa_meanX=ggplot()+
  geom_point(aes(y = zansaX_mean, x = Depth, color = as.factor(Hz)), size=3, alpha=0.5, vector_zansa_mean)+
  geom_line(aes(y = zansaX_mean, x = Depth, color = as.factor(Hz)), vector_zansa_mean)+
  # geom_hline(yintercept = 0, color = "red")+
  # geom_vline(xintercept = 12, color = "blue")+
  # scale_y_continuous(limits=c(-0.025 , 0.025))+
  coord_flip()+
  facet_grid(Density ~ Distance)

zansa_meanY=ggplot()+
  geom_point(aes(y = zansaY_mean, x = Depth, color = as.factor(Hz)), size=3, alpha=0.5, vector_zansa_mean)+
  geom_line(aes(y = zansaY_mean, x = Depth, color = as.factor(Hz)), vector_zansa_mean)+
  # geom_hline(yintercept = 0, color = "red")+
  # geom_vline(xintercept = 12, color = "blue")+
  # scale_y_continuous(limits=c(-0.025 , 0.025))+
  coord_flip()+
  facet_grid(Density ~ Distance)

zansa_meanZ=ggplot()+
  geom_point(aes(y = zansaZ_mean, x = Depth, color = as.factor(Hz)), size=3, alpha=0.5, vector_zansa_mean)+
  geom_line(aes(y = zansaZ_mean, x = Depth, color = as.factor(Hz)), vector_zansa_mean)+
  # geom_hline(yintercept = 0, color = "red")+
  # geom_vline(xintercept = 12, color = "blue")+
  # scale_y_continuous(limits=c(-0.025 , 0.025))+
  coord_flip()+
  facet_grid(Density ~ Distance)

zansa_meanX
zansa_meanY
zansa_meanZ

#藻場内の残渣二乗平均
head(vector_zansa_mean)
upstream_u_mean = vector_mean %>%
  filter(Distance == 26) %>%
  select(Hz, Depth, Density, upstream_u_mean = u_mean) %>% 
  group_by(Hz, Density) %>% 
  summarise(upstream_u_mean = mean(upstream_u_mean))

vector_zansa_mean = vector_zansa_mean %>% left_join(upstream_u_mean, by = c("Hz", "Density"))


ggplot(vector_zansa_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                                    Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(Hz, Depth, Depth2, Density, upstream_u_mean, Distance, zansaX_mean, zansaY_mean, zansaZ_mean, stress_mean, TKE))+
  geom_line(aes(x = upstream_u_mean, y = zansaX_mean, color = Density))+
  geom_point(aes(x = upstream_u_mean, y = zansaX_mean, color = Density))+
  facet_grid(Depth2 ~ Distance)

ggplot(vector_zansa_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                                    Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(Hz, Depth, Depth2, Density, upstream_u_mean, Distance, zansaX_mean, zansaY_mean, zansaZ_mean, stress_mean, TKE))+
  geom_line(aes(x = upstream_u_mean, y = zansaY_mean, color = Density))+
  geom_point(aes(x = upstream_u_mean, y = zansaY_mean, color = Density))+
  facet_grid(Depth2 ~ Distance)

ggplot(vector_zansa_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                                    Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(Hz, Depth, Depth2, Density, upstream_u_mean, Distance, zansaX_mean, zansaY_mean, zansaZ_mean, stress_mean, TKE))+
  geom_line(aes(x = upstream_u_mean, y = zansaZ_mean, color = Density))+
  geom_point(aes(x = upstream_u_mean, y = zansaZ_mean, color = Density))+
  facet_grid(Depth2 ~ Distance)

ggplot(vector_zansa_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                              Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(Hz, Depth, Depth2, Density, upstream_u_mean, Distance, zansaX_mean, zansaY_mean, zansaZ_mean, stress_mean, TKE))+
  geom_line(aes(x = upstream_u_mean, y = TKE, color = Density))+
  geom_point(aes(x = upstream_u_mean, y = TKE, color = Density))+
  facet_grid(Depth2 ~ Distance)

ggplot(vector_zansa_mean %>% mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "edge", Depth),
                                    Depth2 = recode(Depth2, "16" = "above", "6" = "in")) %>% 
         select(Hz, Depth, Depth2, Density, upstream_u_mean, Distance, zansaX_mean, zansaY_mean, zansaZ_mean, stress_mean, TKE))+
  geom_line(aes(x = upstream_u_mean, y = stress_mean, color = Density))+
  geom_point(aes(x = upstream_u_mean, y = stress_mean, color = Density))+
  facet_grid(Depth2 ~ Distance)

#TKE
TKE=ggplot()+
  geom_point(aes(y = TKE, x = Depth, color = as.factor(Hz)), size=3, alpha=0.5, vector_zansa_mean)+
  geom_line(aes(y = TKE, x = Depth, color = as.factor(Hz)), vector_zansa_mean)+
  # geom_hline(yintercept = 0, color = "red")+
  # geom_vline(xintercept = 12, color = "blue")+
  # scale_y_continuous(limits=c(-0.025 , 0.025))+
  coord_flip()+
  facet_grid(Density ~ Distance)

TKE

#stress
stress=ggplot()+
  geom_point(aes(y = stress_mean, x = Depth, color = as.factor(Hz)), size=3, alpha=0.5, vector_zansa_mean)+
  geom_line(aes(y = stress_mean, x = Depth, color = as.factor(Hz)), vector_zansa_mean)+
  # geom_hline(yintercept = 0, color = "red")+
  # geom_vline(xintercept = 12, color = "blue")+
  # scale_y_continuous(limits=c(-0.025 , 0.025))+
  coord_flip()+
  facet_grid(Density ~ Distance)

stress

#X-Z残差相関
Hz05 = vector_zansa %>%
  filter(Hz == 0.5)
Hz1 = vector_zansa %>%
  filter(Hz == 1)
Hz2 = vector_zansa %>%
  filter(Hz == 2)
Hz4 = vector_zansa %>%
  filter(Hz == 4)
Hz6 = vector_zansa %>%
  filter(Hz == 6)
Hz8 = vector_zansa %>%
  filter(Hz == 8)
Hz10 = vector_zansa %>%
  filter(Hz == 10)
Hz15 = vector_zansa %>%
  filter(Hz == 15)

Hz05_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz05)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.2 , 0.2))+
  scale_x_continuous(limits=c(-0.2 , 0.2))+
  facet_grid(Density + Depth ~ Distance)

Hz1_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz1)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.2 , 0.2))+
  scale_x_continuous(limits=c(-0.2 , 0.2))+
  facet_grid(Density + Depth ~ Distance)

Hz2_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz2)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.2 , 0.2))+
  scale_x_continuous(limits=c(-0.2 , 0.2))+
  facet_grid(Density + Depth ~ Distance)

Hz4_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz4)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.2 , 0.2))+
  scale_x_continuous(limits=c(-0.2 , 0.2))+
  facet_grid(Density + Depth ~ Distance)

Hz6_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz6)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.3 , 0.3))+
  scale_x_continuous(limits=c(-0.3 , 0.3))+
  facet_grid(Density + Depth ~ Distance)

Hz8_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz8)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.4 , 0.4))+
  scale_x_continuous(limits=c(-0.4 , 0.4))+
  facet_grid(Density + Depth ~ Distance)

Hz10_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz10)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.7 , 0.7))+
  scale_x_continuous(limits=c(-0.7 , 0.7))+
  facet_grid(Density + Depth ~ Distance)

Hz15_zansa = ggplot()+
  geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz15)+
  geom_hline(yintercept = 0, color = "red")+
  geom_vline(xintercept = 0, color = "red")+
  scale_y_continuous(limits=c(-0.7 , 0.7))+
  scale_x_continuous(limits=c(-0.7 , 0.7))+
  facet_grid(Density + Depth ~ Distance)

Hz05_zansa
Hz1_zansa
Hz2_zansa
Hz4_zansa
Hz6_zansa
Hz8_zansa
Hz10_zansa
Hz15_zansa

# Hz20_zansa = ggplot()+
#   geom_point(aes(y = zansaZ, x = zansaX), size=1, alpha=0.1, Hz20)+
#   geom_hline(yintercept = 0, color = "red")+
#   geom_vline(xintercept = 0, color = "red")+
#   scale_y_continuous(limits=c(-0.7 , 0.7))+
#   scale_x_continuous(limits=c(-0.7 , 0.7))+
#   facet_grid(Density + Depth ~ Distance)

#saving plot-----------------------------------------
# library(gridExtra)
# library(gtable)
# library(grid)
# 
# gp1 = ggplotGrob(X)
# gp2 = ggplotGrob(Y)
# gp3 = ggplotGrob(Z)
# gp4 = ggplotGrob(zansa_meanX)
# gp5 = ggplotGrob(zansa_meanY)
# gp6 = ggplotGrob(zansa_meanZ)
# gp7 = ggplotGrob(TKE)
# gp8 = ggplotGrob(stress)
# gp9 = ggplotGrob(Hz05_zansa)
# gp10 = ggplotGrob(Hz1_zansa)
# gp11 = ggplotGrob(Hz2_zansa)
# gp12 = ggplotGrob(Hz4_zansa)
# gp13 = ggplotGrob(Hz6_zansa)
# gp14 = ggplotGrob(Hz8_zansa)
# gp15 = ggplotGrob(Hz15_zansa)
# gp16 = ggplotGrob(Hz20_zansa)
# 
# maxwidth = unit.pmax(gp1[["widths"]], gp2[["widths"]], gp3[["widths"]], gp4[["widths"]],
#                      gp5[["widths"]], gp6[["widths"]], gp7[["widths"]], gp8[["widths"]],
#                      gp9[["widths"]], gp10[["widths"]], gp11[["widths"]], gp12[["widths"]],
#                      gp13[["widths"]], gp14[["widths"]], gp15[["widths"]], gp16[["widths"]])
# 
# gp1[["widths"]] = maxwidth
# gp2[["widths"]] = maxwidth
# gp3[["widths"]] = maxwidth
# gp4[["widths"]] = maxwidth
# gp5[["widths"]] = maxwidth
# gp6[["widths"]] = maxwidth
# gp7[["widths"]] = maxwidth
# gp8[["widths"]] = maxwidth
# gp9[["widths"]] = maxwidth
# gp10[["widths"]] = maxwidth
# gp11[["widths"]] = maxwidth
# gp12[["widths"]] = maxwidth
# gp13[["widths"]] = maxwidth
# gp14[["widths"]] = maxwidth
# gp15[["widths"]] = maxwidth
# gp16[["widths"]] = maxwidth
# 
# width = 300
# height= 300
# 
# png(file = "fg_xval_mean.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp1)
# dev.off()
# 
# png(file = "fg_yval_mean.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp2)
# dev.off()
# 
# png(file = "fg_zval_mean.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp3)
# dev.off()
# 
# png(file = "fg_zansa_meanX.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp4)
# dev.off()
# 
# png(file = "fg_zansa_meanY.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp5)
# dev.off()
# 
# png(file = "fg_zansa_meanZ.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp6)
# dev.off()
# 
# png(file = "fg_TKE.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp7)
# dev.off()
# 
# png(file = "fg_stress.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp8)
# dev.off()
# 
# width = 300
# height= 300
# 
# png(file = "fg_Hz05_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp9)
# dev.off()
# 
# png(file = "fg_Hz1_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp10)
# dev.off()
# 
# png(file = "fg_Hz2_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp11)
# dev.off()
# 
# png(file = "fg_Hz4_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp12)
# dev.off()
# 
# png(file = "fg_Hz6_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp13)
# dev.off()
# 
# png(file = "fg_Hz8_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp14)
# dev.off()
# 
# png(file = "fg_Hz15_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp15)
# dev.off()
# 
# png(file = "fg_Hz20_zansa.png",
#     type = "cairo",
#     width = width,
#     height = height,
#     units = "mm",
#     res = 600,
#     family = "Times New Roman")
# grid.draw(gp16)
# dev.off()