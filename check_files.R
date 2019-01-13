library(tidyverse)
library(lubridate)

cnames_vectrino=c("time", "counter", "status",
                  "xvel", "yvel", "zvel", "z2vel",
                  "amp1", "amp2", "amp3", "amp4",
                  "snr1", "snr2", "snr3", "snr4",
                  "corr1", "corr2", "corr3", "corr4")

data = data_frame(fnames = dir("Data/Velocity/", pattern = "dat", full =T))
data =
  data %>%
  mutate(data = map(fnames, read_table,
         col_names = cnames_vectrino)) %>%
  mutate(fnames = basename(fnames)) %>%
  separate(fnames, c("Hz", "Depth", "Distance", "Density"))


data =
  data %>% arrange(Hz, Depth, Distance, Density) %>%
  mutate_at(vars(Depth, Distance),
            funs(as.numeric))

data %>%
  filter(Depth > 14) %>%
  mutate(data2 = map(data, function(X) {
    X %>% summarise_at(vars(xvel:z2vel), funs(mean,sd,var))
    })) %>%
  unnest(data2) %>% print(n = Inf)

read_lines("Data/HDR/2_16_41_H_180611.hdr", skip = 1, n_max = 20)






data2 =
  data_frame(fnames = dir("Data/HDR/", pattern = "hdr", full = T)) %>%
  mutate(hdr = map(fnames,
                   function(X) {
                     x = read_lines(X, skip = 1, n_max = 20)
                     t(x[c(4,10, 17)]) %>% as_tibble()
                     }))

data2 %>% unnest() %>% select(4) %>% print(n = Inf)
data2 =
  data2 %>%
  unnest(hdr) %>%
  mutate(measurement = str_extract(V1, "[0-9]{4}/[0-9]{2}/[0-9]{2}"),
           range = str_extract(V2, "[0-9]+\\.[0-9]+"),
         scaling = str_extract(V3, "[0-9].[0-9]?")) %>%
  mutate(scaling = as.numeric(scaling)) %>%
  select(-V1, -V2, -V3)

data2 %>% print(n = Inf)

data2 =
  data2 %>%
  mutate(fnames = basename(fnames)) %>%
  mutate(measurement = as.Date(measurement)) %>%
  separate(fnames, c("Hz", "Depth", "Distance", "Density", "Date")) %>%
  arrange(Hz, Depth, Distance, Density) %>%
  mutate_at(vars(Depth, Distance), funs(as.numeric)) %>%
  mutate(Date = ymd(paste0("20",Date)))

data3 = full_join(data, data2, by = c("Hz", "Depth", "Distance", "Density"))


data3 %>%
  filter(Depth > 14) %>%
  mutate(data2 = map(data, function(X) {
    X %>% summarise_at(vars(xvel), funs(mean))
  })) %>%
  arrange(Hz, Depth, Distance, Density) %>%
  unnest(data2) %>%
  group_by(Hz, Depth, Distance) %>%
  mutate(xvel_x = xvel / mean(xvel)) %>% print(n = Inf)


data3 %>% unnest(data) %>% pull(xvel) %>% range()


sprintf("%02d-%02d-%02d-%s.png", 5, 2, 12, L)

data3 %>%
  slice(1) %>%
  unnest(data) %>%
  ggplot() +
  geom_line(aes(x=time, y=xvel))

data3 %>%
  arrange(Hz, Depth, Distance, Density) %>%
  mutate(n = 1:n()) %>%
  mutate(out = pmap(list(data, Hz, Depth, Distance, Density, n ), function(X, Hz, Depth, Distance, Density, n) {
    xmean = X %>% pull(xvel) %>% mean()
    title = sprintf("%02d-%02d-%02d-%s (%4.2f cm/s)",
                    as.numeric(Hz),
                    as.numeric(Depth),
                    as.numeric(Distance), Density,
                    xmean * 100)
    print(title)
    p1 =
      ggplot(X) +
      geom_line(aes(x = time, y = xvel)) +
      geom_hline(yintercept = xmean, color = "skyblue") +
      ggtitle(title) +
      scale_y_continuous(limits = c(-1.5, 1.1))
    ggsave(p1, file = sprintf("PNG/%04d.jpg", n), width = 180, height = 80, units = "mm")
  }))


system("ffmpeg -r 2 -f image2 -i ./PNG/%04d.jpg -vcodec libx264 -an -vf pad='width=ceil(iw/2)*2:height=ceil(ih/2)*2' velocity_time_series.mp4 -y")




