library(tidyverse)
library(broom)
library(modelr)
library(glue)
library(lemon)

out = read_csv("Modified_Data/vector_all.csv")


get_ejections = function(X) {
    N = nrow(X)
    X %>%
      mutate(event = (H > 0)) %>%
      group_by(quadrant) %>%
      filter(event) %>%
      summarise(S = sum(event)) %>%
      mutate(P = S / N)
}

dset = out %>%
  mutate(data = map(data, function (X) {
    X %>%
      mutate_at(vars(u_out, v_out, w_out, w2_out),
                function(x) {x - mean(x)}) %>%
      rename(up = u_out,
             vp = v_out,
             wp = w_out,
             w2p = w2_out) %>%
      mutate(H = abs(up * wp) / (sqrt(mean(up^2)) * sqrt(mean(wp^2)))) %>%
      mutate(quadrant = rep(NA, length(up))) %>%
      mutate(quadrant = ifelse(wp > 0 & up > 0, 1, quadrant)) %>%
      mutate(quadrant = ifelse(wp < 0 & up > 0, 2, quadrant)) %>%
      mutate(quadrant = ifelse(wp < 0 & up < 0, 3, quadrant)) %>%
      mutate(quadrant = ifelse(wp > 0 & up < 0, 4, quadrant)) %>%
      mutate(quadrant = factor(quadrant, label = c("I", "II", "III", "IV")))
  }))


dset = dset %>%
    mutate(Depth2 = ifelse(Depth >= 7 & Depth <= 13, "Edge", Depth),
           Depth2 = recode(Depth2, "16" = "Above", "6" = "Within")) %>%
    mutate(Depth2 = factor(Depth2, levels = c("Above", "Edge", "Within"), order = TRUE)) %>%
    mutate(Density = factor(Density, levels = c("H", "M", "L"), order = TRUE))


# All Events (H = 0)

dset2 =
  dset %>%
  mutate(data = map(data, get_ejections))

# Quadrant I: Outward interaction
# Quadrant II: Ejection
# Quadrant III: Inward interaction
# Quadrant IV: Sweep

dset2 %>%
  filter(str_detect(Density, "H")) %>%
  unnest() %>%
  ggplot() +
  geom_boxplot(aes(x=quadrant, y = P)) +
  facet_rep_grid(Hz ~ Distance, repeat.tick.labels = TRUE) +
  theme(axis.line = element_line(linetype = 1)) +
  ggtitle("High density")


dset2 %>%
  filter(str_detect(Density, "M")) %>%
  unnest() %>%
  ggplot() +
  geom_boxplot(aes(x=quadrant, y = P)) +
  facet_rep_grid(Hz ~ Distance, repeat.tick.labels = TRUE) +
  theme(axis.line = element_line(linetype = 1))+
  ggtitle("Medium density")


dset2 %>%
  filter(str_detect(Density, "L")) %>%
  unnest() %>%
  ggplot() +
  geom_boxplot(aes(x=quadrant, y = P)) +
  facet_rep_grid(Hz ~ Distance, repeat.tick.labels = TRUE) +
  theme(axis.line = element_line(linetype = 1)) +
  ggtitle("Low density")

dset2 %>%
  unnest() %>%
  ggplot() +
  geom_boxplot(aes(x = quadrant, y = P, fill = Density)) +
  scale_fill_brewer(palette = "Dark2") +
  facet_rep_grid(Depth2 ~ Distance, repeat.tick.labels = TRUE) +
  theme(axis.line = element_line(linetype = 1)) +
  ggtitle("Proportion of events")



dset2 %>%
  unnest() %>%
  glm(P ~ Depth2 * Distance * Density, data = . ,
      family = binomial) %>%
  summary()




