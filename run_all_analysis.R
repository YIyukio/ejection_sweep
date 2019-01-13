library(tidyverse)
library(broom)
library(modelr)
library(glue)
library(lemon)
library(knitr)
library(kableExtra)
set.seed(2019)
options(mc.cores = 4)
col_types = "nnnc_nnnnn____nnnn"
col_names = c("Hz", "Depth", "Distance", "Density",
              "tau","x", "y", "z","z2",
              "u", "v", "w", "w2")
out = read_csv("Modified_Data/vector_all.csv",
               col_types = col_types,
               col_names = col_names, skip = 1)

# A function that makes functions is called a
# function factory.

calculate_TKE = function(up, vp, wp) {
  up*up + vp*vp + wp*wp
}
calculate_stress = function(up, wp) {
  up * wp
}
calculate_mean_stress = function(X) {
  X %>%
    summarise(mean_stress = mean(calculate_stress(up, wp))) %>%
    pull(mean_stress)
}

calculate_mean_velocity = function(X) {
  X %>%
    rowwise() %>%
    mutate(a  = sum(w, w2, na.rm=TRUE)) %>%
    ungroup() %>%
    select(u, v, w=a) %>%
    summarise_all(funs(mean = mean, sd = sd))
}

calculate_total_time = function(X) {
  # In seconds
  X %>%
    summarise(total_time = length(up) * 0.005) %>%
    pull(total_time)
}

identify_quadrant = function(up, wp) {
  case_when(
    (wp > 0 & up > 0) ~ 1,
    (wp > 0 & up < 0) ~ 2,
    (wp < 0 & up < 0) ~ 3,
    (wp < 0 & up > 0) ~ 4
  )
}

quadrant_analysis_function = function(H) {
  function(X, mean_stress) {
    N = nrow(X)
    X %>%
      mutate(hole_size = rep(H, N)) %>%
      mutate(Stress = calculate_stress(up, wp)) %>%
      mutate(in_hole = abs(Stress) / abs(mean_stress) <= H) %>%
      mutate(Quadrant = identify_quadrant(up, wp))
  }
}

calculate_residuals = function(X) {
  X %>%
    rowwise() %>%
    mutate(a  = sum(w, w2, na.rm=TRUE)) %>%
    ungroup() %>%
    mutate_at(vars(u, v, a), function(x) {
      z = mean(x, na.rm=T)
      x - z
    }) %>%
    rename(up = u, vp = v, wp = a) %>%
    ungroup()
}
get_slope_cf = function(X) {
  coef(summary(X)) %>%
    as_tibble(rownames = "Distance") %>%
    filter(str_detect(Distance, "Hz")) %>%
    mutate(Distance=str_extract(Distance, "[0-9]+")) %>%
    mutate(Distance = ifelse(is.na(Distance), 26, Distance),
           Distance = as.numeric(Distance)) %>%
    select(Distance,
           estimate = Estimate,
           std.error = `Std. Error`) %>%
    mutate(l95 = estimate - 1.96*std.error,
           u95 = estimate + 1.96*std.error)
}

out =
  out %>%
  group_by(Hz, Depth, Distance, Density) %>%
  nest()

out =
  out %>%
  mutate(Position = case_when((Depth > 13) ~ 1,
                              (Depth >= 7 & Depth <= 13) ~ 2,
                              (Depth < 7) ~ 3)) %>%
  mutate(Position = factor(Position,
                           labels = c("Above canopy",
                                      "Canopy edge",
                                      "Within canopy"),
                           order = TRUE)) %>%
  mutate(Density = factor(Density,
                          levels = c("H", "M", "L"),
                          labels = c("High", "Medium", "Low"),
                          order = TRUE)) %>%
  mutate(Hz = ifelse(Hz == 5, 0.5, Hz)) %>%
  mutate(Distance = Distance - 28) %>%
  select(-Depth)

run_analysis = function(HSIZE, out) {
  quadrant_analysis = quadrant_analysis_function(H = HSIZE)
  dset =
    out %>% slice(1:2) %>%
    # mutate(data = map(data, function(X){
    #   X %>% drop_na(x, y, z)
    # })) %>%
    mutate(sum = map(data, calculate_mean_velocity)) %>%
    unnest(sum) %>%
    mutate(data = map(data, calculate_residuals)) %>%
    mutate(mean_stress = map_dbl(data, calculate_mean_stress)) %>%
    mutate(total_time = map_dbl(data, calculate_total_time)) %>%
    mutate(data = map2(data, mean_stress,  quadrant_analysis)) %>%
    mutate(data = map(data, function(X) {
      X %>%
        mutate(Event = factor(Quadrant,
                              labels = c("Outward interaction",
                                         "Ejection",
                                         "Inward interaction",
                                         "Sweep"),
                              order = TRUE)) %>%
        mutate(Quadrant = factor(Quadrant,
                                 labels = c("I", "II", "III", "IV"),
                                 order = TRUE))
    })) %>%
    mutate(data = map(data, function(X) {
      X %>%
        select(-(x:z2), -w2, -w)
    }))

  dset %>%
    mutate(summarised_data = map2(data, total_time, function(X, total_time) {
      sampling_frequency = 0.005
      N = X %>% summarise(tmp = sum(!in_hole)) %>% pull(tmp)
      X %>%
        group_by(Quadrant, Event) %>%
        summarise(Stress = sum(Stress) / total_time * sampling_frequency,
                  Duration = sum(!in_hole) / total_time * sampling_frequency,
                  TKE = 0.5 * sum(up^2 + vp^2 + wp^2) / total_time * sampling_frequency,
                  Events = sum(!in_hole),
                  Proportion = Events / N)

    }))
}





# Timing two methods of calculating row sums.
N = 100000
XX = data_frame(x = rnorm(N),
           y = rnorm(N),
           r = runif(N)) %>%
  mutate(x = ifelse(r<0.2, NA, x),
         y = ifelse(r>0.8, NA, y)) %>%
  select(-r)


start_time <- Sys.time()
XX %>%
  mutate(xy = rowSums(cbind(x,y),na.rm=T)) %>%
  summarise_at(vars(xy), funs(mean = mean, sd = sd))
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
XX %>%
  rowwise() %>%
  mutate(xy = sum(x,y,na.rm=T)) %>%
  ungroup() %>%
  summarise_at(vars(xy), funs(mean = mean, sd = sd))
end_time <- Sys.time()
end_time - start_time








