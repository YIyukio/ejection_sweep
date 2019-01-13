library(tidyverse)
library(ggplot2)
library(rstan)
library(loo)
library(broom)


fname = dir("./Modified_Data/", pattern = "all_data", full = T)
all_data = read_csv(fname) %>%
  mutate(Density = as.factor(Density),
         Density = factor(Density, levels = c("High", "Middle", "Low")))
quad_data = read_csv("Modified_Data/all_quadrant_data.csv")
Gross_data = all_data %>% filter(Hz != 0, type == "NP")

# x = Gross_data %>% select(type, Hz, Density, Slope, upstream_u_mean)
# x = x %>% mutate(Density = recode(Density, "Middle" = "Medium"))
#
# y = quad_data %>%
#   filter(hole == 0,
#          str_detect(Event, "Sweep|Ejection"),
#          Distance > 0) %>%
#   group_by(Hz, Density, Position, Event) %>%
#   summarise(TKE = mean(TKE),
#             Stress = mean(Stress),
#             Duration =        mean(Duration),
#             Stress_fraction = mean(Stress_fraction),
#             TKE_fraction =    mean(TKE_fraction))
#
# cdata = full_join(x,y, by = c("Hz", "Density"))
#
# cdata %>%
#   filter(str_detect(type, "NP")) %>%
#   ggplot(aes(x = TKE,
#              y = Slope, color = Event)) +
#   geom_point() +
#   geom_line()+
#   facet_grid(Position~Density)
#
#
# cdata %>%
#   filter(str_detect(type, "NP")) %>%
#   ggplot(aes(x = Stress,
#              y = Slope, color = Event)) +
#   geom_point() +
#   geom_line()+
#   facet_grid(Position~Density)
#
# cdata %>%
#   filter(str_detect(type, "NP")) %>%
#   ggplot(aes(x = Duration,
#              y = Slope, color = Event)) +
#   geom_point() +
#   geom_line()+
#   facet_grid(Position~Density)
#
#
# cdata %>%
#   filter(str_detect(type, "NP")) %>%
#   ggplot(aes(x = TKE_fraction,
#              y = Slope, color = Event)) +
#   geom_point() +
#   geom_line()+
#   facet_grid(Position~Density)
#
# cdata %>%
#   filter(str_detect(type, "NP")) %>%
#   ggplot(aes(x = Stress_fraction,
#              y = Slope, color = Event)) +
#   geom_point() +
#   geom_line()+
#   facet_grid(Position~Density)
#

X = all_data %>%
  group_by(type, Hz, Density) %>%
  summarise(upstream_u_mean = first(upstream_u_mean),
            Slope = first(Slope)) %>%
  filter(Hz!= 0)


X =
  X %>%
  mutate(speed = upstream_u_mean * 100) %>%
  select(speed, type, Slope) %>%
  group_by(type) %>% mutate(n = 1:n()) %>%
  spread(type, Slope)

N = X %>% nrow()
speed = X %>% pull(speed)
net = X %>% pull(NP)
resp = X %>% pull(RP)



M = 500
speed_predict = seq(0, max(speed), length = M)

stanmodel1 = stan_model(file = "nonlinear1.stan")
stanmodel2 = stan_model(file = "nonlinear2.stan")

# offsetrp, slope, offsetnp, pmax, alpha, beta
prior_log_mu    = c( 50, 1,  50,  50, log(5), log(1))
prior_log_sigma = c(100, 1, 100, 100,       1,     1) *
                  c(  5, 5,   5,   5, 1, 1)

rng_seed　= 2019
nchains = 5
iter = 20000
ncores = 5

stanout1 = sampling(stanmodel1,
                    chains = nchains,
                    iter = iter,
                    cores = ncores,
                    seed = rng_seed,
                    control = list(adapt_delta = 0.995,
                                   max_treedepth = 12))

stanout2 = sampling(stanmodel2,
                    chains = nchains,
                    iter = iter,
                    cores = ncores,
                    seed = rng_seed,
                    control = list(adapt_delta = 0.995,
                                   max_treedepth = 12))

#LOOを使って計算
log_lik1 = extract_log_lik(stanout1, merge_chains = FALSE)
log_lik2 = extract_log_lik(stanout2, merge_chains = FALSE)
r_eff_1  = relative_eff(exp(log_lik1), cores = ncores)
r_eff_2  = relative_eff(exp(log_lik2), cores = ncores)
loo1     = loo(log_lik1, r_eff = r_eff_1)
loo2     = loo(log_lik2, r_eff = r_eff_2)
loo1
loo2
print(compare(loo1, loo2), digits = 3)

pars1 = names(extract(stanout1))[c(1,2,3,4,5,7)]
pars2 = names(extract(stanout2))[c(1,2,3,4,5,8)]

print(stanout1, pars1)
print(stanout2, pars2)

preddata1 =
  data_frame(pars = c("yfit_gp", "ypred_np", "ypred_rp")) %>%
  mutate(data = map(pars, function(X) {
     tidy(stanout1,
         pars = X,
         conf.int = TRUE,
         conf.level = 0.95,
         conf.method = "HPDinterval") %>%
      mutate(speed = speed_predict)
  }))

preddata2 =
  data_frame(pars = c("yfit_gp", "ypred_np", "ypred_rp")) %>%
  mutate(data = map(pars, function(X) {
    tidy(stanout2,
         pars = X,
         conf.int = TRUE,
         conf.level = 0.95,
         conf.method = "HPDinterval") %>%
      mutate(speed = speed_predict)
  }))

preddata = bind_rows(preddata1 %>% mutate(model = "Saturating model"),
                     preddata2 %>% mutate(model = "Inhibition model")) %>% unnest()
preddata = preddata %>%
  mutate(pars = recode(pars,
         yfit_gp = "GP",
         ypred_rp = "RP",
         ypred_np = "NP"))

rawdata = all_data %>% filter(Hz != 0)
rawdata = rawdata %>% mutate(speed = upstream_u_mean * 100) %>%
  select(speed, Slope, pars = type)

ggplot() +
  geom_point(aes(x = speed, y = Slope, color = pars), rawdata) +
  geom_ribbon(aes(x = speed,
                  ymin = conf.low,
                  ymax = conf.high, fill = pars),
              data = preddata,
              alpha = 0.5) +
  geom_line(aes(x = speed, y = estimate, color = pars),
            data = preddata) + facet_wrap("model")
