## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.height = 8,
  fig.width = 8,
  comment = "#>"
)

## ----duffy--------------------------------------------------------------------

library(multifunc)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
theme_set(theme_bw(base_size = 14))

data("duffy_2003")
duffyAllVars <- qw(grazer_mass,wkall_chla,tot_algae_mass,
                  Zost_final_mass,sessile_invert_mass,sediment_C)

#re-normalize so that everything is on the same 
#sign-scale (e.g. the maximum level of a function is the "best" function)
#and the dataset is cleaner
duffy <- duffy_2003 %>%
 dplyr::select(id, treatment, 
               diversity, all_of(duffyAllVars)) %>%
 dplyr::mutate(wkall_chla = -1*wkall_chla + 
                 max(wkall_chla, na.rm=T),
               tot_algae_mass = -1*tot_algae_mass + 
                 max(tot_algae_mass, na.rm=T)) 

# Plot functions we want
duffy %>% 
  select(treatment, diversity, all_of(duffyAllVars)) %>%
  pivot_longer(all_of(duffyAllVars)) %>%
  ggplot(aes(x = diversity, y = value, 
             color = treatment)) +
  stat_summary(fun.data = mean_se) +
  stat_smooth(method = "lm", aes(group = name)) +
  facet_wrap(vars(name), scale = "free_y") +
  scale_color_brewer(palette = "Set3") +
  labs(y = "Level of Function", x = "Species Richness",
       color = "Treatment")

## ----std----------------------------------------------------------------------
#first, mean multifunctionality
duffy <- duffy %>%
 cbind(getStdAndMeanFunctions(duffy, duffyAllVars)) %>%
 dplyr::rename(richness=diversity)



duffy <-duffy %>%
  getFuncsMaxed(duffyAllVars,
                threshmin=0.8, threshmax=0.8)

## -----------------------------------------------------------------------------
duffyAllVars.std <- paste0(duffyAllVars, ".std")

#now effective number of functions
duffy <- duffy %>%
  mutate(n_eff_func_1 = eff_num_func(., duffyAllVars.std, q = 1),
         n_eff_func_2 = eff_num_func(., duffyAllVars.std, q = 2),
         
         mf_eff_1 = n_eff_func_1 * meanFunction,
         mf_eff_2 = n_eff_func_2 * meanFunction
         )

## -----------------------------------------------------------------------------
duffy %>%
  select(id, treatment, richness,
         meanFunction, funcMaxed,
         n_eff_func_1, n_eff_func_2,
         mf_eff_1, mf_eff_2) %>%
  pivot_longer(cols = c(meanFunction:mf_eff_2)) %>%
  mutate(name = fct_inorder(name)) %>%
  
  #now a plot
  ggplot(aes(x = richness, y = value,
             color = treatment)) +
  stat_summary(fun.data = mean_se) +
  stat_smooth(method = "lm", color = "black") +
  facet_wrap(vars(name), ncol = 2, scale = "free_y") +
  scale_color_brewer(palette = "Set3") 

## -----------------------------------------------------------------------------
D <- cor_dist(duffy %>% select(all_of(duffyAllVars.std)))

## ----dmean--------------------------------------------------------------------
tau <- dmean(duffy %>% select(all_of(duffyAllVars.std)),
             D)

tau

## ----mf_d---------------------------------------------------------------------
duffy <- duffy %>%
  mutate(mf_eff_1_cor = getMF_eff(., duffyAllVars, q = 1, 
                                  D = D, tau = tau),
         mf_eff_2_cor = getMF_eff(., duffyAllVars, q = 2,
                                  D = D, tau = tau)
         )

## -----------------------------------------------------------------------------
duffy %>%
  select(id, treatment, richness,
         mf_eff_1, mf_eff_2,
         mf_eff_1_cor, mf_eff_2_cor) %>%
  pivot_longer(cols = c(mf_eff_1:mf_eff_2_cor)) %>%
  mutate(name = fct_inorder(name)) %>%
  
  #now a plot
  ggplot(aes(x = richness, y = value,
             color = treatment)) +
  stat_summary(fun.data = mean_se) +
  stat_smooth(method = "lm", color = "black") +
  facet_wrap(vars(name), ncol = 2) +
  scale_color_brewer(palette = "Set3") 


