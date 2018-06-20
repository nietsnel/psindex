
library(tidyverse)
library(rgenoud)
library(ggplot2)
library(GenSA)
library(devtools)




setwd
getwd()
.libPaths()
#needed up update package.

getwd()
setwd("..")


install("psindex")
library(psindex)

debugonce(psindex)


library(devtools)
library(psindex)



# install("lavaan_stability_index")





setwd("..")

library(lavaan)
library(tidyverse)
library(devtools)



model_original <- '
# measurement model
ind60  =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60  + dem60
# residual correlations
y1 ~~ y5
y2 ~~ y4 + y6
y3 ~~ y7
y4 ~~ y8
y6 ~~ y8
'

install("psindex")
library(psindex)


ps_index(model = model_original, data_set = PoliticalDemocracy,
         RMSEA_pert = .01, plot_fpe = TRUE, frac_plot = .3)




ps_index(model = model_original, data_set = PoliticalDemocracy,
         iterations_bin = 5000, GENANNEAL_steps = 100, RMSEA_pert = .03, plot_fpe = FALSE, frac_plot = .3)






View(fpe_wide[1:32,1:4])


ps_index(sem_model = model_original, data_set = PoliticalDemocracy,
         iterations_bin = 5000, GENANNEAL_steps = 100, RMSEA_pert = .04, frac_plot = .3, plot_fpe = TRUE)



dt_speed[1:20, 1:4]

ggplot() +
  geom_point(data=fpe_long, aes(x=variable, y=estimate), stat="identity") +
  geom_point(data = mle_est_long, aes(x=variable, y=estimate, fill="#FF3300", colour = "#FF3300")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
  coord_cartesian(ylim = c(-1, 7.5))





psindex::ps_index(sem_model = model_original, data_set = PoliticalDemocracy, iterations_bin = 5000)



psindex

traceback()


M = matrix(0, nrow = length(x) + 1, ncol = 500)
iters_assign <- 1
DT = as.data.table(M)



# devtools::document()
traceback(lavaan::psindex(data = PoliticalDemocracy, sem_model=model_original, iterations_bin=10000, bandwidth = .02))

model_original <- '
# measurement model
ind60  =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60  + dem60
# residual correlations
y1 ~~ y5
y2 ~~ y4 + y6
y3 ~~ y7
y4 ~~ y8
y6 ~~ y8
'


# 
# data <- data(PoliticalDemocracy, package="lavaan")
# 
# fit <- psindex::sem(model_original, data=PoliticalDemocracy)
# summary(fit, standardized=TRUE)
# rm(fit)
# 
# psindex::sem(model_original, data=PoliticalDemocracy)
# 
# 
# library(lavaan)
# library(lavaan.survey)

library(tidyverse)

######################################################
######################################################
#Graph results

results_array_subset <- results_array[1:3500,]




results_array_subset <- results_array
data_iterations <- as.data.frame(as.matrix(results_array_subset[,]))

data_iterations <- data_iterations %>%
  drop_na()
names <- rbind(names, "discrepancy fx", "count")

names(data_iterations) <- names$new

results_main_long <- data_iterations %>%
  # arrange(`discrepancy fx`)
  # mutate(mle = ifelse(`discrepancy fx`))
  gather(key=variable, value=estimate, -count, -`discrepancy fx`, (1:(dim(data_iterations)[2]-2))) %>%
  arrange(variable)

mle_dat <- results_main_long %>%
  filter(count == 1)
names(results_array_mle)<- names$new[1:(dim(results_array_mle)[2])]

results_array_mle_long <- results_array_mle %>%
  # arrange(`discrepancy fx`)
  # mutate(mle = ifelse(`discrepancy fx`))
  gather(key=variable, value=estimate) %>%
  arrange(variable)

fit_meas <- round(as.data.frame(fitMeasures(fit.initial, c("cfi","rmsea","srmr"))),3)

ggplot() +
  geom_point(data=results_main_long, aes(x=variable, y=estimate), stat="identity") +
  geom_point(data = results_array_mle_long, aes(x=variable, y=estimate, fill="#FF3300", colour = "#FF3300")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
  coord_cartesian(ylim = c(-1, 7.5)) +
  # labs(title="Political Democracy SEM Iteration History", subtitle=paste0("RMSEA perturbation = ", bandwidth, " ", "(CFI = ", fit_meas[1,1], " RMSEA = ", fit_meas[2,1], ")"), x="Parameter", y="Estimate")
  # labs(title="Political Democracy SEM Iteration History", subtitle=paste0("F perturbation  = ", RMSEA_pert, "% ", "(CFI = ", fit_meas[1,1], " RMSEA = ", fit_meas[2,1], ")"), x="Parameter", y="Estimate")




ggsave(filename='C:/Users/jordan/Desktop/MMM project/rmsea_.05.png', height = 4, width = 6.5, dpi=300)




# ggsave(filename='C:/Users/jprendez/Desktop/edms_research_day2018/polidem_2perc.png', height = 4, width = 6.5, dpi=300)


##J.Harring dataset
library(semPlot)


ggplot() +
  geom_point(data=results_main_long, aes(x=variable, y=estimate), stat="identity") +
  geom_point(data = results_array_mle_long, aes(x=variable, y=estimate, fill="#FF3300", colour = "#FF3300")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
  coord_cartesian(ylim = c(-30, 30)) +
  # labs(title="Political Democracy SEM Iteration History", subtitle=paste0("RMSEA perturbation = ", bandwidth, " ", "(CFI = ", fit_meas[1,1], " RMSEA = ", fit_meas[2,1], ")"), x="Parameter", y="Estimate")
  labs(title="Read dataset", subtitle=paste0("F perturbation  = ", bandwidth, "% ", "(CFI = ", fit_meas[1,1], " RMSEA = ", fit_meas[2,1], ")"), x="Parameter", y="Estimate")


ggsave(filename='C:/Users/jprendez/Desktop/edms_research_day2018/polidem_2perc.png', height = 4, width = 6.5, dpi=300)



semPaths(fit.initial)
ggsave(filename='C:/Users/jprendez/Desktop/edms_research_day2018/j.harring5.png', height = 4, width = 6.5, dpi=300)



C:\Users\jprendez\Desktop\edms_research_day2018


# a linear growth model with a time-varying covariate
model <- '
# intercept and slope with fixed coefficients
i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
# regressions
i ~ x1 + x2
s ~ x1 + x2
# time-varying covariates
t1 ~ c1
t2 ~ c2
t3 ~ c3
t4 ~ c4
'
fit <- growth(model, data = Demo.growth)
summary(fit)

#
# ggsave(filename='C:/Users/jordan/Desktop/pictures for prendex/iter_hist.png', height = 4, width = 8, dpi=300)


#
# ggsave(filename='C:/Users/jordan/Desktop/pictures for prendex/iter_hist_last2.png', height = 4, width = 8, dpi=300)
#
#
#
# ggsave(filename='C:/Users/jprendez/Documents/usm/financial_aid_graphs/cost_of_attendance_vs_peers.tiff', height = 4, width = 6.5, dpi=300)
#
#




#
#
