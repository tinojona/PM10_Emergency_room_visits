################################################################################
### Sensitivity Analysis
################################################################################



### Data
#----

# libraries
library(tidyverse); library(lubridate); library(splines); library(dlnm); library(MASS);library(stringr)

# clean env
rm(list = ls())


# load processed data
data <- read.csv("data/emerg_ap_CLEAN.csv") |> mutate(date = as.Date(date))

# monthly aggregates of main variables
data_monthly <- data |>
  group_by(month) |>
  summarise(across(c(all, tmean, pm10), mean, na.rm = TRUE))

data_monthly_long <- data_monthly |>
  pivot_longer(cols = c(pm10, tmean, all),
               names_to = "variable", values_to = "value")

# subpopulations
subpop <- c("all", "all_male", "all_female", "all_below64y", "all_above64y",
            "accid", "renal", "resp", "cerebrov", "cardio")

#----



### Crossbasis
#----

# create the PM10 cross-basis
cb_pm10_2 <- crossbasis(data$pm10_MA_2, lag = 0, argvar = list(fun = "lin"))

# create the temperature crossbasis
cb_temp <- crossbasis(data$tmean,
                      lag=21,
                      argvar=list(fun="ns", knots = quantile(data$tmean, c(.1,.75,.9), na.rm=TRUE)),
                      arglag=list(fun="ns", knots = logknots(21,3)))

# define interaction term for heat and no heat
pm10_heat <- cb_pm10_2 * data$heat_90
pm10_heat_no <- cb_pm10_2 * data$heat_90_no

#----



### PM10 moving average length (5 days)
#----

# empty data frame
model2_subgroups_df = data.frame(all = rep(NA,3),
                                 all_male = rep(NA,3),
                                 all_female = rep(NA,3),
                                 all_below64y = rep(NA,3),
                                 all_above64y = rep(NA,3),
                                 accid = rep(NA,3),
                                 renal = rep(NA,3),
                                 resp = rep(NA,3),
                                 cerebrov = rep(NA,3),
                                 cardio = rep(NA,3))

rownames(model2_subgroups_df) <- c("Estimate", "CI_Low", "CI_High")

# loop through all subpopulations
for(i in 1:ncol(model2_subgroups_df)){

  # extract subpopulation
  subpop_index = colnames(model2_subgroups_df[i])

  # formula
  formula_subpop = as.formula(paste0(subpop_index, " ~ cb_pm10_2 + cb_temp + dow + ns(date, df = ", 19*8, ")"))

  # fir model
  mod_subpop = glm(formula_subpop, data = data, family = quasipoisson)

  # predict
  pred_subpop <- crosspred(cb_pm10_2, mod_subpop, cen = 0, by = 1)

  # save risk and confidence interval
  model2_subgroups_df[1,i] = pred_subpop$matRRfit[50]
  model2_subgroups_df[2,i] = pred_subpop$matRRlow[50]
  model2_subgroups_df[3,i] = pred_subpop$matRRhigh[50]
}

# reproject for plotting
model2_subgroups_df_long <- model2_subgroups_df |>
  t() |>
  as.data.frame() |>
  rownames_to_column("Model") |>
  rename(
    Estimate = Estimate,
    CI_Low = CI_Low,
    CI_High = CI_High
  ) |>
  mutate(Model = c("all", "male", "female", "<64 years", ">64 years", "accident",
                   "renal", "respiratory", "cerebrovascular", "cardiovascular"),
         Model <- factor(Model, levels = c("all", "male", "female", "<64 years", ">64 years", "accident",
                                           "renal", "respiratory", "cerebrovascular", "cardiovascular")))

# save plot
png("plots/X_SENS_model2_subpop_MA5.png", width = 1200, height = 800, res = 100)

ggplot(model2_subgroups_df_long, aes(x = Model, y = Estimate)) +
  geom_point(size = 3, color = "brown1") +
  # lims(y = c(0.48, 1.65)) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, color = "brown1") +
  geom_segment(aes(x = Model, xend = Model, y = CI_Low, yend = CI_High), color = "brown1") +
  labs(title = "",
       y = "relative risk",
       x = "subpopulation") +
  geom_hline(yintercept = 1, color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        # axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  coord_cartesian(ylim = c(0.48, 1.65))

dev.off()

#----



### Temperature lag period 5 days
#----

# define new crossbasis
cb_temp <- crossbasis(data$tmean,
                      lag=5,
                      argvar=list(fun="ns", knots = quantile(data$tmean, c(.1,.75,.9), na.rm=TRUE)),
                      arglag=list(fun="ns"))


# fit Model 2 (including temperature confounding)
mod2 <- glm(all ~ cb_pm10_2 + cb_temp + dow + ns(date, df = 19*8), data = data, family = quasipoisson)

# predict Model 2
pred2 <- crosspred(cb_pm10_2, mod2, cen = 0, by = 1)


# save model
png("plots/X_SENS_model2_all_TempLag5days.png", width = 1800, height = 1500, res = 300)
par(mfrow = c(1,1),
    mar = c(4,4,.5,1),
    mgp = c(2.5, .8, 0)
)

plot(pred2,
     "overall",
     # ylim = c(1,1.5),
     col = "brown1",
     ci.arg = list(col = alpha(colour = "brown1", .15)),
     xlab = "PM10 exposure",
     lwd = 2,
     ylab = "relative risk",
     main = "",
     cex.axis = 2,
     cex.lab = 2
)
abline(v = 50, lty = "dashed")


value = round((  pred2$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.4 , labels = paste0(value  ,"%"), pos = 4, cex = 2)

dev.off()


#----



### Heat percentile threshold
#----

# define interaction term for heat and no heat
pm10_heat <- cb_pm10_2 * data$heat_95
pm10_heat_no <- cb_pm10_2 * data$heat_95_no

# fit mode
mod4_heat <- glm(all ~ cb_pm10_2 + cb_temp + pm10_heat + dow + ns(date, df = 19*8),
                 data = data, family = quasipoisson)
mod4_heat_no <- glm(all ~ cb_pm10_2 + cb_temp + pm10_heat_no + dow + ns(date, df = 19*8),
                    data = data, family = quasipoisson)

# predict the association
pred4_heat <- crosspred(cb_pm10_2, mod4_heat, cen = 0, by = 1)
pred4_heat_no <- crosspred(cb_pm10_2, mod4_heat_no, cen = 0, by = 1)


png("plots/X_SENS_model3_all_95p.png", width = 1550, height = 1200, res = 300)
par(mfrow = c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0)
)


plot(pred4_heat,
     "overall",
     # ylim = c(1,3),
     xlab = "PM10 exposure",
     ylab = "relative risk",
     main = "",
     col = "brown1",
     ci.arg = list(col =alpha(colour = "brown1", .15))
)
lines(pred4_heat_no,           ## cumulative exposure
      "overall",
      col = "grey2",
      ci = "area",
      ci.arg = list(col = alpha(colour = "grey2", .15)),
      lwd = 2)
abline(v = 50, lty = "dashed")
legend("top", ncol = 1, legend = c("Model 4: heat days", "Model 4: non-heat days"), col = c("brown1", "grey2"),
       bty = "n", lwd=c(2,2), cex = 0.7)

value = round((  pred4_heat$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.55 , labels = paste0(value  ,"%"), pos = 4, cex = 1)


dev.off()

#----



### Temperature cummulative relative risk
#----

# fit Model 1
mod1 <- glm(all ~ cb_temp + dow + ns(date, df = 19*8), data = data, family = quasipoisson)

# predict Model 1
pred1 <- crosspred(cb_temp, mod1, cen = 20, by = 1)


png("plots/X_temp_cummulative_RR.png", width = 1500, height = 1400, res = 300)
par(mfrow = c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0)
)
plot(pred1,
     "overall",
     # ylim = c(1,1.5),
     col = "purple3",
     ci.arg = list(col = alpha(colour = "purple3", .15)),
     xlab = "temperature [°C]",
     lwd = 2,
     ylab = "relative risk",
     main = ""
)


dev.off()


# 3d plot
png("plots/X_temp_3d.png", width = 1000, height = 800, res = 100)
plot(pred1,
     "3d",
     # ylim = c(1,1.5),
     col = "purple3",
     # ci.arg = list(col = alpha(colour = "purple3", .15)),
     xlab = "temperature [°C]",
     lwd = 2,
     ylab = "lag",
     zlab = "relative risk",
     zlim = c(0.8,1.1),
     main = "",
     theta = 45,  # horizontal rotation angle (default: 40)
     phi = 20     # vertical angle (default: 30)
)
dev.off()

#----









