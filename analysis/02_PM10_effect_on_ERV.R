################################################################################
### Research Objective 1: the effect of PM10 with/without confounders
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

# create the NO2 cross-basis
cb_no2_2 <- crossbasis(data$no2_MA_2, lag = 0, argvar = list(fun = "lin"))

# create the O3 cross-basis
cb_o3_2 <- crossbasis(data$o3_MA_2, lag = 0, argvar = list(fun = "lin"))

# create the RH cross-basis
cb_relhum_2 <- crossbasis(data$relhum_MA_2, lag = 0, argvar = list(fun = "lin"))

#----



### PM10 Model 1 and 2 for all
#----

# The association between PM10 and emergency room visits for all visits in Model 1
# (without accounting for temperature confounding) and Model 2 (with accounting
# for temperature confounding)

# fit Model 1
mod1 <- glm(all ~ cb_pm10_2 + dow + ns(date, df = 19*8), data = data, family = quasipoisson)

# predict Model 1
pred1 <- crosspred(cb_pm10_2, mod1, cen = 0, by = 1)

# fit Model 2 (including temperature confounding)
mod2 <- glm(all ~ cb_pm10_2 + cb_temp + dow + ns(date, df = 19*8), data = data, family = quasipoisson)

# predict Model 2
pred2 <- crosspred(cb_pm10_2, mod2, cen = 0, by = 1)


# save plot
png("plots/model1_model2_all.png", width = 2000, height = 1200, res = 300)
par(mfrow = c(1,2),
    mar = c(3,3,2,.5),
    mgp = c(1.8, .5, 0)
)
plot(pred1,
     "overall",
     ylim = c(1,1.5),
     col = "grey2",
     ci.arg = list(col = alpha(colour = "grey2", .15)),
     xlab = expression("PM10 ["*mu*"g/"*m^{3}*"]"),
     lwd = 2,
     ylab = "relative risk",
     main = "Model 1"
)
abline(v = 50, lty = "dashed")

# legend("top", ncol = 1, legend = c("Model 1: PM10"), col = c("grey2"),
#        bty = "n", lwd=c(2,2), cex = 1)

value = round((  pred1$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.3 , labels = paste0(value  ,"%"), pos = 4, cex = 1)



plot(pred2,
     "overall",
     ylim = c(1,1.5),
     col = "brown1",
     ci.arg = list(col = alpha(colour = "brown1", .15)),
     xlab = expression("PM10 ["*mu*"g/"*m^{3}*"]"),
     lwd = 2,
     ylab = "relative risk",
     main = "Model 2"
)
abline(v = 50, lty = "dashed")

# legend("top", ncol = 1, legend = c("Model 2: PM10 + Temperature"), col = c("brown1"),
#        bty = "n", lwd=c(2,2), cex = 1)

value = round((  pred2$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.3 , labels = paste0(value  ,"%"), pos = 4, cex = 1)

dev.off()

#----



### PM10 Model 2: subpopulations
#----

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

# loop through subpopulations
for(i in 1:ncol(model2_subgroups_df)){

  # extract subpopulation
  subpop_index = colnames(model2_subgroups_df[i])

  # define formual
  formula_subpop = as.formula(paste0(subpop_index, " ~ cb_pm10_2 + cb_temp + dow + ns(date, df = ", 19*8, ")"))

  # fit model
  mod_subpop = glm(formula_subpop, data = data, family = quasipoisson)

  # predict on the model
  pred_subpop <- crosspred(cb_pm10_2, mod_subpop, cen = 0, by = 1)

  # extract and save estimate and confidence interval
  model2_subgroups_df[1,i] = pred_subpop$matRRfit[50]
  model2_subgroups_df[2,i] = pred_subpop$matRRlow[50]
  model2_subgroups_df[3,i] = pred_subpop$matRRhigh[50]

}

# convert to long format
model2_subgroups_df_long <- model2_subgroups_df |>
  t() |>
  as.data.frame() |>
  rownames_to_column("Model") |>
  rename(
    Estimate = Estimate,
    CI_Low = CI_Low,
    CI_High = CI_High
  ) |>
  mutate(
    # assign prpoper names
    Model = c("all", "male", "female", "<64 years", ">=64 years", "accident",
              "renal", "respiratory", "cerebrovascular", "cardiovascular"),
    Model = factor(Model, levels = c("all", "male", "female", "<64 years", ">=64 years", "accident",
                                     "renal", "respiratory", "cerebrovascular", "cardiovascular")))

# save plot
png("plots/model2_subpop.png", width = 1900, height = 1400, res = 300)

ggplot(model2_subgroups_df_long, aes(x = Model, y = Estimate)) +
  geom_point(size = 3, color = "brown1") +
  lims(y = c(0.48, 1.65)) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.2, color = "brown1") +
  geom_segment(aes(x = Model, xend = Model, y = CI_Low, yend = CI_High), color = "brown1") +
  labs(title = "",
       y = "relative risk",
       x = "Subpopulation") +
  geom_hline(yintercept = 1, color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        # axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

dev.off()

#----



### PM10 Model 2 + additional confounders
#----

# Model2 + NO2
mod3_no2 <- glm(all ~ cb_pm10_2 + cb_temp + cb_no2_2 + dow + ns(date, df = 19*8),
                data = data, family = quasipoisson)

# Model2 + O3
mod3_o3 <- glm(all ~ cb_pm10_2 + cb_temp + cb_o3_2 + dow + ns(date, df = 19*8),
               data = data, family = quasipoisson)

# Model2 + RH
mod3_relhum <- glm(all ~ cb_pm10_2 + cb_temp + cb_relhum_2 + dow + ns(date, df = 19*8),
                   data = data, family = quasipoisson)

# predict on the models
pred3_no2 <- crosspred(cb_pm10_2, mod3_no2, cen = 0, by = 1)
pred3_o3 <- crosspred(cb_pm10_2, mod3_o3, cen = 0, by = 1)
pred3_relhum <- crosspred(cb_pm10_2, mod3_relhum, cen = 0, by = 1)



# save plot
png("plots/model2_plus_additional_confounders.png", width = 2200, height = 800, res = 300)
par(mfrow = c(1,3),
    mar = c(3,3,3,.5),
    mgp = c(1.8, .5, 0)
)

plot(pred2,
     "overall",
     ylim = c(.9,1.5),
     col = "grey2",
     ci.arg = list(col = alpha(colour = "grey2", .15)),
     xlab = expression("PM10 ["*mu*"g/"*m^{3}*"]"),
     lwd = 2,
     ylab = "relative risk",
     main = "Model 3: O3"
)
lines(pred3_o3,
      "overall",
      col = "gold",
      ci = "area",
      ci.arg = list(col = alpha(colour = "gold", .15)),
      lwd = 2,
)
abline(v = 50, lty = "dashed")
value = round((  pred3_o3$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.4 , labels = paste0(value  ,"%"), pos = 4, cex = 1)

# legend("top", ncol = 1, legend = c("Model 3: PM10 + Temperature + O3"), col = c("gold"),
#        bty = "n", lwd=c(2,2), cex = 1)


plot(pred2,
     "overall",
     ylim = c(.9,1.5),
     col = "grey2",
     ci.arg = list(col = alpha(colour = "grey2", .15)),
     xlab = expression("PM10 ["*mu*"g/"*m^{3}*"]"),
     lwd = 2,
     ylab = "relative risk",
     main = "Model 3: RH"
)
lines(pred3_relhum,
      "overall",
      col = "steelblue",
      ci = "area",
      ci.arg = list(col = alpha(colour = "steelblue", .15)),
      lwd = 2,
      lty = 1
)
abline(v = 50, lty = "dashed")
value = round((  pred3_relhum$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.4 , labels = paste0(value  ,"%"), pos = 4, cex = 1)

# legend("top", ncol = 1, legend = c("Model 3: PM10 + Temperature+ RH"), col = c("steelblue"),
#        bty = "n", lwd=c(2,2), cex = 1)


plot(pred2,
     "overall",
     ylim = c(.9,1.5),
     col = "grey2",
     ci.arg = list(col = alpha(colour = "grey2", .15)),
     xlab = expression("PM10 ["*mu*"g/"*m^{3}*"]"),
     lwd = 2,
     ylab = "relative risk",
     main = "Model 3: NO2"
)
lines(pred3_no2,
      "overall",
      col = "green3",
      ci = "area",
      ci.arg = list(col = alpha(colour = "green3", .15)),
      lwd = 2,
      lty = 1
)
abline(v = 50, lty = "dashed")
value = round((  pred3_no2$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.4 , labels = paste0(value  ,"%"), pos = 4, cex = 1)

# legend("top", ncol = 1, legend = c("Model 3: PM10 + Temperature + NO2"), col = c("green3"),
#        bty = "n", lwd=c(2,2), cex = 1)

dev.off()











