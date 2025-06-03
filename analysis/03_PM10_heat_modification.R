################################################################################
### Research Objective 2: heat modification of PM10
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


### PM10 heat modification all
#----

# fit the model for both scenarios
mod4_heat <- glm(all ~ cb_pm10_2 + cb_temp + pm10_heat + dow + ns(date, df = 19*8),
                 data = data, family = quasipoisson)
mod4_heat_no <- glm(all ~ cb_pm10_2 + cb_temp + pm10_heat_no + dow + ns(date, df = 19*8),
                    data = data, family = quasipoisson)

# predict the association
pred4_heat <- crosspred(cb_pm10_2, mod4_heat, cen = 0, by = 1)
pred4_heat_no <- crosspred(cb_pm10_2, mod4_heat_no, cen = 0, by = 1)


# save plot
png("plots/model3_all.png", width = 1550, height = 900, res = 300)
par(mfrow = c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0)
)


plot(pred4_heat,
     "overall",
     ylim = c(.9,2.8),
     xlab = expression("PM10 ["*mu*"g/"*m^{3}*"]"),
     ylab = "relative risk",
     main = "",
     col = "brown1",
     ci.arg = list(col =alpha(colour = "brown1", .15)),
     cex.axis = 0.8,
     cex.lab = 0.8
)

lines(pred4_heat_no,           ## cumulative exposure
      "overall",
      col = "grey2",
      ci = "area",
      ci.arg = list(col = alpha(colour = "grey2", .15)),
      lwd = 2)
abline(v = 50, lty = "dashed")
legend(x = 52, y = 2.9, ncol = 1, legend = c("Model 3: heat days", "Model 3: non-heat days"), col = c("brown1", "grey2"),
       bty = "n", lwd=c(2,2), cex = 0.8)

value = round((  pred4_heat$matRRfit[50] -1 ) * 100, digits = 1)
text(49, 1.8 , labels = paste0(value  ,"%"), pos = 4, cex = 0.8)


dev.off()

#----


### PM10 heat modification subpopulations
#----

# empty data frame
model4_subgroups_df = data.frame(all = rep(NA,6),
                                 all_male = rep(NA,6),
                                 all_female = rep(NA,6),
                                 all_below64y = rep(NA,6),
                                 all_above64y = rep(NA,6),
                                 accid = rep(NA,6),
                                 renal = rep(NA,6),
                                 resp = rep(NA,6),
                                 cerebrov = rep(NA,6),
                                 cardio = rep(NA,6))

rownames(model4_subgroups_df) <- c("heat_Estimate", "heat_CI_Low", "heat_CI_High",
                                   "no_heat_Estimate", "no_heat_CI_Low", "no_heat_CI_High")

# loop through subpopulations
for(i in 1:ncol(model4_subgroups_df)){

  # extract subpopulation
  subpop_index = colnames(model4_subgroups_df[i])

  # formula for both scenarios
  heat_formula_subpop = as.formula(paste0(subpop_index, " ~ cb_pm10_2 + cb_temp + pm10_heat + dow + ns(date, df = ", 19*8, ")"))
  heat_no_formula_subpop = as.formula(paste0(subpop_index, " ~ cb_pm10_2 + cb_temp + pm10_heat_no + dow + ns(date, df = ", 19*8, ")"))

  # models
  heat_mod_subpop = glm(heat_formula_subpop, data = data, family = quasipoisson)
  heat_no_mod_subpop = glm(heat_no_formula_subpop, data = data, family = quasipoisson)

  # prediction
  heat_pred_subpop <- crosspred(cb_pm10_2, heat_mod_subpop, cen = 0, by = 1)
  heat_no_pred_subpop <- crosspred(cb_pm10_2, heat_no_mod_subpop, cen = 0, by = 1)

  # save values in table
  model4_subgroups_df[1,i] = heat_pred_subpop$matRRfit[50]
  model4_subgroups_df[2,i] = heat_pred_subpop$matRRlow[50]
  model4_subgroups_df[3,i] = heat_pred_subpop$matRRhigh[50]
  model4_subgroups_df[4,i] = heat_no_pred_subpop$matRRfit[50]
  model4_subgroups_df[5,i] = heat_no_pred_subpop$matRRlow[50]
  model4_subgroups_df[6,i] = heat_no_pred_subpop$matRRhigh[50]
}


# reproject for plotting
model4_subgroups_df <- model4_subgroups_df |>
  t() |>
  as.data.frame()

colnames(model4_subgroups_df) <- rep("X", 6)

model4_subgroups_df_long <- rbind(model4_subgroups_df[,1:3], model4_subgroups_df[,4:6])

model4_subgroups_df_long <- model4_subgroups_df_long |>
  mutate(scenario = c(rep("heat", 10), rep("no heat", 10)),
         scenario = factor(scenario, levels = c("no heat", "heat")),
         subpopulation = rep(c("all", "male", "female", "<64 years", ">=64 years", "accident", "renal", "respiratory", "cerebrovascular", "cardiovascular"),2),
         subpopulation = factor(subpopulation, levels = c("all", "male", "female", "<64 years", ">=64 years", "accident", "renal", "respiratory", "cerebrovascular", "cardiovascular")))

colnames(model4_subgroups_df_long)[1:3] <- c("Estimate", "CI_Low", "CI_High")



# save plot
png("plots/model3_subpop.png", width = 1900, height = 1400, res = 300)
par(mfrow = c(1,1),
    mar = c(3,3,.5,.5),
    mgp = c(1.8, .5, 0)
)

ggplot(model4_subgroups_df_long, aes(x = subpopulation, y = Estimate, color = scenario)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High),
                position = position_dodge(width = 0.6), width = 0.4) +
  geom_hline(yintercept = 1, linetype = 1, color = "gray2") +
  labs(
    x = "Subpopulation",
    y = "relative risk",
    title = ""
  ) +
  scale_color_manual(values = c("heat" = "brown1", "no heat" = "grey2"),
                     name = "Heat Mode") +
  theme_classic() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.02, 1),  # x and y in [0, 1] — upper left
    legend.justification = c(0, 1),   # anchor legend at top-left
    legend.background = element_rect(fill = alpha("white", 0), color = NA),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    # legend.title = element_text(size = 14),
    # legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    # axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  coord_cartesian(ylim = c(0.5, 2))

dev.off()

#----




