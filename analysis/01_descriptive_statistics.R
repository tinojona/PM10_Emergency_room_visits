################################################################################
### Descriptive Statistics
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


#----


### Seasonal Cycle
#----

# Average seasonal cycle of PM10, temperature and ERV

# scale factor to align ERV with other variables visually
scale_factor <- max(data_monthly$tmean, data_monthly$pm10, na.rm = TRUE) / max(data_monthly$all, na.rm = TRUE)

# save plot
png("plots/seasonal_cycle.png", width = 1500, height = 1200, res = 300)

ggplot() +
  geom_line(data = data_monthly, aes(x = month, y = pm10, color = "PM10"), lwd = 1.2) +
  geom_line(data = data_monthly, aes(x = month, y = tmean, color = "TEMP"), lwd = 1.2) +
  geom_line(data = data_monthly, aes(x = month, y = all * scale_factor, color = "ERV"), lwd = 1.2) +
  scale_color_manual(name = "Variable",
                     values = c("PM10" = "#21908C", "TEMP" = "#440154", "ERV" = "#FDE725")) +
  scale_y_continuous(
    name = expression("Temperature [°C] and PM10 ["*mu*"g/"*m^{3}*"]"),
    sec.axis = sec_axis(~./scale_factor, name = "ERV")
  ) +

  theme_classic()+
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_x_continuous(breaks = 1:12, labels = 1:12, name = "Month")

dev.off()

#----



### Time series
#----

# save plot
png("plots/timeseries.png", width = 1600, height = 1400, res = 300)

# margins
par(mfrow = c(3, 1),
    mar = c(0, 4, 0.5, 0.5),  # Top two plots have zero bottom margin
    oma = c(4, 0, 0, 0),
    mgp = c(2, .5, 0))      # Outer bottom margin for x-axis label

# foehn plot
plot(data$date, data$pm10,
     type = "l",
     xlab = "",
     col = "#21908C",
     ylab = expression("PM10 ["*mu*"g/"*m^{3}*"]"),
     xaxt = "n",
     lwd = 2,
     ylim = c(-2, 200),
     cex.axis = 1.35,
     cex.lab = 1.35
)

text(data$date[nrow(data)],150, labels = "(a)", pos = 2, cex = 1.35)

# temp plot
plot(data$date, data$tmean,
     type = "l",
     xlab = "",
     col = "#440154",
     ylab = "Temperature [°C]",
     xaxt = "n",
     lwd = 2,
     ylim = c(-15,35),
     cex.axis = 1.35,
     cex.lab = 1.35
)

text(data$date[nrow(data)], 32, labels = "(b)", pos = 2, cex = 1.35)

# hosp plot
plot(data$date, data$all,
     type = "l",
     xlab = "",
     ylab = "ERV",
     lwd = 2,
     col = "#FDE725",
     ylim = c(-2,17),
     cex.axis = 1.35,
     cex.lab = 1.35
)

mtext("Time", side = 1, line = 2, outer = TRUE, cex = 1)

text(data$date[nrow(data)], 16, labels = "(c)", pos = 2, cex = 1.35)

dev.off()

#----



### PM10: 2-day moving average
#----

# save plot
png("plots/2-dayMA.png", width = 1400, height = 1000, res = 300)

ggplot(data = data[250:300,], aes(x = date)) +
  geom_line(aes(y = pm10, linetype = "Observation"), color = "black") +
  geom_line(aes(y = pm10_MA_2, linetype = "Moving average"), color = "black") +
  scale_linetype_manual(name = "PM10", values = c("Observation" = "solid", "Moving average" = "dashed")) +
  labs(x = "date", y = expression("PM10 ["*mu*"g/"*m^{3}*"]")) +
  theme_classic() +
  theme(legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1))

dev.off()

#----

