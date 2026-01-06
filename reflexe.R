##########################################################
# Titre : Analyse des réflexes oculomoteurs chez le marlin
# Auteur : Romain Winsback
# Date : 22-12-2025
##########################################################

# packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(circular)

# importation dataset
setwd("~/M2/S1/Analyses spatio-temporelles/Vignon/DM")
ocu <- read.csv2("data.csv")

# mise en forme dataset
str(ocu)
ocu$theta <- as.numeric(ocu$theta)
ocu$Y <- as.numeric(ocu$Y)

ocu <- ocu %>%
  mutate(
    Origin = factor(Origin),
    ID = factor(ID),
    theta_rad = theta * pi/180,  # valeur de theta en radian
    sin_t = sin(theta_rad),     # sinus de theta
    cos_t = sin(theta_rad)    # cosinus de theta
  )

# visualisation
ggplot(ocu, aes(theta, Y, color=Origin)) +
  geom_point(alpha=0.25, size=0.7) +
  geom_smooth(se=FALSE, method="loess", span=0.2) +
  scale_x_continuous(breaks=seq(0,360,60), limits=c(0,360)) +
  theme_bw()

# modèle mixte (harmonique ?)
mod1 <- lmer(Y~Origin*(sin_t + cos_t)+(1|ID), data=ocu)
summary(mod1)
anova(mod1)

mod0 <- lmer(Y~Origin+sin_t + cos_t+(1|ID), data=ocu)
anova(mod0,mod1)


grid <- expand.grid(
  theta = seq(0, 360, by = 1),
  Origin = levels(ocu$Origin)
) %>%
  mutate(
    theta_rad = theta*pi/180,
    sin_t = sin(theta_rad),
    cos_t = cos(theta_rad)
  )

grid$pred <- predict(mod1, newdata = grid, re.form = NA)

ggplot(grid, aes(theta, pred, color = Origin)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(y = "Y prédit (effets fixes)", x = "Angle θ (°)")


