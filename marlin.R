# ============================================================================
# Ocular reflex analysis - Marlins
#
# Author : Romain Winsback
# Date : 27-12-2026
#
# Data not publicly available.
# Dataset provided by the course instructor for educational purposes.
# ============================================================================

set.seed(123)

# ===========================
# Packages
# ===========================
library(dplyr)
library(lme4)
library(circular)
library(ggplot2)
library(mgcv)

# ===========================
# Load data
# ===========================
getwd()
setwd("~/M2/S1/Analyses spatio-temporelles/Vignon/DM")
ocu <- read.csv2("data.csv")

# ===========================
# Data preparation
# ===========================
str(ocu)
summary(ocu)

ocu$theta <- as.numeric(ocu$theta)
ocu$Y <- as.numeric(ocu$Y)
ocu$ID <- as.factor(ocu$ID)
ocu$Origin <- as.factor(ocu$Origin)

ocu <- ocu %>%
  mutate(
    theta_rad = theta * pi/180,  # conversion of angular variable to radians
    sin_t = sin(theta_rad),   
    cos_t = cos(theta_rad)   
  )

# ==========================
# Exploratory data analysis
# ==========================

# 1) Raw signal
ggplot(ocu, aes(theta, Y, color=Origin)) +
  geom_point(alpha=0.25, size=0.7) +
  geom_smooth(se=FALSE, method="loess", span=0.2) +
  scale_x_continuous(breaks=seq(0,360,60), limits=c(0,360)) +
  theme_bw()

# or ...
ggplot(ocu,aes(x=theta_rad,y=Y,color=Origin))+
  geom_point(alpha=0.3) +
  theme_bw()

# 2) Inter-individual variability (some profiles)
sample_ids <- sample(unique(ocu$ID),6)

ggplot(filter(ocu,ID %in% sample_ids),
       aes(x=theta_rad,y=Y,group=ID)) +
  geom_line(colour="black",alpha=0.8) +
  facet_wrap(~ID) +
  scale_x_continuous(
    breaks = seq(0,2*pi,by=pi),
    labels=c("0","π","2π")
  ) +
  labs(x = "Roll angle (radians)",
       y = "Ocular muscle tonus (µV)")

# 3) Visual comparison of average ocular muscle tone
ocu_ind <- ocu %>%
  group_by(ID, Origin) %>%
  summarise(
    mean_Y = mean(Y, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(ocu_ind, aes(x = Origin, y = mean_Y, fill = Origin)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(
    x = "Origin",
    y = "Average muscular tone"
  ) +
  theme_bw()  +
  theme(
    legend.position = "none",
    axis.title=element_text(size=20,face="bold"),
    axis.text=element_text(size=20),
    legend.title=element_text(size=25),
  )


# =============================
# Linear mixed models (LMM)
# =============================

mod1 <- lmer(Y~Origin*(sin_t + cos_t)+(1|ID), data=ocu)
summary(mod1)

mod0 <- lmer(Y~Origin+sin_t + cos_t+(1|ID), data=ocu)
anova(mod0,mod1)

# Model diagnostics 
par(mfrow=c(1,2))
plot(fitted(mod1),resid(mod1))
abline(h = 0, lty=2)

qqnorm(resid(mod1))
qqline(resid(mod1))
par(mfrow=c(1,2))

# Alternative LMM

  # Log-transformation of Y
ocu <- ocu %>%
  mutate(logY = log(Y))

mod2 <- lmer(logY~Origin*(sin_t + cos_t)+(1|ID), data=ocu)
par(mfrow=c(1,2))
plot(fitted(mod2),resid(mod2))
abline(h = 0, lty=2)

qqnorm(resid(mod2))
qqline(resid(mod2))
par(mfrow=c(1,2))

  # Second-order sinusoid
ocu <- ocu %>%
  mutate(sin_2t=sin(2*theta_rad),
         cos_2t=cos(2*theta_rad))

mod_2 <- lmer(Y~Origin*(sin_t+cos_t+sin_2t + cos_2t)+(1|ID), data=ocu)
par(mfrow=c(1,2))
plot(fitted(mod_2),resid(mod_2))
abline(h = 0, lty=2)
qqnorm(resid(mod_2))
qqline(resid(mod_2))

# Third-order sinusoid
ocu <- ocu %>%
  mutate(sin_3t=sin(3*theta_rad),
         cos_3t=cos(3*theta_rad))

mod_3 <- lmer(Y~Origin*(sin_t+cos_t+sin_2t + cos_2t+sin_3t + cos_3t)+(1|ID), data=ocu)
plot(fitted(mod_3),resid(mod_3))
abline(h = 0, lty=2)
qqnorm(resid(mod_3))
qqline(resid(mod_3))

par(mfrow=c(1,2))

anova(mod0,mod1,mod_2,mod_3)

# ===========================
# Modèles GAMM
# ===========================

gamm1 <- gam(
  Y~Origin+
    s(theta_rad,bs="cc",k=120) +   # k adapted so that k-index > 1
    s(ID,bs="re"),
  data=ocu,
  method="REML"
)
gam.check(gamm1)

gamm2 <- gam(
  Y ~ Origin +
    s(theta_rad, by = Origin, bs = "cc",k=120) +
    s(ID, bs = "re"),
  data = ocu,
  method = "REML"
)

anova(gamm1,gamm2,test="Chisq")
gam.check(gamm2)


gamm0 <- gam(
  Y ~ Origin +
    s(theta_rad, by = Origin, bs = "cc", k = 10) +
    s(ID, bs = "re"),
  data = ocu,
  method = "REML"
)
gam.check(gamm0)


gamm3 <- gam(
  Y ~ Origin +
    s(theta_rad, by = Origin, bs = "cc", k = 120) +
    s(ID, bs = "re"),
  data = ocu,
  method = "REML"
)
gam.check(gamm3)

gamm4 <- gam(
  Y ~ Origin +
    s(theta_rad, by = Origin, bs = "cc", k = 120) +
    s(ID, bs = "re"),
  data = ocu,
  method = "REML",
  knots = list(theta_rad = c(0, 2*pi))
)

anova(gamm1,gamm4, test = "Chisq")


# Predictions
theta_grid <- seq(0, 2*pi, length.out = 500)[-500]

newdata <- expand.grid(
  theta_rad = theta_grid,
  Origin = levels(ocu$Origin),
  ID = levels(ocu$ID)[1]   # Arbitrary valid ID
)


newdata$pred <- predict(
  model_gamm_final,
  newdata = newdata,
  exclude = "s(ID)"
)

# Chart
ggplot() +
  geom_point(
    data = ocu,
    aes(x = theta_rad, y = Y, colour = Origin),
    alpha = 0.15,
    size = 0.6
  ) +
  geom_line(
    data = newdata,
    aes(x = theta_rad, y = pred, colour = Origin),
    linewidth = 1.2
  ) +
  scale_x_continuous(
    breaks = seq(0, 2 * pi, by = pi / 2),
    labels = c("0", "π/2", "π", "3π/2", "2π")
  ) +
  labs(
    x = "Roll angle (radians)",
    y = "Ocular muscle tonus",
    colour = "Population",
    caption = "Curves estimated using a GAMM with population-specific circular splines.\nLarge sample size (n ≈ 10,000) allows detection of subtle but consistent differences."
  ) +
  theme_bw() +
  theme(
    text=element_text(size=26),
    axis.title=element_text(size=26,face="bold"),
    axis.text=element_text(size=24),
    legend.title=element_text(size=25),
    legend.text=element_text(size=24),
    plot.caption=element_text(size=22)
  )

