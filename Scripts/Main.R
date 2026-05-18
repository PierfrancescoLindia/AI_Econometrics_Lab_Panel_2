############################################################
# AI ECONOMETRICS LAB PANEL 2
# Dynamic Panel Data: Arellano-Bond Estimation
# Script completo per esercitazione
############################################################

# Clear workspace
rm(list = ls())

# Set working directory
setwd("C:/Users/lindi/Desktop/Econometria/Esercitazione 2")

# Load required packages
library(plm)
library(lmtest)

# Avoid scientific notation
options(scipen = 999)


############################################################
# PART A - PRELIMINARY INSPECTION OF THE PANEL
############################################################

# Import dataset
farms <- read.table("farms1.txt",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)

# Sort data by farm and year
farms <- farms[order(farms$FARM, farms$YEAR), ]

# First inspection
head(farms)
str(farms)
summary(farms)
names(farms)
dim(farms)

# Missing values
colSums(is.na(farms))

# Log-transformed variables
farms$lmilk <- log(farms$MILK)
farms$lcows <- log(farms$COWS)
farms$lland <- log(farms$LAND)
farms$lfeed <- log(farms$FEED)

# Check log variables
head(farms[, c("FARM", "YEAR", "MILK", "COWS", "LAND", "FEED",
               "lmilk", "lcows", "lland", "lfeed")])

summary(farms[, c("lmilk", "lcows", "lland", "lfeed")])

# Declare panel structure
pdata <- pdata.frame(farms, index = c("FARM", "YEAR"))

# Check panel structure
pdim(pdata)

# Number of farms and years
length(unique(farms$FARM))
length(unique(farms$YEAR))
range(farms$YEAR)

# Observations by farm and year
table(farms$FARM)
table(farms$YEAR)

# Descriptive statistics table
desc_stats <- data.frame(
  Variable = c("lmilk", "lcows", "lland", "lfeed"),
  
  Mean = c(
    mean(farms$lmilk, na.rm = TRUE),
    mean(farms$lcows, na.rm = TRUE),
    mean(farms$lland, na.rm = TRUE),
    mean(farms$lfeed, na.rm = TRUE)
  ),
  
  SD = c(
    sd(farms$lmilk, na.rm = TRUE),
    sd(farms$lcows, na.rm = TRUE),
    sd(farms$lland, na.rm = TRUE),
    sd(farms$lfeed, na.rm = TRUE)
  ),
  
  Min = c(
    min(farms$lmilk, na.rm = TRUE),
    min(farms$lcows, na.rm = TRUE),
    min(farms$lland, na.rm = TRUE),
    min(farms$lfeed, na.rm = TRUE)
  ),
  
  Max = c(
    max(farms$lmilk, na.rm = TRUE),
    max(farms$lcows, na.rm = TRUE),
    max(farms$lland, na.rm = TRUE),
    max(farms$lfeed, na.rm = TRUE)
  )
)

desc_stats

# Average lmilk by year
avg_lmilk_year <- aggregate(lmilk ~ YEAR,
                            data = farms,
                            FUN = mean,
                            na.rm = TRUE)

avg_lmilk_year

# Plot average lmilk over time
plot(avg_lmilk_year$YEAR, avg_lmilk_year$lmilk,
     type = "l",
     xlab = "Year",
     ylab = "Average lmilk",
     main = "Average log milk production over time")

# Plot lmilk for selected farms
selected_farms <- unique(farms$FARM)[1:6]

plot(NULL,
     xlim = range(farms$YEAR),
     ylim = range(farms$lmilk, na.rm = TRUE),
     xlab = "Year",
     ylab = "lmilk",
     main = "Log milk production over time for selected farms")

for (f in selected_farms) {
  farm_data <- farms[farms$FARM == f, ]
  lines(farm_data$YEAR, farm_data$lmilk, type = "l")
}

legend("topleft",
       legend = paste("Farm", selected_farms),
       lty = 1,
       cex = 0.8,
       bty = "n")

# Boxplot of lmilk by farm
boxplot(lmilk ~ FARM,
        data = farms,
        xlab = "Farm",
        ylab = "lmilk",
        main = "Distribution of log milk production across farms",
        col = "lightgray",
        las = 2)


############################################################
# PART B - BENCHMARK MODELS
# Pooled OLS and Fixed Effects
############################################################

# Dynamic pooled OLS
pool_mod <- plm(
  lmilk ~ lag(lmilk, 1) + lcows + lland + lfeed,
  data = pdata,
  model = "pooling"
)

summary(pool_mod)

# Dynamic fixed effects
fe_mod <- plm(
  lmilk ~ lag(lmilk, 1) + lcows + lland + lfeed,
  data = pdata,
  model = "within",
  effect = "individual"
)

summary(fe_mod)

# Robust standard errors clustered by farm
pool_robust <- coeftest(pool_mod,
                        vcov = vcovHC(pool_mod,
                                      type = "HC1",
                                      cluster = "group"))

fe_robust <- coeftest(fe_mod,
                      vcov = vcovHC(fe_mod,
                                    type = "HC1",
                                    cluster = "group"))

pool_robust
fe_robust

# Benchmark comparison
benchmark_comparison <- data.frame(
  Variable = c("lagged_lmilk", "lcows", "lland", "lfeed"),
  
  Pooled_OLS = c(
    coef(pool_mod)["lag(lmilk, 1)"],
    coef(pool_mod)["lcows"],
    coef(pool_mod)["lland"],
    coef(pool_mod)["lfeed"]
  ),
  
  Fixed_Effects = c(
    coef(fe_mod)["lag(lmilk, 1)"],
    coef(fe_mod)["lcows"],
    coef(fe_mod)["lland"],
    coef(fe_mod)["lfeed"]
  )
)

benchmark_comparison

nobs(pool_mod)
nobs(fe_mod)


############################################################
# PART C - ARELLANO-BOND DIFFERENCE GMM
# One-step and two-step
############################################################

# Arellano-Bond one-step
ab_mod_1 <- pgmm(
  lmilk ~ lag(lmilk, 1) + lcows + lland + lfeed |
    lag(lmilk, 2:5),
  data = pdata,
  effect = "individual",
  model = "onestep",
  transformation = "d"
)

summary(ab_mod_1, robust = TRUE)

# Arellano-Bond two-step
ab_mod_2 <- pgmm(
  lmilk ~ lag(lmilk, 1) + lcows + lland + lfeed |
    lag(lmilk, 2:5),
  data = pdata,
  effect = "individual",
  model = "twosteps",
  transformation = "d"
)

summary(ab_mod_2, robust = TRUE)

# Comparison one-step vs two-step
ab_comparison <- data.frame(
  Variable = c("lagged_lmilk", "lcows", "lland", "lfeed"),
  
  AB_One_Step = c(
    coef(ab_mod_1)["lag(lmilk, 1)"],
    coef(ab_mod_1)["lcows"],
    coef(ab_mod_1)["lland"],
    coef(ab_mod_1)["lfeed"]
  ),
  
  AB_Two_Step = c(
    coef(ab_mod_2)["lag(lmilk, 1)"],
    coef(ab_mod_2)["lcows"],
    coef(ab_mod_2)["lland"],
    coef(ab_mod_2)["lfeed"]
  )
)

ab_comparison

# Serial correlation tests
mtest(ab_mod_1, order = 1)
mtest(ab_mod_1, order = 2)

mtest(ab_mod_2, order = 1)
mtest(ab_mod_2, order = 2)

# Sargan tests
sargan(ab_mod_1)
sargan(ab_mod_2)


############################################################
# ALTERNATIVE INSTRUMENT STRATEGY
# Reduced instrument set: lag(lmilk, 2:3)
############################################################

ab_mod_small <- pgmm(
  lmilk ~ lag(lmilk, 1) + lcows + lland + lfeed |
    lag(lmilk, 2:3),
  data = pdata,
  effect = "individual",
  model = "twosteps",
  transformation = "d"
)

summary(ab_mod_small, robust = TRUE)

# Comparison baseline vs reduced instruments
instrument_comparison <- data.frame(
  Variable = c("lagged_lmilk", "lcows", "lland", "lfeed"),
  
  AB_lag_2_5 = c(
    coef(ab_mod_2)["lag(lmilk, 1)"],
    coef(ab_mod_2)["lcows"],
    coef(ab_mod_2)["lland"],
    coef(ab_mod_2)["lfeed"]
  ),
  
  AB_lag_2_3 = c(
    coef(ab_mod_small)["lag(lmilk, 1)"],
    coef(ab_mod_small)["lcows"],
    coef(ab_mod_small)["lland"],
    coef(ab_mod_small)["lfeed"]
  )
)

instrument_comparison

# Diagnostics for reduced instrument model
mtest(ab_mod_small, order = 1)
mtest(ab_mod_small, order = 2)
sargan(ab_mod_small)

# Diagnostic comparison
diagnostic_comparison <- data.frame(
  Model = c("AB one-step lag 2:5",
            "AB two-step lag 2:5",
            "AB two-step lag 2:3"),
  
  AR1_pvalue = c(
    mtest(ab_mod_1, order = 1)$p.value,
    mtest(ab_mod_2, order = 1)$p.value,
    mtest(ab_mod_small, order = 1)$p.value
  ),
  
  AR2_pvalue = c(
    mtest(ab_mod_1, order = 2)$p.value,
    mtest(ab_mod_2, order = 2)$p.value,
    mtest(ab_mod_small, order = 2)$p.value
  ),
  
  Sargan_pvalue = c(
    sargan(ab_mod_1)$p.value,
    sargan(ab_mod_2)$p.value,
    sargan(ab_mod_small)$p.value
  )
)

diagnostic_comparison


############################################################
# ALTERNATIVE CLASSIFICATION OF REGRESSORS
# FEED treated as predetermined
############################################################

ab_feed_pred <- pgmm(
  lmilk ~ lag(lmilk, 1) + lcows + lland + lfeed |
    lag(lmilk, 2:3) + lag(lfeed, 2:3),
  data = pdata,
  effect = "individual",
  model = "twosteps",
  transformation = "d"
)

summary(ab_feed_pred, robust = TRUE)

# Comparison with reduced baseline
feed_comparison <- data.frame(
  Variable = c("lagged_lmilk", "lcows", "lland", "lfeed"),
  
  Baseline_lag_2_3 = c(
    coef(ab_mod_small)["lag(lmilk, 1)"],
    coef(ab_mod_small)["lcows"],
    coef(ab_mod_small)["lland"],
    coef(ab_mod_small)["lfeed"]
  ),
  
  Feed_predetermined = c(
    coef(ab_feed_pred)["lag(lmilk, 1)"],
    coef(ab_feed_pred)["lcows"],
    coef(ab_feed_pred)["lland"],
    coef(ab_feed_pred)["lfeed"]
  )
)

feed_comparison

# Diagnostics for feed predetermined model
mtest(ab_feed_pred, order = 1)
mtest(ab_feed_pred, order = 2)
sargan(ab_feed_pred)

# Diagnostic comparison including feed predetermined model
feed_diagnostics <- data.frame(
  Model = c("Baseline lag 2:3",
            "Feed predetermined"),
  
  AR1_pvalue = c(
    mtest(ab_mod_small, order = 1)$p.value,
    mtest(ab_feed_pred, order = 1)$p.value
  ),
  
  AR2_pvalue = c(
    mtest(ab_mod_small, order = 2)$p.value,
    mtest(ab_feed_pred, order = 2)$p.value
  ),
  
  Sargan_pvalue = c(
    sargan(ab_mod_small)$p.value,
    sargan(ab_feed_pred)$p.value
  )
)

feed_diagnostics


############################################################
# FINAL OUTPUTS TO KEEP FOR THE REPORT
############################################################

# Dataset and descriptive statistics
pdim(pdata)
desc_stats
avg_lmilk_year

# Benchmark models
summary(pool_mod)
summary(fe_mod)
pool_robust
fe_robust
benchmark_comparison
nobs(pool_mod)
nobs(fe_mod)

# Arellano-Bond models
summary(ab_mod_1, robust = TRUE)
summary(ab_mod_2, robust = TRUE)
ab_comparison

# Diagnostics
mtest(ab_mod_2, order = 1)
mtest(ab_mod_2, order = 2)
sargan(ab_mod_2)

# Reduced instrument strategy
summary(ab_mod_small, robust = TRUE)
instrument_comparison
diagnostic_comparison

# Feed predetermined
summary(ab_feed_pred, robust = TRUE)
feed_comparison
feed_diagnostics
