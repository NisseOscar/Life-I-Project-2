
########## Task 1 Code ###############
# Run this for different parameter choices to get the results in the report

n.sim <- 1000
max_cens <- 10
set.seed(123)

## all individuals start in state 0

## transition times Tij := i -> j
T02 <- rweibull(n = n.sim, shape = 2, scale = 1)
T01 <- rweibull(n = n.sim, shape = 2, scale = 0.5)
T12 <- rweibull(n = n.sim, shape = 2, scale = 0.5)
Tcens <- runif(n = n.sim, min = 0, max = max_cens)


## -------------------------------
## take care of transition times
## -------------------------------
## first jump
T1 <- pmin(T01, T02)

## first jump potentially censored
T1c <- pmin(T1, Tcens)

## second (potential) jump
T2 <- T01 + T12

## second jump potentially censored
T2c <- pmin(T2, Tcens)

## -------------------------------
## format data
## -------------------------------
## structure of dN = matrix with
## jumps
## dN <- c(0 -> 1, 0 -> 2, 1 -> 2)
dN <- c()

## dY = matrix with changes in
## the at risk process
dY <- c()

## Tj = jump times
Tj <- c()
for (i in (1:n.sim)) {
    Tj <- rbind(Tj, T1c[i])
    if (T1c[i] == T1[i]) {
        if (T1[i] == T01[i]) {
            ## 0 -> 1
            dN <- rbind(dN, c(1, 0, 0))
            dY <- rbind(dY, c(-1, 1, 0))

            ## check if second jump
            Tj <- rbind(Tj, T2c[i])
            if (T2c[i] == T2[i]) {
                ## 1 -> 2
                dN <- rbind(dN, c(0, 0, 1))
                dY <- rbind(dY, c(0, -1, 1))
            } else {
                dN <- rbind(dN, c(0, 0, 0))
                dY <- rbind(dY, c(0, -1, 0))
            }
        } else {
            ## 0 -> 2
            dN <- rbind(dN, c(0, 1, 0))
            dY <- rbind(dY, c(-1, 0, 1))
        }
    } else {
        dN <- rbind(dN, c(0, 0, 0))
        dY <- rbind(dY, c(-1, 0, 0))
    }
}

## sort the unsorted jump times,
## and sort the other matrices
## accordingly
sorted.t <- sort(Tj, index.return = TRUE)
idx.sorted.t <- sorted.t$ix

dN <- dN[idx.sorted.t, ]
dY <- dY[idx.sorted.t, ]
Tj <- Tj[idx.sorted.t]


## add initial states, all individuals
## start in state 0
dN <- rbind(c(0, 0, 0), dN)
Y <- rbind(c(n.sim, 0, 0), dY)
Y <- apply(Y, FUN = cumsum, MARGIN = 2)
Tj <- c(0, Tj)

## initiate transition matrices and matrices
## used for nelson-aalen increments
I <- diag(1, nrow = 3, ncol = 3)
Z <- matrix(0, nrow = 3, ncol = 3)

## array with sequence of transition probabilities
## rows correspond to state 0, 1, 2 from top down
Pt <- array(data = 0, dim = c(3, 3, length(Tj)))
Pt[, , 1] <- I

## function calculating the aalen-johansen estimator
aalen.johansen <- function(dN, Y, Tj) {
    for (i in (2:length(Tj))) {
        dA <- Z
        dA[1, ] <- c(-(dN[i, 1] + dN[i, 2]) / Y[i - 1, 1], dN[i, 1] / Y[i - 1, 1], dN[i, 2] / Y[i - 1, 1])
        dA[2, ] <- c(0, -dN[i, 3] / Y[i - 1, 2], dN[i, 3] / Y[i - 1, 2])
        dA[3, ] <- c(0, 0, 0)
        dA[is.na(dA)] <- 0
        Pt[, , i] <- Pt[, , i - 1] %*% (I + dA)
    }

    return(Pt)
}

## calculate the aalen-johanesen estimator
Pt.aalen.johansen <- aalen.johansen(dN, Y, Tj)

# take out estimated transition probabilities
T11 <- Pt.aalen.johansen[1, 1, ]
T12 <- Pt.aalen.johansen[1, 2, ]
T13 <- Pt.aalen.johansen[1, 3, ]
T22 <- Pt.aalen.johansen[2, 2, ]
T23 <- Pt.aalen.johansen[2, 3, ]

# Plot all in a joint plot
# png(paste("plots/n", n.sim, ".png", sep = ""), width = 400, height = 600)
png(paste("plots/cens", max_cens, ".png", sep = ""), width = 600, height = 500)
plot(Tj, T11, type = "s", ylim = c(0, 1), xlim = c(0, 2.2), ylab = "P_ij(0, t)", xlab = "t")
lines(Tj, T12, type = "s", col = "red")
lines(Tj, T13, type = "s", col = "blue")
lines(Tj, T22, type = "s", col = "green")
lines(Tj, T23, type = "s", col = "purple")
legend("topright", legend = c("P_00", "P_01", "P_02", "P_11", "P_12"), col = c("black", "red", "blue", "green", "purple"), lty = 1)
dev.off()


############## Task 2 Code ################
library("gnm")
library("forecast")
library("StMoMo")
library("demography")
library(ggplot2)
library(dplyr)
library(reshape2)

setwd(".")

MortalityDataSWE <- read.demogdata("SWE_Mx_1x1.txt", "SWE_Exposures_1x1.txt", type = "mortality", label = "Sweden")
StMoMoSWE <- StMoMoData(MortalityDataSWE, series = "female")

# check which model you have defined
mortalityFunc$textFormula

##########  Task 1 ############

fitAges <- 0:90 # ages to be used for model estimation
fitYrs <- 1900:2010 # calendar years to be used for model estimation
totYrs <- 1900:2020 # check upper year limit from data

# nice to have for model evaluation: observed death rates
obsDeathRateSWE <- StMoMoSWE$Dxt[as.character(fitAges), as.character(totYrs)] / StMoMoSWE$Ext[as.character(fitAges), as.character(totYrs)]

# Plot Lexis diagram of survival rate
death_data <- data.frame(obsDeathRateSWE)
death_data$Age <- fitAges
death_data <- melt(death_data, id.vars = "Age", variable.name = "Year", value.name = "DeathRate")
death_data$Year <- as.numeric(gsub("X", "", death_data$Year))
death_data %>% tail(30)

# Plot Lexis diagram
p <- ggplot(death_data, aes(x = Year, y = Age, fill = DeathRate)) +
    geom_tile() + # Fill grid cells with death rates
    scale_fill_gradient(low = "lightblue", high = "blue", name = "Death Rate", trans = "log", breaks = c(0.001, 0.01, 0.1, 1, 10)) +
    geom_abline(intercept = seq(-100, 100, by = 10), slope = 1, color = "grey80", linetype = "dotted") + # Diagonal lines
    labs(
        x = "Calendar Year",
        y = "Age",
        fill = "Death Rate"
    ) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10))
ggsave("./plots/task2/Lexis_diagram.png", plot = p, width = 32, height = 16, units = "cm", dpi = 300)
p


############ Modelling ##############################################

##############################################################################
# Mortality Forecasting
##############################################################################

# Import Swedish mortality data
MortalityDataSWE <- read.demogdata("SWE_Mx_1x1.txt", "SWE_Exposures_1x1.txt", type = "mortality", label = "Sweden")
StMoMoSWE <- StMoMoData(MortalityDataSWE, series = "female")

mortalityFunc <- lc(link = "log") # corresponds to the Poisson log-bilinear model

# check which model you have defined
mortalityFunc$textFormula


fitAges <- 20:90 # ages to be used for model estimation
fitYrs <- 1990:2010 # calendar years to be used for model estimation
totYrs <- 1990:2020 # check upper year limit from data

png(paste("plots/Observed_death_rate_SWE.png", sep = ""), width = 800, height = 600)
# nice to have for model evaluation: observed death rates
obsDeathRateSWE <- StMoMoSWE$Dxt[as.character(fitAges), as.character(totYrs)] / StMoMoSWE$Ext[as.character(fitAges), as.character(totYrs)]
plot(fitAges, obsDeathRateSWE[, 1], type = "l", ylab = "Mortality Rate", xlab = "Age", main = "Observed Mortality Rate in Sweden", col = "blue")
# Add the second column to the same plot
lines(fitAges, obsDeathRateSWE[, 11], col = "red")
lines(fitAges, obsDeathRateSWE[, 21], col = "green")
# Add a legend to differentiate the lines
legend("topleft", legend = c("1990", "2000", "2010"), col = c("blue", "red", "green"), lty = 1)
dev.off()

# estimate the Poisson log-bilinear model
LCfitSWE <- fit(mortalityFunc, data = StMoMoSWE, ages.fit = fitAges, years.fit = fitYrs)

# have a look at the results
png(paste("plots/Model_fit_SWE.png", sep = ""), width = 800, height = 600)
plot(LCfitSWE)
dev.off()

# Make forecasts
LCsim <- simulate(LCfitSWE, nsim = 5000, h = 10)
LCfor <- forecast(LCfitSWE, h = 10)

mxt <- LCfitSWE$Dxt / LCfitSWE$Ext
mxtHat <- fitted(LCfitSWE, type = "rates")
mxtCentral <- LCfor$rates
mxtPred2.5 <- apply(LCsim$rates, c(1, 2), quantile, probs = 0.05)
mxtPred97.5 <- apply(LCsim$rates, c(1, 2), quantile, probs = 0.95)

x <- c("20", "40", "60", "75", "90")
x1 <- c(1, 21, 41, 56, 71)

png(paste("plots/Mortality_rate_SWE.png", sep = ""), width = 600, height = 500)
# Plot the model and forecast
matplot(LCfitSWE$years, t(mxt[x, ]), xlim = range(LCfitSWE$years, LCfor$years), type = "p", xlab = "Year", ylab = "Mortality rate", main = "Mortality rate in Sweden", pch = 20, col = "black", log = "y")
matlines(LCfitSWE$years, t(mxtHat[x, ]), lty = 1, col = "red")
matlines(LCfor$years, t(mxtCentral[x, ]), lty = 4, col = "red")
matlines(LCsim$years, t(mxtPred2.5[x, ]), lty = 3, col = "blue")
matlines(LCsim$years, t(mxtPred97.5[x, ]), lty = 3, col = "blue")

# Add observed values to the plot
for (i in 1:length(x1)) {
    points(LCsim$years, obsDeathRateSWE[x1[i], 22:31], pch = 20, col = "black")
}
text(1990, mxtHat[x, "1995"], labels = c("x=20", "x=40", "x=60", "x=75", "x=90"))
dev.off()

# Plot residuals for the model
png(paste("plots/Residuals_SWE.png", sep = ""), width = 800, height = 600)
res_SWE <- residuals(LCfitSWE)
plot(res_SWE, type = "scatter", reslim = c(-7, 7))
dev.off()


##############################################################################
# Multi-population
##############################################################################

# Import Danish mortality data
MortalityDataDNK <- read.demogdata("DNK_Mx_1x1.txt", "DNK_Exposures_1x1.txt", type = "mortality", label = "Denmark")
StMoMoDNK <- StMoMoData(MortalityDataDNK, series = "male")

# Estimate the Poisson log-bilinear model for Danish male population
LCfitDNK <- fit(mortalityFunc, data = StMoMoDNK, ages.fit = fitAges, years.fit = fitYrs)

# Observed mortality rate
obsDeathRateDNK <- StMoMoDNK$Dxt[as.character(fitAges), as.character(totYrs)] / StMoMoDNK$Ext[as.character(fitAges), as.character(totYrs)]

# Make forecast
LCsimDNK <- simulate(LCfitDNK, nsim = 5000, h = 10)
LCforDNK <- forecast(LCfitDNK, h = 10)

mxtDNK <- LCfitDNK$Dxt / LCfitDNK$Ext
mxtHatDNK <- fitted(LCfitDNK, type = "rates")
mxtCentralDNK <- LCforDNK$rates
mxtPred2.5DNK <- apply(LCsimDNK$rates, c(1, 2), quantile, probs = 0.05)
mxtPred97.5DNK <- apply(LCsimDNK$rates, c(1, 2), quantile, probs = 0.95)

x <- c("20", "40", "60", "75", "90")
x1 <- c(1, 21, 41, 56, 71)

png(paste("plots/Mortality_rate_DK.png", sep = ""), width = 800, height = 600)
# Plot the model and forecast
matplot(LCfitDNK$years, t(mxtDNK[x, ]), xlim = range(LCfitDNK$years, LCforDNK$years), type = "p", xlab = "Year", ylab = "Mortality rate", main = "Mortality rate in Denmark", pch = 20, col = "black", log = "y")
matlines(LCfitDNK$years, t(mxtHatDNK[x, ]), lty = 1, col = "red")
matlines(LCforDNK$years, t(mxtCentralDNK[x, ]), lty = 4, col = "red")
matlines(LCsimDNK$years, t(mxtPred2.5DNK[x, ]), lty = 3, col = "blue")
matlines(LCsimDNK$years, t(mxtPred97.5DNK[x, ]), lty = 3, col = "blue")

# Add observed values
for (i in 1:length(x1)) {
    points(LCsimDNK$years, obsDeathRateDNK[x1[i], 22:31], pch = 20, col = "black")
}
text(1990, mxtHatDNK[x, "1995"], labels = c("x=20", "x=40", "x=60", "x=75", "x=90"))
dev.off()

# Plot residuals
png(paste("plots/Residuals_DK.png", sep = ""), width = 800, height = 600)
res_DNK <- residuals(LCfitDNK)
plot(res_DNK, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

######
# Joint population
######

# Import Norwegian mortality data
MortalityDataNOR <- read.demogdata("NOR_Mx_1x1.txt", "NOR_Exposures_1x1.txt", type = "mortality", label = "Norway")
StMoMoNOR <- StMoMoData(MortalityDataNOR, series = "male")

# estimate the Poisson log-bilinear model for Norwegian male population
LCfitNOR <- fit(mortalityFunc, data = StMoMoNOR, ages.fit = fitAges, years.fit = fitYrs)

# Merging data into a single population
StMoMoTotal <- StMoMoNOR
StMoMoTotal$Ext <- StMoMoTotal$Ext + StMoMoDNK$Ext[, as.character(StMoMoTotal$years)]
StMoMoTotal$Dxt <- StMoMoTotal$Dxt + StMoMoDNK$Dxt[, as.character(StMoMoTotal$years)]

# estimate the Poisson log-bilinear model for joint Danish and Norwegian male population
LCfitTotal <- fit(mortalityFunc, data = StMoMoTotal, ages.fit = fitAges, years.fit = fitYrs)

# Observed mortality rates
obsDeathRateTotal <- StMoMoTotal$Dxt[as.character(fitAges), as.character(totYrs)] / StMoMoTotal$Ext[as.character(fitAges), as.character(totYrs)]

# Make forecast
LCsimTotal <- simulate(LCfitTotal, nsim = 5000, h = 10)
LCforTotal <- forecast(LCfitTotal, h = 10)

mxtTotal <- LCfitTotal$Dxt / LCfitTotal$Ext
mxtHatTotal <- fitted(LCfitTotal, type = "rates")
mxtCentralTotal <- LCforTotal$rates
mxtPred2.5Total <- apply(LCsimTotal$rates, c(1, 2), quantile, probs = 0.05)
mxtPred97.5Total <- apply(LCsimTotal$rates, c(1, 2), quantile, probs = 0.95)

x <- c("20", "40", "60", "75", "90")
x1 <- c(1, 21, 41, 56, 71)

png(paste("plots/Mortality_rate_NOK_DK.png", sep = ""), width = 800, height = 600)
# Plot model and forecast
matplot(LCfitTotal$years, t(mxtTotal[x, ]), xlim = range(LCfitTotal$years, LCforTotal$years), type = "p", xlab = "Year", ylab = "Mortality rate", main = "Mortality rate in Norway and Denmark", pch = 20, col = "black", log = "y")
matlines(LCfitTotal$years, t(mxtHatTotal[x, ]), lty = 1, col = "red")
matlines(LCforTotal$years, t(mxtCentralTotal[x, ]), lty = 4, col = "red")
matlines(LCsimTotal$years, t(mxtPred2.5Total[x, ]), lty = 3, col = "blue")
matlines(LCsimTotal$years, t(mxtPred97.5Total[x, ]), lty = 3, col = "blue")

# Add observed values
for (i in 1:length(x1)) {
    points(LCsimTotal$years, obsDeathRateTotal[x1[i], 22:31], pch = 20, col = "black")
}
text(1990, mxtHatTotal[x, "1995"], labels = c("x=20", "x=40", "x=60", "x=75", "x=90"))
dev.off()

# Plot residuals
res_Tot <- residuals(LCfitTotal)
png(paste("plots/Residuals_NOK_DK.png", sep = ""), width = 800, height = 600)
plot(res_Tot, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()





# Calculate ratio for every age
ages <- as.character(fitAges)
ratiosDK <- sapply(ages, function(age) {
    sum(StMoMoDNK$Dxt[age, as.character(totYrs)]) /
        sum(StMoMoDNK$Ext[age, as.character(totYrs)] * obsDeathRateTotal[age, as.character(totYrs)])
})

# Plot ratios
png(paste("plots/Ratios_DK.png", sep = ""), width = 600, height = 500)
plot(ages, ratiosDK, xlab = "Age", ylab = expression(hat(pi)(x)), main = "Ratios for Danish population", pch = 20)
dev.off()

# Calculate alternative mortality rate for the Danish male population
mu_hat_DK <- matrix(0, nrow = nrow(mxtHatTotal), ncol = ncol(mxtHatTotal))
rownames(mu_hat_DK) <- rownames(mxtHatTotal)

mu_hat_pred_DK <- matrix(0, nrow = nrow(mxtCentralTotal), ncol = ncol(mxtCentralTotal))
rownames(mu_hat_pred_DK) <- rownames(mxtHatTotal)
for (age in 1:71) {
    mu_hat_DK[age, ] <- ratiosDK[age] * mxtHatTotal[age, ]
    mu_hat_pred_DK[age, ] <- ratiosDK[age] * mxtCentralTotal[age, ]
}

png(paste("plots/Comparison_alternative_mortality_rate.png", sep = ""), width = 600, height = 500)
# Plot the model and forecast
matplot(LCfitDNK$years, t(mxtDNK[x, ]), xlim = range(LCfitDNK$years, LCforDNK$years), type = "p", xlab = "Year", ylab = "Mortality rate", main = "Comparison between alternative mortality rate in Denmark", pch = 20, col = "black", log = "y")
matlines(LCfitDNK$years, t(mxtHatDNK[x, ]), lty = 1, col = "red")
matlines(LCforDNK$years, t(mxtCentralDNK[x, ]), lty = 4, col = "red")
matlines(LCfitDNK$years, t(mu_hat_DK[x, ]), lty = 1, col = "green")
matlines(LCforDNK$years, t(mu_hat_pred_DK[x, ]), lty = 4, col = "green")
matlines(LCsimDNK$years, t(mxtPred2.5DNK[x, ]), lty = 3, col = "blue")
matlines(LCsimDNK$years, t(mxtPred97.5DNK[x, ]), lty = 3, col = "blue")

# Add observed values
for (i in 1:length(x1)) {
    points(LCsimDNK$years, obsDeathRateDNK[x1[i], 22:31], pch = 20, col = "black")
}
text(1990, mxtHatDNK[x, "1995"], labels = c("x=20", "x=40", "x=60", "x=75", "x=90"))
dev.off()