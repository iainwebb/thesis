# Setup ----

# replace back in
# ACURE_P3_AOD_Total_jul_level2 for A_O_D
# ERF_PPE_global_jul for ERF
options(scipen=999)

library(readr)
# library(readr, lib="C:/Users/smp22ijw/Desktop/Library/") # if Lenovo
library(DiceKriging)
# library(DiceKriging, lib="C:/Users/smp22ijw/Desktop/Library/") # if Lenovo
# https://colordesigner.io/gradient-generator

# Importing input settings ----
inputs_x_norm_T_matrix <- read.csv("data/ACURE-UKESM1-PPE-par-norm-Design.csv")
p <- ncol(inputs_x_norm_T_matrix)

# Importing observable variable ----

## AOD_Total
library(ncdf4)
# library(ncdf4, lib="C:/Users/smp22ijw/Desktop/Library/") # if Lenovo
ACURE_P3_AOD_Total_jul_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc")
)
ACURE_P3_AOD_Total_jul_array <- ncvar_get(
  ACURE_P3_AOD_Total_jul_nc, names(ACURE_P3_AOD_Total_jul_nc$var)
)
response_A_O_D_array <- ACURE_P3_AOD_Total_jul_array[,,2,]
rm("ACURE_P3_AOD_Total_jul_nc", "ACURE_P3_AOD_Total_jul_array")

# Importing unobservable variable ----

response_ERF_dataframe <- read.table(
  "X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/global_mean_forcing/ERF_PPE_global_jul.dat", 
  col.names = "ERF"
)
# response_ERF_dataframe <- read.table(
#   "data/ERF_PPE_global_jul.dat", 
#   col.names = "ERF"
# ) # if Lenovo

# Importing map-required things ----
ACURE_P3_AOD_Total_jul_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc")
)
# ACURE_P3_AOD_Total_jul_nc <- nc_open(paste("data/ACURE_P3_AOD_Total_jul.nc")) # if Lenovo
longitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"longitude") - 180
latitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"latitude")
rm("ACURE_P3_AOD_Total_jul_nc")

# lm is ~. - plot a GAM smoother for each of the 37 for ERF, then for some rogue AOD, and compare ----
# taken from https://uq4covid.github.io/vignettes/SensitivityAnalysisWithGAM
library(tidyverse)

# ERF ----
## Inputs
X <- inputs_x_norm_T_matrix
parnames = colnames(X)
# dummy input
set.seed(123)
X$dummy <- runif(221)

## Outputs
dat <- response_ERF_dataframe

# X$dat <- dat

X <- cbind(X, dat)

## Plots
X_ordered <- X %>% 
  as_tibble %>% 
  # mutate(y=ERF) %>% 
  gather('parameter', 'value', -ERF)
X_ordered <- transform(X_ordered,
                       parameter=factor(parameter, levels = c(parnames, "dummy")))
X_ordered %>% 
  ggplot(aes(x=value, y=ERF)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~parameter) +
  labs(y='output', x='input') +
  geom_smooth(method = mgcv::gam,
              formula = y ~ s(x, bs = "tp"),
              fill = "red",
              method.args = list(method="GCV.Cp")) +
  ylim(-4,6)

mainEffects <- rep(0, 38)
for(i in 1:38){
  gam1 <- mgcv::gam(X[, 39] ~ s(X[, i]))
  mainEffects[i] <- var(gam1$fitted) / var(X[, 39])
}
round(mainEffects * 100)
# Increase margin size
par(mar=c(1,3,1,0.5))
# par(mar=c(6.5,3,1,0.5))
barplot(mainEffects, 
        # names.arg = "",
        names.arg = c(parnames, "dummy"),
        las=2, cex.names=0.75,
        ylim = c(0,0.7)
        # , ylab = "Main effects for Total ERF (July)"
)
range(mainEffects)
mainEffects[38]

cumsum(sort(round(mainEffects * 100), decreasing=T))
sort(round(mainEffects * 100), decreasing=T)[1:9]
input_names <- colnames(inputs_x_norm_T_matrix)
input_names[sort(round(mainEffects * 100), decreasing=T, index.return=T)$ix[1:9]]
ERF_SAed_indices <- sort(round(mainEffects * 100), decreasing=T, index.return=T)$ix[1:9]
rm(list=setdiff(ls(), "ERF_SAed_indices"))
input_names[ERF_SAed_indices]

# AOD ----
## Inputs
X <- inputs_x_norm_T_matrix
parnames = colnames(X)
# dummy input
set.seed(123)
X$dummy <- runif(221)

## Outputs
# long <- 44; lat <- 89 # 0.49 R2adj, 0.92 AM
# long <- 44; lat <- 90 # 0.25 R2adj, 0.06 AM
long <- 43; lat <- 90 # 0.77 R2adj, 0.37 AM
R2_adj_A_O_D_matrix[long,lat]
sum(!is.na(AM_A_O_D_versus_E_R_F_only_gp_matrix))
192*144 - sum(!is.na(AM_A_O_D_versus_E_R_F_matrix))
AM_A_O_D_versus_E_R_F_matrix[long,lat]

longitude[96+long]
latitude[lat]
dat <- response_A_O_D_array[long, lat,]

# X$dat <- dat

X <- cbind(X, "AOD_July_43_90" = dat)

## Plots
X_ordered <- X %>% 
  as_tibble %>%
  gather('parameter', 'value', -AOD_July_43_90)
X_ordered <- transform(X_ordered,
                       parameter=factor(parameter, levels = c(parnames, "dummy")))

X_ordered %>% 
  ggplot(aes(x=value, y=AOD_July_43_90)) + 
  geom_point(alpha = 0.5) + 
  facet_wrap(~parameter) +
  labs(y='output', x='input') +
  geom_smooth(method = mgcv::gam,
              formula = y ~ s(x, bs = "tp"),
              fill = "red",
              method.args = list(method="GCV.Cp")) +
  ylim(-4,6)

  
range(response_A_O_D_array[long, lat,])

mainEffects <- rep(0, 38)
for(i in 1:38){
  gam1 <- mgcv::gam(X[, 39] ~ s(X[, i]))
  mainEffects[i] <- var(gam1$fitted) / var(X[, 39])
}
round(mainEffects * 100)
# Increase margin size
par(mar=c(1,3,1,0.5))
# par(mar=c(6.5,3,1,0.5))
barplot(mainEffects, 
        # names.arg = "",
        names.arg = c(parnames, "dummy"),
        las=2, cex.names=0.75,
        ylim = c(0,0.7)
        # , ylab = "Main effects for Total ERF (July)"
)
range(mainEffects)
mainEffects[38]
parnames
### import
R2_adj_A_O_D_matrix <-
  matrix(as.numeric(readLines("objects/R2_adj_A_O_D_matrix.txt")),
         nrow = length(longitude))
R2_adj_A_O_D_matrix[43, 90]

R2_adj_E_R_F_matrix <-
  matrix(as.numeric(readLines("objects/R2_adj_E_R_F_matrix.txt")),
         nrow = length(longitude))

AM_A_O_D_versus_E_R_F_matrix <-
  matrix(as.numeric(readLines(paste0("objects/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N, ".txt"))),
         nrow = length(longitude))

## Selecting the two gridboxes

N <- 100000

### Map
AM_A_O_D_versus_E_R_F_only_gp_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_A_O_D_versus_E_R_F_only_gp_matrix[long, lat] <- 
      if (R2_adj_A_O_D_matrix[long, lat] < 0.6) AM_A_O_D_versus_E_R_F_matrix[long, lat]
    else NA
  }
  cat(long, "/192, ", sep = "")
}

library(fields)
library(maps)
image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_only_gp_matrix[97:192,],
                 AM_A_O_D_versus_E_R_F_only_gp_matrix[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("red4","red3","red1","indianred1", "mistyrose",
                   "honeydew1", "darkseagreen1", "palegreen2", "seagreen3", 
                   "springgreen4"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

### Finding indices of the two

threshold <- 0.6
double_gp_indices <- which(R2_adj_A_O_D_matrix < threshold & R2_adj_E_R_F_matrix < threshold, arr.ind = TRUE)
double_gp_indices_dataframe <- as.data.frame(double_gp_indices)
double_gp_indices_sorted_dataframe <- double_gp_indices_dataframe[order(double_gp_indices_dataframe[,1]),]
rownames(double_gp_indices_sorted_dataframe) <- 1:nrow(double_gp_indices_dataframe)
double_gp_indices_sorted_dataframe

range(double_gp_indices_sorted_dataframe[,2])
hist(double_gp_indices_sorted_dataframe[,1])
double_gp_indices_sorted_dataframe[33:47,]

### Map to confirm

AM_A_O_D_versus_E_R_F_two_only_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_A_O_D_versus_E_R_F_two_only_matrix[long, lat] <- 
      if (long == 43 & (lat == 90 | lat == 90)) AM_A_O_D_versus_E_R_F_only_gp_matrix[long, lat]
    else NA
  }
  cat(long, "/192, ", sep = "")
}
table(AM_A_O_D_versus_E_R_F_two_only_matrix)

image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_two_only_matrix[97:192,],
                 AM_A_O_D_versus_E_R_F_two_only_matrix[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("red4","red3","red1","indianred1", "mistyrose",
                   "honeydew1", "darkseagreen1", "palegreen2", "seagreen3", 
                   "springgreen4"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")
