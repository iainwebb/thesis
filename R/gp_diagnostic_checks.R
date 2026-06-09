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

# LOO CV for ERF ----

## km() ----

LOOCVPred <- data.frame()
mean <- lower95 <- upper95 <- rep(NA, 221)
LOOCVPred <- data.frame(cbind(mean, lower95, upper95))
for (obs in 1:221) {
  gp <- km(~., 
           design = inputs_x_norm_T_matrix[-obs,], 
           response = response_ERF_dataframe[-obs,], 
           covtype="gauss", optim.method="BFGS", control=list(maxit=500)
          )
  pred <- predict(object=gp,
                       newdata=data.frame(inputs_x_norm_T_matrix[obs,]),
                       type="UK",checkNames=FALSE,light.return=TRUE
                       )
  LOOCVPred[obs, 1] <- pred$mean
  LOOCVPred[obs, 2] <- pred$lower95
  LOOCVPred[obs, 3] <- pred$upper95
}

df <- data.frame(
  cbind("E_R_F" = response_E_R_F_dataframe, inputs_x_norm_T_matrix)
)
lin_mod <- lm(E_R_F ~., df)
summary(lin_mod)$adj.r.squared

LOOCVPred <- cbind(LOOCVPred, 
                   response_ERF_dataframe, 
                   "lmfittedvalues" = lin_mod$fitted.values
)
minXY <- min(LOOCVPred$lower95, LOOCVPred$ERF, LOOCVPred$lmfittedvalues)
maxXY <- max(LOOCVPred$upper95, LOOCVPred$ERF, LOOCVPred$lmfittedvalues)

par(mar = c(5.1, 4.1, 4.1, 2.1))

em_true_cols <- rep("red", 221)
em_true_cols[LOOCVPred$ERF <= LOOCVPred$upper95 & 
               LOOCVPred$ERF >= LOOCVPred$lower95] <- "black"

plot(mean ~ ERF, LOOCVPred,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = em_true_cols
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:Nval){
  lines(c(LOOCVPred$ERF[obs],LOOCVPred$ERF[obs]),
        c(LOOCVPred$lower95[obs],LOOCVPred$upper95[obs]),lwd=1.2,
        col = em_true_cols[obs]
  )
}
round((221 - sum(LOOCVPred$ERF <= LOOCVPred$upper95 & 
                   LOOCVPred$ERF >= LOOCVPred$lower95)) / 221,
      2)

## lm() ----

# lin_mod <- lm(E_R_F ~., df)
# summary(lin_mod)$adj.r.squared

LOOCVPred_lm <- data.frame()
fittedval <- lower95 <- upper95 <- rep(NA, 221)
LOOCVPred_lm <- data.frame(cbind(fittedval, lower95, upper95))

# df <- data.frame(
#   cbind("E_R_F" = response_E_R_F_dataframe, inputs_x_norm_T_matrix)
# )
# obs <- 220
for (obs in 1:221) {
  linmod <- lm(E_R_F ~., df[-obs,])
  pred <- predict(linmod, 
                  newdata = inputs_x_norm_T_matrix[obs,],
                  interval = "confidence")
  LOOCVPred_lm[obs, 1] <- pred[,1]
  LOOCVPred_lm[obs, 2] <- pred[,2]
  LOOCVPred_lm[obs, 3] <- pred[,3]
  cat(obs, " ")
}

LOOCVPred_lm <- cbind(LOOCVPred_lm, 
                   "ERF" = response_ERF_dataframe
)
minXY <- min(LOOCVPred_lm$lower95, LOOCVPred_lm$ERF, LOOCVPred_lm$ERF)
maxXY <- max(LOOCVPred_lm$upper95, LOOCVPred_lm$ERF, LOOCVPred_lm$ERF)

par(mar = c(5.1, 4.1, 4.1, 2.1))

lmfitted_true_cols <- rep("red", 221)
lmfitted_true_cols[LOOCVPred_lm$ERF <= LOOCVPred_lm$upper95 & 
               LOOCVPred_lm$ERF >= LOOCVPred_lm$lower95] <- "black"

plot(fittedval ~ ERF, LOOCVPred_lm,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = lmfitted_true_cols
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:Nval){
  lines(c(LOOCVPred_lm$ERF[obs],LOOCVPred_lm$ERF[obs]),
        c(LOOCVPred_lm$lower95[obs],LOOCVPred_lm$upper95[obs]),lwd=1.2,
        col = lmfitted_true_cols[obs]
  )
}
round((221 - sum(LOOCVPred_lm$ERF <= LOOCVPred_lm$upper95 & 
                   LOOCVPred_lm$ERF >= LOOCVPred_lm$lower95)) / 221,
      2)

# LOO CV for AOD ----

long <- 44; lat <- 89 # 0.49 R2adj, 0.92 AM
# long <- 44; lat <- 90 # 0.25 R2adj, 0.06 AM
# long <- 43; lat <- 90 # 0.77 R2adj, 0.37 AM
response_AOD_July_gridbox <- data.frame(response_A_O_D_array[long, lat,])

AM_A_O_D_versus_E_R_F_matrix[long, lat]

library(lhs)
x_star_matrix <- t(randomLHS(N, p))
x_norm_matrix <- t(inputs_x_norm_T_matrix)
x_star_T_dataframe <- data.frame(t(x_star_matrix))
colnames(x_star_T_dataframe) <- colnames(inputs_x_norm_T_matrix)
# make H_matrix
n <- nrow(inputs_x_norm_T_matrix)
H_matrix <- cbind(c(rep(1, n)), t(unname(x_norm_matrix)))

q <- 1 + p

library(fields)
# library(fields, lib="C:/Users/smp22ijw/Desktop/Library/") # if Lenovo
corGaussian <- function(inputs, inputs2, phi) {
  
  if (missing(inputs2) || is.null(inputs2))
    return(corGaussianSquare(inputs, phi))
  
  delta <- (phi)
  exp(-(rdist(inputs / rep(delta, each = nrow(inputs)), inputs2 / rep(delta, each = nrow(inputs2))) ^ 2))
}

assign(paste0("gp_d_dxi_normalised_E_R_F_matrix"),
       matrix(as.numeric(readLines(paste0("objects/gp_d_dxi_normalised_E_R_F_matrix", "_", N, "_2.txt"))),
              nrow = p)
)

for (rerun in 1:1) {
  gp <- km(~.,
           design = inputs_x_norm_T_matrix,
           response = response_A_O_D_array[long, lat,],
           covtype="gauss", optim.method="BFGS", control=list(maxit=500))
  betas_hat_matrix <- matrix(gp@trend.coef,
                             nrow = q)
  sigma_sq_hat <- gp@covariance@sd2
  l_hat_vector <- gp@covariance@range.val
  l_hat_matrix <- matrix(gp@covariance@range.val,
                         nrow = p)
  l_hat_diag_matrix <- diag(as.vector(l_hat_matrix))
#   # assign(paste0("betas_with_", N, "_gb_", g, "_attempt_", c), gp@trend.coef)
#   # assign(paste0("l_hat_with_", N, "_gb_", g, "_attempt_", c), gp@covariance@range.val)
  x_star_predictions_list <- predict(gp,
                                     newdata = x_star_T_dataframe,
                                     type="SK"
  )
  A_matrix <-
    corGaussian(t(x_norm_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
  # make A_inv_matrix
  A_inv_matrix <- solve(A_matrix)
# 
  # partial derivatives
  partial_derivatives_dataframe <- data.frame(rep(NA, N))
  for (i in 2:p) {
    partial_derivatives_dataframe <- cbind(partial_derivatives_dataframe,
                                           data.frame(rep(NA, N)))
  }
  colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx_",p), 1:p)
#   
  # make t(x_star)^T
  t_x_star_T_matrix <- corGaussian(t(x_star_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
#   
  for (i in 1:p) {
    d_dxi_t_x_star_T_matrix <- -2 / l_hat_matrix[i,]^2 * (matrix(rep(x_star_matrix[i,], n), ncol = n) - matrix(rep(x_norm_matrix[i,], N), ncol = n, byrow = T)) * t_x_star_T_matrix[,]
    for (k in 1:N) {
      partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dxi_t_x_star_T_matrix[k,] %*% A_inv_matrix %*% (response - H_matrix %*% betas_hat_matrix)))
    }
    cat(long, "-", lat, " (", i, "/", p, ") ", sep = "")
  }
#   
  partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
#   
  partial_derivatives_normalised_transposed_matrix <- t(data.matrix(partial_derivatives_normalised_dataframe))
#   
  AM <- sum(abs(colSums(partial_derivatives_normalised_transposed_matrix *
                                                  gp_d_dxi_normalised_E_R_F_matrix))
                                        ) / N
  print(AM)
}

gp
betas_hat_matrix
sigma_sq_hat
l_hat_vector

LOOCVPred_km_AOD_44_89 <- leaveOneOut.km(gp, 
                                         "SK", 
                                         trend.reestim=FALSE)

LOOCVPred_km_AOD_44_89_for_plotting <- cbind("mean" = LOOCVPred_km_AOD_44_89$mean,
                                             "lower95" = LOOCVPred_km_AOD_44_89$mean -
                                               qnorm(0.975) * LOOCVPred_km_AOD_44_89$sd,
                                             "upper95" = LOOCVPred_km_AOD_44_89$mean +
                                               qnorm(0.975) * LOOCVPred_km_AOD_44_89$sd,
                                             response_AOD_July_gridbox
                                             )
colnames(LOOCVPred_km_AOD_44_89_for_plotting)[4] <- "AOD"
minXY <- min(LOOCVPred_km_AOD_44_89_for_plotting$lower95, 
             LOOCVPred_km_AOD_44_89_for_plotting$AOD)
maxXY <- max(LOOCVPred_km_AOD_44_89_for_plotting$upper95, 
             LOOCVPred_km_AOD_44_89_for_plotting$AOD)

par(mar = c(5.1, 4.1, 4.1, 2.1))

em_true_cols <- rep("red", 221)
em_true_cols[LOOCVPred_km_AOD_44_89_for_plotting$AOD <= LOOCVPred_km_AOD_44_89_for_plotting$upper95 & 
               LOOCVPred_km_AOD_44_89_for_plotting$AOD >= LOOCVPred_km_AOD_44_89_for_plotting$lower95] <- "black"

plot(mean ~ AOD, LOOCVPred_km_AOD_44_89_for_plotting,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = em_true_cols
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:Nval){
  lines(c(LOOCVPred_km_AOD_44_89_for_plotting$AOD[obs],
          LOOCVPred_km_AOD_44_89_for_plotting$AOD[obs]),
        c(LOOCVPred_km_AOD_44_89_for_plotting$lower95[obs],
          LOOCVPred_km_AOD_44_89_for_plotting$upper95[obs]),
        lwd=1.2, col = em_true_cols[obs]
  )
}
round((221 - sum(LOOCVPred_km_AOD_44_89_for_plotting$ERF <= LOOCVPred_km_AOD_44_89_for_plotting$upper95 & 
                   LOOCVPred_km_AOD_44_89_for_plotting$ERF >= LOOCVPred_km_AOD_44_89_for_plotting$lower95)) / 221,
      2)

round((221 - sum(LOOCVPred_km_AOD_44_89_for_plotting$AOD <= LOOCVPred_km_AOD_44_89_for_plotting$upper95 &
                   LOOCVPred_km_AOD_44_89_for_plotting$AOD >= LOOCVPred_km_AOD_44_89_for_plotting$lower95)) / 221, 
      2)







LOOCVPred <- data.frame()
mean <- lower95 <- upper95 <- rep(NA, 221)
LOOCVPred <- data.frame(cbind(mean, lower95, upper95))
for (obs in 1:221) {
  gp <- km(~., 
           design = inputs_x_norm_T_matrix[-obs,], 
           response = response_AOD_July_gridbox[-obs,], 
           covtype="gauss", optim.method="BFGS", control=list(maxit=500)
  )
  pred <- predict(object=gp,
                  newdata=data.frame(inputs_x_norm_T_matrix[obs,]),
                  type="UK",checkNames=FALSE,light.return=TRUE
  )
  LOOCVPred[obs, 1] <- pred$mean
  LOOCVPred[obs, 2] <- pred$lower95
  LOOCVPred[obs, 3] <- pred$upper95
}

df <- data.frame(
  cbind("A_O_D" = response_AOD_July_gridbox, inputs_x_norm_T_matrix)
)
colnames(df)[1] <- "A_O_D"
lin_mod <- lm(A_O_D ~., df)
# summary(lin_mod)$adj.r.squared

LOOCVPred <- cbind(LOOCVPred, 
                   "AOD" = response_AOD_July_gridbox, 
                   "lmfittedvalues" = lin_mod$fitted.values
)
minXY <- min(LOOCVPred$lower95, LOOCVPred$AOD, LOOCVPred$lmfittedvalues)
maxXY <- max(LOOCVPred$upper95, LOOCVPred$AOD, LOOCVPred$lmfittedvalues)

par(mar = c(5.1, 4.1, 4.1, 2.1))

em_true_cols <- rep("red", 221)
em_true_cols[LOOCVPred$AOD <= LOOCVPred$upper95 & 
               LOOCVPred$AOD >= LOOCVPred$lower95] <- "black"

plot(mean ~ AOD, LOOCVPred,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = em_true_cols
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:Nval){
  lines(c(LOOCVPred$AOD[obs],LOOCVPred$AOD[obs]),
        c(LOOCVPred$lower95[obs],LOOCVPred$upper95[obs]),lwd=1.2,
        col = em_true_cols[obs]
  )
}
round((221 - sum(LOOCVPred$ERF <= LOOCVPred$upper95 & 
                   LOOCVPred$ERF >= LOOCVPred$lower95)) / 221,
      2)
