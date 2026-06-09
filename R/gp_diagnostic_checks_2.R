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

## Importing input settings
inputs_x_norm_T_matrix <- read.csv("data/ACURE-UKESM1-PPE-par-norm-Design.csv")
p <- ncol(inputs_x_norm_T_matrix)

## Importing observable variable

### AOD_Total
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

## Importing unobservable variable

response_ERF_dataframe <- read.table(
  "X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/global_mean_forcing/ERF_PPE_global_jul.dat", 
  col.names = "ERF"
  )

## Importing map-required things
ACURE_P3_AOD_Total_jul_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc")
)
# ACURE_P3_AOD_Total_jul_nc <- nc_open(paste("data/ACURE_P3_AOD_Total_jul.nc")) # if Lenovo
longitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"longitude") - 180
latitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"latitude")
rm("ACURE_P3_AOD_Total_jul_nc")

# LOO CV for ERF ----

## km() ----

gp_ERF <- km(~., 
             design = inputs_x_norm_T_matrix, 
             response = response_ERF_dataframe,
             covtype="gauss", optim.method="BFGS", 
             control=list(maxit=500)
)

LOOCVPred_km_ERF <- leaveOneOut.km(gp_ERF, 
                                         "SK", 
                                         trend.reestim=FALSE)

LOOCVPred_km_ERF_for_plotting <- cbind("mean" = LOOCVPred_km_ERF$mean,
                                       "lower95" = LOOCVPred_km_ERF$mean -
                                               qnorm(0.975) * LOOCVPred_km_ERF$sd,
                                       "upper95" = LOOCVPred_km_ERF$mean +
                                               qnorm(0.975) * LOOCVPred_km_ERF$sd,
                                       response_ERF_dataframe
)
minXY <- min(LOOCVPred_km_ERF_for_plotting$lower95, 
             LOOCVPred_km_ERF_for_plotting$ERF)
maxXY <- max(LOOCVPred_km_ERF_for_plotting$upper95, 
             LOOCVPred_km_ERF_for_plotting$ERF)

par(mar = c(5.1, 4.1, 4.1, 2.1))

cols_LOOCVPred_km_ERF <- rep("red", 221)
cols_LOOCVPred_km_ERF[LOOCVPred_km_ERF_for_plotting$ERF <= LOOCVPred_km_ERF_for_plotting$upper95 & 
               LOOCVPred_km_ERF_for_plotting$ERF >= LOOCVPred_km_ERF_for_plotting$lower95] <- "black"

plot(mean ~ ERF, LOOCVPred_km_ERF_for_plotting,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = cols_LOOCVPred_km_ERF
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:221){
  lines(c(LOOCVPred_km_ERF_for_plotting$ERF[obs],
          LOOCVPred_km_ERF_for_plotting$ERF[obs]),
        c(LOOCVPred_km_ERF_for_plotting$lower95[obs],
          LOOCVPred_km_ERF_for_plotting$upper95[obs]),
        lwd=1.2, col = cols_LOOCVPred_km_ERF[obs]
  )
}
round((221 - sum(LOOCVPred_km_ERF_for_plotting$ERF <= LOOCVPred_km_ERF_for_plotting$upper95 & 
                   LOOCVPred_km_ERF_for_plotting$ERF >= LOOCVPred_km_ERF_for_plotting$lower95)) / 221,
      2)

## lm() ----

ERF_df <- data.frame(
  cbind("E_R_F" = response_ERF_dataframe, inputs_x_norm_T_matrix)
)
ERF_lin_mod <- lm(ERF ~., ERF_df)
summary(ERF_lin_mod)$adj.r.squared

LOOCVPred_lm_ERF <- data.frame()
fittedval <- lower95 <- upper95 <- rep(NA, 221)
LOOCVPred_lm_ERF <- data.frame(cbind(fittedval, lower95, upper95))

for (obs in 1:221) {
  linmod <- lm(ERF ~., ERF_df[-obs,])
  pred <- predict(linmod, 
                  newdata = inputs_x_norm_T_matrix[obs,],
                  interval = "confidence")
  LOOCVPred_lm_ERF[obs, 1] <- pred[,1]
  LOOCVPred_lm_ERF[obs, 2] <- pred[,2]
  LOOCVPred_lm_ERF[obs, 3] <- pred[,3]
  cat(obs, " ")
}

LOOCVPred_lm_ERF_for_plotting <- cbind(LOOCVPred_lm_ERF,
                                       "ERF" = response_ERF_dataframe
)
minXY <- min(LOOCVPred_lm_ERF_for_plotting$lower95, 
             LOOCVPred_lm_ERF_for_plotting$ERF)
maxXY <- max(LOOCVPred_lm_ERF_for_plotting$upper95,
             LOOCVPred_lm_ERF_for_plotting$ERF)

par(mar = c(5.1, 4.1, 4.1, 2.1))

cols_LOOCVPred_lm_ERF <- rep("red", 221)
cols_LOOCVPred_lm_ERF[LOOCVPred_lm_ERF_for_plotting$ERF <= 
                        LOOCVPred_lm_ERF_for_plotting$upper95 & 
                        LOOCVPred_lm_ERF_for_plotting$ERF >= 
                        LOOCVPred_lm_ERF_for_plotting$lower95] <- "black"

plot(fittedval ~ ERF, LOOCVPred_lm_ERF_for_plotting,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = cols_LOOCVPred_lm_ERF
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:221){
  lines(c(LOOCVPred_lm_ERF_for_plotting$ERF[obs],
          LOOCVPred_lm_ERF_for_plotting$ERF[obs]),
        c(LOOCVPred_lm_ERF_for_plotting$lower95[obs],
          LOOCVPred_lm_ERF_for_plotting$upper95[obs]),
        lwd=1.2, col = cols_LOOCVPred_lm_ERF[obs]
  )
}
round((221 - sum(LOOCVPred_lm_ERF_for_plotting$ERF <= 
                   LOOCVPred_lm_ERF_for_plotting$upper95 & 
                   LOOCVPred_lm_ERF_for_plotting$ERF >= 
                   LOOCVPred_lm_ERF_for_plotting$lower95)) / 221,
      2)


# LOO CV for AOD ----

## getting ERF ready for AM ----
q <- 1 + p
betas_hat_matrix <- matrix(gp_ERF@trend.coef,
                           nrow = q)
sigma_sq_hat <- gp_ERF@covariance@sd2
l_hat_vector <- gp_ERF@covariance@range.val
l_hat_matrix <- matrix(gp_ERF@covariance@range.val,
                       nrow = p)
l_hat_diag_matrix <- diag(as.vector(l_hat_matrix))
library(lhs)
N <- 100000
x_star_matrix <- t(randomLHS(N, p))
x_norm_matrix <- t(inputs_x_norm_T_matrix)
x_star_T_dataframe <- data.frame(t(x_star_matrix))
colnames(x_star_T_dataframe) <- colnames(inputs_x_norm_T_matrix)
n <- nrow(inputs_x_norm_T_matrix)
H_matrix <- cbind(c(rep(1, n)), t(unname(x_norm_matrix)))
x_star_predictions_list <- predict(gp_ERF,
                                   newdata = x_star_T_dataframe,
                                   type="SK"
)
library(fields)
corGaussian <- function(inputs, inputs2, phi) {
  
  if (missing(inputs2) || is.null(inputs2))
    return(corGaussianSquare(inputs, phi))
  
  delta <- (phi)
  exp(-(rdist(inputs / rep(delta, each = nrow(inputs)), inputs2 / rep(delta, each = nrow(inputs2))) ^ 2))
}
A_matrix <-
  corGaussian(t(x_norm_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
A_inv_matrix <- solve(A_matrix)
partial_derivatives_dataframe <- data.frame(rep(NA, N))
for (i in 2:p) {
  partial_derivatives_dataframe <- cbind(partial_derivatives_dataframe, 
                                         data.frame(rep(NA, N)))
}
colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx_",p), 1:p)
t_x_star_T_matrix <- corGaussian(t(x_star_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
for (i in 1:p) {
  d_dxi_t_x_star_T_matrix <- -2 / l_hat_matrix[i,]^2 * (matrix(rep(x_star_matrix[i,], n), ncol = n) - matrix(rep(x_norm_matrix[i,], N), ncol = n, byrow = T)) * t_x_star_T_matrix[,]
  for (k in 1:N) {
    partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dxi_t_x_star_T_matrix[k,] %*% A_inv_matrix %*% (response_ERF_dataframe[,1] - H_matrix %*% betas_hat_matrix)))
  }
  cat(i, "/", p, " ", sep = "")
}
ERF_partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
ERF_partial_derivatives_normalised_transposed_matrix <- t(data.matrix(ERF_partial_derivatives_normalised_dataframe))

## km() ----

# long <- 44; lat <- 89 # 0.49 R2adj, 0.92 AM
long <- 44; lat <- 90 # 0.25 R2adj, 0.06 AM
# long <- 43; lat <- 90 # 0.77 R2adj, 0.37 AM
response_AOD_July_44_90 <- data.frame(response_A_O_D_array[long, lat,])

### check AM ----

threshold <- 0.6
assign(paste0("AM_A_O_D_versus_E_R_F_matrix"),
       matrix(as.numeric(readLines(paste0("objects/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N, ".txt"))),
              nrow = length(longitude))
)
AM_A_O_D_versus_E_R_F_matrix[long, lat]

gp_AOD_44_90 <- km(~.,
         design = inputs_x_norm_T_matrix,
         response = response_A_O_D_array[long, lat,],
         covtype="gauss", optim.method="BFGS", control=list(maxit=500))
betas_hat_matrix <- matrix(gp_AOD_44_90@trend.coef,
                           nrow = q)
sigma_sq_hat <- gp_AOD_44_90@covariance@sd2
l_hat_vector <- gp_AOD_44_90@covariance@range.val
l_hat_matrix <- matrix(gp_AOD_44_90@covariance@range.val,
                       nrow = p)
l_hat_diag_matrix <- diag(as.vector(l_hat_matrix))
x_star_predictions_list <- predict(gp_AOD_44_90,
                                   newdata = x_star_T_dataframe,
                                   type="SK"
)
A_matrix <-
  corGaussian(t(x_norm_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
A_inv_matrix <- solve(A_matrix)
partial_derivatives_dataframe <- data.frame(rep(NA, N))
for (i in 2:p) {
  partial_derivatives_dataframe <- cbind(partial_derivatives_dataframe,
                                         data.frame(rep(NA, N)))
}
colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx_",p), 1:p)
t_x_star_T_matrix <- corGaussian(t(x_star_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
for (i in 1:p) {
  d_dxi_t_x_star_T_matrix <- -2 / l_hat_matrix[i,]^2 * (matrix(rep(x_star_matrix[i,], n), ncol = n) - matrix(rep(x_norm_matrix[i,], N), ncol = n, byrow = T)) * t_x_star_T_matrix[,]
  for (k in 1:N) {
    partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dxi_t_x_star_T_matrix[k,] %*% A_inv_matrix %*% (response_AOD_July_44_90[,1] - H_matrix %*% betas_hat_matrix)))
  }
  cat(long, "-", lat, " (", i, "/", p, ") ", sep = "")
}
partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
partial_derivatives_normalised_transposed_matrix <- t(data.matrix(partial_derivatives_normalised_dataframe))
AM <- sum(abs(colSums(partial_derivatives_normalised_transposed_matrix *
                        ERF_partial_derivatives_normalised_transposed_matrix))
) / N
AM


### do the LOO CV ----

LOOCVPred_km_AOD_44_90 <- leaveOneOut.km(gp_AOD_44_90, 
                                   "SK", 
                                   trend.reestim=FALSE)

LOOCVPred_km_AOD_44_90_for_plotting <- data.frame(cbind("mean" = LOOCVPred_km_AOD_44_90$mean,
                                       "lower95" = LOOCVPred_km_AOD_44_90$mean -
                                         qnorm(0.975) * LOOCVPred_km_AOD_44_90$sd,
                                       "upper95" = LOOCVPred_km_AOD_44_90$mean +
                                         qnorm(0.975) * LOOCVPred_km_AOD_44_90$sd,
                                       "AOD" = response_A_O_D_array[long, lat,])
)
minXY <- min(LOOCVPred_km_AOD_44_90_for_plotting$lower95, 
             LOOCVPred_km_AOD_44_90_for_plotting$AOD)
maxXY <- max(LOOCVPred_km_AOD_44_90_for_plotting$upper95, 
             LOOCVPred_km_AOD_44_90_for_plotting$AOD)

par(mar = c(5.1, 4.1, 4.1, 2.1))

cols_LOOCVPred_km_AOD_44_90 <- rep("red", 221)
cols_LOOCVPred_km_AOD_44_90[
  LOOCVPred_km_AOD_44_90_for_plotting$AOD <= LOOCVPred_km_AOD_44_90_for_plotting$upper95 &
    LOOCVPred_km_AOD_44_90_for_plotting$AOD >= LOOCVPred_km_AOD_44_90_for_plotting$lower95] <- "black"

plot(mean ~ AOD, LOOCVPred_km_AOD_44_90_for_plotting,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = cols_LOOCVPred_km_AOD_44_90
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:221){
  lines(c(LOOCVPred_km_AOD_44_90_for_plotting$AOD[obs],
          LOOCVPred_km_AOD_44_90_for_plotting$AOD[obs]),
        c(LOOCVPred_km_AOD_44_90_for_plotting$lower95[obs],
          LOOCVPred_km_AOD_44_90_for_plotting$upper95[obs]),
        lwd=1.2, col = cols_LOOCVPred_km_AOD_44_90[obs]
  )
}
round((221 - sum(
  LOOCVPred_km_AOD_44_90_for_plotting$AOD <= LOOCVPred_km_AOD_44_90_for_plotting$upper95 &
    LOOCVPred_km_AOD_44_90_for_plotting$AOD >= LOOCVPred_km_AOD_44_90_for_plotting$lower95)) / 221,
      2)

## lm() ----

AOD_44_90_df <- data.frame(
  cbind("AOD" = response_A_O_D_array[long, lat,], inputs_x_norm_T_matrix)
)
AOD_44_90_lin_mod <- lm(AOD ~., AOD_44_90_df)
summary(AOD_44_90_lin_mod)$adj.r.squared

LOOCVPred_lm_AOD_44_90 <- data.frame()
fittedval <- lower95 <- upper95 <- rep(NA, 221)
LOOCVPred_lm_AOD_44_90 <- data.frame(cbind(fittedval, lower95, upper95))

for (obs in 1:221) {
  linmod <- lm(AOD ~., AOD_44_90_df[-obs,])
  pred <- predict(AOD_44_90_lin_mod, 
                  newdata = inputs_x_norm_T_matrix[obs,],
                  interval = "confidence")
  LOOCVPred_lm_AOD_44_90[obs, 1] <- pred[,1]
  LOOCVPred_lm_AOD_44_90[obs, 2] <- pred[,2]
  LOOCVPred_lm_AOD_44_90[obs, 3] <- pred[,3]
  cat(obs, " ")
}

LOOCVPred_lm_AOD_44_90_for_plotting <- cbind(LOOCVPred_lm_AOD_44_90,
                                       "AOD" = response_A_O_D_array[long, lat,]
                                       )
minXY <- min(LOOCVPred_lm_AOD_44_90_for_plotting$lower95, 
             LOOCVPred_lm_AOD_44_90_for_plotting$AOD)
maxXY <- max(LOOCVPred_lm_AOD_44_90_for_plotting$upper95,
             LOOCVPred_lm_AOD_44_90_for_plotting$AOD)

par(mar = c(5.1, 4.1, 4.1, 2.1))

cols_LOOCVPred_lm_AOD_44_90 <- rep("red", 221)
cols_LOOCVPred_lm_AOD_44_90[LOOCVPred_lm_AOD_44_90_for_plotting$AOD <= 
                        LOOCVPred_lm_AOD_44_90_for_plotting$upper95 & 
                        LOOCVPred_lm_AOD_44_90_for_plotting$AOD >= 
                        LOOCVPred_lm_AOD_44_90_for_plotting$lower95] <- "black"

plot(fittedval ~ AOD, LOOCVPred_lm_AOD_44_90_for_plotting,
     xlim = c(minXY, maxXY), ylim = c(minXY, maxXY),
     pch=20, col = cols_LOOCVPred_lm_AOD_44_90
)
abline(0,1, col = "blue", lwd=2)
for(obs in 1:221){
  lines(c(LOOCVPred_lm_AOD_44_90_for_plotting$AOD[obs],
          LOOCVPred_lm_AOD_44_90_for_plotting$AOD[obs]),
        c(LOOCVPred_lm_AOD_44_90_for_plotting$lower95[obs],
          LOOCVPred_lm_AOD_44_90_for_plotting$upper95[obs]),
        lwd=1.2, col = cols_LOOCVPred_lm_AOD_44_90[obs]
  )
}
round((221 - sum(LOOCVPred_lm_AOD_44_90_for_plotting$AOD <= 
                   LOOCVPred_lm_AOD_44_90_for_plotting$upper95 & 
                   LOOCVPred_lm_AOD_44_90_for_plotting$AOD >= 
                   LOOCVPred_lm_AOD_44_90_for_plotting$lower95)) / 221,
      2)


