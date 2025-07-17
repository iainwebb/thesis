# Setup ----

# replace back in
# ACURE_P3_AOD_Total_jul_level2 for A_O_D
# ERF_PPE_global_jul for E_R_F
options(scipen=999)

library(readr)
library(DiceKriging)
# https://colordesigner.io/gradient-generator

# Importing input settings ----
inputs_x_norm_T_matrix <- read.csv("data/ACURE-UKESM1-PPE-par-norm-Design.csv")
p <- ncol(inputs_x_norm_T_matrix)

# Importing observable variable ----

## AOD_Total
library(ncdf4)
ACURE_P3_AOD_Total_jul_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc")
  )
ACURE_P3_AOD_Total_jul_array <- ncvar_get(
  ACURE_P3_AOD_Total_jul_nc, names(ACURE_P3_AOD_Total_jul_nc$var)
  )
response_A_O_D_array <- ACURE_P3_AOD_Total_jul_array[,,2,]
rm("ACURE_P3_AOD_Total_jul_nc", "ACURE_P3_AOD_Total_jul_array")

# Importing unobservable variable ----

response_E_R_F_dataframe <- read.table(
  "X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/global_mean_forcing/ERF_PPE_global_jul.dat", 
  col.names = "E_R_F"
  )

# Importing map-required things ----
ACURE_P3_AOD_Total_jul_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc")
  )
longitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"longitude") - 180
latitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"latitude")
rm("ACURE_P3_AOD_Total_jul_nc")

# Calculating R2_adj and lm partial deriv's for all gridboxes, then lm partial deriv's for lm'able, for both variables ----

threshold <- 0.60

## AOD_Total

### from scratch

# R2_adj_A_O_D_matrix <- matrix(NA, nrow = 192, ncol = 144)
# lm_all_d_dxi_normalised_A_O_D_array <- array(NA, dim = c(length(longitude), length(latitude), p))
# lm_best_d_dxi_normalised_A_O_D_array <- array(NA, dim = c(length(longitude), length(latitude), p))
# 
# for (long in 1:192) {
#   for (lat in 1:144) {
#     df <- data.frame(
#       cbind("A_O_D" = response_A_O_D_array[long, lat, ], inputs_x_norm_T_matrix)
#       )
#     lin_mod <- lm(A_O_D ~., df)
#     R2_adj_A_O_D_matrix[long, lat] <-
#       summary(lin_mod)$adj.r.squared
#     lm_all_d_dxi_normalised_A_O_D_array[long, lat,] <-
#       as.vector(summary(lin_mod)$coefficients[2:38, 1]) / sqrt(sum((as.vector(summary(lin_mod)$coefficients[2:38, 1]))^2))
#     lm_best_d_dxi_normalised_A_O_D_array[long, lat,] <-
#       if (summary(lin_mod)$adj.r.squared >= threshold) lm_all_d_dxi_normalised_A_O_D_array[long, lat,] 
#       else rep(NA, p)
#   }
#   cat(long, "/192 ", sep = "")
# }
# rm("long", "lat", "df", "lin_mod")
# write_lines(as.vector(R2_adj_A_O_D_matrix),
#             file="object/R2_adj_A_O_D_matrix.txt")
# write_lines(as.vector(lm_all_d_dxi_normalised_A_O_D_array),
#             file="object/lm_all_d_dxi_normalised_A_O_D_array.txt")
# write_lines(as.vector(lm_best_d_dxi_normalised_A_O_D_array),
#             file=paste0("object/lm_best_d_dxi_normalised_A_O_D_array_0", threshold*100, ".txt"))

### from before

R2_adj_A_O_D_matrix <-
  matrix(as.numeric(readLines("object/R2_adj_A_O_D_matrix.txt")),
         nrow = length(longitude))

lm_all_d_dxi_normalised_A_O_D_array <-
  array(as.numeric(readLines("object/lm_all_d_dxi_normalised_A_O_D_array.txt")),
        dim = c(length(longitude), length(latitude), p))

assign(paste0("lm_best_d_dxi_normalised_A_O_D_array"),
       array(as.numeric(readLines(paste0("object/lm_best_d_dxi_normalised_A_O_D_array_0", threshold*100, ".txt"))),
             dim = c(length(longitude), length(latitude), p))
       )
sum(is.na(lm_best_d_dxi_normalised_A_O_D_array)) / 37 # check

### Frequency table
table(
  ifelse(floor(R2_adj_A_O_D_matrix * 10) / 10 < 0.6, "<0.6",
         floor(R2_adj_A_O_D_matrix * 10) / 10)
)
### Table with percentages
100*round(
  table(
    ifelse(floor(R2_adj_A_O_D_matrix * 10) / 10 < 0.6, "<0.6",
           floor(R2_adj_A_O_D_matrix * 10) / 10)
  ) / (192*144), 3)
hist(R2_adj_A_O_D_matrix)
### Map
library(fields)
library(maps)
image.plot(longitude, latitude,
           rbind(R2_adj_A_O_D_matrix[97:192,],
                 R2_adj_A_O_D_matrix[1:96,]),
           breaks = c(0,0.6,0.7,0.75,0.8,1),
           col = c("blue","green","orange","yellow", "white"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.6,0.7,0.75,0.8,1),labels=as.character(c(0,0.6,0.7,0.75,0.8,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1)

## ERF_global

### from scratch

# R2_adj_E_R_F_matrix <- matrix(NA, nrow = 192, ncol = 144)
# lm_d_dxi_normalised_E_R_F_array <- array(NA, dim = c(length(longitude), length(latitude), p))
# 
# for (long in 1:192) {
#   for (lat in 1:144) {
#     df <- data.frame(
#       cbind("E_R_F" = response_E_R_F_dataframe, inputs_x_norm_T_matrix)
#     )
#     lin_mod <- lm(E_R_F ~., df)
#     R2_adj_E_R_F_matrix[long, lat] <- summary(lin_mod)$adj.r.squared
#     lm_d_dxi_normalised_E_R_F_array[long, lat,] <-
#       if (summary(lin_mod)$adj.r.squared >= threshold) as.vector(summary(lin_mod)$coefficients[2:38, 1]) / sqrt(sum((as.vector(summary(lin_mod)$coefficients[2:38, 1]))^2)) else NA
#   }
#   cat(long, "/192 ", sep = "")
# }
# rm("long", "lat", "df", "lin_mod")
# write_lines(as.vector(R2_adj_E_R_F_matrix),
#             file="object/R2_adj_E_R_F_matrix.txt")
# write_lines(as.vector(lm_d_dxi_normalised_E_R_F_array),
#             file="object/lm_d_dxi_normalised_E_R_F_array.txt")
# 
# sum(is.na(lm_d_dxi_normalised_E_R_F_array)) / 37 # check
# 192 * 144

### from before

R2_adj_E_R_F_matrix <-
  matrix(as.numeric(readLines("object/R2_adj_E_R_F_matrix.txt")),
         nrow = length(longitude))

lm_d_dxi_normalised_E_R_F_array <-
  array(as.numeric(readLines("object/lm_d_dxi_normalised_E_R_F_array.txt")),
        dim = c(length(longitude), length(latitude), p))

sum(is.na(lm_d_dxi_normalised_E_R_F_array)) / 37 # check

# GP work for use in AM work ----

N <- 100000
q <- 1 + p
n <- nrow(inputs_x_norm_T_matrix)

library(lhs)
x_star_matrix <- t(randomLHS(N, p))
x_norm_matrix <- t(inputs_x_norm_T_matrix)
x_star_T_dataframe <- data.frame(t(x_star_matrix))
colnames(x_star_T_dataframe) <- colnames(inputs_x_norm_T_matrix)
# make H_matrix
H_matrix <- cbind(c(rep(1, n)), t(unname(x_norm_matrix)))

library(fields)
corGaussian <- function(inputs, inputs2, phi) {
  
  if (missing(inputs2) || is.null(inputs2))
    return(corGaussianSquare(inputs, phi))
  
  delta <- (phi)
  exp(-(rdist(inputs / rep(delta, each = nrow(inputs)), inputs2 / rep(delta, each = nrow(inputs2))) ^ 2))
}

gp_d_dxi_normalised_function <- function(response) {
  gp <- km(~., 
           design = inputs_x_norm_T_matrix, 
           response = response, 
           covtype="gauss", optim.method="BFGS", control=list(maxit=500))
  betas_hat_matrix <- matrix(gp@trend.coef,
                             nrow = q)
  sigma_sq_hat <- gp@covariance@sd2
  l_hat_vector <- gp@covariance@range.val
  l_hat_matrix <- matrix(gp@covariance@range.val,
                         nrow = p)
  l_hat_diag_matrix <- diag(as.vector(l_hat_matrix))
  # assign(paste0("betas_with_", N, "_gb_", g, "_attempt_", c), gp@trend.coef)
  # assign(paste0("l_hat_with_", N, "_gb_", g, "_attempt_", c), gp@covariance@range.val)
  x_star_predictions_list <- predict(gp,
                                     newdata = x_star_T_dataframe,
                                     type="SK"
  )
  A_matrix <-
    corGaussian(t(x_norm_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
  # make A_inv_matrix
  A_inv_matrix <- solve(A_matrix)
  
  # partial derivatives
  partial_derivatives_dataframe <- data.frame(rep(NA, N))
  for (i in 2:p) {
    partial_derivatives_dataframe <- cbind(partial_derivatives_dataframe, 
                                           data.frame(rep(NA, N)))
  }
  colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx_",p), 1:p)
  
  # make t(x_star)^T
  t_x_star_T_matrix <- corGaussian(t(x_star_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))
  
  for (i in 1:p) {
    d_dxi_t_x_star_T_matrix <- -2 / l_hat_matrix[i,]^2 * (matrix(rep(x_star_matrix[i,], n), ncol = n) - matrix(rep(x_norm_matrix[i,], N), ncol = n, byrow = T)) * t_x_star_T_matrix[,]
    for (k in 1:N) {
      partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dxi_t_x_star_T_matrix[k,] %*% A_inv_matrix %*% (response - H_matrix %*% betas_hat_matrix)))
    }
    cat(long, "-", lat, " (", i, "/", p, ") ", sep = "")
  }
  
  partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
  
  partial_derivatives_normalised_transposed_matrix <- t(data.matrix(partial_derivatives_normalised_dataframe))
  
  return(partial_derivatives_normalised_transposed_matrix)
}

# 4-case AM calculation ----

# from scratch

# AM_A_O_D_versus_E_R_F_matrix_case0 <- matrix(NA, nrow = 192, ncol = 144)

# export
# write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix_case0),
#             file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_case0_0", threshold*100, "_", N, ".txt"))

# import
assign(paste0("AM_A_O_D_versus_E_R_F_matrix_case0"),
       matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_case0_0", threshold*100, "_", N, ".txt"))),
              nrow = length(longitude))
)
sum(is.na(AM_A_O_D_versus_E_R_F_matrix_case0))

## Case 1: R2_adj >= threshold for both variables

### from scratch

# AM_A_O_D_versus_E_R_F_matrix_case1 <- AM_A_O_D_versus_E_R_F_matrix_case0

# sum(R2_adj_A_O_D_matrix >= threshold & R2_adj_E_R_F_matrix >= threshold)
 
# for (long in 1:192) {
#   for (lat in 1:144) {
#     AM_A_O_D_versus_E_R_F_matrix_case1[long, lat] <-
#       if (R2_adj_A_O_D_matrix[long, lat] >= threshold & R2_adj_E_R_F_matrix[long, lat] >= threshold)
#         abs(sum(lm_best_d_dxi_normalised_A_O_D_array[long, lat,] *
#                           lm_d_dxi_normalised_E_R_F_array[long, lat,]))
#       else NA
#     if (R2_adj_A_O_D_matrix[long, lat] >= threshold & R2_adj_E_R_F_matrix[long, lat] >= threshold)
#       write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix_case1),
#                   file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_case1_0", threshold*100, "_", N, ".txt"))
#   }
#   cat(long, "/192 ", sep = "")
# }

# export
# write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix_case1),
#             file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_case1_0", threshold*100, "_", N, ".txt"))

# import
assign(paste0("AM_A_O_D_versus_E_R_F_matrix_case1"),
       matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_case1_0", threshold*100, "_", N, ".txt"))),
              nrow = length(longitude))
)

sum(is.na(AM_A_O_D_versus_E_R_F_matrix_case1))
image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_matrix_case1[97:192,],
                 AM_A_O_D_versus_E_R_F_matrix_case1[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

## Case 2: R2_adj < threshold for observable, >= threshold for unobservable

### from scratch

# AM_A_O_D_versus_E_R_F_matrix_case2 <- AM_A_O_D_versus_E_R_F_matrix_case1

# sum(R2_adj_A_O_D_matrix < threshold & R2_adj_E_R_F_matrix >= threshold)

# for (long in 1:192) {
#   for (lat in 1:144) {
#     AM_A_O_D_versus_E_R_F_matrix_case2[long, lat] <- 
#       if (R2_adj_A_O_D_matrix[long, lat] < threshold & R2_adj_E_R_F_matrix[long, lat] >= threshold)
#         sum(abs(colSums(gp_d_dxi_normalised_function(response_A_O_D_array[long, lat,]) * 
#                   matrix(rep(lm_d_dxi_normalised_E_R_F_array[long, lat,], N), nrow = p, byrow = F)))
#         ) / N
#       else AM_A_O_D_versus_E_R_F_matrix_case2[long, lat]
#     if (R2_adj_A_O_D_matrix[long, lat] < threshold & R2_adj_E_R_F_matrix[long, lat] >= threshold)
#       write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix_case2),
#                   file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_case2_0", threshold*100, "_", N, ".txt"))
#   }
#   cat(long, "/192 ", sep = "")
# }

# export
# write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix_case2),
#             file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_case2_0", threshold*100, "_", N, ".txt"))

# import
assign(paste0("AM_A_O_D_versus_E_R_F_matrix_case2"),
       matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_case2_0", threshold*100, "_", N, ".txt"))),
              nrow = length(longitude))
)
sum(is.na(AM_A_O_D_versus_E_R_F_matrix_case2))

image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_matrix_case2[97:192,],
                 AM_A_O_D_versus_E_R_F_matrix_case2[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

## Case 3: R2_adj >= threshold for observable, < threshold for unobservable

### from scratch

# AM_A_O_D_versus_E_R_F_matrix_case3<- AM_A_O_D_versus_E_R_F_matrix_case2

# sum(R2_adj_A_O_D_matrix >= threshold & R2_adj_E_R_F_matrix < threshold)

### from scratch

# gp_d_dxi_normalised_E_R_F_matrix <- gp_d_dxi_normalised_function(response_E_R_F_dataframe[,1])

# export
# write_lines(as.vector(gp_d_dxi_normalised_E_R_F_matrix),
            # file=paste0("object/gp_d_dxi_normalised_E_R_F_matrix", "_", N, ".txt"))

### import
assign(paste0("gp_d_dxi_normalised_E_R_F_matrix"),
       matrix(as.numeric(readLines(paste0("object/gp_d_dxi_normalised_E_R_F_matrix", "_", N, ".txt"))),
              nrow = p)
)

# for (long in 1:192) {
#   for (lat in 1:144) {
#     AM_A_O_D_versus_E_R_F_matrix_case3[long, lat] <- 
#       if (R2_adj_A_O_D_matrix[long, lat] >= threshold & R2_adj_E_R_F_matrix[long, lat] < threshold)
#         sum(abs(colSums(matrix(rep(lm_best_d_dxi_normalised_A_O_D_array[long, lat,], N), nrow = p, byrow = F) *
#                           gp_d_dxi_normalised_E_R_F_matrix))
#         ) / N
#       else AM_A_O_D_versus_E_R_F_matrix_case3[long, lat]
#     }
#   cat(long, "/192 ", sep = "")
# }

# export
# write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix_case3),
#             file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_case3_0", threshold*100, "_", N, ".txt"))


### import
assign(paste0("AM_A_O_D_versus_E_R_F_matrix_case3"),
       matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_case3_0", threshold*100, "_", N, ".txt"))),
              nrow = length(longitude))
)

sum(is.na(AM_A_O_D_versus_E_R_F_matrix_case3)) # check

image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_matrix_case3[97:192,],
                 AM_A_O_D_versus_E_R_F_matrix_case3[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

## Case 4: R2_adj < threshold for observable, < threshold for unobservable

### from scratch

# AM_A_O_D_versus_E_R_F_matrix <- AM_A_O_D_versus_E_R_F_matrix_case3

# sum(R2_adj_A_O_D_matrix < threshold & R2_adj_E_R_F_matrix < threshold)

# index_counter <- 0

# for (long in 1:192) {
#   for (lat in 1:144) {
#     AM_A_O_D_versus_E_R_F_matrix[long, lat] <-
#       if (R2_adj_A_O_D_matrix[long, lat] < threshold & R2_adj_E_R_F_matrix[long, lat] < threshold)
#         sum(abs(colSums(gp_d_dxi_normalised_function(response_A_O_D_array[long, lat,]) *
#                           gp_d_dxi_normalised_E_R_F_matrix))
#         ) / N
#     else AM_A_O_D_versus_E_R_F_matrix[long, lat]
#     index_counter <-
#       if (R2_adj_A_O_D_matrix[long, lat] < threshold & R2_adj_E_R_F_matrix[long, lat] < threshold)
#         index_counter + 1
#     else index_counter
#     if (R2_adj_A_O_D_matrix[long, lat] < threshold & R2_adj_E_R_F_matrix[long, lat] < threshold)
#       write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix),
#                   file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N, ".txt"))
#   }
# }
# rm("long", "lat", "p", "q", "x_star_matrix", "x_norm_matrix", "x_star_T_dataframe", 
#    "n", "H_matrix", "corGaussian", "gp_d_dxi_normalised_function", "index_counter")

# export
# write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix),
#             file=paste0("object/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N, ".txt"))

### import
AM_A_O_D_versus_E_R_F_matrix <-
  matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N, ".txt"))),
         nrow = length(longitude))
sum(is.na(AM_A_O_D_versus_E_R_F_matrix)) # check

range(AM_A_O_D_versus_E_R_F_matrix)
hist(AM_A_O_D_versus_E_R_F_matrix)

image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_matrix[97:192,],
                 AM_A_O_D_versus_E_R_F_matrix[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

# First scatterplots

# AM_A_O_D_versus_E_R_F_all_lm_matrix <- matrix(NA, nrow = 192, ncol = 144)
# assign(paste0("gp_d_dxi_normalised_E_R_F_matrix"),
#        matrix(as.numeric(readLines(paste0("object/gp_d_dxi_normalised_E_R_F_matrix", "_", N, ".txt"))),
#               nrow = p)
# )
# for (long in 1:192) {
#   for (lat in 1:144) {
#     AM_A_O_D_versus_E_R_F_all_lm_matrix[long, lat] <-
#         sum(abs(colSums(matrix(rep(lm_all_d_dxi_normalised_A_O_D_array[long, lat,], N), nrow = p, byrow = F) *
#                           gp_d_dxi_normalised_E_R_F_matrix))
#         ) / N
#   }
#   cat(long, "/192, ", sep = "")
# }

# # export
# write_lines(as.vector(AM_A_O_D_versus_E_R_F_all_lm_matrix),
#             file=paste0("object/AM_A_O_D_versus_E_R_F_all_lm_matrix_", N, ".txt"))

# import
assign(paste0("AM_A_O_D_versus_E_R_F_all_lm_matrix"),
       matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_all_lm_matrix.txt"))),
              nrow = length(longitude))
)

R2_adj_AM_all_lm <- cbind("AOD_Total_R2_adj" 
                          = as.vector(R2_adj_A_O_D_matrix), 
                          "AM_A_O_D_versus_E_R_F_all_lm_matrix" 
                          = as.vector(AM_A_O_D_versus_E_R_F_all_lm_matrix)
)
blue_black_colours <- ifelse(as.vector(R2_adj_A_O_D_matrix) < 0.60, "blue", "black")
plot(R2_adj_AM_all_lm, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours,
     xlab = "R2_adj for AOD_Total",
     ylab = "AM AOD_Total versus ERF")
R2_adj_AM <- cbind("AOD_Total_R2_adj" 
                   = as.vector(R2_adj_A_O_D_matrix), 
                   "AM_for_AOD_Total_versus_ERF" 
                   = as.vector(AM_A_O_D_versus_E_R_F_matrix)
)
plot(R2_adj_AM, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours,
     xlab = "R2_adj for AOD_Total",
     ylab = "AM AOD_Total versus ERF")

rm(AM_A_O_D_versus_E_R_F_matrix)

# Tracking two below-threshold gridboxes

N <- 100000

AM_A_O_D_versus_E_R_F_matrix <-
  matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N, ".txt"))),
         nrow = length(longitude))
sum(is.na(AM_A_O_D_versus_E_R_F_matrix)) # check

## Map
AM_A_O_D_versus_E_R_F_only_gp_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_A_O_D_versus_E_R_F_only_gp_matrix[long, lat] <- 
      if (R2_adj_A_O_D_matrix[long, lat] < 0.6) AM_A_O_D_versus_E_R_F_matrix[long, lat]
    else NA
  }
  cat(long, "/192, ", sep = "")
}

image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_only_gp_matrix[97:192,],
                 AM_A_O_D_versus_E_R_F_only_gp_matrix[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

# Selecting the two gridboxes
double_gp_indices <- which(R2_adj_A_O_D_matrix < threshold & R2_adj_E_R_F_matrix < threshold, arr.ind = TRUE)
double_gp_indices_dataframe <- as.data.frame(double_gp_indices)
double_gp_indices_sorted_dataframe <- double_gp_indices_dataframe[order(double_gp_indices_dataframe[,1]),]
rownames(double_gp_indices_sorted_dataframe) <- 1:nrow(double_gp_indices_dataframe)
double_gp_indices_sorted_dataframe

range(double_gp_indices_sorted_dataframe[,2])
hist(double_gp_indices_sorted_dataframe[,1])
double_gp_indices_sorted_dataframe[34:47,]

AM_A_O_D_versus_E_R_F_only_gp_2_gbs_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_A_O_D_versus_E_R_F_only_gp_2_gbs_matrix[long, lat] <- 
      if (R2_adj_A_O_D_matrix[long, lat] < 0.6) AM_A_O_D_versus_E_R_F_only_gp_matrix[long, lat]
    else NA
  }
  cat(long, "/192, ", sep = "")
}

image.plot(longitude, latitude,
           rbind(AM_A_O_D_versus_E_R_F_only_gp_matrix[97:192,],
                 AM_A_O_D_versus_E_R_F_only_gp_matrix[1:96,]),
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

# Comparison between AM's obtaining with 100 versus 100,000

assign(paste0("AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N),
       matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", N, ".txt"))),
              nrow = length(longitude))
)
assign(paste0("AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", 100000),
       matrix(as.numeric(readLines(paste0("object/AM_A_O_D_versus_E_R_F_matrix_0", threshold*100, "_", 100000, ".txt"))),
              nrow = length(longitude))
)
plot(as.vector(AM_A_O_D_versus_E_R_F_matrix_060_100000)[which(as.vector(R2_adj_A_O_D_matrix) >= 0.60)] ~
       as.vector(AM_A_O_D_versus_E_R_F_matrix_060_100)[which(as.vector(R2_adj_A_O_D_matrix) >= 0.60)],
     col = "black", pch = 15, xlim = c(0,1), ylim = c(0,1),
     xlab = "AM for Total AOD (level 2) and Total ERF (both July), N = 100, threshold = 0.60",
     ylab = "AM for Total AOD (level 2) and Total ERF (both July), N = 100,000, threshold = 0.60")
points(as.vector(AM_A_O_D_versus_E_R_F_matrix_060_100000)[which(as.vector(R2_adj_A_O_D_matrix) < 0.60)] ~
         as.vector(AM_A_O_D_versus_E_R_F_matrix_060_100)[which(as.vector(R2_adj_A_O_D_matrix) < 0.60)],
       col = "blue", pch = 15)
range(AM_A_O_D_versus_E_R_F_matrix_060_100); range(AM_A_O_D_versus_E_R_F_matrix_060_100000)
sum(AM_A_O_D_versus_E_R_F_matrix_060_100 - AM_A_O_D_versus_E_R_F_matrix_060_100000 > 0.05)


