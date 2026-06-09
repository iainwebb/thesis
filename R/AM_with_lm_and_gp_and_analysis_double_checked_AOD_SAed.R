# Preamble #####################################################################

## This file aims to test the alignment measure with two variables,one global 
## and by-gridbox. They are:
## ACURE_P3_AOD_Total_jul_level2, referred to throughout as "AOD"
## and
## ERF_PPE_global_jul, referred to as "ERF"
## These are taken from the A-CURE monthly mean PPE from the first version of 
## the U.K. Earth System Model (UKESM1)

## Various objects are exported to /objects to avoid having to create them from
## scratch each time. To help with this, we specify the following:
options(scipen=999)

# Packages #####################################################################

## We load the following packages. If running locally, use
# library(readr, lib="C:/Users/smp22ijw/Desktop/Library/")
# library(DiceKriging, lib="C:/Users/smp22ijw/Desktop/Library/")
## or otherwise use
library(readr)
library(DiceKriging)

# Importing the inputs  ########################################################

## The normalised input settings from the PPE are read in as a 221 × 37 
## (transposed) matrix
inputs_x_norm_T_matrix_SAed <- 
  as.matrix(
    read.csv("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")
  )[,ERF_SAed_indices]
inputs_x_norm_matrix_SAed <- t(inputs_x_norm_T_matrix_SAed)

## There are p variables in all, and n runs
p_SAed <- ncol(inputs_x_norm_T_matrix_SAed)
n <- nrow(inputs_x_norm_T_matrix_SAed)

# Importing the outputs ########################################################

## AOD is the observable variable. Since it is by-gridbox, it's saved as a
## NetCDF file, so we first load the required package, with the local code first
# library(ncdf4, lib="C:/Users/smp22ijw/Desktop/Library/")
library(ncdf4)
response_AOD_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc")
)
response_AOD_all_levels_array <- ncvar_get(response_AOD_nc, names(response_AOD_nc$var))

## We only want level 2, so extract that as a 192 × 144 × 221 array
response_AOD_array <- response_AOD_all_levels_array[,,2,]

## Remove all unneeded objects
rm("response_AOD_nc", "response_AOD_all_levels_array")

## ERF is the unobserved variable, imported as 221 × 1 dataframe, with the 
## local version first
# response_ERF_dataframe <- read.table("data/ERF_PPE_global_jul.dat", col.names = "ERF")
response_ERF_dataframe <- read.table(
  "X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/global_mean_forcing/ERF_PPE_global_jul.dat", 
  col.names = "ERF"
)

# Importing map details ########################################################

## We import the same NetCDF as previously, then extract longitudes and 
## latitudes, adjusting the longitudes ready for producing maps
response_AOD_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc")
)
longitude = ncvar_get(response_AOD_nc,"longitude") - 180
latitude = ncvar_get(response_AOD_nc,"latitude")

## We also load required packages, with local code first
# library(readr, lib="C:/Users/smp22ijw/Desktop/Library/")
library(fields)
# library(maps, lib="C:/Users/smp22ijw/Desktop/Library/")
library(maps)

## Remove all unneeded objects
rm("response_AOD_nc")

# lm work ######################################################################

## In case linear models are fine, and there's no need to fit GPs, we check 
# adjusted R^2 values

## For ERF this only needs to be done once
R2adj_ERF_SAed <- summary(
  lm(ERF ~., data.frame(cbind("ERF" = response_ERF_dataframe, inputs_x_norm_T_matrix_SAed)))
  )$adj.r.squared

## Next for AOD, by gridbox
R2adj_AOD_matrix_SAed <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    df <- data.frame(
      cbind("AOD" = response_AOD_array[long, lat, ], inputs_x_norm_T_matrix_SAed)
      )
    R2adj_AOD_matrix_SAed[long, lat] <- summary(lm(AOD ~., df))$adj.r.squared
  }
  cat(long, "/192 ", sep = "")
}
rm("long", "lat", "df")

## We export both objects to /objects, with the import code below
write_lines(as.vector(R2adj_AOD_matrix_SAed),
            file="objects/R2adj_AOD_matrix_SAed.txt")
# R2adj_AOD_matrix_SAed <-
#   matrix(as.numeric(readLines("objects/R2adj_AOD_matrix_SAed.txt")),
#          nrow = length(longitude))
write_lines(R2adj_ERF_SAed,
            file="objects/R2adj_ERF_SAed.txt")
# R2adj_ERF_SAed <-
#   as.numeric(readLines("objects/R2adj_ERF_SAed.txt"))

## ERF's adjusted R^2 is low, so we will fit a GP. For AOD, we can produce a
## table and plots to understand how many and which gridboxes require a GP

table(
  ifelse(floor(R2adj_AOD_matrix_SAed * 10) / 10 < 0.1, "<0.1",
         floor(R2adj_AOD_matrix_SAed * 10) / 10)
)

100*round(
  table(
    ifelse(floor(R2adj_AOD_matrix_SAed * 10) / 10 < 0.5, "<0.5",
           floor(R2adj_AOD_matrix_SAed * 10) / 10)
  ) / (192*144), 3)

hist(R2adj_AOD_matrix_SAed, xlim = c(0,1))

image.plot(longitude, latitude,
           rbind(R2adj_AOD_matrix_SAed[97:192,],
                 R2adj_AOD_matrix_SAed[1:96,]),
           breaks = c(0,0.1,0.6,0.7,0.75,0.8,1),
           col = c("magenta","#0000ff","orange","green","yellow", "white"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.6,0.7,0.75,0.8,1),
             labels=as.character(c(0,0.1,0.6,0.7,0.75,0.8,1)),
             mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1)

## Let's fit an lm for every gridbox, for use with most later on, and also for 
## comparison when fitting GPs. We'll also set a threshold for the adjusted R^2
threshold <- 0.10

pdsnormed_AOD_lmall_array_SAed <- 
  pdsnormed_AOD_lmbestNArest_array_SAed <- 
  array(NA, dim = c(length(longitude), length(latitude), p_SAed))
for (long in 1:192) {
  for (lat in 1:144) {
    df <- data.frame(
      cbind("AOD" = response_AOD_array[long, lat, ], inputs_x_norm_T_matrix_SAed)
    )
    lin_mod <- lm(AOD ~., df)
    R2adj <- summary(lin_mod)$adj.r.squared
    pdslm <- as.vector(summary(lin_mod)$coefficients[2:10, 1])
    pdsnormed_AOD_lmall_array_SAed[long, lat,] <- pdslm / sqrt(sum(pdslm^2))
    pdsnormed_AOD_lmbestNArest_array_SAed[long, lat,] <-
      if (R2adj >= threshold) pdsnormed_AOD_lmall_array_SAed[long, lat,]
    else rep(NA, p_SAed)
  }
  cat(long, "/192 ", sep = "")
}
rm("long", "lat", "df", "lin_mod", "R2adj", "pdslm")

# We'll export both objects, with import code below
write_lines(as.vector(pdsnormed_AOD_lmall_array_SAed),
            file="objects/pdsnormed_AOD_lmall_array_SAed.txt")
# pdsnormed_AOD_lmall_array_SAed <-
#   array(as.numeric(readLines("objects/pdsnormed_AOD_lmall_array_SAed.txt")),
#         dim = c(length(longitude), length(latitude), p))
write_lines(as.vector(pdsnormed_AOD_lmbestNArest_array_SAed),
            file=paste0("objects/pdsnormed_AOD_lmbestNArest_array_0", threshold*100, "_SAed.txt"))
# assign(paste0("pdsnormed_AOD_lmbestNArest_array_SAed"),
#        array(as.numeric(readLines(paste0("objects/pdsnormed_AOD_lmbestNArest_array_0", threshold*100, "_SAed.txt"))),
#              dim = c(length(longitude), length(latitude), p))
# )
# sum(is.na(pdsnormed_AOD_lmbestNArest_array)) / 37 # check

# GP set-up ####################################################################

## We specify some variables, namely N (the size of the sample from the input
## parameter space) and q (the number of betas that will be estimated) 

N <- 100000
q_SAed <- 1 + p_SAed

## We load the required packages for sampling
library(lhs)

## We create the sample
xstar_matrix_SAed <- t(randomLHS(N, p_SAed))
# TO DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# xstar_T_dataframe <- data.frame(t(xstar_matrix))
# colnames(xstar_T_dataframe) <- colnames(inputs_x_norm_T_matrix)
# TO DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Our estimated posterior derivative requires H, which doesn't depend on the 
# response (and thus on the gridbox in the case of AOD)
H_matrix_SAed <- cbind(c(rep(1, n)), t(unname(inputs_x_norm_matrix_SAed)))

## We define a function called gp_d_dxi_normalised_function which returns the
## normalised GP estimates of the partial derivatives. It takes as its single
## argument a 221-long vector or 221 x 1 dataframe of responses

gp_d_dxi_normalised_function <- function(response, loglik_index) {
  
  gp <- km(~., 
           design = inputs_x_norm_T_matrix_SAed, 
           response = response, 
           covtype="gauss", optim.method="BFGS", control=list(maxit=500)
           )
  
  ## Using the posterior estimates of the range parameters, we make the matrix 
  ## A, its inverse, and t(xstar)^T
  lhat_vector_SAed <- gp@covariance@range.val
  lhat_matrix_SAed <- matrix(lhat_vector_SAed, nrow = p_SAed)
  lhatinv_diag_matrix_SAed <- diag(1/lhat_vector_SAed)
  
  A_matrix_SAed <- matrix(NA, nrow = n, ncol = n)
  for (r in 1:n) {
    for (c in 1:n) {
      if (is.na(A_matrix_SAed[r,c])) {
        A_matrix_SAed[r,c] <- A_matrix_SAed[c,r] <- 
          exp(-0.5 *
                matrix(inputs_x_norm_matrix_SAed[,r] - inputs_x_norm_matrix_SAed[,c], nrow = 1) %*%
                (lhatinv_diag_matrix_SAed^2) %*%
                matrix(inputs_x_norm_matrix_SAed[,r] - inputs_x_norm_matrix_SAed[,c], ncol = 1)
              )
      }
    }
  }
  Ainv_matrix_SAed <- solve(A_matrix_SAed)
  
  txstar_T_matrix_SAed <- matrix(NA, nrow = N, ncol = n)
  for (r in 1:N) {
    for (c in 1:n) {
      txstar_T_matrix_SAed[r,c] <- exp(-0.5 *
                                      matrix(xstar_matrix_SAed[,r] - inputs_x_norm_matrix_SAed[,c], nrow = 1) %*% 
                                      (lhatinv_diag_matrix_SAed^2) %*%
                                      matrix(xstar_matrix_SAed[,r] - inputs_x_norm_matrix_SAed[,c], ncol = 1)
      )
    }
    cat("txstarT: ", r, "/", N, ",", sep = "")
  }
  
  ## Secondly, we use the trend estimates to estimate the 37 posterior partial
  ## derivatives at each of the N points, then normalise them
  betahat_matrix_SAed <- matrix(gp@trend.coef, nrow = q_SAed)
  
  pd_T_dataframe_SAed <- data.frame(matrix(NA, nrow = p_SAed, ncol = N))
  
  for (r in 1:p_SAed) {
    pd_T_dataframe_SAed[r,] <-
      matrix(betahat_matrix_SAed[r+1,], nrow = N) + 
      (((matrix(rep(inputs_x_norm_matrix_SAed[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix_SAed[r,], n), nrow = N, byrow = F)) / (lhat_vector_SAed[r]^2)) *
         txstar_T_matrix_SAed) %*% 
      Ainv_matrix_SAed %*% (as.matrix(response, ncol = n) - H_matrix_SAed %*% betahat_matrix_SAed)
    cat("pd_t:", r, "/", p_SAed, ",", sep = "")
  }
  
  pd_dataframe_SAed <- t(pd_T_dataframe_SAed)
  
  pdnormed_dataframe_SAed <- pd_dataframe_SAed / sqrt(rowSums(pd_dataframe_SAed^2))
  
  pdnormed_T_matrix_SAed <- t(data.matrix(pdnormed_dataframe_SAed))
  
  ## just when repeating runs
  # assign(paste0("loglik_long142lat88_", loglik_index), logLik.km(gp), 
  #        envir = .GlobalEnv)
  
  return(pdnormed_T_matrix_SAed)
}

## First we'll fit the GP for ERF and calculate the normalised partial
## derivatives
gppdnormalised_t_ERF_matrix_SAed <- gp_d_dxi_normalised_function(response_ERF_dataframe)

## This produces the following AMs (with the lm-estimated partial derivatives
## for AOD):
AM_ERFgp_AODlmall_matrix_SAed <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_AODlmall_matrix_SAed[long,lat] <-
      mean(abs(colSums(gppdnormalised_t_ERF_matrix_SAed * matrix(rep(pdsnormed_AOD_lmall_array_SAed[long,lat,], N), ncol = N))))
    # AM_A_O_D_versus_E_R_F_matrix_case1[long, lat] <-
    #   if (R2_adj_A_O_D_matrix[long, lat] >= threshold & R2_adj_E_R_F_matrix[long, lat] >= threshold)
    #     abs(sum(lm_best_d_dxi_normalised_A_O_D_array[long, lat,] *
    #                       lm_d_dxi_normalised_E_R_F_array[long, lat,]))
    #   else NA
    # if (R2_adj_A_O_D_matrix[long, lat] >= threshold & R2_adj_E_R_F_matrix[long, lat] >= threshold)
    #   write_lines(as.vector(AM_A_O_D_versus_E_R_F_matrix_case1),
    #               file=paste0("objects/AM_A_O_D_versus_E_R_F_matrix_case1_0", threshold*100, "_", N, ".txt"))
  }
  cat(long, "/192 ", sep = "")
}

AM_ERFgp_AODlmbest_matrix_SAed <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_AODlmbest_matrix_SAed[long, lat] <-
      if (R2adj_AOD_matrix_SAed[long, lat] < threshold)
        -0.05
    else AM_ERFgp_AODlmall_matrix_SAed[long, lat]
  }
}
sum(is.na(AM_ERFgp_AODlmbest_matrix_SAed))

## We produce a map to show these AMs, excluding the one for which the adjusted
## R^2 was too low
image.plot(longitude, latitude,
           rbind(AM_ERFgp_AODlmbest_matrix_SAed[97:192,],
                 AM_ERFgp_AODlmbest_matrix_SAed[1:96,]),
           breaks = c(-0.1,0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           col = c("blue","red4","red3","red1","indianred1", "mistyrose",
                   "honeydew1", "darkseagreen1", "palegreen2", "seagreen3", 
                   "springgreen4"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(-0.06,0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c("NA",0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

## as well as one with them all in
image.plot(longitude, latitude,
           rbind(AM_ERFgp_AODlmall_matrix_SAed[97:192,],
                 AM_ERFgp_AODlmall_matrix_SAed[1:96,]),
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

## Next we'll fit GPs for AOD for those gridboxes where the adjusted R^2 was
## deemed too low. It gets exported each new gridbox that gets added
AM_ERFgp_AODlmbestgprest_matrix_SAed <- matrix(NA, nrow = 192, ncol = 144)

index_counter <- 0

for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_AODlmbestgprest_matrix_SAed[long, lat] <-
      if (R2adj_AOD_matrix_SAed[long, lat] < threshold)
      mean(abs(colSums(gppdnormalised_t_ERF_matrix_SAed * gp_d_dxi_normalised_function(response_AOD_array[long, lat,]))))
    else AM_ERFgp_AODlmall_matrix_SAed[long, lat]
    index_counter <-
      if (R2adj_AOD_matrix_SAed[long, lat] < threshold)
        index_counter + 1
    else index_counter
    if (R2adj_AOD_matrix_SAed[long, lat] < threshold)
      write_lines(as.vector(AM_ERFgp_AODlmbestgprest_matrix_SAed),
                  file=paste0("objects/AM_ERFgp_AODlmbestgprest_matrix_0", threshold*100, "_", N, "_SAed.txt"))
  }
  write_lines(as.vector(AM_ERFgp_AODlmbestgprest_matrix_SAed),
              file=paste0("objects/AM_ERFgp_AODlmbestgprest_matrix_0", threshold*100, "_", N, "_SAed.txt")) 
}
# rm("long", "lat", "p", "q", "x_star_matrix", "inputs_x_norm_matrix", "x_star_T_dataframe",
#    "n", "H_matrix", "corGaussian", "gp_d_dxi_normalised_function", "index_counter")

sum(is.na(AM_ERFgp_AODlmbestgprest_matrix_SAed))

range(AM_ERFgp_AODlmbestgprest_matrix_SAed)

table(floor(AM_ERFgp_AODlmall_matrix_SAed * 10) / 10)
table(floor(AM_ERFgp_AODlmbest_matrix_SAed * 10) / 10)
table(floor(AM_ERFgp_AODlmbestgprest_matrix_SAed * 10) / 10)

100*round(
  table(
    ifelse(floor(R2adj_AOD_matrix_SAed * 10) / 10 < 0.6, "<0.6",
           floor(R2adj_AOD_matrix_SAed * 10) / 10)
  ) / (192*144), 3)

hist(R2adj_AOD_matrix_SAed, xlim = c(0,1))

image.plot(longitude, latitude,
           rbind(AM_ERFgp_AODlmbestgprest_matrix_SAed[97:192,],
                 AM_ERFgp_AODlmbestgprest_matrix_SAed[1:96,]),
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

# R2_adj_A_O_D_matrix <-
#   matrix(as.numeric(readLines("objects/R2_adj_A_O_D_matrix.txt")),
#          nrow = length(longitude))
R2adj_AMalllm_SAed <- cbind("AOD_R2adj" 
                          = as.vector(R2adj_AOD_matrix_SAed), 
                          "AM_AOD_versus_ERFalllm_matrix" 
                          = as.vector(AM_ERFgp_AODlmall_matrix_SAed)
)
blue_black_colours_SAed <- ifelse(as.vector(R2adj_AOD_matrix_SAed) < 0.10, "blue", "black")
plot(R2adj_AMalllm_SAed, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours_SAed,
     xlab = "R2_adj for AOD_Total",
     ylab = "AM AOD_Total versus ERF")
R2adj_AM_SAed <- cbind("AOD_Total_R2_adj" 
                   = as.vector(R2adj_AOD_matrix_SAed), 
                   "AM_for_AOD_Total_versus_ERF" 
                   = as.vector(AM_ERFgp_AODlmbestgprest_matrix_SAed)
)
plot(R2adj_AM_SAed, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours_SAed,
     xlab = "R2_adj for AOD_Total",
     ylab = "AM AOD_Total versus ERF")
sum(abs(AM_ERFgp_AODlmall_matrix_SAed[R2adj_AOD_matrix_SAed < 0.10] - 
  AM_ERFgp_AODlmbestgprest_matrix_SAed[R2adj_AOD_matrix_SAed < 0.10]) < 0.00001)

############################
gped <- 1
which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)
AM_gped <- rep(NA, 113)
for (gped in 1:113) {
  AM_gped[gped] <- 
    AM_ERFgp_AODlmbestgprest_matrix[which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[gped,1],
                                    which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[gped,2]]
}
table(floor(AM_gped * 10) / 10)
which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[which(AM_gped == max(AM_gped)),]
AM_ERFgp_AODlmbestgprest_matrix[16,51]
AM_ERFgp_AODlmall_matrix[16,51]

gped_lmAMtoAMgpdiff <- rep(NA, 113)
for (gped in 1:113) {
  gped_lmAMtoAMgpdiff[gped] <- 
    AM_ERFgp_AODlmall_matrix[which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[gped,1],
                             which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[gped,2]] -
    AM_ERFgp_AODlmbestgprest_matrix[which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[gped,1],
                                    which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[gped,2]]
}
which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[which(gped_lmAMtoAMgpdiff == max(gped_lmAMtoAMgpdiff))]
AM_ERFgp_AODlmall_matrix[which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[59,1],
                         which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[59,2]]
AM_ERFgp_AODlmbestgprest_matrix[which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[59,1],
                                which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[59,2]]
which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[59,]

long <- 142; lat <- 88
R2adj_AOD_matrix[long, lat] < threshold
AM_ERFgp_AODlmall_matrix[long, lat]
AM_ERFgp_AODlmbestgprest_matrix[long, lat]
# AM_long142lat88 <- rep(NA, 100)
# loglik_long142lat88 <- rep(NA, 100)
for (reps in 2:100) {
  AM_long142lat88[reps] <-
    mean(abs(colSums(gppdnormalised_t_ERF_matrix * gp_d_dxi_normalised_function(response_AOD_array[long, lat,], reps))))
  loglik_long142lat88[reps] <- eval(parse(text = paste0("loglik_long142lat88_", reps)))
}
AM_long142lat88
loglik_long142lat88
stripchart(AM_long142lat88, method = "jitter", pch = 20, col = "blue", main = "One-Dimensional Scatterplot")
abline(v = AM_ERFgp_AODlmall_matrix[which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[59,1],
                                    which(R2adj_AOD_matrix < threshold, arr.ind = TRUE)[59,2]])
