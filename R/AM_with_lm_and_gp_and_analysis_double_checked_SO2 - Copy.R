# Preamble #####################################################################

## This file aims to test the alignment measure with two variables,one global 
## and by-gridbox. They are:
## ACURE_P3_SO2_Total_jul_level2, referred to throughout as "SO2"
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
inputs_x_norm_T_matrix <- 
  as.matrix(
    read.csv("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")
  )
inputs_x_norm_matrix <- t(inputs_x_norm_T_matrix)

## There are p variables in all, and n runs
p <- ncol(inputs_x_norm_T_matrix)
n <- nrow(inputs_x_norm_T_matrix)

# Importing the outputs ########################################################

## SO2 is the observable variable. Since it is by-gridbox, it's saved as a
## NetCDF file, so we first load the required package, with the local code first
# library(ncdf4, lib="C:/Users/smp22ijw/Desktop/Library/")
library(ncdf4)
response_SO2_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/SO2/ACURE_P3_SO2_PD_221_jul.nc")
)
response_SO2_array <- ncvar_get(response_SO2_nc, names(response_SO2_nc$var))

## Remove all unneeded objects
rm("response_SO2_nc")

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
response_SO2_nc <- nc_open(
  paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/SO2/ACURE_P3_SO2_PD_221_jul.nc")
)
longitude = ncvar_get(response_SO2_nc,"longitude") - 180
latitude = ncvar_get(response_SO2_nc,"latitude")

## We also load required packages, with local code first
# library(readr, lib="C:/Users/smp22ijw/Desktop/Library/")
library(fields)
# library(maps, lib="C:/Users/smp22ijw/Desktop/Library/")
library(maps)

## Remove all unneeded objects
rm("response_SO2_nc")

# lm work ######################################################################

## In case linear models are fine, and there's no need to fit GPs, we check 
# adjusted R^2 values

## For ERF this only needs to be done once
R2adj_ERF <- summary(
  lm(ERF ~., data.frame(cbind("ERF" = response_ERF_dataframe, inputs_x_norm_T_matrix)))
)$adj.r.squared

## Next for SO2, by gridbox
R2adj_SO2_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    df <- data.frame(
      cbind("SO2" = response_SO2_array[long, lat, ], inputs_x_norm_T_matrix)
    )
    R2adj_SO2_matrix[long, lat] <- summary(lm(SO2 ~., df))$adj.r.squared
  }
  cat(long, "/192 ", sep = "")
}
rm("long", "lat", "df")

## We export both objects to /objects, with the import code below
write_lines(as.vector(R2adj_SO2_matrix),
            file="objects/R2adj_SO2_matrix.txt")
# R2adj_SO2_matrix <-
#   matrix(as.numeric(readLines("objects/R2adj_SO2_matrix.txt")),
#          nrow = length(longitude))
write_lines(R2adj_ERF,
            file="objects/R2adj_ERF.txt")
# R2adj_ERF <-
#   as.numeric(readLines("objects/R2adj_ERF.txt"))

## ERF's adjusted R^2 is low, so we will fit a GP. For SO2, we can produce a
## table and plots to understand how many and which gridboxes require a GP

table(
  ifelse(floor(R2adj_SO2_matrix * 10) / 10 < 0.5, "<0.5",
         floor(R2adj_SO2_matrix * 10) / 10)
)

100*round(
  table(
    ifelse(floor(R2adj_SO2_matrix * 10) / 10 < 0.5, "<0.5",
           floor(R2adj_SO2_matrix * 10) / 10)
  ) / (192*144), 3)

hist(R2adj_SO2_matrix, xlim = c(0,1))

image.plot(longitude, latitude,
           rbind(R2adj_SO2_matrix[97:192,],
                 R2adj_SO2_matrix[1:96,]),
           breaks = c(0,0.5,0.6,0.7,0.75,0.8,1),
           col = c("magenta","#0000ff","orange","green","yellow", "white"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.5,0.6,0.7,0.75,0.8,1),
             labels=as.character(c(0,0.5,0.6,0.7,0.75,0.8,1)),
             mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1)

## Let's fit an lm for every gridbox, for use with most later on, and also for 
## comparison when fitting GPs. We'll also set a threshold for the adjusted R^2
threshold <- 0.50

pdsnormed_SO2_lmall_array <- 
  pdsnormed_SO2_lmbestNArest_array <- 
  array(NA, dim = c(length(longitude), length(latitude), p))
for (long in 1:192) {
  for (lat in 1:144) {
    df <- data.frame(
      cbind("SO2" = response_SO2_array[long, lat, ], inputs_x_norm_T_matrix)
    )
    lin_mod <- lm(SO2 ~., df)
    R2adj <- summary(lin_mod)$adj.r.squared
    pdslm <- as.vector(summary(lin_mod)$coefficients[2:38, 1])
    pdsnormed_SO2_lmall_array[long, lat,] <- pdslm / sqrt(sum(pdslm^2))
    pdsnormed_SO2_lmbestNArest_array[long, lat,] <-
      if (R2adj >= threshold) pdsnormed_SO2_lmall_array[long, lat,]
    else rep(NA, p)
  }
  cat(long, "/192 ", sep = "")
}
rm("long", "lat", "df", "lin_mod", "R2adj", "pdslm")

# We'll export both objects, with import code below
write_lines(as.vector(pdsnormed_SO2_lmall_array),
            file="objects/pdsnormed_SO2_lmall_array.txt")
# pdsnormed_SO2_lmall_array <-
#   array(as.numeric(readLines("objects/pdsnormed_SO2_lmall_array.txt")),
#         dim = c(length(longitude), length(latitude), p))
write_lines(as.vector(pdsnormed_SO2_lmbestNArest_array),
            file=paste0("objects/pdsnormed_SO2_lmbestNArest_array_0", threshold*100, ".txt"))
# assign(paste0("pdsnormed_SO2_lmbestNArest_array"),
#        array(as.numeric(readLines(paste0("objects/pdsnormed_SO2_lmbestNArest_array_0", threshold*100, ".txt"))),
#              dim = c(length(longitude), length(latitude), p))
# )
# sum(is.na(pdsnormed_SO2_lmbestNArest_array)) / 37 # check

# GP set-up ####################################################################

## We specify some variables, namely N (the size of the sample from the input
## parameter space) and q (the number of betas that will be estimated) 

N <- 100000
q <- 1 + p

## We load the required packages for sampling
library(lhs)

## We create the sample
xstar_matrix <- t(randomLHS(N, p))
## Export this to /objects, with the import code below
# write_lines(as.vector(xstar_matrix),
#             file="objects/xstar_matrix.txt")
xstar_matrix <-
  matrix(as.numeric(readLines("objects/xstar_matrix.txt")),
         nrow = length(longitude))
# TO DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# xstar_T_dataframe <- data.frame(t(xstar_matrix))
# colnames(xstar_T_dataframe) <- colnames(inputs_x_norm_T_matrix)
# TO DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Our estimated posterior derivative requires H, which doesn't depend on the 
# response (and thus on the gridbox in the case of SO2)
H_matrix <- cbind(c(rep(1, n)), t(unname(inputs_x_norm_matrix)))

## We define a function called gp_d_dxi_normalised_function which returns the
## normalised GP estimates of the partial derivatives. It takes as its single
## argument a 221-long vector or 221 x 1 dataframe of responses

response <- response_ERF_dataframe

gp_d_dxi_normalised_function <- function(response, loglik_index) {
  
  gp <- km(~., 
           design = inputs_x_norm_T_matrix, 
           response = response, 
           covtype="gauss", optim.method="BFGS", control=list(maxit=500)
  )
  
  ## Using the posterior estimates of the range parameters, we make the matrix 
  ## A, its inverse, and t(xstar)^T
  lhat_vector <- gp@covariance@range.val
  lhat_matrix <- matrix(lhat_vector, nrow = p)
  lhatinv_diag_matrix <- diag(1/lhat_vector)
  
  A_matrix <- matrix(NA, nrow = n, ncol = n)
  for (r in 1:n) {
    for (c in 1:n) {
      if (is.na(A_matrix[r,c])) {
        A_matrix[r,c] <- A_matrix[c,r] <- 
          exp(-0.5 *
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*%
                (lhatinv_diag_matrix^2) %*%
                matrix(inputs_x_norm_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
          )
      }
    }
  }
  Ainv_matrix <- solve(A_matrix)
  
  txstar_T_matrix <- matrix(NA, nrow = N, ncol = n)
  for (r in 1:N) {
    for (c in 1:n) {
      txstar_T_matrix[r,c] <- exp(-0.5 *
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], nrow = 1) %*% 
                                    (lhatinv_diag_matrix^2) %*%
                                    matrix(xstar_matrix[,r] - inputs_x_norm_matrix[,c], ncol = 1)
      )
    }
    cat("txstarT: ", r, "/", N, ",", sep = "")
  }
  
  ## Secondly, we use the trend estimates to estimate the 37 posterior partial
  ## derivatives at each of the N points, then normalise them
  betahat_matrix <- matrix(gp@trend.coef, nrow = q)
  
  pd_T_dataframe <- data.frame(matrix(NA, nrow = p, ncol = N))
  
  for (r in 1:p) {
    pd_T_dataframe[r,] <-
      matrix(betahat_matrix[r+1,], nrow = N) + 
      (((matrix(rep(inputs_x_norm_matrix[r,], N), nrow = N, byrow = T) - matrix(rep(xstar_matrix[r,], n), nrow = N, byrow = F)) / (lhat_vector[r]^2)) *
         txstar_T_matrix) %*% 
      Ainv_matrix %*% (as.matrix(response, ncol = n) - H_matrix %*% betahat_matrix)
    cat("pd_t:", r, "/", p, ",", sep = "")
  }
  
  pd_dataframe <- t(pd_T_dataframe)
  
  pdnormed_dataframe <- pd_dataframe / sqrt(rowSums(pd_dataframe^2))
  
  pdnormed_T_matrix <- t(data.matrix(pdnormed_dataframe))
  
  ## just when repeating runs
  # assign(paste0("loglik_long142lat88_", loglik_index), logLik.km(gp), 
  #        envir = .GlobalEnv)
  
  return(pdnormed_T_matrix)
}

## First we'll fit the GP for ERF and calculate the normalised partial
## derivatives
gppdnormalised_t_ERF_matrix <- gp_d_dxi_normalised_function(response_ERF_dataframe)

## This produces the following AMs (with the lm-estimated partial derivatives
## for SO2):
AM_ERFgp_SO2lmall_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_SO2lmall_matrix[long,lat] <-
      mean(abs(colSums(gppdnormalised_t_ERF_matrix * matrix(rep(pdsnormed_SO2_lmall_array[long,lat,], N), ncol = N))))
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

AM_ERFgp_SO2lmbest_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_SO2lmbest_matrix[long, lat] <-
      if (R2adj_SO2_matrix[long, lat] < threshold)
        -0.05
    else AM_ERFgp_SO2lmall_matrix[long, lat]
  }
}
sum(is.na(AM_ERFgp_SO2lmbest_matrix))

## We produce a map to show these AMs, excluding the one for which the adjusted
## R^2 was too low
image.plot(longitude, latitude,
           rbind(AM_ERFgp_SO2lmbest_matrix[97:192,],
                 AM_ERFgp_SO2lmbest_matrix[1:96,]),
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
           rbind(AM_ERFgp_SO2lmall_matrix[97:192,],
                 AM_ERFgp_SO2lmall_matrix[1:96,]),
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

## Next we'll fit GPs for SO2 for those gridboxes where the adjusted R^2 was
## deemed too low. It gets exported each new gridbox that gets added
AM_ERFgp_SO2lmbestgprest_matrix <- matrix(NA, nrow = 192, ncol = 144)

index_counter <- 0

for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_SO2lmbestgprest_matrix[long, lat] <-
      if (R2adj_SO2_matrix[long, lat] < threshold)
        mean(abs(colSums(gppdnormalised_t_ERF_matrix * gp_d_dxi_normalised_function(response_SO2_array[long, lat,]))))
    else AM_ERFgp_SO2lmall_matrix[long, lat]
    index_counter <-
      if (R2adj_SO2_matrix[long, lat] < threshold)
        index_counter + 1
    else index_counter
    if (R2adj_SO2_matrix[long, lat] < threshold)
      write_lines(as.vector(AM_ERFgp_SO2lmbestgprest_matrix),
                  file=paste0("objects/AM_ERFgp_SO2lmbestgprest_matrix_0", threshold*100, "_", N, ".txt"))
  }
  write_lines(as.vector(AM_ERFgp_SO2lmbestgprest_matrix),
              file=paste0("objects/AM_ERFgp_SO2lmbestgprest_matrix_0", threshold*100, "_", N, ".txt")) 
}
# rm("long", "lat", "p", "q", "x_star_matrix", "inputs_x_norm_matrix", "x_star_T_dataframe",
#    "n", "H_matrix", "corGaussian", "gp_d_dxi_normalised_function", "index_counter")

sum(is.na(AM_ERFgp_SO2lmbestgprest_matrix))

range(AM_ERFgp_SO2lmbestgprest_matrix)

table(floor(AM_ERFgp_SO2lmall_matrix * 10) / 10)
table(floor(AM_ERFgp_SO2lmbest_matrix * 10) / 10)
table(floor(AM_ERFgp_SO2lmbestgprest_matrix * 10) / 10)

100*round(
  table(
    ifelse(floor(R2adj_SO2_matrix * 10) / 10 < 0.6, "<0.6",
           floor(R2adj_SO2_matrix * 10) / 10)
  ) / (192*144), 3)

hist(R2adj_SO2_matrix, xlim = c(0,1))

image.plot(longitude, latitude,
           rbind(AM_ERFgp_SO2lmbestgprest_matrix[97:192,],
                 AM_ERFgp_SO2lmbestgprest_matrix[1:96,]),
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
R2adj_AMalllm <- cbind("SO2_R2adj" 
                       = as.vector(R2adj_SO2_matrix), 
                       "AM_SO2_versus_ERFalllm_matrix" 
                       = as.vector(AM_ERFgp_SO2lmall_matrix)
)
blue_black_colours <- ifelse(as.vector(R2adj_SO2_matrix) < 0.60, "blue", "black")
plot(R2adj_AMalllm, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours,
     xlab = "R2_adj for SO2_Total",
     ylab = "AM SO2_Total versus ERF")
R2adj_AM <- cbind("SO2_Total_R2_adj" 
                  = as.vector(R2adj_SO2_matrix), 
                  "AM_for_SO2_Total_versus_ERF" 
                  = as.vector(AM_ERFgp_SO2lmbestgprest_matrix)
)
plot(R2adj_AM, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours,
     xlab = "R2_adj for SO2_Total",
     ylab = "AM SO2_Total versus ERF")


############################
gped <- 1
which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)
AM_gped <- rep(NA, 113)
for (gped in 1:113) {
  AM_gped[gped] <- 
    AM_ERFgp_SO2lmbestgprest_matrix[which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[gped,1],
                                    which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[gped,2]]
}
table(floor(AM_gped * 10) / 10)
which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[which(AM_gped == max(AM_gped)),]
AM_ERFgp_SO2lmbestgprest_matrix[16,51]
AM_ERFgp_SO2lmall_matrix[16,51]

gped_lmAMtoAMgpdiff <- rep(NA, 113)
for (gped in 1:113) {
  gped_lmAMtoAMgpdiff[gped] <- 
    AM_ERFgp_SO2lmall_matrix[which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[gped,1],
                             which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[gped,2]] -
    AM_ERFgp_SO2lmbestgprest_matrix[which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[gped,1],
                                    which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[gped,2]]
}
which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[which(gped_lmAMtoAMgpdiff == max(gped_lmAMtoAMgpdiff))]
AM_ERFgp_SO2lmall_matrix[which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[59,1],
                         which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[59,2]]
AM_ERFgp_SO2lmbestgprest_matrix[which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[59,1],
                                which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[59,2]]
which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[59,]

long <- 142; lat <- 88
R2adj_SO2_matrix[long, lat] < threshold
AM_ERFgp_SO2lmall_matrix[long, lat]
AM_ERFgp_SO2lmbestgprest_matrix[long, lat]
# AM_long142lat88 <- rep(NA, 100)
# loglik_long142lat88 <- rep(NA, 100)
for (reps in 2:100) {
  AM_long142lat88[reps] <-
    mean(abs(colSums(gppdnormalised_t_ERF_matrix * gp_d_dxi_normalised_function(response_SO2_array[long, lat,], reps))))
  loglik_long142lat88[reps] <- eval(parse(text = paste0("loglik_long142lat88_", reps)))
}
AM_long142lat88
loglik_long142lat88
stripchart(AM_long142lat88, method = "jitter", pch = 20, col = "blue", main = "One-Dimensional Scatterplot")
abline(v = AM_ERFgp_SO2lmall_matrix[which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[59,1],
                                    which(R2adj_SO2_matrix < threshold, arr.ind = TRUE)[59,2]])
