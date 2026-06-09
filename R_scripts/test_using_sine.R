setwd("U:/ManWin/My Documents/thesis")

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
n <- 5
inputs_x_norm_T_matrix <- matrix(seq(0.1, 0.9, length.out=5))
inputs_x_norm_matrix <- t(inputs_x_norm_T_matrix)

## There are p variables in all, and n runs
p <- ncol(inputs_x_norm_T_matrix)
n <- nrow(inputs_x_norm_T_matrix)

# Importing the outputs ########################################################

## ERF is the unobserved variable, imported as 221 × 1 dataframe, with the 
## local version first
# response_ERF_dataframe <- read.table("data/ERF_PPE_global_jul.dat", col.names = "ERF")

response_test_dataframe <- data.frame(y = matrix(2*sin(3*inputs_x_norm_T_matrix - 0.5)))
plot(response_test_dataframe$y ~ as.vector(inputs_x_norm_matrix), xlim = c(0,1))
# GP set-up ####################################################################

## We specify some variables, namely N (the size of the sample from the input
## parameter space) and q (the number of betas that will be estimated) 

N <- 100
q <- 1 + p

## We load the required packages for sampling
library(lhs)

## We create the sample
xstar_matrix <- t(matrix(seq(0, 1, length.out=N)))

## Our estimated posterior derivative requires H, which doesn't depend on the 
# response (and thus on the gridbox in the case of AOD)
H_matrix <- cbind(c(rep(1, n)), t(unname(inputs_x_norm_matrix)))

## We define a function called gp_d_dxi_normalised_function which returns the
## normalised GP estimates of the partial derivatives. It takes as its single
## argument a 221-long vector or 221 x 1 dataframe of responses

source("./R_scripts/1_gp_d_dxi_normalised_function.R")

## First we'll fit the GP for ERF and calculate the normalised partial
## derivatives
gppdnormalised_t_test_matrix <- gp_d_dxi_normalised_function(response_test_dataframe)

## This produces the following AMs (with the lm-estimated partial derivatives
## for AOD):
AM_ERFgp_AODlmall_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_AODlmall_matrix[long,lat] <-
      mean(abs(colSums(gppdnormalised_t_ERF_matrix * matrix(rep(pdsnormed_AOD_lmall_array[long,lat,], N), ncol = N))))
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

AM_ERFgp_AODlmbest_matrix <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_AODlmbest_matrix[long, lat] <-
      if (R2adj_AOD_matrix[long, lat] < threshold)
        -0.05
    else AM_ERFgp_AODlmall_matrix[long, lat]
  }
}
sum(is.na(AM_ERFgp_AODlmbest_matrix))

## We produce a map to show these AMs, excluding the one for which the adjusted
## R^2 was too low
image.plot(longitude, latitude,
           rbind(AM_ERFgp_AODlmbest_matrix[97:192,],
                 AM_ERFgp_AODlmbest_matrix[1:96,]),
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
           rbind(AM_ERFgp_AODlmall_matrix[97:192,],
                 AM_ERFgp_AODlmall_matrix[1:96,]),
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

# AM_ERFgp_AODlmbestgprest_matrix <- matrix(NA, nrow = 192, ncol = 144)

index_counter <- 0

for (long in 1:192) {
  for (lat in 1:144) {
    AM_ERFgp_AODlmbestgprest_matrix[long, lat] <-
      if (R2adj_AOD_matrix[long, lat] < threshold)
        mean(abs(colSums(gppdnormalised_t_ERF_matrix * gp_d_dxi_normalised_function(response_AOD_array[long, lat,]))))
    else AM_ERFgp_AODlmall_matrix[long, lat]
    index_counter <-
      if (R2adj_AOD_matrix[long, lat] < threshold)
        index_counter + 1
    else index_counter
    if (R2adj_AOD_matrix[long, lat] < threshold)
      write_lines(as.vector(AM_ERFgp_AODlmbestgprest_matrix),
                  file=paste0("objects/AM_ERFgp_AODlmbestgprest_matrix_0", threshold*100, "_", N, ".txt"))
  }
  write_lines(as.vector(AM_ERFgp_AODlmbestgprest_matrix),
              file=paste0("objects/AM_ERFgp_AODlmbestgprest_matrix_0", threshold*100, "_", N, ".txt")) 
}
# rm("long", "lat", "p", "q", "x_star_matrix", "inputs_x_norm_matrix", "x_star_T_dataframe",
#    "n", "H_matrix", "corGaussian", "gp_d_dxi_normalised_function", "index_counter")

sum(is.na(AM_ERFgp_AODlmbestgprest_matrix))

range(AM_ERFgp_AODlmbestgprest_matrix)

table(floor(AM_ERFgp_AODlmall_matrix * 10) / 10)
table(floor(AM_ERFgp_AODlmbest_matrix * 10) / 10)
table(floor(AM_ERFgp_AODlmbestgprest_matrix * 10) / 10)

100*round(
  table(
    ifelse(floor(R2adj_AOD_matrix * 10) / 10 < 0.6, "<0.6",
           floor(R2adj_AOD_matrix * 10) / 10)
  ) / (192*144), 3)

hist(R2adj_AOD_matrix, xlim = c(0,1))

image.plot(longitude, latitude,
           rbind(AM_ERFgp_AODlmbestgprest_matrix[97:192,],
                 AM_ERFgp_AODlmbestgprest_matrix[1:96,]),
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
R2adj_AMalllm <- cbind("AOD_R2adj" 
                       = as.vector(R2adj_AOD_matrix), 
                       "AM_AOD_versus_ERFalllm_matrix" 
                       = as.vector(AM_ERFgp_AODlmall_matrix)
)
blue_black_colours <- ifelse(as.vector(R2adj_AOD_matrix) < 0.60, "blue", "black")
plot(R2adj_AMalllm, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours,
     xlab = "R2_adj for AOD_Total",
     ylab = "AM AOD_Total versus ERF")
R2adj_AM <- cbind("AOD_Total_R2_adj" 
                  = as.vector(R2adj_AOD_matrix), 
                  "AM_for_AOD_Total_versus_ERF" 
                  = as.vector(AM_ERFgp_AODlmbestgprest_matrix)
)
plot(R2adj_AM, xlim = c(0,1), ylim = c(0,1), pch=15, cex=0.5,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours,
     xlab = "R2_adj for AOD_Total",
     ylab = "AM AOD_Total versus ERF")


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
