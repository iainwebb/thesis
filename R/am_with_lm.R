# HIDE -------
# i <- 36
# k <- 9
# for (i in 1:p) {
#   for (k in 1:N) {
#     d_dx_i_t_x_star_T_matrix <- matrix(-2 / l_hat_matrix[i,]^2 * (rep(x_star_matrix[i,k],n) - x_norm_matrix[i,]) * t_x_star_T_matrix[k,], ncol = n)
#     partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix %*% A_inv_matrix %*% (AODs_Total_jul_level2_gb - H_matrix %*% betas_hat_matrix)))
#   }
# }
# 
# d_dx_i_t_x_star_T_matrix <- matrix(-2 / l_hat_matrix[i,]^2 * (rep(x_star_matrix[i,k],n) - x_norm_matrix[i,]) * t_x_star_T_matrix[k,], ncol = n)
# 
# 
# partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe
# 
# dim(abs(partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe))

# START ------
# setwd("U:/ManWin/My Documents/thesis")
setwd("C:/Users/smp22ijw/Desktop")
rm(list = ls())
calcs <- 1

library(ncdf4)
library(DiceKriging)
library(lhs)

# AOD_Total work ----

JDO_jul_nc <- nc_open(paste("data/ACURE_P3_AOD_Total_jul.nc"))
JDO_jul <- ncvar_get(JDO_jul_nc, names(JDO_jul_nc$var))
AODs_Total_jul_level2 <- JDO_jul[,,2,]

# R2_adj_matrix <- matrix(NA, nrow = 192, ncol = 144)

for (g1 in 1:192) {
  for (g2 in 1:144) {
    df <- data.frame(cbind("AOD_Total_jul" = AODs_Total_jul_level2[g1,g2,], x_norm_T_matrix))
    # R2[g1,g2] <- summary(lm(AOD_Total_jul ~., df))$r.squared
    R2_adj[g1,g2] <- summary(lm(AOD_Total_jul ~., df))$adj.r.squared
  }
  cat(g1, "/192, ", sep = "")
}
range(R2_adj_matrix)
hist(R2_adj_matrix)
library(gplots)
heatmap.2(t(R2_adj_matrix), Rowv = NA, Colv = NA)
heatmap(t(R2_adj_matrix), Rowv = NA, Colv = NA)
R2_adj_rounded_matrix <- round(R2_adj_matrix,1)
R2_adj_factor_matrix <- factor(R2_adj_rounded_matrix)
AOD_Total_df <- expand.grid(1:192, 1:144)
AOD_Total_df$R2_adj <- R2_adj_factor_matrix
table(R2_adj_factor_matrix)
library(ggplot2)
ggplot(AOD_Total_df, aes(Var1, Var2, fill= R2_adj)) + 
  geom_tile(color = "white",
            lwd = 0.5,
            linetype = 1)+
  theme_bw() +
  scale_fill_manual(values = c('yellow', 'yellow', 'yellow', 'yellow', 
                               'yellow', 'yellow', 'yellow', 
                               'magenta', 
                               'black', 'black'))


library(fields)
library(maps)

par(mfrow = c(1,1))

longitude = ncvar_get(JDO_jul_nc,"longitude") - 180
latitude = ncvar_get(JDO_jul_nc,"latitude")

image.plot(longitude, latitude,
           rbind(R2_adj_matrix[97:192,],
                 R2_adj_matrix[1:96,]),
           xlab = "longitude", ylab = "latitude")
map("world",lwd=1.2,add=TRUE)

colour_breaks_vector <- c(0, 0.8, 1)
breaks_number <- length(colour_breaks_vector)
midpoints_vector <- seq(1.5, 1.5 + breaks_number - 2, 1)
colours_number <- length(midpoints_vector)
colours_vector <- as.vector(tim.colors(colours_number))
R2_adj_matrix_for_plotting <- R2_adj_matrix
R2_adj_matrix_for_plotting <- 
  ifelse(R2_adj_matrix_for_plotting < colour_breaks_vector[2], midpoints_vector[1], 
         ifelse(R2_adj_matrix_for_plotting < colour_breaks_vector[3], midpoints_vector[2]
                ))

100 * (sum(R2_adj_matrix < 0.7) - sum(R2_adj_matrix < 0.6)) / (192*144)
image.plot(longitude, latitude,
           rbind(R2_adj_matrix[97:192,],
                 R2_adj_matrix[1:96,]),
           breaks = c(0,0.6,0.7,0.8,1),
           col = c("red", "pink", "lavenderblush", "white"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.6,0.7,0.8,1),labels=as.character(c(0,0.6,0.7,0.8,1)),mgp=c(3,0.5,0)
           )
           )
map("world",lwd=1.2,add=TRUE, lty=3)

R2 <- round(R2_adj_matrix, 1)
R2 <- ifelse(R2 <= 0.6, 0.6, R2)
table(R2)
round(table(R2) / sum(table(R2)),2)

image.plot(longitude,
           latitude,
           PlotData,
           xaxs='i',
           ylab="Latitude", xlab="Longitude", cex.lab=1.3,
           zlim=c(1,length(ColBreakVec)),
           axis.args=list(
             at=seq(1,length(ColBreakVec),1),labels=as.character(ColBreakVec),mgp=c(3,0.5,0)
           ),
           breaks=seq(1,length(ColBreakVec),1),
           col=IWColoursVec,
           cex.axis=1,
           main=PlotTitleIn,cex.main=1.5,
           legend.args=list( text=LegendTitleIn, col=1, cex=1.2, side=4, line=2.8))
map("world",lwd=1.2,add=TRUE)


ColBreakVec <- c(0, 10, 20, 50, 100, 200, 500, 1000, 3000, 5000, 8000)
NumBreaks <- length(ColBreakVec)
MidsVec <- seq(1.5, 1.5+NumBreaks-2, 1)
NumCols <- length(MidsVec)

#### Conversion of PlotData from CCN to MidsVec elements

PlotData <- 
  ifelse(PlotData < ColBreakVec[2], MidsVec[1], 
         ifelse(PlotData < ColBreakVec[3], MidsVec[2], 
                ifelse(PlotData < ColBreakVec[4], MidsVec[3], 
                       ifelse(PlotData < ColBreakVec[5], MidsVec[4], 
                              ifelse(PlotData < ColBreakVec[6], MidsVec[5], 
                                     ifelse(PlotData < ColBreakVec[7], MidsVec[6], 
                                            ifelse(PlotData < ColBreakVec[8], MidsVec[7], 
                                                   ifelse(PlotData < ColBreakVec[9], MidsVec[8], 
                                                          ifelse(PlotData < ColBreakVec[10], MidsVec[9], 
                                                                 ifelse(PlotData < ColBreakVec[11], MidsVec[10]
                                                                 ))))))))))
# table(PlotData)
# nrow(PlotData) * ncol(PlotData)
# length(PlotData)
# hist(PlotData)
```

```{r median-ccn-map-plot, fig.cap = "Annual Mean CCN for the median simulator run.", message=F}

#### Set up labels

LegendTitleIn = expression("Annual mean CCN (cm"^-3*")")
# "ANN CCN Conc. (cm^-3)"
# PlotTitleIn = paste("Annual Mean CCN: Run ",EnsembleMemberNo,sep="")
PlotTitleIn=""

#### Make the colour scale:

IWColoursVec<-as.vector(tim.colors(NumCols))


# ERF work ----


erf_jul <- read.table("data/ERF_PPE_global_jul.dat", col.names = "ERF_jul")

# heatmap(t(AODs_Total_jul_level2[,,125]), Rowv = NA, Colv = NA)

# SETUP
# specify the inputs ----------------------
# x_norm_T_matrix <- read.csv("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")
x_norm_T_matrix <- read.csv("data/ACURE-UKESM1-PPE-par-norm-Design.csv")
#                                                                       [,1:4]

forcing_dataframe <- data.frame(cbind(erf_jul, x_norm_T_matrix))
# R2[g1,g2] <- summary(lm(AOD_Total_jul ~., df))$r.squared
forcing_R2 <- summary(lm(ERF_jul ~., forcing_dataframe))$adj.r.squared

x_norm_matrix <- t(x_norm_T_matrix)
p <- ncol(x_norm_T_matrix)
colnames(x_norm_T_matrix) <- paste0(rep("x.",p), 1:p)

# number of observations ------------------
n <- nrow(x_norm_T_matrix)

# PRIOR MEAN FUNCTION
q <- 1 + p
# and corresponding formula
h_formula <- ~.

# DERIVATIVE WORK
# defining covariance matrix making function
sq_exp_cov_function <- function(matrix_1, matrix_2, l = l_hat_diag_matrix){
  cov_matrix <- matrix(rep(NA, ncol(matrix_1) * ncol(matrix_2)), nrow = ncol(matrix_1))
  for (a in 1:nrow(cov_matrix)) {
    for (b in 1:ncol(cov_matrix)) {
      cov_matrix[a,b] <- exp(-matrix((matrix_1[,a] - matrix_2[,b]), nrow = 1) %*% l_hat_diag_matrix %*% matrix((matrix_1[,a] - matrix_2[,b]), ncol = 1))
    }
  }
  return(cov_matrix)
}
library(fields)
corGaussian <- function(inputs, inputs2, phi) {
  
  if (missing(inputs2) || is.null(inputs2))
    return(corGaussianSquare(inputs, phi))
  
  delta <- (phi)
  exp(-(rdist(inputs / rep(delta, each = nrow(inputs)), inputs2 / rep(delta, each = nrow(inputs2))) ^ 2))
}
# define x_star  ---------------------------
N <- 200000
x_star_matrix <- t(randomLHS(N, p))
# x_star_matrix <- t(read.csv("data/x_norm_matrix.csv", header = FALSE, sep = ""))[,sample(1:200000, 10000)]
N <- ncol(x_star_matrix)
x_star_T_dataframe <- data.frame(t(x_star_matrix))
colnames(x_star_T_dataframe) <- paste0(rep("x.",p), 1:p)
# make H_matrix
H_matrix <- cbind(c(rep(1, n)), t(unname(x_norm_matrix)))

# erf gp ----

# fit the GP ------------------------------
f_GP <- km(~., 
           design = x_norm_T_matrix, 
           response = erf_jul, 
           covtype="gauss", optim.method="BFGS", control=list(maxit=500))
# extract estimates
betas_hat_matrix <- matrix(f_GP@trend.coef,
                           nrow = q)
sigma_sq_hat <- f_GP@covariance@sd2
l_hat_vector <- f_GP@covariance@range.val
l_hat_matrix <- matrix(f_GP@covariance@range.val,
                       nrow = p)
l_hat_diag_matrix <- diag(as.vector(l_hat_matrix))
assign(paste0("betas_with_", N, "_gb_", g, "_attempt_", c), f_GP@trend.coef)
assign(paste0("l_hat_with_", N, "_gb_", g, "_attempt_", c), f_GP@covariance@range.val)

# DERIVATIVE WORK
# predict at x_star values using GP --------
x_star_predictions_list <- predict(f_GP,
                                   newdata = x_star_T_dataframe,
                                   type="SK"
)
# create A_matrix
# A_matrix <- 
#   sq_exp_cov_function(x_norm_matrix, x_norm_matrix) 
# + diag(0.001, nrow(x_norm_T_matrix))
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
colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx.",p), 1:p)

# make t(x_star)^T
# t_x_star_T_matrix <- sq_exp_cov_function(x_star_matrix, x_norm_matrix)
t_x_star_T_matrix <- corGaussian(t(x_star_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))

i <- sample(1:37, 1)
k <- N/2
for (i in 1:p) {
  d_dx_i_t_x_star_T_matrix <- -2 / l_hat_matrix[i,]^2 * (matrix(rep(x_star_matrix[i,], n), ncol = n) - matrix(rep(x_norm_matrix[i,], N), ncol = n, byrow = T)) * t_x_star_T_matrix[,]
  for (k in 1:N) {
    partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix[k,] %*% A_inv_matrix %*% (AODs_Total_jul_level2_gb - H_matrix %*% betas_hat_matrix)))
  }
}

partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
assign(paste0("partial_derivatives_normalised_response_surface_erf_jul_dataframe"), partial_derivatives_normalised_dataframe)

# AOD_jul_part_der_array <- array(NA, dim = c(192, 144, 37))

for (g1 in 1:192) {
  for (g2 in 1:144) {
    df <- data.frame(cbind("AOD_Total_jul" = AODs_Total_jul_level2[g1,g2,], x_norm_T_matrix))
    AOD_jul_part_der_array[g1,g2,] <- summary(lm(AOD_Total_jul ~., df))$coefficients[2:(p+1),1]
  }
  cat(g1, "/192, ", sep = "")
}

sink("AOD_jul_part_der_vector.txt") # all output goes into this file
a <- as.vector(AOD_jul_part_der_array)
a
sink()

sink("AOD_jul_part_der_array.txt") # all output goes into this file
a <- AOD_jul_part_der_array
a
sink() # turn off sink
test_array <- readLines("AOD_jul_part_der_array.txt")

# AM_AOD_erf_matrix <- matrix(NA, nrow = 192, ncol = 144)

for (g1 in 1:192) {
  for (g2 in 1:144) {
    AM_AOD_erf_matrix[g1,g2] <- 
      sum(abs(rowSums(partial_derivatives_normalised_response_surface_erf_jul_dataframe * 
                        AOD_jul_part_der_array[g1,g2,]))
      ) / N
  }
  cat(g1, "/192, ", sep = "")
}

AM_AOD_erf_matrix_for_plotting <- rbind(AM_AOD_erf_matrix[97:192,],
                                        AM_AOD_erf_matrix[1:96,])
par(mfrow = c(1, 1))
image.plot(longitude, latitude,
           AM_AOD_erf_matrix_for_plotting,
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
map("world",lwd=1.2,add=TRUE, lty=3, col = "white")

hist(AM_AOD_erf_matrix_for_plotting)

image.plot(longitude, latitude,
           AM_AOD_erf_matrix_for_plotting,
           breaks = c(0,0.01,0.02,0.03,0.04,0.05,
                      0.06,0.07,0.08,0.09,0.2,1),
           col = c("red4","red3","red1","indianred1", "mistyrose",
                   "honeydew1", "darkseagreen1", "palegreen2", "seagreen3", 
                   "springgreen3","springgreen4"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.01,0.02,0.03,0.04,0.05,
                  0.06,0.07,0.08,0.09,0.2,1),labels=as.character(c(0,0.01,0.02,0.03,0.04,0.05,
                                                               0.06,0.07,0.08,0.09,0.2,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=3, col = "white")


# AOD gp ----

gridboxes <- matrix(c(20, 120, 20, 120), nrow = 2, byrow = T)
# gridboxes <- matrix(c(20, 120, 20, 121), nrow = 2, byrow = T)
# gridboxes <- matrix(c(20, 120, 100, 85), nrow = 2, byrow = T)

c <- 1
g <- 1

for (c in 1:calcs) {
  
  start_time <- Sys.time()
  
  for (g in 1:2) {
    # SETUP
    # observations ----------------------------
    AODs_Total_jul_level2_gb <- AODs_Total_jul_level2[gridboxes[g,1], gridboxes[g,2],]
    
    # GP
    # fit the GP ------------------------------
    f_GP <- km(~., 
               design = x_norm_T_matrix, 
               response = AODs_Total_jul_level2_gb, 
               covtype="gauss", optim.method="BFGS", control=list(maxit=500))
    # extract estimates
    betas_hat_matrix <- matrix(f_GP@trend.coef,
                               nrow = q)
    sigma_sq_hat <- f_GP@covariance@sd2
    l_hat_vector <- f_GP@covariance@range.val
    l_hat_matrix <- matrix(f_GP@covariance@range.val,
                           nrow = p)
    l_hat_diag_matrix <- diag(as.vector(l_hat_matrix))
    assign(paste0("betas_with_", N, "_gb_", g, "_attempt_", c), f_GP@trend.coef)
    assign(paste0("l_hat_with_", N, "_gb_", g, "_attempt_", c), f_GP@covariance@range.val)
    
    # DERIVATIVE WORK
    # predict at x_star values using GP --------
    x_star_predictions_list <- predict(f_GP,
                                       newdata = x_star_T_dataframe,
                                       type="SK"
    )
    # create A_matrix
    # A_matrix <- 
    #   sq_exp_cov_function(x_norm_matrix, x_norm_matrix) 
    # + diag(0.001, nrow(x_norm_T_matrix))
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
    colnames(partial_derivatives_dataframe) <- paste0(rep("d_dx.",p), 1:p)
    
    # make t(x_star)^T
    # t_x_star_T_matrix <- sq_exp_cov_function(x_star_matrix, x_norm_matrix)
    t_x_star_T_matrix <- corGaussian(t(x_star_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))

    i <- sample(1:37, 1)
    k <- N/2
    for (i in 1:p) {
      d_dx_i_t_x_star_T_matrix <- -2 / l_hat_matrix[i,]^2 * (matrix(rep(x_star_matrix[i,], n), ncol = n) - matrix(rep(x_norm_matrix[i,], N), ncol = n, byrow = T)) * t_x_star_T_matrix[,]
      for (k in 1:N) {
        partial_derivatives_dataframe[k,i] <- as.vector(unlist(x_star_matrix[i,k] * betas_hat_matrix[i+1,] - 2 / l_hat_matrix[i,] * d_dx_i_t_x_star_T_matrix[k,] %*% A_inv_matrix %*% (AODs_Total_jul_level2_gb - H_matrix %*% betas_hat_matrix)))
      }
    }
    
    partial_derivatives_normalised_dataframe <- partial_derivatives_dataframe / sqrt(rowSums(partial_derivatives_dataframe^2))
    assign(paste0("partial_derivatives_normalised_response_surface_", g, "_dataframe"), partial_derivatives_normalised_dataframe)
    
  }
  
  AM <- sum(abs(rowSums(partial_derivatives_normalised_response_surface_1_dataframe * partial_derivatives_normalised_response_surface_2_dataframe))) / N
  
  # assign(paste0("AM_with_", N, "_attempt_", c), AM)
  write(AM, 
        file=paste0("R/", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_AM_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(AM, 
        file=paste0("R/", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("betas_with_", N, "_gb_", 1, "_attempt_", c)))), 
        file=paste0("R/", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("betas_with_", N, "_gb_", 2, "_attempt_", c)))), 
        file=paste0("R/", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 1, "_attempt_", c)))), 
        file=paste0("R/", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  write(print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 2, "_attempt_", c)))), 
        file=paste0("R/", 
                    gridboxes[1,1], "_", gridboxes[1,2], "_", 
                    gridboxes[2,1], "_", gridboxes[2,2],
                    "_all_", N, "_pts_", p, "_pars.txt"),append=TRUE)
  
  end_time <- Sys.time()
  
  # assign(paste0("time_with_", N, "_attempt_", c), end_time - start_time)
  
}

# for (c in 1:calcs) {
#   print(eval(parse(text = paste0("AM_with_", N, "_attempt_", c))))
#   print(eval(parse(text = paste0("betas_with_", N, "_gb_", 1, "_attempt_", c))))
#   print(eval(parse(text = paste0("betas_with_", N, "_gb_", 2, "_attempt_", c))))
#   print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 1, "_attempt_", c))))
#   print(eval(parse(text = paste0("l_hat_with_", N, "_gb_", 2, "_attempt_", c))))
# }
# 
# for (c in 1:calcs) {
#   print(eval(parse(text = paste0("time_with_", N, "_attempt_", c))))
# }
# 
# print(c(p, N))

corGaussian <- function(inputs, inputs2, phi) {
  
  if (missing(inputs2) || is.null(inputs2))
    return(corGaussianSquare(inputs, phi))
  
  delta <- (phi)
  exp(-(rdist(inputs / rep(delta, each = nrow(inputs)), inputs2 / rep(delta, each = nrow(inputs2))) ^ 2))
}

sq_exp_cov_function <- function(matrix_1, matrix_2, l = l_hat_diag_matrix){
  cov_matrix <- matrix(rep(NA, ncol(matrix_1) * ncol(matrix_2)), nrow = ncol(matrix_1))
  for (a in 1:nrow(cov_matrix)) {
    for (b in 1:ncol(cov_matrix)) {
      cov_matrix[a,b] <- exp(-matrix((matrix_1[,a] - matrix_2[,b]), nrow = 1) %*% l_hat_diag_matrix %*% matrix((matrix_1[,a] - matrix_2[,b]), ncol = 1))
    }
  }
  return(cov_matrix)
}
t_x_star_T_matrix <- sq_exp_cov_function(x_star_matrix, x_norm_matrix)

library(fields)
t_x_star_T_matrix_2 <- corGaussian(t(x_star_matrix), t(x_norm_matrix), 1/sqrt(l_hat_vector))

sum(abs(t_x_star_T_matrix - t_x_star_T_matrix_2) > 0.0001)
