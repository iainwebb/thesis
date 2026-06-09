library(RobustGaSP)
library(knitr)
library(readr)
library(lhs)
library(tidyverse)
library(ncdf4)
library(fields)
library(maps)
library(kableExtra)
library(ggplot2)

setwd("U:/ManWin/My Documents/thesis/AM")

source("H_maker.R")
source("estimated_parameters.R")
source("estimated_parameters_with_nugget.R")
source("A_maker.R")
source("t_maker.R")
source("loo_cv.R")
source("loo_validation_plot.R")
source("byhand_preds.R")
source("byhand_plots.R")

### Set-up #####################################################################

# option selection ####

output_variable1 <- "AOD"

# if looking at gridbox with itself
gridbox_with_another <- "no"
output_variable2 <- "ACI"
long1 <- long2 <- 5; lat1 <- lat2 <- 114

## if looking at 2 different gridboxes
# gridbox_with_another <- "yes"
# output_variable2 <- output_variable1
# long1 <- 5; lat1 <- 114
# long2 <- 6; lat2 <- 114
# long2 <- 45; lat2 <- 55

## month of interest
month <- 7

## whether or not to remove the 25 problematic runs for which prim_so4_diam < 0.1
prim_so4_diam_removal <- "no"

## whether or not to perform SA, and required minimum sum of main effects
perform_SA <- "no"
SA_main_effects_min <- if(perform_SA == "yes"){0.95}
# 1.204
if(perform_SA == "no"){rm(SA_main_effects_min)}
nugget_estimation <- ifelse(perform_SA == "yes", TRUE, FALSE)
nugget_value <- if(perform_SA == "no"){1e-5}
if(perform_SA == "yes"){rm(nugget_value)}

## parameter estimate method
em <- "post_mode"
# em <- "mle"

# map things ####

response_AOD_nc <- nc_open(
  paste0("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_", 
         tolower(month.abb[month]),
         ".nc")
)
## Extract longitudes and latitudes, adjusting the longitudes ready for producing maps
longitude = ncvar_get(response_AOD_nc,"longitude") - 180
latitude = ncvar_get(response_AOD_nc,"latitude")
if(output_variable1 != "AOD" & output_variable2 != "AOD"){rm(response_AOD_nc)}

# inputs ####

input <- 
  as.matrix(read.csv(
    "X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")
  )
prim_so4_diam_low_indices <- which(input[,7] >= 0.1)
if(prim_so4_diam_removal == "yes"){input <- input[prim_so4_diam_low_indices,]}
input_names <- colnames(input)
dim_inputs <- length(input_names)
num_obs <- nrow(input)
H <- H_maker(input)

# test inputs ####

num_testing_input <- 10000
# generate points to be predicted
# testing_input <- maximinLHS(n=num_testing_input, k=dim_inputs)
# testing_input <- matrix(runif(10000*37), ncol = 37)
## Export 
# write_lines(as.vector(testing_input),
#             file=
#               paste0(
#                 "U:/ManWin/My Documents/thesis/Rmd/objects_for_Rmd/testing_input_",
#                 dim_inputs,
#                 "D_", num_testing_input, "newpts.txt")
#             )
## Import
testing_input <-
  matrix(as.numeric(readLines(paste0("testing_input_", dim_inputs, "D_", num_testing_input, "newpts.txt"))),
         ncol = dim_inputs)
if(prim_so4_diam_removal == "yes") {testing_input[,7] <- 0.1 + testing_input[,7] * 0.9}
H_testing <- H_maker(testing_input)

# outputs ####

if(output_variable1 == "AOD"){response_AOD_array <- ncvar_get(response_AOD_nc, names(response_AOD_nc$var))[,,2,]; output1 <- matrix(response_AOD_array[long1,lat1,], nrow = length(response_AOD_array[long1,lat1,]))}
if(output_variable2 == "AOD"){response_AOD_array <- ncvar_get(response_AOD_nc, names(response_AOD_nc$var))[,,2,]; output2 <- matrix(response_AOD_array[long2,lat2,], nrow = length(response_AOD_array[long2,lat2,]))}
if(output_variable2 == "ACI"){output2 <- as.matrix(
  read.table("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Optimal-Constraint-Data/data_for_emulation/jul/ACI_PPE_global_jul.dat",
             col.names = "ACI"
  )
)
}
if(output_variable2 == "ARI"){output2 <- as.matrix(
  read.table("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Optimal-Constraint-Data/data_for_emulation/jul/ARI_clear_PPE_global_jul.dat",
             col.names = "ARI"
  )
)
}
colnames(output1) <- output_variable1
colnames(output2) <- output_variable2

if(prim_so4_diam_removal == "yes"){output1 <- matrix(output1[prim_so4_diam_low_indices,], nrow=num_obs); output2 <- matrix(output2[prim_so4_diam_low_indices,], nrow=num_obs)}

rm(response_AOD_nc, response_AOD_array)

# type of analysis ####

type_of_analysis_statement <- ifelse(gridbox_with_another == "yes",
                                     paste(
                                       output_variable1, "at the pair of gridboxes at longitude", longitude[long1], "and latitude", latitude[lat1], "and longitude", longitude[long2], "and latitude", latitude[lat2], "(shown in green and blue respectively in Figure 1) are considered, for", month.name[month]
                                     ),
                                     paste(
                                       output_variable1, "and", output_variable2, "at the gridbox at longitude", longitude[long1], "and latitude", latitude[lat1], "(shown in green in Figure 1) are considered, for", month.name[month]
                                     )
)
rm(month, output_variable1, output_variable2)
type_of_analysis_statement

# prim_so4_diam removals ####

prim_so4_diam_statement <- ifelse(prim_so4_diam_removal == "yes",
                                  paste("The", 221 - length(prim_so4_diam_low_indices), "runs for which prim_so4_diam < 0.1 have been removed."),
                                  paste("The", 221 - length(prim_so4_diam_low_indices), "runs for which prim_so4_diam < 0.1 have not been removed.")
)
prim_so4_diam_statement

# remove objects ####

rm(prim_so4_diam_statement, type_of_analysis_statement)

# map of selected gridboxes ####

map_matrix <- matrix(NA, nrow = 192, ncol = 144)
map_matrix[long1, lat1] <- 0
if(gridbox_with_another == "yes"){map_matrix[long2, lat2] <- 1}

pdf('gridboxes.pdf')
image(longitude, latitude, 
      rbind(map_matrix[97:192,],
            map_matrix[1:96,]),
      breaks = c(0,0.5,1),
      col = c("seagreen3","blue"),
      xlab = "longitude", ylab = "latitude"
)
maps::map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

### Sensitivity analysis #######################################################

# SA for output 1 ####

if(perform_SA == "yes"){
  input_output1 <- cbind(input, output1)
  ## Plots if needed
  # input_output_ordered <- input_output %>% 
  #   as_tibble %>%
  #   gather('parameter', 'value', -AOD)
  # input_output_ordered <- transform(input_output_ordered,
  #                        parameter=factor(parameter, levels = c(input_names)))
  # input_output_ordered %>%
  #   ggplot(aes(x=value, y=AOD)) +
  #   geom_point(alpha = 0.5) +
  #   facet_wrap(~parameter) +
  #   labs(y='output', x='input') +
  #   geom_smooth(method = mgcv::gam,
  #               formula = y ~ s(x, bs = "tp"),
  #               fill = "red",
  #               method.args = list(method="GCV.Cp"))
  main_effects1 <- rep(0, dim_inputs)
  for(i in 1:dim_inputs){
    gam1 <- mgcv::gam(input_output1[, dim_inputs+1] ~ s(input_output1[, i]))
    main_effects1[i] <- var(gam1$fitted) / var(input_output1[, dim_inputs+1])
  }
  # par(mar=c(6.5,3,1,0.5))
  # barplot(main_effects1,
  #         # names.arg = "",
  #         names.arg = c(input_names),
  #         las=2, cex.names=0.75,
  #         ylim = c(0,max(main_effects1)+0.1)
  # )
  main_effects_sorted1 <- sort(main_effects1, decreasing = T)
  main_effects_sum1 <- 0
  main_effects_num_vars_1 <- 0
  for (i in 1:dim_inputs) {
    if(main_effects_sum1 < SA_main_effects_min)
    {main_effects_sum1 <- main_effects_sum1 + main_effects_sorted1[i]; 
    main_effects_num_vars_1 <- main_effects_num_vars_1 + 1}
  }
  main_effects_vars1 <- data.frame(cbind(input_names, main_effects1))[order(data.frame(cbind(input_names, main_effects1))$main_effects1, decreasing = T)[1:main_effects_num_vars_1],1]
  
  rm(gam1, input_output1)
}
  
# SA for output 2 ####
  
if(perform_SA == "yes"){
  input_output2 <- cbind(input, output2)
  ## Plots if needed
  # input_output_ordered <- input_output %>% 
  #   as_tibble %>%
  #   gather('parameter', 'value', -AOD)
  # input_output_ordered <- transform(input_output_ordered,
  #                        parameter=factor(parameter, levels = c(input_names)))
  # input_output_ordered %>%
  #   ggplot(aes(x=value, y=AOD)) +
  #   geom_point(alpha = 0.5) +
  #   facet_wrap(~parameter) +
  #   labs(y='output', x='input') +
  #   geom_smooth(method = mgcv::gam,
  #               formula = y ~ s(x, bs = "tp"),
  #               fill = "red",
  #               method.args = list(method="GCV.Cp"))
  main_effects2 <- rep(0, dim_inputs)
  for(i in 1:dim_inputs){
    gam2 <- mgcv::gam(input_output2[, dim_inputs+1] ~ s(input_output2[, i]))
    main_effects2[i] <- var(gam2$fitted) / var(input_output2[, dim_inputs+1])
  }
  # par(mar=c(6.5,3,1,0.5))
  # barplot(main_effects2,
  #         # names.arg = "",
  #         names.arg = c(input_names),
  #         las=2, cex.names=0.75,
  #         ylim = c(0,max(main_effects2)+0.1)
  # )
  main_effects_sorted2 <- sort(main_effects2, decreasing = T)
  main_effects_sum2 <- 0
  main_effects_num_vars_2 <- 0
  for (i in 1:dim_inputs) {
    if(main_effects_sum2 < SA_main_effects_min)
    {main_effects_sum2 <- main_effects_sum2 + main_effects_sorted2[i]; 
    main_effects_num_vars_2 <- main_effects_num_vars_2 + 1}
  }
  main_effects_vars2 <- data.frame(cbind(input_names, main_effects2))[order(data.frame(cbind(input_names, main_effects2))$main_effects2, decreasing = T)[1:main_effects_num_vars_2],1]
  
  rm(gam2, input_output2)
}
  
# SA statements ####
  
SA_statement <- ifelse(perform_SA == "yes",
                       paste0(
                         "GAM-based sensitivity analysis is performed for the two outputs. To reach ", 100*SA_main_effects_min, "% of the variance of the two outputs requires ", main_effects_num_vars_1, " and ", main_effects_num_vars_2, " inputs respectively (with main effects summing to ", round(main_effects_sum1 * 100), "% and ", round(main_effects_sum2 * 100), "%); see Table 1 and Figure 2."
                       ),
                       paste0("Sensitivity analysis is not carried out")
)
SA_union_statement <- ifelse(perform_SA == "yes",
                             paste0(
                               "An emulator will be fitted for the union of inputs from Table 1, namely the ", length(union(main_effects_vars1, main_effects_vars2)), " inputs:"),
                             ""
)
SA_statement

# SA tables ####

if(perform_SA == "yes"){
  k1 <- data.frame(cbind(main_effects_vars1, 100*round(main_effects_sorted1[1:main_effects_num_vars_1],2)))
  k2 <- data.frame(cbind(main_effects_vars2, 100*round(main_effects_sorted2[1:main_effects_num_vars_2],2)))
  p1 <- data.frame(cbind(main_effects_vars1, round(main_effects_sorted1[1:main_effects_num_vars_1],2)))
  p2 <- data.frame(cbind(main_effects_vars2, round(main_effects_sorted2[1:main_effects_num_vars_2],2)))
  colnames(p1) <- colnames(p2) <- c("Input", "Main effect")
  
  kable_styling(
    kable(
      list(k1, k2),
      format = "latex",
      caption = 'Post-SA main effects for the first (left) and second (right) outputs.',
      col.names = c("Input", "Main effect (%)"),
      align=c("l","c")
    ),
    latex_options = "HOLD_position")
}

# post-SA inputs and objects ####

if(perform_SA == "yes"){
  SA_indices <- which(main_effects1 >= main_effects_sorted1[main_effects_num_vars_1] | main_effects2 >= main_effects_sorted2[main_effects_num_vars_2])
  
  input <- input[,SA_indices]
  
  dim_inputs <- ncol(input)
  
  testing_input <- testing_input[,SA_indices]
  
  H <- H_maker(input)
  H_testing <- H_maker(testing_input)
  
  input_names <- union(main_effects_vars1, main_effects_vars2)
}

# SA pie charts ####

if(perform_SA == "yes"){
  c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown"
  )
  
  # p <- merge(p1, p2, by="Input", all=T)
  # colnames(p) <- c("Input", "Main effect 1", "Main effect 2")
  # p
  
  p <- rbind(p1, p2)
  p$Output <- c(rep("1", main_effects_num_vars_1),
                rep("2", main_effects_num_vars_2))
  p$Colour <- rep(NA, nrow(p))
  for(i in 1:length(input_names)){
    p$Colour[p$Input == input_names[i]] <- c25[i]
  }
  
  pie(as.numeric(p$`Main effect`[1:main_effects_num_vars_1]),
      p$Input[1:main_effects_num_vars_1],
      col = p$Colour[1:main_effects_num_vars_1],
      clockwise = T, init.angle = 180)
  pie(as.numeric(p$`Main effect`[(main_effects_num_vars_1+1):(main_effects_num_vars_1+main_effects_num_vars_2)]),
      p$Input[(main_effects_num_vars_1+1):(main_effects_num_vars_1+main_effects_num_vars_2)],
      col = p$Colour[(main_effects_num_vars_1+1):(main_effects_num_vars_1+main_effects_num_vars_2)],
      clockwise = T, init.angle = 180)
}

# SA union statement ####

if(perform_SA == "yes"){
  for (i in 1:length(union(main_effects_vars1, main_effects_vars2))) {
    cat(i, ". ", union(main_effects_vars1, main_effects_vars2)[i], "\n", sep="")
  }
}

### Emulator fit and validation ################################################

# estimation-method ####

estimation_method <- ifelse(em == "post_mode", "maximum posterior mode estimation", "maximum likelihood estimation")

# remove objects ####

if(perform_SA == "yes"){
  rm(main_effects_num_vars_1, main_effects_num_vars_2, main_effects_sorted1, main_effects_sorted2, main_effects_sum1, main_effects_sum2, main_effects_vars1, main_effects_vars2, main_effects1, main_effects2)
}

# emulator fit ####

if(perform_SA == "yes"){
  emulator1 <- rgasp(design = input, response = output1, trend = H, lower_bound=FALSE, zero.mean = "No",
                     kernel_type = 'pow_exp', alpha = rep(2, ncol(input)),        # squared exponential
                     method = em,
                     nugget.est=TRUE
  )
  emulator2 <- rgasp(design = input, response = output2, trend = H, lower_bound=FALSE, zero.mean = "No",
                     kernel_type = 'pow_exp', alpha = rep(2, ncol(input)),        # squared exponential
                     method = em,
                     nugget.est=TRUE
  )
}
if(perform_SA == "no"){
  # start.time <- Sys.time()
  emulator1 <- rgasp(design = input, response = output1, trend = H, lower_bound=FALSE, zero.mean = "No",
                     kernel_type = 'pow_exp', alpha = rep(2, ncol(input)),        # squared exponential
                     method = em,
                     nugget = nugget_value,
                     nugget.est=FALSE
  )
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # time.taken
  # start.time <- Sys.time()
  emulator2 <- rgasp(design = input, response = output2, trend = H, lower_bound=FALSE, zero.mean = "No",
                     kernel_type = 'pow_exp', alpha = rep(2, ncol(input)),        # squared exponential
                     method = em,
                     nugget = nugget_value,
                     nugget.est=FALSE
  )
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # time.taken
}

if(perform_SA == "yes"){
  estimated_parameters_with_nugget(emulator1, 1)
  estimated_parameters_with_nugget(emulator2, 2)
  estimated_parameters1 <- c(estimated_ranges1, estimated_trends1, estimated_var1, estimated_nugget1)
  estimated_parameters2 <- c(estimated_ranges2, estimated_trends2, estimated_var2, estimated_nugget2)
}
if(perform_SA == "no"){
  estimated_parameters(emulator1, 1)
  estimated_parameters(emulator2, 2)
  estimated_parameters1 <- c(estimated_ranges1, estimated_trends1, estimated_var1, nugget_value)
  estimated_parameters2 <- c(estimated_ranges2, estimated_trends2, estimated_var2, nugget_value)
}
rm(nugget_estimation)

# parameter estimates ####

options(digits = 1, scipen = 0)
kable(
  cbind(
    c(1:(dim_inputs+1), 1:dim_inputs, 1, 1),
    round(estimated_parameters1),
    round(estimated_parameters2)
  ),
  align = c("c", "c", "c"),
  caption = "Estimated parameters",
  col.names = c("Ranges, trends, variance and nugget", "Output 1", "Output 2")
)
options(scipen = 999)
if(perform_SA == "yes"){
  estimated_nugget1
  estimated_nugget2
}

# posterior predictions ####

predictions1 <- RobustGaSP::predict(emulator1, testing_input, testing_trend = H_testing)
predictions2 <- RobustGaSP::predict(emulator2, testing_input, testing_trend = H_testing)

# A, A_inv and t ####

if(perform_SA == "yes"){
  A1 <- A_maker(input, estimated_ranges1) + diag(estimated_nugget1, num_obs)
  A2 <- A_maker(input, estimated_ranges2) + diag(estimated_nugget2, num_obs)
}
if(perform_SA == "no"){
  A1 <- A_maker(input, estimated_ranges1) + diag(nugget_value, num_obs)
  A2 <- A_maker(input, estimated_ranges2) + diag(nugget_value, num_obs)
}
A_inv1 <- solve(A1)
A_inv2 <- solve(A2)

t1 <- t_maker(input, testing_input, estimated_ranges1)
t2 <- t_maker(input, testing_input, estimated_ranges2)

# loo CV ####

loo_cv(emulator1, output1, 1)
loo_cv(emulator2, output2, 2)

# validation checks ####

kable(cbind(
  c("1", "2"),
  c(round(mean(predictions1$upper95 - predictions1$lower95),2),
    round(mean(predictions2$upper95 - predictions2$lower95),2)),
  c(round(mean(loo_cv1$upper95 - loo_cv1$lower95),2),
    round(mean(loo_cv2$upper95 - loo_cv2$lower95),2)),
  c(round(100 * sum(loo_cv1$validated) / nrow(loo_cv1)),
    round(100 * sum(loo_cv2$validated) / nrow(loo_cv1)))
), 
      align=c("l","c","c","c"),
      caption = "Validation checks",
      col.names = c("Output", "Prediction 95% CI mean width", "LOO CV 95% CI mean width", "% validated in LOO CV")
      )

# loo CV plots ####

par(mar = c(4, 4, .1, .1))

line_type <- rep(2, 221)
line_type[prim_so4_diam_low_indices] <- 1

num_obs <- nrow(loo_cv1)
cols <- rep("red", num_obs)
cols[loo_cv1$validated] <- "black"

pdf('plots/loo_cv1.pdf')
plot(loo_cv1$loo_prediction ~ loo_cv1$actual,
     xlim = c(min(loo_cv1$lower95, loo_cv1$actual), 
              max(loo_cv1$upper95, loo_cv1$actual)),
     ylim = c(min(loo_cv1$lower95, loo_cv1$actual), 
              max(loo_cv1$upper95, loo_cv1$actual)),
     col = cols, pch = 20,
     xlab = "Model", ylab = "Emulator")
abline(0,1, col="blue")
segments(loo_cv1$actual,
         loo_cv1$lower95,
         loo_cv1$actual,
         loo_cv1$upper95,
         col = cols,
         lty = line_type)

num_obs <- nrow(loo_cv2)
cols <- rep("red", num_obs)
cols[loo_cv2$validated] <- "black"

pdf('plots/loo_cv2.pdf')
plot(loo_cv2$loo_prediction ~ loo_cv2$actual,
     xlim = c(min(loo_cv2$lower95, loo_cv2$actual), 
              max(loo_cv2$upper95, loo_cv2$actual)),
     ylim = c(min(loo_cv2$lower95, loo_cv2$actual), 
              max(loo_cv2$upper95, loo_cv2$actual)),
     col = cols, pch = 20,
     xlab = "Model", ylab = "Emulator")
abline(0,1, col="blue")
segments(loo_cv2$actual,
         loo_cv2$lower95,
         loo_cv2$actual,
         loo_cv2$upper95,
         col = cols,
         lty = line_type)

### Check of by-hand calculations ##############################################

# by-hand checks ####

par(mar = c(4, 4, .1, .1))

byhand_predictions1 <- byhand_preds(testing_input,
                                   estimated_trends1,
                                   t1,
                                   A_inv1,
                                   output1,
                                   H)
byhand_predictions2 <- byhand_preds(testing_input,
                                   estimated_trends2,
                                   t2,
                                   A_inv2,
                                   output2,
                                   H)

byhand_score1 <- sum(abs(predictions1$mean - as.vector(byhand_predictions1))) / num_testing_input
# / (10^(-8)*num_testing_input)
byhand_score2 <- sum(abs(predictions2$mean - as.vector(byhand_predictions2))) / num_testing_input
# / (10^(-8)*num_testing_input)

byhand_score1
byhand_score2

byhand_plots(input, output, testing_input, predictions1, byhand_predictions1)
byhand_plots(input, output, testing_input, predictions2, byhand_predictions2)

### Calculation of the alignment measure #######################################

# partial derivatives ####

pds1 <- data.frame(matrix(NA, nrow = num_testing_input, ncol = dim_inputs))
for (r in 1:dim_inputs) {
  pds1[,r] <-
    matrix(estimated_trends1[r+1], nrow = num_testing_input) + 
    ((2*(matrix(rep(input[,r], num_testing_input), nrow = num_testing_input, byrow = T) - matrix(rep(testing_input[,r], num_obs), nrow = num_testing_input, byrow = F)) / (estimated_ranges1[r]^2)) *
       t1) %*% 
    A_inv1 %*% (as.matrix(output1, ncol = n) - H %*% estimated_trends1)
  # cat("pd_t:", r, "/", dim_inputs, ",", sep = "")
}
pds_normed1 <- pds1 / sqrt(rowSums(pds1^2))

pds2 <- data.frame(matrix(NA, nrow = num_testing_input, ncol = dim_inputs))
for (r in 1:dim_inputs) {
  pds2[,r] <-
    matrix(estimated_trends2[r+1], nrow = num_testing_input) + 
    ((2*(matrix(rep(input[,r], num_testing_input), nrow = num_testing_input, byrow = T) - matrix(rep(testing_input[,r], num_obs), nrow = num_testing_input, byrow = F)) / (estimated_ranges2[r]^2)) *
       t2) %*% 
    A_inv2 %*% (as.matrix(output2, ncol = n) - H %*% estimated_trends2)
  # cat("pd_t:", r, "/", dim_inputs, ",", sep = "")
}
pds_normed2 <- pds2 / sqrt(rowSums(pds2^2))

# alignment measure ####

AM <- mean(abs(rowSums(pds_normed1 * pds_normed2)))

plot(output2 ~ output1)
abline(0,1, col="blue")