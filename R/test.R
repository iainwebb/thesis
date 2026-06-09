# --------------------------------------------------------------------------------------------------
# indices of low R2_adj gridboxes
# --------------------------------------------------------------------------------------------------
library(ncdf4)
ACURE_P3_AOD_Total_jul_nc <- nc_open(paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc"))
longitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"longitude") - 180
latitude = ncvar_get(ACURE_P3_AOD_Total_jul_nc,"latitude")
rm("ACURE_P3_AOD_Total_jul_nc")
ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix <-
  matrix(as.numeric(readLines("objects/R2_adj_ACURE_P3_AOD_Total_jul_level2_matrix.txt")),
         nrow = length(longitude))
threshold <- 0.60
ERF_PPE_global_jul_R2_adj_matrix <-
  matrix(as.numeric(readLines("objects/R2_adj_ERF_PPE_global_jul_matrix.txt")),
         nrow = length(longitude))
double_gp_indices <- which(ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix < threshold & ERF_PPE_global_jul_R2_adj_matrix < threshold, arr.ind = TRUE)
double_gp_indices_dataframe <- as.data.frame(double_gp_indices)
double_gp_indices_sorted_dataframe <- double_gp_indices_dataframe[order(double_gp_indices_dataframe[,1]),]
rownames(double_gp_indices_sorted_dataframe) <- 1:nrow(double_gp_indices_dataframe)
double_gp_indices_sorted_dataframe
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# recording ALL lm'ed partial derivatives, not just those above the threshold
# --------------------------------------------------------------------------------------------------
inputs_x_norm_T_matrix <- read.csv("data/ACURE-UKESM1-PPE-par-norm-Design.csv")
library(ncdf4)
ACURE_P3_AOD_Total_jul_nc <- nc_open(paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_jul.nc"))
ACURE_P3_AOD_Total_jul_array <- ncvar_get(ACURE_P3_AOD_Total_jul_nc, names(ACURE_P3_AOD_Total_jul_nc$var))
ACURE_P3_AOD_Total_jul_level2_array <- ACURE_P3_AOD_Total_jul_array[,,2,]
rm("ACURE_P3_AOD_Total_jul_nc", "ACURE_P3_AOD_Total_jul_array")

lm_all_d_dxi_normalised_ACURE_P3_AOD_Total_jul_level2_array <-
  array(as.numeric(readLines("objects/lm_all_d_dxi_normalised_ACURE_P3_AOD_Total_jul_level2_array.txt")),
        dim = c(length(longitude), length(latitude), ncol(inputs_x_norm_T_matrix)))

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# AM for when AOD_Total has been lm'ed for ALL gridboxes (not just those above the threshold)
# --------------------------------------------------------------------------------------------------
p <- ncol(inputs_x_norm_T_matrix)
N <- 100000
assign(paste0("gp_d_dxi_normalised_ERF_PPE_global_jul_matrix"),
       matrix(as.numeric(readLines(paste0("objects/gp_d_dxi_normalised_ERF_PPE_global_jul_matrix", "_", format(N, scientific = F), ".txt"))),
              nrow = p)
)

sum(is.na(lm_all_d_dxi_normalised_ACURE_P3_AOD_Total_jul_level2_array))

# --------------------------------------------------------------------------------------------------
# AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix <- matrix(NA, nrow = 192, ncol = 144)
# --------------------------------------------------------------------------------------------------
# CREATING
# for (long in 1:192) {
#   for (lat in 1:144) {
#     AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix[long, lat] <- 
#         sum(abs(colSums(matrix(rep(lm_all_d_dxi_normalised_ACURE_P3_AOD_Total_jul_level2_array[long, lat,], N), nrow = p, byrow = F) *
#                           gp_d_dxi_normalised_ERF_PPE_global_jul_matrix))
#         ) / N
#   }
#   cat(long, "/192, ", sep = "")
# }
# library(readr)
# write_lines(as.vector(AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix),
#             file=paste0("objects/AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix_", threshold*100,
#             "_", format(N, scientific = F), ".txt"))

# from before

assign(paste0("AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix"),
       matrix(as.numeric(readLines(paste0("objects/AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix_", threshold*100,
                                          "_", format(N, scientific = F), ".txt"))),
              nrow = length(longitude))
)

# MAP
AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix_for_plotting <- rbind(AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix[97:192,],
                                                                                        AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix[1:96,])
library(fields)
library(maps)
image.plot(longitude, latitude,
           AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix_for_plotting,
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           # col = c("grey90","grey80","grey70", "grey60", "grey50",
           #         "grey40", "grey30", "grey20", "grey10", "grey0"),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0),
             legend.lab = "AM for AOD Total at level 2 for July versus ERF for July"
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

range(AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix_for_plotting)
hist(AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix_for_plotting, xlim = c(0,1))
image.plot(longitude, latitude,
           rbind(ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix[97:192,],
                 ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix[1:96,]),
           breaks = c(0,0.6,0.7,0.75,0.8,1),
           col = c("blue","green","orange","yellow", "white"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.6,0.7,0.75,0.8,1),labels=as.character(c(0,0.6,0.7,0.75,0.8,1)),mgp=c(3,0.5,0)
           )
)
map("world",lwd=1.2,add=TRUE, lty=1)
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# Now with 113 gp'ed
# --------------------------------------------------------------------------------------------------
AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix <- 
  matrix(as.numeric(readLines("objects/AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix.txt")),
         nrow = length(longitude))
AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_2 <- matrix(NA, nrow = 192, ncol = 144)
for (long in 1:192) {
  for (lat in 1:144) {
    AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_2[long, lat] <- 
      if (ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix[long, lat] < 0.6) AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix[long, lat]
      else NA
  }
  cat(long, "/192, ", sep = "")
}
AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_2_for_plotting <- rbind(AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_2[97:192,],
                                                                                          AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_2[1:96,])
image.plot(longitude, latitude,
           AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_2_for_plotting,
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           # col = c("grey90","grey80","grey70", "grey60", "grey50",
           #         "grey40", "grey30", "grey20", "grey10", "grey0"),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0),
             legend.lab = "AM for AOD Total at level 2 for July versus ERF for July"
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_for_plotting <- rbind(AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix[97:192,],
                                                                                          AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix[1:96,])
image.plot(longitude, latitude,
           AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_for_plotting,
           breaks = c(0,0.1,0.2,0.3,0.4,0.5,
                      0.6,0.7,0.8,0.9,1),
           # col = c("grey90","grey80","grey70", "grey60", "grey50",
           #         "grey40", "grey30", "grey20", "grey10", "grey0"),
           col = c("#ffcdcd","#ffbcbb","#ffaba8", "#ff9b95", "#ff8a82",
                   "#ff796d", "#ff6556", "#ff4d3d", "#ff3324", "#ff0000"),
           xlab = "longitude", ylab = "latitude",
           axis.args=list(
             at=c(0,0.1,0.2,0.3,0.4,0.5,
                  0.6,0.7,0.8,0.9,1),labels=as.character(c(0,0.1,0.2,0.3,0.4,0.5,
                                                           0.6,0.7,0.8,0.9,1)),mgp=c(3,0.5,0),
             legend.lab = "AM for AOD Total at level 2 for July versus ERF for July"
           )
)
map("world",lwd=1.2,add=TRUE, lty=1, col = "black")

hist(AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix_for_plotting, xlim = c(0,1))

# SCATTERPLOT
blue_black_colours <- ifelse(as.vector(ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix) < 0.60, "blue", "black")
R2_adj_AM_all_lm <- cbind("AOD_Total_R2_adj" 
                   = as.vector(ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix), 
                   "AM_for_AOD_Total_versus_ERF" 
                   = as.vector(AM_ACURE_P3_AOD_Total_jul_level2_all_lm_versus_ERF_PPE_global_jul_matrix)
)
plot(R2_adj_AM_all_lm, xlim = c(0,1), ylim = c(0,1), pch=15,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours)
R2_adj_AM <- cbind("AOD_Total_R2_adj" 
                          = as.vector(ACURE_P3_AOD_Total_jul_level2_R2_adj_matrix), 
                          "AM_for_AOD_Total_versus_ERF" 
                          = as.vector(AM_ACURE_P3_AOD_Total_jul_level2_versus_ERF_PPE_global_jul_matrix)
)
plot(R2_adj_AM, xlim = c(0,1), ylim = c(0,1), pch=15,
     cex.axis = 1.5, cex.lab = 1.5, col = blue_black_colours)
