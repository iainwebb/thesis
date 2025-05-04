# library(ncdf4)

month <- 3

for (i in month:month){
  print(tolower(month.abb[i]))
  # assign(paste("AODs_Total_", tolower(month.abb[1]), sep=""), array()
  file = paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_", tolower(month.abb[i]), "_PD_221.nc", sep="")
  monthly_netCDF <- nc_open(file)
  var_names <- names(monthly_netCDF$var)[-(2:4)]
  AODs_Totals_monthly <- array(0, dim(ncvar_get(monthly_netCDF, var_names[1])))
  for (j in 1:7) {
    print(paste(j, ": ", var_names[j], sep=""))
    AODs_Totals_monthly <- AODs_Totals_monthly + ncvar_get(monthly_netCDF, var_names[j])
    rm(j)
  }
  assign(paste("AODs_Total_", tolower(month.abb[i]), sep=""), AODs_Totals_monthly)
  rm(AODs_Totals_monthly, i, file, var_names)
}
dim(AODs_Total_jan)
monthly_netCDF$dim

# names(monthly_netCDF$var)
# monthly_netCDF$var$atmosphere_optical_thickness_due_to_soluble_aitken_mode_ambient_aerosol

## Check my March totals same as JDO's
JDO_mar <- nc_open(paste("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD/ACURE_P3_AOD_Total_", tolower(month.abb[month]), ".nc", sep=""))
JDO_mar <- ncvar_get(JDO_mar, names(JDO_mar$var))
sum(AODs_Total_mar - JDO_mar != 0) # 0 means all entries same, anything else indicates some difference

AODs_Total_mar_level2 <- AODs_Total_mar[,,2,]
plot(AODs_Total_mar_level2[])

# Matrix
m <- matrix(rnorm(100), ncol = 10)
colnames(m) <- paste("Col", 1:10)
rownames(m) <- paste("Row", 1:10)

heatmap(m) 

heatmap(t(AODs_Total_mar_level2[,,125]), Rowv = NA, Colv = NA)

X_norm <- read.csv("X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Design/ACURE-UKESM1-PPE-par-norm-Design.csv")
library(DiceKriging)
fit <- km(~., design = X_norm, response = AODs_Total_mar_level2[192/2, 144/2, ], covtype="gauss", optim.method="BFGS", control=list(maxit=500))

file = paste("C:/Users/smp22ijw/Downloads/ACURE_P3_AOD_Total_jul.nc")
monthly_netCDF <- nc_open(file)
