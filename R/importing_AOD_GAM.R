# Setup ----
mainDir <- "c:/path/to/main/dir"
subDir <- "outputDirectory"
Dir <- "X:/johnson_group/Aerosol-MFR/A-CURE-UKESM1-PPE-Data/Outputs/AOD_GAM/"

# Code ----
latseq <- seq(from = 34.375, by = 1.25, length.out = 31)
lonseq <- seq(from = -10.3125, by = 1.875, length.out = 26)
options(timeout = 10000)
options(timeout = max(10000, getOption("timeout")))
# library(beepr)
beep_on_error(
for (mon in tolower(month.abb)) {
  for (i in latseq) {
  if (file.exists(paste0(Dir, mon, "/lat", i))){
    setwd(file.path(paste0(Dir, mon, "/lat", i)))
  } else {
      dir.create(file.path(paste0(Dir, mon, "/lat", i)), recursive = T)
      setwd(file.path(paste0(Dir, mon, "/lat", i)))
  }
    for (j in lonseq) {
    url <- paste0(
      "https://gws-access.jasmin.ac.uk/public/acure/Iain/AOD_GAM/", mon,
      "/lat", i, "/GAM_gradient_signs_AOD_", mon, "_898610_ilat_", i, "_ilon_", j, ".dat")
    destfile <- paste0(
      "./GAM_gradient_signs_AOD_", mon, "_898610_ilat_", i, "_ilon_", j, ".dat")
    download.file(url, destfile)
    url <- paste0(
      "https://gws-access.jasmin.ac.uk/public/acure/Iain/AOD_GAM/", mon,
      "/lat", i, "/GAM_variances_AOD_", mon, "_898610_ilat_", i, "_ilon_", j, ".dat")
    destfile <- paste0(
      "./GAM_variances_AOD_", mon, "_898610_ilat_", i, "_ilon_", j, ".dat")
    download.file(url, destfile)
    }
  }
}
)
