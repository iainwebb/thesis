# we recommend running this is a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

library(cmdstanr)
library(posterior)
install.packages("bayesplot")
library(bayesplot)
color_scheme_set("brightblue")

check_cmdstan_toolchain()

install_cmdstan(cores = 2)
cmdstan_path()
