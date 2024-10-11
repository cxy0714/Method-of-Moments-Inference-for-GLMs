# Function to install and load multiple packages
install_and_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

# List of packages to install or load
packages_to_load <- c('MASS',
                      'mvtnorm',
                      'lsbclust',
                      'rje',
                      'nleqslv',
                      'rootSolve',
                      'SMUT')

# Install or load the packages
install_and_load(packages_to_load)