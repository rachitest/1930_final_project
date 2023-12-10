library(devtools)
library(fs)
library(usethis)

if (.Platform$OS.type == "windows") {
  get_base_path = path_split(path_wd())[[1]][1]
  
  if (grepl("UTH/Biostats/Fall 2023", get_base_path)) {
    pkg_path = path(get_base_path, "mycglassoG2")
  } else {
    pkg_path = path(get_base_path, "UTH", "Biostats", "Fall 2023", "mycglassoG2")
  }
  
} else {
  get_base_path = path_join(path_split(path_wd())[[1]][-c(7:9)])
  
  if (grepl("UTH/Biostats/Fall 2023", get_base_path)) {
    pkg_path = path(get_base_path, "mycglassoG2")
  } else {
    pkg_path = path(get_base_path, "UTH", "Biostats", "Fall 2023", "mycglassoG2")
  }
}

create_package(pkg_path,
               fields = list(Title = "PH 1930 Final Project",
                             Version = "0.1.0",
                             Description = "This is a package that implements C++ functions for Coordinate Gradient Descent and ADMM in solving LASSO and Elastic Net Regularization problems.",
                             `Authors@R` = c(person('Rachit', 'Sabharwal', email = 'rachit.sabharwal@uth.tmc.edu', role = c('aut', 'cre', 'cph')),
                                             person('Caroline', 'Schaefer', email = 'caroline.m.schaefer@uth.tmc.edu', role = c('aut', 'cph')),
                                             person('Yu Bin', 'Chen', email = 'yu.bin.chen@uth.tmc.edu', role = c('aut', 'cph')),
                                             person('Pagna', 'Sok', email = 'pagna.sok.1@uth.tmc.edu', role = c('aut', 'cph'))),
                             License = "GPL-3"),
               open = FALSE)

# Add the C++ files to the src directory of your package
if (grepl("UTH/Biostats/Fall 2023", get_base_path)) {
  
dir_create(path(path_join(path_split(path_wd())[[1]][-c(7:9)]), "mycglassoG2", "src"))

file_copy(path(get_base_path, "1930_final_project", "coordinate_descent_algo_alt", ext = "cpp"), 
          path(path_join(path_split(path_wd())[[1]][-c(7:9)]), "mycglassoG2", "src", "coordinate_descent_algo_alt", ext = "cpp"))

file_copy(path(get_base_path, "1930_final_project", "admm_cpp_v3", ext = "cpp"), 
          path(path_join(path_split(path_wd())[[1]][-c(7:9)]), "mycglassoG2", "src", "admm_cpp_v3", ext = "cpp"))
} else {
  
dir_create(path(path_join(path_split(path_wd())[[1]][-c(7:9)]), "UTH", "Biostats", "Fall 2023", "mycglassoG2", "src"))

file_copy(path(get_base_path, "UTH", "Biostats", "Fall 2023", "1930_final_project", "coordinate_descent_algo_alt", ext = "cpp"), 
          path(path_join(path_split(path_wd())[[1]][-c(7:9)]), "UTH", "Biostats", "Fall 2023", "mycglassoG2", "src", "coordinate_descent_algo_alt", ext = "cpp"))

file_copy(path(get_base_path, "UTH", "Biostats", "Fall 2023", "1930_final_project", "admm_cpp_v3", ext = "cpp"), 
          path(path_join(path_split(path_wd())[[1]][-c(7:9)]), "UTH", "Biostats", "Fall 2023", "mycglassoG2", "src", "admm_cpp_v3", ext = "cpp"))
}

.rs.restartR()

setwd(pkg_path)

use_rcpp_armadillo("coordinate_descent_algo_alt")
use_rcpp_armadillo("admm_cpp_v3")

devtools::document(pkg_path)

check(pkg_path, args = "--as-cran")

build(pkg_path)
