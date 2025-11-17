if (!require("remotes")){
  install.packages("remotes")
}

if (!require("renv")){
  remotes::install_github("rstudio/renv")
}

if (!file.exists("renv.lock")){
    renv::init()
} else {
    renv::restore(prompt = FALSE)
}

rmarkdown::render("Main_Analysis.Rmd", output_dir="../../../../html_local")