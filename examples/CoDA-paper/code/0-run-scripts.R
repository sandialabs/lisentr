# Run script to clean and prepare folders and run all CoDA paper analyses

# Clear data, figures, output folders
unlink("figs/", recursive = TRUE)
unlink("output/", recursive = TRUE)

# Create needed folders
dir.create("figs/")
dir.create("output/")

# Run all scripts
rmarkdown::render("examples/CoDA-paper/code/1-simulation_rho_phi.Rmd")
rmarkdown::render("examples/CoDA-paper/code/2-simulation_sigma.Rmd")
rmarkdown::render("examples/CoDA-paper/code/3-simulation_results.Rmd")
rmarkdown::render("examples/CoDA-paper/code/4-merra2_data_prep.Rmd")
rmarkdown::render("examples/CoDA-paper/code/5-merra2_EDA_plot.Rmd")
rmarkdown::render("examples/CoDA-paper/code/6-merra2_ESN_tuning.Rmd")
rmarkdown::render("examples/CoDA-paper/code/7-merra2_ESN_and_FI.Rmd")
rmarkdown::render("examples/CoDA-paper/code/8-merra2_results.Rmd")
rmarkdown::render("examples/CoDA-paper/code/9-merra2_rmse_plot.Rmd")
