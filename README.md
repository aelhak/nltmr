## Using linear and natural cubic splines, SITAR, and latent trajectory models to characterise nonlinear longitudinal growth trajectories in cohort studies 
Elhakeem et al.

This page contains the R analysis code used for the paper and the example analysis dataset, the files are described below:

data_setup.R: load and combine the three cohort datasets and harmonise the variable names for analysis.

polynomials_example.R: fit LME models with linear and polynomial functions and produce figure 1.

descriptive_analysis.R: run descriptive analysis, including figure 2 (scatterplots and lineplots) and supplementary tables on numbers of participants 

linear_spline_lme_models.R: fit LME models with linear spline functions, calculate model BIC, caluclate residuals for best models, calculate predicted trajectories, calculate estimated velocities (for table 3), plot trajectories (figure 3).

natural_cubic_spline_lme_models.R: fit LME models with natural cubic spline functions, calculate model BIC, caluclate residuals for best models, calculate predicted trajectories, plot estimated trajectories and velocities (figure 4), calculate mean age at peak velocity (for table 4).

sitar_models.R: fit SITAR models, calculate model BIC, caluclate variance explained, calculate predicted trajectories, plot estimated trajectories and velocities (figure 5), calculate mean age at peak velocity (for table 4).

overlayed_trajectories_plot.R: create figure 6 showing trajectories from the LME spline and SITAR models in the same plot.

growth_mixture_models.R: fit 1-5 CLclass GMM's and calculate fit statistics and plot trajectories for each model (for supplementary file), plot trajectories from the selected models (figure 7) and posterior probabilities (for supplement).

lme_velocity_curves_example.R: fit LME natural cubic spline models with increasing number of knots for the random spline and plot velocity curves (figure 8).

multicohort_sitar_models.R: fit SITAR models to all cohorts combined and plot trajectories and calcualte mean age at peak vlocity (figure 9).

synth_cohort.csv: example dataset simulated from LME natural cubic spline model predcitions in PBMAS males and females. 

getting_started_ex_dat.R: shows how to load the example dataset for practicing plots and models.
