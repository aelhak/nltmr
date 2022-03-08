## Using linear and natural cubic (regression) splines, SITAR, and latent trajectory models to characterise nonlinear longitudinal growth trajectories in cohort studies 
Ahmed Elhakeem, Rachael Hughes, Kate Tilling, Diana Cousminer, Stefan Jackowski, Tim Cole, Alex Kwong, Zheyuan Li, Struan Grant, Adam Baxter-Jones, Babette Zemel, Deborah Lawlor: BMC Med Res Methodol. 2022: doi: 10.1186/s12874-022-01542-8

This page contains the R analysis code used for the paper and the example analysis dataset, the files are described below:

**data_setup.R**: script used to load and combine the three cohort datasets and harmonise the variable names for analysis.

**polynomials_example.R**: fit LME models with linear and polynomial functions and create Fig 1.

**descriptive_analysis.R**: run descriptive analysis, including Fig 2 (scatterplots and lineplots) and supplementary tables on numbers of participants 

**linear_spline_lme_models.R**: fit LME models with linear spline functions, calculate model BIC, calculate residuals for best models, calculate predicted trajectories, calculate estimated velocities (for Table 3), plot trajectories (Fig 3).

**natural_cubic_spline_lme_models.R**: fit LME models with natural cubic spline functions, calculate model BIC, calculate residuals for best models, calculate predicted trajectories, plot estimated trajectories and velocities (Fig 4), calculate mean age at peak velocity (for Table 4).

**sitar_models.R**: fit SITAR models, calculate model BIC, calculate variance explained, calculate predicted trajectories, plot estimated trajectories and velocities (Fig 5), calculate mean age at peak velocity (for Table 4).

**overlayed_trajectories_plot.R**: create plot showing trajectories from the LME spline and SITAR models in the same plot (Fig 6).

**growth_mixture_models.R**: fit 1-5 CLclass GMM's and calculate fit statistics and plot trajectories for each model (for supplementary file), plot trajectories from the selected models (Fig 7) and posterior probabilities (for supplement).

**gamm_example.R**: fit GAMM to PBMAS cohort and plot mean trajectory and velocity curves (Fig 8)

**lme_velocity_curves_example.R**: fit LME natural cubic spline models with increasing number of knots for the random spline and plot velocity curves (Fig 9).

**multicohort_sitar_models.R**: fit SITAR models to all cohorts combined and plot trajectories and calculate mean age at peak vlocity (Fig 10).

**synth_cohort.csv**: example dataset simulated from LME natural cubic spline model predcitions in PBMAS males and females. 

**getting_started_ex_dat.R**: shows how to load the example dataset for practicing plots and models.
