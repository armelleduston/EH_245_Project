Files for simulations:
- run_sim.R runs a single iteration of simulation. sbatch used to run multiple iterations in parallel.
- combineResults.R merges results of simulation results.
- summarize_sims.R Takes merged simulation results and produces tables/figures

Files for HELIX analysis: 
- Make_Covariates.R preprocesses covariate data
- Make_Exposures.R preprocesses exposure data
- HELIXanalysis.R runs data analysis and produces results/figures
- HELIXanalysis_supplementary.R runs supplementary analysis with separate metals indices

