# Confidence intervals for the Cox model test error from cross-validation

This repo contains all the code for the method introduced in *Confidence intervals for the Cox model test error from cross-validation* by Min Woo Sun, Robert Tibshirani.
The arXiv preprint is availabe on https://arxiv.org/abs/2201.10770.

The core nested cross validation algorithm functions are in `R/ncv.R` and the base functions used for all simulation experiements are in `R/simulation.R`. You can run the function `ncv_repeated()` from  `R/ncv.R` if you want to get the mean squared error estimate of your cox model's cross-validation c-index error. 

The code for all simulations discussed in the paper is in `experiments/`. All the R scripts for the different simulations can be run on a single compute with 1 cpu, but we highly recommend using multiple cores to parallelize the repeated nested cross validation algorithm to reduce the compute time. For our work, we used computers with 32 cpu and at least 64GB ram on Stanford University's computing cluster Sherlock. All scripts used to submit sbatch jobs in Sherlock is in `experiments/sherlock`.


## Setting up the conda environment
Before running any of the scripts make sure to use `conda` to install all the necessary R packages using the `environment.yml` file by running the following commands:
```
conda env create -f environment.yml
conda activate nestedcv_r_env
```

## Unit tests

### Simulation n=100 p=10
`unittest/unittest_n100p10.R` can be used to check if the simulation function runs as expected for a simple test case scenario with 100 simulated samples and 10 features. You should expect to get lower bound miscoverage of 0.125 for both the standard cross-validation and nested cross-validation confidence interval approaches. This will result in the script printing "PASS". Do NOT run the ncv functions in parallel, i.e. keep the variable MC_CORES as 1. To run this unit test:
1. Activate the conda environment as shown above.
2. Eexcute the script by running the command `Rscript ./unittest_n100p10.R`. If you run into permission issues, make sure to make the script executable using `chmod`, e.g. `chmod u+x unittest_n100p10.R `. 
