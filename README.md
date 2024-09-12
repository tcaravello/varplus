# varplus
 Replication Files for "Evaluating Policy Counterfactuals: A VAR-Plus Approach" by Caravello, McKay & Wolf

Tested in Matlab R2022b on a Dell Inspiron 15.
# Important
In order to produce most Figures and Tables, the posterior mode and posterior draws for policy causal effect matrices are required.
Due to its size (~7.5 GB), these are not included in this repository. They must be downloaded from [here](https://www.dropbox.com/scl/fo/zi78j833q1py3w312dm0x/ADee6QczKcb7dqTgQXxRTwY?rlkey=ryb5ilyx0ywjftk2j1lmsz882&e=1&dl=0). Then, the /suff_stats folder in this repository must be replaced with the one downloaded above. 

To ensure that all codes run, the variable "path"---located near the top of the various m-files---needs to be changed to reflect the local machine.

# Contents
[readme.pdf](https://github.com/tcaravello/varplus/blob/main/varplus_readme.pdf): contains detailed instructions for using this repository. The main points are summarized below.

[main_varplus.m](https://github.com/tcaravello/varplus/blob/main/main_varplus.m): produces all figures and tables in the paper.

[applications](https://github.com/tcaravello/varplus/tree/main/applications): creates Figures 5-9, D.1-D.4
* [second_moments](https://github.com/tcaravello/varplus/tree/main/applications/second_moments) creates Figures 5, 6, D.1 and D.2
* [hist_evol](https://github.com/tcaravello/varplus/tree/main/applications/hist_evol) creates Figures 7 and D.3
* [hist_scenario](https://github.com/tcaravello/varplus/tree/main/applications/hist_scenario) creates Figures 8, 9 and D.4

[invertibility](https://github.com/tcaravello/varplus/tree/main/invertibility) : Functions used to create Figure 1 of the paper.

[model_estim](https://github.com/tcaravello/varplus/tree/main/model_estim) : obtains the posterior for the causal effect matrices for all models, as well as the posterior model probabilities.
* [get_param_post.m](https://github.com/tcaravello/varplus/tree/main/model_estim/get_param_post.m): obtain descriptive statistics for posterior of strucutral parameters.
* [get_posterior_mode.m](https://github.com/tcaravello/varplus/tree/main/model_estim/get_posterior_mode.m): obtains the posterior mode.
* [get_posterior_probs.m](https://github.com/tcaravello/varplus/tree/main/model_estim/get_posterior_probs.m): get posterior model probabilities, and create posterior draws across different subsets of models.
* [plot_model_irfs.m](https://github.com/tcaravello/varplus/tree/main/model_estim/plot_model_irfs.m): plot fitted impulse responses, as well as impulse responses for different models and horizons.
* [sample_posterior.m](https://github.com/tcaravello/varplus/tree/main/model_estim/sample_posterior.m): obtain posterior draws for each model.

[var_inputs](https://github.com/tcaravello/varplus/tree/main/var_inputs) : estimate 1) Wold Decomposition, 2) IRF to monetary policy shock using VARs

* [run_var_fcst_evol.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_fcst_evol.m) creates VAR forecasts for each date, used when constructing the historical evolution.
* [run_var_fcst_scenario.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_fcst_scenario.m) creates VAR forecasts at a given point in time, used when constructing the historical scenario.
* [run_var_mbc.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_mbc.m) estimates the impulse response to the Main Business Cycle shock by Angeletos, Collard & Dellas (2020, AER)
* [run_var_mp_ad.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_mp_ad.m) estimates the impulse response to a monetary policy shock, using the shock series constructed by Aruoba & Drechsel (2022)
* [run_var_spf_fcst_compare.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_spf_fcst_compare.m) compares the forecasting performance of VARs and the Survey of Profesional Forecasters.
* [run_var_swfactors.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_swfactors.m) compares the forecasting performance of VARs with and without the Stock & Watson (2016) factors.
* [run_var_wold.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_wold.m) estimates Wold IRFs for the sample 1960Q1 - 2019Q4
* [run_var_wold_early.m](https://github.com/tcaravello/varplus/tree/main/var_inputs/run_var_wold_early.m) estimates Wold IRFs for the sample 1960Q1 - 2007Q1

[_auxiliary_functions](https://github.com/tcaravello/varplus/tree/main/_auxiliary_functions): auxiliary functions and routines.

# Acknowledgements
* Our solution of the HANK model uses: first, several files from the replication codes of the article Ahn et al. (2017); second, the CompEcon toolbox of Miranda and Fackler, available here: www4.ncsu.edu/~pfackler/compecon; and third, the ergodicdist.m
function, written by Marco Maffezzoli. The codes closely build on those used by one of the authors in Wolf (2024).
* The files jbfill.m and winsorize.m are available on Mathworks file exchange.
* The file regcyc.m is taken from the replication codes for the article Hamilton (2018).


