# Code for "[Gaussian Process Nowcasting: Application to COVID-19 Mortality Reporting](https://arxiv.org/abs/2102.11249)"

https://github.com/ihawryluk/GP_nowcasting



#### Abstract

Updating observations of a signal due to the delays in the measurement process is a common problem in signal processing, with prominent examples in a wide range of fields. An important example of this problem is the nowcasting of COVID-19 mortality: given a stream of reported counts of daily deaths, can we correct for the delays in reporting to paint an accurate picture of the present, with uncertainty? Without this correction, raw data will often mislead by suggesting an improving situation. We present a flexible approach using a latent Gaussian process that is capable of describing the changing auto-correlation structure present in the reporting time-delay surface. This approach also yields robust estimates of uncertainty for the estimated nowcasted numbers of deaths. We test assumptions in model specification such as the choice of kernel or hyper priors, and evaluate model performance on a challenging real dataset from Brazil. Our experiments show that Gaussian process nowcasting performs favourably against both comparable methods, and a small sample of expert human predictions. Our approach has substantial practical utility in disease modelling --- by applying our approach in COVID-19 mortality data from Brazil, where reporting delays are large, we can make informative predictions on important epidemiological quantities such as the current effective reproduction number.

### Contents

#### data

All data needed to reproduce the analysis is contained in the *data* folder. 

- *df_SIVEP_nowcast_Brazil_31-05-2021, df_SIVEP_nowcast_allStates_08-02-2021.csv* -- the SIVEP-Gripe datasets up to the 31st May 2021 and 8th Feb 2021 releases, which were downloaded from <https://opendatasus.saude.gov.br/dataset/bd-srag-2020> website and further processed
- *experts_response.csv* -- nowcasts of the group of anonymous infectious epidemiology experts, as described in section 4.3 and presented in Fig.6
- *get_states_data.R* -- script for downloading the SIVEP-Gripe datasets and pre-processing them to get *df_SIVEP_nowcast_allStates_08-02-2021.csv*

#### tests_results 

Contains compact results of nowcasts done in the retrospective tests for 9 models. Each file corresponds to a single model and a single date and contains a table with the epidemiological week, mean and credible intervals. 

#### stan_models

Contains stan models for each of the models used in retrospective tests and described in the manuscript. The GP models are described in the manuscript in detail, and the NobBS model is based on the  model from https://doi.org/10.1371/journal.pcbi.1007735

#### sensitivity

Contains scripts and stan models for the sensitivity testing described in the manuscript. 

#### example_Rt

Contains scripts, stan models and data required to generate R_t estimates for raw data, ground truth and nowcasted data, as shown in Figure 1. The nowcasted data was generated by *Daily_nowcasts_for_Rt_example_fig1.ipynb*. 

#### other files

The remaining files were used to carry out model runs, plotting and analysis of the results. Specifically *nowcasting_functions.py* contains most of the functions required to pre-process the data and carry out nowcasting, using a specified model.




### Requirements:

- Python, R, Jupyter notebooks

- Python packages: pandas, matplotlib, numpy, pystan, seaborn, arviz, sklearn, scipy, pickle, datetime

- R packages: rstan, matrixStats, data, lubridate, gdata, dplyr, tidyr, EnvStats, scales, tidyverse, abind, xtable, ggplot2, gridExtra, ggpubr, bayesplot, cowplot, openxlsx, loo, patchwork, properscoring

