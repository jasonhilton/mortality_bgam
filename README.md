---
title: "README"
output: html_document
---



# Projecting UK Mortality using Bayesian Generalised Additive Models

This repository contains the code need to reproduce the results in the paper "Projecting UK Mortality using Bayesian Generalised Additive Models". 

# Requirements

To run the basic model, you need R (version 3+ should be fine) and the ability to install R libraries, and `rstan` in particular. You will probably need to install Rtools, following the `rstan` [installation instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows). 

To run the model repeatedly at different cutpoints and produce a 'stacked' forecast, a High Performance Computing platform, ([`iridis4`](https://www.southampton.ac.uk/isolutions/staff/iridis.page) based at the University of Southampton) was used to run simulations in parallel and thus minimise runtime. The portable batch system (PBS) scripts needed to conduct the runs on this platform are included, and it should be possible to run these on similar systems with minimal adaptation.

To produce the final manuscript, Rmarkdown, knitr, latex and pandoc are used.

## Storage
The complete set of MCMC results takes up about 20gb.

---

# Reproducing the results.

The results presented in the paper can be reproduced by following the steps below. All code chucks in this document should be run from the terminal from the base directory of this repository (that is, the one that contains this file).

## Installing R packages
- The first step is to install a number of R packages. The command below does so programmatically, and can be run from the command line.

```bash
Rscript install_required_packages.R
```
At the time of writing, the development version of the [`loo`](https://github.com/stan-dev/loo) is also needed, as some features relating to stacking are not yet incorporated into the full release. The script above installs the latest version of `loo` from github using the devtools package.

## Getting the data
UK mortality data from the human mortality database is required. This data can also be obtained using an script. Running this script will create a `data` subfolder within this repository and download files from the [Human Mortality Database](http://www.mortality.org/) (HMD) to it. However, an account and password is needed to access this data. Included in the repository is the file `dummy_hmd_credentials.yaml`. Editing this to refer to your username and password, and resaving it as `hmd_credentials.yaml` in the same location will the data to be obtained with the following command: 


```bash
Rscript download_mortality_data.R
```
This script uses the third-party `HMDHFDplus` package to obtain the data. Please ensure you use a unique password for HMD and don't reuse one from something else, and take appropriate precautions with your credential information. I am not a security expert, and I cannot guarantee the security of these credentials. If you prefer, you can download the UK 1x1 deaths and exposure data manually from the hmd website and save them in `data/hmd/deaths.rds` and `data/hmd/exposures.rds` respectively.

## Sampling for a single model

The script `run_single_stan_model.R` produces samples for a single model. 
Exactly what data and which model is used is determined by the configuration file.


```bash
Rscript run_single_stan_model.R
```
Outputs are saved to a directory determined by the config file - by default this is in the `results`, in subfolders relating to the model type and sex.
The `config` folder contains configurations needed for the analyses in the paper.
The models themselves are located in the `stan` folder. 

A configuration file and sex can be optionally specified at the command line as follows:

```bash
Rscript run_single_stan_model.R config/single_sex.yaml male
Rscript run_single_stan_model.R config/two_sex.yaml
```


## Running separate model for each cutpoint and stacking the results

_Requires HPC environment_

Start by cloning or copying this repository into your working directory on your HPC resource.

A bash script is provided that submits the relevant jobs to the queue management system via the `qsub` command. A config file and sex must be provided as arguments as shown. This scripts submits one job for each cutpoint, each requiring 4 processors (1 for each MCMC chain). The script also specifies that an additional R script `run_stacking.R` should be run when all jobs are finished. This script calculates the model weights and produces the weighted forecasts. All results are saved to the `results` folder by default. The actual batch scripts used can be found in the `pbs` folder.


```bash
./send_jobs.sh config/single_sex.yaml male
./send_jobs.sh config/single_sex.yaml female
./send_jobs.sh config/two_sex.yaml
```

## Running stacked models on held-back data subset
For informal model assessment, the process is repeated on a subset of the data excluding the latest ten years. A similar command as before is needed here.


```bash
./send_jobs_two_sex.sh config/single_sex_hb.yaml male
./send_jobs_two_sex.sh config/single_sex_hb.yaml female
./send_jobs_two_sex.sh config/two_sex_hb.yaml
```

## Create plot data
This script processes the raw model objects and produces the data required to create the plots in the paper. It can be submitted as a batch HPC job as follows:


```bash
qsub PBS/create_plot_data.pbs
```

## Compile the paper
To produce the `tex` files needed to produce the manuscript, the below commands run knitr and pandoc to process the R markdown files `gam_mortality.Rmd` and `Appendix.Rmd`.


```bash
cd paper
Rscript process_md.R
```
The resulting tex file can be `\include`-d in a tex template and compiled as normal.


## R session information

### Local Machine

```r
library(dplyr)
library(tidyr)
library(lavyeval)
```

```
## Error in library(lavyeval): there is no package called 'lavyeval'
```

```r
library(tidyselect)
```

```
## Warning: package 'tidyselect' was built under R version 3.4.4
```

```r
library(magrittr)
library(rstan)
library(purrr)
library(ggfan)
library(HMDHFDplus)
```

```
## Warning: package 'HMDHFDplus' was built under R version 3.4.4
```

```r
library(boot)
library(abind)
library(tibble)
library(loo)
```

```
## This is loo version 2.0.0.
## **NOTE: As of version 2.0.0 loo defaults to 1 core but we recommend using as many as possible. Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session. Visit mc-stan.org/loo/news for details on other changes.
```

```r
library(XLConnect)
```

```
## Warning: package 'XLConnect' was built under R version 3.4.4
```

```
## Loading required package: XLConnectJars
```

```
## Warning: package 'XLConnectJars' was built under R version 3.4.4
```

```
## XLConnect 0.2-15 by Mirai Solutions GmbH [aut],
##   Martin Studer [cre],
##   The Apache Software Foundation [ctb, cph] (Apache POI),
##   Graph Builder [ctb, cph] (Curvesapi Java library)
```

```
## http://www.mirai-solutions.com
## https://github.com/miraisolutions/xlconnect
```

```r
library(curl)
```

```
## Warning: package 'curl' was built under R version 3.4.4
```

```r
library(yaml)
```

```
## Warning: package 'yaml' was built under R version 3.4.4
```

```r
library(devtools)
```

```
## Warning: package 'devtools' was built under R version 3.4.3
```

```r
devtools::session_info()
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.2 (2017-09-28)
##  system   x86_64, mingw32             
##  ui       RTerm                       
##  language (EN)                        
##  collate  English_United Kingdom.1252 
##  tz       Europe/London               
##  date     2018-06-14
```

```
## Packages -----------------------------------------------------------------
```

```
##  package       * version    date       source                            
##  abind         * 1.4-5      2016-07-21 CRAN (R 3.3.1)                    
##  assertthat      0.2.0      2017-04-11 CRAN (R 3.3.3)                    
##  backports       1.1.2      2017-12-13 CRAN (R 3.4.3)                    
##  base          * 3.4.2      2017-09-28 local                             
##  bindr           0.1.1      2018-03-13 CRAN (R 3.4.4)                    
##  bindrcpp      * 0.2.2      2018-03-29 CRAN (R 3.4.4)                    
##  bitops          1.0-6      2013-08-17 CRAN (R 3.2.0)                    
##  boot          * 1.3-20     2017-08-06 CRAN (R 3.4.2)                    
##  colorspace      1.4-0      2017-12-23 R-Forge (R 3.4.3)                 
##  compiler        3.4.2      2017-09-28 local                             
##  curl          * 3.2        2018-03-28 CRAN (R 3.4.4)                    
##  datasets      * 3.4.2      2017-09-28 local                             
##  devtools      * 1.13.5     2018-02-18 CRAN (R 3.4.3)                    
##  digest          0.6.15     2018-01-28 CRAN (R 3.4.3)                    
##  dplyr         * 0.7.5      2018-05-19 CRAN (R 3.4.4)                    
##  evaluate        0.10.1     2017-06-24 CRAN (R 3.4.1)                    
##  ggfan         * 0.1.2      2018-06-14 CRAN (R 3.4.4)                    
##  ggplot2       * 2.2.1.9000 2018-06-13 Github (tidyverse/ggplot2@4db5122)
##  glue            1.2.0      2017-10-29 CRAN (R 3.4.2)                    
##  graphics      * 3.4.2      2017-09-28 local                             
##  grDevices     * 3.4.2      2017-09-28 local                             
##  grid          * 3.4.2      2017-09-28 local                             
##  gridExtra       2.3        2017-09-09 CRAN (R 3.4.2)                    
##  gtable          0.2.0      2016-02-26 CRAN (R 3.2.4)                    
##  highr           0.7        2018-06-09 CRAN (R 3.4.4)                    
##  HMDHFDplus    * 1.1.8      2018-05-22 local                             
##  htmltools       0.3.6      2017-04-28 CRAN (R 3.4.0)                    
##  inline          0.3.15     2018-05-18 CRAN (R 3.4.4)                    
##  knitr         * 1.20       2018-02-20 CRAN (R 3.4.4)                    
##  labeling        0.3        2014-08-23 CRAN (R 3.2.0)                    
##  lazyeval        0.2.1      2017-10-29 CRAN (R 3.4.2)                    
##  loo           * 2.0.0      2018-04-11 CRAN (R 3.4.2)                    
##  magrittr      * 1.5        2014-11-22 CRAN (R 3.2.0)                    
##  matrixStats     0.53.1     2018-02-11 CRAN (R 3.4.4)                    
##  memoise         1.1.0      2017-04-21 CRAN (R 3.4.0)                    
##  methods         3.4.2      2017-09-28 local                             
##  munsell         0.5.0      2018-06-12 CRAN (R 3.4.4)                    
##  parallel        3.4.2      2017-09-28 local                             
##  pillar          1.2.3      2018-05-25 CRAN (R 3.4.4)                    
##  pkgconfig       2.0.1      2017-03-21 CRAN (R 3.4.1)                    
##  plyr            1.8.4      2016-06-08 CRAN (R 3.3.1)                    
##  purrr         * 0.2.5      2018-05-29 CRAN (R 3.4.4)                    
##  R6              2.2.2      2017-06-17 CRAN (R 3.4.1)                    
##  Rcpp            0.12.17    2018-05-18 CRAN (R 3.4.4)                    
##  RCurl           1.95-4.10  2018-01-04 CRAN (R 3.4.3)                    
##  rJava           0.9-10     2018-05-29 CRAN (R 3.4.4)                    
##  rlang           0.2.1      2018-05-30 CRAN (R 3.4.4)                    
##  rmarkdown     * 1.10       2018-06-11 CRAN (R 3.4.4)                    
##  rprojroot       1.3-2      2018-01-03 CRAN (R 3.4.3)                    
##  rstan         * 2.17.3     2018-01-20 CRAN (R 3.4.4)                    
##  scales        * 0.5.0.9000 2018-06-13 Github (hadley/scales@c28066a)    
##  StanHeaders   * 2.17.2     2018-01-20 CRAN (R 3.4.4)                    
##  stats         * 3.4.2      2017-09-28 local                             
##  stats4          3.4.2      2017-09-28 local                             
##  stringi         1.2.3      2018-06-12 CRAN (R 3.4.4)                    
##  stringr         1.3.1      2018-05-10 CRAN (R 3.4.4)                    
##  tibble        * 1.4.2      2018-01-22 CRAN (R 3.4.4)                    
##  tidyr         * 0.8.1      2018-05-18 CRAN (R 3.4.4)                    
##  tidyselect    * 0.2.4      2018-02-26 CRAN (R 3.4.4)                    
##  tools           3.4.2      2017-09-28 local                             
##  utils         * 3.4.2      2017-09-28 local                             
##  withr           2.1.2      2018-06-13 Github (jimhester/withr@dbcd7cd)  
##  XLConnect     * 0.2-15     2018-04-05 CRAN (R 3.4.4)                    
##  XLConnectJars * 0.2-15     2018-04-05 CRAN (R 3.4.4)                    
##  XML             3.98-1.11  2018-04-16 CRAN (R 3.4.4)                    
##  yaml          * 2.1.19     2018-05-01 CRAN (R 3.4.4)
```


----

### Iridis


```r
print(readRDS(file.path("results", "processed", "iridis_sesh.rds")))
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.2 (2017-09-28)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  C                           
##  tz       <NA>                        
##  date     2018-06-14
```

```
## Packages -----------------------------------------------------------------
```

```
##  package     * version date       source        
##  abind       * 1.4-5   2016-07-21 CRAN (R 3.4.2)
##  assertthat    0.2.0   2017-04-11 CRAN (R 3.3.2)
##  base        * 3.4.2   2017-10-09 local         
##  bindr         0.1     2016-11-13 cran (@0.1)   
##  bindrcpp    * 0.2     2017-06-17 cran (@0.2)   
##  colorspace    1.3-2   2016-12-14 CRAN (R 3.4.2)
##  compiler      3.4.2   2017-10-09 local         
##  datasets    * 3.4.2   2017-10-09 local         
##  devtools    * 1.13.1  2017-05-13 CRAN (R 3.3.2)
##  digest        0.6.12  2017-01-27 CRAN (R 3.3.2)
##  dplyr       * 0.7.4   2017-09-28 CRAN (R 3.4.2)
##  ggplot2     * 2.2.1   2016-12-30 CRAN (R 3.3.2)
##  glue          1.1.1   2017-06-21 cran (@1.1.1) 
##  graphics    * 3.4.2   2017-10-09 local         
##  grDevices   * 3.4.2   2017-10-09 local         
##  grid          3.4.2   2017-10-09 local         
##  gridExtra     2.2.1   2016-02-29 CRAN (R 3.3.2)
##  gtable        0.2.0   2016-02-26 CRAN (R 3.3.2)
##  inline        0.3.14  2015-04-13 CRAN (R 3.3.2)
##  lazyeval      0.2.1   2017-10-29 CRAN (R 3.4.2)
##  loo         * 2.0.0   2018-04-11 CRAN (R 3.4.2)
##  magrittr    * 1.5     2014-11-22 CRAN (R 3.4.2)
##  matrixStats   0.52.2  2017-04-14 cran (@0.52.2)
##  memoise       1.1.0   2017-04-21 CRAN (R 3.3.2)
##  methods       3.4.2   2017-10-09 local         
##  munsell       0.4.3   2016-02-13 CRAN (R 3.4.2)
##  parallel      3.4.2   2017-10-09 local         
##  pillar        1.1.0   2018-01-14 CRAN (R 3.4.2)
##  pkgconfig     2.0.1   2017-03-21 CRAN (R 3.4.2)
##  plyr          1.8.4   2016-06-08 CRAN (R 3.3.2)
##  purrr       * 0.2.4   2017-10-18 CRAN (R 3.4.2)
##  R6            2.2.2   2017-06-17 cran (@2.2.2) 
##  Rcpp          0.12.11 2017-05-22 CRAN (R 3.3.2)
##  rlang         0.1.6   2017-12-21 CRAN (R 3.4.2)
##  rstan       * 2.17.2  2017-12-21 CRAN (R 3.4.2)
##  scales        0.4.1   2016-11-09 CRAN (R 3.3.2)
##  StanHeaders * 2.17.1  2017-12-20 CRAN (R 3.4.2)
##  stats       * 3.4.2   2017-10-09 local         
##  stats4        3.4.2   2017-10-09 local         
##  tibble      * 1.4.1   2017-12-25 CRAN (R 3.4.2)
##  tidyr       * 0.7.2   2017-10-16 CRAN (R 3.4.2)
##  tidyselect    0.2.3   2017-11-06 CRAN (R 3.4.2)
##  utils       * 3.4.2   2017-10-09 local         
##  withr         1.0.2   2016-06-20 CRAN (R 3.3.2)
##  yaml        * 2.1.16  2017-12-12 CRAN (R 3.4.2)
```
