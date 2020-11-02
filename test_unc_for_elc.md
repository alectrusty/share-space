Estimating Uncertainty at a given spatial scale
================
Alec
11/2/2020

    ##                 geotag                               plunitid      aglb
    ## 1: W122687810N41639962 {dd4f0bb2-fa42-4d1a-b809-2a208fb62101}  79.43349
    ## 2: W122687811N41639782 {dd4f0bb2-fa42-4d1a-b809-2a208fb62101} 101.20421
    ## 3: W122687811N41639872 {dd4f0bb2-fa42-4d1a-b809-2a208fb62101}  88.30994
    ## 4: W122687812N41639602 {dd4f0bb2-fa42-4d1a-b809-2a208fb62101} 120.73101
    ## 5: W122687812N41639692 {dd4f0bb2-fa42-4d1a-b809-2a208fb62101}  91.80973
    ## 6: W122687813N41639422 {dd4f0bb2-fa42-4d1a-b809-2a208fb62101}  80.10650
    ##    SizeClass sc_mean_err sc_sd_err
    ## 1:  (75,150]   -8.709237  65.60805
    ## 2:  (75,150]   -8.709237  65.60805
    ## 3:  (75,150]   -8.709237  65.60805
    ## 4:  (75,150]   -8.709237  65.60805
    ## 5:  (75,150]   -8.709237  65.60805
    ## 6:  (75,150]   -8.709237  65.60805

## Estimating Uncertainty at the Stand Level

Given a grid-level dataset with corresponding planning unit ID and error
bins:

``` r
# set up cluster for parallel computing
cl <- makeCluster(detectCores())

# iterate through all unique planning unit ID's, applying a bootstrapping function
uncertainty.dt <- (error.grid.dt$plunitid %>% unique) %>% parLapply(cl, ., function(x, error.grid.dt){
  library(data.table)
  library(magrittr)
  
  # grab all grid points for the given planning unit
  grid.dt <- error.grid.dt[plunitid %in% x,]

  # iterate through all grid cells, taking 100 bootstrap samples from their corresponding error distribution
  # adding that error (+/-) to the predicted AGLB
  bootstrap.dt <- lapply(1:NROW(grid.dt),function(x){
    grid.dt$aglb[x]+rnorm(100,grid.dt$sc_mean_err[x],grid.dt$sc_sd_err[x])
    }) %>%
    do.call(rbind,.) %>%
    data.table
  
  # matrix rotation so that we can do cumulative sum along the bootstrap columns
  # This is an area based measurement (MT AGLB /acre) -- that's why we take the mean
  unc.predicts <- bootstrap.dt[ , (names(bootstrap.dt)) := lapply(.SD, function(x){mean(x)}), .SDcols = names(bootstrap.dt)][1,] %>%
    unlist
  
  # summarize spatial and model uncertainty - mean, sd, standard error, sampling error
  data.table(
    plunitid=x,
    unc_sp_mean=mean(grid.dt$aglb),
    unc_sp_sd=sd(grid.dt$aglb),
    unc_sp_stderr=sqrt(var(grid.dt$aglb)/length(grid.dt$aglb)),
    unc_sp_smperr=(sqrt(var(grid.dt$aglb)/length(grid.dt$aglb))/mean(grid.dt$aglb))*100*1.645,
    unc_mod_mean=mean(unc.predicts),
    unc_mod_sd=sd(unc.predicts),
    unc_mod_stderr=sqrt(var(unc.predicts)/length(unc.predicts)),
    unc_mod_smperr=(sqrt(var(unc.predicts)/length(unc.predicts))/mean(unc.predicts))*100*1.645)
}, error.grid.dt) %>%
  rbindlist

stopCluster(cl)
```

Each planning unit will be populated with the following summary table:

    ##                                  plunitid unc_sp_mean unc_sp_sd
    ## 1: {dd4f0bb2-fa42-4d1a-b809-2a208fb62101}    99.41618  37.47442
    ##    unc_sp_stderr unc_sp_smperr unc_mod_mean unc_mod_sd unc_mod_stderr
    ## 1:      1.334971      2.208923     93.81618   2.142725      0.2142725
    ##    unc_mod_smperr
    ## 1:      0.3757117
