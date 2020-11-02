Estimating Uncertainty at a given spatial scale
================

Dataset needed: a grid-level dataset related to the spatial boundaries of interest. In
this case, each grid cell is associated with a given planning unit ID.
You will also need the AGLB prediction and mean/sd of error associated
with the AGLB size class of that grid cell:

<table class="table" style="width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;">

geotag

</th>

<th style="text-align:left;">

plunitid

</th>

<th style="text-align:right;">

aglb

</th>

<th style="text-align:left;">

SizeClass

</th>

<th style="text-align:right;">

sc\_mean\_err

</th>

<th style="text-align:right;">

sc\_sd\_err

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

W122687810N41639962

</td>

<td style="text-align:left;">

{dd4f0bb2-fa42-4d1a-b809-2a208fb62101}

</td>

<td style="text-align:right;">

79.43349

</td>

<td style="text-align:left;">

(75,150\]

</td>

<td style="text-align:right;">

\-8.709237

</td>

<td style="text-align:right;">

65.60805

</td>

</tr>

<tr>

<td style="text-align:left;">

W122687811N41639782

</td>

<td style="text-align:left;">

{dd4f0bb2-fa42-4d1a-b809-2a208fb62101}

</td>

<td style="text-align:right;">

101.20421

</td>

<td style="text-align:left;">

(75,150\]

</td>

<td style="text-align:right;">

\-8.709237

</td>

<td style="text-align:right;">

65.60805

</td>

</tr>

<tr>

<td style="text-align:left;">

W122687811N41639872

</td>

<td style="text-align:left;">

{dd4f0bb2-fa42-4d1a-b809-2a208fb62101}

</td>

<td style="text-align:right;">

88.30994

</td>

<td style="text-align:left;">

(75,150\]

</td>

<td style="text-align:right;">

\-8.709237

</td>

<td style="text-align:right;">

65.60805

</td>

</tr>

<tr>

<td style="text-align:left;">

W122687812N41639602

</td>

<td style="text-align:left;">

{dd4f0bb2-fa42-4d1a-b809-2a208fb62101}

</td>

<td style="text-align:right;">

120.73101

</td>

<td style="text-align:left;">

(75,150\]

</td>

<td style="text-align:right;">

\-8.709237

</td>

<td style="text-align:right;">

65.60805

</td>

</tr>

<tr>

<td style="text-align:left;">

W122687812N41639692

</td>

<td style="text-align:left;">

{dd4f0bb2-fa42-4d1a-b809-2a208fb62101}

</td>

<td style="text-align:right;">

91.80973

</td>

<td style="text-align:left;">

(75,150\]

</td>

<td style="text-align:right;">

\-8.709237

</td>

<td style="text-align:right;">

65.60805

</td>

</tr>

<tr>

<td style="text-align:left;">

W122687813N41639422

</td>

<td style="text-align:left;">

{dd4f0bb2-fa42-4d1a-b809-2a208fb62101}

</td>

<td style="text-align:right;">

80.10650

</td>

<td style="text-align:left;">

(75,150\]

</td>

<td style="text-align:right;">

\-8.709237

</td>

<td style="text-align:right;">

65.60805

</td>

</tr>

</tbody>

</table>

Given the above dataset, run the following boostrapping function across
all planning units of interest:

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
  unc.predicts <- bootstrap.dt[, (names(bootstrap.dt)) := lapply(.SD, function(x){mean(x)}), .SDcols = names(bootstrap.dt)][1,] %>%
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

<table class="table" style="width: auto !important; ">

<thead>

<tr>

<th style="text-align:left;">

plunitid

</th>

<th style="text-align:right;">

unc\_sp\_mean

</th>

<th style="text-align:right;">

unc\_sp\_sd

</th>

<th style="text-align:right;">

unc\_sp\_stderr

</th>

<th style="text-align:right;">

unc\_sp\_smperr

</th>

<th style="text-align:right;">

unc\_mod\_mean

</th>

<th style="text-align:right;">

unc\_mod\_sd

</th>

<th style="text-align:right;">

unc\_mod\_stderr

</th>

<th style="text-align:right;">

unc\_mod\_smperr

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

{dd4f0bb2-fa42-4d1a-b809-2a208fb62101}

</td>

<td style="text-align:right;">

99.41618

</td>

<td style="text-align:right;">

37.47442

</td>

<td style="text-align:right;">

1.334971

</td>

<td style="text-align:right;">

2.208923

</td>

<td style="text-align:right;">

94.05518

</td>

<td style="text-align:right;">

2.093305

</td>

<td style="text-align:right;">

0.2093305

</td>

<td style="text-align:right;">

0.3661135

</td>

</tr>

</tbody>

</table>
