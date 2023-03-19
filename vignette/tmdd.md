Model: tmdd
================

- [Reference](#reference)
- [Example](#example)
- [Code](#code)

# Reference

**Note**: this is the general, fully-parameterized TMDD model.

    Peletier LA, Gabrielsson J. Dynamics of target-mediated drug
    disposition: characteristic profiles and parameter identification. 
    J Pharmacokinet Pharmacodyn. 2012 Oct;39(5):429-51. 
    doi: 10.1007/s10928-012-9260-6. Epub 2012 Aug 1. PMID: 22851162; 
    PMCID: PMC3446204.

# Example

This simulation replicates Figure 3 in the paper.

``` r
library(depot)
library(dplyr)
```

``` r
mod <- depot("tmdd", end = 600, delta = 0.1) 
```

Doses

- 1.5 mg/kg
- 5.0 mg/kg
- 15 mg/kg
- 45 mg/kg

``` r
e <- expand.ev(amt = c(1.5, 5, 15, 45))
```

Note that these are mg/kg doses. Model parameters are also
weight-normalized.

``` r
mod %>% mrgsim_e(e) %>% plot(CP ~time, logy = TRUE)
```

![](tmdd_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Code

``` c
[ param ]

Vt = 0.1
CLd = 0.003
CL = 0.001
kon = 0.091
koff = 0.001
kin = 0.11
kout = 0.009
keRL = 0.003
R0 = 12
Vc = 0.05;

[ cmt ] L R RL T

[ pk ]

R_0 = R0;

double keL = CL/Vc;

F_L = 1/Vc;

[ des ]

dxdt_L = - kon * L * R + koff * RL - keL * L - (CLd/Vc) * L + (CLd/Vt) * T;

dxdt_R = kin - kout * R - kon * L * R + koff * RL;

dxdt_RL = kon * L * R - (koff + keRL) * RL;

dxdt_T = (CLd/Vc) * L - (CLd/Vt) * T;

[ capture ] CP = L
```
