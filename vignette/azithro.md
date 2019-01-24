Model: azithro
================

# Reference

    Zhao Q, Tensfeldt TG, Chandra R, Mould DR. Population pharmacokinetics of
    azithromycin and chloroquine in healthy adults and paediatric malaria subjects
    following oral administration of fixed-dose azithromycin and chloroquine
    combination tablets. Malar J. 2014 Jan 29;13:36. doi: 10.1186/1475-2875-13-36.
    PubMed PMID: 24472224; PubMed Central PMCID: PMC3909452.

# Example

``` r
library(depot)
library(dplyr)
```

``` r
mod <- depot("azithro", end = 168, delta = 0.1)
```

  - 500 mg x1; then 250 mg QD x 4

<!-- end list -->

``` r
e <- ev_rx("500 then 250 q24 x 4 after 24")

e
```

    . Events:
    .   time cmt amt evid ii addl
    . 1    0   1 500    1  0    0
    . 2   24   1 250    1 24    3

``` r
mod %>% 
  mrgsim_e(e, add = 0.05) %>% 
  plot(CP + PER2 ~time, logy=TRUE)
```

![](/Users/kyleb/git/mrgsolve/depot/vignette/azithro_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

# Code

``` c
[PROB]
1: Zhao Q, Tensfeldt TG, Chandra R, Mould DR. Population pharmacokinetics of
azithromycin and chloroquine in healthy adults and paediatric malaria subjects
following oral administration of fixed-dose azithromycin and chloroquine
combination tablets. Malar J. 2014 Jan 29;13:36. doi: 10.1186/1475-2875-13-36.
PubMed PMID: 24472224; PubMed Central PMCID: PMC3909452.

https://malariajournal.biomedcentral.com/articles/10.1186/1475-2875-13-36

[SET] end=10*24, delta=0.1

[PARAM] @annotated
TVCL :   100 : Typical value of clearance (L/h)
TVV1 :   186 : Typical central volume of distribution (L)
TVQ2 :   180 : Typical intercomp clearance 1 (L/h)
TVV2 :  2890 : Typical peripheral volume of distribution 1 (L)
Q3   :  10.6 : Intercomp clearance 2 (L/h)
V3   :  2610 : Peripheral volume of distirbution 2 (L)
KA   : 0.259 : Absorption rate constant (1/h)
WT   :    70 : Patient weight (kg)

[CMT] @annotated
GUT  : Dosing compartment (mg)
CENT : Central comaprtment (mg)
PER2 : First peripheral compartment (mg)
PER3 : Second peripheral compartment (mg)
  
[OMEGA] @annotated @block
ETACL : 0.097969 : ETA on clearance
ETAV1 : 0.280000 1.2769 : ETA on V1


[MAIN]
double CL = TVCL*pow(WT/70.0,0.75)*exp(ETACL);
double V1 = TVV1*(WT/70)*exp(ETAV1);
double Q2 = TVQ2*pow(WT/70.0,0.75);
double V2 = TVV2*(WT/70.0);

[SIGMA] @annotated
RUV : 0 : Residual unexplained variability

[ODE]
dxdt_GUT  = -KA*GUT;
dxdt_CENT =  KA*GUT - (CL+Q2+Q3)*CENT/V1 + Q2*PER2/V2 + Q3*PER3/V3;
dxdt_PER2 =  Q2*(CENT/V1 - PER2/V2);
dxdt_PER3 =  Q3*(CENT/V1 - PER3/V3);

[TABLE]
capture CP = CENT/(V1/1000.0)*exp(RUV);
```
