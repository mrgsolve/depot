Model: moxi
================

# Reference

    Wicha SG, Haak T, Zink K, Kees F, Kloft C, Kees MG. Population
    pharmacokinetics and target attainment analysis of moxifloxacin in patients with 
    diabetic foot infections. J Clin Pharmacol. 2015 Jun;55(6):639-46.
    doi:10.1002/jcph.464. Epub 2015 Feb 13. PubMed PMID: 25600294.

# Example

``` r
library(depot)
library(dplyr)
```

``` r
mod <- depot("moxi", end = 96, delta = 0.1) %>% zero_re()
```

  - 400 mg IV over 1 hour daily for 3 days

<!-- end list -->

``` r
e <- ev_rx("400 over 1 q 24 x 3")

e
```

    . Events:
    .   time cmt amt evid ii addl rate
    . 1    0   1 400    1 24    2  400

``` r
out <- mod %>% mrgsim_e(e)

plot(out, fDV ~ time)
```

![](/Users/kyleb/git/mrgsolve/depot/vignette/moxi_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Code

``` c
[ PROB ]
1: Wicha SG, Haak T, Zink K, Kees F, Kloft C, Kees MG. Population
pharmacokinetics and target attainment analysis of moxifloxacin in patients with 
diabetic foot infections. J Clin Pharmacol. 2015 Jun;55(6):639-46.
doi:10.1002/jcph.464. Epub 2015 Feb 13. PubMed PMID: 25600294.

- https://www.ncbi.nlm.nih.gov/pubmed/25600294

[ PARAM ] WGT = 70, IBW = 70, fu = 0.6

[ PKMODEL ] cmt="CENT PERIPH"

[ MAIN ]
double TVCL  = THETA1;
double ASCCL = TVCL*pow(IBW/70,0.75);
double CL    = ASCCL*exp(ETA(1));

double TVV1  = THETA2;
double ASCV1 = TVV1*WGT/70.0;
double V1    = ASCV1*exp(ETA(2));

double Q = THETA3;

double TVV2  = THETA4;
double ASCV2 = TVV2*WGT/70.0;
double V2    = ASCV2;

[ TABLE ]
capture IPRED = CENT/V1;
double  DV = IPRED+IPRED*EPS(1)+EPS(2);

capture fIPRED = fu*IPRED;
capture fDV = fu*DV;

[ CAPTURE ] DV

[ THETA ] @annotated
1.21E+01 : typical clearance (L/h)
6.81E+01 : typical volume - cent (L)
2.03E+01 : intercompartmental clearance (L/h) 
4.46E+01 : volume - periph (L)

[ OMEGA ] @annotated
ETA_CL : 6.33E-02 : ETA on CL
ETA_V1 : 7.25E-02 : ETA on V1

[ SIGMA ]
3.40E-02
4.65E-03

[ SET ] delta=1, end=24*3
```
