
[PROB]
# Secukinumab PK/PD 
  - FDA Clin Pharm Review 125504Orig1s000ClinPharmR.pdf
  - Clinical Pharmacology Reviewer: Jie Wang, Ph.D.
  - Pharmcometrics Reviewer: Jee Eun Lee, Ph.D.
  - Pharmacometrics Team Leader: Jeffry Florian, Ph.D.
  - Clinical Pharmacology Team Leader: Yow-Ming Wang, Ph.D.
  
https://www.accessdata.fda.gov/drugsatfda_docs/nda/2015/125504Orig1s000ClinPharmR.pdf

[SET] req = "RESP", end=84, delta=1

[CMT] @annotated
SC     : SC dosing compartment (mg)
CENT   : Central compartment (mg) 
PERIPH : Peripheral compartment (mg)
RESP   : PASI score

[PARAM] @annotated
WT   : 90 : Patient weight (kg)
BASE : 21 : Baseline PASI

[GLOBAL]
inline double antilogit(const double a) {
  return exp(a)/(1.0+exp(a));
}

#define CP  (CENT/(V2/1000))

[THETA] @annotated
0.19 : Typical value of clearance (L/day)
3.61 : Typical value of central volume (L)
0.39 : Typical value of intercomp clearance (L/day)
2.87 : Typical value of peripheral volume (L)
0.18 : Absorption rate constant <sc dosing> (1/day)
1.00 : Weight effect on CL
0.81 : Weight effect on V2
0.68 : Weight effect on Q
0.56 : Weight effect on V3
0.99 : Bioav from SC dose, logit scale

[THETA] @name theta @annotated
0.0012   : Typical value of PASI kout
19150.5  : Typical value of EC50 (ng/ml)
56.38    : Typical value of PASI Emax
2.07     : Hyperbolic parameter for PASI Emax model


[MAIN]
double  CL = THETA1*pow(WT/90.0,THETA6)*exp(ECL);
double  V2 = THETA2*pow(WT/90.0,THETA7)*exp(EVC);
double  Q  = THETA3*pow(WT/90.0,THETA8); 
double  V3 = THETA4*pow(WT/90.0,THETA9)*exp(EVP);
double  KA = THETA5*exp(EKA);

double KOUT = theta1*exp(EKOUT);
double EC50 = theta2*exp(EEC);
double EMAX = theta3;
double GAM  = theta4;

F_SC = antilogit(THETA10);

double KIN = BASE*KOUT;
RESP_0 = BASE;

[ODE]
dxdt_SC     = -KA*SC;
dxdt_CENT   =  KA*SC - (CL+Q)*CENT/V2 + Q*PERIPH/V3;
dxdt_PERIPH =  Q*(CENT/V2 - PERIPH/V3);

double DE   =  EMAX*pow(CP,GAM)/(pow(CP,GAM) + pow(EC50,GAM));
dxdt_RESP   =  KIN - KOUT*(1.0+DE)*RESP;


[OMEGA] 
@labels ECL EVC EVP EKA
@cor 
0.1024
0.7 0.09 
0.14 0.72 0.0324
0 0 0 0.1225


[OMEGA]
@labels EKOUT EEC ESL
0.5776 1 0.000196

[SIGMA] @labels PROP ADD 
0.0289 137641

[SIGMA] @labels PDADD PDPROP
0.14822 0.1177862

[TABLE]
double IPRED = CP/1000;
double DV =  (CP*(1.0+PROP)+ADD)/1000;
double PASI = RESP*(1.0+PDPROP)+PDADD;

if(NEWIND <=1) {
  double base = PASI; 
}
double dPASI = PASI/base;

[CAPTURE] @annotated
IPRED : Individual predicted secuk concentration (mg/L)
DV    : IPRED with residual variability  (mg/L)
PASI  : PASI score with residual variability (.)
dPASI : PASI score, change from baseline (.)
base  : Real baseline PASI
  
