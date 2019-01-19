[PROB]
1: Li C, Kuti JL, Nightingale CH, Nicolau DP. Population pharmacokinetic analysis
and dosing regimen optimization of meropenem in adult patients. J Clin Pharmacol.
2006 Oct;46(10):1171-8. PubMed PMID: 16988206.

https://www.ncbi.nlm.nih.gov/pubmed/16988206

[SET] delta=0.1, end=8, req=""

[PKMODEL] cmt = "CENT, PERIPH"  

[PARAM] @annotated
WT   : 70 : Weight (kg)
CLCR : 83 : Creatinine clearance (ml/min)
AGE  : 35 : Age (years)
 
[THETA] @annotated
 1.50E+01 : Typical value of clearance (L/h)
 1.27E+01 : Typical value of volume 1 (L)
 1.52E+01 : Intercompartmental clearance (L/h) 
 1.24E+01 : Typical value of volume 2 (L) 
-4.47E-01 : AGE on CL
 8.20E-01 : WT on V1
 1.88E-01 : Proportional error standard deviation
 4.76E-01 : Additive error standard deviation
 6.20E-01 : CLCR on CL 

[MAIN]

double TVCL     = THETA1;
double TVV1     = THETA2;
double TVQ      = THETA3;
double TVV2     = THETA4;
double CL_AGE   = THETA5;
double V1_WT    = THETA6;
double RUV_PROP = THETA7;
double RUV_ADD  = THETA8;
double CL_CLCR  = THETA9;

double LOGTWT = log((WT/70.0)); 
  
double LOGTAGE = log((AGE/35.0));
  
double LOGTCLCR = log((CLCR/83.0));
  
double MU_1 = log(TVCL) + CL_AGE * LOGTAGE + CL_CLCR * LOGTCLCR;

double CL =  exp(MU_1 +  ETA(1)) ;

double MU_2 = log(TVV1) + V1_WT * LOGTWT;
double V1 =  exp(MU_2 +  ETA(2)) ;

double MU_3 = log(TVQ);
double Q =  exp(MU_3 +  ETA(3)) ;

double MU_4 = log(TVV2);
double V2 =  exp(MU_4 +  ETA(4));

[OMEGA] @annotated
ECL : 8.84E-02 : ETA on CL
EV1 : 9.76E-02 : ETA on V1
EQ  : 1.03E-01 : ETA on Q
EV2 : 7.26E-02 : ETA on V2

[SIGMA] @annotated
EP : 1 : Not used

[TABLE] 
capture CC = (CENT/V1);
double IPRED = CC;
double W = sqrt((RUV_ADD*RUV_ADD)+ (RUV_PROP*RUV_PROP*IPRED*IPRED));
capture Y = IPRED+W*EPS(1);
