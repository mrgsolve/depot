[ PROB ]
1: Yan X, Lowe PJ, Fink M, Berghout A, Balser S, Krzyzanski W. Population
pharmacokinetic and pharmacodynamic model-based comparability assessment of a
recombinant human Epoetin Alfa and the Biosimilar HX575. J Clin Pharmacol. 2012
Nov;52(11):1624-44. doi: 10.1177/0091270011421911. Epub 2011 Dec 12. 
PubMed PMID: 22162538; PubMed Central PMCID: PMC3483448.

https://www.ncbi.nlm.nih.gov/pubmed/22162538

$CMT @annotated
A1 : EPO SC
A2 : EPO SERUM
A3 : EPO TISSUE
A4 : PRECURSOR 1
A5 : PRECURSOR 2
A6 : PRECURSOR 3
A7 : RET
A8 : RBC
A9 : HGB


$OMEGA @annotated
ED1   : 4.20E-01 : .
EKA   : 1.11E-01 : .
EVC   : 9.59E-02 : .
ECL   : 2.92E-01 : .
EC0   : 4.17E-02 : .
ERET0 : 7.60E-02 : .
ETP1  : 0.00E+00 : .
ESC50 : 0.00E+00 : .
ESMAX : 8.44E-01 : .
EKD   : 0.00E+00 : .
EKSAI : 1.80E-01 : .
ERBC0 : 2.84E-03 : .
EMCH  : 1.86E-03 : .

$SIGMA
5.24E-02 0 7.70E-02  1.11E-03

$THETA @annotated
2.64E+00 : D1
3.01E-02 : KA
4.51E+00 : VC
3.74E-01 : CL
5.16E-01 : F1
0.00E+00 : ALAG1  
1.25E+01 : C0 
7.24E-02 : K23  
1.70E-01 : K32
4.90E+01 : RET0  
0.00E+00 : TD3
5.80E+01 : TP1
0.00E+00 : TP2  
0.00E+00 : TP3  
6.68E+00 : SC50
7.56E-02 : SMAX  
4.12E+02 : KD  
1.32E+00 : KSAI  
4.63E+00 : RBO  
3.03E+01 : MCH


$MAIN
D_A1        = THETA1*exp(ED1);
double KA   = THETA2*exp(EKA);
double VC   = THETA3*exp(EVC);
double CL   = THETA4*exp(ECL);
F_A1        = THETA5;
ALAG_A1     = THETA6;
double C0   = THETA7*exp(EC0);
double K23  = THETA8;
double K32  = THETA9;
double RR0  = THETA10*exp(ERET0);
double TD3  = THETA11;
double TP1  = THETA12*exp(ETP1);
double TP2  = THETA13;
double TP3  = THETA14;
double SC50 = THETA15*exp(ESC50);
double SMAX = THETA16*exp(ESMAX);
double KD   = THETA17*exp(EKD);
double KSAI = THETA18*exp(EKSAI);
double RB0  = THETA19*exp(ERBC0);
double MCH  = THETA20*exp(EMCH);
double TR   = TP1;
TP3  = TP1;
TP2  = TP1;
double TB   = (RB0-RR0/1000.0)*TR*1000.0/RR0;

A1_0 = 0;
A2_0 = C0*VC;
A3_0 = C0*VC*K23/K32;
double RTOT0 = RR0*TP2/TR;
double RC0 = RTOT0*C0/(KD+C0);
double SC0 = SMAX*RC0/(SC50+RC0);

A4_0 = RR0*TP1/(TR*SC0)*(1.0/pow(2.0,5.0));
A5_0 = RR0*TP2/TR;
A6_0 = RR0*TP3/TR;
A7_0 = RR0;
A8_0 = RB0 - RR0/1000.0;

$ODE
double RTOT = A5;
double KEPO = CL*C0+RTOT0*KSAI*VC*C0/(KD+C0);

dxdt_A1 = -KA*A1;
dxdt_A2 = KEPO + KA*A1 - CL*A2/VC - K23*A2 + K32*A3 - 
  RTOT*KSAI*VC*(A2/VC)/(KD+A2/VC);
dxdt_A3 = K23*A2 - K32*A3;

double C2 = A2/VC;
double RC = RTOT*C2/(KD+C2);
double SC = SMAX*RC/(SC50+RC);
double KIN0 = RR0/TR*(1.0/pow(2.0,5.0));

dxdt_A4 = KIN0 - (1.0/TP1)*SC*A4;
dxdt_A5 = pow(2.0,5.0)*(1.0/TP1)*SC*A4 - (1.0/TP2)*A5;
dxdt_A6 = (1.0/TP2)*A5 - (1.0/TP3)*A6;
dxdt_A7 = (1.0/TP3)*A6 - A7/TR;
dxdt_A8 = A7/(TR*1000.0) - A8/TB;
dxdt_A9 = 0;

$TABLE
capture EPOi  = A2/VC;
double EPOdv = EPOi*(1+EPS(1))+EPS(2);

capture RETi = A7;
double RETdv = RETi*(1+EPS(3));

capture RBCi = A8+A7/1000.0;
double RBCdv = RBCi*(1+EPS(4));

capture HGBi = MCH*RBCi/10.0;
double HGBdv = HGBi*(1+EPS(4));

$CAPTURE RETdv

$SET req = ""
  


