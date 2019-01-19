[PROB]
1: Krzyzanski W, Wiczling P, Lowe P, Pigeolet E, Fink M, Berghout A, Balser S.
Population modeling of filgrastim PK-PD in healthy adults following intravenous
and subcutaneous administrations. J Clin Pharmacol. 
2010 Sep;50(9 Suppl):101S-112S. doi: 10.1177/0091270010376966. 
PubMed PMID: 20881223.

https://www.ncbi.nlm.nih.gov/pubmed/20881223


PKPD NEUPOGEN BASELINE(BAS) at 0.0246 ng/ml


[SET] end=120, delta=0.1

[PARAM] ROUT  = 0, BAS = 0.0246

[CMT]
ABS
CENT
FRR
DR
B1 B2 B3 B4 B5 B6 B7 B8 B9
NB

[THETA]
6.26E-01  //1
6.42E-01  //2
1.00E+00  //3
6.77E+00  //4
1.48E-01  //5
2.56E+00  //6
1.27E+00  //7
1.01E-01  //8
2.11E-01  //9
0.00E+00  //10
7.23E-02  //11
0.00E+00  //12
1.02E-02  //13
1.65E+00  //14
3.21E+00  //15
3.43E+01  //16
3.23E+01  //17
0.00E+00  //18

[MAIN]
double FF    = THETA1;
double KA    = THETA2*exp(ETA(4));
double FR    = THETA3;
double D2    = THETA4;
double KEL   = THETA5*exp(ETA(2));
double VD    = THETA6*exp(ETA(3));
double KD    = THETA7;
double KINT  = THETA8;
double KSI   = THETA9*exp(ETA(5));
double KOFF  = THETA10;
double KMT   = THETA11;
double KBB1  = THETA12;
double KTT   = THETA13;
double NB0   = THETA14*exp(ETA(1));
double SC1   = THETA15*exp(ETA(6));
double SM1   = THETA16*exp(ETA(7));
double SM2   = THETA17*exp(ETA(8));
double SM3   = THETA18;
double CP0   = BAS;


D_CENT = D2;

F_ABS  = FF;
F_CENT = 0;

if(ROUT == 1) {
  F_ABS  = 0;
  F_CENT = 1;
};

double KON   = KOFF / KD;

// HILL FUNCTION AT BASELINE
double H10     = 1.0 + SM1 * CP0 / (SC1 + CP0);
double H20     = 1.0 + SM2 * CP0 / (SC1 + CP0);
double H30     = 1.0 + SM3 * CP0 / (SC1 + CP0);

// SECONDARY RELATIONSHIPS
double KINB = KMT * NB0 / H10;
double BM10 = KINB * H10 / (KBB1 * H30 + KTT * H20);
double BM20 = KTT * H20 * BM10/(KBB1 * H30 + KTT * H20);
double BM30 = KTT * H20 * BM20/(KBB1 * H30 + KTT * H20);
double BM40 = KTT * H20 * BM30/(KBB1 * H30 + KTT * H20);
double BM50 = KTT * H20 * BM40/(KBB1 * H30 + KTT * H20);
double BM60 = KTT * H20 * BM50/(KBB1 * H30 + KTT * H20);
double BM70 = KTT * H20 * BM60/(KBB1 * H30 + KTT * H20);
double BM80 = KTT * H20 * BM70/(KBB1 * H30 + KTT * H20);
double BM90 = KTT * H20 * BM80/(KBB1 * H30 + KTT * H20);

double NT0   = BM10+BM20+BM30+BM40+BM50+BM60+BM70+BM80+BM90+NB0;
double AC0   = CP0 * VD;
double RTOT0 = KSI*NT0;
double ADR0  = RTOT0*AC0/(KD+CP0);
double KIN   = KEL * AC0 + KINT * ADR0;

//;----- INITIAL CONDITIONS
ABS_0  = 0;
CENT_0 = AC0+ADR0;
FRR_0  = 0;
DR_0   = 0;
B1_0   = BM10;
B2_0   = BM20;
B3_0   = BM30;
B4_0   = BM40;
B5_0   = BM50;
B6_0   = BM60;
B7_0   = BM70;
B8_0   = BM80;
B9_0   = BM90;
NB_0   = NB0;

[ODE]

double ATOT = CENT;
double ADR  = DR;
double NT = B1+B2+B3+B4+B5+B6+B7+B8+B9+NB;
double RTOT = KSI*NT;
double BB = RTOT-ATOT/VD+KD;
double CP = 0.5*(-BB+sqrt(pow(BB,2.0)+4*KD*ATOT/VD));
double AC = CP*VD;
ADR = RTOT*AC/(KD+CP);
double H1 = 1.0 + SM1 * CP / (SC1 + CP);
double H2 = 1.0 + SM2 * CP / (SC1 + CP);
double H3 = 1.0 + SM3 * CP / (SC1 + CP);


//;-------- SET OF DIFFERENTIAL EQUATIONS:
dxdt_ABS  = - KA * ABS;
dxdt_CENT =   KA * ABS  + KIN - KEL * AC - KINT * ADR;
dxdt_FRR  =   0;
dxdt_DR   =   0;
dxdt_B1   =  KINB* H1*(B1/BM10) - KBB1 * H3 * B1 - KTT * H2 * B1;
dxdt_B2   =  KTT * H2 * B1 - KBB1 * H3 * B2 - KTT * H2 * B2;
dxdt_B3   =  KTT * H2 * B2 - KBB1 * H3 * B3 - KTT * H2 * B3;
dxdt_B4   =  KTT * H2 * B3 - KBB1 * H3 * B4 - KTT * H2 * B4;
dxdt_B5   =  KTT * H2 * B4 - KBB1 * H3 * B5 - KTT * H2 * B5;
dxdt_B6   =  KTT * H2 * B5 - KBB1 * H3 * B6 - KTT * H2 * B6;
dxdt_B7   =  KTT * H2 * B6 - KBB1 * H3 * B7 - KTT * H2 * B7;
dxdt_B8   =  KTT * H2 * B7 - KBB1 * H3 * B8 - KTT * H2 * B8;
dxdt_B9   =  KTT * H2 * B8 - KBB1 * H3 * B9 - KTT * H2 * B9;
dxdt_NB   =  KTT * H2 * B9 + KBB1 * H3 *(B1+B2+B3+B4+B5+B6+B7+B8+B9)-KMT*NB;

[TABLE]
  // ZNB   = A(14)
  // ZNT=A(5)+A(6)+A(7)+A(8)+A(9)+A(10)+A(11)+A(12)+A(13)+A(14)
  // RRTOT    = KSI*ZNT
  // BBB=RRTOT-A(2)/VD+KD
  // ZCP=0.5*(-BBB+SQRT(BBB**2+4*KD*A(2)/VD))
//;----------SOLUTIONS:
double ZNB   = NB;
double ZNT   = B1+B2+B3+B4+B5+B6+B7+B8+B9+NB;
double RRTOT = KSI*ZNT;
double BBB   = RRTOT - CENT/VD + KD;
double ZCP   = 0.5*(-BBB+sqrt(pow(BBB,2.0) + 4*KD*CENT/VD));

// ;----- OUTPUT:
CP = ZCP*(1+EPS(2))+EPS(1);
capture RESP = ZNB*(1+EPS(4))+EPS(3);
capture DV = 0;
if(self.cmt==2) DV = CP;
if(self.cmt==4) DV = RESP;

[CAPTURE] CP


[OMEGA]
2.98E-01
3.12E-01
3.28E-01
0.00E+00
2.24E-01
8.03E-01
1.28E-02

[SIGMA]
0 0 0 0
