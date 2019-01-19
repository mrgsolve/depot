$PROB
1: Sadiq MW, Nielsen EI, Khachman D, Conil JM, Georges B, Houin G, Laffont CM,
Karlsson MO, Friberg LE. A whole-body physiologically based pharmacokinetic
(WB-PBPK) model of ciprofloxacin: a step towards predicting bacterial killing at 
sites of infection. J Pharmacokinet Pharmacodyn. 2017 Apr;44(2):69-79. doi:
10.1007/s10928-016-9486-9. Epub 2016 Aug 30. PubMed PMID: 27578330; PubMed
Central PMCID: PMC5376394.

$SET atol = 1E-9, rtol = 1E-9

$GLOBAL
#define A1   ART
#define A2   VEN
#define A3   LUN
#define A4   BRA
#define A5   HRT
#define A6   SKN
#define A7   MUS
#define A8   ADI
#define A9   SPL
#define A10  GIO
#define A11  LIV
#define A12  KID
#define A13  RES
#define A14  ELI

$CMT @annotated
ART : arteries
VEN : veins
LUN : lungs
BRA : brain
HRT : heart
SKN : skin
MUS : muscle
ADI : adipose
SPL : spleeen
GIO : GIO
LIV : liver
KID : kidneys
RES : rest
ELI : elimiated drug

$PARAM @annotated
WT  :  77  : weight (kg)
SEX :   0  : sex
CRCL:  82  : creatinine clearance (ml/min/1.73m2)
FUP : 0.65 : unbound fraction

$MAIN
if(SEX!=1.0) {
  double VART = 0.017083*(WT)   ; // Arterial volume (L)
  double VVEN = 0.051248*(WT)   ; // Venous volume (L)
  double VLUN = 0.00643836*(WT) ; // Lungs weight (Kg) ICRP
  double VBRA = 0.0191781*(WT)  ; // Brain volume (L)
  double VHRT = 0.004167*(WT)   ; // Heart volume (L) ICRP
  double VSKN = 0.0383*(WT/1.18); // Skin volume (L)
  double VMUS = 0.2916*(WT)     ; // Muscle volume (L) ICRP
  double VADI = 0.3*(WT/0.916)  ; // Adipose volume (L) ICRP
  double VSPL = 0.00247*(WT)    ; // Spleen volume (L) ICRP
  double VGIO = 0.01644*(WT)    ; // GIO volume (L) stomach+SI colon ICRP
  double VHEP = 0.020724*(WT)   ; // Hepatic volume (L) ICRP
  double VKID = 0.0042466*(WT)  ; // Kidney volume (L) ICRP
  double VRES = 0.23220194*(WT) ; // Rest of body volume bones ICRP
}
if(SEX==1) {
  VART = 0.01918*(WT)       ; // Arterial volume (L)
  VVEN = 0.05753*(WT)       ; // Venous volume (L)
  VLUN = 0.00643836*(WT)    ; // Lungs weight (Kg) ICRP
  VBRA = 0.01918*(WT)       ; // Brain volume (L)
  VHRT = 0.004521*(WT)      ; // Heart volume (L) ICRP
  VSKN = 0.045205*(WT/1.18) ; // Skin volume (L)
  VMUS = 0.3973*(WT)        ; // Muscle volume (L) ICRP
  VADI = 0.171233*(WT/0.916); // Adipose volume (L) ICRP
  VSPL = 0.00247*(WT)       ; // Spleen volume (L) ICRP
  VGIO = 0.01644*(WT)       ; // GIO volume (L) stomach+SI colon ICRP
  VHEP = 0.02466*(WT)       ; // Hepatic volume (L) ICRP
  VKID = 0.004247*(WT)      ; // Kidney volume (L) ICRP
  VRES = 0.23612164*(WT)    ; // Rest of body volume bones ICRP
}

double CO = (15*pow(WT,0.74)); // Cardiac output in L/h
double QART = CO             ; // Artery Blood flow
double QVEN = CO             ; // Venous blood flow
double QLUN = CO             ; // Lung blood flow

if(SEX!=1.0) {
  double QBRA = CO*(0.12) ; // Brain blood flow
  double QHRT = CO*(0.05) ; // Heart blood flow
  double QSKN = CO*(0.05) ; // Skin blood flow
  double QMUS = CO*(0.12) ; // Muscle blood flow
  double QADI = CO*(0.085); // Adipose blood flow
  double QSPL = CO*(0.03) ; // Spleen blood flow
  double QGIO = CO*(0.16) ; // GIO blood flow
  double QKID = CO*(0.17) ; // Kidney blood flow
  double QHEPA= CO*(0.065); // Hepatic arterial blood flow
  double QRES = CO*(0.15) ; // Rest of Body blood flow
} 
if(SEX==1.0) { 
  QBRA = CO*(0.12) ;  // Brain blood flow
  QHRT = CO*(0.04) ;  // Heart blood flow
  QSKN = CO*(0.05) ;  // Skin blood flow
  QMUS = CO*(0.17) ;  // Muscle blood flow
  QADI = CO*(0.05) ;  // Adipose blood flow
  QSPL = CO*(0.03) ;  // Spleen blood flow
  QGIO = CO*(0.15) ;  // GIO blood flow
  QKID = CO*(0.19) ;  // Kidney blood flow
  QHEPA= CO*(0.065);  // Hepatic arterial blood flow
  QRES = CO*(0.135);  // Rest of Body blood flow
}

double QHEPT= QHEPA + QSPL + QGIO ; // Total hepatic bloodflow HEPA+SPL+GIO)

double TVCLH = clh ; //Hepatic or Non-renal CL
double CLH = TVCLH*exp(ETA(1));
double TVRSEC = resch ; //Factor renal secretion
double RSEC = TVRSEC;
double CRCL2 = CRCL > 150.0 ? 150.0 : CRCL;
double CLR = (((CRCL2*60/1000.0)*FUP)*(1+RSEC))*exp(ETA(1));
double CL = CLR+CLH;

double KLUN = klun*exp(ETA(2));
double KBRA = kbra*exp(ETA(2));
double KHRT = khrt*exp(ETA(2));
double KSKN = kskn*exp(ETA(2));
double KMUS = kmus*exp(ETA(2));
double KADI = kadi*exp(ETA(2));
double KSPL = kspl*exp(ETA(2));
double KGIO = kgio*exp(ETA(2));
double KHEP = khep*exp(ETA(2));
double KKID = kkid*exp(ETA(2));
double KRES = kres*exp(ETA(2));

$PARAM @annotated
clh   : exp(2.60)   : Hepatic clearance
resch : exp(0.674)  : Renal secretion
klun  : exp(1.20)   : KLUN
kbra  : exp(-0.257) : KBRA
khrt  : exp(1.30)   : KHRT
kskn  : exp(-0.335) : KSKN
kmus  : exp(0.0229) : KMUS
kadi  : exp(-0.885) : KADI
kspl  : exp(0.668)  : KSPL
kgio  : exp(1.21)   : KGIO
khep  : exp(1.27)   : KHEP
kkid  : exp(2.09)   : KKID
kres  : exp(1.35)   : KRES

$ODE

double C1 = A1/VART ;  // AMT/Arterial volume
double C2 = A2/VVEN ;  // AMT/Venous volume
double C3 = A3/VLUN ;  // AMT/Lung volume
double C4 = A4/VBRA ;  // AMT/Brain volume
double C5 = A5/VHRT ;  // AMT/Heart volume
double C6 = A6/VSKN ;  // AMT/Skin volume
double C7 = A7/VMUS ;  // AMT/Muscle volume
double C8 = A8/VADI ;  // AMT/Adipose volume
double C9 = A9/VSPL ;  // AMT/Spleen volume
double C10= A10/VGIO;  // AMT/GIO volume
double C11= A11/VHEP;  // AMT/Hepatic volume
double C12= A12/VKID;  // AMT/Kidney volume
double C13= A13/VRES;  // AMT/Rest volume

double VENIN1 = (C4*QBRA/KBRA)+(C5*QHRT/KHRT)+(C6*QSKN/KSKN)+(C7*QMUS/KMUS);
double VENIN2 = (C8*QADI/KADI)+(C11*QHEPT/KHEP)+(C12*QKID/KKID)+(C13*QRES/KRES);
double VENOUT = (QLUN*C2);
double HEPIN  = (C1*QHEPA)+(C9*QSPL/KSPL)+(C10*QGIO/KGIO);
double HEPOUT = (C11*QHEPT/KHEP)+(C11*CLH*FUP/KHEP);
double CAB    = C3*0.79 ; 

dxdt_ART = (QLUN*C3/KLUN)-(QART*C1);
dxdt_VEN = VENIN1+VENIN2-VENOUT;
dxdt_LUN = (QVEN*C2)-(QLUN*C3/KLUN);
dxdt_BRA = (QBRA*C1)-(QBRA*C4/KBRA);
dxdt_HRT = (QHRT*C1)-(QHRT*C5/KHRT);
dxdt_SKN = (QSKN*C1)-(QSKN*C6/KSKN);
dxdt_MUS = (QMUS*C1)-(QMUS*C7/KMUS);
dxdt_ADI = (QADI*C1)-(QADI*C8/KADI);
dxdt_SPL = (QSPL*C1)-(QSPL*C9/KSPL);
dxdt_GIO = (QGIO*C1)-(QGIO*C10/KGIO);
dxdt_LIV = HEPIN-HEPOUT;
dxdt_KID = (QKID*C1)-(QKID*C12/KKID)-C1*CLR;
dxdt_RES = (QRES*C1)-(QRES*C13/KRES);
dxdt_ELI = C1*CLR+(C11*CLH*FUP/KHEP);

$TABLE
double CC1  = A1/VART*FUP  ; // AMT/Arterial volume
double CC2  = A2/VVEN*FUP  ; // AMT/Venous volume
double CC3  = A3/VLUN*0.79 ; // AMT/Lung volume
double CC4  = A4/VBRA*0.79 ; // AMT/Brain volume
double CC5  = A5/VHRT*0.79 ; // AMT/Heart volume
double CC6  = A6/VSKN*0.65 ; // AMT/Skin volume
double CC7  = A7/VMUS*0.79 ; // AMT/Muscle volume
double CC8  = A8/VADI*0.79 ; // AMT/Adipose volume
double CC9  = A9/VSPL*0.79 ; // AMT/Spleen volume
double CC10 = A10/VGIO*0.67; // AMT/GIO volume
double CC11 = A11/VHEP*0.79; // AMT/Hepatic volume
double CC12 = A12/VKID*0.79; // AMT/Kidney volume
double CC13 = A13/VRES*0.79; // AMT/Rest volume

$CAPTURE CP = A2/VVEN Ckid = CC12 Clung = CC3
