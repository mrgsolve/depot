
[PARAM] @annotated
TVCL : 0.168  : Clearance (L/d)
TVVC : 2.88   : Central volume of distribution (L)
TVQ  : 0.384  : Intercompartmental clearance (L/d)
TVVP : 2.85   : Peripheral volume of distribution (L)
TVVMAX : 0.114  : Maximum reaction velocity (mg/d)
KM   : 0.0784 : Michaelis constant (mg/L)
BASE : 2.09   : Baseline 
IMAX : 0.961  : Maximal inhibitory activity 
IC50 : 0.535  : Concentration required for 50% inhibition
WT   : 78     : Patient weight (kg)
  
[CMT] @annotated
CENT    : Central compartment (mg)
PERIPH  : Peripheral compartment (mg)
  
[MAIN]
double CL = TVCL*pow(WT/78,0.75);
double VC = TVVC*(WT/78);
double Q  = TVQ *pow(WT/78,0.75);
double VP = TVVP*(WT/78);
double VMAX = TVVMAX*exp(EVMAX);
                       
[OMEGA] @annotated
EVMAX : 0 : IIV on VMAX

[SIGMA] @annotated
EPS_PK : 0.0876 : RUV PK
EPS_PD : 0.209  : RUV PD

  
[ODE]
double CP = CENT/VC;
double CLNL = VMAX/(KM+CP);

dxdt_CENT = -(CL+Q+CLNL)*CP + Q*PERIPH/VP;
dxdt_PERIPH = Q*(CP - PERIPH/VP);
  
[TABLE]

CP = CENT/VC;
//EMAX-(EMAX-BASE)*(CP**GAM/(EC50**GAM+CP**GAM))
capture PD = IMAX-(IMAX-BASE)*(CP/(IC50+CP));


[CAPTURE] CP CL CLNL