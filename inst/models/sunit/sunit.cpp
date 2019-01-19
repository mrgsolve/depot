[ PROB ]
1: Houk BE, Bello CL, Kang D, Amantea M. A population pharmacokinetic
meta-analysis of sunitinib malate (SU11248) and its primary metabolite (SU12662) 
in healthy volunteers and oncology patients. Clin Cancer Res. 2009 Apr
1;15(7):2497-506. doi: 10.1158/1078-0432.CCR-08-1893. Epub 2009 Mar 3. 
PubMed PMID: 19258444.

$PARAM
TVCL = 51.8
TVVC = 2030
TVKA = 0.195
TVQ = 7.22
TVVP = 583
WTVC = 0.459
SEXCL = -0.0876
ASIANCL = -0.130
GISTCL = -0.285
SOLIDCL = -0.269
MRCCCL = -0.258
SEX = 0, ASIAN = 0, GIST = 0
SOLID = 0, MRCC = 0, WT = 76.9

$MAIN
double CL  = TVCL * (1+SEXCL*SEX) * (1+ASIANCL*ASIAN) * 
  (1+GISTCL*GIST) * (1+SOLIDCL*SOLID) * (1+MRCCCL*MRCC) * exp(ETA(1));

double V2 = TVVC*pow(WT/76.9, WTVC)*exp(ETA(2));
double KA = TVKA*exp(ETA(3));
double Q  = TVQ;
double V3 = TVVP;

$OMEGA 0.14 0.18 0.64

$SIGMA 0.146

$PKMODEL cmt = "GUT CENT, PERIPH", depot = TRUE

$POST
capture CP = (1000*CENT/V2);

