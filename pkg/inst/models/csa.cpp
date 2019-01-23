[PROB]
# Yoshikado et al. (2016)

1: Yoshikado T, Yoshida K, Kotani N, Nakada T, Asaumi R, Toshimoto K, Maeda K,
Kusuhara H, Sugiyama Y. Quantitative Analyses of Hepatic OATP-Mediated
Interactions Between Statins and Inhibitors Using PBPK Modeling With a Parameter 
Optimization Method. Clin Pharmacol Ther. 2016 Nov;100(5):513-523. 
doi: 10.1002/cpt.391. Epub 2016 Jul 28. PubMed PMID: 27170342.

https://www.ncbi.nlm.nih.gov/pubmed/27170342
  
[ CMT ] @annotated
gut : Gut compartment
cent : CsA central compartment
me : Muscle extracellular
se : Skin extracullular
ae : Adipose extracellular
mc : Muscle intracellular
sc : Skin intracellular
ac : Adipose intracellular
liv1 : First  liver compartment 
liv2 : Second liver compartment
liv3 : Third  liver compartment
liv4 : Fourth liver compartment
liv5 : Fifth  liver compartment

[ PARAM ] 
Kp_mus = 2.98
Kp_adi = 17.3
Kp_ski = 13.6
Kp_liv = 16.7
fb     = 0.06
PSmus = 245/60
PSski = 37.4/60
PSadi = 10.2/60
fhCLint = 0.587/60

fafg = 0.572
Clr = 0
ka = 0.999
tlag =  0.254 
wt = 60

[PARAM] @annotated
Qh   : 1.200 : Hepatic blood flow (L/hr/kg)
Qmus : 0.642 : Muscle blood flow (L/hr/kg)
Qski : 0.257 : Skin blood flow (L/hr/kg)
Qadi : 0.223 : Adipose blood flow (L/hr/kg)

Vcent : 0.075 : Volume central (L/kg)
Vliv : 0.0241 : Volume liver (L/kg)
Vmus : 0.4290 : Volume muscle (L/kg)
Vski : 0.1110 : Volume skin (L/kg)
Vadi : 0.1430 : Volume adipose tissue (L/kg)

exFliv : 0.278 : Liver fraction extracellular
exFmus : 0.146 : Muscle fraction extracellular
exFski : 0.321 : Skin fraction extracellular
exFadi : 0.145 : Adipose fraction extracellular


[MAIN]

if(NEWIND <=1) {
  double Vme = Vmus*exFmus;
  double Vae = Vadi*exFadi;
  double Vse = Vski*exFski;
  double Vmc = Vmus-Vme;
  double Vac = Vadi-Vae;
  double Vsc = Vski-Vse;
  double dVliv = Vliv/5.0;
}

[ODE]
// CsA concentrations
double Ccent = cent/Vcent;

double Cme = me/Vme;
double Cse = se/Vse;
double Cae = ae/Vae;

double Cmc = mc/Vmc;
double Csc = sc/Vsc;
double Cac = ac/Vac;

double Cliv1 = liv1/dVliv;
double Cliv2 = liv2/dVliv;
double Cliv3 = liv3/dVliv;
double Cliv4 = liv4/dVliv;
double Cliv5 = liv5/dVliv;

dxdt_gut = -ka/fafg*gut;

dxdt_cent = 
  Qh*Cliv5/Kp_liv 
  - Qh*Ccent 
  - Clr*Ccent 
  - Qmus*(Ccent-Cme) 
  - Qski*(Ccent-Cse) 
  - Qadi*(Ccent-Cae);
  
dxdt_me = Qmus*(Ccent-Cme) - PSmus*fb*(Cme-Cmc/Kp_mus);
dxdt_se = Qski*(Ccent-Cse) - PSski*fb*(Cse-Csc/Kp_ski);
dxdt_ae = Qadi*(Ccent-Cae) - PSadi*fb*(Cae-Cac/Kp_adi);
  
dxdt_mc = PSmus*fb*(Cme-Cmc/Kp_mus);
dxdt_sc = PSski*fb*(Cse-Csc/Kp_ski);
dxdt_ac = PSadi*fb*(Cae-Cac/Kp_adi);
  
dxdt_liv1 = Qh*(Ccent-Cliv1/Kp_liv) - (fhCLint/5.0)*Cliv1 + ka*gut;
dxdt_liv2 = Qh*(Cliv1-Cliv2)/Kp_liv - (fhCLint/5.0)*Cliv2;
dxdt_liv3 = Qh*(Cliv2-Cliv3)/Kp_liv - (fhCLint/5.0)*Cliv3;
dxdt_liv4 = Qh*(Cliv3-Cliv4)/Kp_liv - (fhCLint/5.0)*Cliv4;
dxdt_liv5 = Qh*(Cliv4-Cliv5)/Kp_liv - (fhCLint/5.0)*Cliv5;
  


[TABLE]

capture CSA = cent/Vcent;
capture CSAliv = liv1/dVliv;

