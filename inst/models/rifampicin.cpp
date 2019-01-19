[ PROB ] 
1: Asaumi R, Toshimoto K, Tobe Y, Hashizume K, Nunoya KI, Imawaka H, Lee W,
Sugiyama Y. Comprehensive PBPK Model of Rifampicin for Quantitative Prediction of
Complex Drug-Drug Interactions: CYP3A/2C9 Induction and OATP Inhibition Effects. 
CPT Pharmacometrics Syst Pharmacol. 2018 Mar;7(3):186-196. doi:
10.1002/psp4.12275. Epub 2018 Feb 5. PubMed PMID: 29368402; PubMed Central PMCID:
PMC5869557.


[ PARAM ]

Rdif = 0.129
beta = 0.2 // 0.2/0.5/0.8
gamma = 0.778

Km_u_uptake = 0.146 // 1.23 // ug/mL

SFKp = 6.65 // 1 // KpScaleRif

Emax_UGT_RIF = 1.34 // 1.00
EC50_u_UGT_RIF = 0.0526 // ug/mL

Emax_CYP3A4_RIF = 4.57 // 12.3
EC50_u_CYP3A4_RIF = 0.05 // ug/mL

kdeg_UGT_liver = 0.0158 // per h
kdeg_CYP3A4_liver = 0.0158 // per h
kdeg_UGT_ent = 0.0288 // per h

fm_UGT_liver = 0.759 // same
fm_UGT_ent = 0.759 // same

[ PARAM ]
fB = 0.0778
fH = 0.0814
fE = 0.115

Fa = 1.000
Fg = 0.943

[ PARAM ]
Kp_skin = 0.326
Kp_muscle = 0.0947
Kp_adipose = 0.0629
Kp_serosa = 0.200

[PARAM]
fBCLint_all_kg = 0.251 // 0.204 // L/h/kg
PSdif_E_kg = 0.161 // 0.143 // L/h/kg
CLrenal_kg = 0.011 // L/h/kg

[ PARAM ]
Qvilli_kg = 0.257 // L/hr/kg
Qh_kg = 1.24 // L/hr/kg
Qmuscle_kg = 0.642 // L/h/kg
Qskin_kg = 0.257 // L/h/kg
Qadipose_kg = 0.223 // L/h/kg
Qserosa_kg = 0.274 // L/h/kg

[ PARAM ]
VHE_kg = 0.0067 // L/kg (Vi)
VHC_kg = 0.0174 // L/kg (Vh)
Vcentral_kg = 0.0743 // L/kg (VbRif)
Vskin_kg = 0.111 // L/kg
Vadipose_kg = 0.143 // L/kg
Vmuscle_kg = 0.429 // L/kg
Vserosa_kg = 0.00893 // L/kg
Vent_kg = 0.00739 // L/kg
Vmucblood_kg = 0.00099 // L/kg

[PARAM]
ka = 37.6 // 3.26 // per hr
WT  = 80


[CAPTURE] Qh VHE VHC Vmax_uptake Ccentral

[ MAIN ]
if(NEWIND <= 1) {
  // -------------------------------------------
  double fBCLint_all = fBCLint_all_kg*WT;
  double CLint_all   = fBCLint_all / fB;
  double CLrenal     = CLrenal_kg*WT;
  double PSdif_E     = PSdif_E_kg*WT;
  // -------------------------------------------
  double Qvilli    = Qvilli_kg*WT;
  double Qh        = Qh_kg*WT;
  double Qmuscle   = Qmuscle_kg*WT;
  double Qskin     = Qskin_kg*WT;
  double Qadipose  = Qadipose_kg*WT;
  double Qserosa   = Qserosa_kg*WT;
  double Qhart     = Qh - Qserosa - Qvilli;
  // -------------------------------------------
  double VHE = VHE_kg*WT;
  double VHC = VHC_kg*WT;
  double Vcentral = Vcentral_kg*WT;
  double Vskin = Vskin_kg*WT;
  double Vadipose = Vadipose_kg*WT;
  double Vmuscle = Vmuscle_kg*WT;
  double Vserosa = Vserosa_kg*WT;
  double Vent = Vent_kg*WT;
  double Vmucblood = Vmucblood_kg*WT;
  // -------------------------------------------
  double Vmax_uptake = 1.0 / (1 + Rdif) * CLint_all / beta * Km_u_uptake;
  double PSdif_inf = Rdif /  (1 + Rdif) * CLint_all / beta;
  double PSdif_eff = Rdif /  (1 + Rdif) * CLint_all / beta / gamma;
  double CLint = Rdif / (1 + Rdif) * CLint_all / (1 - beta) / gamma;
  double Qgut = fE * PSdif_E * Qvilli / (Qvilli + fB * PSdif_E);
  double CLint_E = (Qgut * (1.0 / Fg - 1.0) - (1.0- Fa) * fE * PSdif_E * 20.0) / fE;

}

[CMT]
Xgutlumen
central Cmuscle Cskin Cadipose Cserosa Cmucblood Cent
CHE1 CHE2 CHE3 CHE4 CHE5
CHC1 CHC2 CHC3 CHC4 CHC5

[INIT]
UGT_ratio_HC1 = 1
UGT_ratio_HC2 = 1
UGT_ratio_HC3 = 1
UGT_ratio_HC4 = 1
UGT_ratio_HC5 = 1
UGT_ratio_ent = 1

// CYP3A4_ratio_HC1 = 1
// CYP3A4_ratio_HC2 = 1
// CYP3A4_ratio_HC3 = 1
// CYP3A4_ratio_HC4 = 1
// CYP3A4_ratio_HC5 = 1

[ ODE ]

double Ccentral = central/Vcentral;

// Model Equations for Rifampicin
//   Central compartment

// Muscle: 63
// Skin: 64
// Adipose: 65
// Central: 51
// CHE5: 60

//Vcentral * dCcentral / dt = 
// y51' = 1 / VbRif * ( 
//   Qh * y60 + 
//   Qm * ( y63 / ( KpScaleRif * KpmRif ) - y51 ) + 
//   Qs * ( y64 / ( KpScaleRif * KpsRif ) - y51 ) + 
//   Qa * ( y65 / ( KpScaleRif * KpaRif ) - y51 ) - 
//   ( Qhart + Qser + Qvilli ) * y51 - CLrRif * y51 )
//  )
dxdt_central = 
  Qh       * CHE5 - 
  Qhart    * Ccentral - 
  Qserosa  * Ccentral - 
  Qvilli   * Ccentral - 
  CLrenal  * Ccentral + 
  Qmuscle  * (Cmuscle  / (SFKp * Kp_muscle)  - Ccentral) + 
  Qskin    * (Cskin    / (SFKp * Kp_skin)    - Ccentral) + 
  Qadipose * (Cadipose / (SFKp * Kp_adipose) - Ccentral);

// Distribution compartments (tissue; muscle_ skin_ adipose and serosa)
// Vtissue * dCmuscle / dt = 
dxdt_Cmuscle = 
  (1.0/Vmuscle) * Qmuscle * (Ccentral - Cmuscle / (SFKp * Kp_muscle));

//  Vtissue * dCskin / dt = 
dxdt_Cskin = 
  (1.0/Vskin) * Qskin * (Ccentral - Cskin / (SFKp * Kp_skin));

// Vtissue * dCadipose / dt = 
// y65' = 1 / Va * Qa * ( y51 - y65 / ( KpScaleRif * KpaRif ) )
dxdt_Cadipose = 
  (1.0/Vadipose) * Qadipose * (Ccentral - Cadipose / (SFKp * Kp_adipose));

// Vtissue * dCserosa / dt = 
dxdt_Cserosa = 
  (1.0/Vserosa) * Qserosa * (Ccentral - Cserosa / (SFKp * Kp_serosa));

// Mucosal blood compartment
//  Vmucosal blood * dCmucosal blood / dt = 
// y68' = 1 / Vmucb * ( Qvilli * ( y51 - y68 ) + 
//  fgRif * PSdifentRif * y67 - fbRif * PSdifentRif * y68 ),
dxdt_Cmucblood = 
  Qvilli * (Ccentral - Cmucblood) + 
  fE * PSdif_E * Cent - fB * PSdif_E * Cmucblood;

dxdt_Cmucblood = dxdt_Cmucblood * (1/Vmucblood);

//  Gut lumen compartment
//  dXgut lumen / dt = - 
// y66' = - KaRif / FaRif * y66 + fgRif * PSdifentRif * 20 * y67,
dxdt_Xgutlumen = 
  - ka / Fa * Xgutlumen + fE * PSdif_E * 20 * Cent;
  
//  Enterocyte compartment
//  Vent * dCent / dt = 
// 
// y67' = 1 / Vent * ( KaRif * y66 + fbRif * PSdifentRif * y68 - 
// fgRif * ( PSdifentRif * 21 + CLmetgRif * ( 1 + fmUGTgRif * ( y86 - 1 ) ) ) * y67 ), 
dxdt_Cent = 
  ka * Xgutlumen + 
  fB * PSdif_E * Cmucblood - 
  fE * (PSdif_E * 21 + 
  CLint_E * (1 + fm_UGT_ent * (UGT_ratio_ent - 1))) * Cent;

dxdt_Cent = dxdt_Cent * (1/Vent);

// Extent of enzyme induction in hepatocytes (enzyme; UGT, CYP3A or CYP2C9. HC_i; i = 1~5)

// UGT
//  dEnzymeratio_HC_i / dt = 
dxdt_UGT_ratio_HC1 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC1 / (fH * CHC1 + EC50_u_UGT_RIF) - UGT_ratio_HC1);

dxdt_UGT_ratio_HC2 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC2 / (fH * CHC2 + EC50_u_UGT_RIF) - UGT_ratio_HC2);

dxdt_UGT_ratio_HC3 = 
  kdeg_UGT_liver *
  (1 + Emax_UGT_RIF * fH * CHC3 / (fH * CHC3 + EC50_u_UGT_RIF) - UGT_ratio_HC3);

dxdt_UGT_ratio_HC4 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC4 / (fH * CHC4 + EC50_u_UGT_RIF) - UGT_ratio_HC4);

dxdt_UGT_ratio_HC5 = 
  kdeg_UGT_liver * 
  (1 + Emax_UGT_RIF * fH * CHC5 / (fH * CHC5 + EC50_u_UGT_RIF) - UGT_ratio_HC5);

//  Extent of enzyme induction in ent (enzyme; UGT or CYP3A)
// UGT
// dEnzymeratio_ent / dt = 
dxdt_UGT_ratio_ent = 
  kdeg_UGT_ent * 
  (1 + Emax_UGT_RIF * fE * Cent / (fE * Cent + EC50_u_UGT_RIF) - UGT_ratio_ent);


// CYP3A4
// dxdt_CYP3A4_ratio_HC1 = 
//   kdeg_CYP3A4_liver * 
//   (1 + Emax_CYP3A4_RIF * fH * CHC1 / (fH * CHC1 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC1);
// dxdt_CYP3A4_ratio_HC2 = 
//   kdeg_CYP3A4_liver * 
//   (1 + Emax_CYP3A4_RIF * fH * CHC2 / (fH * CHC2 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC2);
// dxdt_CYP3A4_ratio_HC3 = 
//   kdeg_CYP3A4_liver *
//   (1 + Emax_CYP3A4_RIF * fH * CHC3 / (fH * CHC3 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC3);
// dxdt_CYP3A4_ratio_HC4 = 
//   kdeg_CYP3A4_liver * 
//   (1 + Emax_CYP3A4_RIF * fH * CHC4 / (fH * CHC4 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC4);
// dxdt_CYP3A4_ratio_HC5 = 
//   kdeg_CYP3A4_liver * 
//   (1 + Emax_CYP3A4_RIF * fH * CHC5 / (fH * CHC5 + EC50_u_CYP3A4_RIF) - CYP3A4_ratio_HC5); 



//  Hepatic extracellular compartments (HE_i; i = 1~5)
//  (VHE / 5) * dCHE_i / dt = 
// y52' = 5 / Vi * ( 
//  Qhart  * y51 + 
//  Qvilli * y68 + 
//  Qser * y62 / ( KpScaleRif * KpserRif ) - 
//  Qh * y52 + 
//  ( fhRif * PSdifeffRif * y53 - fbRif * 
//  ( VmaxOATPRif / ( freeKmOATPRif + fbRif * y52 ) + PSdifinfRif ) * y52 ) / 5 ),
dxdt_CHE1 = 
  Qhart  * Ccentral + 
  Qvilli * Cmucblood + 
  Qserosa * Cserosa / (SFKp * Kp_serosa) - 
  Qh * CHE1 + 
  (fH * PSdif_eff * CHC1 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE1) + PSdif_inf) * CHE1) / 5.0;   

dxdt_CHE1 = dxdt_CHE1 * (5.0/VHE);

//  (VHE / 5) * dCHE_i / dt = 
dxdt_CHE2 = 
  Qh * (CHE1 - CHE2) + 
  (fH * PSdif_eff * CHC2 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE2) + PSdif_inf) * CHE2) / 5.0;// (i = 2~5)

dxdt_CHE2 = dxdt_CHE2 * (5.0/VHE);

dxdt_CHE3 = 
  Qh * (CHE2 - CHE3) + 
  (fH * PSdif_eff * CHC3 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE3) + PSdif_inf) * CHE3) / 5.0;// (i = 2~5)

dxdt_CHE3 = dxdt_CHE3 * (5.0/VHE);

dxdt_CHE4 = 
  Qh * (CHE3 - CHE4) + 
  (fH * PSdif_eff * CHC4 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE4) + PSdif_inf) * CHE4) / 5.0;// (i = 2~5)

dxdt_CHE4 = dxdt_CHE4 * (5.0/VHE);

// y60' = 5 / Vi * ( Qh * ( y58 - y60 ) + 
// ( fhRif * PSdifeffRif * y61 - fbRif * 
// ( VmaxOATPRif / ( freeKmOATPRif + fbRif * y60 ) + PSdifinfRif ) * y60 ) / 5 ),
dxdt_CHE5 = 
  Qh * (CHE4 - CHE5) + 
  (fH * PSdif_eff * CHC5 - 
   fB * (Vmax_uptake / (Km_u_uptake + fB * CHE5) + PSdif_inf) * CHE5) / 5.0;// (i = 2~5)

dxdt_CHE5 = dxdt_CHE5 * (5.0/VHE);

// Hepatocytes compartments (HC_i; i = 1~5) 
//  (VHC / 5) * dCHC_i / dt = 
// y53' = 5 / Vh * ( 
//  fbRif * ( VmaxOATPRif / ( freeKmOATPRif + fbRif * y52 ) + 
// PSdifinfRif ) * y52 - fhRif * PSdifeffRif * y53 - 
// fhRif * CLintRif * ( 1 + fmUGThRif * ( y81 - 1 ) ) * y53 ) / 5,

dxdt_CHC1 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE1) + PSdif_inf) * CHE1 - 
   fH * PSdif_eff * CHC1 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC1 - 1)) * CHC1) / 5.0;

dxdt_CHC2 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE2) + PSdif_inf) * CHE2 - 
   fH * PSdif_eff * CHC2 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC2 - 1)) * CHC2) / 5.0;

dxdt_CHC3 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE3) + PSdif_inf) * CHE3 - 
   fH * PSdif_eff * CHC3 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC3 - 1)) * CHC3) / 5.0;

dxdt_CHC4 = 
  (5.0/VHC) * 
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE4) + PSdif_inf) * CHE4 - 
   fH * PSdif_eff * CHC4 -
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC4 - 1)) * CHC4) / 5.0;

dxdt_CHC5 = 
  (5.0/VHC) *
  (fB * (Vmax_uptake / (Km_u_uptake + fB * CHE5) + PSdif_inf) * CHE5 - 
   fH * PSdif_eff * CHC5 - 
   fH * CLint * (1 + fm_UGT_liver * (UGT_ratio_HC5 - 1)) * CHC5) / 5.0;

