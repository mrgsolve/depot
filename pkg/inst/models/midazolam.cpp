[ PROB ]

1: Asaumi R, Toshimoto K, Tobe Y, Hashizume K, Nunoya KI, Imawaka H, Lee W,
Sugiyama Y. Comprehensive PBPK Model of Rifampicin for Quantitative Prediction of
Complex Drug-Drug Interactions: CYP3A/2C9 Induction and OATP Inhibition Effects. 
CPT Pharmacometrics Syst Pharmacol. 2018 Jan 25. doi: 10.1002/psp4.12275. [Epub
ahead of print] PubMed PMID: 29368402.

https://www.ncbi.nlm.nih.gov/pubmed/29368402
http://onlinelibrary.wiley.com/doi/10.1002/psp4.12275/abstract

[ PARAM ]
mSFKp          = 0.201 // 1
mCLperm_gut_kg = 0.151


[ PARAM ]

mfB = 0.0545
mFa = 1.000

[ PARAM ]

mKp_liver   = 6.96
mKp_muscle  = 4.00
mKp_skin    = 20.4
mKp_adipose = 34.4

[PARAM]
mfBCLint_kg    = 0.469 // 0.528 // L/h/kg
mfECLint_E_kg  = 0.107 // L/h/kg
mCLrenal       = 0.000

[ PARAM ]
Qvilli_kg   = 0.257 // L/h/kg
Qh_kg       = 1.240 // L/h/kg
Qmuscle_kg  = 0.642 // L/h/kg
Qskin_kg    = 0.257 // L/h/kg
Qadipose_kg = 0.223 // L/h/kg
Qserosa_kg  = 0.274 // L/h/kg
Qportal_kg  = 0.531 // L/h/kg

[ PARAM ]
VHE_kg       = 0.0067  // L/kg (Vi)
VHC_kg       = 0.0174  // L/kg (Vh)
Vcentral_kg  = 0.0743  // L/kg (VbRif)
mVcentral_kg = 0.571   // L/kg (estimated)
Vskin_kg     = 0.111   // L/kg
Vadipose_kg  = 0.143   // L/kg
Vmuscle_kg   = 0.429   // L/kg
Vserosa_kg   = 0.00893 // L/kg
Vent_kg      = 0.00739 // L/kg
Vmucblood_kg = 0.00099 // L/kg
Vportal_kg   = 0.001   // L/kg

[PARAM]
mka = 1.29  // 5.51 per hr
WT  = 80    // Not sure 

[ TABLE ] 
capture Cmidazolam = 1000*mCcentral;


[ MAIN ]
if(NEWIND <= 1) {
  // -------------------------------------------
  double mfBCLint    = mfBCLint_kg*WT;
  double mCLperm_gut = mCLperm_gut_kg*WT;
  double mfECLint_E  = mfECLint_E_kg*WT;
  double Qvilli      = Qvilli_kg*WT;
  double Qh          = Qh_kg*WT;
  double Qmuscle     = Qmuscle_kg*WT;
  double Qskin       = Qskin_kg*WT;
  double Qadipose    = Qadipose_kg*WT;
  double Qserosa     = Qserosa_kg*WT;
  double Qhart       = Qh - Qserosa - Qvilli;
  double Qportal     = Qportal_kg*WT;
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
  double mVcentral = mVcentral_kg*WT;
  double Vportal = Vportal_kg*WT;
  // -------------------------------------------
  double mQgut = Qvilli * mCLperm_gut / (Qvilli + mCLperm_gut);
}

[CMT]
Mgutlumen
mcentral mCmuscle mCskin mCadipose
CLIV1 CLIV2 CLIV3 CLIV4 CLIV5  
Cportal 

[CAPTURE] mCcentral = (1000*mcentral/mVcentral)
  
[ ODE ]

double mCcentral = mcentral/mVcentral;

dxdt_mcentral = 
  Qh * (CLIV5 / (mSFKp * mKp_liver)) - 
  (Qh-Qportal) * mCcentral +
  Qmuscle      * (mCmuscle  / (mSFKp * mKp_muscle)  - mCcentral) + 
  Qskin        * (mCskin    / (mSFKp * mKp_skin)    - mCcentral) +
  Qadipose     * (mCadipose / (mSFKp * mKp_adipose) - mCcentral) -
  Qportal      * mCcentral - 
  mCLrenal     * mCcentral;

dxdt_CLIV1 = 
  (Qh-Qportal) * mCcentral + 
  Qportal * Cportal - 
  Qh * CLIV1 / (mSFKp * mKp_liver) - 
  mfBCLint / 5.0 * CLIV1 / 
  (mSFKp * mKp_liver);

dxdt_CLIV1 = dxdt_CLIV1 * (5/(VHE+VHC));
  
dxdt_CLIV2 = 
  (Qh * (CLIV1 - CLIV2) - 
   mfBCLint / 5.0 * CLIV2) / 
   (mSFKp * mKp_liver); 
  
dxdt_CLIV2 = dxdt_CLIV2 * (5/(VHE+VHC));   
 
dxdt_CLIV3 = 
  (Qh * (CLIV2 - CLIV3) - 
  mfBCLint / 5.0 * CLIV3) / 
  (mSFKp * mKp_liver); 
 
dxdt_CLIV3 = dxdt_CLIV3 * (5/(VHE+VHC));   
 
dxdt_CLIV4 = 
  (Qh * (CLIV3 - CLIV4) - 
  mfBCLint / 5.0 * CLIV4) / 
  (mSFKp * mKp_liver); 
 
dxdt_CLIV4 = dxdt_CLIV4 * (5.0/(VHE+VHC));    

dxdt_CLIV5 = 
  (Qh * (CLIV4 - CLIV5) - 
  mfBCLint / 5.0 * CLIV5) / 
  (mSFKp * mKp_liver); 

dxdt_CLIV5 = dxdt_CLIV5 * (5.0/(VHE+VHC));    

dxdt_Cportal = 
  Qportal * (mCcentral - Cportal) + 
  mka * mQgut / (mQgut + mfECLint_E) * Mgutlumen; 

dxdt_Cportal = dxdt_Cportal * (1.0/Vportal); 
  
dxdt_mCmuscle = 
  (1.0/Vmuscle) * Qmuscle * (mCcentral - mCmuscle / (mSFKp * mKp_muscle));

dxdt_mCskin = 
  (1.0/Vskin) * Qskin * (mCcentral - mCskin / (mSFKp * mKp_skin));

dxdt_mCadipose = 
  (1.0/Vadipose) * Qadipose * (mCcentral - mCadipose / (mSFKp * mKp_adipose));

dxdt_Mgutlumen = -mka/mFa * Mgutlumen;

