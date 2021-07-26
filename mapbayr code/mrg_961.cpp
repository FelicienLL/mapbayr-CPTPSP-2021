$PROB Voriconazole, Liu et al, AAC 2014 [ORAL !]
//unit: vol (L), qte (microg), time (h)


$PARAM @annotated
TVKM      : 1150  : KM constant (microg/L) // original 1.15 microg/mL
TVVMAX    : 113   : VMAX constant (microg/h/70kg)
VMAXINH   : 1.50  : Maximum fract of Vmax inhibition ()
VMAXSCALE : 0.583 : Scaling on random effect ()
T50       : 2.42  : Time of half inhibition (h)
TVCL      : 5.30  : Clearance (L/h/70kg)
TVV2      : 77.6  : Central volume (L/70kg)
TVV3      : 89.5  : Peripheral volume (L/70kg)
TVQ       : 15.9  : Intercompartmental clearance (L/h/70kg)
TVF1      : 0.595 : Typical bioavailability ()
BCF1      : 0.367 : Box-Cox parameter F1 ()
KA        : 1.2   : Absorption rate (h-1)
ALAG1     : 1     : Lag time (h)
TVRATE2   : 12800  : Typical rate (microg/h)
// rate; we assume prior treatment by voriconazole


ETA1 : 0 : VMAX_KM
ETA2 : 0 : CL
ETA3 : 0 : V2
ETA4 : 0 : V3
ETA5 : 0 : Q
ETA6 : 0 : F1
ETA7 : 0 : R2

$PARAM @annotated @covariates
WT : 70 : Body weight (kg)

$OMEGA
1.91  // VMAX_KM
0.634 // CL
0.139 // V2
0.831 // V3
0.459 // Q
0.713 // F1
0.910 // R2

$SIGMA
0 // err prop
0.017424   //  err log additive 0.132 au carre

$CMT @annotated
DEPOT   : Depot () [ADM]
CENTRAL : Central () [OBS]
PERIPH  : Peripheral () []

$TABLE
double DV = CONC * exp(EPS(2)) ; //dose in microg, vol in Liter = PRED in microg/L

$MAIN
double KM = TVKM * exp(ETA1 + ETA(1)) ;
double iVMAX = TVVMAX * exp(VMAXSCALE * (ETA1 + ETA(1))) * pow(WT/70, 0.75) ;
double iVMAXINH  = exp(VMAXINH) / (1 + exp(VMAXINH)) ;
double CL = TVCL * exp(ETA2 + ETA(2)) * pow(WT/70, 0.75) ;
double V2 = TVV2 * exp(ETA3 + ETA(3)) * WT/70 ;
double V3 = TVV3 * exp(ETA4 + ETA(4)) * WT/70 ;
double Q  = TVQ  * exp(ETA5 + ETA(5)) * pow(WT/70, 0.75) ;
double ETATR = (exp((ETA(6) + ETA6) * BCF1) - 1) / BCF1 ;
double LTVF1 = log(TVF1 / (1 - TVF1)) ;
double F1 = exp(LTVF1 + ETATR) / (1 + exp(LTVF1 + ETATR)) ;
double RATE2 = TVRATE2 * exp(ETA7 + ETA(7)) ;

double K20 = CL / V2 ;
double K23 = Q / V2 ;
double K32 = Q / V3 ;

F_DEPOT = F1 ;
ALAG_DEPOT = ALAG1 ;

$ODE
double VMAX = iVMAX * (1 - iVMAXINH * ((SOLVERTIME - 1)/ ((SOLVERTIME-1) + (T50 - 1)))) ;
double CONC  = CENTRAL / V2 ;

dxdt_DEPOT   = - KA * DEPOT ;
dxdt_CENTRAL =   KA * DEPOT - (K20 + K23) * CENTRAL + K32 * PERIPH - (VMAX * CONC)/(CONC + KM) + RATE2;
dxdt_PERIPH  = K23 * CENTRAL - K32 * PERIPH ;

$CAPTURE DV
